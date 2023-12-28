program main_CDR3d
  use param; use geometry; use biunit; use inputInfo
  use library; use boundVal; use timeInt; use sourceTerm

implicit none

  external :: fdate
  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  character(len=19)                                 :: name_inputFile
  character(len=13)                                 :: geometry_File
  double precision, allocatable, dimension(:,:)     :: A_K, A_C, A_F, presc !,Jsource
  double precision, allocatable, dimension(:,:)     :: basfun, dN_dxi, dN_deta
  double precision, allocatable, dimension(:,:)     :: hes_xixi, hes_xieta, hes_etaeta
  !double precision,              dimension(3,4) :: Hesxieta
  double precision, allocatable, dimension(:,:)     :: Bvs
  integer,          allocatable, dimension(:,:)     :: condition, ifpre
  integer,          allocatable, dimension(:)       :: nofix
  double precision                                              :: start, finish, k_y
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
  double precision, allocatable, dimension(:,:,:)   :: grad_u_sol
  double precision, allocatable, dimension(:,:)     :: AK_LU, u_sol
  double precision, allocatable, dimension(:)       :: Ex_field, ky
  external                                          :: dgbtrf, dgbtrs, dgbrfs
  integer,          allocatable, dimension(:)       :: S_ipiv!, S_iwork
  character(len=1)                                  :: S_trans
  character(len=24)           :: date
  integer                                           :: S_m, S_n, S_nrhs, info, S_ldSol, ii
  
  !-----Input File
  call cpu_time(start)
  name_inputFile = 'tMaxwelinputCDR.c' 
  
  !--------------- Input Data ---------------!
  call inputData(name_inputFile, geometry_file)
  
  !--------------- Geometry -----------------!
  call readMesh(geometry_File)
  call GeneralInfo(name_inputFile, geometry_File)

  !---------- Shape Functions ---------------!
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta)

  !------- Computing half bandwidth  --------!
  call BandWidth

  !------- Setting Boundary Conditions ------!
  !Esta funcion crea el archivo BVs.dat
  call SetBoundVal( nBVs, nBVscol)
  !Designo memoria para la matriz de nodos con valor en la frontera
  allocate( condition(nBvs,nBVscol-ndofn), BVs(nBVs,ndofn) ) 
  call ReadFile(nBVs, nBVscol, condition, BVs)!Reading the boundary values file and store it at BVs
  !memoria para el numero de nodo, si esta preescrito o no y su valor correspondiente
  allocate( nofix(nBVs), ifpre(ndofn,nBVs), presc(ndofn,nBVs) )
  call VinculBVs( condition,  BVs, nofix, ifpre, presc )

  !----- Setting the wave numbers range ------!
  call WaveNumbers(1.0d-7,4.0d-2,10,ky) !k**2 = i*mu*sigma*2pi*f
  k_y = ky(idk_y)

  !-------- Problem Type Definition ----------!
  select case(ProbType)
    case('TIME')
      
      print*, ' '
      print*, "!=============== Time Integration =============! "
      call Timeintegration(k_y, basfun, dN_dxi, dN_deta,hes_xixi,hes_xieta,hes_etaeta,&
      &                    nofix, ifpre, presc, Ex_field)
     
      ! call Res_Matlab(Ex_field)
     
      !---------- Memory Relase -----------!
      deallocate( basfun, dN_dxi, dN_deta, BVs, nofix, ifpre, presc)
      deallocate(hes_xixi, hes_xieta, hes_etaeta) 
    case('STAT')
      !---------- Global Matrix and Vector ------!
      print*,'  •BUILDING GLOBAL MATRIX (K) AND VECTOR (F).....'
      !the allocate of A_K and A_F are inside of GlobalSystem
      call GlobalSystem(k_y,basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      call ApplyBVs(nofix,ifpre,presc,A_K, A_F)
      
      !Applying the source term for DC simulation  j = I*δ(x-x0)δ(y-y0)*δ(z-z0)
      if(BCsProb.eq.5)then
        continue
      else
        if(((ndofn.eq.1).or.(ndofn.eq.3)).and.((exacSol.eq.5).or.(exacSol.eq.2)))then
          print*,'delta source'
          if(ndofn.eq.1)then
            do ii=1,nodalSrc
              if(ii.eq.2)Icurr(1) = -1.0*Icurr(1)
              !print*,Icurr
              A_F((srcLoc(ii)-1)*ndofn+1,1) = -Icurr(1)
            end do
          else
            do ii=1,nodalSrc
              if(ii.eq.2)Icurr(1) = -1.0*Icurr(1)
              !print*,Icurr
              A_F((srcLoc(ii)-1)*ndofn+1,1) = -Icurr(1)
              A_F((srcLoc(ii)-1)*ndofn+2,1) = -Icurr(2)
            end do
          endif
        end if
      endif
      
      !----- Setting MKL-Solver Parammeters -----!
      S_m     = size(A_K,2)  !antes ntotv
      S_n     = size(A_K,2)  !antes ntotv
      S_ldSol = max(1,S_n)
      S_trans = 'N'
      S_nrhs  = 1
      
      allocate( AK_LU(ldAKban,ntotv), u_sol(S_ldSol,1) ) 
      allocate( S_ipiv(max(1,min(S_m, S_n)) ) )
      
      print*, ' '
      print*, '!================ MKL Solver ==================!'
      print*,'  •SETTING VARIABLES FOR MKL LIBRARY.....'
      AK_LU = A_K     !AK_band(ldab,*) The array AK_band contains the matrix A_K in band storage
      u_sol = A_F     !Sol_vec will be rewrited by LAPACK solution avoiding lose A_F
      
      !---------- Solving System of Equations by INTEL MKL solver -----------!
      print*,'  •PERFORMING LU BAND DECOMPOSITION.....'
      call dgbtrf(S_m, S_n, lowban, upban, AK_LU, ldAKban, S_ipiv, info)  
      call MKLfactoResult('dgbtrf',info) 
      print*,'  •SOLVING SYSTEM OF EQUATIONS.....'
      
      !---------- Computing nodal unknowns ux_i, uy_i and p_i ---------------!
      call dgbtrs( S_trans, S_n, lowban, upban, S_nrhs, AK_LU, ldAKban, S_ipiv, u_sol, S_ldSol, info )
      call MKLsolverResult('dgbtrs',info)
   
      !---------- Computing Inverse Fourier Transform------------------------!
      ! call iFT(ky, u_sol, 3D_sol)

      !---------- Computing E=-∇u or -∂B/∂t=∇xE -----------------------------!
      call PostPro_EMfield(basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,u_sol,grad_u_sol)
      
      !---------- Print and write results -----------------------------------!
      call GID_results(u_sol, grad_u_sol) 
      !call Res_Matlab(u_sol)
      
      write(*,*)
      !---------- Memory Relase ---------------------------------------------!
      DEALLOCATE( A_F, A_K, AK_LU, u_sol)
      deallocate( basfun, dN_dxi, dN_deta, BVs,  nofix, ifpre, presc)
      deallocate(hes_xixi, hes_xieta, hes_etaeta) 
      
    case DEFAULT
      PRINT*, ' ' 
      print*, 'In Problem Type Transient or Static'
      write(*,*) 'Invalid Problem Type.'
      stop
  end select

  call cpu_time(finish); call fdate(date); 
  write(*,'(A11,f14.5,A,f5.2,A,A)   ')' &
    &CPU-Time =', finish,' ~',((finish-start)/60.0)/10.0,' minutes. Finished on: ', date
  print*,' ' 

end program main_CDR3d


! print*, 'Antes de BV'
! write(*,"(I2,1x, 16F8.5)") (i, (A_K(i,j), j=1,ntotv), i=1,ldAKban)
! do i =1, ntotv
!   print '(I2, 1x, E13.5)',i, A_F(i,1)
! end do
   
! print*, 'Despues de BV'
! write(*,"(I2,1x, 16F8.5)") (i, (AK_LU(i,j), j=1,ntotv), i=1,ldAKban)
! do i =1, ntotv
!   print '(I2,1X,F9.5)', i, u_sol(i,1)
! end do
