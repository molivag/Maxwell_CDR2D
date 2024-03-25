program main_CDR3d
  use param; use geometry; use biunit; use inputInfo
  use library; use boundVal; use timeInt; use sourceTerm

  implicit none

  external :: fdate
  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  character(len=19)                                 :: name_inputFile
  character(len=:), allocatable                     :: geometry_File
  double precision, allocatable, dimension(:,:)     :: A_K, A_C, A_F, presc
  double precision, allocatable,  dimension(:,:)    :: Jsource
  double precision, allocatable, dimension(:,:)     :: basfun, dN_dxi, dN_deta
  double precision, allocatable, dimension(:,:)     :: hes_xixi, hes_xieta, hes_etaeta
  !double precision,              dimension(3,4) :: Hesxieta
  double precision, allocatable, dimension(:,:)     :: Bvs
  integer,          allocatable, dimension(:,:)     :: condition, ifpre
  integer,          allocatable, dimension(:)       :: nofix
  double precision                                              :: start, finish
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
  double precision, allocatable, dimension(:,:,:)   :: grad_u_sol
  double precision, allocatable, dimension(:,:)     :: AK_LU, u_sol, E_3D
  double precision, allocatable, dimension(:)       :: Ex_field
  external                                          :: dgbtrf, dgbtrs, dgbrfs
  integer,          allocatable, dimension(:)       :: S_ipiv!, S_iwork
  character(len=1)                                  :: S_trans
  character(len=24)           :: date
  integer                                           :: S_m, S_n, S_nrhs, info, S_ldSol, ii_ky, ii
  double precision, allocatable, dimension(:,:)     :: store_Spec
  
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


  !====================== Hasta aqui iria la parte general de todas las simulaciones para cada ky
  !-------- Problem Type Definition ----------!
  print*, ' '
  print'(A)', " !=============== Computation of 2D problems for each wave number =============! "
  select case(ProbType)
    case('TIME')
      ! varPrimerSet = 10 -5 = 5
      ! varSeconSet = varPrimerSet+1 = 6
      if(splits == 'N')then
        ! do ii_ky = 1,tot_ky-5 !First running of the code with tot_ky = 5
        do ii_ky = 1,tot_ky 
        call TimeIntegration(ii_ky, basfun, dN_dxi, dN_deta,hes_xixi,hes_xieta,hes_etaeta, nofix, ifpre, presc, Ex_field)
        end do
        !---------- Memory Relase -----------!
        deallocate( basfun, dN_dxi, dN_deta, BVs, nofix, ifpre, presc)
        deallocate(hes_xixi, hes_xieta, hes_etaeta) 
        if((ii_ky==size(nodal_ky)+1) .and. (TwoHalf=='Y')) call invDFT(E_3D) !if para test cortos pero sin dividir la ejecucion
        
      elseif(splits == 'Y')then !If the problem split run as kind of parallel
        print*,'dividido'
        do ii_ky = tot_ky-4, tot_ky! this do for the case of 10 ky, will go from 6 to 10
          call TimeIntegration(ii_ky, basfun, dN_dxi, dN_deta,hes_xixi,hes_xieta,hes_etaeta, nofix, ifpre, presc, Ex_field)
        end do
        !---------- Memory Relase -----------!
        deallocate( basfun, dN_dxi, dN_deta, BVs, nofix, ifpre, presc)
        deallocate(hes_xixi, hes_xieta, hes_etaeta) 
        call invDFT(E_3D) 
        
      endif
      
    case('STAT')
      !---------- Global Matrix and Vector ------!
      ! if(ProbType == 'STAT') t_steps=0
      t_steps=0
      print*,'  •BUILDING GLOBAL MATRIX (K) AND VECTOR (F).....'
      if_run_parallel: if(splits == 'N')then
        loop_each_ky: do ii_ky = 1,tot_ky !First running of the code with tot_ky = 5
          if(TwoHalf == 'Y')then !if it is a 2.5D problem, these variable gonna be updating 
            ky_id = nodal_ky(ii_ky)
            k_y = WaveNumbers(ii_ky) !Value of wave number for current problem (used in reama)
            print'(A17,I0,A3,E11.4)', " !-------> for ky",ii_ky,"=",k_y
          else
            continue
          endif
          !the allocation of A_K and A_F are inside of GlobalSystem
          call GlobalSystem(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
          call ApplyBVs(nofix,ifpre,presc,A_K, A_F)
          call currDensity(Jsource) 
          Jsource = Jsource+A_F
          !----- Setting MKL-Solver Parammeters -----!
          S_m     = size(A_K,2)  !antes ntotv
          S_n     = size(A_K,2)  !antes ntotv
          S_ldSol = max(1,S_n)
          S_trans = 'N'
          S_nrhs  = 1
          allocate( store_Spec(S_ldSol,t_steps+1) )
          allocate( AK_LU(ldAKban,ntotv), u_sol(S_ldSol,1) ) 
          allocate( S_ipiv(max(1,min(S_m, S_n)) ) )
          print*, ' '
          print*, '!================ MKL Solver ==================!'
          print*,'  •SETTING VARIABLES FOR MKL LIBRARY.....'
          AK_LU = A_K      !AK_band(ldab,*) contains the matrix A_K in band storage
          u_sol = Jsource  !u_sol gonna be rewrited by LAPACK avoiding lose Jsource
          !---------- Solving System of Equations by INTEL MKL solver -----------!
          print*,'  •PERFORMING LU BAND DECOMPOSITION.....'
          call dgbtrf(S_m, S_n, lowban, upban, AK_LU, ldAKban, S_ipiv, info)  
          call MKLfactoResult('dgbtrf',info) 
          print*,'  •SOLVING SYSTEM OF EQUATIONS.....'
          !---------- Computing nodal unknowns ux_i, uy_i and p_i ---------------!
          call dgbtrs(S_trans,S_n,lowban,upban,S_nrhs,AK_LU,ldAKban,S_ipiv,u_sol,S_ldSol,info )
          call MKLsolverResult('dgbtrs',info)
          
          !!---------- Computing E=-∇u or -∂B/∂t=∇xE -----------------------------!
          !call PostPro_EMfield(basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,u_sol,grad_u_sol)
          
          if(TwoHalf == 'Y')then
            do ii = 1, S_ldSol
              store_Spec(ii,1) = u_sol(ii,1) 
            end do
            call storeSpectrum(ii_ky, store_Spec)
          endif 
          
          !---------- Print and write results -----------------------------------!
          ! u_sol = log10(u_sol)
          call GID_results(u_sol) 
          !---------- Memory Relase ---------------------------------------------!
          DEALLOCATE( A_F, A_K, AK_LU)
          if(TwoHalf=='Y')DEALLOCATE( u_sol)
          DEALLOCATE( store_Spec, S_ipiv )
          
        end do loop_each_ky
        deallocate(BVs,  nofix, ifpre, presc)
        !---------- Performing the Inverse Fourier Transform -----------------------------!
        if(TwoHalf=='Y')call invDFT(E_3D) 
        
      endif if_run_parallel
      if(TwoHalf =='N')goto 10
      !---------- Computing E=-∇u or -∂B/∂t=∇xE -----------------------------!
      call PostPro_EMfield(basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,E_3D,grad_u_sol)
      deallocate(hes_xixi, hes_xieta, hes_etaeta) 
      deallocate( basfun, dN_dxi, dN_deta)
      call spatialProfile_BubbleSort(E_3D)
      TwoHalf = 'N'
      File_Nodal_Vals = File_3DNodal_Vals
      ! E_3D = log10(E_3D)
      call GID_results(E_3D, grad_u_sol) 
      !call Res_Matlab(u_sol)
      DEALLOCATE( E_3D)
      
    case DEFAULT
      PRINT*, ' ' 
      print*, 'In Problem Type Transient or Static'
      write(*,*) 'Invalid Problem Type.'
      stop
    end select
    
  10 write(*,*)
  ! u_sol = 10.0**(u_sol)
  print*, 'Writing spatial profile at z=0'
  ! call spatialProfile_BubbleSort(u_sol)

  call cpu_time(finish); call fdate(date)
  print*, ' '
  print'(A)',' - - -' 
  write(*,'(A11,f14.5,A,f5.2,A,A)   ')' &
    &CPU-Time =', finish,' ~',((finish-start)/60.0)/10.0,' minutes. Finished on: ', date
  ! primero entre 10 p[ara pasar de cientos de segundos a decenas de segundos y luego dividir esas decenas
  ! en segmentos de sesenta  segundos que serian un minuto, por eso de ambas divisiones
  write(*,*) 'Computation finished. Run-time is ', finish-start, 'seconds.'
  write(*,*) 'Esto es en minutos',  ((finish-start)/60.0)
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
