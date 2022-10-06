program main_CDR3d
  use library; use param ; use biunit ; use boundVal ; use timeInt

implicit none

  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  double precision, allocatable, dimension(:,:) :: A_K, A_C, A_F, presc
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision, allocatable, dimension(:,:) :: hes_xixi, hes_xieta, hes_etaeta
  !double precision,              dimension(3,4) :: Hesxieta
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
  external :: dgbtrf, dgbtrs, dgbrfs
  double precision, allocatable, dimension(:,:) :: AK_LU, u_sol
  !double precision, allocatable, dimension(:)   :: S_ferr, S_berr, S_work
  integer, allocatable,          dimension(:)   :: S_ipiv!, S_iwork
  character(len=1) :: S_trans
  integer          :: S_m, S_n, S_nrhs, info, S_ldSol, ans!, workdim
  

  !---------- Input Data -----------!
  call cpu_time(start)
  call inputData( )
  call GeneralInfo( )

  !---------- Shape Functions ---------------!
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta)

  !------- Computing half bandwidth  --------!
  call BandWidth( )
  
  !---------- Global Matrix and Vector ------!
  call GlobalSystem(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
  !the allocate of A_K and A_F are inside of GlobalSystem

  !------- Setting Boundary Conditions ------!
  call SetBoundVal( nBVs, nBVscol) !Esta funcion crea el archivo BVs.dat
  allocate( BVs(nBVs, nBVscol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadFile("BVs.dat", nBVs, nBVscol, BVs)!Reading the boundary values file and store it at BVs
  !memoria para el numero de nodo, si esta preescrito o no y su valor correspondiente
  allocate( nofix(nBVs), ifpre(ndofn,nBVs), presc(ndofn,nBVs) )
  call VinculBVs(  BVs, nofix, ifpre, presc )

  !----- Setting MKL-Solver Parammeters -----!
  S_m     = size(A_K,2)  !antes ntotv
  S_n     = size(A_K,2)  !antes ntotv
  S_ldSol = max(1,S_n)
  S_trans = 'N'
  S_nrhs  = 1
  !workdim = max(1,3*S_n)


  !-------- Problem Type Definition --------!
  if(ProbType .eq. 'trans')then
    
    call TimeIntegration(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta,&
    &                    time_ini, time_fin, max_time, nofix, ifpre, presc, &
    &                    S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol)
   
    !---------- Memory Relase -----------!
    DEALLOCATE( N, dN_dxi, dN_deta, BVs, nofix, ifpre, presc)
    DEALLOCATE(hes_xixi, hes_xieta, hes_etaeta) 
  
  else
    
    call ApplyBVs(nofix,ifpre,presc,A_K, A_F)
    
    !---------- Memory Relase -----------!
    DEALLOCATE( N, dN_dxi, dN_deta, BVs,  nofix, ifpre, presc)
    DEALLOCATE(hes_xixi, hes_xieta, hes_etaeta) 
    
    allocate( AK_LU(ldAKban,ntotv), u_sol(S_ldSol,1)) 
    allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
    
    !allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )

    print*, ' '
    print*, '!================ MKL Solver ==============!'
    AK_LU = A_K     !AK_band(ldab,*) The array AK_band contains the matrix A_K in band storage
    u_sol = A_F     !Sol_vec will be rewrited by LAPACK solution avoiding lose A_F
   
    print*, ' '
    
    !---------- Solving System of Equations by retpla solver -----------!
    print*,'  •INITIALIZING BAND LU DECOMPOSITION.....'                                                    
    call dgbtrf(S_m, S_n, lowban, upban, AK_LU, ldAKban, S_ipiv, info)  
    call MKLfactoResult('dgbtrf',info) 
    print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
    call dgbtrs( S_trans, S_n, lowban, upban, S_nrhs, AK_LU, ldAKban, S_ipiv, u_sol, S_ldSol, info )
    call MKLsolverResult('dgbtrs',info)
    print*, ' '
    
    !---------- Print and write results -----------!
    print*, '!====== Name of Postprocess file'
    read(*,*) File_PostProcess 
    call PosProcess(u_sol, 'msh') 
    call PosProcess(u_sol, 'res')
    write(*,*)'Matlab File?. Y=1, N=2'; read(*,*) ans
    if(ans.eq.1)call Res_Matlab(u_sol)
    
    !print*, ' '
    !print*, 'Shape of Global K: ',shape(A_K)
    !print*, 'Shape of Global F: ',shape(A_F)
    !print*, 'Shape of Solution: ',shape(u_sol)
    write(*,*)
    !---------- Memory Relase -----------!
    DEALLOCATE( A_F, A_K, AK_LU, u_sol)
    
  endif

  call cpu_time(finish)
  write(*,'(A11,f9.2,A8)')' CPU-Time =', finish-start, ' Seconds', ' '

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
