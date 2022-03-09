program main_CDR3d
  use library; use param ; use biunit ; use boundVal ; use timeInt
  !use solver

implicit none

  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  double precision, allocatable, dimension(:,:) :: A_K, A_F, presc
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision,              dimension(3,4) :: Hesxieta
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish, time_ini, time_fin, max_time, u0_cond
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
  external :: dgbtrf, dgbtrs, dpbtrs, dpbtrf
  double precision, allocatable, dimension(:,:) :: AK_LU, u_sol
  integer, allocatable, dimension(:) :: S_ipiv
  integer :: S_m, S_n, S_nrhs, info, S_ldSol
  character(len=1) :: S_trans, S_uplo

  !---------- Input Data -----------!
  call cpu_time(start)
  call inputData( )
  call GeneralInfo( )

  !---------- Shape Functions ---------------!
  !print*, '!=============== INFO DURING EXECUTION ===============!'
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, nne, N, dN_dxi, dN_deta, Hesxieta)

  !------- Computing half bandwidth  --------!
  call BandWidth( )
  
  !---------- Global Matrix and Vector ------!
  !the allocate of A_K is inside of GlobalSystem, first compute the semi bandwidth, then allocate A_K
  call GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_K, A_F)
  !print*, 'Antes de BV'
  !write(*,"(I2,1x, 16F13.9)") (i, (A_K(i,j), j=1,ntotv), i=1,ldAKban)
  !print*,'!=============== Output Files ================!'
  !call writeMatrix(A_K, 10, 'A_K.dat', A_F, 20, 'A_F.dat')
  !call writeMatrixAKbLU,60,'-', Sols, 70, 'SolMKL_LU.dat')

  !------- Setting Boundary Conditions ------!
  call SetBoundVal( nBVs, nBVscol) !Esta funcion crea el archivo BVs.dat
  allocate( BVs(nBVs, nBVscol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadIntegerFile(30,"BVs.dat", nBVs, nBVscol, BVs)!Reading the boundary values file and it will store at BVs
  allocate( nofix(nBVs), ifpre(ndofn,nBVs), presc(ndofn,nBVs) ) !memoria para el numero de nodo, si esta preescrito y su valor
  call VinculBVs(  BVs, nofix, ifpre, presc )

  !----- Setting MKL-Solver Parammeters -----!
  S_m     = size(A_K,2)  !antes ntotv
  S_n     = size(A_K,2)  !antes ntotv
  S_ldSol = max(1,S_n)
  S_trans = 'N'
  S_uplo  = 'U'
  S_nrhs  = 1

  !-------- Problem Type Definition --------!
  if(ProbType .eq. 'trans')then
    time_ini = 0.0   !Estos valores
    time_fin = 0.5       !Deben ser leidos en
    max_time = 10.0           !el archivo de entrada
    u0_cond  = 0.0  
    
    call BackwardEuler(N, dN_dxi, dN_deta, Hesxieta, time_ini, time_fin, max_time, u0_cond,&
    &                      nofix, ifpre, presc, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol )
   
    !---------- Memory Relase -----------!
    DEALLOCATE( N, dN_dxi, dN_deta, BVs,  nofix, ifpre, presc)
    
  else
    call ApplyBVs(nofix,ifpre,presc,A_K, A_F)
    !print*, 'Despues de BV'
    !write(*,"(I2,1x, 16F13.9)") (i, (A_K(i,j), j=1,ntotv), i=1,ldAKban)
    
    !---------- Memory Relase -----------!
    DEALLOCATE( N, dN_dxi, dN_deta, BVs,  nofix, ifpre, presc)
    allocate( AK_LU(ldAKban,ntotv), u_sol(S_ldSol,1)) 
    allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
    
    print*, ' '
    print*, '!================ MKL Solver ==============!'
    AK_LU = A_K                 !AK_band(ldab,*) The array AK_band contains the matrix A_K in band storage
    u_sol = A_F                 !Sol_vec will be rewrited by LAPACK solution avoiding lose A_F
    !---------- Solving System of Equations by retpla solver -----------!
    !  print*, 'Cholesky-decomposition of band-storaged Matrix'          
    !  call solsistem(A_K,A_K,1,A_F,Sols)                                                             
    print*,'  •INITIALIZING BAND LU DECOMPOSITION.....'                                                    
    !print*, 'Cholesky-decomposition of band-storaged Matrix'     
    !call dpbtrf( uplo, S_n, nband, AKb_LU, ldAKban, info )
    call dgbtrf(S_m, S_n, lowban, upban, AK_LU, ldAKban, S_ipiv, info)  
    call MKLfactoResult('dgbtrf',info)  !Aqui agregar el paso del tiempo para en cada tiempo indicar el info de ejecucion
    print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
    !call dpbtrs( uplo, S_n, nband, 1, AKb_LU, ldAKbn, b, ldb, info )                                     
    call dgbtrs( S_trans, S_n, lowban, upban, S_nrhs, AK_LU, ldAKban, S_ipiv, u_sol, S_ldSol, info )
    call MKLsolverResult('dgbtrs',info) !Aqui agregar el paso del tiempo para en cada tiempo indicar el info de ejecucion
    !---------- Print and write results -----------!
    call PosProcess(u_sol, File_PostMsh, 'msh') !se debe agregar el nt como dummyvariable
    call PosProcess(u_sol, File_PostRes, 'res')
    !call GID_PostProcess(u_sol, File_PostMsh, File_PostRes,'msh', 0) !se debe agregar el nt como dummyvariable
    !call GID_PostProcess(u_sol, File_PostMsh, File_PostRes,'res', 0) !se debe agregar el nt como dummyvariable
    
    print*, ' '
    print*, 'Shape of Global K: ',shape(A_K)
    print*, 'Shape of Global F: ',shape(A_F)
    print*, 'Shape of Solution: ',shape(u_sol)
    write(*,*)
    !---------- Memory Relase -----------!
    DEALLOCATE( A_F, A_K, AK_LU, u_sol)
    
  endif

  call cpu_time(finish)
  write(*,'(A11,f9.2,A8)')' CPU-Time =', finish-start, ' Seconds', ' '

end program main_CDR3d
