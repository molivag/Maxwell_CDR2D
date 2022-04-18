program main_CDR3d
  use library; use param ; use biunit ; use boundVal ; use timeInt

implicit none

  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  double precision, allocatable, dimension(:,:) :: A_K, A_F, presc
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision,              dimension(3,4) :: Hesxieta
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
  external :: dgbtrf, dgbtrs, dgbrfs
  double precision, allocatable, dimension(:,:) :: AK_LU, u_sol
  double precision, allocatable, dimension(:)   :: S_ferr, S_berr, S_work
  integer, allocatable, dimension(:) :: S_ipiv, S_iwork
  integer :: S_m, S_n, S_nrhs, info, S_ldSol, workdim, S_info,  i, j
  character(len=1) :: S_trans
  

  !---------- Input Data -----------!
  call cpu_time(start)
  call inputData( )
  call GeneralInfo( )

  !---------- Shape Functions ---------------!
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, nne, N, dN_dxi, dN_deta, Hesxieta)

  !------- Computing half bandwidth  --------!
  call BandWidth( )
  
  !---------- Global Matrix and Vector ------!
  call GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_K, A_F) !the allocate of A_K and A_F are inside of GlobalSystem
 ! print*,'!=============== Output Files ================!'
 ! call writeMatrix(A_K, 10, 'A_K.dat', A_F, 20, 'A_F.dat')
 ! call writeMatrix(AKbLU,60,'-', Sols, 70, 'SolMKL_LU.dat')

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
  S_nrhs  = 1
  workdim = max(1,3*S_n)


  !-------- Problem Type Definition --------!
  if(ProbType .eq. 'trans')then
   ! time_ini = 0.0   !Estos valores
   ! time_fin = 1.0       !Deben ser leidos en
   ! max_time = 35           !el archivo de entrada
   ! u0_cond  = 0.0  
    
    call BackwardEuler(N, dN_dxi, dN_deta, Hesxieta, time_ini, time_fin, max_time, u0_cond,&
    &                      nofix, ifpre, presc, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol, workdim )
   
    !---------- Memory Relase -----------!
    DEALLOCATE( N, dN_dxi, dN_deta, BVs,  nofix, ifpre, presc)
    
  else
   ! print*, 'Antes de BV'
   ! write(*,"(I2,1x, 16F8.5)") (i, (A_K(i,j), j=1,ntotv), i=1,ldAKban)
   ! do i =1, ntotv
   !   print '(I2, 1x, E13.5)',i, A_F(i,1)
   ! end do
    call ApplyBVs(nofix,ifpre,presc,A_K, A_F)
    
    !---------- Memory Relase -----------!
    DEALLOCATE( N, dN_dxi, dN_deta, BVs,  nofix, ifpre, presc)
    allocate( AK_LU(ldAKban,ntotv), u_sol(S_ldSol,1)) 
    allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
    
    allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )

    print*, ' '
    print*, '!================ MKL Solver ==============!'
    AK_LU = A_K                 !AK_band(ldab,*) The array AK_band contains the matrix A_K in band storage
    u_sol = A_F                 !Sol_vec will be rewrited by LAPACK solution avoiding lose A_F
   ! print*, 'Despues de BV'
   ! write(*,"(I2,1x, 16F8.5)") (i, (AK_LU(i,j), j=1,ntotv), i=1,ldAKban)
   ! do i =1, ntotv
   !   print '(I2,1X,F9.5)', i, u_sol(i,1)
   ! end do
   ! print*, ' '
    
    !---------- Solving System of Equations by retpla solver -----------!
    print*,'  •INITIALIZING BAND LU DECOMPOSITION.....'                                                    
    call dgbtrf(S_m, S_n, lowban, upban, AK_LU, ldAKban, S_ipiv, info)  
    call MKLfactoResult('dgbtrf',info)  !Aqui agregar el paso del tiempo para en cada tiempo indicar el info de ejecucion
    print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
    call dgbtrs( S_trans, S_n, lowban, upban, S_nrhs, AK_LU, ldAKban, S_ipiv, u_sol, S_ldSol, info )
    call MKLsolverResult('dgbtrs',info)
    print*, ' '
   ! print*,"Solution vector"
   !   do i =1, ntotv
   !     print '(I2, 1x, F13.5)',i, u_sol(i,1)
   !   end do
   ! print*, ' '
   ! 
   !  print*,'  •REFINING SOLUTION..... '
   !  call dgbrfs(S_trans, S_n, lowban, upban, S_nrhs, A_K, ldAKban, AK_LU, ldAKban, S_ipiv,&
   ! &           A_F, S_ldSol, u_sol, S_ldSol, S_ferr, S_berr, S_work, S_iwork, S_info )
   !  call MKLsolverResult('dgbrfs',S_info) 
   !  print*, 'Refined Solution'
   !    do i =1, ntotv
   !      print '(9F13.5)', u_sol(i,1)
   !    end do
   !  print*, ' '
   !  Write (*,*) 'Backward errors (machine-dependent)'
   !  Write (*,99) S_berr(1:S_nrhs)
   !  Write (*,*) 'Estimated forward error bounds (machine-dependent)'
   !  Write (*,99) S_ferr(1:S_nrhs)
    
   
   ! 99 Format ((3X,1P,7E11.1))
    
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
