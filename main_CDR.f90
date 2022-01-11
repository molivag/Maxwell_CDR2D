program main_CDR3d
  use library
  use param
  use biunit  

  implicit none
  
  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  double precision, allocatable, dimension(:,:) :: A_K, A_F, presc
  double precision, allocatable, dimension(:,:) :: AK_LU, Solution
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision,              dimension(3,4) :: Hesxieta
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
  external                             :: mkl_dgetrfnp, dgetrf, dgetrs
  integer                              :: S_m, S_n, S_lda, S_ldb, S_infoLU, S_nrhs , S_infoSOL
  integer, allocatable, dimension(:)   :: S_ipiv
  character(len=1)                     :: S_trans
  
  !---------- Input Data -----------!
  call cpu_time(start)
  call inputData ( )
  call GeneralInfo( )

  !---------- Shape Functions -----------!
  !print*, '!=============== INFO DURING EXECUTION ===============!'
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, nne, N, dN_dxi, dN_deta, Hesxieta)  
  
  !---------- Global Matrix and Vector -----------!
  !the allocate of A_K is inside of GlobalSystem, first compute the semi bandwidth, then allocate A_K
  call GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_K, A_F)
  call writeMatrix(A_K, 10, 'A_K.dat', A_F, 20, 'A_F.dat')

  !---------- Boundary Conditions -----------!
  call SetBoundVal( nBVs, nBVscol) !Esta funcion crea el archivo BVs.dat
  allocate( BVs(nBVs, nBVscol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadIntegerFile(30,"BVs.dat", nBVs, nBVscol, BVs)!Lectura del archivo de boundary values will stored at BVs
  allocate( nofix(nBVs), ifpre(ndofn,nBVs), presc(ndofn,nBVs) ) 
  call VinculBVs(  BVs, nofix, ifpre, presc )
  call ApplyBVs(nofix,ifpre,presc,A_K,A_F)
  call writeMatrix(A_K, 40, 'A_Kbv.dat', A_F, 50, 'A_Fbv.dat')
  
  !---------- Memory Relase -----------!
  DEALLOCATE(N)
  DEALLOCATE(dN_dxi)
  DEALLOCATE(dN_deta)
  DEALLOCATE(BVs)
  DEALLOCATE(nofix, ifpre, presc)
  
  !---------- Solving System of Equations -----------!
  print*,'!================ SOLVER (LAPACK) ================!'
  allocate(AK_LU(nband+1,ntotv))
  allocate(Solution(ntotv, 1))
  S_m   = size(AK_LU,1)
  S_n   = size(AK_LU,2)
  S_lda = max(1,size(AK_LU,2)) ! lda ≥ max(1, n).
  S_trans  = 'N'
  S_nrhs   = 1
  S_ldb    = max(1,size(Solution,1))
  Solution = A_F !Solucion sera reescrito por la solucion de lapack asi no reescribo el vector global.
  AK_LU    = A_K
  allocate( S_ipiv(max(1,min(S_m, S_n)) ) )
  print*,'  •INITIALIZING LU FACTORIZATION A = P*L*U.....'
  call dgetrf( S_m, S_n, AK_LU, S_lda, S_ipiv, S_infoLU )
  call MKLfactoResult( S_infoLU )
  print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
  call dgetrs( S_trans, S_n, S_nrhs, AK_LU, S_lda, S_ipiv, Solution, S_ldb, S_infoSOL )
  call MKLsolverResult( S_infoSOL )
  
  !---------- Print and write results -----------!
  !print*, 'Writing postprocesses files.....'
  call writeMatrix(AK_LU, 60, 'AKLU.dat', Solution, 70, 'Sol.dat')
  print*, 'Shape of Global K: ',shape(A_K)
  print*, 'Shape of Global F: ',shape(A_F)
  !call PosProcess(Solution, File_PostMsh, 'msh')
  !call PosProcess(Solution, File_PostRes, 'res')

  !---------- Memory Relase -----------!
  DEALLOCATE( S_ipiv)
  DEALLOCATE( A_F)
  DEALLOCATE( A_K)        
  DEALLOCATE( AK_LU)
  DEALLOCATE( Solution)
  call cpu_time(finish)
  write(*,'(A11,f9.2,A8)'),' CPU-Time =', finish-start, ' Seconds', ' '

end program main_CDR3d
