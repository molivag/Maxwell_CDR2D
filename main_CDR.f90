program main_CDR3d
  use library
  use param
  use param

  implicit none
  
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  double precision, allocatable, dimension(:,:) :: A_K, rhsgl, AK_LU, Solution, N, dN_dxi, dN_deta
  integer, allocatable, dimension(:,:) :: BVs
  integer                              :: nBVs, nBVscol
  real                                 :: start, finish
 
  !=============== S O L V E R ============================ 
  external                             :: mkl_dgetrfnp, dgetrf, dgetrs
  integer                              :: S_m, S_n, S_lda, S_ldb, S_infoLU, S_nrhs , S_infoSOL
  integer, allocatable, dimension(:)   :: S_ipiv
  character(len=1)                     :: S_trans
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  call cpu_time(start)

  call GeneralInfo( )
  call ReadIntegerFile(10,File_element, nelem, nne + 1, lnods)  
  call ReadRealFile(20,File_coord, nnodes,3, coord) !Para dreducir el numero de subrutinas, usar la sentencia option par
  call ReadTensors(30, File_tensors, difma, conma, reama, force) !Para dreducir el numero de subrutinas, usar la sentencia option par

  print*, ' '
  print*, '!=============== INFO DURING EXECUTION ===============!'
  
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, nne, N, dN_dxi, dN_deta)  
  
  allocate(A_K(2*nnodes, 2*nnodes), AK_LU(2*nnodes, 2*nnodes) )
  call SetBounCond( nBVs, nBVscol) !Esta funcion crea el archivo bcsVP.dat
  allocate( BVs(nBVs, nBVscol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadIntegerFile(60,"BVs.dat", nBVs, nBVscol, BVs)!Llamo el archivo de valores en la frontera y lo guardo en BVs
  
  call GlobalK( A_K, dN_dxi, dN_deta)

  allocate(rhsgl(2*nnodes, 1), Solution(2*nnodes, 1))
  rhsgl = 0.0 !initializing source vector (rhsgl) 
  call ApplyBoundCond(nBVs, BVs, A_K, rhsgl )
  
  Solution = rhsgl !Solucion sera reescrito por la solucion de lapack asi no reescribo el vector global.
  AK_LU    = A_K
  DEALLOCATE(N)
  DEALLOCATE(dN_dxi)
  DEALLOCATE(dN_deta)
  DEALLOCATE(BVs)
  DEALLOCATE(rhsgl )
  
  print*,' '
  print*,'!=============== SOLVER (LAPACK) ===============!'
  S_m   = size(AK_LU,1)
  S_n   = size(AK_LU,2)
  S_lda = max(1,size(AK_LU,1)) ! lda ≥ max(1, n).
  S_trans = 'N'
  S_nrhs  = 1
  S_ldb   = max(1,size(Solution,1))
  allocate( S_ipiv(max(1,min(S_m, S_n)) ) )
 
  print*,'  •INITIALIZING LU FACTORIZATION A = P*L*U.....'
  call dgetrf( S_m, S_n, AK_LU, S_lda, S_ipiv, S_infoLU )
  call MKLfactoResult( S_infoLU )

  print*, ' '  
  print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
  call dgetrs( S_trans, S_n, S_nrhs, AK_LU, S_lda, S_ipiv, Solution, S_ldb, S_infoSOL )
  call MKLsolverResult( S_infoSOL )
  
  print*, 'Writing postprocesses files.....'
  DEALLOCATE( S_ipiv)
  
  !call writeMatrix(A_K, 111, 'GlobalK.dat', Solution, 444, 'Sol.dat')
  call PosProcess(Solution, File_PostMsh, 'msh')
  call PosProcess(Solution, File_PostRes, 'res')

  DEALLOCATE( A_K)        
  DEALLOCATE( AK_LU)
  DEALLOCATE( Solution)
  print*,' '  
  print"(A6,A19, A38)", ' File ',File_PostMsh,' written succesfully in Pos/ . . . . .'
  print"(A6,A19, A38)", ' File ',File_PostRes, 'written succesfully in Pos/ . . . . .'
  print*, ' ' 
  call cpu_time(finish)
  print '(A11,f9.2,A8)',' CPU-Time =', finish-start, ' Seconds'
  print"(A)", ' ' 

end program main_CDR3d
