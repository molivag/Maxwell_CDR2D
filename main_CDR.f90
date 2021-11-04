program main_CDR3d
  use library
  use param
  use param

  implicit none
  
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  double precision, allocatable, dimension(:,:) :: A_K, Sv, AK_LU, Solution, N, dN_dxi, dN_deta
  integer, allocatable, dimension(:,:) :: Fbcsvp
  integer                              :: NoBV, NoBVcol
  real                                 :: start, finish
 
  !=============== S O L V E R ============================ 
  external                             :: mkl_dgetrfnp, dgetrf, dgetrs
  integer                              :: S_m, S_n, S_lda, S_ldb, S_infoLU, S_nrhs , S_infoSOL
  integer, allocatable, dimension(:)   :: S_ipiv
  character(len=1)                     :: S_trans
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  call cpu_time(start)

  call GeneralInfo( )
  call ReadIntegerFile(10,File_element, Nelem, nUne + 1, elements)  
  call ReadRealFile(20,File_nodes, n_nodes,3, nodes) !Para dreducir el numero de subrutinas, usar la sentencia option par
  !call ReadReal(30,File_material, materials)    !Para dreducir el numero de subrutinas, usar la sentencia option para      
  !call ReadIntegerFile(40,File_pnodes, n_nodes,2, pnodes)
  !call ReadIntegerFile(50,File_pelement, Nelem,nPne + 1, pelements)
  print*, ' '
  print*, '!=============== INFO DURING EXECUTION ===============!'
  
  call GaussQuadrature(gauss_points, gauss_weights)
  call ShapeFunctions(gauss_points, nUne, N, dN_dxi, dN_deta)  
  
  allocate(A_K(2*n_nodes, 2*n_nodes), AK_LU(2*n_nodes, 2*n_nodes) )
  call SetBounCond( NoBV, NoBVcol) !Esta funcion crea el archivo bcsVP.dat
  allocate( Fbcsvp(NoBV, NoBVcol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadIntegerFile(60,"Fbcsvp.dat", NoBV, NoBVcol, Fbcsvp)!Llamo el archivo de valores en la frontera y lo guardo en Fbcsvp
  
  call GlobalK( A_K, dN_dxi, dN_deta)

  allocate(Sv(2*n_nodes, 1), Solution(2*n_nodes, 1))
  Sv = 0.0 !initializing source vector (Sv) 
  call ApplyBoundCond(NoBV, Fbcsvp, A_K, Sv )
  
  Solution = Sv !Solucion sera reescrito por la solucion de lapack asi no reescribo el vector global.
  AK_LU    = A_K
  DEALLOCATE(N)
  DEALLOCATE(dN_dxi)
  DEALLOCATE(dN_deta)
  DEALLOCATE(Fbcsvp)
  DEALLOCATE(Sv )
  
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
