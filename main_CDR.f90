program main_CDR3d
  use library
  use param

  implicit none
  
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  double precision, allocatable, dimension(:,:) :: A_K, A_F, presc
  double precision, allocatable, dimension(:,:) :: AK_LU, Solution
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision,            dimension(3,nne) :: Hesxieta
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish

  !integer                             :: nBVs, nBVscol      esta en el modulo como variable global
 
  !=============== S O L V E R ============================ 
  !external                             :: mkl_dgetrfnp, dgetrf, dgetrs
  !integer                              :: S_m, S_n, S_lda, S_ldb, S_infoLU, S_nrhs , S_infoSOL
  !integer, allocatable, dimension(:)   :: S_ipiv
  !character(len=1)                     :: S_trans
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  call cpu_time(start)

  call GeneralInfo( )
  call ReadIntegerFile(10,File_element, nelem, nne + 1, lnods)  
  call ReadRealFile(20,File_coord, nnodes,3, coord) 
  call ReadTensors(30, File_tensors, difma, conma, reama, force)

  print*, '!=============== INFO DURING EXECUTION ===============!'
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, nne, N, dN_dxi, dN_deta, Hesxieta)  
 
  call SetBoundVal( nBVs, nBVscol) !Esta funcion crea el archivo BVs.dat
  allocate( BVs(nBVs, nBVscol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  print*, 'numero de BV from main', nBVs 
  print*, 'numero de BV columns from main', nBVscol 
  print*, ' '
  call ReadIntegerFile(40,"BVs.dat", nBVs, nBVscol, BVs)!Lectura del archivo de valores en la frontera y lo guardo en BVs
  allocate(A_K(nband+1,ntotv))
  allocate( nofix(nBVs), ifpre(ndofn,nBVs), presc(ndofn,nBVs) ) 
  allocate(AK_LU(nband+1,ntotv))
  allocate(A_F(ntotv, 1), Solution(ntotv, 1))

  call VinculBVs(  BVs, nofix, ifpre, presc )
  call Prevop(A_K)

  ! Despues de la llamada de prevop, se guarda el valor de nband que esta global en el modulo y se puede asignar dimension a A_K
  !call GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_K, A_F)
  !call ApplyBVs(nofix,ifpre,presc,rigid,gload)
  
  print*, ' '
  print*, 'Shape of Global K: ',shape(A_K)
  print*, 'Shape of Global F: ',shape(A_F)

  Solution = A_F !Solucion sera reescrito por la solucion de lapack asi no reescribo el vector global.
  AK_LU    = A_K
  DEALLOCATE(N)
  DEALLOCATE(dN_dxi)
  DEALLOCATE(dN_deta)
  DEALLOCATE(BVs)
  
  print*,' '
  print*,'!=============== SOLVER (LAPACK) ===============!'
  !S_m   = size(AK_LU,1)
  !S_n   = size(AK_LU,2)
  !S_lda = max(1,size(AK_LU,1)) ! lda ≥ max(1, n).
  !S_trans = 'N'
  !S_nrhs  = 1
  !S_ldb   = max(1,size(Solution,1))
  !allocate( S_ipiv(max(1,min(S_m, S_n)) ) )
 
  !print*,'  •INITIALIZING LU FACTORIZATION A = P*L*U.....'
  !call dgetrf( S_m, S_n, AK_LU, S_lda, S_ipiv, S_infoLU )
  !call MKLfactoResult( S_infoLU )

  !print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
  !call dgetrs( S_trans, S_n, S_nrhs, AK_LU, S_lda, S_ipiv, Solution, S_ldb, S_infoSOL )
  !call MKLsolverResult( S_infoSOL )
  !
  !
  !print*, 'Writing postprocesses files.....'
  !call writeMatrix(A_K, 100, 'A_K.dat', A_F, 200, 'A_F.dat')
  !call writeMatrix(AK_LU, 300, 'AKLU.dat', Solution, 400, 'Sol.dat')
  !!call PosProcess(Solution, File_PostMsh, 'msh')
  !!call PosProcess(Solution, File_PostRes, 'res')

  !DEALLOCATE( S_ipiv)
  DEALLOCATE( A_F)
  DEALLOCATE( A_K)        
  !DEALLOCATE( AK_LU)
  !DEALLOCATE( Solution)
  call cpu_time(finish)
  print '(A11,f9.2,A8)',' CPU-Time =', finish-start, ' Seconds'
  print*, ' ' 

end program main_CDR3d
