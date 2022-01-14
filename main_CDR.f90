program main_CDR3d
  use library
  use param
  use biunit 
  use solver
  !use boundVal

  implicit none
  
  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  double precision, allocatable, dimension(:,:) :: A_K, A_F, Solution, presc
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision,              dimension(3,4) :: Hesxieta
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
 ! external                             :: 
 ! integer                              ::
 ! integer, allocatable, dimension(:)   :: 
 ! character(len=1)                     :: 
  
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
  print*, 'nbvs', nbvs
  print*, 'nbvscol', nbvscol
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

  allocate(Solution(ntotv, 1))
  Solution = 0.0
  call solsistem(A_K,A_K,1,A_F,Solution)
  
  
  !---------- Print and write results -----------!
  !print*, 'Writing postprocesses files.....'
  call writeMatrix(A_K, 60, '-.dat', Solution, 70, 'SolCholesky.dat')
  print*, ' ' 
  print*, 'Shape of Global K: ',shape(A_K)
  print*, 'Shape of Global F: ',shape(A_F)
  print*,' Shape of Solution: ',shape(Solution)
  !call PosProcess(Solution, File_PostMsh, 'msh')
  !call PosProcess(Solution, File_PostRes, 'res')

  !---------- Memory Relase -----------!
  DEALLOCATE( A_F)
  DEALLOCATE( A_K)        
  DEALLOCATE( Solution)
  call cpu_time(finish)
  write(*,'(A11,f9.2,A8)'),' CPU-Time =', finish-start, ' Seconds', ' '

end program main_CDR3d
