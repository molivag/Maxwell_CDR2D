program main_CDR3d
  use library
  use param
  use biunit 
  use boundVal
  !use solver
!
  implicit none
  
  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  double precision, allocatable, dimension(:,:) :: A_K, A_F,  presc
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision,              dimension(3,4) :: Hesxieta
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish
  ! - - - - - - - - * * * Variable declaration (SOLVER) * * * * * * * - - - - - - - !
  external :: dgbtrf, dgbtrs, dpbtrs, dpbtrf
  integer, allocatable, dimension(:) :: S_ipiv
  integer :: S_m, S_n, S_nrhs, info, S_ldSol
  double precision, allocatable, dimension(:,:) :: AKbLU, Sols
  character(len=1) :: S_trans, S_uplo

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

  !---------- Boundary Conditions -----------!
  call SetBoundVal( nBVs, nBVscol) !Esta funcion crea el archivo BVs.dat
  allocate( BVs(nBVs, nBVscol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadIntegerFile(30,"BVs.dat", nBVs, nBVscol, BVs)!Reading of boundary values file will store at BVs
  allocate( nofix(nBVs), ifpre(ndofn,nBVs), presc(ndofn,nBVs) ) 
  call VinculBVs(  BVs, nofix, ifpre, presc )
  call ApplyBVs(nofix,ifpre,presc,A_K,A_F )

  !---------- Memory Relase -----------!
  DEALLOCATE(N)
  DEALLOCATE(dN_dxi)
  DEALLOCATE(dN_deta)
  DEALLOCATE(BVs)
  DEALLOCATE(nofix, ifpre, presc)



  write(*,*) ''
  print*,'!================= MKL <S>OLVER ===============!'
  S_m     = ntotv       !size(AK,1)
  S_n     = ntotv       !size(AK,2)
  S_ldSol = max(1,S_n)
  S_trans = 'N'
  S_uplo  = 'U'
  S_nrhs  = 1


  allocate( AKbLU(ldAKban,ntotv))
  allocate( Sols(S_ldSol, 1))
  allocate( S_ipiv(max(1,min(S_m, S_n)) ) ) !size (min(m,n))
  
  AKbLU = A_K         !AK_band(ldab,*) The array AK_band contains the matrix A_K in band storage
  Sols  = A_F         !Sol_vec will be rewrited by LAPACK solution avoiding lose A_F

  !---------- Solving System of Equations by retpla solver -----------!
  !  print*, 'Cholesky-decomposition of band-storaged Matrix'
  !  call solsistem(A_K,A_K,1,A_F,Sols)
  
  print*,'  •INITIALIZING BAND LU DECOMPOSITION.....'
  !print*, 'Cholesky-decomposition of band-storaged Matrix'
  !call dpbtrf( uplo, S_n, nband, AKb_LU, ldAKban, info )
  call dgbtrf(S_m, S_n, lowban, upban, AKbLU, ldAKban, S_ipiv, info)
  call MKLfactoResult('dgbtrf',info)
  
  print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
  !call dpbtrs( uplo, S_n, nband, 1, AKb_LU, ldAKban, b, ldb, info )
  call dgbtrs( S_trans, S_n, lowban, upban, S_nrhs, AKbLU, ldAKban, S_ipiv, Sols, S_ldSol, info )
  call MKLsolverResult('dgbtrs',info)
  

  !---------- Print and write results -----------!
  print*,'!=============== Output Files ================!'
  !call writeMatrix(A_K, 10, 'A_K.dat', A_F, 20, 'A_F.dat')
  call writeMatrix(AKbLU,60,'-', Sols, 70, 'SolMKL_LU.dat')
  print*, ' ' 
  print*, 'Shape of Global K: ',shape(A_K)
  print*, 'Shape of Global F: ',shape(A_F)
  print*, 'Shape of Solution: ',shape(Sols)
  write(*,*)
  call PosProcess(Sols, File_PostMsh, 'msh')
  call PosProcess(Sols, File_PostRes, 'res')
  write(*,*)

  !---------- Memory Relase -----------!
  DEALLOCATE( A_F)
  DEALLOCATE( A_K)
  DEALLOCATE( AKbLU)
  DEALLOCATE( Sols)
  call cpu_time(finish)
  write(*,'(A11,f9.2,A8)')' CPU-Time =', finish-start, ' Seconds', ' '

end program main_CDR3d
