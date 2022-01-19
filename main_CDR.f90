program main_CDR3d
  use library
  use param
  use biunit 
  use solver
  !use boundVal

  implicit none
  
  ! - - - - - - - - - - * * * Variable declaration * * * * * * * - - - - - - - - - -!
  double precision, allocatable, dimension(:,:) :: A_K, A_Kbv, A_F, Solution, presc
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
  
    external :: dgbtrf, dgbtrs, dpbtrs, dpbtrf
    integer, allocatable :: ipiv(:)
    integer :: i, j, k,l, ii,jj, kl, ku, ldab, S_m, S_n, nrhs, info, ldb
    double precision, allocatable :: ab(:, :), ablap(:, :), A_Ktot(:,:), AKb_LU(:,:), b(:,:)
    character(len=1) :: trans, uplo
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


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
  allocate(A_Kbv(nband+1,ntotv))
  A_Kbv = A_K
  call VinculBVs(  BVs, nofix, ifpre, presc )
  call ApplyBVs(nofix,ifpre,presc,A_Kbv,A_F )
  call writeMatrix(A_Kbv, 40, 'A_Kbv.dat', A_F, 50, 'A_Fbv.dat')
  !write(*,"(4F10.5)") ((A_K(i,j), j=1,S_n), i=1,nband+1)  
 print*, 'nabnd',nband 
  !---------- Memory Relase -----------!
  DEALLOCATE(N)
  DEALLOCATE(dN_dxi)
  DEALLOCATE(dN_deta)
  DEALLOCATE(BVs)
  DEALLOCATE(nofix, ifpre, presc)
  
  !---------- Solving System of Equations -----------!

  allocate(Solution(ntotv, 1))
  Solution = 0.0
  call solsistem(A_Kbv,A_Kbv,1,A_F,Solution)
  
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  write(*,*)
  write(*,*)
  S_m=9
  S_n=9
  kl =0
  ku =4
  ldab = 2*kl + ku + 1
  trans= 'N'
  uplo = 'U'
  nrhs = 1
  ldb  = max(1, S_n)
  
  allocate(A_Ktot(S_m,S_n), ab(ldab,S_n), ablap(ldab,S_n), ipiv(S_n), AKb_LU(ldab,S_n), b(ldb,1))
  !    Read A from data file
  ab     = 0.0
  ablap  = 0.0
  AKb_LU = 0.0
  ipiv   = 0 
  
  open(unit=1010, file= 'globalK.dsc', STATUS="old", ACTION="read") 
  read(1010,*) ((A_Ktot(ii,jj), jj=1,S_n), ii=1,S_m) 
  
  close(1010)
  print*, ''
  print("(A,2(I2,2x))"), 'The general banded square matrix A_Ktot of shape:', shape(A_Ktot)
  write(*,"(9F10.5)") ((A_Ktot(i,j), j=1,S_n), i=1,S_m)  

  print*, ' '
print("(A,2(I2,2x))"),'Transforming the general banded matrix to general-band storage mode of shape: ', shape(ab)
  k   = kl + ku + 1
  do j=1,S_n
     do  i=max(1, j-ku),min(S_n,j+kl)!i=max(j-ku,1),min(j+kl,S_n)
       ablap(i-j+k,j)=a_ktot(i,j)
     end do
   end do
  write(*,"(9F10.5)") ((ablap(i,j), j=1,S_n), i=1,ldab)  

  print*, ' '
  print("(A,2(I2,2x))"), 'The banded matrix from retpla.f whit shape:', shape(A_K)
  write(*,"(9F10.5)") ((A_K(i,j), j=1,ntotv), i=1,nband+1)  
  !read(10,*)((ab(k+i-j,j),j=max(i-kl,1),min(i+ku,S_n)), i=1, S_m)

  print*, ' '
  print*, 'Reordering the band matrix from retpla to LAPACK format'
  i=0
  j=0
  do k = nband+1, 1,-1
    i = i+1
    j = 1
    do l = ntotv,1,-1
      ab(i,j) = A_K(k,l)
      j=j+1
    end do
  end do
  write(*,"(9F10.5)") ((ab(i,j), j=1,S_n), i=1,ldab)  



  AKb_LU = ab

  write(*,*)
  !call dpbtrf( uplo, S_n, nband, AKb_LU, ldab, info )
  call dgbtrf(S_m, S_n, kl, ku, AKb_LU, ldab, ipiv, info)

  !print*, 'info Cholesky decomposition' , info
  !print*, 'Cholesky-decomposition of band-storaged Matrix'
  print*, 'info LU decomposition' , info
  print*, 'LU-decomposition of band-storaged Matrix'

  write(*,"(9F10.5)") ((AKb_LU(i,j), j=1,S_n), i=1,ldab)  
  Write (*, *)
  Write (*, *) 'IPIV: Permutation'
  Write (*, 100) ipiv(1:min(S_m,S_n))
  If (info/=0) Then
    Write (*, *) 'The factor U is singular'
  End If
  
  print*, 'Iniciando solver' 
  b = A_F
  !call dpbtrs( uplo, S_n, nband, 1, AKb_LU, ldab, b, ldb, info )
  call dgbtrs( trans, S_n, kl, ku, nrhs, AKb_LU, ldab, ipiv, b, ldb, info )
  print*, 'info solver' , info
  print*, 'Solution'
  write(*,"(1F15.5)") (b(i,1), i=1,S_n)  

  100 Format((3X,9I11))










  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !---------- Print and write results -----------!
  !print*, 'Writing postprocesses files.....'
  !call writeMatrix(A_K, 60, '-.dat', Solution, 70, 'SolCholesky.dat')
  !print*, ' ' 
  print*, 'Shape of Global K: ',shape(A_K)
  print*, 'Shape of Global F: ',shape(A_F)
  print*,' Shape of Solution: ',shape(Solution)
  !call PosProcess(Solution, File_PostMsh, 'msh')
  !call PosProcess(Solution, File_PostRes, 'res')

  !---------- Memory Relase -----------!
  DEALLOCATE( A_F)
  DEALLOCATE( A_K)        
  DEALLOCATE( A_Kbv)        
  DEALLOCATE( Solution)
  call cpu_time(finish)
  write(*,'(A11,f9.2,A8)'),' CPU-Time =', finish-start, ' Seconds', ' '

end program main_CDR3d
