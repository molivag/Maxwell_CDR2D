program main_CDR3d
  use library
  use param

  implicit none
  
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  double precision, allocatable, dimension(:,:) :: A_K, A_F, presc
  double precision, allocatable, dimension(:,:) :: AK_LU, Solution
  double precision, allocatable, dimension(:,:) :: N, dN_dxi, dN_deta
  double precision,            dimension(3,4) :: Hesxieta ! Puede quedar como allocatable si se allocate en ShapeFunctions subroutine.
  integer, allocatable, dimension(:,:)          :: BVs, ifpre
  integer, allocatable, dimension(:)            :: nofix
  real                                          :: start, finish

  !=============== S O L V E R ============================ 
  external                             :: mkl_dgetrfnp, dgetrf, dgetrs
  integer                              :: S_m, S_n, S_lda, S_ldb, S_infoLU, S_nrhs , S_infoSOL, i,j,k,l
  integer, allocatable, dimension(:)   :: S_ipiv
  character(len=1)                     :: S_trans
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  call cpu_time(start)
  
  call inputData ( )
  call GeneralInfo( )
  call ReadTensors(30, difma, conma, reama, force)

  print*,'size of difma', size(difma)
  print*,'size of conma', size(conma)
  print*,'size of reama', size(reama)


  !print*,'nelem ', nelem
  !print*,'nnodes', nnodes
  !print*,'nne   ', nne
  !print*,'ndofn ', ndofn
  !print*,'totGp ', totGp
  !print*,'kstab ', kstab
  !print*,'ktaum ', ktaum
  !print*,'ntotv ', ntotv
  !print*,'nevab ', nevab
  !print*,'maxban', maxband
  !! print*,'#........'
  !print*,'hnatu', hnatu
  !print*,'patau', patau
  !print*,'ElemType ',ElemType
  !print*,' '
  !
  !do i =1, nelem
  !  print*, (lnods(i,j), j=1,nne+1)
  !enddo

  !print*, ' '

  !do i =1, nnodes
  !  print*, (coord(i,j), j=1,DimPr+1)
  !enddo
  print*, '!=============== INFO DURING EXECUTION ===============!'
  call GaussQuadrature(ngaus, weigp)
  call ShapeFunctions(ngaus, nne, N, dN_dxi, dN_deta, Hesxieta)  

  print*, ' '
  print*, 'Diffusion matrix'
  do i = 1,dimPr
    do j = 1,DimPr
      print"(A,2I1)", 'k_',i,j
      do k = 1,ndofn
        print*,( difma(k,l,i,j), l=1,ndofn)
      end do
      print*,' '
    end do
  end do
 print*, 'Convection matrix'
  do k = 1, DimPr
    print"(A,2I1)",'A_',k
    do i = 1, ndofn
      write(*, *)( conma(i,j,k) ,j=1, ndofn)
    end do
    print*,''
  end do
  print*,'Reaction'
  do i=1,ndofn
    write(*, *)( reama(i,j) ,j=1,ndofn)
  end do
  print*,'force'
  do i =1, ndofn
    print*, force(i)
  end do

  call SetBoundVal( nBVs, nBVscol) !Esta funcion crea el archivo BVs.dat
  allocate( BVs(nBVs, nBVscol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadIntegerFile(40,"BVs.dat", nBVs, nBVscol, BVs)!Lectura del archivo de valores en la frontera y lo guardo en BVs
  allocate( nofix(nBVs), ifpre(ndofn,nBVs), presc(ndofn,nBVs) ) 
  allocate(A_F(ntotv, 1), Solution(ntotv, 1))

  call VinculBVs(  BVs, nofix, ifpre, presc )

  !el allocate de A_K esta dentro de GlobalSystem, primero calcula el semiancho de banda, luego aloja A_K y luego calcula Ke y rhslo
  call GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_K, A_F)
  print*, 'Shape of Global K: ',shape(A_K)
  call writeMatrix(A_K, 100, 'A_K.dat', A_F, 200, 'A_F.dat')
  call ApplyBVs(nofix,ifpre,presc,A_K,A_F)
  call writeMatrix(A_K, 111, 'A_Kbv.dat', A_F, 222, 'A_Fbv.dat')
  
  print*, ' '
  print*, 'Shape of Global F: ',shape(A_F)
  allocate(AK_LU(nband+1,ntotv))
  Solution = A_F !Solucion sera reescrito por la solucion de lapack asi no reescribo el vector global.
  AK_LU    = A_K
  DEALLOCATE(N)
  DEALLOCATE(dN_dxi)
  DEALLOCATE(dN_deta)
  DEALLOCATE(BVs)
  deallocate(nofix, ifpre, presc)
  
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

  print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
  call dgetrs( S_trans, S_n, S_nrhs, AK_LU, S_lda, S_ipiv, Solution, S_ldb, S_infoSOL )
  call MKLsolverResult( S_infoSOL )
  
  
  print*, 'Writing postprocesses files.....'
  call writeMatrix(AK_LU, 300, 'AKLU.dat', Solution, 400, 'Sol.dat')
  !call PosProcess(Solution, File_PostMsh, 'msh')
  !call PosProcess(Solution, File_PostRes, 'res')

  DEALLOCATE( S_ipiv)
  DEALLOCATE( A_F)
  DEALLOCATE( A_K)        
  DEALLOCATE( AK_LU)
  DEALLOCATE( Solution)
  call cpu_time(finish)
  print '(A11,f9.2,A8)',' CPU-Time =', finish-start, ' Seconds'
  print*, ' ' 

end program main_CDR3d
