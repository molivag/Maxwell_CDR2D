module param
  
  implicit none

  character(len=14) :: ElemType
  integer           :: nevab, ntotv
  integer           :: upban, lowban, totban, ldAKban !variables defined in GlobalSystem
  integer           :: DimPr, nelem, nnodes, nne, ndofn, totGp, kstab, ktaum, maxband, nBVs, nBVscol, nband
  real              :: hnatu, patau
  integer,          allocatable, dimension(:,:)     :: lnods
  real,             allocatable, dimension(:,:)     :: coord
  double precision, allocatable, dimension(:,:)     :: ngaus, weigp

  double precision, allocatable, dimension(:,:,:,:) :: difma
  double precision, allocatable, dimension(:,:,:)   :: conma
  double precision, allocatable, dimension(:,:)     :: reama !Tensor materials
  double precision, allocatable, dimension(:)       :: force !Force vector 
  
  !character(len=20), parameter :: File_element  = 'lnods.dat'
  !character(len=20), parameter :: File_coord    = 'coord.dat'
  !character(len=20), parameter :: File_tensors  = 'tensors.dat'
  character(len=15), parameter :: File_PostMsh  = 'CDR3d.post.msh'
  character(len=15), parameter :: File_PostRes  = 'CDR3d.post.res'
  contains
    
    subroutine inputData( ) 
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !
      ! subrutina que lee todos los parametros de entrada para la simulacion, 
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion
      !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      implicit none
      
      integer :: i,j, stat! k, l, 
      character(len=80) :: msg
      character(len=*), parameter  :: fileplace = "./"
      ! double precision, allocatable, dimension(:,:,:,:) :: difma
      ! double precision, allocatable, dimension(:,:,:)   :: conma
      ! double precision, allocatable, dimension(:,:)     :: reama 
      ! double precision, allocatable, dimension(:)       :: force 
      open(5, file=fileplace//'inputCDR.dsc',status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      ! read(5, 100)  ElemType, DimPr, nelem, nnodes, nne, ndofn, totGp, maxband, kstab, ktaum, patau, hnatu
      read(5, 100)  ElemType, DimPr, nelem, nnodes, nne, & 
      ndofn, totGp, maxband, kstab, ktaum, patau, hnatu

      allocate(lnods(nelem,nne+1))
      allocate(coord(nnodes,Dimpr+1))
      allocate( difma(ndofn,ndofn,DimPr,DimPr) )
      allocate( conma(ndofn,ndofn,DimPr) )
      allocate( reama(ndofn,ndofn) )
      allocate( force(ndofn) )
      
      difma = 0.0
      conma = 0.0
      reama = 0.0
      force = 0.0
      
      if(ndofn.eq.1)then
        read(5,101) difma(1,1,1,1), difma(1,1,1,2), difma(1,1,2,2)
        ! read(5,101) difma(1,1,1,2)
        ! read(5,101) difma(1,1,2,2)
    
        read(5,101) conma(1,1,1)
        read(5,101) conma(1,1,2)
    
        read(5,101) reama(1,1)
        
        read(5,105) force(1)
        
      elseif(ndofn.eq.2) then
        read(5,102) difma(1,1,1,1), difma(1,2,1,1), difma(2,1,1,1),difma(2,2,1,1)
        read(5,102) difma(1,1,1,2), difma(1,2,1,2), difma(2,1,1,2),difma(2,2,1,2)
        read(5,102) difma(1,1,2,2), difma(1,2,2,2), difma(2,1,2,2),difma(2,2,2,2)
    
        read(5,102) conma(1,1,1), conma(1,2,1), conma(2,1,1), conma(2,2,1)
        read(5,102) conma(1,1,2), conma(1,2,2), conma(2,1,2), conma(2,2,2)
    
        read(5,102) reama(1,1), reama(1,2), reama(2,1), reama(2,2)
        
        read(5,106) force(1), force(2)
        
      elseif(ndofn.eq.3)then
    
        read(5,103) &
        difma(1,1,1,1), difma(1,2,1,1), difma(1,3,1,1), &
        difma(2,1,1,1), difma(2,2,1,1), difma(2,3,1,1), &
        difma(3,1,1,1), difma(3,2,1,1), difma(3,3,1,1)
        read(5,103) &
        difma(1,1,1,2), difma(1,2,1,2), difma(1,3,1,2), &
        difma(2,1,1,2), difma(2,2,1,2), difma(2,3,1,2), &
        difma(3,1,1,2), difma(3,2,1,2), difma(3,3,1,2)
        read(5,103) &
        difma(1,1,2,2), difma(1,2,2,2), difma(1,3,2,2), &
        difma(2,1,2,2), difma(2,2,2,2), difma(2,3,2,2), &
        difma(3,1,2,2), difma(3,2,2,2), difma(3,3,2,2)
           
        read(5,103) &
        conma(1,1,1), conma(1,2,1), conma(1,3,1), &
        conma(2,1,1), conma(2,2,1), conma(2,3,1), &
        conma(3,1,1), conma(3,2,1), conma(3,3,1)
        read(5,103) &
        conma(1,1,2), conma(1,2,2), conma(1,3,2), &
        conma(2,1,2), conma(2,2,2), conma(2,3,2), &
        conma(3,1,2), conma(3,2,2), conma(3,3,2)
        
        read(5,103) &
        reama(1,1), reama(1,2), reama(1,3), &
        reama(2,1), reama(2,2), reama(2,3), &
        reama(3,1), reama(3,2), reama(3,3)
        
       read(5,107) force(1), force(2), force(3)
      
      end if

      do i=1,nelem
        read(5,*,iostat=stat,iomsg=msg) (lnods(i,j), j =1,nne+1)
        IF ( stat /= 0 )then
          print*,stat
          print*, msg
        end if
      end do
      do i=1,nnodes
        read(5,*,iostat=stat,iomsg=msg) (coord(i,j), j =1,DimPr+1 )
        IF ( stat /= 0 )then
          print*,stat
          print*, msg
        end if
      end do
      
      close(5)
     
      nevab = ndofn*nne   
      ntotv = ndofn*nnodes

      if(ndofn.eq.1)then
        difma(1,1,2,1) = difma(1,1,1,2)
      elseif(ndofn.eq.2)then
        difma(1,1,2,1) = difma(1,1,1,2)
        difma(1,2,2,1) = difma(2,1,1,2)
        difma(2,1,2,1) = difma(1,2,1,2)
        difma(2,2,2,1) = difma(2,2,1,2)
      elseif(ndofn.eq.3) then
        difma(1,3,2,1) = difma(3,1,1,2)
        difma(2,3,2,1) = difma(3,2,1,2)
        difma(3,1,2,1) = difma(1,3,1,2)
        difma(3,2,2,1) = difma(2,3,1,2)
        difma(3,3,2,1) = difma(3,3,1,2)
      end if
      
      100 format( 7/, 11x, A14, /, 7(11x,I5,/), 2/, 2(11x,I5,/),11x,f3.4,/,11x,F3.42,2/)    
      101 format(1/,e15.5,2/)!,3/,e9.2,3/,e9.2,3/,e9.2)!,/,e9.2,/,e9.2)
      102 format(1/,e15.5, e15.5,/, e15.5,e15.5,/)
      103 format(1/,3(e15.5))
      105 format(1/,e15.5,2/) !e15.5, e15.5,/)
      106 format(1/,e15.5,e15.5,2/) !e15.5, e15.5,/)
      107 format(1/,e15.5,e15.5,e15.5,2/) !e15.5, e15.5,/)
      ! print*, ' '
      ! print*, 'Diffusion matrix'
      ! do i = 1,dimPr
      !   do j = 1,DimPr
      !     print"(A,2I1)", 'k_',i,j
      !     do k = 1,ndofn
      !       print"(e10.3,1x,e10.3, 1x, e10.3)",( difma(k,l,i,j), l=1,ndofn)
      !     end do
      !     print*,' '
      !   end do
      ! end do
      ! print*, ' '  
      ! print*, 'Convection matrix'
      ! do k = 1, DimPr
      !   print"(A,2I1)",'A_',k
      !   do i = 1, ndofn
      !     write(*, *)( conma(i,j,k) ,j=1, ndofn)
      !   end do
      !   print*,''
      ! end do
      ! print*,'Reaction'
      ! do i=1,ndofn
      !   write(*, *)( reama(i,j) ,j=1,ndofn)
      ! end do
      ! print*,'force'
      ! do i =1, ndofn
      !   print*, force(i)
      ! end do
      ! print*, ' '
    end subroutine inputData
    
    
end module param
