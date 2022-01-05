module param
  
  implicit none

  character(len=14) :: ElemType
  integer           :: nevab, ntotv
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
  !character(len=17), parameter :: File_PostMsh  = 'CDR3d.post.msh'
  !character(len=17), parameter :: File_PostRes  = 'CDR3d.post.res'
  contains
    
    subroutine inputData( ) 
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !
      ! subrutina que lee todos los parametros de entrada para la simulacion, 
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion
      !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      implicit none
      
      integer :: i,j, stat
      CHARACTER(len=80) :: msg
      character(len=*), parameter  :: fileplace = "~/Dropbox/1.Doctorado/1.Research/1.Computing/Fortran/2.ConDifRea/Geo/"
      
      open(5, file=fileplace//'parameters.dat',status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      read(5, 102)  ElemType, DimPr, nelem, nnodes, nne, ndofn, totGp, maxband, kstab, ktaum, patau, hnatu
      
      allocate(lnods(nelem,nne+1))
      allocate(coord(nnodes,Dimpr+1))
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
      
      102 format( 7/, 11x, A14, /, 7(11x,I5,/), 2/, 2(11x,I5,/),11x,f3.4,/,11x,F3.4,2/)
      
    end subroutine inputData
    
    subroutine ReadTensors(nr, difma, conma, reama, force)
      
      !character(len=*), parameter   :: fileplace = "~/Dropbox/1.Doctorado/1.Research/1.Computing/Fortran/2.ConDifRea/Geo/"
      !character(len=*), intent (in) :: FileName
      integer :: nr
      double precision, allocatable, dimension(:,:,:,:) :: difma
      double precision, allocatable, dimension(:,:,:)   :: conma
      double precision, allocatable, dimension(:,:)     :: reama 
      double precision, allocatable, dimension(:)       :: force 
      
      !open (unit = nr, file = fileplace//FileName, status='old', iostat = status)
      allocate( difma(ndofn,ndofn,DimPr,DimPr) )
      allocate( conma(ndofn,ndofn,DimPr) )
      allocate( reama(ndofn,ndofn) )
      allocate( force(ndofn) )
      
      difma = 0.0
      conma = 0.0
      reama = 0.0
      force = 0.0
      
      2 format(39x,2(e15.5),/,39x,2(e15.5))
      !5 format(39x,3(E15.5),/,39x,3(E15.5),/,39x,3(E15.5))
      !read(nr,1) npoin,nelem,nnode,ngaut,ndofn
      
      if(ndofn.eq.1) then
        difma(1,1,1,1)= 1.0e-4
        difma(1,1,1,2)= 0.0e-0
        difma(1,1,2,2)= 1.0e-4
        conma(1,1,1)  = 1.0e+1
        conma(1,1,2)  = 0.0e+0
        reama(1,1)    = 0.0
        force(1)      = 1.0e+0
      elseif(ndofn.eq.2)then
        difma(1,1,1,1) = 1.0e-0 
        difma(1,2,1,1) = 0.0e-0
        difma(2,1,1,1) = 0.0e-0
        difma(2,2,1,1) = 1.0e-0
        difma(1,1,1,2) = 0.0e-0
        difma(1,2,1,2) = 0.0e-0
        difma(2,1,1,2) = 0.0e-0
        difma(2,2,1,2) = 0.0e-0
        difma(1,1,2,2) = 1.0e-0
        difma(1,2,2,2) = 0.0e-0
        difma(2,1,2,2) = 0.0e-0
        difma(2,2,2,2) = 1.0e-0
        !read(nr,2) difma(1,1,1,1),difma(1,2,1,1),difma(2,1,1,1),difma(2,2,1,1)
        !read(nr,2) difma(1,1,1,2),difma(1,2,1,2),difma(2,1,1,2),difma(2,2,1,2)
        !read(nr,2) difma(1,1,2,2),difma(1,2,2,2),difma(2,1,2,2),difma(2,2,2,2)
        
        conma(1,1,1) = 0.0e+2 
        conma(1,2,1) = 0.0e-0
        conma(2,1,1) = 0.0e-0 
        conma(2,2,1) = 0.0e+2
        conma(1,1,2) = 0.0e+0
        conma(1,2,2) = 0.0e-0
        conma(2,1,2) = 0.0e-0
        conma(2,2,2) = 0.0e+1
        !read(nr,2) conma(1,1,1), conma(1,2,1),conma(2,1,1), conma(2,2,1)
        !read(nr,2) conma(1,1,2), conma(1,2,2),conma(2,1,2), conma(2,2,2)
       
        reama(1,1) = 8.0e+4
        reama(1,2) = 0.0e-0
        reama(2,1) = 0.0e-0
        reama(2,2) = 8.0e+4
        !read(nr,2) reama(1,1), reama(1,2), reama(2,1), reama(2,2)
        
        force(1) = 1.0e+0
        force(2) = 1.0e+0
        !read(nr,3) force(1), force(2)
        
        !print*, force 
        
      else if(ndofn.eq.3) then                              
        !print*, 'test if'
       ! read(nr,5) difma(1,1,1,1),difma(1,2,1,1),difma(1,3,1,1)
       ! read(nr,5) difma(2,1,1,1),difma(2,2,1,1),difma(2,3,1,1)
       ! read(nr,5) difma(3,1,1,1),difma(3,2,1,1),difma(3,3,1,1)
       ! 
       ! read(nr,5) difma(1,1,1,2),difma(1,2,1,2),difma(1,3,1,2)
       ! read(nr,5) difma(2,1,1,2),difma(2,2,1,2),difma(2,3,1,2) 
       ! read(nr,5) difma(3,1,1,2),difma(3,2,1,2),difma(3,3,1,2)
       ! 
       ! read(nr,5) difma(1,1,2,2),difma(1,2,2,2),difma(1,3,2,2)
       ! read(nr,5) difma(2,1,2,2),difma(2,2,2,2),difma(2,3,2,2)
       ! read(nr,5) difma(3,1,2,2),difma(3,2,2,2),difma(3,3,2,2)
       ! 
       ! read(nr,5) conma(1,1,1), conma(1,2,1), conma(1,3,1)
       ! read(nr,5) conma(2,1,1), conma(2,2,1), conma(2,3,1)
       ! read(nr,5) conma(3,1,1), conma(3,2,1), conma(3,3,1)
       ! 
       ! read(nr,5) conma(1,1,2), conma(1,2,2), conma(1,3,2)
       ! read(nr,5) conma(2,1,2), conma(2,2,2), conma(2,3,2)
       ! read(nr,5) conma(3,1,2), conma(3,2,2), conma(3,3,2)
       !
       ! read(nr,5) reama(1,1), reama(1,2), reama(1,3)
       ! read(nr,5) reama(2,1), reama(2,2), reama(2,3)
       ! read(nr,5) reama(3,1), reama(3,2), reama(3,3)
       ! 
       ! read(nr,3) force(1), force(2), force(3)
        
        
      end if
      
      close (nr)
      
      !The slash / descriptor 
      !begins a new line (record) on output and 
      !skips to the next line on input, 
      !and ignoring any unread information on the current record format(6/) o 6/
      3 format(39x,3(e15.5))
      5 format(39x,3(E15.5),/,39x,3(E15.5),/,39x,3(E15.5))
      
      if(ndofn.eq.2)then 
        difma(1,1,2,1) = difma(1,1,1,2)
        difma(1,2,2,1) = difma(2,1,1,2)
        difma(2,1,2,1) = difma(1,2,1,2)
        difma(2,2,2,1) = difma(2,2,1,2)
      endif
      if(ndofn.eq.3) then
        difma(1,3,2,1)=difma(3,1,1,2)
        difma(2,3,2,1)=difma(3,2,1,2)
        difma(3,1,2,1)=difma(1,3,1,2)
        difma(3,2,2,1)=difma(2,3,1,2)
        difma(3,3,2,1)=difma(3,3,1,2)
      end if
      
    end subroutine ReadTensors
    
end module param
