module geometry
use param
  
  implicit none

  double precision, allocatable, dimension(:,:)     :: coord !, coordRef
  integer,          allocatable, dimension(:,:)     :: lnods !, lnodsRef
  character(len=4) :: ElemType
  integer           :: nelem, nnodes, nevab, ntotv

  !common/contr/nin,nou,DimPr,nelem,nne,nnodes
  integer, parameter :: mxnod=20, mxelm=20000, mxpoi=20000 
  integer, parameter :: mxnow=20, mxelw=40000, mxpow=30000
  
  contains
    
    subroutine readMesh(file_mesh)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !                                                                                   !
      ! subrutina que lee todos los parametros de entrada para la simulacion,             !
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion          !
      !                                                                                   !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      implicit none
      
      character(len=*), parameter  :: fileplace = "Msh/"
      character(len=19)            :: file_mesh
      character(len=180)           :: msg
      double precision :: coorw(DimPr,mxpow), tempo(DimPr, mxpow)
      integer          :: ielem, jpoin, idime, i,j, stat
      integer          :: lnodw(mxelw,mxpow), lnod_add(mxelw,mxpow), tlnod(mxelw,mxpow)
      integer          :: npoiw,nelew,nnodw, npoif
      
      
      open(5, file=fileplace//file_mesh,status='old', action='read',IOSTAT=stat, IOMSG=msg)
      !open(5, file=fileplace//'inputCDR.dsc',status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      nelem = initElem
      nnodes= initNodes
      
      allocate( lnods(initElem,nne))
      allocate( coord(Dimpr,initNodes))
      lnods = 0.0
      coord = 0.0
      
      !do i=1,skipline
      !read(5,*) !se salta todas las lineas del input file hasta donde comienza la malla
      !end do
      do i=1,2
        read(5,*) !se salta todas las lineas del input file hasta donde comienza la malla
        if ( stat /= 0 )then
          print*, ' ' 
          print*, 'error in read mesh module geometry' 
          print'(A9,I3)','iostat= ',stat
          print'(A8,1x,A180)','iomsg= ',msg
          print'(A55,A)', 'error >>>> Something wrong during reading of mesh file ',file_mesh
          print*, ' ' 
          stop
        end if
      end do
      do i=1,nnodes !number of total nodes
        read(5,*,iostat=stat,iomsg=msg) jpoin,(coord(idime,jpoin), idime =1,DimPr )
        IF ( stat /= 0 )then
          print*,'iostat= ',stat
          print*, msg
        end if
      end do
      do i=1,3
        read(5,*) !se salta todas las lineas del input file hasta donde comienza la malla
        if ( stat /= 0 )then
          print*,'iostat= ',stat
          print*, msg
        endif
      end do
      do i=1,nelem
        read(5,*,iostat=stat,iomsg=msg) ielem,(lnods(ielem,j), j =1,nne)
        IF ( stat /= 0 )then
          print*,'iostat= ',stat
          print*, msg
        END IF
      end do
      
      close(5)
      
      if(refiType.eq.'NO')then
        ElemType = InitElemType
        nevab    = initnevab
        ntotv    = initntotv
        
        goto 101
        
      elseif(refiType.eq.'PS'.or.refitype.eq.'CB')then
        !
        !***  Undertakes the mesh change
        !
        call AddNodes(refiType,npoiw,nnodw,coorw,lnod_add)
        !
        !***  Checks if there are repeated nodes and output of results
        ! 
        call SplitElem(refiType,lnod_add,nelew,lnodw) 
        !
        !***  Checks if there are repeated nodes and reallocate coord and lnods
        ! 
        call checkMesh(coorw,lnodw,nnodw,nelew,npoiw,npoif,tempo,tlnod)
        
        !
        !* Recounting of nodes and elements after the refination
        !
        nnodes = npoif
        nelem  = nelew 
        nne    = nnodw
        
      end if
     
      101 continue
      
      nevab = ndofn*nne
      ntotv = ndofn*nnodes
      
    end subroutine readMesh
    !
    !*************************************************************************     
    !
    subroutine AddNodes(refiType,npoiw,nnodw,coorw,lnodw)
      !******************************************************************************
      !
      !**** This program changes the element type of a finite element mesh.
      !**** The possibilities are the following:      
      !
      !**** For NDIME = 2
      !      
      !****     NNODE = 3     ---->    NNODE = 7    (more nodes are created)
      !****     NNODE = 9     ---->    NNODE = 3    (same number of nodes, more elements)
      !
      !****     NNODE = 4     ---->    NNODE = 5    (more nodes are created)
      !****     NNODE = 5     ---->    NNODE = 3    (triangular elements of 3nodes)
      !
      !******************************************************************************
      !       
      character(len=2),intent(in)      :: refitype
      integer                                               :: idime, inode, ielem, ipoin 
      double precision, intent(out), dimension(DimPr,mxpow) :: coorw
      integer         , intent(out), dimension(mxelw,mxnow) :: lnodw
      integer         , intent(out) :: nnodw ,npoiw
      
      !***  Initializations
      do idime=1,DimPr
        do ipoin=1,nnodes
          coorw(idime,ipoin)=coord(idime,ipoin)
        end do
      end do
      
      npoiw=nnodes
      nnodw=3 !is a fix value because both PS & CB, has as a final element a 3-node ∆ element
      
      do ielem=1,nelem
        !*** 2D: NNODE = 3 --> NNODE = 4 6 & 7
        if(refiType.eq.'PS')then
          if(nne.ne.3)then 
            write(*,'(a)') 'PS refinment not compatible with quadrilateral element'
            write(*,'(a)') '>>> Verify PS must be nne=3'
            stop
          else
            continue
            !if (nnodw.ge.6) then 
              !Para construir el elemento de 7 nodos. Primero hay que construir el de 6
              !por eso la codicion dice mayor o igual, si pongo 7 entonces 
              !entrara al de 6 y luego al de 7
              do idime=1,2
                coorw(idime,npoiw+1)= 0.5*(coord(idime,lnods(ielem,1))+coord(idime,lnods(ielem,2)))
                coorw(idime,npoiw+2)= 0.5*(coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
                coorw(idime,npoiw+3)= 0.5*(coord(idime,lnods(ielem,3))+coord(idime,lnods(ielem,1)))
              end do
              do inode=1,3
                lnodw(ielem,  inode)=lnods(ielem,inode)
                lnodw(ielem,3+inode)=npoiw+inode
              end do
              npoiw=npoiw+3
            !if(nnodw.eq.7) then
              do idime=1,2
                coorw(idime,npoiw+1)= (1.0/3.0)*(coord(idime,lnods(ielem,1)) + &
                &                     coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
              end do
              lnodw(ielem,7)=npoiw+1
              npoiw=npoiw+1
            !endif
          endif
          
        elseif(refiType.eq.'CB')then
          if(nne.ne.4)then
            
            write(*,'(a)') 'CB refinment not compatible with triangular element'
            write(*,'(a)') '>>> Verify CB must be nne=4'
            stop
          else
            continue
            !  !Creation of one more node at the center of 4 node quadrilateral
            do idime=1,2
              coorw(idime,npoiw+1)= 0.25*(coord(idime,lnods(ielem,1))+coord(idime,lnods(ielem,2))&
              &                         +coord(idime,lnods(ielem,3))+coord(idime,lnods(ielem,4)))
            end do
            
            lnodw(ielem  ,1)=lnods(ielem,1)
            lnodw(ielem  ,2)=lnods(ielem,2)
            lnodw(ielem  ,3)=lnods(ielem,3)
            lnodw(ielem  ,4)=lnods(ielem,4)
            lnodw(ielem  ,5)=npoiw+1
            npoiw=npoiw+1
            
          endif
        else
          write(*,'(A)')' --No refinment selected-- '
          goto 101
        end if
        
      end do
      
      if(npoiw.gt.mxpow) then
        write(*,'(a,i5)') 'Increase MXPOW to ' , npoiw
        stop
      end if
      
      101 continue
      
    end subroutine AddNodes 
    !
    !*************************************************************************     
    !
    subroutine SplitElem(refiType,lnod_add,nelew,lnodw) 
      
      implicit none 
      
      character(len=2), intent(in) :: refitype
      integer         , intent(in), dimension(mxelw,mxnow) :: lnod_add !add nodes
      integer                                               :: ielem
      integer         , intent(out), dimension(mxelw,mxnow) :: lnodw
      integer         , intent(out)                         ::  nelew 
      
      nelew=nelem
      do ielem=1,nelem
        !*** 2D: NNODE = 3 --> NNODE = 4 6 & 7
        if(refiType.eq.'PS')then
          if(nne.ne.3)then 
            write(*,'(a)') 'PS refinment not compatible with quadrilateral element'
            write(*,'(a)') '>>> Verify PS must be nne=3'
            stop
          else
            continue
            !Split of a 7-node ∆ element into six 3-nodes ∆ elements. 
            !  ^
            !  |        3
            !  |        o
            !  |       / \
            !  |      /   \
            !  Y    6o  7  o5
            !  |    /   o   \
            !  |   /         \
            !  |  o-----o-----o
            !  |  1     4     2
            !  |
            !  +--------X-------->
            lnodw(ielem  ,1)=lnod_add(ielem,1)
            lnodw(ielem  ,2)=lnod_add(ielem,4)
            lnodw(ielem  ,3)=lnod_add(ielem,7)
            lnodw(nelew+1,1)=lnod_add(ielem,1)
            lnodw(nelew+1,2)=lnod_add(ielem,7)
            lnodw(nelew+1,3)=lnod_add(ielem,6)
            lnodw(nelew+2,1)=lnod_add(ielem,7)
            lnodw(nelew+2,2)=lnod_add(ielem,3)
            lnodw(nelew+2,3)=lnod_add(ielem,6)
            lnodw(nelew+3,1)=lnod_add(ielem,4)
            lnodw(nelew+3,2)=lnod_add(ielem,2)
            lnodw(nelew+3,3)=lnod_add(ielem,7)
            lnodw(nelew+4,1)=lnod_add(ielem,2)
            lnodw(nelew+4,2)=lnod_add(ielem,5)
            lnodw(nelew+4,3)=lnod_add(ielem,7)
            lnodw(nelew+5,1)=lnod_add(ielem,5)
            lnodw(nelew+5,2)=lnod_add(ielem,3)
            lnodw(nelew+5,3)=lnod_add(ielem,7)
            nelew=nelew+5
            
            ElemType = InitElemType
           
          endif
          
        elseif(refiType.eq.'CB')then
          if(nne.ne.4)then
            
            write(*,'(a)') 'CB refinment not compatible with triangular element'
            write(*,'(a)') '>>> Verify CB must be nne=4'
            stop
          else
            continue
            !  !Creation of one more node at the center of 4 node quadrilateral
            !  
            !  !Split of a 5-node element into four 3 nodes ∆ elements. 
            !  !  ^
            !  !  | 4       3
            !  !  | o-------o
            !  !  | | \   / |
            !  !  | |  \5/  |
            !  !  Y |   o   |        Is a square compund by 4 three-nodes triangle      
            !  !  | |  / \  |
            !  !  | | /   \ |
            !  !  | o-------o
            !  !  | 1       2
            !  !  +------X------->
            !  
            lnodw(ielem  ,1)=lnod_add(ielem,1)
            lnodw(ielem  ,2)=lnod_add(ielem,2)
            lnodw(ielem  ,3)=lnod_add(ielem,5)
            lnodw(nelew+1,1)=lnod_add(ielem,1)
            lnodw(nelew+1,2)=lnod_add(ielem,5)
            lnodw(nelew+1,3)=lnod_add(ielem,4)
            lnodw(nelew+2,1)=lnod_add(ielem,5)
            lnodw(nelew+2,2)=lnod_add(ielem,3)
            lnodw(nelew+2,3)=lnod_add(ielem,4)
            lnodw(nelew+3,1)=lnod_add(ielem,2)
            lnodw(nelew+3,2)=lnod_add(ielem,3)
            lnodw(nelew+3,3)=lnod_add(ielem,5)
            nelew=nelew+3
            
            ElemType = 'TRIA'
            
          endif
        else
          write(*,'(A)')' --No refinment selected-- '
          goto 101
        end if                                            ! nne
        
        
      end do                                              ! ielem=1,nelem
      
      if(nelew.gt.mxelw) then
        write(*,'(a,i5)') 'Increase MXELW to ' , nelew
        stop
      end if
      
      101 continue
      
    end subroutine SplitElem
    !
    !*************************************************************************     
    !
    subroutine checkMesh(coorw,lnodw,nnodw,nelew,npoiw,npoif,tempo,tlnod)
      
      !implicit     real*8 (a-h,o-z)
      
      implicit none !integer :: (i,j)
      
      integer         , intent(in) :: lnodw(mxelw,mxnow)
      double precision, intent(in)  :: coorw(DimPr,mxpow)
      integer         , intent(in)  :: nnodw, nelew, npoiw
      integer                       :: ipoiw, ielew, inodw, idime, kpoiw, jpoiw, lpoiw(npoiw)
      double precision              :: x0, x1, y0, y1, z1, z0, dist
      double precision, intent(out) :: tempo(DimPr,npoiw) !Total EleMent POint
      integer         , intent(out) :: tlnod(mxelw,mxnow) !Total List of NODes
      integer         , intent(out) :: npoif
      
      integer :: ielem, ipoin, inode
      
      !dimension :: coorw(DimPr,mxpow),  
      !double precision, intent(out), dimension :: tempo(DimPr,npoiw)
      !integer, intent(out), dimension          :: lnodw(mxelw,mxnow)
     
      !allocate(lnodw(mxelw,mxnow))

      do ipoiw=1,npoiw
        lpoiw(ipoiw)=ipoiw
      end do
     
      z1=0.0
      z0=0.0
      do ipoiw=nnodes+2,npoiw
        x1=coorw(1,ipoiw)
        y1=coorw(2,ipoiw)
        if(DimPr.eq.3) z1=coorw(3,ipoiw)
        do jpoiw=nnodes+1,ipoiw-1
          x0=coorw(1,jpoiw)
          y0=coorw(2,jpoiw)
          if(DimPr.eq.3) z0=coorw(3,jpoiw)
          dist=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
          if(dist.lt.1.0e-7) then
            lpoiw(ipoiw)=-abs(lpoiw(jpoiw))
            do kpoiw=ipoiw+1,npoiw
              lpoiw(kpoiw)=lpoiw(kpoiw)-1
            end do
          end if
        end do
      end do
      
      npoif=0
      do ipoiw=1,npoiw
        if(lpoiw(ipoiw).gt.0) then
          npoif=npoif+1
          do idime=1,DimPr
            tempo(idime,npoif)=coorw(idime,ipoiw)
          end do
        end if
      end do
      
      do ielew=1,nelew
        do inodw=1,nnodw
          tlnod(ielew,inodw)=abs(lpoiw(lnodw(ielew,inodw)))
        end do
      end do
      
      !* Setting the refinment mesh in previous coord and lnods (global)) variables 
      deallocate(coord, lnods)
      allocate(coord(DimPr,npoif),lnods(nelew,nnodw))
     
      do ielem=1,nelew
        do inode=1,nnodw
          lnods(ielem,inode)=tlnod(ielem,inode)
        end do
      end do
      
      do ipoin=1,npoif
          do idime=1,DimPr
            coord(idime,ipoin)=tempo(idime,ipoin)
          end do
      end do
      
      !write(*,'(a)') 'ELEMENTS'
      !do ielem=1,nelew
      !  write(*,10) ielem,(lnods(ielem,inode),inode=1,nnodw)
      !end do
      !write(*,'(a)') 'END_ELEMENTS'
     
      !write(*,'(a)') 'COORDINATES'
      !do ipoin=1,npoif
      !  write(*,20) ipoin,(coord(idime,ipoin),idime=1,DimPr)
      !end do
      !write(*,'(a)') 'END_COORDINATES'
      ! 
      !10 format(1x,i5,10(1x,i5))
      !20 format(5x,i6,3(2x,f16.9))
      
    end subroutine checkMesh
    
  !end contains 
    
end module geometry
