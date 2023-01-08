module refine
  use param
  implicit none
  !common/contr/nin,nou,DimPr,nelem,nne,nnodes
  integer, parameter :: mxnod=20, mxelm=20000, mxpoi=20000 
  integer, parameter :: mxdim=3, mxnow=20, mxelw=40000, mxpow=30000
  integer, :: nin, nou
  
  
  contains
    
    
    subroutine MacroElement(reftype, coordRef, lnodsRef)
      !******************************************************************************
      !
      !**** This program changes the element type of a finite element mesh.
      !**** The possibilities are the following:      
      !
      !**** For NDIME = 2
      !      
      !****     NNODE = 3     ---->    NNODE = 4, 6 & 7    (more nodes are created)
      !****     NNODE = 9     ---->    NNODE = 3           (same number of nodes)
      !
      !**** The data from the initial mesh is supposed to be read from the .geo
      !**** file of FANTOM. The results are also presented in this format.
      !
      !******************************************************************************
      
      double precision, intent(in), dimension :: coord(mxdim,mxpoi)
      integer,          intent(in), dimension :: lnods(mxelm,mxnod)
      double precision, dimension :: tempo(mxdim,mxpow)
      integer, dimension          :: lpoiw(mxpow)
      double precision, intent(out), dimension :: coorw(mxdim,mxpow)
      integer,          intent(out), dimension :: lnodw(mxelw,mxnow)
      
      write(6,'(a,$)') ' >>> Original mesh file:'
      read(5,'(a)') file_mesh
      write(6,'(a,$)') ' >>> Final mesh file:'
      read(5,'(a)') file_resu
      nin=7
      nou=8
      open(nin,file=file_mesh,status='old')
      open(nou,file=file_resu,status='unknown')
      
      !Esta rutina la usare para leer la malla como yo la tengo a como la necesita este modulo, es decir
      !
      !      local                     DREMm
      !lnods(nelem,nne) --------> lnods(nelem,nne)
      !coord(DimPr,nnodes) -----> coord(nnodes,DimPr)
      
      read(nin,*) !esta linea esta puest apara leer el el encabezado ELEMENTS
      do ielem = 1,nelem
        read(nin,*) jelem,(lnods(jelem,inode),inode=1,nne)
      enddo
      read(nin,*)! esta linea lee END ELEMENTS
      
      read(nin,*)! esta linea lee COORDINATES
      do ipoin = 1,nnodes
        read(nin,*) jpoin,(coord(idime,jpoin),idime=1,DimPr)
      enddo
      
      !***  Open files
      !
      !     call opfile
      !
      !***  Identify control parameters
      !
      !     call contro
      !
      !***  Read geometry (initial mesh)
      !
      !     call reageo(lnods,coord)
      !
      !***  Undertakes the mesh change
      !
      call mesdiv(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,npoiw,nelew,nnodw)
      !
      !***  Checks if there are repeated nodes and output of results
      !
      call mescek(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,tempo,lpoiw,npoiw,nelew,nnodw)
      stop
    end program MacroElement
   
    !
    !*************************************************************************     
    !
    subroutine contro
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,DimPr,nelem,nne,nnodes
      character*80 string
     
      DimPr=-1
      i_chk= 1
      nelem= 0
      nnodes= 0
      icoun= 0
      nne=-1
      read(nin,*)
      20 jcoun= 1
      
      read(nin,'(a80)') string
      10 do while(string(jcoun:jcoun).eq.' '.and.jcoun.le.80)
           jcoun = jcoun+1
         enddo
         kcoun = 1
         do while(string(jcoun:jcoun).ne.' '.and.jcoun.le.80)
           if(string(jcoun:jcoun).eq.'/') then
             nelem = nelem+1
             i_chk = i_chk+1
             goto 20
           endif
           if(kcoun.eq.1) nne = nne+1
           kcoun = kcoun+1
           jcoun = jcoun+1
         enddo
         if(jcoun.lt.80) goto 10
         do while(string(1:7).ne.'END_ELE')
           read(nin,'(a80)') string
           nelem=nelem+1
         enddo
      read(nin,*)
      read(nin,'(a80)') string
      jcoun= 1
      30 do while(string(jcoun:jcoun).eq.' '.and.jcoun.le.80)
           jcoun = jcoun+1
         enddo
         kcoun = 1
         do while(string(jcoun:jcoun).ne.' '.and.jcoun.le.80)
           if(kcoun.eq.1) DimPr = ndime+1
           kcoun = kcoun+1
           jcoun = jcoun+1
         enddo
      if(jcoun.lt.80) goto 30
        do while(string(1:7).ne.'END_COO')
          read(nin,'(a80)') string
          nnodes=nnodes+1  !npoin=npoin+1
        enddo
      end if
      
      nelem=nelem/i_chk
      rewind(nin)
     
    end subroutine contro 
    !
    !*************************************************************************     
    !
    subroutine mesdiv( mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,npoiw,nelew,nnodw)
      implicit     real*8 (a-h,o-z)
      !common/contr/nin,nou,DimPr,nelem,nne,nnodes
      dimension    coord(DimPr,nnodes), lnods(nelem,nne), coorw(ndime,mxpow), lnodw(mxelw,mxnow)
      
      !***  Initializations
      do idime=1,DimPr
        do ipoin=1,nnodes
          coorw(idime,ipoin)=coord(idime,ipoin)
        end do
      end do
      
      nelew=nelem
      npoiw=nnodes
      
      !***  Splits the elements
      write(6,'(a,$)') ' >>> Nodes of the final elements:'
      read(5,*) nnodw
      
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
          
          !***  2D: NNODE = 7 --> NNODE = 3
          !Este If deberia salir del condicional para, despues de tener el elemento de 7 nodos
          !dividirlo en 6 elementos de 3 nodos cada uno (powell Sabin)
          
          !Division of a 7 node element into six 3 nodes elements. 
          !  ^
          !  |        3
          !  |        o
          !  |       / \
          !  |      /   \
          !  Y   6 o  7  o 5
          !  |    /   o   \
          !  |   /    |    \
          !  |  o-----o-----o
          !  |  1     4     2
          !  |
          !  +--------X-------->
          lnodw(ielem  ,1)=lnods(ielem,1)
          lnodw(ielem  ,2)=lnods(ielem,4)
          lnodw(ielem  ,3)=lnods(ielem,7)
          lnodw(nelew+1,1)=lnods(ielem,1)
          lnodw(nelew+1,2)=lnods(ielem,7)
          lnodw(nelew+1,3)=lnods(ielem,6)
          lnodw(nelew+2,1)=lnods(ielem,7)
          lnodw(nelew+2,2)=lnods(ielem,3)
          lnodw(nelew+2,3)=lnods(ielem,6)
          lnodw(nelew+3,1)=lnods(ielem,4)
          lnodw(nelew+3,2)=lnods(ielem,2)
          lnodw(nelew+3,3)=lnods(ielem,7)
          lnodw(nelew+4,1)=lnods(ielem,2)
          lnodw(nelew+4,2)=lnods(ielem,5)
          lnodw(nelew+4,3)=lnods(ielem,7)
          lnodw(nelew+5,1)=lnods(ielem,5)
          lnodw(nelew+5,2)=lnods(ielem,3)
          lnodw(nelew+5,3)=lnods(ielem,7)
          nelew=nelew+5
          
        elseif(refiType.eq.'CB')then
          if(nne.ne.4)then
            
            write(*,'(a)') 'CB refinment not compatible with triangular element'
            write(*,'(a)') '>>> Verify CB must be nne=4'
            stop
          else
            continue
            !  !Creation of one more node at the center of 4 node quadrilateral
            !  do idime=1,2
            !    coorw(idime,npoiw+1)= 0.25*(coord(idime,lnods(ielem,1))+coord(idime,lnods(ielem,2))&
            !    &                         +coord(idime,lnods(ielem,3))+coord(idime,lnods(ielem,4)))
            !  end do
            !  
            !  
            !  !Division of a 5 node element into four 3 nodes elements. 
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
            !  lnodw(ielem  ,1)=lnods(ielem,1)
            !  lnodw(ielem  ,2)=lnods(ielem,2)
            !  lnodw(ielem  ,3)=lnods(ielem,5)
            !  lnodw(nelew+1,1)=lnods(ielem,1)
            !  lnodw(nelew+1,2)=lnods(ielem,5)
            !  lnodw(nelew+1,3)=lnods(ielem,4)
            !  lnodw(nelew+2,1)=lnods(ielem,5)
            !  lnodw(nelew+2,2)=lnods(ielem,3)
            !  lnodw(nelew+2,3)=lnods(ielem,4)
            !  lnodw(nelew+3,1)=lnods(ielem,2)
            !  lnodw(nelew+3,2)=lnods(ielem,3)
            !  lnodw(nelew+3,3)=lnods(ielem,5)
            !  nelew=nelew+3
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
        end if                                            ! nne
        
        
      end do                                              ! ielem=1,nelem
      
      if(nelew.gt.mxelw) then
        write(6,'(a,i5)') 'Increase MXELW to ' , nelew
        stop
      end if
      if(npoiw.gt.mxpow) then
        write(6,'(a,i5)') 'Increase MXPOW to ' , npoiw
        stop
      end if
      
      101 print*, 'no_refinment'
      continue
      
    end subroutine mesdiv 
    !
    !*************************************************************************     
    !
    subroutine mescek(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,tempo,lpoiw,npoiw,nelew,nnodw)
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,DimPr,nelem,nne,nnodes
      dimension :: coord(DimPr,nnodes), lnods(nelem,nne), coorw(ndime,mxpow)
      dimension :: lnodw(mxelw,mxnow), tempo(DimPr,npoiw), lpoiw(npoiw)
       
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
          lnodw(ielew,inodw)=abs(lpoiw(lnodw(ielew,inodw)))
        end do
      end do
    
      write(nou,'(a)') 'ELEMENTS'
      do ielem=1,nelew
        !if(nnodw.eq.16) then
        !  write(nou,11) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        !else if (nnodw.gt.16) then
        !  write(nou,12) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        !else
        !  write(nou,10) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        !end if
        write(nou,10) ielem,(lnodw(ielem,inode),inode=1,nnodw)
      end do
      write(nou,'(a)') 'END_ELEMENTS'
     
      write(nou,'(a)') 'COORDINATES'
      do ipoin=1,npoif
        write(nou,20) ipoin,(tempo(idime,ipoin),idime=1,DimPr)
      end do
      write(nou,'(a)') 'END_COORDINATES'
      
      10 format(1x,i5,10(1x,i5))
      !11 format(1x,i3,16(1x,i3))      
      !12 format(1x,i5,10(1x,i5)," /",/,6x,10(1x,i5))
      20 format(5x,i6,3(2x,f16.9))
    end mescek
    
  !end contains 
 
end module refine
