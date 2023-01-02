module refin
  use param

  contains
    !program Changenod
    !    c******************************************************************************
    !    c
    !    c**** This program changes the element type of a finite element mesh.
    !    c**** The possibilities are the following:      
    !    c
    !    c**** For NDIME = 2
    !    c      
    !    c****     NNODE = 3     ---->    NNODE = 4, 6 & 7    (more nodes are created)
    !    c****     NNODE = 9     ---->    NNODE = 3           (same number of nodes)
    !    c
    !    c**** The data from the initial mesh is supposed to be read from the .geo
    !    c**** file of FANTOM. The results are also presented in this format.
    !    c
    !    c******************************************************************************
    !   
    !    implicit     real*8 (a-h,o-z)
    !    common/contr/nin,nou,ndime,nelem,nnode,npoin
    !    parameter
    !    .  (mxnod=20, mxelm=20000, mxpoi=20000, mxdim=3,
    !    .   mxnow=20, mxelw=40000, mxpow=30000)
    !     dimension
    !    .  lnods(mxelm,mxnod), lnodw(mxelw,mxnow),
    !    .  coord(mxdim,mxpoi), coorw(mxdim,mxpow),
    !    .  tempo(mxdim,mxpow), lpoiw(      mxpow)
    !    
    !    !***  Open files
    !    !
    !          call opfile
    !    !
    !    !***  Identify control parameters
    !    !
    !          call contro
    !    !
    !    !***  Read geometry (initial mesh)
    !    !
    !          call reageo(lnods,coord)
    !    !
    !    !***  Undertakes the mesh change
    !    !
    !          call mesdiv(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,npoiw,nelew,nnodw)
    !    !
    !    !***  Checks if there are repeated nodes and output of results
    !    !
    !          call mescek(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,tempo,lpoiw,npoiw,nelew,nnodw)
    !          stop
    !end program Changenod
   
    ! 
    !*************************************************************************     
    !
    subroutine opfile
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      character*20 file_mesh, file_resu
      
      write(6,'(a,$)') ' >>> Original mesh file:'
      read(5,'(a)') file_mesh
      write(6,'(a,$)') ' >>> Final mesh file:'
      read(5,'(a)') file_resu
      nin=7
      nou=8
      open(nin,file=file_mesh,status='old')
      open(nou,file=file_resu,status='unknown')
    end opfile
    !
    !*************************************************************************     
    !
    subroutine contro
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      character*80 string
     
      ndime=-1
      i_chk= 1
      nelem= 0
      npoin= 0
      icoun= 0
      nnode=-1
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
           if(kcoun.eq.1) nnode = nnode+1
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
           if(kcoun.eq.1) ndime = ndime+1
           kcoun = kcoun+1
           jcoun = jcoun+1
         enddo
      if(jcoun.lt.80) goto 30
        do while(string(1:7).ne.'END_COO')
          read(nin,'(a80)') string
          npoin=npoin+1
        enddo
      end if
      
      nelem=nelem/i_chk
      rewind(nin)
     
    end subroutine contro 
    !
    !*************************************************************************     
    !
    subroutine reageo(lnods,coord)
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      dimension    coord(ndime,npoin), lnods(nelem,nnode)
      
      read(nin,*)
      do ielem = 1,nelem
        read(nin,*) jelem,(lnods(jelem,inode),inode=1,nnode)
      enddo
      read(nin,*)
      
      read(nin,*)
      do ipoin = 1,npoin
        read(nin,*) jpoin,(coord(idime,jpoin),idime=1,ndime)
      enddo
      
    end subroutine contro
    !
    !*************************************************************************     
    !
    subroutine mesdiv( mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,npoiw,nelew,nnodw)
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      dimension    coord(ndime,npoin), lnods(nelem,nnode), coorw(ndime,mxpow), lnodw(mxelw,mxnow)
      
      !***  Initializations
      do idime=1,ndime
        do ipoin=1,npoin
          coorw(idime,ipoin)=coord(idime,ipoin)
        end do
      end do
      
      nelew=nelem
      npoiw=npoin
      
      !***  Splits the elements
      write(6,'(a,$)') ' >>> Nodes of the final elements:'
      read(5,*) nnodw
      
      do ielem=1,nelem
        if(ndime.eq.2) then
        !*** 2D: NNODE = 3 --> NNODE = 4 6 & 7
          if(nnode.eq.3) then
            if (nnodw.ge.6) then !Para construir el elemento de 7 nodos. Primero hay que construir el de 6
              ! por eso la codicion dice mayor o igual, si pongo 7 entonces entrara al de 6 y luego al de 7
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
              if(nnodw.eq.7) then
                do idime=1,2
                  coorw(idime,npoiw+1)= (1.0/3.0)*(coord(idime,lnods(ielem,1)) +
                                        coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
                end do
                lnodw(ielem,7)=npoiw+1
                npoiw=npoiw+1
            end if
          !   *** 2D:  NNODE = 6 --> NNODE = 7
          else if (nnode.eq.6) then             !mismo numero de nodos pero mas elementos
            if (nnodw.eq.3) then
              lnodw(ielem,1)=lnods(ielem,1)
              lnodw(ielem,2)=lnods(ielem,4)
              lnodw(ielem,3)=lnods(ielem,6)
              nelew=nelew+1
              lnodw(nelew,1)=lnods(ielem,4)
              lnodw(nelew,2)=lnods(ielem,2)
              lnodw(nelew,3)=lnods(ielem,5)
              nelew=nelew+1
              lnodw(nelew,1)=lnods(ielem,4)
              lnodw(nelew,2)=lnods(ielem,5)
              lnodw(nelew,3)=lnods(ielem,6)
              nelew=nelew+1
              lnodw(nelew,1)=lnods(ielem,6)
              lnodw(nelew,2)=lnods(ielem,5)
              lnodw(nelew,3)=lnods(ielem,3)
            else if (nnodw.eq.7) then
              do idime=1,2
                coorw(idime,npoiw+1)= (1.0/3.0)*(coord(idime,lnods(ielem,1))+&
                 &                     coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
              end do
              do inode=1,6
                lnodw(ielem,inode)=lnods(ielem,inode)
              end do
              lnodw(ielem,7)=npoiw+1
              npoiw=npoiw+1
            end if
        !***  2D: NNODE = 9 --> NNODE = 3
            else if(nnode.eq.9) then
              
              if(nnodw.eq.3) then
                lnodw(ielem  ,1)=lnods(ielem,1)
                lnodw(ielem  ,2)=lnods(ielem,5)
                lnodw(ielem  ,3)=lnods(ielem,9)
                lnodw(nelew+1,1)=lnods(ielem,1)
                lnodw(nelew+1,2)=lnods(ielem,9)
                lnodw(nelew+1,3)=lnods(ielem,8)
                lnodw(nelew+2,1)=lnods(ielem,8)
                lnodw(nelew+2,2)=lnods(ielem,9)
                lnodw(nelew+2,3)=lnods(ielem,7)
                lnodw(nelew+3,1)=lnods(ielem,8)
                lnodw(nelew+3,2)=lnods(ielem,7)
                lnodw(nelew+3,3)=lnods(ielem,4)
                lnodw(nelew+4,1)=lnods(ielem,5)
                lnodw(nelew+4,2)=lnods(ielem,2)
                lnodw(nelew+4,3)=lnods(ielem,6)
                lnodw(nelew+5,1)=lnods(ielem,5)
                lnodw(nelew+5,2)=lnods(ielem,6)
                lnodw(nelew+5,3)=lnods(ielem,9)
                lnodw(nelew+6,1)=lnods(ielem,9)
                lnodw(nelew+6,2)=lnods(ielem,6)
                lnodw(nelew+6,3)=lnods(ielem,3)
                lnodw(nelew+7,1)=lnods(ielem,9)
                lnodw(nelew+7,2)=lnods(ielem,3)
                lnodw(nelew+7,3)=lnods(ielem,7)
                nelew=nelew+7
              endif
              
            end if
           
          end if                                            ! nnode
        end if                                              ! ndime
      end do                                                ! ielem=1,nelem
      
      if(nelew.gt.mxelw) then
        write(6,'(a,i5)') 'Increase MXELW to ' , nelew
        stop
      end if
      if(npoiw.gt.mxpow) then
        write(6,'(a,i5)') 'Increase MXPOW to ' , npoiw
        stop
      end if
    end subroutine mesdiv 
    !
    !*************************************************************************     
    !
    subroutine mescek(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,tempo,lpoiw,npoiw,nelew,nnodw)
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      dimension :: coord(ndime,npoin), lnods(nelem,nnode), coorw(ndime,mxpow)
      dimension :: lnodw(mxelw,mxnow), tempo(ndime,npoiw), lpoiw(npoiw)
       
      do ipoiw=1,npoiw
        lpoiw(ipoiw)=ipoiw
      end do
     
      z1=0.0
      z0=0.0
      do ipoiw=npoin+2,npoiw
        x1=coorw(1,ipoiw)
        y1=coorw(2,ipoiw)
        if(ndime.eq.3) z1=coorw(3,ipoiw)
        do jpoiw=npoin+1,ipoiw-1
          x0=coorw(1,jpoiw)
          y0=coorw(2,jpoiw)
          if(ndime.eq.3) z0=coorw(3,jpoiw)
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
          do idime=1,ndime
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
        if(nnodw.eq.16) then
          write(nou,11) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        else if (nnodw.gt.16) then
          write(nou,12) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        else
          write(nou,10) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        end if
      end do
      write(nou,'(a)') 'END_ELEMENTS'
    
      write(nou,'(a)') 'COORDINATES'
      do ipoin=1,npoif
        write(nou,20) ipoin,(tempo(idime,ipoin),idime=1,ndime)
      end do
      write(nou,'(a)') 'END_COORDINATES'
      
      10 format(1x,i5,10(1x,i5))
      11 format(1x,i3,16(1x,i3))      
      12 format(1x,i5,10(1x,i5)," /",/,6x,10(1x,i5))
      20 format(5x,i6,3(2x,e16.9))
    end mescek
    
  !end contains 
 
end module refin 
