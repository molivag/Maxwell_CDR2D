module solver
  use param


  contains
    
    subroutine solsistem(amatr,cmatr,nsist,bvect,xvect)
      !*****************************************************************************
      !
      !     Resolucio de sistemes lineals amb matriu en banda, simetriques i  
      !     definides positives emprant la descomposicio de Cholesky.
      !
      !             A . x = b                  L . z = b
      !                             -----> 
      !             L . Lt= A                   Lt. x = z
      !
      !     nequa = num. d'equacions = num.d'incognites ----> es ntotv y es variable global por eso no es dummy
      !     nsban = semiample de banda   ----> es nband y es variable global por eso no es dummy
      !     amatr = triangle superior de la matriu de coeficients emmagatzemat
      !             en banda per files. Si A es la matriu original,
      !             amatr(i,j) = A(j,i+j-1)
      !     cmatr = matriu L transposta de la descomposicio de Cholesky, 
      !             emmagatzemada en banda per files
      !     nsist = num. de sistemes que cal resoldre
      !     bvect = vector de termes independents
      !     xvect = vector d'incognites
      !
      !*****************************************************************************
      implicit real*8(a-h,o-z)
      double precision, intent(in) :: amatr(nband+1,ntotv),  bvect(ntotv,1)
      double precision, intent(inout) :: cmatr(nband+1,ntotv)
      integer                                   :: nsist, nequa, nsban
      double precision, allocatable, intent(out):: xvect(:,:)
      
      allocate ( xvect(ntotv,nsist) )
      
      nequa = ntotv
      nsban = nband
      
      call cholesky(nequa,nsban,amatr,cmatr)
      call sub1 (nequa,nsban,cmatr,nsist,bvect,xvect)
      call sub2 (nequa,nsban,cmatr,nsist,xvect,xvect)
      return
    end subroutine solsistem
    
    !   Descomposicio de Cholesky 
    subroutine cholesky(n,l,a,c)
      
      implicit double precision(a-h,o-z)
      double precision, intent(in) :: a(nband+1,ntotv) 
      integer                                   :: n, l, i, j, j2, k, k1
      double precision, intent(out):: c(l+1,n)
      
      
      do i = 1,n
        c(1,i) = a(1,i)
        k1 = i-l
        if (k1.lt.1) k1 = 1
        if (k1.le.i-1) then 
          do k = k1,i-1
            c(1,i) = c(1,i) - c(i-k+1,k)*c(i-k+1,k)
          end do
        end if
        if (c(1,i).le.0.) call error(' Element no positiu a la fila ',i)
        c(1,i) = dsqrt(c(1,i))
        j2 = i + l
        if (j2.gt.n) j2 = n
        if (i+1.le.j2) then 
          do j = i+1,j2
            c(j-i+1,i) = a(j-i+1,i)
            k1 = j - l
            if (k1.lt.1) k1 = 1
            if (k1.le.i-1) then 
              do k = k1,i-1
                c(j-i+1,i) = c(j-i+1,i) - c(i-k+1,k)*c(j-k+1,k)
              end do
            end if
            c(j-i+1,i) = c(j-i+1,i)/c(1,i)
          end do
        end if
      end do
      
      return
    end subroutine cholesky
    
    ! Substitucio endavant
    subroutine sub1(n,l,c,m,b,z)
      
      implicit real*8(a-h,o-z)
      real*8 c(l+1,n),b(n,m),z(n,m)
      
      do is = 1,m
        z(1,is) = b(1,is)/c(1,1)
        do i = 2,n
          z(i,is) = b(i,is)
          k1 = i-l
          if (k1.lt.1) k1 = 1
          if (k1.le.i-1) then 
            do k = k1,i-1
              z(i,is) = z(i,is)-c(i-k+1,k)*z(k,is)
            end do
          end if
          if (c(1,i).le.0.) call error(' Pivot no positiu a la fila ',i)
          z(i,is) = z(i,is)/c(1,i)
        end do
      end do
      
      return
    end subroutine sub1
    
    !   Substitucio enrere
    subroutine sub2(n,l,c,m,z,x)
      
      implicit real*8(a-h,o-z)
      real*8  c(l+1,n),z(n,m),x(n,m)
      
      do is = 1,m
        x(n,is) = z(n,is)/c(1,n)
        do i = n-1,1,-1
          x(i,is) = z(i,is)
          k2 = i+l
          if (k2.gt.n) k2 = n
          if (i+1.le.k2) then 
            do k = i+1,k2
              x(i,is) = x(i,is)-c(k-i+1,i)*x(k,is)
            end do
          end if
          if (c(1,i).le.0.) call error(' Pivot no positiu a la fila ',i)
          x(i,is) = x(i,is)/c(1,i)
        end do
      end do
      return
    end subroutine sub2
   !
    ! Sortida de possibles errors
    subroutine error(alfa,i)
      
      implicit real*8(a-h,o-z)
      integer :: i
      character alfa*(*)
      write(*,'(a,i5)') alfa,i
      stop
      
    end subroutine error
    !
    !
    !subroutine outdes(tdisp)
    !  !*****************************************************************************
    !  !
    !  !   Escriptura dels moviments de l'estructura a l'arxiu de resultats
    !  !
    !  !*****************************************************************************
    !  implicit real*8(a-h,o-z)
    !  character c_coor*12,c_prop*14,c_dofn*9,c_load*13,c_stre*10
    !  common /contrl/ npoin,nelem,nmats,nvfix,nload,nband,ntotv
    !  common /arxius/ n_lect,n_escr,n_kltt,n_tran
    !  common /titols/ c_coor(2),c_prop(3),c_dofn(3),c_load(3),c_stre(3)
    !  real*8 tdisp(ntotv)
    !  
    !  write(n_escr,100) (c_dofn(i),i=1,3)
    !  100 format(//' MOVIMENTS DELS NODES DE L''ESTRUCTURA',/,
    !           ' --------- ---- ----- -- ------------',//,
    !           '        Num. de node',3(5x,a9))
    !  do ipoin=1,npoin
    !    itotv=3*(ipoin-1)
    !    write(n_escr,10) ipoin,(tdisp(itotv+i),i=1,3)
    !  end do
    !  10  format(i20,3(f14.8))
    !  
    !  return
    !end subroutine outdes
    !  
    !  
    !subroutine calesf(tdisp,treac,etmat,tmatr,lnods,matno, nofix,ifpre)
    !  !*****************************************************************************
    !  !
    !  !   Calcula i escriu els esforcos d'extrem de  barra i les reaccions
    !  !
    !  !*****************************************************************************
    !  implicit real*8(a-h,o-z)
    !  character c_coor*12,c_prop*14,c_dofn*9,c_load*13,c_stre*10
    !  common /contrl/ npoin,nelem,nmats,nvfix,nload,nband,ntotv
    !  common /arxius/ n_lect,n_escr,n_kltt,n_tran
    !  common /titols/ c_coor(2),c_prop(3),c_dofn(3),c_load(3),c_stre(3)
    !  real*8    tdisp(ntotv), treac(3,nvfix), etmat(6,6), tmatr(3,3), sextr(6), sglob(6) 
    !  integer*4 lnods(2,nelem),matno(nelem),nofix(nvfix), ifpre(3,nvfix)
    !  
    !  !***  Rebobina els arxius de Kl.Tt i T
    !  rewind(n_kltt)
    !  rewind(n_tran)
    !  
    !  !***  Titols per escriure els esforcos d'extrem de barra
    !  write(n_escr,100) (c_stre(i),i=1,3)
    !  100 format(//,' ESFORCOS EN EXTREMS DE BARRA',/,' -------- -- ------- -- -----',&
    !  //,'  Barra','            ',3(5x,a10) )
    !  
    !  !***  Llac sobre les barres
    !  do ielem=1,nelem
    !    !***  Calcul dels esforcos d'extrem de barra (P_l = K_l . T^t . d_g): SEXTR
    !    read(n_kltt) ((etmat(i,j),i=1,6),j=1,6)
    !    do inode=1,2
    !      do idofn=1,3
    !        stres=0.0
    !        ievab=(inode-1)*3+idofn
    !        do jnode=1,2
    !          jpoin=lnods(jnode,ielem)
    !          do jdofn=1,3
    !            jevab=(jnode-1)*3+jdofn
    !            jtotv=(jpoin-1)*3+jdofn
    !            stres=stres+etmat(ievab,jevab)*tdisp(jtotv)
    !          end do
    !        end do
    !        sextr(ievab)=stres
    !      end do
    !    end do
    !    ! 
    !    !***  Escriptura dels esforcos d'extrem de barra en coordenades locals
    !    !
    !    ipoi1=lnods(1,ielem)
    !    ipoi2=lnods(2,ielem)
    !    write(n_escr,10) ielem,'Node 1',ipoi1,(sextr(i),i=1,3)
    !    write(n_escr,11)       'Node 2',ipoi2,(sextr(i),i=4,6)
    !    !
    !    !***  Calcul de reaccions en el cas de barres coaccionades
    !    !
    !    if (matno(ielem).lt.0) then
    !      read(n_tran) ((tmatr(i,j),i=1,3),j=1,3)
    !    !
    !    !***  Esforcos d'extrem de barra en coordenades globals: SGLOB
    !    !
    !      do i=1,3
    !        sglob(i  )=0.0
    !        sglob(i+3)=0.0
    !        do j=1,3
    !          sglob(i  )=sglob(i  )+tmatr(i,j)*sextr(j  )
    !          sglob(i+3)=sglob(i+3)+tmatr(i,j)*sextr(j+3)
    !        end do
    !      end do 
    !    !
    !    !***  Identificacio del g.d.ll. prescrit i obtencio de la reaccio
    !    !
    !      do ivfix=1,nvfix
    !        if(nofix(ivfix).eq.ipoi1) then
    !          do i=1,3
    !            if(ifpre(i,ivfix).eq.1)&
    !            treac(i,ivfix)=treac(i,ivfix)+sglob(i  )
    !          end do
    !        else if(nofix(ivfix).eq.ipoi2) then
    !          do i=1,3
    !            if(ifpre(i,ivfix).eq.1)&
    !            treac(i,ivfix)=treac(i,ivfix)+sglob(i+3)
    !          end do
    !        end if
    !      end do
    !    end if
    !  end do
    !  !
    !  !***  Escriptura de les reaccions
    !  !
    !  write(n_escr,101) (c_load(i),i=1,3)
    !  101 format(//' REACCIONS EN NODES COACCIONATS',/,
    !    &         ' --------- -- ----- -----------',//,
    !    &         '    Node prescrit',3(7x,a13))
    !    
    !    do ivfix=1,nvfix
    !      ipoin=nofix(ivfix)
    !      write(n_escr,12) ipoin,(treac(i,ivfix),i=1,3)
    !    end do
    !  !
    !  !***  Formats
    !  !
    !  10  format(i7,2x,a6,i4,3f15.8)
    !  11  format(7x,2x,a6,i4,3f15.8)
    !  12  format(i17,3f20.8)
    !  
    !  return
    !end subroutine calesf
    !
    !subroutine tanarx
    !  !*****************************************************************************
    !  !
    !  !   Tancament dels arxius emprats en el programa
    !  !
    !  !*****************************************************************************
    !  implicit real*8(a-h,o-z)
    !  common /arxius/ n_lect,n_escr,n_kltt,n_tran
    !  close(unit=n_lect,status='keep')
    !  close(unit=n_escr,status='keep')
    !  close(unit=n_kltt,status='delete')
    !  close(unit=n_tran,status='delete')
    !  
    !  return
    !end subroutine tenarx
    
  ! end contains      
  
end module solver
