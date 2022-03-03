module library
  use param
  use biunit

  contains

    subroutine GeneralInfo( )
      external :: fdate
      character(len=24) :: date
      call fdate(date)
      print*, ' '
      print*, '- - - - 2D Convetion-Diffusion-Reaction Simulation - - - - '
      print*, ' '
      print*,' ',date
      print*,'!================= GENERAL INFO ===============!'
      write(*,"(A19,4x,a13,3X,A1)") ' -Element type:           ', ElemType,''
      write(*,"(A19,4X,I6,1X,A10)") ' -Problem dimension:      ', DimPr, '  '
      write(*,"(A19,4X,I6,1X,A10)") ' -Elements:               ', nelem,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' -Nodal points:           ', nnodes, ' '
      write(*,"(A19,4X,I6,1X,A10)") ' -DoF per node:           ', ndofn, '  '
      write(*,"(A19,4X,I6,1X,A10)") ' -Nodes per element:      ', nne, '    '
      write(*,"(A19,4X,I6,1X,A10)") ' -Total Gauss points:     ', totGp,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' -Element variabless:     ', nevab,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' -Total unknowns:         ', ntotv,'   '

    endsubroutine GeneralInfo

    subroutine ReadIntegerFile(UnitNum, FileName, NumRows, NumCols, IntegerArray)

      integer :: i, j, status
      integer, intent(in)            :: UnitNum, NumRows, NumCols
      character(len=*), parameter    :: fileplace = "./"
      character (len=*), intent (in) :: FileName
      integer, dimension (1:NumRows, 1:NumCols), intent (out) :: IntegerArray


      open (unit = UnitNum, file =fileplace//FileName, status='old', action='read' , iostat = status)

      read(UnitNum,*) ((IntegerArray(i,j), j=1,NumCols), i=1,NumRows)
      if (status.ne.0) then
        print *, "Status_Int_File  ", status
      else
        continue
      end if
      close (UnitNum)

    end subroutine ReadIntegerFile

    subroutine SetElementNodes(elm_num, element_nodes, nodeIDmap)

      implicit none

      integer,intent(in)                      :: elm_num ! number of element for each elemental integral in do of K global
      real,dimension(nne,DimPr), intent(out)  :: element_nodes
      integer,dimension(nne), intent(out)     :: nodeIDmap
      integer                                 :: i,j,global_node_id


      element_nodes = 0.0
      nodeIDmap = 0

      do i = 1, nne
        global_node_id = lnods(elm_num,i+1)
        do j=1 ,DimPr
          element_nodes(i,j) = coord(global_node_id,j+1)
        end do
        nodeIDmap(i) = global_node_id
      end do

    end subroutine SetElementNodes

    function J2D( element_nodes, dN_dxi, dN_deta, Gp)
      implicit none

      integer, intent(in)                      :: Gp !esta variable se usara en el lazo principal con el punto de Gauss
      real, dimension(nne,DimPr), intent(in)  :: element_nodes
      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(DimPr,nne)  :: Basis2D
      double precision, dimension(1,nne)      :: Nxi, Neta
      double precision, dimension(DimPr,DimPr) :: J2D

      !con estas instrucciones extraigo la columna de Nx como renglon y lo guardo en Nxi, Gp se
      !ira moviendo conforme la funcion J2D sea llamada en el lazo principal para cada elemento lo mismo para Neta con dN_deta
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)

      !Las siguientes tres lineas realizan de forma implicita el calculo de las derivadas
      !espaciales es decir dN/dx and dN/dy (eq. 5.114 - 5.117). Las derivadas espaciales
      !no se calcula explicitamente, en su lugar se usa:

      !            d/dy = (d/deta)J^-1      (ver eq. 5.77 y 5.109)

      Basis2D(1,:) = Nxi(1,:)
      Basis2D(2,:) = Neta(1,:)
      J2D = matmul(Basis2D,element_nodes)

      return

    end function J2D

    subroutine DerivativesXY(Gp, InvJaco, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)

      implicit none

      double precision, dimension(DimPr,DimPr),intent(in):: InvJaco
      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(3,nne), intent(in)     :: Hesxieta
      integer, intent(in)                                :: Gp !esta variable se usara en el lazo principal
      double precision, dimension(1,nne)                 :: Nxi, Neta
      double precision, dimension(2,nne)                 :: derst
      integer                                            :: idime, inode, jdime
      double precision, dimension(2,nne), intent(out)    :: dN_dxy
      double precision, dimension(3,nne), intent(out)    :: HesXY


      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)


      do idime=1,DimPr
        do inode=1,nne
          dN_dxy(idime,inode) = 0.0
          derst(1,inode) = Nxi(1,inode)
          derst(2,inode) = Neta(1,inode)
          do jdime=1,2
            dN_dxy(idime,inode)  = dN_dxy(idime,inode) + InvJaco(idime,jdime) * derst(jdime,inode)
          end do
        end do
      end do

      !The Hessian matrix
      HesXY = 0.0
      do inode=1,nne
        HesXY(1,inode) = InvJaco(1,1)*InvJaco(1,1)*Hesxieta(1,inode)+&
          2.0*InvJaco(1,1)*InvJaco(2,1)*Hesxieta(2,inode)+InvJaco(2,1)*InvJaco(2,1)*Hesxieta(3,inode)

        HesXY(2,inode) = InvJaco(1,1)*InvJaco(1,2)*Hesxieta(1,inode) + (InvJaco(1,1)*InvJaco(2,2)+&
          InvJaco(2,1)*InvJaco(1,2)) * Hesxieta(2,inode) + InvJaco(2,1)*InvJaco(2,2)*Hesxieta(3,inode)

        HesXY(3,inode) = InvJaco(2,2)*InvJaco(2,2)*Hesxieta(3,inode)+&
          2.0*InvJaco(2,2)*InvJaco(1,2)*Hesxieta(2,inode)+InvJaco(1,2)*InvJaco(1,2)*Hesxieta(1,inode)
      end do


    end subroutine DerivativesXY

    function inv2x2(A)

      implicit none

      double precision, dimension(DimPr, DimPr), intent(in) :: A
      double precision, dimension(DimPr, DimPr)             :: inv2x2
      double precision, dimension(DimPr,DimPr)              :: cofactor
      double precision, parameter :: EPS = 1.0E-10
      double precision            :: det


      det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      if (abs(det) .le. EPS) then
        inv2x2 = 0.0D0
        return
      end if

      cofactor(1,1) = +A(2,2)
      cofactor(1,2) = -A(2,1)
      cofactor(2,1) = -A(1,2)
      cofactor(2,2) = +A(1,1)

      inv2x2 = transpose(cofactor) / det

      return

    end function inv2x2

    subroutine invmtx(a,deter,b)

      ! This routine inverts a square matrix A -> Mat(ndofn,ndofn). The
      ! inverse is stored in B. Its determinant is DETER

      implicit none

      double precision, intent(in) :: a(3,3)
      double precision :: deter, t1, t2, t3, denom
      double precision, intent(out) :: b(3,3)

      !nvers of a 1*1 matrix

      if(ndofn.eq.1) then
        deter=a(1,1)
        if(deter.eq.0.0) return
        b(1,1) = 1.0/a(1,1)
        return
      endif

      !invers of a 2*2 matrix

      if(ndofn.eq.2) then
        deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
        if(deter.eq.0.) return
        denom=1.0/deter
        b(1,1) = a(2,2)*denom
        b(2,2) = a(1,1)*denom
        b(2,1) =-a(2,1)*denom
        b(1,2) =-a(1,2)*denom
        return
      endif

      !inverse of a 3*3 matrix

      if(ndofn.eq.3) then
        t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
        t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
        t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
        deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
        if(deter.eq.0.0) return
        denom = 1./deter
        b(1,1) = t1*denom
        b(2,1) = t2*denom
        b(3,1) = t3*denom
        b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
        b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
        b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
        b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
        b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
        b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom
        return
      endif

    end subroutine invmtx

    subroutine sqrtma(mainp,maout)

      ! Square root of matrix. In the case NDOFN = 3, it is assumed that this
      ! matrix has the form diag(A,a), where A is a 2 x 2 matrix.

      implicit none

      double precision, intent(inout) :: mainp(3,3)
      integer :: i
      double precision, intent(out) :: maout(3,3)

      do i = 1,2
        mainp(i,i) = abs(mainp(i,i))
      end do

      call sqrtm2(mainp,maout)

      if(ndofn.eq.3)then
        maout(3,3) = sqrt(abs(mainp(3,3)))
      endif

    end subroutine sqrtma

    subroutine sqrtm2(mainp,maout)

      ! Square root of a 2 x 2 matrix (from a 3 x 3 matrix)

      implicit none
      double precision, intent(in) :: mainp(3,3)
      double precision :: a, b, c, d, aux1, aux2, vap1, vap2, det, sq1, sq2
      double precision, intent(out) :: maout(3,3)

      a = mainp(1,1)
      b = mainp(1,2)
      c = mainp(2,1)
      d = mainp(2,2)
      aux1 =  0.5*(a+d)
      aux2 = 0.25*(a-d)*(a-d) + b*c
      if(aux2.lt.0.0) then
        !call runend('SQRTMA: Non real eigenvalue in A')
        print*, 'SQRTMA: Non real eigenvalue in A'
        stop
      else if(aux2.lt.1.0e-10) then                         ! b or c = 0, a = d
        maout(1,1) = sqrt(a)
        maout(1,2) = 0.0
        maout(2,1) = 0.0
        maout(2,2) = sqrt(d)
        if(abs(b).gt.1.0e-10) then
          maout(1,2) = b/(sqrt(a) + sqrt(d))
        else if(abs(c).gt.1.0e-10) then
          maout(2,1) = c/(sqrt(a) + sqrt(d))
        end if
      else
        vap1 = aux1 + sqrt(aux2)                            ! vep1 = ( b ,vap1 -a )
        vap2 = aux1 - sqrt(aux2)                            ! vep2 = ( vap2 -d, c )
        if(abs(b)+abs(vap1-a).lt.1.0e-15) then
          sq1  = vap1
          vap1 = vap2
          vap2 = sq1
        end if
        sq1 = sqrt(vap1)
        sq2 = sqrt(vap2)
        vap1 = vap1 - a
        vap2 = vap2 - d
        det = b*c - vap1*vap2
        maout(1,1) = (b*c*sq1-vap1*vap2*sq2)/det
        maout(1,2) =  b*vap2*(-sq1+sq2)/det
        maout(2,1) =  c*vap1*( sq1-sq2)/det
        maout(2,2) =(-vap1*vap2*sq1+b*c*sq2)/det
      end if

    end subroutine sqrtm2

    function m22det(A)

      implicit none
      double precision :: m22det
      double precision, dimension(2,2), intent(in)  :: A



      m22det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

      return

    end function m22det


    subroutine gather(lnods, vecgl, veclo)
      !        gather(vecgl,veclo,lnods,ndofn,nnode)
      !   call gather(coord,elcod,lnods(1,ielem),2,nnode)
      !*****************************************************************************
      !
      !     Gather operations: Recupera los nodos elementales (veclo) del vector global vecgl
      !
      !*****************************************************************************
      !
      !veclo - vector global
      !vecgl - vector global
      
      implicit none
      
      double precision, dimension(*), intent(in) :: vecgl(*)
      integer, intent(in) ::   lnods(nne)
      integer   inode,idofn,ipoin,ievab,itotv
      double precision, dimension(ndofn*nne), intent(out) :: veclo(ndofn*nne)
      
      do inode=1,nne
        ipoin=lnods(inode)
        do idofn=1,ndofn
          ievab=(inode-1)*ndofn+idofn
          itotv=(ipoin-1)*ndofn+idofn
          veclo(ievab)=vecgl(itotv)
        end do
      end do
      
    end subroutine gather








    function elemSize(InvJacobian)
      implicit none

      double precision, dimension(DimPr,DimPr), intent(in) :: InvJacobian
      double precision :: hx, hy, elemSize
     ! hx    = sqrt(xjaci(1,1)**2+xjaci(2,1)**2)
     ! hy    = sqrt(xjaci(1,2)**2+xjaci(2,2)**2)

      hx    = sqrt(InvJacobian(1,1)**2+InvJacobian(2,1)**2)
      hy    = sqrt(InvJacobian(1,2)**2+InvJacobian(2,2)**2)

      elemSize = hnatu/(min(hx,hy))     !hnatu = Reference element length en mod_param

      return

    end function elemsize

    subroutine Galerkin(dvol, basis, dNdxy, Ke, Ce, rhslo)

      implicit none

      double precision, intent(in) :: basis(nne), dNdxy(DimPr,nne)
      double precision, intent(in) :: dvol
      integer :: inode, idofn, ievab, jevab, jnode, jdofn, i, j
      double precision ::  diff, convec, reac, cpcty
      double precision, intent(out) :: Ke(nevab,nevab), rhslo(nevab), Ce(nevab,nevab)
      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              diff=0.0
              do i=1,2
                do j=1,2                      !conductivity tensor
                  diff=diff+ dNdxy(i,inode) * difma(idofn,jdofn,i,j)* dNdxy(j,jnode)
                end do
              end do
              convec=0.0
              do i=1,2
                convec = convec + basis(inode) * conma(idofn,jdofn,i) * dNdxy(i,jnode)
              end do
              reac = basis(inode) * reama(idofn,jdofn) * basis(jnode)
              cpcty = basis(inode) * basis(jnode)
              Ke(ievab,jevab) = Ke(ievab,jevab) + (diff + convec + reac) * dvol
              Ce(ievab,jevab) = Ce(ievab,jevab) + cpcty * dvol                                 !element Capacity (Mass) matrix
            end do
          end do
          rhslo(ievab) = rhslo(ievab) + basis(inode) * force(idofn) * dvol
        end do
      end do

    end subroutine Galerkin

    subroutine pertur( idofn, jdofn, workm, derxy, basis, pertu )
      
      !***************************************************************************
      !
      ! Perturbation of the test function according to the type of method:
      !
      !   SUPG    :                 A_i   V,i           (kstab=1)
      !   GLS     : -(K_ij V,j),i + A_i   V,i + S   V   (kstab=2)
      !   SGS, TG :  (K_ij V,j),i + A^t_i V,i - S^t V   (kstab=3,5)
      !   CG      :            diag(A_i)  V,i           (kstab=4)
      !***************************************************************************
      
      implicit none
      
      
      double precision, intent(in)     :: workm(2,2),derxy(2),basis
      integer                          :: idofn, jdofn, k, l
      double precision                 :: prod1, prod2, prod3
      double precision, intent(in out) :: pertu
      
      !common/proper/difma,conma,reama,force
      
      
      ! SUPG
      if(kstab.eq.1) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        pertu=prod1
        
        ! Galerkin least squares
      else if(kstab.eq.2) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            prod2=prod2+difma(jdofn,idofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(jdofn,idofn)*basis
        pertu=-prod2+prod1+prod3
        
        ! Subgrid scale & Taylor Galerkin
      else if((kstab.eq.3).or.(kstab.eq.5)) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(idofn,jdofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            prod2=prod2+difma(idofn,jdofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(idofn,jdofn)*basis
        pertu=prod2+prod1-prod3
        
        ! Characteristic Galerkin
      else if(kstab.eq.4) then
        prod1=0.0
        if(idofn.eq.jdofn) then
          do k=1,2
            prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
          end do
        end if
        pertu=prod1
      end if
      
    end subroutine pertur
        
    subroutine Stabilization(dvolu, basis, derxy,hesxy,tauma,Ke,rhslo)
      !subroutine Stabilization(dvolu, basis, derxy,hesxy,tauma,Ke,rhslo,pertu,workm,resid)

      ! Contribution to the system matrix and RHS from the stabilization term
      
      implicit none

      double precision, intent(in)  :: basis(nne), derxy(DimPr,nne), hesxy(3,nne), tauma(3,3)
      double precision, intent(in)  :: dvolu
      double precision              :: pertu(nevab,ndofn), workm(2,2),  resid(ndofn,nevab)
      double precision              :: prod1, prod2, prod3
      integer                       :: ievab, inode, idofn, jdofn, jevab, jnode, k, l
      double precision, intent(out) :: Ke(nevab,nevab), rhslo(nevab)

      ! integer :: nnode,ndofn,nevab,kstab,n_ini
      !difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3)
      !common/proper/difma,conma,reama,force

      ievab = 0
      ! n_ini = ndofn*nevab
      ! v_ini = 0.0
      ! call initia(pertu,n_ini,v_ini)
      pertu =  0.0

      do inode=1,nne
        workm(1,1)=hesxy(1,inode)
        workm(2,2)=hesxy(3,inode)
        workm(1,2)=hesxy(2,inode)
        workm(2,1)=workm(1,2)
        do idofn=1,ndofn

          ievab=ievab+1
          do jdofn=1,ndofn
            prod1 = reama(jdofn,idofn)*basis(inode)
            prod2=0.0
            do k=1,2
              prod2 = prod2 + conma(jdofn,idofn,k)*derxy(k,inode)
            end do
            prod3=0.0
            do k=1,2
              do l=1,2
                prod3 = prod3 + difma(jdofn,idofn,k,l)*workm(k,l)
              end do
            end do

            resid(jdofn,ievab) = prod1 + prod2 - prod3
            call pertur( idofn, jdofn, workm, derxy(1,inode), basis(inode), pertu(ievab,jdofn) )
          end do
        end do
      end do

      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              prod1=0.0
              do k=1,ndofn
                do l=1,ndofn
                  prod1 = prod1 + pertu(ievab,k)*tauma(k,l)*resid(l,jevab)
                end do
              end do
              Ke(ievab,jevab) = Ke(ievab,jevab) + prod1 * dvolu
            end do
          end do

          prod1=0.0
          do k=1,ndofn
            do l=1,ndofn
              prod1 = prod1 + pertu(ievab,k) * tauma(k,l) * force(l)
            end do
          end do
          rhslo(ievab) = rhslo(ievab) + prod1 * dvolu
        end do
      end do

    end subroutine Stabilization

    subroutine TauMat(hmaxi,tauma)
      !
      !     Matrix of intrinsic time scales, computed as
      !
      !     TAU = PATAU * [ 4 K / h^2 + 2 A / h + S ]^{-1}

      implicit none

      double precision, intent(in) :: hmaxi
      integer :: i, j, k
      !double precision :: difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3) !Declaradas en parameters
      double precision :: chadi(3,3), chaco(3,3), chare(3,3), tauin(3,3)
      double precision :: a, b, c, tau, det            !hnatu -> declarado en parameters
      double precision, intent(out) :: tauma(3,3)               !ndofn -> en parameters

      !common/numert/hnatu,patau,ksoty,kprec,kstab,ktaum
      !common/proper/difma,conma,reama,force

      !v_ini = 0.0
      !call initia(tauma,9,v_ini)
      tauma = 0.0
      if(kstab.eq.0) return
      !call initia(tauin,9,v_ini)
      !call initia(chaco,9,v_ini)
      !call initia(chadi,9,v_ini)
      tauin = 0.0
      chaco = 0.0
      chadi = 0.0

      !  Characteristic convection matrix: A = sqrt | A_i A_i |
      do i=1,ndofn
        do j=1,ndofn
          chaco(i,j)=0.0
          do k=1,ndofn
            chaco(i,j) = chaco(i,j) + conma(i,k,1)*conma(k,j,1) + conma(i,k,2)*conma(k,j,2)
          end do
        end do
      end do
      call sqrtma(chaco,chaco)

      !  Characteristic diffusion matrix: K = sqrt( K_ij K_ij )
      do i=1,ndofn
        do j=1,ndofn
          chadi(i,j)=0.0
          do k=1,ndofn
            chadi(i,j) = chadi(i,j) + difma(i,k,1,1) * difma(k,j,1,1) + &
              &difma(i,k,1,2)*difma(k,j,1,2)*2.0 + difma(i,k,2,2)*difma(k,j,2,2)
          end do
        end do
      end do

      call sqrtma(chadi,chadi)

      !  Characteristic reaction matrix: S = | S |
      do i=1,ndofn
        do j=1,ndofn
          chare(i,j)=0.0
          do k=1,ndofn
            chare(i,j)=chare(i,j) + reama(i,k) * reama(k,j)
          end do
        end do
      end do
      call sqrtma(chare,chare)

      ! Invers of the matrix of characteristic times
      do i=1,ndofn
        do j=1,ndofn
          tauin(i,j) = 4.0*chadi(i,j)/(hmaxi*hmaxi) + 2.0*chaco(i,j)/hmaxi + chare(i,j)
        end do
      end do

      !  Matrix tau, corresponding to:
      !     KTAUM = 0: T = t I, where t is the minimum of all the admissible tau's
      !           = 1: T = diag(t1,t2,t3), where ti is the minimum of the admissible tau's for the i-th row (equation)
      !           = 2: T = [ 4 K / h^2 + 2 A / h + S ]^{-1}

      if(ktaum.eq.0) then
        tau = 0.0
        do i=1,ndofn
          do j=1,ndofn
            tau = max(tau,abs(tauin(i,j)))
          end do
        end do
        tau = patau/tau
        do i=1,ndofn
          tauma(i,i) = tau
        end do

      else if(ktaum.eq.1) then
        a = 0.0
        b = 0.0
        c = 0.0
        do j=1,ndofn
          a = max(a,abs(tauin(    1,j)))
          b = max(b,abs(tauin(    2,j)))
          c = max(c,abs(tauin(ndofn,j)))
        end do
        a = patau/a
        b = patau/b
        c = patau/c
        tauma(    1,    1) = a
        tauma(    2,    2) = b
        tauma(ndofn,ndofn) = c

      else if(ktaum.eq.2) then
        call invmtx(tauin,det,tauma)
        do i = 1,ndofn
          do j = 1,ndofn
            tauma(i,j) = tauma(i,j)*patau
          end do
        end do
        tauma(ndofn,ndofn) = 0.0

      else if(ktaum.eq.3) then
        a = 1.0/(patau*difma(1,1,1,1) /(hmaxi*hmaxi) + reama(1,1))
        tauma(1,1) = a
        tauma(2,2) = a
        a = (hmaxi*hmaxi*hmaxi*hmaxi)/(patau*patau)
        a = a*(patau/(hmaxi*hmaxi*reama(1,1)) + 1.0d0/(difma(1,1,1,1)))
        tauma(3,3) = a
      end if

    end subroutine TauMat

    subroutine VinculBVs(  BVs, nofix, ifpre, presc )

      implicit none

      !integer, intent(in)              :: nBvs, nBVscol ya no se ponen estan en el modulo parameter y se comunica el valor
      integer, intent(in) :: BVs( nBvs, nBVscol)
      integer             :: i, j
      double precision, intent(out) :: presc(ndofn,nBVs)
      integer, intent(out)          :: ifpre(ndofn,nBVs)
      integer, intent(out)          :: nofix(nBVs)


      select case(ndofn)
        case(1)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = BVs(j,1)
              ifpre(i,j) = BVs(j,2) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,3)
            end do
          end do

        case(2)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = BVs(j,1)
              ifpre(i,j) = BVs(j,i+1) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,i+3)
            end do
          end do

        case(3)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = BVs(j,1)
              ifpre(i,j) = BVs(j,i+1) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,i+4)
            end do
          end do

        case DEFAULT
          write(*,*) 'Exceeded DoF'
        end select
    end subroutine VinculBVs

    subroutine BandWidth( )
      !use stdlib_linalg, only: diag

      implicit none

      integer :: iband, ielem, inode, ipoin, jnode, jpoin, nband       ! , C, D
      !integer :: i,j,k
      !real, allocatable, dimension(:) :: A, B, BB
      !real, allocatable :: AA(:,:)

      iband=0
      do ielem =1, nelem
        do inode = 1, nne
          ipoin = lnods(ielem,inode+1) !Este +1 es para que comience en los nodos (columna 2) y no del numeor de elemento
          do jnode = 1, nne
            jpoin = lnods(ielem,jnode+1)
            iband = max(iband,abs(jpoin-ipoin))
          end do
        end do
        !if (iband.gt.nband) nband=iband !nband pasa como variable global por que se usa en  ApplyBVal y otros
        nband = iband
      end do
      nband=(nband+1)*ndofn-1
      if(nband.ge.maxband) then      !Puedo poner a maxband como variable local (solo se usa aqui) si lo hago debe ir como dummy var
        !                            en subroutine globaK y en bandwidth. Tambien ldakban puede pasar a ldA y ponerla como out pues
        !                            se usa en global assemb y aout no vale la pena ponerla como global.
        write(*,'(a,i5,a)') ' >>> Hay que aumentar MAXBAND a ',nband+1,' !!!'
        stop
      end if
      upban  = nband
      lowban = nband
      totban = lowban + upban + 1
      ldAKban= 2*lowban + upban + 1

      !AA = diag([(0.5*i+1*0.851,i=1,10)]) ! creates a 10 by 10 identity matrix

      !do k =1,10
      !  write(*,"(10(1x,f10.2))") (AA(k,j),j=1,10)
      !end do
      !allocate(BB(size(AA,1)))
      !allocate(A(size(conma,1)))
      !allocate(B(size(reama,1)))

      !BB = diag(AA)
      !A = diag(conma(ndofn,ndofn,1))
      !B = diag(reama)
      !select case(ndofn)
      !
      !case(1)
      !  C = abs(sum(conma))
      !  D = abs(sum(reama))
      !  if((C.EQ.0) .OR. (D.EQ.0))then
      !    write(*,*) '-Global Matrix is symmetric'
      !    upban  = nband
      !    lowban = nband
      !    totban = lowban + upban + 1
      !    ldAKban= 2*lowban + upban + 1
      !  else
      !    write(*,*) '-Global Matrix is structural-or-non symmetric'
      !    upban  = nband
      !    lowban = nband
      !    totban = lowban + upban + 1
      !    ldAKban= 2*lowban + upban + 1
      !  endif
      !
      !case(2)
      !  C = abs(sum(conma))
      !  D = abs(sum(reama))
      !  if((C.EQ.0) .OR. (D.EQ.0))then
      !    write(*,*) '-Global Matrix is symmetric'
      !    upban  = nband
      !    lowban = 0
      !    totban = lowban + upban + 1
      !    ldAKban= 2*lowban + upban + 1
      !  else
      !    write(*,*) '-Global Matrix is structural-or-non symmetric'
      !    upban  = nband
      !    lowban = nband
      !    totban = lowban + upban + 1
      !    ldAKban= 2*lowban + upban + 1
      !  endif
      !
      !case(3)
      !  C = abs(sum(conma) )
      !  D = abs(sum(reama) )
      !  if((C.EQ.0) .OR. (D.EQ.0))then
      !    write(*,*) '-Global Matrix is symmetric'
      !    upban  = nband
      !    lowban = 0
      !    totban = lowban + upban + 1
      !    ldAKban= 2*lowban + upban + 1
      !  else
      !    write(*,*) '-Global Matrix is structural-or-non symmetric'
      !    upban  = nband
      !    lowban = nband
      !    totban = lowban + upban + 1
      !    ldAKban= 2*lowban + upban + 1
      !  endif
      !
      !end select


      write(*,*) ''
      print*,'!================ Bandwidth Info ==============!'

      write(*,"(A15,9X,I6,1X,A9)")'-UpBand:      ', upban,'   '
      write(*,"(A15,9X,I6,1X,A9)")'-LowBand:     ', lowban,'  '
      write(*,"(A15,9X,I6,1X,A9)")'-TotBand:     ', totban,'  '
      write(*,"(A15,9X,I6,1X,A9)")'-ledimAK:     ', ldAKban,' '


    end subroutine BandWidth

    subroutine GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_K, A_C, A_F)

      implicit none

      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(3,nne), intent(in)     :: Hesxieta
      double precision, dimension(nne)          :: basis
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce
      double precision, dimension(nevab)        :: rhslo
      double precision, dimension(3,3)          :: tauma
      real, dimension(nne,DimPr)                :: element_nodes
      integer, dimension( nne + 1, nelem)       :: lnods2
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_K, A_C, A_F
      
      call BandWidth( )
      allocate( A_K(ldAKban,ntotv), A_C(ldAKban,ntotv), A_F(ntotv, 1) )
      
      !duda rhslo se declara como a(n) y en la rutina assembleF como a(n,1), pero compila y ejecuta bien. ¿Poooor?
      A_K = 0.0
      A_F = 0.0
      do ielem = 1, nelem 
        !gather
        Ke = 0.0       !Esto es amate
        rhslo = 0.0    !rhslo(nevab)
        call SetElementNodes(ielem, element_nodes, nodeIDmap)
        !do-loop: compute element stiffness matrix Ke
        do igaus = 1, TotGp
          Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
          detJ = m22det(Jaco)
          Jinv = inv2x2(Jaco)
          dvol = detJ *  weigp(igaus,1)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          call Galerkin(dvol, basis, dN_dxy, Ke, Ce, rhslo) !amate lo llame Ke
          call TauMat(hmaxi,tauma)
          !!call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, rhslo, pertu,workm,resid)
          call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, rhslo)
        end do
        lnods2=transpose(lnods)
        call Assemble_K(nodeIDmap, Ke, A_K)     !Assemble Global Conductivity Matrix K
        call Assemble_K(nodeIDmap, Ce, A_C)     !Assemble Global Capacity Matrix C          cambiar nombre por AssemGlobalMat
        call AssembleF(nodeIDmap, rhslo, A_F)   !Assemble Global Source vector F
      end do

      
      !print*, 'shape of tauma', shape(tauma)
      !print*, ' '
      !do i = 1,3
      !  print'(4F10.3)', (tauma(i,j), j=1,3)
      !end do
      !
      !print*, 'shape of HesXY', shape(HesXY)
      !print*, ' '
      !do i = 1,3
      !  print'(4F10.3)', (HesXY(i,j), j=1,nne)
      !end do
      !
      !print*, 'shape of Hesxieta', shape(Hesxieta)
      !print*, ' '
      !do i = 1,3
      !  print'(4F10.3)', (Hesxieta(i,j), j=1,nne)
      !end do

    end subroutine GlobalSystem

    subroutine Assemble_K(lnods,Ke,A_K)
      !subroutine Assemble_K(ielem,lnods,Ke,A_K)
      !*****************************************************************************
      !
      !    Fa l'assembly de les matrius de CDR de cada elemento en la matriu global
      !
      !*****************************************************************************

      implicit none
      !common /contrl/ npoin,nelem,nmats,nvfix,nload,nband,ntotv
      double precision, intent(in) :: Ke(nevab,nevab)
      integer, intent(in) :: lnods(nne)
      integer :: inode, ipoin, idofn, ievab, itotv, jnode, jpoin, jdofn, jevab, jtotv, jband !, i,j,k,l
      double precision, intent(in out) :: A_K(ldAKban,ntotv)
     !                                        ldAKban= 2*lowban + upban + 1
     !                                          totban = lowban + upban + 1
      !if(lowban.EQ.0)then
      !  write(*,*) 'Symmetric case'
      !  do inode=1,nne    !nne = number of node in the element
      !    ipoin=lnods(inode)
      !    !print*,'ipoin',ipoin
      !    do idofn=1,ndofn
      !      ievab=(inode-1)*ndofn+idofn
      !      itotv=(ipoin-1)*ndofn+idofn
      !      do jnode=1,nne
      !        jpoin=lnods(jnode)
      !        do jdofn=1,ndofn
      !          jevab=(jnode-1)*ndofn+jdofn
      !          jtotv=(jpoin-1)*ndofn+jdofn

      !          jband=jtotv-itotv+1   !   ------> Original de retpla
      !          if (jband.ge.1)then
      !          A_K(jband,itotv)=A_K(jband,itotv) +Ke(ievab,jevab)
      !          endif
      !        end do
      !      end do
      !    end do
      !  end do
      !  !Reorder the global band matrix into LAPACK format
      !  !i=0
      !  !do k = nband+1, 1,-1
      !  !  i = i+1
      !  !  j = 1
      !  !  do l = ntotv,1,-1
      !  !    A_K(i,j) = A_K(k,l)
      !  !    j=j+1
      !  !  end do
      !  !end do
      !else
      !  !Non or structural symmetric case
      !  do inode=1,nne    !nne = number of node in the element
      !    ipoin=lnods(inode)
      !    !print*,'ipoin',ipoin
      !    do idofn=1,ndofn
      !      ievab=(inode-1)*ndofn+idofn
      !      itotv=(ipoin-1)*ndofn+idofn
      !      do jnode=1,nne
      !        jpoin=lnods(jnode)
      !        do jdofn=1,ndofn
      !          jevab=(jnode-1)*ndofn+jdofn
      !          jtotv=(jpoin-1)*ndofn+jdofn

      !          !jband=jtotv-itotv+1      ------> Original de retpla
      !          !A_K(jband,itotv)=A_K(jband,itotv) + ...

      !          !jband= upban +1 +itotv-jtotv    !Algoritmo de recuperacion para LAPACK
      !          jband = itotv-jtotv + totban
      !          if (jband.ge.1)then
      !            A_K(jband,jtotv)=A_K(jband,jtotv)+Ke(ievab,jevab)
      !          endif
      !        end do
      !      end do
      !    end do
      !  end do
      !endif

      do inode=1,nne    !nne = number of node in the element
        ipoin=lnods(inode)
        !print*,'ipoin',ipoin
        do idofn=1,ndofn
          ievab=(inode-1)*ndofn+idofn
          itotv=(ipoin-1)*ndofn+idofn
          do jnode=1,nne
            jpoin=lnods(jnode)
            do jdofn=1,ndofn
              jevab=(jnode-1)*ndofn+jdofn
              jtotv=(jpoin-1)*ndofn+jdofn

              !jband=jtotv-itotv+1      ------> Original de retpla
              !A_K(jband,itotv)=A_K(jband,itotv) + ...

              !jband= upban +1 +itotv-jtotv    !Algoritmo de recuperacion para LAPACK
              jband = itotv-jtotv + totban
              if (jband.ge.1)then
                A_K(jband,jtotv)=A_K(jband,jtotv)+Ke(ievab,jevab)
              endif
            end do
          end do
        end do
      end do

      return
    end subroutine Assemble_K

    subroutine AssembleF( nodeIDmap, fe, F_global)

      implicit none

      double precision, dimension(nevab,1), intent(in)   :: fe  !nevab = Num of element variable
      integer, dimension(nne), intent(in)                :: nodeIDmap
      integer :: i, rowNode, row
      double precision, dimension(ntotv,1),intent(inout):: F_global

      do i = 1, nne
        rowNode = nodeIDmap(i)                !global id for row nodes
        row =  ndofn * rowNode - (ndofn-1)    !row number in the global F
        F_global(row:row+ndofn-1,1) =  F_global(row:row+ndofn-1,1) + fe( (i-1)*ndofn+1:i*ndofn, 1)

      end do

    end subroutine AssembleF

    subroutine ApplyBVs(nofix,ifpre,presc,rigid,gload)
      !                                   ,A_K ,A_F
      !        vincul(rigid,gload,treac,nofix,ifpre,presc)
      !*****************************************************************************
      !
      !   Imposa les condicions de contorn
      !
      !*****************************************************************************

      implicit none !real*8(a-h,o-z)

      !Agregar un common a BVs para guardar nBVs y nBVscol asi como nband.
      !Ya se hizo y se uso el modulo mod_param para guardarlos ahi y el valor se comparte
      double precision,intent(in) :: presc(ndofn,nBVs)
      integer, intent(in)         :: nofix(nBVs), ifpre(ndofn,nBVs)
      !                                    nvfix             ,nvfix
      !common /contrl/ npoin,nelem,nmats,nvfix,nload,nband,ntotv
      integer :: ivfix, idofn, itotv, jpoin, jdofn, jtotv, itot1, jband, itot2, nvfix
      double precision,intent(inout)  :: rigid(ldAKban,ntotv), gload(ntotv)

      nvfix = nBVs
      !***  Inicialitzacio de les reaccions per als graus de llibertat prescrits
      do ivfix=1,nvfix
        do idofn=1,ndofn!3
          if (ifpre(idofn,ivfix).eq.1) then
            itotv              = (nofix(ivfix)-1)*ndofn+idofn   !3+idofn
            !treac(idofn,ivfix) = -gload(itotv) !Yo no uso treac, solo gload. Preguntar a Ramon como definir gload
          end if
        end do
      end do

      !***  Llac sobre els nodes coaccionats
      do ivfix=1,nvfix
        jpoin=nofix(ivfix)
        do jdofn=1,ndofn !3
          if (ifpre(jdofn,ivfix).eq.1) then
            
            !***  Ca de grau de llibertat prescrit
            jtotv=(jpoin-1)*ndofn+jdofn
            
            !***  Modificacio de les equacions anteriors a la del g.d.ll. prescrit (arriba de diagonal)
            if (jtotv.gt.1) then
              itot1=jtotv-nband
              if (itot1.lt.1) itot1=1
              do itotv=itot1,jtotv-1
                jband=itotv-jtotv+totban               !Algoritmo de recuperacion
                gload(itotv) = gload(itotv)-rigid(jband,jtotv)*presc(jdofn,ivfix)
                rigid(jband,jtotv)=0.0
                
              end do
            end if
            
            !***  Modificacio de les equacions posteriors a la del g.d.ll. prescrit (abajo de la diagonal)
            if (jtotv.lt.ntotv) then
              itot2=jtotv+nband
              if (itot2.gt.ntotv) itot2=ntotv
              do itotv=jtotv+1,itot2
                jband=jtotv-itotv+totban                !Algoritmo de recuperacion
                gload(itotv) = gload(itotv)-rigid(jband,itotv)*presc(jdofn,ivfix)
                rigid(jband,itotv)=0.0
              end do
            end if
            
            !***  Equacio trivial per al grau de llibertat prescrit (en la diagonal)
            rigid(totban,jtotv)=1.0
            gload(jtotv)=presc(jdofn,ivfix)
          end if
        end do
      end do
      return
      
    end subroutine ApplyBVs
    
    subroutine MKLfactoResult( routine_name, num )
      implicit none
      
      character(*), intent(in)  :: routine_name
      integer,      intent(in)  :: num
      character(len=34) :: text
      character(len=44) :: text2
      integer           :: val
      external          :: xerbla
      
      text  = '   *FACTORIZATION DONE WITH STATUS'
      text2 = '   *FACTORIZATION HAS BEEN COMPLETED, BUT U('
      if ( num .eq. 0 ) then
        !print*, ' '
        write(*, 101) text, num, ', THE EXECUTION IS SUCCESSFUL.'
      elseif(num .lt. 0 )then
        val = abs(num)
        print*, ' '
        write(*, 102) '    THE',val,'-TH PARAMETER HAD AN ILLEGAL VALUE.'
        call xerbla( routine_name, num )
      elseif(num .gt. 0 )then
        print*, ' '
        write(*, 103) text2,num,',',num,') = 0.0. Then U IS EXACTLY SINGULAR.'
        print"(A)",'   DIVISION BY 0 WILL OCCUR IF USE THE FACTOR U FOR SOLVING A SYSTEM OF LINEAR EQUATIONS.'
        print*, ' '
        call xerbla( routine_name, num )
        print*, ' ~ ~ ~ Stopping the execution'
        print*, ' '
        !stop
      endif
      print*, ' '

      101 format (A, 1x, I1, A)
      102 format (A, I4, A)
      103 format (A, I3, A, I3, A)
    end subroutine MKLfactoResult

    subroutine MKLsolverResult(routine_name, num )
      implicit none
      character(*), intent(in)  :: routine_name
      integer,      intent(in)  :: num
      integer :: val
      character(len=30) :: text
      character(len=35) :: text2
      character(len=30) :: text3
      external :: xerbla
      text =  '   *SYSTEM SOLVED WITH STATUS'
      text2 = '-TH PARAMETER HAD AN ILLEGAL VALUE.'
      text3 = '   *THE LEADING MINOR OF ORDER'
      if ( num .eq. 0 ) then
        write(*,101) text, num, ', THE EXECUTION IS SUCCESSFUL.'
      elseif(num .lt. 0 )then
        val = abs(num)
        write(*,102) '    THE',val, text2
        call xerbla( routine_name, num )
      elseif(num .gt. 0 )then
        print*, ' '
        write(*, 103) text3,num,' (THEREFORE THE MATRIX A ITSELF) IS NOT'
        print*,'   POSITIVE-DEFINITE. THE FACTORIZATION COULD NOT BE COMPLETED. '
        print*,'   THIS MAY INDICATE AN ERROR IN FORMING THE MATRIX A.'
        print*, ' '
        call xerbla( routine_name, num )
        print*, ' ~ ~ ~ Stopping the execution'
        print*, ' '
        stop
      endif
      101 format (A, 1x, I1, A)
      102 format (A, I3, A)
      103 format (A30, I3, A)
      print*,' '
    end subroutine MKLsolverResult

    subroutine writeMatrix(Matrix, unit1, name1, Vector, unit2, name2)
      implicit none

      character(len=*), parameter    :: fileplace = "Res/"
      character(*) :: name1, name2
      integer :: i, j, mrow, ncol, unit1, unit2
      double precision, dimension(ldAKban ,ntotv ), intent(in) :: Matrix
      double precision, dimension(ntotv ,1), intent(in) :: Vector

      100 format (900E15.5)

      mrow = size(Matrix,1)
      ncol = size(Matrix,2)
      open(unit=unit1, file= fileplace//name1, ACTION="write", STATUS="replace")

      do i=1,mrow
        write(unit1, 100)( Matrix(i,j) ,j=1,ncol)
      end do
      close(unit1)

      open(unit=unit2, file= fileplace//name2, ACTION="write", STATUS="replace")
      do i=1,ncol
        write(unit2, 100) Vector(i,1)
      end do
      close(unit2)
      write(*,*) 'files: ', name1,' and ', name2, ' written succesfully on Res/'

    end subroutine writeMatrix

    subroutine PosProcess(solution, nameFile1, activity)

      implicit none

      character(len=*), parameter    :: fileplace = "Pos/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      character(*), intent(in)                :: nameFile1, activity
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcor, ycor
      integer                                 :: ipoin, ii

      solution_T = transpose(solution)
      xcor  = spread(coord(:,2),dim = 1, ncopies= 1)
      ycor  = spread(coord(:,3),dim = 1, ncopies= 1)

      open(unit=555, file= fileplace//nameFile1, ACTION="write", STATUS="replace")

      if(activity == "msh")then !quitar este if y acomodar el numero de unidad

        write(555,902) 'MESH', '"Domain"', 'dimension', DimPr, 'ElemType', ElemType, 'Nnode', nne
        write(555,"(A)") '#2D Convection-Diffusion-Reaction'
        write(555,900) '#Element tipe: ', ElemType,'/',ElemType
        write(555,"(A)")'Coordinates'
        write(555,"(A)") '#   No        X           Y'
        do ipoin = 1, nnodes
          write(555,906) ipoin, xcor(1,ipoin), ycor(1,ipoin)
        end do
        write(555,"(A)") 'End Coordinates'
        write(555,"(A)") 'Elements'
        do ipoin = 1, nelem
          write(555,908) lnods(ipoin,:)
        end do
        write(555,"(A)") 'End Elements'
        close(555)
        print"(A6,A19,A30)", ' File ',File_PostMsh,'written succesfully in Pos/ '

      elseif(activity == "res")then
        write(555,"(A)") 'GiD Post Results File 1.0'
        write(555,"(A)") '#2D Convection-Diffusion-Reaction'

        ! se escribe el res de las componentes de la velocidad
        select case(ndofn)
          case(1)
            write(555,"(A)") 'Result "DoF" "Concentration" 0 Scalar OnNodes'
            write(555,"(A)") 'ComponentNames "" '
            write(555,"(A)") 'Values'
            write(555,*) '#',   'No    ','             ux '
            !  se escribe el res para el caso escalar de un grado de libertad
            write(555,914)
            do ipoin = 1, nnodes
              write(555,916) ipoin, solution(ipoin, 1)
            end do
            !An alternative work around is to explicitly designate the elements to be read using an io-implied-do.
            !Something like
            !read (unit=10, fmt=*, iostat=iostat) (mat(pcnt,i),i=1,m)
            write(555,"(A)") 'End Values'
          case(2)
            write(555,"(A)") 'Result "DoF" "Concentration" 0 Vector OnNodes'
            write(555,"(A)") 'ComponentNames "u" "v" "--" "" '
            write(555,"(A)") 'Values'
            write(555,*) '#',   'No    ','             ux ','               uy '
            do ipoin = 1, nnodes
              write(555,918) ipoin, solution_T(1, ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
            end do
            write(555,"(A)") 'End Values'
          case(3)
            write(555,"(A)") 'Result "DoF" "Concentration" 0 Vector OnNodes'
            write(555,"(A)") 'ComponentNames "u" "v" "w" "" '
            write(555,"(A)") 'Values'
            write(555,*) '#',   'No    ','             ux ','               uy'
           ! do ipoin = 1, nnodes
           !   write(555,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
           ! end do
            do ipoin = 1, nnodes
              write(555,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1)
            end do
            write(555,"(A)") 'End Values'
        end select
        write(555,"(A)") 'Result "P" "Preassure" 0 Scalar OnNodes'
        write(555,"(A)") 'ComponentNames "" '
        write(555,"(A)") 'Values'
        write(555,*) '#',   'No    ','             P '
        !  se escribe el res para el caso escalar de un grado de libertad
        write(555,914)
        ii=1
        do ipoin = 3, nnodes*3,3
          write(555,916) ii, solution_T(1,ipoin)
          ii=ii+1
        end do
        write(555,"(A)") 'End Values'
        print"(A6,A19,A30)", ' File ',File_PostRes, 'written succesfully in Pos/ '
        
        close(555)
      else
        write(*,"(A)") ' < < Error > > Postprocess activity must be "msh" or "res" '
        close(555)
        stop
      end if


      900 format(A15, A13, A1, A13)
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      906 format(I7,2(3x,f9.4)) !format for msh
      908 format(9(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,E12.5)  !format for scalar case
      918 format(I7,3x,f15.5,3x,f15.5) !format for res velocity
      919 format(I7,3(3x,E15.5)) !format for res velocity

    end subroutine PosProcess



  !Fin de contains


end module library
