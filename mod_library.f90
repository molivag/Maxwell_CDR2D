module library
  use param
  use biunit

  contains
    
    subroutine GeneralInfo( )
      external :: fdate
      character(len=24) :: date
      character(len=4) :: aaaa 
      
      integer :: i,j, k, l
     
      
      if(kstab.eq.3 .or. kstab.eq.5)then
        aaaa = 'SGS'
      elseif(kstab.eq.2)then
        aaaa = 'GLS'
      elseif(kstab.eq.4)then
        aaaa = 'CG'
      elseif(kstab.eq.1)then
        aaaa = 'SUPG'
      elseif(kstab.eq.0)then
        aaaa = 'NONE'
      else
        write(*,'(A)') '> > >Error in stabilization method'
      endif
      
      
      call fdate(date)
      print*, ' '
      print*, '- - - - 2D Convetion-Diffusion-Reaction Simulation - - - - '
      print*, ' '
      print*,' ',date
      print*,'!================= GENERAL INFO ===============!'
      write(*,"(A19,4x,a13,3X,A1)") ' - Element type:           ', ElemType,''
      write(*,"(A19,7x,a5,3X,A1)")  ' - Problem Type:           ', ProbType,''
      write(*,"(A19,4X,I6,1X,A10)") ' - Problem dimension:      ', DimPr, '  '
      write(*,"(A19,4X,I6,1X,A10)") ' - Elements:               ', nelem,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' - Nodal points:           ', nnodes, ' '
      write(*,"(A19,4X,I6,1X,A10)") ' - DoF per node:           ', ndofn, '  '
      write(*,"(A19,4X,I6,1X,A10)") ' - Nodes per element:      ', nne, '    '
      write(*,"(A19,4X,I6,1X,A10)") ' - Total Gauss points:     ', totGp,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' - Element variabless:     ', nevab,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' - Total unknowns:         ', ntotv,'   '
      
      print*, ' '
      print*,'!========== STABILIZATION PARAMETERS ==========!'
      write(*,"(A26,3x,a4,3X,A1)") ' - Stabilization method:   ', aaaa,''
      write(*,"(A26,3x,I2,3X,A1)")  ' - Type of Tau matrix:    ', ktaum,''
      write(*,"(A26,3X,f3.1,1X,A10)") ' - Param. to obtain TAU:   ', patau, '  '
      write(*,"(A26,3X,f3.1,1X,A10)") ' - Lenght ref. element:    ', hnatu,'   '
      write(*,"(A26,2X,f5.2,1X,A10)") ' - Algorithmic constant:   ', Cu, ' '
      write(*,"(A26,1X,e13.5,1X,A10)") ' - Magnetic Permeability:  ', mu, '  '
      write(*,"(A26,1X,f5.1,1X,A10)") ' - Constante of lenght:   ', ell, '    '
      write(*,"(A26,3X,f3.1,1X,A10)") ' - Exponent of mesh size:    ', i_exp,'   '
      write(*,"(A26,1X,f8.4,1X,A10)") ' - Mesh size 2^-i:        ', 10*2**(-i_exp),'   '
      write(*,"(A26,1X,e13.5,1X,A10)") ' - Stab. param.1 (Cu):    ', Cu*mu*(helem**2/ell**2),'   '
      write(*,"(A26,2X,e14.5,2X,A10)") ' - Stab. param.2 (ℓ):       ', ell**2 / mu,'   '
      
      print*, ' '
      print*,'!============ TENSOR COEFFICIENTS  ============!'
      print*, 'Diffusion'
      do i = 1,dimPr
        do j = 1,DimPr
          print"(A,2I1)", 'k_',i,j
          do k = 1,ndofn
            print"(e15.5,1x,e15.5, 1x, e15.5)",( difma(k,l,i,j), l=1,ndofn)
          end do
          !print*,' '
        end do
      end do
      print*, ' '  
      print*, 'Convection'
      do k = 1, DimPr
        print"(A,2I1)",'A_',k
        do i = 1, ndofn
          write(*, "(f10.3, 1x, f10.3, 1x, f10.3)")( conma(i,j,k) ,j=1, ndofn)
        end do
        print*,' '
      end do
      print*,'Reaction'
      do i=1,ndofn
        write(*,"(f10.3, 1x, f10.3, 1x, f10.3)" )( reama(i,j) ,j=1,ndofn)
      end do
        print*,' '
      print*,'External Forces'
      do i =1, ndofn
        print"(f10.3)", force(i)
      end do
      print*, ' '
      
      
    endsubroutine GeneralInfo
    
    subroutine ReadFile(FileName, NumRows, NumCols, IntegerArray)
      
      integer :: i, j, status
      integer, intent(in)            :: NumRows, NumCols
      character(len=*), parameter    :: fileplace = "./"
      character (len=*), intent (in) :: FileName
      integer, dimension (1:NumRows, 1:NumCols), intent (out) :: IntegerArray
      
      
      open (unit=99, file=fileplace//FileName, status='old', action='read' , iostat = status)
      
      read(99,*) ((IntegerArray(i,j), j=1,NumCols), i=1,NumRows)
      if (status.ne.0) then
        print *, "Status while reading BVs file is:", status
      else
        continue
      end if
      close (99)
      
    end subroutine ReadFile
    
    subroutine SetElementNodes(elm_num, element_nodes, nodeIDmap)
      
      implicit none
      
      ! number of element for each elem integral in loop of K global
      integer,intent(in)                                   :: elm_num 
      integer                                              :: i,j,global_node_id
      double precision, dimension(nne,DimPr), intent(out)  :: element_nodes
      integer,dimension(nne), intent(out)                  :: nodeIDmap
      
      
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
      
      integer, intent(in)      :: Gp ! se usara en el lazo principal con el punto de Gauss
      double precision, dimension(nne,DimPr), intent(in) :: element_nodes
      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(DimPr,nne)  :: Basis2D!, derst
      !double precision, dimension(DimPr,nne)  :: elcod
      double precision, dimension(1,nne)      :: Nxi, Neta
      double precision, dimension(DimPr,DimPr):: J2D!, xjacm
      !integer                                 :: ideime, jdime 
      !con estas instrucciones extraigo la columna de Nx como renglon y lo guardo en Nxi, Gp se
      !ira moviendo conforme la funcion J2D sea llamada en el lazo principal para 
      !cada elemento lo mismo para Neta con dN_deta
      
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)
      
      !Las siguientes tres lineas realizan de forma implicita el calculo de las derivadas
      !espaciales es decir dN/dx and dN/dy (eq. 5.114 - 5.117). Las derivadas espaciales
      !no se calcula explicitamente, en su lugar se usa:
      
      !            d/dy = (d/deta)J^-1      (ver eq. 5.77 y 5.109)
      
      Basis2D(1,:) = Nxi(1,:)
      Basis2D(2,:) = Neta(1,:)
      
      J2D = matmul(Basis2D,element_nodes)
      
     
      !do inode = 1,nne
      !  derst(1,inode) = dN_dxi(i,Gp) 
      !  derst(2,inode) = dN_deta(i,Gp) 
      !end do
      !
      !do jdime = 1,2
      !  do inode = 1,nne
      !   elcod(jdime,inode) = element_nodes(inode, jdime)
      !  end do
      !end do
     
      !
      !do idime=1,2
      !  do jdime=1,2
      !    xjacm(idime,jdime)=0.0
      !    do inode=1,nnode
      !      xjacm(idime,jdime) = xjacm(idime,jdime) + derst(idime,inode)*elcod(jdime,inode)
      !    end do
      !  end do
      !end do
      
      return
      
    end function J2D
    
    function m22det(A)
      
      implicit none
      
      double precision :: m22det
      double precision, dimension(2,2), intent(in)  :: A
      
      m22det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
      
      !djacb  = xjacm(1,1)*xjacm(2,2) - xjacm(1,2)*xjacm(2,1)
      
      if(m22det.lt.1.e-8)then
        write(*,*)'Element with non-positive Jacobian'
      end if
     
      return
      
    end function m22det
    
    function inv2x2(A)
      
      implicit none
      
      double precision, parameter :: EPS = 1.0E-10
      double precision, dimension(DimPr, DimPr), intent(in) :: A
      double precision, dimension(DimPr, DimPr)             :: inv2x2
      double precision, dimension(DimPr,DimPr)              :: cofactor
      double precision            :: det
      
      
      det = A(1,1)*A(2,2) - A(2,1)*A(1,2)
      
      if (abs(det) .le. EPS) then
        inv2x2 = 0.0D0
        return
      end if
      
      cofactor(1,1) = +A(2,2)
      cofactor(1,2) = -A(2,1)
      cofactor(2,1) = -A(1,2)
      cofactor(2,2) = +A(1,1)
      
      inv2x2 = transpose(cofactor) / det
      
      
      !xjaci(1,1)= xjacm(2,2)/djacb
      !xjaci(2,2)= xjacm(1,1)/djacb
      !xjaci(1,2)=-xjacm(1,2)/djacb
      !xjaci(2,1)=-xjacm(2,1)/djacb
      
      return
      
    end function inv2x2
    
    subroutine DerivativesXY(Gp,InvJaco,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,dN_dxy, HesXY)
      
      implicit none
      
      double precision, dimension(DimPr,DimPr),intent(in)    :: InvJaco
      double precision, dimension(nne,totGp), intent(in)     :: dN_dxi, dN_deta
      double precision, dimension(nne,totGp), intent(in)     :: Hes_xixi, Hes_xieta, Hes_etaeta 
      integer, intent(in)                                    :: Gp !Variable en el lazo principal
      double precision, dimension(3,nne)                     :: Hesxieta
      double precision, dimension(2,nne)                     :: derst
      double precision, dimension(1,nne)                     :: Nxi, Neta
      double precision, dimension(1,nne)                     :: Hes_1, Hes_2, Hes_3 
      integer                                                :: idime, inode, jdime
      double precision, dimension(DimPr,nne), intent(out)    :: dN_dxy
      double precision, dimension(DimPr+1,nne), intent(out)  :: HesXY
      
      
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)
      
      Hes_1 = spread(Hes_xixi(:,Gp),dim = 1, ncopies= 1)
      Hes_2 = spread(Hes_xieta(:,Gp),dim = 1, ncopies= 1)
      Hes_3 = spread(Hes_etaeta(:,Gp),dim = 1, ncopies= 1)
     
      do inode = 1, nne
        Hesxieta(1,inode) = Hes_1(1,inode)
        Hesxieta(2,inode) = Hes_2(1,inode)
        Hesxieta(3,inode) = Hes_3(1,inode)
      end do
      
      
      do idime=1,DimPr
        do inode=1,nne
          dN_dxy(idime,inode) = 0.0
          derst(1,inode) = Nxi(1,inode)
          derst(2,inode) = Neta(1,inode)
          do jdime=1,2
            dN_dxy(idime,inode) = dN_dxy(idime,inode) + InvJaco(idime,jdime) * derst(jdime,inode)
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
      double precision, dimension(nevab), intent(out) :: veclo
      
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

    subroutine Galerkin(dvol, basis, dNdxy, source, Ke, Ce, Fe)
      
      implicit none
      
      double precision, intent(in) :: basis(nne), source(ndofn), dNdxy(DimPr,nne)
      double precision, intent(in) :: dvol
      integer :: inode, idofn, ievab, jevab, jnode, jdofn, i, j
      double precision ::  diff, convec, reac, cpcty
      double precision, intent(out) :: Ke(nevab,nevab), Fe(nevab), Ce(nevab,nevab)
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
                do j=1,2   
                  !write(*,"(A6,I2,A,I2,A,I2,A,I2,A3,f12.5)")'difma(',idofn,',',jdofn,',',i,',',j,') = ' ,difma(idofn,jdofn,i,j)
                  !call param_stab(idofn, jdofn, i, j, cte)         !conductivity tensor
                  !diff = diff+ dNdxy(i,inode) * cte * difma(idofn,jdofn,i,j)* dNdxy(j,jnode)
                  diff = diff+ dNdxy(i,inode) * difma(idofn,jdofn,i,j)* dNdxy(j,jnode)
                  !print"(A8, f10.5)",'Product ', cte * difma(idofn,jdofn,i,j)
                  !print*, '- - - - - - - - - - - - - - - - - - -'
                end do
              end do
              convec=0.0
              do i=1,2
                !print*,conma(idofn,jdofn,i)
                convec = convec + basis(inode) * conma(idofn,jdofn,i) * dNdxy(i,jnode)
              end do
              reac = basis(inode) * reama(idofn,jdofn) * basis(jnode)
              cpcty = basis(inode) * basis(jnode)
              Ke(ievab,jevab) = Ke(ievab,jevab) + (diff + convec + reac) * dvol
              Ce(ievab,jevab) = Ce(ievab,jevab) + cpcty * dvol     !element Capacity (Mass) matrix
            end do
          end do
          !Fe(ievab) = Fe(ievab) + basis(inode) * source(idofn) * dvol
          Fe(ievab) = Fe(ievab) + basis(inode) * force(idofn) * source(idofn) * dvol
          !Fe(ievab) = Fe(ievab) + basis(inode) * force(idofn) * dvol
        end do
      end do
      
    end subroutine Galerkin
   
    !subroutine param_stab(idofn, jdofn, i, j, cte)       
    !  !***********************************************************!
    !  !                                                           !
    !  ! Subroutine which check the dofn, x and y position in the  !
    !  ! diffusion tensor and take the coefficient to multiply     !
    !  ! the term of the PDE to its corresponding coeff.           !
    !  !                                                           !
    !  ! coeficients:                                              !
    !  !             Cuλ(h^2/ell^2),  λ  and   ell^2/λ             !
    !  !                                                           !
    !  !  λ represents the magnetic permeability µ                 !
    !  !  (call it mu in the code)                                 !
    !  !                                                           !
    !  ! h = 2^-i ; computed as in the paper                       !
    !  !***********************************************************!
    !  
    !  implicit none
    !  
    !  integer, intent(in) :: idofn, jdofn, i, j
    !  double precision, intent(out) :: cte
    !  
    !  cte = 0.0 
    !  
    !  if(idofn.eq.1)then
    !    if(jdofn.eq.1)then
    !      if(i==1 .and. j==1)then
    !        !cte = Cu*mu*(helem**2/ell**2)
    !      end if
    !     
    !      if(i==2 .and. j==2)then
    !        cte = mu
    !      endif
    !      
    !    elseif(jdofn==2)then
    !      if(i==1 .and. j==2)then
    !        !cte = Cu*mu*(helem**2/ell**2)
    !      end if
    !      
    !      if(i==2.and.j==1)then
    !        cte = mu
    !      end if
    !    end if
    !    
    !  elseif(idofn==2)then
    !    if(jdofn.eq.1)then
    !      if(i==1 .and. j==2)then
    !        cte = mu
    !      end if
    !      
    !      if(i==2 .and. j==1)then
    !        !cte = Cu*mu*(helem**2/ell**2)
    !      endif
    !     
    !    elseif(jdofn==2)then
    !      if(i==1 .and. j==1)then
    !        cte = mu
    !      end if
    !      
    !      if(i==2.and.j==2)then
    !        !cte = Cu*mu*(helem**2/ell**2)
    !      end if
    !    end if
    !    
    !  elseif(idofn==3 .and. jdofn==3)then
    !    if( i==j )then
    !      cte = ell**2/mu
    !    endif
    !    
    !  else
    !    continue
    !  end if
    !  !close(10)
    !  
    !  
    !  !9 format(A20,A6,I1,A1,I1,A1,I1,A1,I1,A1,I1,A1)
    !  !Next lines are to taste the 
    !  !print*, 'hmaxi,', h
    !  !print*, 'Cu µ h^2/ell^2', Cu*mu*(h**2/ell**2)
    !  !print*, 'ell^2/µ', ell**2/mu
    !  
    !end subroutine param_stab
    
    
    ! subroutine source_term_orig(element_nodes, source)
    !  !         source_term(idofn, source)
    !   implicit none
    !
    !   !***********************************************************!
    !   !The source term is given by:                               !
    !   !                                                           !
    !   !              u = grad(r^{2n/3}*sin(2ntheta/3))            !
    !   ! where:                                                    !
    !   ! r = sqrt(x^2 + y^2)   ;   theta = atan(y/x)               !
    !   !                                                           !
    !   !***********************************************************!
    
    !   !integer, intent(in) :: ievab
    !   integer :: i, j
    !   real    :: n
    !   real, dimension(nne,DimPr), intent(in)        :: element_nodes
    !   double precision, allocatable, dimension(:,:) :: r_coor, theta_coor, x_coor, y_coor
    !   double precision, dimension(nevab,1), intent(out)  :: source

    !   if(idofn.eq.3)goto 101

    !   allocate(r_coor(nne,1), theta_coor(nne,1))
    !   allocate(x_coor(nne,1), y_coor(nne,1))
      
    !   do j= 1, 3
    !     do i =1, nne
    !       print*, element_nodes(i,j)!LAS OPERACIONES TIENEN QUE SER ELEMENTALES Y AQUI DEBE X y Y IGUALARSE 
    !     end do
    !   end do


    !   do i =1, nne
    !     x_coor(i,1) = element_nodes(i,2)!LAS OPERACIONES TIENEN QUE SER ELEMENTALES Y AQUI DEBE X y Y IGUALARSE 
    !     y_coor(i,1) = element_nodes(i,3)!A NODEIMAP
    !   end do

    !   !ESTAS VARIABLES QUEDARIAN DE 4X4
    !   r_coor = sqrt(x_coor**2 + y_coor**2) 
    !   theta_coor = atan(y_coor/x_coor)

     

    !   !Si se construye el vector global F elemento a elemento (elemental)
    !   if(idofn.eq.1)then
    !     source = (2*n/3) * r_coor**(n/3) * sin(2*n*theta_coor/3)
    !   else
    !     source = (2*n/3*r_coor) * r_coor**(2*n/3) * cos(2*n*theta_coor/3)
    !   end if

    !   !Si se construye el vector global F directamente (sin proyección elemental)
    !   ! i = 1
    !   ! do j =1, ndofn - 3
    !   !     source(i,1) = (2*n/3) * r_coor**(n/3) * sin(2*n*theta_coor/3)
    !   !     source(i+1,1) = (2*n/3*r_coor) * r_coor**(2*n/3) * cos(2*n*theta_coor/3)
    !   !     source(i+2,1) = A_F(j*3,1)
    !   !     i=i+3
    !   ! end do
    
    !   101 continue
    
    ! end subroutine source_term_orig
    
    subroutine source_term(igaus, source)
      
      implicit none
      
      !***********************************************************!
      !The source term is given by:                               !
      !                                                           !
      !              u = (fi*psi',-fi'*psi)                       !
      !                                                           !
      !                                                           !
      !   f = Lu       ;   where L is the diferential operator    !
      !                                                           !
      !***********************************************************!
      
      integer, intent(in) :: igaus
      double precision, dimension(totGp) :: x_coor, y_coor
      double precision :: dey_dydx, dex_dy2, dey_dx2, dex_dxdy
      double precision :: x, y
      double precision, dimension(ndofn), intent(out)  :: source
      
      
      x_coor = ngaus(:,1)
      y_coor = ngaus(:,2)
      
      x = x_coor(igaus)               ! xi-coordinate of point j 
      y = y_coor(igaus)               ! eta-coordinate of point j 
      
      
      !Derivatives in x-direction
      dey_dydx = 2-4*y
      dex_dy2  = 0.0
      
      !Derivatives in y-direction
      dey_dx2  = 0.0
      dex_dxdy = 4*x-2
      
      source(1) = mu*(dey_dydx - dex_dy2)
      source(2) = mu*(-dey_dx2  + dex_dxdy)
      source(3) = 0.0 !force(ndofn)
      
    end subroutine source_term
    
    !subroutine source_term(igaus, source)
    !  implicit none
    !  
    !  !***********************************************************!
    !  !The source term is given by:                               !
    !  !                                                           !
    !  !              u = grad(r^{2n/3}*sin(2ntheta/3))            !
    !  ! where:                                                    !
    !  ! r = sqrt(x^2 + y^2)   ;   theta = atan(y/x)    
    !  !                                                           !
    !  ! and                                                       !
    !  !                                                           !
    !  !   f = Lu       ;   where L is the diferential operator    !
    !  !                                                           !
    !  !***********************************************************!
    !  
    !  !integer, intent(in) :: ievab
    !  integer, intent(in) :: igaus
    !  double precision, dimension(totGp) :: x_coor, y_coor
    !  !double precision, dimension(totGp) :: x, y
    !  real    :: n
    !  !integer :: i
    !  double precision :: dey_dydx, dex_dy2, dex_dx2, dey_dxdy, dey_dx2, dex_dxdy, dex_dydx, dey_dy2 
    !  double precision :: x, y, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, mm, exp_1, exp_2
    !  double precision, dimension(ndofn), intent(out)  :: source
    !  
    !  
    !  x_coor = ngaus(:,1)
    !  y_coor = ngaus(:,2)
    !  
    !  x = x_coor(igaus)               ! xi-coordinate of point j 
    !  y = y_coor(igaus)               ! eta-coordinate of point j 
    !  
    !  !terms for derivatives
    !  n  = n_val
    !  aa = (2.0*n**2)/27.0
    !  bb = (2.0*n)/27.0
    !  cc = x**2.0*(4.0*n + 3.0) - y**2.0*(n+3.0)
    !  dd = atan(y/x)
    !  ee = sin(2.0*n/3.0 * dd)
    !  ff = cos(2.0*n/3.0 * dd)
    !  gg = (x**2 + y**2)
    !  exp_1 = -(2.0 + n/6.0)
    !  exp_2 = (n/3.0 - 5.0/2.0)
    !  hh = (2.0*n**2 - 9.0*n + 9.0)
    !  ii = (4.0*n**2 - 6.0*n + 9.0)
    !  jj = (x**2 - y**2)
    !  kk = x**2.0*(n+3.0) - y**2.0*(4.0*n+3.0)
    !  ll = (8.0*n**2.0 - 24.0*n +27)
    !  mm = 8.0*n*x*y*(n-3.0)
    !  
    !  !Derivatives in x-direction
    !  dey_dydx = bb * gg**exp_2 *( x*y* ll * ff + 4.0*n*(n-3.0)*jj * ee )  
    !  dex_dy2  = aa * gg**exp_1 *( cc * ee - 4*x*y*(n+3.0) * ff)
    !  dex_dx2  = aa * gg**exp_1 *( kk * ee - 4*x*y*(n+3.0) * ff)
    !  dey_dxdy = bb * gg**exp_2 *( x*y* ll * ff + 4.0*n*(n-3.0)*jj * ee )
    !  
    !  !Derivatives in y-direction
    !  dey_dx2  = bb * gg**exp_2 * ( (2.0*x**2 * hh - y**2 * ii)*ff - mm * ee )
    !  dex_dxdy = aa * gg**exp_1 * ( 2*(n+3)* gg * ff + x*y*(5*n+ 6) * ee )
    !  dex_dydx = aa * gg**exp_1 * ( 2*(n+3)* gg * ff + x*y*(5*n+ 6) * ee )
    !  dey_dy2  = bb * gg**exp_2 * ( (x**2 * ii - 2.0*y**2 * hh)*ff - mm * ee ) 
    !  
    !  source(1) = mu*dey_dydx + mu*dex_dy2 +  Cu*mu*(helem**2/ell**2) * (dex_dx2 + dey_dxdy )
    !  source(2) =-mu*dey_dx2 + mu*dex_dxdy +  Cu*mu*(helem**2/ell**2) * (dex_dydx - dey_dy2 )
    !  if(ndofn.eq.3)then
    !    source(3) = force(ndofn)
    !  elseif(ndofn.eq.2)then
    !    continue
    !  else
    !    print*, 'Source term is a bidimensional field, not enough DoF'
    !  end if
    !  
    !  ! print*, ' Se imprime el termino de fuente '
    !  ! do i =1,ndofn
    !  !   print*, source(i)
    !  ! end do
    !  
    !end subroutine source_term
    
    
    subroutine Galerkin_Init_Cond(dvol, basis, u0_cond, C0e, u0e)
     
      !Esta rutina proyecta la condicion inicial al dominio de elementos mediante Galerkin        
      implicit none
      
      double precision, intent(in) :: basis(nne), u0_cond(nne)
      double precision, intent(in) :: dvol
      !real, intent(in)             :: u0_cond
      integer :: inode, idofn, ievab, jevab, jnode, jdofn
      double precision :: cpcty
      
      double precision, intent(out) :: u0e(nevab), C0e(nevab,nevab)
      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              cpcty = basis(inode) * basis(jnode)
              C0e(ievab,jevab) = C0e(ievab,jevab) + cpcty * dvol       !element Capacity (Mass) matrix
            end do
          end do
          u0e(ievab) = u0e(ievab) + basis(inode) * u0_cond(inode) * dvol
          !u0e(ievab) = u0e(ievab) + basis(inode) * u0_cond * dvol
        end do
      end do
      
    end subroutine Galerkin_Init_Cond
    
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
      
      
      ! SUPG
      if(kstab.eq.1) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        pertu=prod1
        
        ! galerkin least squares
      else if(kstab.eq.2)then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            !call param_stab(jdofn,idofn,k,l,cte1)
            !prod2=prod2+cte1*difma(jdofn,idofn,k,l)*workm(k,l)
            prod2=prod2+difma(jdofn,idofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(jdofn,idofn)*basis
        pertu=-prod2+prod1+prod3
        
        ! subgrid scale & taylor galerkin
      else if((kstab.eq.3).or.(kstab.eq.5)) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(idofn,jdofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            !call param_stab(idofn,jdofn,k,l,cte2)
            !prod2=prod2+cte2*difma(idofn,jdofn,k,l)*workm(k,l)
            prod2=prod2+difma(idofn,jdofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(idofn,jdofn)*basis
        pertu=prod2+prod1-prod3
        
        ! characteristic galerkin
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
    
    
    subroutine TauMat(hmaxi,tauma)
      
      !Matrix of intrinsic time scales, computed as
      ! TAU = PATAU * [ 4 K / h^2 + 2 A / h + S ]^{-1}
     
      implicit none
      
      double precision, intent(in) :: hmaxi
      integer :: i, j, k
      double precision :: chadi(3,3), chaco(3,3), chare(3,3), tauin(3,3)
      double precision :: a, b, c, tau, det                     !hnatu -> declarado en parameters
      double precision, intent(out) :: tauma(3,3)               !ndofn -> en parameters
      
      tauma = 0.0
      if(kstab.eq.0) return
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
              &difma(i,k,1,2)*difma(k,j,1,2)*2.0 + difma(i,k,2,1)*difma(k,j,2,1)*2.0 + &
              &difma(i,k,2,2)*difma(k,j,2,2)
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
      
      ! Matrix tau, corresponding to:
      ! KTAUM = 0: T = t I, where t is the minimum of all the admissible tau's
      !       = 1: T = diag(t1,t2,t3), ti is the minimum of the admissible tau's for the i-th row (equation)
      !       = 2: T = [ 4 K / h^2 + 2 A / h + S ]^{-1}
      
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
        !call param_stab(1,1,1,1,cte)
        !a = 1.0/(patau*cte*difma(1,1,1,1) /(hmaxi*hmaxi) + reama(1,1))
        a = 1.0/(patau*difma(1,1,1,1) /(hmaxi*hmaxi) + reama(1,1))
        tauma(1,1) = a
        tauma(2,2) = a
        a = (hmaxi*hmaxi*hmaxi*hmaxi)/(patau*patau)
        a = a*(patau/(hmaxi*hmaxi*reama(1,1)) + 1.0d0/(difma(1,1,1,1))) !cte*(difma(1,1,1,1)))
        tauma(3,3) = a
      end if
      
    end subroutine TauMat
    
    
    subroutine Stabilization(dvolu, basis, derxy,HesXY,source, tauma, Ke,Fe)
      !subroutine Stabilization(dvolu, basis, derxy,HesXY,tauma,Ke,Fe,pertu,workm,resid)
      
      ! Contribution to the system matrix and RHS from the stabilization term
      
      implicit none
      
      double precision, intent(in)  :: basis(nne), derxy(DimPr,nne), HesXY(3,nne), tauma(3,3)
      double precision, intent(in)  :: dvolu, source(ndofn)
      double precision              :: pertu(nevab,ndofn), workm(2,2),  resid(ndofn,nevab)
      double precision              :: prod1, prod2, prod3
      integer                       :: ievab, inode, idofn, jdofn, jevab, jnode, k, l
      double precision, intent(out) :: Ke(nevab,nevab), Fe(nevab)
      
      ! integer :: nnode,ndofn,nevab,kstab,n_ini
      !difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3)
      !common/proper/difma,conma,reama,force
      
      ievab = 0
      print*, 'in'
      pertu =  0.0
      
      do inode=1,nne
        workm(1,1)=HesXY(1,inode)
        workm(2,2)=HesXY(3,inode)
        workm(1,2)=HesXY(2,inode)
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
                !call param_stab(jdofn,idofn,k,l,cte)
                !prod3 = prod3 + cte*difma(jdofn,idofn,k,l)*workm(k,l)
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
              prod1 = prod1 + pertu(ievab,k) * tauma(k,l) * force(l) * source(l)
            end do
          end do
          Fe(ievab) = Fe(ievab) + prod1 * dvolu
        end do
      end do
      
    end subroutine Stabilization
    
    
    subroutine VinculBVs(  BVs, nofix, ifpre, presc )
      
      implicit none
      
      !integer, intent(in)      :: nBvs, nBVscol ya no se ponen estan en el modulo parameter y se comunica el valor
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
      
      print*,'!================ Bandwidth Info ==============!'
      
      write(*,"(A15,9X,I6,1X,A9)")' - UpBand:      ', upban,'   '
      write(*,"(A15,9X,I6,1X,A9)")' - LowBand:     ', lowban,'  '
      write(*,"(A15,9X,I6,1X,A9)")' - TotBand:     ', totban,'  '
      write(*,"(A15,9X,I6,1X,A9)")' - ledimAK:     ', ldAKban,' '
      
    end subroutine BandWidth
    
    subroutine TimeContribution(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, delta_t, ugl_pre, A_F)
      
      implicit none
      
      double precision, allocatable, dimension(:,:), intent(in out) :: ugl_pre
      double precision, dimension(nne,TotGp), intent(in):: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in):: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn)        :: source
      double precision, dimension(nne)          :: basis
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce, rhs_CN
      double precision, dimension(nevab)        :: Fe, Fe_time, ue_pre, time_cont
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)    :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ, delta_t
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_F
      
      allocate(A_F(ntotv, 1) )
      
      A_F = 0.0
      do ielem = 1, nelem 
        Ke = 0.0; Fe = 0.0; Ce = 0.0
        call SetElementNodes(ielem, element_nodes, nodeIDmap)
        
        call gather(nodeIDmap, ugl_pre, ue_pre)
        time_cont = ue_pre * 1.0/delta_t
        !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
        do igaus = 1, TotGp
          Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
          detJ = m22det(Jaco)
          Jinv = inv2x2(Jaco)
          dvol = detJ *  weigp(igaus,1)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          call source_term(igaus, source)
          call Galerkin(dvol, basis, dN_dxy, source, Ke, Ce, Fe) !amate lo llame Ke
          call TauMat(hmaxi,tauma)
          !call Stabilization(dvol, basis, dN_dxy, HesXY, source, tauma, Ke, Fe)

          select case(theta)
            case(2)
              Fe_time = Fe + matmul(Ce,time_cont)
            case(3)
              rhs_CN  = (1.0/delta_t)*Ce - 0.5*Ke
              Fe_time = 0.5*Fe + matmul(rhs_CN,ue_pre)
          end select

          !Fe_time = Fe + matmul(Ce,time_cont)
        end do
        
        call Assemb_Glob_Vec(nodeIDmap, Fe_time, A_F) !Assemble Global Source vector F
        
      end do
      
    end subroutine TimeContribution
    
    
    subroutine GlobalSystem(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      
      implicit none
      
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn)        :: source
      double precision, dimension(nne)          :: basis
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce
      double precision, dimension(nevab)        :: Fe
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)    :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ, aaa
      
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_K, A_C, A_F
      
      allocate(A_K(ldAKban,ntotv), A_C(ldAKban,ntotv), A_F(ntotv, 1))
      
      !duda: Fe se declaró como a(n) y en la rutina assembleF como a(n,1)
      !      pero compila y ejecuta bien. ¿Poooor?
      
      A_K = 0.0
      A_F = 0.0
      do ielem = 1, nelem 
        !gather
        Ke = 0.0    !Esto es amate
        Fe = 0.0    !Fe(nevab)
        Ce = 0.0    !elemental capacity matrix (not used in static case)
        call SetElementNodes(ielem, element_nodes, nodeIDmap)
        !do-loop: compute element stiffness matrix Ke
        do igaus = 1, TotGp
          Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
          detJ = m22det(Jaco)
          Jinv = inv2x2(Jaco)
          dvol = detJ *  weigp(igaus,1)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          call TauMat(hmaxi,tauma)
          call source_term(igaus, source)
          call Galerkin(dvol, basis, dN_dxy, source, Ke, Ce, Fe) !amate lo llame Ke
          !call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe, pertu,workm,resid)
          if(kstab.ne.0)call Stabilization(dvol, basis, dN_dxy, HesXY, source, tauma, Ke, Fe)
          !call Stabilization(dvol, basis, dN_dxy, HesXY, source, tauma, Ke, Fe)
        end do
        
        call Assemb_Glob_Mat(nodeIDmap, Ce, A_C)     !Assemble Global Conductivity Matrix K
        call Assemb_Glob_Mat(nodeIDmap, Ke, A_K)     !Assemble Global Conductivity Matrix K
        call Assemb_Glob_Vec(nodeIDmap, Fe, A_F)     !Assemble Global Source vector F
      end do
      
      aaa = maxval(coord(:,2))*2**(-i_exp) 
      if(aaa.ne.hmaxi) write(*,'(A)') '> > >Element size does not match'
      
    end subroutine GlobalSystem
    
    subroutine Assemb_Glob_Mat(lnods,Ke,A_K)
      !subroutine Assemb_Glob_Mat(ielem,lnods,Ke,A_K)
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
    end subroutine Assemb_Glob_Mat

    subroutine Assemb_Glob_Vec( nodeIDmap, fe, F_global)
      
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
      
    end subroutine Assemb_Glob_Vec
    
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
      integer :: ivfix, idofn, itotv, jband, nvfix, ktotv, ipoin
      double precision,intent(inout)  :: rigid(ldAKban,ntotv), gload(ntotv)
      
      nvfix = nBVs
      
      !***  Llac sobre els nodes coaccionats (Lazo sobre los nodos preescritos)
      do ivfix=1,nvfix
        ipoin=nofix(ivfix)
        do idofn=1,ndofn !3
          if (ifpre(idofn,ivfix).eq.1) then
            
            !***  Ca de grau de llibertat prescrit
            itotv=(ipoin-1)*ndofn+idofn
            do ktotv=1,ntotv
              jband = ktotv-itotv + totban
              if((lowban +1.le. jband).and.(jband.le.ldAKban))then
                gload(ktotv) = gload(ktotv)-rigid(jband,itotv)*presc(idofn,ivfix)
                rigid(jband,itotv)=0!7.777
              end if
            end do
            do ktotv= 1, ntotv
              jband = itotv-ktotv + totban
              if((lowban +1.le. jband).and.(jband.le.ldAKban))then
                rigid(jband,ktotv)=0!5.555 !columna K
              end if
            end do
            rigid(totban,itotv)=1.0
            gload(itotv)=presc(idofn,ivfix)
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
        stop
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
    
    subroutine Convergence(solution, x, y, exact_x, exact_y, error)
      implicit none
      
      double precision, dimension(ntotv,1), intent(in) :: solution
      double precision, dimension(nnodes), intent(in)  :: x, y
      !declaracion de variables relacionadas con la solucion exacta
      !double precision                                 :: aa, bb, cc, dd, ee, exp1, exp2
      double precision                                 :: fi, psi, der_fi, der_psi 
      integer                                          :: i
      !real                                             :: n
      !end variables for exact solution
      
      
      double precision, dimension(1, ntotv)   :: solution_T
      double precision                        :: x_FEM, y_FEM, ex_x, ex_y
      
      double precision, dimension(nnodes), intent(out) :: exact_y, exact_x
      double precision, intent(out)                    :: error
      
      solution_T = transpose(solution)
      error = 0.0 
      !! * * * * * * * * * Exact solution * * * * * * * * * * * * *
      ! L- domain w/singular solution 
      !n    = n_val
      !aa   = (2.0*n)/3.0
      !exp1 = -(n/6.0)
      !exp2 = n/3.0 - 1.0/2.0
      
      !do i = 1, nnodes
      !  bb = atan(y(i)/x(i))
      !  cc = sin(2.0*n/3.0 * bb)
      !  dd = cos(2.0*n/3.0 * bb)
      !  ee = (x(i)**2 + y(i)**2)
      !  
      !  exact_x(i) = aa*ee**exp1*cc 
      !  exact_y(i) = aa*ee**exp2*dd
      !  
      !end do
      
      !A simple test u(fi*psi',fi'*psi) 
      !write(*,*) '       FEMx','            Ex_x', '            FEMy','           Ex_y'
      do i = 1, nnodes  !simple Function
        fi  = x(i)*(1.0-x(i)) 
        psi = y(i)*(1.0-y(i))
        der_fi = (1.0-2*x(i))
        der_psi = (1.0-2*y(i))
        
        ex_x = (fi * der_psi)
        ex_y = -(der_fi * psi)
        x_FEM = solution_T(1,ndofn*i-2)
        y_FEM = solution_T(1,ndofn*i-1)
        !write(*,"(4(f15.5,1x))") x_FEM, ex_x, y_FEM, ex_y
        error = error + (x_FEM - ex_x)**2 + (y_FEM - ex_y)**2 
        
        exact_x(i) = (fi * der_psi)
        exact_y(i) = -(der_fi * psi)
        
      end do
      error = sqrt(error/float(nnodes))
      !! * * * * * * * * * (end) Exact solution * * * * * * * * * * 
      
      
    end subroutine Convergence
    
    
    subroutine Res_Matlab(solution)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Res/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcoor, ycoor
      double precision, dimension(nnodes)     :: exact_y, exact_x
      double precision, dimension(nnodes)     :: x, y
      character(len=12)                       :: extension
      character(len=12)                        :: coord_name, conec_name
      character(len=10)                        :: error_name
      character(len=10)                        :: identy
      integer                                 :: i,j, ipoin, ans
      double precision                        :: error, xmax
      
      
      extension = "_matlab.txt"
      
      solution_T = transpose(solution)
      xcoor = spread(coord(:,2),dim = 1, ncopies= 1)
      ycoor = spread(coord(:,3),dim = 1, ncopies= 1)
      
      x  = xcoor(1,:)
      y  = ycoor(1,:)
      
      call Convergence(solution, x, y, exact_x, exact_y,  error)
      
      !write(*,*) '       FEMx','            Ex_x', '            FEMy','           Ex_y'
      
      open(unit=111, file= fileplace//File_PostProcess//extension, ACTION="write", STATUS="replace")
      do ipoin = 1, nnodes  !   uh_x    uh_y    uex_x   uex_y
        write(111,906) solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1),&
        &              exact_x(ipoin), exact_y(ipoin)
        !write(*,"(4(f15.5,1x))") solution_T(1, ndofn*ipoin-2), exact_x(ipoin),&
        !&                      solution_T(1,ndofn*ipoin-1), exact_y(ipoin)
      end do
      !print*, ' '
      !print*, '!====== Matlab file ======'
      print"(A6,A24,A30)", ' File ',File_PostProcess//extension, 'written succesfully in Res/ '
      print*, ' '
      
      close(111)
      
      write(*,*) '!====== Element size'
      read(*,*) identy
      
      write(*,*)'Mesh File?. Y=1, N=2'; read(*,*) ans
      if(ans.eq.1)then
        
        write(*,*) '!====== Name of coordinates file'
        read(*,*) coord_name
        write(*,*) '!====== Name of connectivity file'
        read(*,*) conec_name
        
        open(unit=333, file= fileplace//conec_name//extension, ACTION="write", STATUS="replace")
        open(unit=444, file= fileplace//coord_name//extension, ACTION="write", STATUS="replace")
        
        
        do i=1,nelem
          write(333,902) (lnods(i,j+1), j =1,nne)
        end do
        
        do ipoin = 1, nnodes
          write(444,904) ipoin, x(ipoin), y(ipoin)
        end do
        print*, ' '
        print*, '!====== Mesh file ======'
        print"(A6,A12,A3,A12,A30)", ' File ',coord_name, ' and ', conec_name, 'written succesfully in Res/ '
        print*, ' '
        
        close(444)
      else
        continue
      end if
      
      error_name = 'err_'//identy
      open(unit=777, file= fileplace//error_name//extension, ACTION="write", STATUS="replace")
      write(777,"(1x,E15.5)") error
      close(777)

      xmax = maxval(coord(:,2)) !the greatest number in x column
      print*, ' '
      write(*,"(A9,f7.3,A25,E15.5)")' For h = ', xmax*2**(-i_exp), 'the error estimation is ', error 
      print*, ' '
     
      902 format(10(1x,I7) ) !format for msh
      904 format(I7,2(3x,f9.4) ) !format for msh
      906 format(4(F25.10, 3x))
      
    end subroutine Res_Matlab
    
    subroutine PostProcess(solution, activity)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Pos/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      character(*), intent(in)                :: activity
      character(len=10)                       :: extension 
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcor, ycor
      integer                                 :: ipoin, ii
      
      solution_T = transpose(solution)
      xcor  = spread(coord(:,2),dim = 1, ncopies= 1)
      ycor  = spread(coord(:,3),dim = 1, ncopies= 1)
      
      
      if(activity == "msh")then !quitar este if y acomodar el numero de unidad
        extension = ".post.msh"
        open(unit=555, file= fileplace//File_PostProcess//extension, ACTION="write", STATUS="replace")
        
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
        print"(A6,A24,A30)", ' File ',File_PostProcess//'.post.msh','written succesfully in Pos/ '
       
      elseif(activity == "res")then
        extension = ".post.res"
        open(unit=555, file= fileplace//File_PostProcess//extension, ACTION="write", STATUS="replace")
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
            !An alternative work around is by explicitly designate the 
            !elements to be read using an io-implied-do
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
            do ipoin = 1, nnodes
              write(555,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1)
            end do
            write(555,"(A)") 'End Values'
            write(555,"(A)") 'Result "P" "Preassure" 0 Scalar OnNodes'
            write(555,"(A)") 'ComponentNames "" '
            write(555,"(A)") 'Values'
            write(555,*) '#',   'No    ','             p '
            !  se escribe el res para el caso escalar de un grado de libertad
            write(555,914)
            ii=1
            do ipoin = 3, nnodes*3,3
              write(555,916) ii, solution_T(1,ipoin)
              ii=ii+1
            end do
            write(555,"(A)") 'End Values'
            print"(A6,A24,A30)", ' File ',File_PostProcess//'.post.res', 'written succesfully in Pos/ '
        end select
        
        close(555)
      else
        write(*,"(A)") ' < < Error > > Postprocess activity must be "msh" or "res" '
        close(555)
        stop
      end if
    
      
      900 format(A15, A13, A1, A13)
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      906 format(I7,2(3x,f9.4)) !format for msh
      908 format(10(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,E12.5)  !format for scalar case
      918 format(I7,3x,E12.5,3x,E12.5) !format for res velocity
      919 format(I7,3(4x,F25.10)) !format for res velocity
    
    end subroutine PostProcess
    
    subroutine GID_PostProcess(solution, activity, step_value, interval, time_final)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Pos/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      character(*), intent(in)                :: activity
      integer, intent(in)                     :: step_value
      real, intent(in)                        :: interval, time_final
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcor, ycor
      integer                                 :: ipoin, ii
      
      solution_T = transpose(solution)
      xcor  = spread(coord(:,2),dim = 1, ncopies= 1)
      ycor  = spread(coord(:,3),dim = 1, ncopies= 1)
      
      
      if(activity == "msh")then !quitar este if y acomodar el numero de unidad
        open(unit=100, file= fileplace//File_PostProcess, ACTION="write", STATUS="replace")
        
        write(100,902) 'MESH', '"Domain"', 'dimension', DimPr, 'ElemType', ElemType, 'Nnode', nne
        write(100,"(A)") '#2D Convection-Diffusion-Reaction'
        write(100,900) '#Element tipe: ', ElemType,'/',ElemType
        write(100,"(A)")'Coordinates'
        write(100,"(A)") '#   No        X           Y'
        do ipoin = 1, nnodes
          write(100,906) ipoin, xcor(1,ipoin), ycor(1,ipoin)
        end do
        write(100,"(A)") 'End Coordinates'
        write(100,"(A)") 'Elements'
        do ipoin = 1, nelem
          write(100,908) lnods(ipoin,:)
        end do
        write(100,"(A)") 'End Elements'
        close(100)
        print"(A11,A19,A30)", ' Mesh file ',File_PostProcess//'.post.msh', 'written succesfully in Pos/ '
      ! if(status.eq.0)then continue 
      elseif(activity == "res")then
       
        if(step_value == 0)then
          open(unit=200, file= fileplace//File_PostProcess, ACTION="write", STATUS="replace")
          write(200,"(A)") 'GiD Post Results File 1.0'
          write(200,"(A)") '#2D Convection-Diffusion-Reaction'
        else
          continue
        endif
        open(unit=200, file= fileplace//File_PostProcess, ACTION="write", STATUS="old", position="append")
        
        ! se escribe el res de las componentes de la velocidad
        select case(ndofn)
          case(1)
            write(200,"(A29, I3, A)") 'Result "DoF" "Concentration" ', step_value,' Scalar OnNodes'
            write(200,"(A)") 'ComponentNames "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             Ex '
            !  se escribe el res para el caso escalar de un grado de libertad
            write(200,914)
            do ipoin = 1, nnodes
              write(200,916) ipoin, solution(ipoin, 1)
            end do
            !An alternative work around is to explicitly designate the elements to be read using an io-implied-do.
            !Something like
            !read (unit=10, fmt=*, iostat=iostat) (mat(pcnt,i),i=1,m)
            write(200,"(A)") 'End Values'
          case(2)
            write(200,"(A29, I3, A)") 'Result "DoF" "Concentration" ', step_value,' Vector OnNodes'
            write(200,"(A)") 'ComponentNames "u" "v" "--" "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             Ex ','               Ey '
            do ipoin = 1, nnodes
              write(200,918) ipoin, solution_T(1, ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
            end do
            write(200,"(A)") 'End Values'
          case(3)
            write(200,"(A29, I3, A)") 'Result "DoF" "Concentration" ', step_value,' Vector OnNodes'
            write(200,"(A)") 'ComponentNames "u" "v" "w" "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             Ex ','               Ey'
           ! do ipoin = 1, nnodes
           !   write(200,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
           ! end do
            do ipoin = 1, nnodes
              write(200,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1)
            end do
            write(200,"(A)") 'End Values'
            write(200,"(A22, I3, A)") 'Result "P" "Preassure"', step_value,' Scalar OnNodes'
            write(200,"(A)") 'ComponentNames "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','     p '
            !  se escribe el res para el caso escalar de un grado de libertad
            write(200,914)
            ii=1
            do ipoin = 3, nnodes*3,3
              write(200,916) ii, solution_T(1,ipoin)
              ii=ii+1
            end do
            write(200,"(A)") 'End Values'
        end select
      else
        write(*,"(A)") ' < < Error > > Postprocess activity must be "msh" or "res" non ', activity
        close(200)
        stop
      end if
      close(200)
      !la siguiente instruccion debe usarse con nt no con time pero solo es para avanzar
      if(interval == time_final) then
        print*, ' '
        print"(1x, A26,A30)", File_PostProcess//'.post.res', 'written succesfully in Pos/ '
      endif
      
      
      900 format(A15, A13, A1, A13)
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      906 format(I7,2(3x,f9.4)) !format for msh
      908 format(9(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,F12.5)  !format for scalar case
      918 format(I7,3x,E15.5,3x,E15.5) !format for res velocity
      919 format(I7,3(4x,F20.10)) !format for res velocity
      
    end subroutine GID_PostProcess
   
    
    
  !Fin de contains


end module library
