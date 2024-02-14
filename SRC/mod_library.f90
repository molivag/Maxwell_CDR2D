module library
  use param
  use geometry
  use biunit

  contains
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine ReadFile(NumRows, NumCols, condition, condition_value)
      
      integer :: i, j, stat1, stat2
      integer, intent(in)            :: NumRows, NumCols
      character(len=*), parameter    :: fileplace = "./"
      character (len=9)              :: FileName1, FileName2
      character(len=180) :: msg
      integer, dimension(NumRows,NumCols-ndofn), intent (out) :: condition
      double precision, dimension(NumRows,ndofn),intent (out) :: condition_value
      
      FileName1 ='ifpre.dat' 
      FileName2 ='BoVal.dat'
      
      open(unit=10,file=fileplace//FileName1, status='old', action='read', iostat=stat1, iomsg=msg)
      open(unit=20,file=fileplace//FileName2, status='old', action='read', iostat=stat2, iomsg=msg)
      
      read(10,*,iostat=stat1,iomsg=msg) ((condition(i,j), j=1,NumCols-ndofn), i=1,NumRows)
      if (stat1.ne.0) then
        print*, ' '
        print*,'!============ STATUS READING BOUND VAL. ============!'
        print "(A38,I2)", "- Status while reading ifPre file is: ", stat1
        print*, ' '
        print'(A8,1x,A180)','iomsg= ',msg
      else
        continue
      end if
      
      read(20,*,iostat=stat2,iomsg=msg) ((condition_value(i,j), j=1,ndofn), i=1,NumRows)
      if (stat2.ne.0) then
        print "(A38,I2)", "- Status while reading BoVal file is: ", stat2
        print*, ' '
        print'(A8,1x,A180)','iomsg= ',msg
      else
        continue
      end if
      
      
      close (10)
      close (20)
      
    end subroutine ReadFile
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine SetElementNodes(elm_num, element_nodes, nodeIDmap, xi_cor, yi_cor)
      
      implicit none
      
      ! number of element for each elem integral in loop of K global
      integer,intent(in)                                   :: elm_num 
      integer                                              :: i,j,global_node_id
      double precision, dimension(DimPr,nne)               :: aaa
      double precision, dimension(nne,DimPr), intent(out)  :: element_nodes
      double precision,dimension(nne), intent(out)         :: xi_cor, yi_cor
      integer,dimension(nne), intent(out)                  :: nodeIDmap
      
      element_nodes = 0.0
      nodeIDmap = 0
      
      do i = 1, nne
        global_node_id = lnods(elm_num,i)
        nodeIDmap(i)   = global_node_id
        do j=1 ,DimPr
          aaa(j,i) = coord(j,global_node_id)
        end do
        xi_cor(i) = aaa(1,i)
        yi_cor(i) = aaa(2,i)
      end do
      
      element_nodes = transpose(aaa)
      
    end subroutine SetElementNodes
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine Jacobian( element_nodes, dN_dxi, dN_deta, Gp, xjacm, djacb, xjaci)
      implicit none
      
      integer, intent(in)      :: Gp ! se usara en el lazo principal con el punto de Gauss
      double precision, dimension(nne,DimPr), intent(in) :: element_nodes
      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(DimPr,nne)  :: derst!, Basis2D, derst
      double precision, dimension(DimPr,nne)  :: elcod
      double precision, dimension(DimPr,DimPr), intent(out) :: xjacm, xjaci
      double precision                        , intent(out) :: djacb
      double precision, parameter                           :: EPS = 1.0E-10
      integer                                               :: idime, jdime, inode
      
      !double precision, dimension(1,nne)      :: Nxi, Neta
      !double precision, dimension(DimPr,DimPr):: J2D
      
      !con estas instrucciones extraigo la columna de Nx como renglon y lo guardo en Nxi, Gp se
      !ira moviendo conforme la funcion J2D sea llamada en el lazo principal para 
      !cada elemento lo mismo para Neta con dN_deta
      
      !Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      !Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)
      
      !Las siguientes tres lineas realizan de forma implicita el calculo de las derivadas
      !espaciales es decir dN/dx and dN/dy (eq. 5.114 - 5.117). Las derivadas espaciales
      !no se calcula explicitamente, en su lugar se usa:
      
      !            d/dy = (d/deta)J^-1      (ver eq. 5.77 y 5.109)
      
      !Basis2D(1,:) = Nxi(1,:)
      !Basis2D(2,:) = Neta(1,:)
      !
      !J2D = matmul(Basis2D,element_nodes)
      
      
      do jdime = 1,2
        do inode = 1,nne
          elcod(jdime,inode) = element_nodes(inode, jdime)
        end do
      end do
     
      do inode = 1,nne
        derst(1,inode) = dN_dxi(inode,Gp) 
        derst(2,inode) = dN_deta(inode,Gp) 
      end do
      !print*," " 
      !print"(A8, 9(1x,f9.5))","dN_dxi: ", (derst(1,inode), inode=1,nne)
      !print"(A9, 9(1x,f9.5))","dN_deta: ", (derst(2,inode), inode=1,nne)
     
      ! Jacobian Matrix
      do idime=1,2
        do jdime=1,2
          xjacm(idime,jdime)=0.0
          do inode=1,nne
            !xjacm(idime,jdime) = xjacm(idime,jdime) + derst(idime,inode)*element_nodes(inode,jdime)
            xjacm(idime,jdime) = xjacm(idime,jdime) + derst(idime,inode)*elcod(jdime,inode)
          end do
        end do
      end do
      
      !do idime=1,2
      !  write(*,"(f10.5, 1x, f10.5)" )( xjacm(idime,jdime) ,jdime=1,2)
      !end do
      
      ! Determinant of Jacobian Matrix
      djacb  = xjacm(1,1)*xjacm(2,2) - xjacm(1,2)*xjacm(2,1)
      !print"(A8,f15.5)", 'det Jaco: ',djacb 
      if(djacb.lt.1.e-8)then
        write(*,*)'Element with non-positive Jacobian'
      end if
      
      ! Inverse of Jacobian
      if (abs(djacb) .le. EPS) then
        xjaci = 0.0! 0.0D0
        return
      end if
      xjaci(1,1)= xjacm(2,2)/djacb
      xjaci(2,2)= xjacm(1,1)/djacb
      xjaci(1,2)=-xjacm(1,2)/djacb
      xjaci(2,1)=-xjacm(2,1)/djacb
      
      return
      
    end subroutine Jacobian
    
    !function djacb(xjacm)
    !  
    !  implicit none
    !  
    !  double precision :: djacb
    !  double precision, dimension(2,2), intent(in)  :: xjacm
    !  integer :: idime, jdime 
    !  !m22det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
    !  do idime=1,2
    !    write(*,"(f10.5, 1x, f10.5)" )( xjacm(idime,jdime) ,jdime=1,2)
    !  end do
    !  
    !  djacb  = xjacm(1,1)*xjacm(2,2) - xjacm(1,2)*xjacm(2,1)
    !  print"(A8,f15.5)", 'det Jaco: ',djacb 
    !  if(djacb.lt.1.e-8)then
    !    write(*,*)'Element with non-positive Jacobian'
    !  end if
    ! 
    !  return
    !  
    !end function djacb
    !
    !function xjaci(detJ,xjacm)
    !  
    !  implicit none
    !  
    !  double precision, parameter :: EPS = 1.0E-10
    !  double precision, dimension(DimPr, DimPr), intent(in) :: xjacm
    !  double precision, intent(in)                          :: detJ
    !  
    !  
    !  if (abs(detJ) .le. EPS) then
    !    xjaci = 0.0! 0.0D0
    !    return
    !  end if
    !  
    !  xjaci(1,1)= xjacm(2,2)/detJ
    !  xjaci(2,2)= xjacm(1,1)/detJ
    !  xjaci(1,2)=-xjacm(1,2)/detJ
    !  xjaci(2,1)=-xjacm(2,1)/detJ
    !  
    !  return
    !  
    !end function xjaci
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine Galerkin(hmaxi, dvol, basis, dNdxy, EMsource, Ke, Ce, Fe)
      
      implicit none
      
      double precision, intent(in) :: basis(nne), EMsource(ndofn), dNdxy(DimPr,nne)
      double precision, intent(in) :: dvol, hmaxi
      integer :: inode, idofn, ievab, jevab, jnode, jdofn, i, j
      double precision ::  diff, convec, reac, cpcty, coef
      double precision, intent(in out) :: Ke(nevab,nevab), Fe(nevab), Ce(nevab,nevab)
      
      ievab=0
      do inode=1,nne
        reac = 0.0
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              diff=0.0
              x_difma_loop: do i=1,2
                z_difma_loop: do j=1,2   
                  !write(*,"(A6,I2,A,I2,A,I2,A,I2,A3,e12.5)")&
                  !&'difma(',idofn,',',jdofn,',',i,',',j,') = ',difma(idofn,jdofn,i,j)
                  call param_stab(idofn, jdofn, i, j, hmaxi, coef) !conductivity tensor
                  !print"(A8, e12.4)",'Product ', coef * difma(idofn,jdofn,i,j)
                  diff = diff + dNdxy(i,inode) * coef * difma(idofn,jdofn,i,j)* dNdxy(j,jnode)
                  !print"(A8, e12.5)",'diff ', diff
                  !print*, '- - - - - - - - - - - - - - - - - - -'
                  !print*, ' '
                end do z_difma_loop
              end do x_difma_loop
              convec=0.0
              conma_loop: do i=1,2
                !print*,conma(idofn,jdofn,i)
                convec = convec + basis(inode) * conma(idofn,jdofn,i) * dNdxy(i,jnode)
              end do conma_loop
              
              !LO IDEAL ES QUE NO SE EVALUE EL MISMO IF EN CADA GRADO DE LIBERATD SINO QUE SE
              !HAGAN VARIOS GALERKIN PARA EVITAR LA MAYOR CANTIDAD DE EVALUACIONES IF POSIBLES
              if(TwoHalf =='Y')then
                !if it is dealing with a 2.5D Problem
                if(oper == 'LAPL')then
                  ! print*,'Entra a Laplacian'
                  ! print*,k_y
                  ! print*,sigma
                  reac = basis(inode) * k_y**2 *sigma* reama(idofn,jdofn) * basis(jnode)
                elseif(oper == 'MAXW')then
                  reac = basis(inode) * lambda*k_y**2 * reama(idofn,jdofn) * basis(jnode)
                endif
              else
                print*,'!if it is NOT dealing with a 2.5D Problem'
                reac = basis(inode) * reama(idofn,jdofn) * basis(jnode)
              endif
              
              Ke(ievab,jevab) = Ke(ievab,jevab) + (diff + convec + reac) * dvol
              
              cpcty = sigma * (basis(inode) * basis(jnode) )
              Ce(ievab,jevab) = Ce(ievab,jevab) + cpcty * dvol     !element Capacity (Mass) matrix
            end do
          end do
          Fe(ievab) = Fe(ievab) + basis(inode) * EMsource(idofn) * dvol
          !Fe(ievab) = Fe(ievab) + basis(inode) * force(idofn) * dvol
        end do
      end do
      
    end subroutine Galerkin
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine param_stab(idofn, jdofn, i, j, elem_size_h, coeff)       
      !***********************************************************!
      !                                                           !
      ! Subroutine which check the dofn, x and y position in the  !
      ! diffusion tensor and take the coefficient to multiply     !
      ! the term of the PDE to its corresponding coefficient.     !
      !                                                           !
      ! coeficients:                                              !
      !             Cuλ(h^2/ell^2),  λ  and   ell^2/λ             !
      !                                                           !
      !  λ = 1/µ  = reluctivity of the medium                     !
      !                                                           !
      !  h = elem_size_h 
      !                                                           !
      !***********************************************************!
      
      implicit none
      
      integer, intent(in) :: idofn, jdofn, i, j
      double precision, intent(in) :: elem_size_h
      double precision, intent(out) :: coeff
      
      
      !print*, 'getting into param_stab'
      
      coeff = 0.0 
      !print'(A9,F10.5)', 'elem_size_h  : ', elem_size_h
      !print'(A9,F10.5)', 'elem_size_h^2: ', elem_size_h**2
      if(kstab.eq.6)then !coefficients for MVAF 
        ! print*,'!coefficients for MVAF'
        !print*,'stabi',kstab
        if(idofn.eq.1)then
          if(jdofn.eq.1)then                      !difma(idofn,jdofn,i,j)
            if(i==1 .and. j==1)then                      !difma(1,1,1,1)
              coeff = Cu*lambda*(elem_size_h**2/ell**2)
              !print'(A9,F10.5)', 'elem_size_h  : ', elem_size_h
              !print'(A2,e12.5)','Su', coeff
            end if
           
            if(i==2 .and. j==2)then                      !difma(1,1,2,2)
              coeff = lambda
              !print'(A2,e12.5)','λ ', coeff
            endif
            
          elseif(jdofn==2)then                           !difma(1,2,1,2)
            if(i==1 .and. j==2)then
              coeff = Cu*lambda*(elem_size_h**2/ell**2)
              !print'(A9,F10.5)', 'elem_size_h  : ', elem_size_h
              !print'(A2,e12.5)','Su', coeff
            end if
            
            if(i==2.and.j==1)then                        !difma(1,2,2,1)
              coeff = lambda
              !print'(A3,e12.5)','λ ', coeff
            end if
          end if
          
        elseif(idofn==2)then
          if(jdofn.eq.1)then
            if(i==1 .and. j==2)then                      !difma(2,1,1,2)
              coeff = lambda
              !print'(A3,e12.5)', 'λ ', coeff
            end if
            
            if(i==2 .and. j==1)then                      !difma(2,1,2,1)
              coeff = Cu*lambda*(elem_size_h**2/ell**2)
              !print'(A9,F10.5)', 'elem_size_h  : ', elem_size_h
              !print'(A2,e12.5)','Su', coeff
            endif
           
          elseif(jdofn==2)then
            if(i==1 .and. j==1)then
              coeff = lambda
              !print'(A3,e12.5)','λ ', coeff
            end if
            
            if(i==2.and.j==2)then                        !difma(2,2,2,2)
              coeff = Cu*lambda*(elem_size_h**2/ell**2)
              !print'(A9,F10.5)', 'elem_size_h  : ', elem_size_h
              !print'(A2,e12.5)','Su', coeff
            end if
          end if
          
        elseif(idofn==3 .and. jdofn==3)then              !difma(3,3,1,1) or difma(3,3,2,2)
          if( i==j )then
            coeff = ell**2/lambda
              !print'(A2,e12.5)','Sp', coeff
          endif
          
        else
          continue
        end if
        !close(10)
      else
        if((idofn == jdofn).and.(i == j))then
            ! write(*,"(A6,I2,A,I2,A,I2,A,I2,A3,e12.5)")&
            ! &'difma(',idofn,',',jdofn,',',i,',',j,') = ',difma(idofn,jdofn,i,j)
            !Aqui tengo que poner un identificador o LGO QUE dependiendo
            !el problema, ponga una propiedad fisica u otra
            if(kstab == 0)then
              ! print*,'!The coeficients for Laplacian operator or 2nd derivatives respect to itslefs'
              coeff = sigma
            else
              coeff = lambda
            endif
           
        end if
      end if
      
    end subroutine param_stab
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine pertur(hmaxi, idofn, jdofn, workm, derxy, basis, pertu )
      
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
      
      double precision, intent(in)     :: workm(2,2),derxy(2),basis, hmaxi
      integer                          :: idofn, jdofn, k, l
      double precision                 :: prod1, prod2, prod3, cte1, cte2
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
            call param_stab(jdofn, idofn, k, l, hmaxi, cte1) 
            prod2=prod2+cte1*difma(jdofn,idofn,k,l)*workm(k,l)
            !prod2=prod2+difma(jdofn,idofn,k,l)*workm(k,l)
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
            call param_stab(idofn, jdofn, k, l, hmaxi, cte2) 
            prod2=prod2+cte2*difma(idofn,jdofn,k,l)*workm(k,l)
            !prod2=prod2+difma(idofn,jdofn,k,l)*workm(k,l)
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine TauMat(hmaxi,tauma)
      
      !Matrix of intrinsic time scales, computed as
      ! TAU = PATAU * [ 4 K / h^2 + 2 A / h + S ]^{-1}
     
      implicit none
      
      double precision, intent(in) :: hmaxi
      integer :: i, j, k
      double precision :: chadi(3,3), chaco(3,3), chare(3,3), tauin(3,3)
      double precision :: a, b, c, tau, det             !hnatu -> declarado en parameters
      double precision, intent(out) :: tauma(3,3)       !ndofn -> en parameters
      double precision :: cte, cte1, cte2, cte3, cte4, cte5, cte6, cte7, cte8
      
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
            call param_stab(i, k, 1, 1, hmaxi, cte1)
            call param_stab(k, j, 1, 1, hmaxi, cte2)
            call param_stab(i, k, 1, 2, hmaxi, cte3)
            call param_stab(k, j, 1, 2, hmaxi, cte4)
            call param_stab(i, k, 2, 1, hmaxi, cte5)
            call param_stab(k, j, 2, 1, hmaxi, cte6)
            
            call param_stab(i, k, 2, 2, hmaxi, cte7)
            call param_stab(k, j, 2, 2, hmaxi, cte8)
            
            chadi(i,j) = chadi(i,j) + cte1*difma(i,k,1,1) * cte2*difma(k,j,1,1) + &
              &cte3*difma(i,k,1,2)*cte4*difma(k,j,1,2)*2.0 + &
              &cte5*difma(i,k,2,1)*cte5*difma(k,j,2,1)*2.0 + &
              &cte7*difma(i,k,2,2)*cte8*difma(k,j,2,2)
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
        call param_stab(1,1,1,1, hmaxi, cte)
        a = 1.0/(patau*cte*difma(1,1,1,1) /(hmaxi*hmaxi) + reama(1,1))
        !a = 1.0/(patau*difma(1,1,1,1) /(hmaxi*hmaxi) + reama(1,1))
        tauma(1,1) = a
        tauma(2,2) = a
        a = (hmaxi*hmaxi*hmaxi*hmaxi)/(patau*patau)
        !a = a*(patau/(hmaxi*hmaxi*reama(1,1)) + 1.0d0/(difma(1,1,1,1))) !cte*(difma(1,1,1,1)))
        a = a*(patau/(hmaxi*hmaxi*reama(1,1)) + 1.0d0/cte*(difma(1,1,1,1)))
        tauma(3,3) = a
      end if
      
    end subroutine TauMat
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    !
    subroutine Stabilization(hmaxi, dvolu, basis, derxy,HesXY,source, tauma, Ke,Fe)
      !subroutine Stabilization(dvolu, basis, derxy,HesXY,tauma,Ke,Fe,pertu,workm,resid)
      
      ! Contribution to the system matrix and RHS from the stabilization term
      
      implicit none
      
      double precision, intent(in)  :: basis(nne), derxy(DimPr,nne), HesXY(3,nne), tauma(3,3)
      double precision, intent(in)  :: dvolu, source(ndofn), hmaxi
      double precision              :: pertu(nevab,ndofn), workm(2,2),  resid(ndofn,nevab)
      double precision              :: prod1, prod2, prod3, cte
      integer                       :: ievab, inode, idofn, jdofn, jevab, jnode, k, l
      double precision, intent(out) :: Ke(nevab,nevab), Fe(nevab)
      
      ! integer :: nnode,ndofn,nevab,kstab,n_ini
      !difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3)
      !common/proper/difma,conma,reama,force
      
      ievab = 0
      !print*, 'in'
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
                call param_stab(jdofn, idofn, k, l, hmaxi, cte) 
                prod3 = prod3 + cte*difma(jdofn,idofn,k,l)*workm(k,l)
                !prod3 = prod3 + difma(jdofn,idofn,k,l)*workm(k,l)
              end do
            end do
            
            resid(jdofn,ievab) = prod1 + prod2 - prod3
            call pertur(hmaxi,idofn,jdofn,workm,derxy(1,inode),basis(inode),pertu(ievab,jdofn) )
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine VinculBVs(condition,  BVs, nofix, ifpre, presc )
      
      implicit none
      
      !integer, intent(in)      :: nBvs, nBVscol ya no se ponen estan en el modulo parameter y se comunica el valor
      integer, intent(in) :: condition( nBvs, nBVscol-ndofn)
      double precision, intent(in) :: BVs( nBvs, ndofn)
      integer             :: i, j
      double precision, intent(out) :: presc(ndofn,nBVs)
      integer, intent(out)          :: ifpre(ndofn,nBVs)
      integer, intent(out)          :: nofix(nBVs)
      
      
      select case(ndofn)
        case(1)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = condition(j,1)
              ifpre(i,j) = condition(j,2) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,1)
            end do
          end do
          
        case(2)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = condition(j,1)
              ifpre(i,j) = condition(j,i+1) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,i)
            end do
          end do
         
        case(3)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = condition(j,1)
              ifpre(i,j) = condition(j,i+1) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,i)
            end do
          end do
          
        case DEFAULT
          write(*,*) 'Exceeded DoF'
        end select
    end subroutine VinculBVs
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine BandWidth( )
      !this routine computes the bandwidth 
      
      implicit none
      
      integer :: iband, ielem, inode, ipoin, jnode, jpoin, nband       ! , C, D
      
      iband=0
      do ielem =1, nelem
        do inode = 1, nne
          !Este +1 es para que comience en los nodos (columna 2) y no del numeor de elemento
          ipoin = lnods(ielem,inode) 
          do jnode = 1, nne
            jpoin = lnods(ielem,jnode)
            iband = max(iband,abs(jpoin-ipoin))
          end do
        end do
        !if (iband.gt.nband) nband=iband !nband pasa como variable global por que se usa en  ApplyBVal y otros
        nband = iband
      end do
      nband=(nband+1)*ndofn-1
      
      upban  = nband
      lowban = nband
      totban = lowban + upban + 1
      ldAKban= 2*lowban + upban + 1   
      
      ! ldAKban represents the leading dimension of the Banded Global Stiffnes matrix AK_band . 
      ! This variable specifies how the AK_band matrix is stored in memory and is used to correctly 
      ! access and manipulate the AK_band matrix during the Solver process. 
      ! In Fortran, the leading dimension specifies the number of rows used for storage. 
      
      print*, ' '
      print*,'!================ Bandwidth Info ==============!'
      write(*,"(A30,2X,I6,1X,A9      )")' - Uppper Band              : ', upban,'   '
      write(*,"(A30,2X,I6,1X,A9      )")' - Lower Band               : ', lowban,'  '
      write(*,"(A30,2X,I6,1X,A9      )")' - Bandwidth                : ', totban,'  '
      write(*,"(A30,2X,I6,1X,A9      )")' - Leading dimension of K   : ', ldAKban,' '
      
    end subroutine BandWidth
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine GlobalSystem(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      
      use sourceTerm
      
      implicit none
      
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      ! double precision,                       intent(in) :: k_y
      double precision, dimension(ndofn)        :: EMsource
      double precision, dimension(nne)          :: basis, xi_cor, yi_cor
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce
      double precision, dimension(nevab)        :: Fe
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)    :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_K, A_C, A_F
      
      allocate(A_K(ldAKban,ntotv), A_C(ldAKban,ntotv), A_F(ntotv, 1))
      
      !duda: Fe se declaró como a(n) y en la rutina assembleF como a(n,1)
      !      pero compila y ejecuta bien. ¿Poooor?
      
      A_K = 0.0
      A_F = 0.0
      A_C = 0.0
      
      do ielem = 1, nelem 
        !gather
        Ke = 0.0    !Esto es amate
        Fe = 0.0    !Fe(nevab)
        Ce = 0.0    !elemental capacity matrix (not used in static case)
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        !do-loop: compute element stiffness matrix Ke
        do igaus = 1, TotGp
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
          call DerivativesXY(igaus,Jinv,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,dN_dxy,HesXY)
          dvol  = detJ *  weigp(igaus,1)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          
          call source_term(ielem, basis, xi_cor, yi_cor, EMsource)
          ! print"(A12,I3,A26,3f15.5)",'for element',ielem,' the RHS contribution is: ', EMsource
          call Galerkin(hmaxi, dvol, basis, dN_dxy, EMsource, Ke, Ce, Fe) !amate lo llame Ke
          !call Galerkin(dvol, basis, dN_dxy, EMsource, Ke, Ce, Fe) !amate lo llame Ke
          !call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe, pertu,workm,resid)
          if(kstab.eq.6.or.kstab.eq.0)then
            !print*, kstab, ' sin estabi'
            continue
          else
            !print*, 'entra en stabi GlobalSystem'
            call TauMat(hmaxi,tauma)
            call Stabilization(hmaxi, dvol, basis, dN_dxy, HesXY, EMsource, tauma, Ke, Fe)
          endif
          
        end do
        !stop
        
        call Assemb_Glob_Mat(nodeIDmap, Ce, A_C)     !Assemble Global Conductivity Matrix K
        call Assemb_Glob_Mat(nodeIDmap, Ke, A_K)     !Assemble Global Conductivity Matrix K
        call Assemb_Glob_Vec(nodeIDmap, Fe, A_F)     !Assemble Global Source vector F
        
      end do
      
      !aaa = maxval(coord(1,:))*2**(-i_exp) 
      !if(aaa.ne.hmaxi) write(*,'(A)') '> > >Element size does not match'
      
    end subroutine GlobalSystem
    
    !subroutine AssembleK(K, ke, node_id_map, ndDOF)

    !  implicit none
    !  real(8), dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(in out)  :: K 
    !  !  Global Stiffnes matrix debe 
    !  !  llevar inout por que entra como variable (IN) 
    !  !  pero en esta funcion se modifica (out)
    !  real(8), dimension(2*nUne, 2*nUne), intent(in)   :: ke
    !  integer, dimension(nUne,1), intent(in)           :: node_id_map
    !  integer, intent(in)                              :: ndDOF 
    !  integer :: i, j, row_node, row, col_node, col !nodal Degrees of Freedom
    !  
    !  do i = 1, nUne
    !    row_node = node_id_map(i,1)
    !    row = ndDOF*row_node - (ndDOF-1)
    !    
    !    do j = 1, nUne
    !      col_node = node_id_map(j,1)
    !      col = ndDOF*col_node - (ndDOF-1)
    !      K(row:row+ndDOF-1, col:col+ndDOF-1) =  K(row:row+ndDOF-1, col:col+ndDOF-1) + &
    !      ke((i-1)*ndDOF+1:i*ndDOF,(j-1)*ndDOF+1:j*ndDOF)
    !    enddo
    !    
    !  enddo
    !  
    !  return
    !  
    !end subroutine AssembleK
   
    !subroutine ApplyBoundCond( NoBV, Fbcsvp, A_K, Sv )
    !  ! - - - - - - - - - - * * * * * * * * * * - - - - - - - 
    !  ! Set velocity (u) and preasure (p) boundary condition by penalty method
    !  ! - - - - - - - - - - * * * * * * * * * * - - - - - - - - - -
    !  implicit none
    !                      !Dof
    !  integer , dimension(NoBV,Dof), intent(in) :: Fbcsvp
    !  double precision, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(in out) :: A_K  !Global Stiffnes matrix
    !  double precision, dimension(2*n_nodes+n_pnodes, 1), intent(in out) :: Sv
    !  double precision :: param, coeff
    !  integer          :: preasure_row, NoBV, i, component, node_id, pnode_id
    !  
    !  !Esencialmente la siguiente instruccion hace: A_K(1*2-1,:) = A_K(1,:) Es decir, obtene el valor maximo de
    !  !la primera fila de la matriz global K (A_K). No le veo el caso pero lo dejamos asi.
    !  param = maxval(A_K(int(Fbcsvp(1,1))*2-1,:))
    !  coeff = abs(param) * 1.0E7

    !  print*, 'param', param
    !  print*, 'coeff', coeff
    !  
    !  preasure_row = 2*n_nodes

    !  do i =1, NoBV
    !    !print*, 'iteration', i
    !    node_id   = Fbcsvp(i,1) !se pone este int() pq la 1a y 2a col de Fbcsvp esta leida como integer pero 
    !    !print*, 'node_id', node_id
    !    component = Fbcsvp(i,2)!la matriz completa esta declarada como real en esta subroutina y en el main.
    !    !print*, 'component', component
    !    !print*, shape(Fbcsvp)
    !    !print*, Fbcsvp(i,:)
    !    if ( component .le. 2 ) then
    !      !print*, 'component of Boundary value', component
    !      !print*,'La pausa', 2*node_id-2+component,' ', 2*node_id-2 +component
    !      !read(*,*)
    !      A_K(2*node_id-2+component, 2*node_id-2 +component) = coeff
    !      Sv( 2*node_id-2+component, 1) = Fbcsvp(i,3)*coeff 
    !    else                                                     
    !      pnode_id = pnodes(node_id,2)
    !      !print*, 'pnode_id', pnode_id 
    !      !print*, 'preasure_row', preasure_row
    !      A_K(preasure_row+pnode_id, preasure_row + pnode_id) = coeff
    !      !print*, preasure_row+pnode_id, preasure_row + pnode_id
    !      Sv(preasure_row+pnode_id,1) = Fbcsvp(i,3)*coeff 
    !      !el tres es por que en la columna 3 esta el valor de la condicon de forntera
    !    end if
    !  end do

    !end subroutine ApplyBoundCond
   
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    
    !subroutine source_term(basis, xi_cor, yi_cor, source)
    !  
    !  implicit none     !interpolating X and Y
    !  
    !  !***********************************************************!
    !  !The source term is given by:                               !
    !  !                                                           !
    !  !              u = (fi*psi',-fi'*psi)                       !
    !  !                                                           !
    !  !   f = Lu       ;   where L is the diferential operator    !
    !  !                                                           !
    !  !***********************************************************!
    !  
    !  double precision,dimension(nne), intent(in) :: basis, xi_cor, yi_cor
    !  double precision :: dey_dydx, dex_dy2, dey_dx2, dex_dxdy
    !  double precision :: x, y
    !  integer :: ii
    !  double precision, dimension(ndofn), intent(out)  :: source
    !  
    !  !En cada llamada de esta funcion entrara basis con diferente gauss point
    !  
    !  x = 0.0
    !  y = 0.0
    !  
    !  !x= matmul(basis,xi_cor)
    !  do ii = 1, nne 
    !    x = x + basis(ii)*xi_cor(ii)
    !    y = y + basis(ii)*yi_cor(ii)
    !    !print"(3(1x,f10.7))",  basis(ii), yi_cor(ii), basis(ii)*yi_cor(ii)
    !  end do
    !  
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N1', basis(1),' xi1', xi_cor(1), 'x1', basis(1)*xi_cor(1)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N2', basis(2),' xi2', xi_cor(2), 'x2', basis(2)*xi_cor(2)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N3', basis(3),' xi3', xi_cor(3), 'x3', basis(3)*xi_cor(3)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N4', basis(4),' xi4', xi_cor(4), 'x4', basis(4)*xi_cor(4)
    !  !print"(A9,f10.5)",'x_gaus = ', x
    !  !
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N1', basis(1),' yi1', yi_cor(1), 'y1', basis(1)*yi_cor(1)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N2', basis(2),' yi2', yi_cor(2), 'y2', basis(2)*yi_cor(2)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N3', basis(3),' yi3', yi_cor(3), 'y3', basis(3)*yi_cor(3)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N4', basis(4),' yi4', yi_cor(4), 'y4', basis(4)*yi_cor(4)
    !  !print"(A9,f10.5)",'y_gaus = ', y
    !  
    !  !Derivatives in x-direction
    !  dey_dydx = 2.0-(4.0*y)
    !  dex_dy2  = 0.0 
    !  
    !  !Derivatives in y-direction
    !  dey_dx2  = 0.0
    !  dex_dxdy = (4.0*x)-2
    !  
    !  if(ndofn.eq.3)then
    !    source(1) = lambda*(dey_dydx - dex_dy2)
    !    source(2) = lambda*(-dey_dx2  + dex_dxdy)
    !    source(3) = 0.0
    !  elseif(ndofn.eq.2)then
    !    source(2) = lambda*(-dey_dx2  + dex_dxdy)
    !  else
    !    source(1) = lambda*(dey_dydx - dex_dy2)
    !  endif
    !  
    !end subroutine source_term
    
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine checkMKL(typpe, time,info)
      
      implicit none
      
      character(1), intent(in) :: typpe
      integer     , intent(in) :: info, time
      
      select case(typpe)
        case('f')
          if(info.ne.0)then
            print'(A34,I0)', '<<<Error in factorization at time: ', time
            call MKLfactoResult('dgbtrf',info) 
          endif
        case('s')
          if(info.ne.0)then
            print'(A48,I0)', '<<<Error in solving system of equation at time: ', time
            call MKLsolverResult('dgbtrs',info) 
          endif
      end select
      
    end subroutine checkMKL
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
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
        print*, ' ~ ~ ~ '
        print*, ' ~ ~ ~ Stopping the execution'
        print*, ' '
        stop
      endif
      print*, ' '

      101 format (A, 1x, I1, A)
      102 format (A, I4, A)
      103 format (A, I3, A, I3, A)
    end subroutine MKLfactoResult
    
    !subroutine MKLCondNumber(norm, ab, ipiv)
    !  implicit none
    !  
    !  double precision, dimension(ldAKban ,ntotv ), intent(in) :: ab
    !  integer, dimension(ntotv), intent(in) :: ipiv
    !  character(*), intent(in)              :: norm         !must be 1 or I
    !  !double precision          :: work(*)
    !  double precision, allocatable, dimension(:) :: work
    !  integer, dimension(ntotv) :: iwork
    !  double precision          :: rcond, val
    !  integer                   :: orderMatrix
    !  external                  :: xerbla, dgbcon
    !  integer                   :: info
    !  
    !  allocate(work(max(1,3*ntotv)))
    !  orderMatrix = ldAKban*ntotv   
    !  
    !  !val = dlangb(norm, orderMatrix, lowban, upban, ab, ldakban, work )
    !  val = 23550.9023
    !  call dgbcon( norm, orderMatrix, upban, lowban, ab, ldakban, ipiv, val, rcond, work, iwork, info )
    !  
    !  call xerbla('dgbcon',info)

    !  print*,'Reciprocal of Condition Number', rcond 
    !  print*,'Condition Number', 1./rcond 
    !  
    !  
    !  
    !end subroutine MKLCondNumber
    
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
      !print*,' '
    end subroutine MKLsolverResult
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine writeMatrix(Matrix, name1, Vector, name2)
      implicit none
      
      character(len=*), parameter    :: fileplace = "Res/"
      character(*) :: name1, name2
      integer :: i, j, mrow, ncol
      double precision, dimension(ldAKban ,ntotv ), intent(in) :: Matrix
      !double precision, dimension(ntotv), intent(in) :: Vector
      integer, dimension(ntotv), intent(in) :: Vector
     
     
      mrow = size(Matrix,1)
      ncol = size(Matrix,2)
      open(unit=10, file= fileplace//name1, ACTION="write", STATUS="replace")
     
      do i=1,mrow
        write(10, 100)( Matrix(i,j) ,j=1,ncol)
      end do
      close(10)
     
      open(unit=20, file= fileplace//name2, ACTION="write", STATUS="replace")
      do i=1,ncol
        !write(unit2, 100) Vector(i,1)
        write(20, 200) Vector(i)
      end do
      close(20)
      write(*,*) 'files: ', name1,' and ', name2, ' written succesfully on Res/'
      
      
      100 format (900E15.5)
      200 format (900I7)
      
    end subroutine writeMatrix 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine PostPro_EMfield(N,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,glob_potential,grad_sol)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Res/"
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ntotv), intent(in) :: glob_potential
      double precision, dimension(nevab)        :: elem_potential
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(nne)          :: xi_cor, yi_cor, basis
      double precision, dimension(nne,DimPr)    :: element_nodes
      double precision, dimension(nelem)        :: du_dx, du_dy 
      integer         , dimension(nne)          :: nodeIDmap
      double precision                          :: detJ, x, y
      integer                                   :: ielem, igaus, inode! ,idime, ipoin, ii, jj
      double precision, allocatable, dimension(:,:,:), intent(out) :: grad_sol
      
      allocate(grad_sol(nelem,totGp,DimPr))
      
      if(postpro.eq.1)then
        write(*,'(A)') 'Running Post Process'
      elseif(postpro.eq.2)then
        write(*,'(A)') ' '
        write(*,'(A)') ' -No Post Process'
        goto 110
      else
        write(*,'(A)') 'No postrocess option defined'
        goto 110
      end if
     
      open(unit=100, file= fileplace//'electric_field'//'.dat', ACTION="write", STATUS="replace")
      
      write(100,50) 'ielem','igaus','ex','ey','x', 'y'
      
      do ielem = 1,nelem
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        call gather(nodeIDmap, glob_potential, elem_potential) 
        
        !grad_sol = 0.0
       do igaus = 1, TotGp
          du_dx = 0.0 
          du_dy = 0.0 
          x = 0.0
          y = 0.0
          
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
          call DerivativesXY(igaus,Jinv,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,dN_dxy,HesXY)
          do inode = 1, nne
            basis(inode) = N(inode,igaus)
          end do
          
          !gradient computation (electric field)
          do inode = 1, nne
            du_dx(ielem) = du_dx(ielem) + dN_dxy(1,inode) * elem_potential(inode)
            du_dy(ielem) = du_dy(ielem) + dN_dxy(2,inode) * elem_potential(inode)
          end do
          
          grad_sol(ielem,igaus,1) = -du_dx(ielem)
          grad_sol(ielem,igaus,2) = -du_dy(ielem) 
          
          do inode = 1, nne 
            x = x + basis(inode)*xi_cor(inode)
            y = y + basis(inode)*yi_cor(inode)
            !print"(3(1x,f10.7))",  basis(ibase), yi_cor(ibase), basis(ibase)*yi_cor(ibase)
          end do
          
          write(100,60) ielem, igaus, x, y, grad_sol(ielem,igaus,1), grad_sol(ielem,igaus,2)
        end do
        !print*,' '
      end do
      
      close(100)
      110 continue 
      !do ielem =1,nelem
      !  do igaus = 1,TotGp
      !    write(*,'(2(1x,e15.5))') ( grad_sol(ielem, igaus, jj), jj=1,2 )
      !  end do
      !end do
      50 format(2x,A5,8x,A5,10x,A2,15x,A2,16x,A,15x,A)
      60 format(2(I5,1x),1x,2(1x,F15.5),3x,2(E15.5))
      
    end subroutine PostPro_EMfield
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine storeSpectrum( indexx, spectrum)
      
      implicit none
      
      character(len=*), parameter  :: path1 = "Pos/Plots/Spectrums/"
      character(len=8)                                :: id_file
      ! double precision, dimension(ntotv, 1),intent(in) :: u_pre
      double precision, dimension(ntotv, 1),intent(in) :: spectrum
      !double precision, dimension(ntotv, t_steps+1)            :: spectrum
      integer                              ,intent(in) :: indexx
      integer                                          :: ipoin, ii, jj
      
      id_file = files_ky(indexx)
      
      ! do ii = 1, S_ldSol
      !   store_Spec(ii,time+1) = u_pre(ii,1) 
      ! end do
      
      ! Se guardan las componentes del campo transformado para cada tiempo:
      !   t1  t2  t3
      !----- ----  --- 
      !  Êx1  Êx1  Êx1
      !  Êy1  Êy1  Êy1
      !  Êz1  Êz1  Êz1
      !  Êx2  Êx2  Êx2
      !  Êy2  Êy2  Êy2
      !  Êz2  Êz2  Êz2
      !  Êx3  Êx3  Êx3
      !  Êy3  Êy3  Êy3
      !  Êz3  Êz3  Êz3
     
      open(unit=200, file=path1//'spectrumE_field_'//id_file, ACTION="write", STATUS="replace")
      write(200,"(A)") '#Transformed Voltage (total variables, time steps)'
      do ii = 1, ntotv
        write(200,917) (spectrum(ii,jj), jj = 1, t_steps+1)
      end do
      close(200)
      
      
      917 format(999(3x,E15.5)) !format for store transformed E-field
      
    end subroutine storeSpectrum
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine invDFT
      
      use, intrinsic                                  :: iso_c_binding
      implicit none
      
      external                                        :: fdate 
      character(len=*), parameter                     :: path1 = "Pos/Plots/Spectrums/"
      character(len=*), parameter                     :: path2 = "Pos/Plots/"
      double precision, parameter                     :: pi = 4*atan(1.d0)
      !Hacer una rutina que contenga los casos posibles de cantidad de numeros de onda a emplear
      !Que sean 10, 15, 20 
      character(len=8)                                :: id_file
      character(len=24)                               :: date
      character(len=180)                              :: msg
      double precision                                :: ky, delta_ky, dt
      double precision, allocatable, dimension(:,:)   :: E_xyzt, E_3D
      double precision, dimension(t_steps+1)          :: dummy
      double precision, allocatable, dimension(:,:,:) :: E_hat_ky
      integer                                         :: nt,ii,jj,kk,ll,stat,iwn,totv,idDoF
      
      ! Declara la función sleep de C
      interface
        subroutine usleep(useconds) bind(c, name="usleep")
          import :: c_int
          integer(c_int), value :: useconds
        end subroutine usleep
      end interface
      print*, ' ' 
      print'(A)', " !=============== Performing the Inverse Fourier Transform =============! "
      ! Pausa durante 1 segundo (1,000,000 microsegundos)
      call usleep(1000000)
      !Solo si el numero de onda actual es igual al numero total de numeros de onda, se 
      !ejecutara el siguiente algoritmo que realiza la Transformada Inversa de Fourier 
      !en los siguiente pasos:
      !       1) Leer los resultados 2D para cada numero de onda y guardarlos en una matriz de matrices Ê(ky,
      !       2) 
      !       3)
      !       4)
      
      ! Se guardan las componentes del campo transformado para cada tiempo:
      !   t1  t2  t3
      !----- ----  --- 
      !  Êx1  Êx1  Êx1
      !  Êy1  Êy1  Êy1
      !  Êz1  Êz1  Êz1
      !  Êx2  Êx2  Êx2
      !  Êy2  Êy2  Êy2
      !  Êz2  Êz2  Êz2
      !  Êx3  Êx3  Êx3
      !  Êy3  Êy3  Êy3
      !  Êz3  Êz3  Êz3
      !Como estoy dividiendo los problemas, deberia  agregar un check para
      !ver que todos los archivos .dat existen, si no existen, esperar unn tiempo 
      !y volver a revisar, y si existen, ejecutar la transformada
      
      allocate( E_hat_ky(tot_ky,t_steps+1,ntotv))
      allocate( E_xyzt(ntotv,t_steps+1))
      allocate( E_3D(ntotv,1))
      
      idDoF = (ndofn-1) != 3-1=2 = vecor or 1-1 = 0 = scalar problem 
      reading_2D_results: do iwn = 1,tot_ky
        id_file = files_ky(iwn)
        open(5,file=path1//'spectrumE_field_'//id_file, status='old',action='read',IOSTAT=stat, IOMSG=msg)
        IF ( stat /= 0 )then
          print'(53A,I0)', 'ioStat for OPPENING 2D electric field spectrum files ', stat
          print*, msg
        end if
        read(5,*,iostat=stat,iomsg=msg) !se salta los encabezados
        do jj = 1,ntotv
          read(5,*,iostat=stat,iomsg=msg) (E_hat_ky(iwn,nt,jj), nt =1,t_steps+1 )
          IF ( stat /= 0 )then
            print'(53A,I0)', 'ioStat for READING 2D electric field spectrum files ', stat
            print*, msg
          end if
        end do
        close(5) 
      end do reading_2D_results
      
      !!Este print es para imprimir en pantalla las matrices que conforman la matriz de resultados 2D
      !do iwn = 1, tot_ky
      !  print'(A3,I0)','Ky',iwn
      !  do ii = 1,t_steps+1
      !    print'(99(E11.3))', (E_hat_ky(iwn,ii,jj), jj=1,ntotv)
      !  end do
      !end do
     
      call fdate(date) 
      open(unit=300, file=path2//shape_spec_file, ACTION="write", STATUS="replace")
      write(300,"(A,1x,A)") '%2D-CDR-EM Simulation: Ê-field vs ky   ', date
      if(ProbType=='TIME')then
        write(300,"(A13, I0)") '%At the time ', t_steps-2 
        ky_loop1: do ll = 1,tot_ky
          ky = WaveNumbers(ll)
          ! Ex_hat = (ndofn-1) = 3-1=2
          ! Ey_hat = (ndofn-2) = 3-2=1
          ! Ez_hat = (ndofn-3) = 3-3=0
          ! write(300,904) ll, ky, Ex_hat, Ey_hat, Ez_hat 
          write(300,904) ll, ky, (E_hat_ky(ll,t_steps-2,ndofn*receivers(jj)-2), jj=1,nodalRec),&
            &(E_hat_ky(ll,t_steps-2,ndofn*receivers(jj)-1), jj=1,nodalRec),&
            &(E_hat_ky(ll,t_steps-2,ndofn*receivers(jj)-0), jj=1,nodalRec)
          !Ex_fieldi(time) = (Sol_T(1,(ndofn*receivers(ipoin)+1)), ipoin=1,nodalRec)
        end do ky_loop1
      else
        !This is for a static and scalar problem
        if(ProbType=='STAT')t_steps=3 
        write(300,"(A)") '% No           ky                    Êx'
        write(300,"(A12, I0)") '%Receivers: ', nodalRec
        write(300,903) (coord(1,receivers(jj)), jj=1,nodalRec) 
        ky_loop2: do ll = 1,tot_ky
          ky = WaveNumbers(ll)
          write(300,904) ky, (E_hat_ky(ll,t_steps-2,ndofn*receivers(jj)), jj=1,nodalRec)
          !Ex_fieldi(time) = (Sol_T(1,(ndofn*receivers(ipoin)+1)), ipoin=1,nodalRec)
        end do ky_loop2
        if(ProbType=='STAT')t_steps=0
      endif
      close(300)
     
     
      
      ! = = = = = = = = = = = = = PERFORMING INVERSE FOURIER TRANSFORM = = = = = = = = = = = = = = 
      !                            1  __∞_
      !             Ê(x,y,z,t) = ___  \     E(x,ky,z,t) * exp(-i*ky*y) dky   !How to choose dky?
      !                           2π  /___
      !                               ky=0
      !
      print*,'Esto es t_steps antes de comenzar la TDF', t_steps
      delta_ky = (ky_max-ky_min)/ tot_ky
      E_xyzt = 0.0
      time_loop: do ii =1,t_steps+1  
        nodes_loop: do jj =1,nnodes
          totv = jj*ndofn
          DoF_loop: do kk = idDoF,0,-1
            ky_loop3: do ll =1, tot_ky
              ky=WaveNumbers(ll)
              E_xyzt(totv-kk,ii) = E_xyzt(totv-kk,ii) + E_hat_ky(ll,ii,totv-kk)*exp(-cmplx(0,1)*y_iFT*ky) * delta_ky
              ! E_xyzt(totv-kk,ii) = E_xyzt(totv-kk,ii) + E_hat_ky(ll,ii,totv-kk)*cos(ky*y_iFT) * delta_ky
            end do ky_loop3
          end do DoF_loop
        end do nodes_loop
      end do time_loop
      E_xyzt = (1.0/2.0*pi) * E_xyzt   !Transformada Exponencial
      ! E_xyzt = (1.0/pi) * E_xyzt       !Transformada coseno Moghaddam et al. 1991
      ! E_xyzt = (2./pi) * E_xyzt          !Transformada coseno Queralt et al. 1989 
      
      ! print*,' '
      ! print*,' '
      ! print'(A)','E_xyzt'
      ! do jj=1,ntotv
      !   print'(99(E11.3))', (E_xyzt(jj,ii), ii = 1,t_steps+1)
      ! end do

      TwoHalf = 'N'
      File_Nodal_Vals = File_3DNodal_Vals
      E_3D = 0.0
      if(ProbType == 'TIME')then
        dt = time_ini
        time_loop2: do ii =1,t_steps+1  
          dt = dt + delta_t!,time_fin,delta_t
          nodes_loop2: do jj =1,nnodes
            totv = jj*ndofn
            DoF_loop2: do kk = idDoF,0,-1
              E_3D(totv-kk,1) = E_xyzt(totv-kk,ii)
              end do DoF_loop2
          end do nodes_loop2
          call GID_PostProcess(1, E_3D, 'res', ii-1, dt, time_fin, dummy)
        end do time_loop2
        call GID_PostProcess(1, E_3D, 'msh', 0, dt, time_fin, dummy)
      else
        nodes_loop3: do jj =1,nnodes
        totv = jj*ndofn
          DoF_loop3: do kk = idDoF,0,-1
            E_3D(totv-kk,1) = E_xyzt(totv-kk,1)
          end do DoF_loop3
        end do nodes_loop3
        call GID_results(E_3D) 
      endif
      open(unit=10, file=path2//"spatial_profile.dat", ACTION="write", STATUS="replace")
      do jj =1,nodalRec
          ! print*, ndofn*receivers(jj)
          write(10,902) coord(1,ndofn*receivers(jj)), E_3D(ndofn*receivers(jj),1)
          ! El error esta aqui por que no se imprime lo que debe ser el valor de E3D
      end do
      close(10)
        
        
      !  NN = 19
      !  w0 = 2.0*pi/NN
      !  reAll = 0.0
      !  imagi = 0.
      !  ! Esto deberia ser el campo electrico transformado que debo leer desde un inputfile
      !  do i = 1, size(t)
      !    arg = 2*pi*12.5*t(i)
      !    f_t(i)= cos(arg)
      !  end do
      !  do k =0, NN-1
      !    do n = 0, NN-1
      !      angle = k*w0*n
      !      reall(k) = reall(k) + f_t(n)*cos(angle) / NN
      !      imagi(k) = imagi(k) - f_t(n)*sin(angle) / NN
      !      fhat(k) = fhat(k) + f_t(n)*exp(-cmplx(0,1)*angle)
      !    end do
      !    write(*,'(I5,2x,5(F7.3,1x))') k, f_t(k), reall(k), imagi(k), fhat(k)
      !  end do
      902 format(5(E18.6))
      903 format(99f12.5)    !format to print the profile file
      904 format(99(e18.6))    !format to print the profile file
      ! 904 format(I5,2x,e15.6,5x,99(e17.6))    !format to print the profile file
    end subroutine invDFT
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine GID_results(solution, grad_sol)
    ! subroutine GID_results(solution, grad_sol,profile)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Pos/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      double precision, dimension(nelem,totGp,DimPr), intent(in), optional :: grad_sol
      ! integer                                       , intent(in), optional :: profile
      character(len=10)                       :: ext1, ext2
      character(len=15)                       :: Elem_Type
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcor, ycor
      integer                                 :: ipoin, ii, ielem, inode, jj, igaus, icomp, RESconma
      
      solution_T = transpose(solution)
      xcor  = spread(coord(1,:),dim = 1, ncopies= 1)
      ycor  = spread(coord(2,:),dim = 1, ncopies= 1)
      if(TwoHalf=='Y')File_Nodal_Vals=File_Nodal_Vals_ky//ky_id
     
      print*, ' '
      print*, '!============== Output files ==================!'
      
      ext1 = ".post.msh"
      ext2 = ".post.res"
      if(ElemType.eq.'QUAD')then
        Elem_Type = 'Quadrilateral'
      elseif(ElemType.eq.'TRIA')then
        Elem_Type = 'Triangle'
      endif
      open(unit=555, file= fileplace//File_Nodal_Vals//ext1, ACTION="write", STATUS="replace")
      
      write(555,902) 'MESH', '"Domain"', 'dimension', DimPr, 'ElemType', Elem_Type, 'Nnode', nne
      write(555,"(A)") '#2D Convection-Diffusion-Reaction'
      write(555,900) '#Element tipe: ', ElemType,'/',ElemType
      
      
      write(555,"(A)")'Coordinates'
      write(555,"(A)") '#   No        X           Y'
      do ipoin = 1, nnodes
        ! write(555,906) ipoin, xcor(1,ipoin),  ycor(1,ipoin)
        write(555,906) ipoin, xcor(1,ipoin), y_iFT,  ycor(1,ipoin)
      end do
      write(555,"(A)") 'End Coordinates'
      
      
      
      write(555,"(A)") 'Elements'
      do ielem=1,nelem
        write(555,908) ielem,(lnods(ielem,inode),inode=1,nne)
      end do
      write(555,"(A)") 'End Elements'
      
      
      close(555)
      write(*,"(A7,A21,A28)") ' -File ',File_Nodal_Vals//'.post.msh','written succesfully in Pos/'
      
      
      open(unit=555, file= fileplace//File_Nodal_Vals//ext2, ACTION="write", STATUS="replace")
      write(555,"(A)") 'GiD Post Results File 1.0'
      write(555,"(A)") '#2D Convection-Diffusion-Reaction'
      
      ! se escribe el res de las componentes de la velocidad
      select case(ndofn)
        case(1)
          print*, 'Case 1 es un PROBLEMA ESCALAR'
          write(555,"(A)") 'Result "phi" "Electric Potential" 0 Scalar OnNodes'
          write(555,"(A)") 'ComponentNames "" '
          write(555,"(A)") 'Values'
          write(555,*) '#',   'No    ','             ux '
          !  se escribe el res para el caso escalar de un grado de libertad
          write(555,914)
          do ipoin = 1, nnodes
            write(555,916) ipoin, solution(ipoin, 1)
          end do
          write(555,"(A)") 'End Values'
          
        case(2)
          write(555,"(A)") 'Result "EM Field" "EM field" 0 Vectsur OnNodes'
          write(555,"(A)") 'ComponentNames "ex" "ey" "--" "" '
          write(555,"(A)") 'Values'
          write(555,*) '#',   'No    ','             ex ','               ey '
          do ipoin = 1, nnodes
            write(555,919) ipoin, solution_T(1, ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
          end do
          write(555,"(A)") 'End Values'
          
        case(3)
          RESconma = sum(sum(conma, DIM=1))
          if(RESconma == 0)then
            write(555,"(A)") 'Result "EM field" "EM  field" 0 Vector OnNodes'
            write(555,"(A)") 'ComponentNames "Ex" "Ey" "Ez" "" '
            write(555,"(A)") 'Values'
            write(555,"(A)") '#No                 ex                ey                 ez'
            do ipoin = 1, nnodes
              write(555,918) ipoin,&
              &       solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
            end do
            write(555,"(A)") 'End Values'
            
          else
            write(555,"(A)") 'Result "EM field" "EM  field" 0 Vector OnNodes'
            write(555,"(A)") 'ComponentNames "ex" "ey" "--" "" '
            write(555,"(A)") 'Values'
            write(555,*) '#',   'No    ','             ex ','               ey'
            do ipoin = 1, nnodes
              write(555,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1)
            end do
            write(555,"(A)") 'End Values'
            write(555,"(A)") 'Result "Multiplier" "Multiplier" 0 Scalar OnNodes'
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
            write(*,"(A7,A21,A28)") ' -File ',File_Nodal_Vals//'.post.res','written succesfully in Pos/'
          endif
      end select

      ! if(present(profile)then
      !   open(unit=10, file=path2//"spatial_profile.dat", ACTION="write", STATUS="replace")
      !   do ipoin =1,nodalRec
      !       write(10,903) coord(1,receivers(ipoin)), solution_T(1,receivers(ipoin))
      !   end do
      !   close(10)
      ! else
      !   write(*,'(A)') 'No spatial profile required'
      ! endif
      
      if (present(grad_sol) )then
        if(postpro.eq.1)then
          write(*,'(A)') 'Writing Post Process.....'
          write(555,"(A,A)") 'GaussPoints "GP_1" ElemType ', Elem_Type 
          write(555,"(A24,I1)") 'Number Of Gauss Points: ', totGp
          write(555,"(A,A)") 'Natural Coordinates: ', "Given"
          do igaus = 1, totGp
          write(555,901) (ngaus(igaus,jj), jj=1,DimPr)
          end do
          write(555,"(A)") 'End GaussPoints'

          write(555,"(A)") 'Result "Electric field" "ANALYSIS" 0 Vector OnGaussPoints "GP_1" '
          write(555,"(A)") 'ComponentNames "Ex" "Ey" '
          write(555,"(A)") 'Values'
          do ielem = 1, nelem
          write(555,'(I0)',Advance='NO') ielem
          do igaus = 1, totGp
          write(555,903) (grad_sol(ielem,igaus,icomp), icomp=1,2)
          end do
          end do
          write(555,"(A)") 'End Values'
        elseif(postpro.eq.2)then
          !write(*,'(A)') 'None Post Process'
          110 close(555)
        else
          write(*,'(A)') 'No postrocess option defined'
        end if
      endif
      
      
      900 format(A15, A13, A1, A13)
      901 format(2(f10.5))
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      903 format(2(E15.5))
      906 format(I7,4x,3(f15.5,3x)) !format for msh
      908 format(10(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,E15.5)  !format for scalar case
      918 format(I7,3(4x,E15.5)) !format for res velocity
      919 format(I7,2(4x,E15.5)) !format for res velocity
      
    end subroutine GID_results 
   
    
    
    
    
    !- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- -
    ! Agregamos las rutinas de tiempo del commit 625fbee: May2 2022
    
    subroutine GID_PostProcess(id,solution, activity, time, timeStep, time_final, Ex_field)
      
      ! use E0field
      
      implicit none
      external                                             :: fdate 
      
      character(len=*), parameter  :: fileplace = "Pos/"
      character(len=*), parameter  :: fileplace2 = "Res/FEM_TEM/"
      character(len=*), parameter  :: fileplace3 = "Exact_Sol_TEM/2D_DoubleLine_WholeSpace/"
      character(len=24)                                :: date
      double precision, dimension(ntotv, 1),intent(in) :: solution
      character(*)                         ,intent(in) :: activity
      integer                              ,intent(in) :: time
      double precision                     ,intent(in) :: timeStep, time_final
      character(len=10)                                :: ext1, ext2
      character(len=4)                                 :: ext3
      character(len=15)                                :: Elem_Type
      double precision, dimension(1, ntotv)            :: Sol_T
      double precision                                 :: Ez_r(ntotv), tEz
      double precision, dimension(1,nnodes)            :: xcor, ycor
      integer                                          :: ipoin, ii, ielem, inode, time2,id, RESconma
      double precision, dimension(t_steps+1), intent(out) :: Ex_field
      double precision :: x_profile
      
      !double precision :: delta_t, timeStep2
      
      
      Sol_T = transpose(solution)
      xcor  = spread(coord(1,:),dim = 1, ncopies= 1)
      ycor  = spread(coord(2,:),dim = 1, ncopies= 1)
      if(TwoHalf=='Y')File_Nodal_Vals=File_Nodal_Vals_ky//ky_id
      
      !delta_t 1e-3 ( time_fin - time_ini ) / (t_steps + 1.0)
      
      if(id.eq.1)then 
        ext1 = ".post.msh"
        ext2 = ".post.res"
        ext3 = ".dat"
      elseif(id.eq.2)then
        ext1 = ".ppost.msh"
        ext2 = ".ppost.res"
      endif
      if(ElemType.eq.'QUAD')then
        Elem_Type = 'Quadrilateral'
      elseif(ElemType.eq.'TRIA')then
        Elem_Type = 'Triangle'
      endif
      !open(unit=555, file= fileplace//File_Nodal_Vals//ext1, ACTION="write", STATUS="replace")
      
      if(activity == "msh")then !quitar este if y acomodar el numero de unidad
        open(unit=100, file= fileplace//File_Nodal_Vals//ext1, ACTION="write", STATUS="replace")
        
        write(100,902) 'MESH', '"Domain"', 'dimension', DimPr, 'ElemType', Elem_Type, 'Nnode', nne
        write(100,"(A)") '#2D Convection-Diffusion-Reaction'
        write(100,900) '#Element tipe: ', ElemType,'/',ElemType
        write(100,"(A)")'Coordinates'
        write(100,"(A)") '#   No        X           Y'
        do ipoin = 1, nnodes
          write(100,906) ipoin, xcor(1,ipoin), y_iFT ,ycor(1,ipoin)
        end do
       
        write(100,"(A)") 'End Coordinates'
        write(100,"(A)") 'Elements'
        do ielem=1,nelem
          write(100,908) ielem,(lnods(ielem,inode),inode=1,nne)
        end do
        write(100,"(A)") 'End Elements'
        close(100)
        print"(A11,A21,A30)", ' Mesh file ',File_Nodal_Vals//ext1, 'written succesfully in Pos/ '
       
      elseif(activity == "res")then
       
        if(time == 0)then
          open(unit=200, file= fileplace//File_Nodal_Vals//ext2, ACTION="write", STATUS="replace")
          write(200,"(A)") 'GiD Post Results File 1.0'
          write(200,"(A)") '#2D Convection-Diffusion-Reaction'
        else
          continue
        endif
        open(unit=200,file=fileplace//File_Nodal_Vals//ext2,ACTION="write",STATUS="old",position="append")
        
        ! se escribe el res de las componentes de la velocidad
        select case(ndofn)
          case(1)
            write(200,"(A21, I0, A)") 'Result "E" "E-field" ', time,' Scalar OnNodes'
            write(200,"(A)") 'ComponentNames "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ex '
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
            write(200,"(A21, I0, A)") 'Result "E" "E-field" ', time,' Vector OnNodes'
            write(200,"(A)") 'ComponentNames "u" "v" "--" "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ex ','               ey '
            do ipoin = 1, nnodes
              write(200,918) ipoin, Sol_T(1, ndofn*ipoin-1), Sol_T(1,ndofn*ipoin)
            end do
            write(200,"(A)") 'End Values'
          case(3)
            RESconma = sum(sum(conma, DIM=1))
            if(RESconma == 0)then
              write(200,"(A21, I0, A)") 'Result "E" "E-field" ', time,' Vector OnNodes'
              write(200,"(A)") 'ComponentNames "Ex" "Ey" "Ez" "" '
              write(200,"(A)") 'Values'
              write(200,*) '#No      ','             ex ','               ey','               ez'
              do ipoin = 1, nnodes
                write(200,919) ipoin,&
                &       Sol_T(1, ndofn*ipoin-2), Sol_T(1,ndofn*ipoin-1), Sol_T(1,ndofn*ipoin)
              end do
              write(200,"(A)") 'End Values'
              
            else
              write(200,"(A21, I0, A)") 'Result "E" "E-field" ', time,' Vector OnNodes'
              write(200,"(A)") 'ComponentNames "Ex" "Ey" "Ez" "" '
              write(200,"(A)") 'Values'
              write(200,*) '#',   'No    ','             ex ','               ey'
           !   do ipoin = 1, nnodes
           !     write(200,919) ipoin, Sol_T(1, ndofn*ipoin-2), Sol_T(1,ndofn*ipoin-1), Sol_T(1,ndofn*ipoin)
           !   end do
              do ipoin = 1, nnodes
                write(200,919) ipoin, Sol_T(1, ndofn*ipoin-2), Sol_T(1,ndofn*ipoin-1)
              end do
              write(200,"(A)") 'End Values'
              write(200,"(A24, I0, A)") 'Result "P" "Multiplier" ', time,' Scalar OnNodes'
              write(200,"(A)") 'ComponentNames "" '
              write(200,"(A)") 'Values'
              write(200,*) '#',   'No    ','     P '
              !  se escribe el res para el caso escalar de un grado de libertad
              write(200,914)
              ii=1
              do ipoin = 3, nnodes*3,3
                write(200,916) ii, Sol_T(1,ipoin)
                ii=ii+1
              end do
              write(200,"(A)") 'End Values'
            end if
        end select
        
        
        
      elseif(activity == "profile")then
        Ex_field = 0.0 
        if(time == 0)then
          call fdate(date) 
          open(unit=300, file=fileplace2//profile_name//ext3, ACTION="write", STATUS="replace")
          write(300,"(A,1x,A)") '%2dCDREM simulation: E-field vs time   ', date
          write(300,"(A)") ' '
          write(300,"(A)") '% - - - - - - - Component ex'
          if(time == 0) write(300,"(A)") '% time     timeStep        receiver1'
        else
          continue
        endif
        open(unit=300, file=fileplace2//profile_name//ext3, ACTION="write",STATUS="old",position="append")
        
        write(300,904) time, timeStep, (Sol_T(1,(ndofn*receivers(ipoin))), ipoin=1,nodalRec)
        !Ex_fieldi(time) = (Sol_T(1,(ndofn*receivers(ipoin)+1)), ipoin=1,nodalRec)
        !Creo que aqui en lugar del +1 debio ser -2 para que el grado de libertad del receptor
        !corresponda a la componente x

        !if( time == t_steps+1 ) then
        !  write(300,"(A)") ' '
        !  write(300,"(A)") '% - - - - - - - Component ey'
        !  if(time == 0) write(300,"(A)") '% time       timeStep          receiver1'
        !  timeStep2 = 0.0
        !  do time2 = 1, t_steps+1
        !    timeStep2 = timeStep2 + delta_t
        !    write(300,904) time2, timeStep2, (Sol_T(1,(ndofn*receivers(ipoin)+2)), ipoin=1,nodalRec)
        !  end do
        !endif
        
      elseif(activity == "spatial")then
        !id_poin = 113
        !open(unit=10, file= fileplace3//"Id_spatial_profile.dat", ACTION="write", STATUS="replace")
        !do ipoin =1,nnodes
        !  if(ycor(1,ipoin).eq.0.0)then
        !    id_poin = id_poin+1
        !    write(10,906) ipoin, xcor(1,ipoin)
        !  else
        !    continue
        !  endif
        !end do
        !close(10)



        !if(time == 1)then
        !  open(unit=6, file=fileplace3//"test_for_commit.dat", STATUS="replace", ACTION="write")
        !else
        !  open(unit=6, file=fileplace3//"test_for_commit.dat", ACTION="write", STATUS="old", position="append")
        !endif
        
        !open(unit=5, file= fileplace3//"Id_spatial_profile.dat", status='old', action='read')
        !!allocate(efile_profile(id))
        
        !tEz = timeStep 
        !call Efield_WholeSpace(time, tEz, Ez_r)
        
        !write(6,'(A7,I0,A,e10.3,A)') ' #t',time,'=',timeStep
        !write(6,'(A)') "index xcor FEM Exact"
        !if(ndofn.eq.1)then
        !  do ii = 1,id_poin                !id_poin viene del modulo E0field como variable global 
        !    read(5,*) ipoin, x_profile
        !    write(6,918) ipoin, x_profile, Sol_T(1, ndofn*ipoin), Ez_r(ndofn*ipoin)
        !  end do
        !else
        !  do ii = 1,id_poin
        !    read(5,*) ipoin, x_profile
        !    write(6,918) ipoin, x_profile, Sol_T(1, ndofn*ipoin-2), Ez_r(ndofn*ipoin)
        !  end do
        !endif
        !write(6,*)' '
        !write(6,*)' '
        !close(5)
       !close(6)
        
      else
        write(*,"(A)") ' < < Error > > Postprocess activity must be "msh", "res" or "profile" non ', activity
        close(200)
        close(300)
        stop
      end if
      
      close(200)        !Estos close van aqui y el position append permite abrir una y otra vez el
      close(300)        !mismo archivo escribiendo a continuacion de donde se quedo el archivo anterior
      
      !la siguiente instruccion debe usarse con timeStep no con time pero solo es para avanzar
      !if(time == t_steps+1.and.activity.eq."profile") then
      if((time == t_steps).and.(activity.ne."profile")) then
        print*, ' '
        print"(1x, A21,A30)", File_Nodal_Vals//'.post.res', 'written succesfully in Pos/ '
        print*, ' '
      endif
      
      !if(simul.eq.6.and.ndofn.1)then
      !  do ipoin = 1, nnodes
      !    write(555,906) ipoin, xcor(1,ipoin), ycor(1,ipoin)
      !  end do
      !else
      !  write(*,'A')
      !
      
      900 format(A15, A13, A1, A13)
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      904 format(I5,2x,e15.6,2x,99(e15.6))    !format to print the profile file
      906 format(I7,3(f15.5,3x)) !format for msh
      908 format(9(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,E12.5)  !format for scalar case
      918 format(I7,3x,5(E14.6,3x)) !format for res velocity
      919 format(I7,3(3x,E15.5)) !format for res velocity
      
    end subroutine GID_PostProcess
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine infoTime(time)
      
      implicit none 
      integer, intent(in) :: time
      double precision    :: check1, check2, check3, check4, check5
      double precision    :: check6, check7, check8, check9, check10
      
      check1=((t_steps/20)*100)/t_steps
      check2=((t_steps/8)*100)/t_steps
      check3=((t_steps/3)*100)/t_steps
      check4=((t_steps/2.2)*100)/t_steps
      check5=((t_steps/1.9)*100)/t_steps
      check6=((t_steps/1.6)*100)/t_steps
      check7=((t_steps/1.4)*100)/t_steps
      check8=((t_steps/1.2)*100)/t_steps
      check9=((t_steps/1.06)*100)/t_steps
      check10=((t_steps/1.01)*100)/t_steps
      check1= ceiling(check1)
      check2= ceiling(check2)
      check3= ceiling(check3)
      check4= ceiling(check4)
      check5= ceiling(check5)
      check6= ceiling(check6)
      check7= ceiling(check7)
      check8= ceiling(check8)
      check9= ceiling(check9)
      check10= ceiling(check10)

      ! print*,'1',check1
      ! print*,'2',check2
      ! print*,'3',check3
      ! print*,'4',check4
      ! print*,'5',check5
      ! print*,'6',check6
      ! print*,'7',check7
      ! print*,'8',check8
      ! print*,'9',check9
      ! print*,'10',check10
      
      
      if((time==floor(t_steps*check1*0.01)))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check1), '%'
      elseif(time==floor(t_steps*check2*0.01))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check2), '%'
      elseif(time==floor(t_steps*check3*0.01))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check3), '%'
      elseif(time==floor(t_steps*check4*0.01))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check4), '%'
      ! elseif(time==floor(t_steps*check5))then
        ! print*, ' -Completed5',check5, '%'
      elseif(time==floor(t_steps*check6*0.01))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check6), '%'
      ! elseif(time==floor(t_steps*check7))then
        ! print*, ' -Completed',check7, '%'
      elseif(time==floor(t_steps*check8*0.01))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check8), '%'
      elseif(time==floor(t_steps*check9*0.01))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check9), '%'
      elseif((time==floor(t_steps*check10*0.01)))then
        print'(A12,i0,1x,A)', ' -Completed ',int(check10), '%'
      endif
    end subroutine infoTime
    
    
    !subroutine GlobalSystem_Time(N,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldSol,delta_t,ugl_pre,A_F)
    subroutine prevTime(N,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldSol,ugl_pre,A_M)
      !        BDF1_u_prev
      
      use sourceTerm
      
      implicit none
      
      double precision, allocatable, dimension(:,:), intent(in out) :: ugl_pre
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      ! double precision,                       intent(in) :: k_y
      !double precision, dimension(3,nne), intent(in)     :: Hesxieta
      integer                               , intent(in) :: S_ldSol
      double precision, dimension(ndofn)        :: EMsource
      double precision, dimension(nne)          :: basis, xi_cor, yi_cor
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce, rhs_CN
      double precision, dimension(nevab)        :: Fe, Mu_time, ue_pre, time_cont
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)    :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ!, delta_t
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_M
      
      allocate( A_M(ntotv, 1) )
      
      A_M = 0.0
      do ielem = 1, nelem 
       
        Ke = 0.0; Ce = 0.0; Fe = 0.0; Mu_time = 0.0 
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        !gather
        call gather(nodeIDmap, ugl_pre, ue_pre)
        time_cont = ue_pre
        
        !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
        do igaus = 1, TotGp
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
          dvol = detJ *  weigp(igaus,1)
          !call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)
          call DerivativesXY(igaus,Jinv,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,dN_dxy,HesXY)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          
          call source_term(ielem, basis, xi_cor, yi_cor, EMsource)
          !En este Galerkin deberia quitar la construccion de la matriz Ke para evitar Calculamos
          !inecesarios pues solo se calcula Ce y Fe y Ke ya no, REVISA RUTINA Galerkin_prevTime 
          !y ver si esa se puede implementar tal cual aqui en lugar de Galerkin completo
          call Galerkin(hmaxi, dvol, basis, dN_dxy, EMsource, Ke, Ce, Fe) !amate lo llame Ke
          !call Galerkin(dvol, basis, dN_dxy, Ke, Ce, Fe) 
          !!call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe, pertu,workm,resid)
          
          if(kstab.eq.6.or.kstab.eq.0)then
            !print*, kstab, ' sin estabi'
            continue
          else
            !print*, 'entra en stabi Prev Time'
            call TauMat(hmaxi,tauma)
            call Stabilization(hmaxi, dvol, basis, dN_dxy, HesXY, EMsource, tauma, Ke, Fe)
          endif
          !estas multiplicaciones deberias ser globales pero por la matriz en banda se deben 
          !hacer locales, es lo mismo.
          select case(theta)
            case(2)
              Mu_time = matmul(Ce,time_cont) 
            case(3) !Este caso se debe quitar de aqi. Esta rutina no se llama en Crank-Nicholson
              rhs_CN  = (1.0/delta_t)*Ce - 0.5*Ke
              Mu_time = 0.5*Fe + matmul(rhs_CN,ue_pre)
          endselect
          
        end do
        
        !call Assemb_Glob_Mat(nodeIDmap, Ce, A_C)      !Assemble Global Capacity Matrix C 
        call Assemb_Glob_Vec(nodeIDmap, Mu_time, A_M) !Assemble Global Source vector F
        
      end do
      
      
    end subroutine prevTime




    subroutine Galerkin_prevTime(dvol, basis, Ce, Fe)
      
      implicit none
      
      double precision, intent(in) :: basis(nne)
      double precision, intent(in) :: dvol
      integer :: inode, idofn, ievab, jevab, jnode, jdofn
      double precision ::  cpcty
      double precision, intent(out) :: Fe(nevab), Ce(nevab,nevab)
      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              cpcty = basis(inode) * basis(jnode)
              Ce(ievab,jevab) = Ce(ievab,jevab) + cpcty * dvol                                 !element Capacity (Mass) matrix
            end do
          end do
          Fe(ievab) = Fe(ievab) + basis(inode) * force(idofn) * dvol
        end do
      end do
    
    end subroutine Galerkin_prevTime
   
    
    
    
    !subroutine TimeLevels(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, delta_t, Ucurr, Uprev, F_plus_MU)
    !  
    !  use sourceTerm
    !  
    !  
    !  implicit none
    !  
    !  double precision, allocatable, dimension(:,:), intent(in out) :: Uprev, Ucurr
    !  double precision, dimension(nne,TotGp), intent(in):: N, dN_dxi, dN_deta
    !  double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
    !  !double precision, dimension(3,nne), intent(in)    :: Hesxieta
    !  double precision, dimension(ndofn)        :: EMsource
    !  double precision, dimension(nne)          :: basis, xi_cor, yi_cor
    !  double precision, dimension(DimPr,nne)    :: dN_dxy
    !  double precision, dimension(3,nne)        :: HesXY
    !  double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
    !  double precision, dimension(nevab, nevab) :: Ke, Ce
    !  double precision, dimension(nevab)        :: Fe, Fe_time, ue_prev, ue_curr, AvrgeTime
    !  !double precision, dimension(3,3)          :: tauma
    !  double precision, dimension(nne,DimPr)                :: element_nodes
    !  integer, dimension(nne)                   :: nodeIDmap
    !  double precision                          :: dvol, hmaxi, detJ, delta_t
    !  integer                                   :: igaus, ibase, ielem
    !  double precision, dimension(ntotv, 1), intent(out)  :: F_plus_MU
    !  
    !  !allocate( F_plus_MU )
    !  F_plus_MU = 0.0
    !
    !  do ielem = 1, nelem 
    !    Ke = 0.0; Fe = 0.0; Ce = 0.0
    !    call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
    !    call gather(nodeIDmap, Uprev, ue_prev)
    !    call gather(nodeIDmap, Ucurr, ue_curr)
    !    AvrgeTime =  (2*ue_curr + 0.5*ue_prev) / delta_t 
    !    
    !    !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
    !    do igaus = 1, TotGp
    !      call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
    !      !Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
    !      !detJ = m22det(Jaco)
    !      !Jinv = inv2x2(Jaco)
    !      !dvol = detJ *  weigp(igaus,1)
    !      !call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)
    !      call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
    !      hmaxi = elemSize(Jinv)
    !      do ibase = 1, nne
    !        basis(ibase) = N(ibase,igaus)
    !      end do
    !      
    !      
    !      !call Galerkin(dvol, basis, dN_dxy, Ke, Ce, Fe) 
    !      call source_term(ielem, basis, xi_cor, yi_cor, EMsource)
    !      call Galerkin_prevTime(dvol, basis, Ce, Fe)
    !      !call TauMat(hmaxi,tauma)
    !      !!call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe, pertu,workm,resid)
    !      !if(kstab.ne.6.or.kstab.ne.0)call Stabilization(dvol, basis, dN_dxy, HesXY, EMsource, tauma, Ke, Fe)
    !      !call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe)
    !      Fe_time = Fe - matmul(Ce,AvrgeTime)
    !    end do
    !
    !    call Assemb_Glob_Vec(nodeIDmap, Fe_time, F_plus_MU) !Assemble Global Source vector F
    !    
    !  end do
    !  
    !end subroutine TimeLevels



    
    
    
    
    
    
    
    
  !Fin de contains


end module library
