module biunit 
  use geometry
 
  contains
    
    ! # # # # # # # # # # # # # # # # # # # #  
    
    subroutine GaussQuadrature(ngaus, weigp)
      implicit none
      
      double precision, allocatable, dimension(:,:),intent(out) :: ngaus, weigp
      
      allocate(ngaus(totGp,2))
      allocate(weigp(totGp,1))
      
      ngaus = 0.0
      weigp = 0.0
      
      Select case(ElemType) 
        case('Quadrilateral')
          if(totGp.eq.1)then 
            ngaus(1,1) = 0.0
            ngaus(1,2) = 0.0
            weigp(1,1) = 4.0
          else if(totGp.eq.4) then
            ngaus(1,1) =  sqrt(1.0/3.0)
            ngaus(2,1) = -sqrt(1.0/3.0)
            ngaus(3,1) = -sqrt(1.0/3.0)
            ngaus(4,1) =  sqrt(1.0/3.0)
            
            ngaus(1,2) =  sqrt(1.0/3.0)
            ngaus(2,2) =  sqrt(1.0/3.0)
            ngaus(3,2) = -sqrt(1.0/3.0)
            ngaus(4,2) = -sqrt(1.0/3.0)
           
            weigp(1,1) = 1.0
            weigp(2,1) = 1.0
            weigp(3,1) = 1.0
            weigp(4,1) = 1.0
          else if(totGp.eq.9) then
            ngaus(1,1) =  sqrt(3.0/5.0)
            ngaus(2,1) =  0.0
            ngaus(3,1) = -sqrt(3.0/5.0)
            ngaus(4,1) = -sqrt(3.0/5.0)
            ngaus(5,1) = -sqrt(3.0/5.0) 
            ngaus(6,1) =  0.0
            ngaus(7,1) =  sqrt(3.0/5.0)
            ngaus(8,1) =  sqrt(3.0/5.0)
            ngaus(9,1) =  0.0
            
            ngaus(1,2) =  sqrt(3.0/5.0)
            ngaus(2,2) =  sqrt(3.0/5.0)
            ngaus(3,2) =  sqrt(3.0/5.0)
            ngaus(4,2) =  0.0
            ngaus(5,2) = -sqrt(3.0/5.0)
            ngaus(6,2) = -sqrt(3.0/5.0)
            ngaus(7,2) = -sqrt(3.0/5.0)
            ngaus(8,2) = 0.0
            ngaus(9,2) = 0.0
            
            weigp(1,1) = 25.0/81.0
            weigp(2,1) = 40.0/81.0
            weigp(3,1) = 25.0/81.0
            weigp(4,1) = 40.0/81.0
            weigp(5,1) = 25.0/81.0
            weigp(6,1) = 40.0/81.0
            weigp(7,1) = 25.0/81.0
            weigp(8,1) = 40.0/81.0
            weigp(9,1) = 64.0/81.0
          else
            write(*,*) 'Invalid number of Gauss points for this element (Q1)'   
            stop
          end if
        case("Triangle") 
          if(totGp.eq.1) then
            ngaus(1,1) = 1.0/3.0
            ngaus(1,2) = 1.0/3.0
            weigp(1,1) = 1.0 
          else if(totGp.eq.3) then
            ngaus(1,1) = 1.0/6.0
            ngaus(2,1) = 2.0/3.0
            ngaus(3,1) = 1.0/6.0
            
            ngaus(1,2) = 1.0/6.0
            ngaus(2,2) = 1.0/6.0
            ngaus(3,2) = 2.0/3.0
            
            weigp(1,1) = 1.0/3.0
            weigp(2,1) = 1.0/3.0
            weigp(3,1) = 1.0/3.0
          else if(totGp.eq.7) then
            ngaus(1,1) = 0.1012865073235
            ngaus(2,1) = 0.7974269853531
            ngaus(3,1) = 0.1012865073235
            ngaus(4,1) = 0.4701420641051
            ngaus(5,1) = 0.4701420641051
            ngaus(6,1) = 0.0597158717898
            ngaus(7,1) = 0.3333333333333
            
            ngaus(1,2) = 0.1012865073235
            ngaus(2,2) = 0.1012865073235
            ngaus(3,2) = 0.7974269853531
            ngaus(4,2) = 0.0597158717898
            ngaus(5,2) = 0.4701420641051
            ngaus(6,2) = 0.4701420641051
            ngaus(7,2) = 0.3333333333333
            
            weigp(1,1) = 0.1259391805448 
            weigp(2,1) = 0.1259391805448 
            weigp(3,1) = 0.1259391805448 
            weigp(4,1) = 0.1323941527885 
            weigp(5,1) = 0.1323941527885 
            weigp(6,1) = 0.1323941527885 
            weigp(7,1) = 0.225
          else
            write(*,*) 'Invalid number of Gauss poooints for this element (P1)'   
            stop
          end if
        case DEFAULT
          print*, 'En GaussQuadrature'
          write(*,*) 'Invalid type of element.'
          stop
      end select
      
      
      
    end subroutine GaussQuadrature
    
   
   
    ! subroutine ShapeFunctions(s,t,nnode,shape,derst,hesst)
    subroutine  ShapeFunctions(ngaus, basfun, dN_dxi, dN_deta, hes_ss, hes_st, hes_tt)
      !***********************************************************************
      !
      !     Evaluate the shape functions and its natural first derivatives 
      !     and Hessian matrix
      !
      !***********************************************************************
      implicit None
      
      double precision, dimension(:,:), intent(in) :: ngaus
      double precision, dimension(totGp) :: xi_vector, eta_vector
      double precision                   :: xi, eta,l1,l2,l3
      integer                            :: j
      double precision, allocatable, dimension(:,:), intent(out) :: basfun, dN_dxi, dN_deta
      double precision, allocatable, dimension(:,:), intent(out) :: hes_ss, hes_st, hes_tt
      
      allocate( basfun(Nne,totGp), dN_dxi(Nne,totGp), dN_deta(Nne,totGp) )
      allocate( hes_ss(Nne,totGp), hes_st(Nne,totGp), hes_tt(Nne,totGp) )
      
      
      !         s = xi    t = eta 
      
      
      xi_vector  = ngaus(:,1)     ! xi-coordinate of point j
      eta_vector = ngaus(:,2)
      
      hes_ss = 0.0 
      hes_st = 0.0
      hes_tt = 0.0
      basfun = 0.0 
      dN_dxi = 0.0 
      dN_deta = 0.0
      
      select case(ElemType)
        case('Quadrilateral')
          if(nne.eq.4) then
            
            do j=1,totGp 
              xi = xi_vector(j)    ! xi-coordinate of point j 
              eta= eta_vector(j)   ! eta-coordinate of point j 
              
              basfun(1,j) = 0.25*(1+xi)*(1+eta)
              basfun(2,j) = 0.25*(1-xi)*(1+eta)
              basfun(3,j) = 0.25*(1-xi)*(1-eta)
              basfun(4,j) = 0.25*(1+xi)*(1-eta)
              
              dN_dxi(1,j) = 0.25*(1+eta)
              dN_dxi(2,j) =-0.25*(1+eta)
              dN_dxi(3,j) =-0.25*(1-eta)
              dN_dxi(4,j) = 0.25*(1-eta)
              
              dN_deta(1,j) = 0.25*(1+xi)
              dN_deta(2,j) = 0.25*(1-xi)
              dN_deta(3,j) =-0.25*(1-xi)
              dN_deta(4,j) =-0.25*(1+xi)
              
              hes_st(1,j) = 0.25
              hes_st(2,j) =-0.25
              hes_st(3,j) = 0.25
              hes_st(4,j) =-0.25
            end do
          !else if(nne.eq.8) then
          !  do j = 1, totGP
          !   
          !    xi = xi_vector(j)    ! xi-coordinate of point j 
          !    eta= eta_vector(j)   ! eta-coordinate of point j 
          !    
          !    N(1,j) = 0.25*s9*st*t9
          !    N(2,j) = 0.25*s1*st*t9
          !    N(3,j) = 0.25*s1*st*t1  
          !    N(4,j) = 0.25*s9*st*t1
          !    N(5,j) = 0.5*(1.0-ss)*t*t9
          !    N(6,j) = 0.5*s*s1*(1.0-tt)
          !    N(7,j) = 0.5*(1.0-ss)*t*t1
          !    N(8,j) = 0.5*s*s9*(1.0-tt)
          !    
          !    dN_dxi(1,j) = 
          !    dN_dxi(2,j) = 
          !    dN_dxi(3,j) = 
          !    dN_dxi(4,j) = 
          !    dN_dxi(5,j) = 
          !    dN_dxi(6,j) = 
          !    dN_dxi(7,j) = 
          !    dN_dxi(8,j) = 
          !    
          !    dN_deta(1,j) = 
          !    dN_deta(2,j) = 
          !    dN_deta(3,j) = 
          !    dN_deta(4,j) = 
          !    dN_deta(5,j) = 
          !    dN_deta(6,j) = 
          !    dN_deta(7,j) = 
          !    dN_deta(8,j) = 
          !    
          !    hes_ss(1,j) = 
          !    hes_ss(2,j) = 
          !    hes_ss(3,j) = 
          !    hes_ss(4,j) = 
          !    hes_ss(5,j) = 
          !    hes_ss(6,j) = 
          !    hes_ss(7,j) = 
          !    hes_ss(8,j) = 
          !    
          !    hes_st(1,j) =
          !    hes_st(2,j) =
          !    hes_st(3,j) =
          !    hes_st(4,j) =
          !    hes_st(5,j) =
          !    hes_st(6,j) =
          !    hes_st(7,j) =
          !    hes_st(8,j) =
          !    
          !    hes_tt(1,j) =
          !    hes_tt(2,j) =
          !    hes_tt(3,j) =
          !    hes_tt(4,j) =
          !    hes_tt(5,j) =
          !    hes_tt(6,j) =
          !    hes_tt(7,j) =
          !    hes_tt(8,j) =
          !  end do
          else if(nne.eq.9) then
            do j = 1, totGP
              xi = xi_vector(j)    ! xi-coordinate of point j 
              eta= eta_vector(j)   ! eta-coordinate of point j 
              
              basfun(1,j) = 0.25*xi*(1.+xi)*eta*(1.+eta)
              basfun(2,j) =-0.25*xi*(1.-xi)*eta*(1.+eta)
              basfun(3,j) = 0.25*xi*(1.-xi)*eta*(1.-eta)
              basfun(4,j) =-0.25*xi*(1.+xi)*eta*(1.-eta)
              basfun(5,j) = 0.50*(1.0-xi**2)*eta*(1.0+eta)
              basfun(6,j) =-0.50*xi*(1.-xi)*(1.0-eta**2)
              basfun(7,j) =-0.50*(1.-xi**2)*eta*(1.-eta)
              basfun(8,j) = 0.50*xi*(1.+xi)*(1.0-eta**2)
              basfun(9,j) = (1.0-xi**2)*(1.0-eta**2) 
              
              !basfun(1,j) = 0.25*(xi**2-xi)*(eta**2-eta);
              !basfun(2,j) = 0.25*(xi**2+xi)*(eta**2-eta);
              !basfun(3,j) = 0.25*(xi**2+xi)*(eta**2+eta);
              !basfun(4,j) = 0.25*(xi**2-xi)*(eta**2+eta);
              !basfun(5,j) = 0.5*(1-xi**2)*(eta**2-eta);
              !basfun(6,j) = 0.5*(xi**2+xi)*(1-eta**2);
              !basfun(7,j) = 0.5*(1-xi**2)*(eta**2+eta);
              !basfun(8,j) = 0.5*(xi**2-xi)*(1-eta**2);
              !basfun(9,j) = (1-xi**2)*(1-eta**2);
              
              dN_dxi(1,j) = 0.5*(0.5+xi)*eta*(1.0+eta)
              dN_dxi(2,j) = 0.5*(xi-0.5)*eta*(1.0+eta)
              dN_dxi(3,j) = eta*(xi*(0.5*eta-0.5)-0.25*eta+0.25)
              dN_dxi(4,j) = eta*(xi*(0.5*eta-0.5)+0.25*eta-0.25)
              dN_dxi(5,j) =-xi*eta*(1.+eta) 
              dN_dxi(6,j) = xi*(1.0-eta**2)+0.5*(eta**2)-0.5
              dN_dxi(7,j) =-xi*(eta-1.0)*eta 
              dN_dxi(8,j) = xi*(1.0-eta**2)-0.5*(eta**2)+0.5
              dN_dxi(9,j) = 2*xi*(-1.0+eta**2)
              
              !! Derivatives of S w.r.t xi
              !dN_dxi(1,j) = 0.25*(2*xi-1.0)*(eta**2-eta);
              !dN_dxi(2,j) = 0.25*(2*xi+1.0)*(eta**2-eta);
              !dN_dxi(3,j) = 0.25*(2*xi+1.0)*(eta**2+eta);
              !dN_dxi(4,j) = 0.25*(2*xi-1)*(eta**2+eta);
              !dN_dxi(5,j) = -(xi)*(eta**2-eta);
              !dN_dxi(6,j) = 0.5*(2*xi+1)*(1-eta**2);
              !dN_dxi(7,j) = -(xi)*(eta**2+eta);
              !dN_dxi(8,j) = 0.5*(2*xi-1)*(1-eta**2);
              !dN_dxi(9,j) = (-2*xi)*(1-eta**2);
              
              
              dN_deta(1,j) = 0.25*xi*(1.0+xi)*(2*eta+1) 
              dN_deta(2,j) = xi*(xi*(0.5*eta+0.25)-0.5*eta-0.25)
              dN_deta(3,j) = xi*(xi*(0.5*eta-0.25)-0.5*eta+0.25)
              dN_deta(4,j) = 0.5*xi*(xi+1.0)*(eta-0.5)
              dN_deta(5,j) = xi*xi*(-eta-0.5)+eta+0.5
              dN_deta(6,j) =-(xi-1.0)*xi*eta 
              dN_deta(7,j) = xi**2*(0.5-eta)+eta-0.5
              dN_deta(8,j) =-xi*(xi+1.0)*eta 
              dN_deta(9,j) = 2.0*(-1.0+xi**2)*eta
              
              !! Derivatives of S w.r.t eta
              !dN_deta(1,j) = 0.25*(xi**2-xi)*(2.*eta-1);
              !dN_deta(2,j) = 0.25*(xi**2+xi)*(2.*eta-1);
              !dN_deta(3,j) = 0.25*(xi**2+xi)*(2.*eta+1);
              !dN_deta(4,j) = 0.25*(xi**2-xi)*(2.*eta+1);
              !dN_deta(5,j) = 0.5*(1-xi**2)*(2*eta-1);
              !dN_deta(6,j) = 0.5*(xi**2+xi)*(-2.*eta);
              !dN_deta(7,j) = 0.5*(1-xi**2)*(2.*eta+1);
              !dN_deta(8,j) = 0.5*(xi**2-xi)*(-2.*eta);
              !dN_deta(9,j) = (1-xi**2)*(-2.*eta);
              

              hes_ss(1,j) = 0.5*eta*(1.0+eta) 
              hes_ss(2,j) = 0.5*eta*(eta+1.0) 
              hes_ss(3,j) = 0.5*eta*(eta-1.0) 
              hes_ss(4,j) = 0.5*eta*(eta-1.0) 
              hes_ss(5,j) =-eta*(1.0+eta) 
              hes_ss(6,j) = (1.0-eta**2) 
              hes_ss(7,j) =-eta*(eta-1.0)
              hes_ss(8,j) = (1.0-eta**2)
              hes_ss(9,j) = 2.0*(-1.0+eta**2)
              
              hes_st(1,j) = xi*(eta+0.5)+0.5*eta+0.25 
              hes_st(2,j) = xi*(eta+0.5)-0.5*eta-0.25 
              hes_st(3,j) = xi*(eta-0.5)-0.5*eta+0.25  
              hes_st(4,j) = xi*(eta-0.5)+0.5*eta-0.25  
              hes_st(5,j) =-xi*(2.0*eta+1.0)
              hes_st(6,j) = eta*(1.0-2.0*xi)  
              hes_st(7,j) = xi*(1.0-2.0*eta) 
              hes_st(8,j) =-eta*(2.0*xi+1.0)
              hes_st(9,j) = 4.0*xi*eta
              
              hes_tt(1,j) = 0.5*xi*(1.0+xi) 
              hes_tt(2,j) = 0.5*xi*(xi-1.0) 
              hes_tt(3,j) = 0.5*xi*(xi-1.0) 
              hes_tt(4,j) = 0.5*xi*(1.0+xi) 
              hes_tt(5,j) = 1.0-xi**2 
              hes_tt(6,j) =-xi*(xi-1.0)
              hes_tt(7,j) = 1.0-xi**2
              hes_tt(8,j) =-xi*(1.0+xi)
              hes_tt(9,j) = 2.0*(-1.0+xi**2) 
            end do
          else
          write(*,*) 'Invalid number of nodes in element. Q must be 4 or 9'
          stop
          end if
        case('Triangle')
          if(nne.eq.3) then
            do j = 1, totGP
             
              xi = xi_vector(j)    ! xi-coordinate of point j 
              eta= eta_vector(j)   ! eta-coordinate of point j 
              
              l1 = 1.0-xi-eta
              l2 = xi
              l3 = eta
              
              basfun(1,j) = l1
              basfun(2,j) = l2
              basfun(3,j) = l3
              
              dN_dxi(1,j) =-1.0
              dN_dxi(2,j) = 1.0
              dN_dxi(3,j) = 0.0
              
              dN_deta(1,j) =-1.0
              dN_deta(2,j) = 0.0
              dN_deta(3,j) = 1.0
            end do
          else if(nne.eq.6) then
            do j = 1, totGP
              xi = xi_vector(j)    ! xi-coordinate of point j 
              eta= eta_vector(j)   ! eta-coordinate of point j 
              
              l1 = 1.0-xi-eta
              l2 = xi
              l3 = eta
              
              basfun(1,j) = 2.0*l1*(l1-0.5)
              basfun(2,j) = 2.0*l2*(l2-0.5)
              basfun(3,j) = 2.0*l3*(l3-0.5)
              basfun(4,j) = 4.0*l1*l2
              basfun(5,j) = 4.0*l2*l3
              basfun(6,j) = 4.0*l3*l1 
              
              dN_dxi(1,j) = 4.0*(xi+eta-0.75)
              dN_dxi(2,j) = 4.0*(xi-0.25)
              dN_dxi(3,j) = 0.0000000000000
              dN_dxi(4,j) =-4.0*((2.0*xi)+eta-1.0)
              dN_dxi(5,j) = 4.0*eta
              dN_dxi(6,j) =-4.0*eta
              
              dN_deta(1,j) = 4.0*(xi+eta-0.75) 
              dN_deta(2,j) = 0.0000000000000
              dN_deta(3,j) = 4*(eta-0.25)
              dN_deta(4,j) =-4.0*xi 
              dN_deta(5,j) = 4.0*xi
              dN_deta(6,j) =-4.0*(xi+(2.0*eta)-1.0)
              
              hes_ss(1,j) = 4.0 
              hes_ss(2,j) = 4.0
              hes_ss(4,j) =-8.0
              
              hes_st(1,j) = 4.0
              hes_st(4,j) =-4.0
              hes_st(5,j) = 4.0
              hes_st(6,j) =-4.0
              
              hes_tt(1,j) = 4.0 
              hes_tt(3,j) = 4.0
              hes_tt(6,j) =-8.0
              
            end do
          else
          write(*,*) 'Invalid number of nodes in element. P must be 3 or 6'
          stop
          end if
        case DEFAULT
          print*, 'In ShapeFucntions'
          write(*,*) 'Invalid type of element. Must be Q or P'
          stop
      end select
     
    end subroutine ShapeFunctions
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !subroutine ShapeFunctions(ngaus, Nne,  N, dN_dxi, dN_deta, Hesxieta) 
    !  
    !  implicit None
    !  
    !  double precision, dimension(:,:), intent(in) :: ngaus
    !  double precision, dimension(totGp) :: xi_vector, eta_vector
    !  double precision, dimension(3,Nne) :: Hesxieta
    !  integer, dimension(Nne,DimPr)      :: master_nodes
    !  double precision                   :: xi, eta, mn_xi, mn_eta
    !  integer                            :: i, j, jj, k, Nne
    !  double precision, allocatable, dimension(:,:), intent(out) :: N
    !  double precision, allocatable, dimension(:,:), intent(out), optional :: dN_dxi, dN_deta
    !  
    !  ! = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !  
    !  
    !  allocate( N(Nne,totGp) )
    !  N = 0.0
    !  xi_vector  = ngaus(:,1)     ! xi-coordinate of point j
    !  eta_vector = ngaus(:,2)

    !  select case(ElemType)

    !    CASE ('Quadrilateral')
    !     select case(Nne)
    !      case  (8)
    !        ! Shape functions for square (quadrilaters) linear elements
    !        !
    !        !  |
    !        !  |
    !        !  | o- - o - -o
    !        !  Y |         |
    !        !  | o         o
    !        !  | |         |
    !        !  | o- - o - -o
    !        !  |
    !        !  +--------X-------->
    !        !write(*,100) ' Element type: ', ElemType, 'whit', Nne, 'nodes per element'
    !        !coordinates of the nodes of the master element
    !        master_nodes = reshape([1, -1, -1, 1, 0, -1, 0, 1, 1, 1, -1, -1, 1, 0, -1, 0], [Nne,DimPr])
    !        !NOTA ** Para que el reshape funcione correctamente, o que produzca el par de valores deseado, 
    !        !primero se deben colocar todos los valores en x, luego todos los de y y luego, si hubiera todos 
    !        !los de z para que al acomodarse salga el par suponiendo un reshape de 3,2 debe acomodarse:
    !        !x1, x2, x3, y1, y2, y3 DUDA *Siempre es asi*
    !        
    !        
    !        ! dN(xi,eta)/dx = dN/dxi(dxi/dx) + dN/deta(deta/dx)
    !        ! dN(xi,eta)/dy = dN/dxi(dxi/dy) + dN/deta(deta/dy)
    !        ! Aqui se calculan as funciones de forma N y parte de las derivadas dN/dxi and dN_deta
    !        ! mas no las derivadas dN/dx and dN/dy completas
    !        if (present(dN_dxi) .and. present(dN_deta))then
    !          !write(*,*) 'Using derivatives of shape functions -', Nne
    !          allocate(dN_dxi(Nne,totGp) )
    !          allocate(dN_deta(Nne,totGp))
    !          dN_dxi  = 0.0
    !          dN_deta = 0.0
    !          do j = 1, totGp
    !            
    !            dN_dxi(5,j)=-xi*(1+eta)
    !            dN_dxi(6,j)=-1.0/2*(1-eta**2)
    !            dN_dxi(7,j)=-xi*(1-eta)
    !            dN_dxi(8,j)=1.0/2*(1-eta**2)
    !            dN_deta(5,j)=1.0/2*(1-xi**2)
    !            dN_deta(6,j)=(1-xi)*(-eta)
    !            dN_deta(7,j)=-1.0/2*(1-xi**2)
    !            dN_deta(8,j)=(1+xi)*(-eta)
    !            
    !            do i = 1, 4
    !              mn_xi = master_nodes(i,1)
    !              mn_eta= master_nodes(i,2)
    !              if (i==1) then
    !                jj=8
    !              else
    !                jj=i+3
    !              end if
    !              k=i+4
    !              dN_dxi(i,j)= mn_xi*(1.0 + mn_eta*eta)/4.0 - 1.0/2*(dN_dxi(jj,j)+dN_dxi(k,j))
    !              dN_deta(i,j)= mn_eta*(1.0 + mn_xi*xi)/4.0 - 1.0/2*(dN_deta(jj,j)+dN_deta(k,j))
    !            end do
    !            
    !          end do
    !        else
    !          !write(*,*) 'Using only shape functions'
    !          continue
    !        end if
    !        !Despues de evaluar si estan las derivadas como variable dummy 
    !        !en la llamada construye las
    !        !funciones de forma.
    !        do j = 1, totGp
    !          xi  = xi_vector(j) ! xi-cor of point j
    !          eta = eta_vector(j)! eta-cor of point j
    !          N(5,j)=1.0/2*(1-xi**2)*(1+eta)
    !          N(6,j)=1.0/2*(1-xi)*(1-eta**2)
    !          N(7,j)=1.0/2*(1-xi**2)*(1-eta)
    !          N(8,j)=1.0/2*(1+xi)*(1-eta**2)
    !          do i = 1, 4
    !            mn_xi = master_nodes(i,1)
    !            mn_eta= master_nodes(i,2)
    !            if (i==1) then
    !              jj=8
    !            else
    !              jj=i+3
    !            end if
    !            k=i+4
    !            N(i,j)=(1.0 + mn_xi*xi)*(1.0 + mn_eta*eta)/4.0 - 1.0/2*(N(jj,j)+N(k,j))
    !          end do
    !        end do
    !        
    !      case  (4)
    !        ! Shape functions for square (quadrilaters) linear elements
    !        !
    !        !  |
    !        !  |
    !        !  | o- - - -o
    !        !  Y |       |
    !        !  | |       |
    !        !  | o- - - -o
    !        !  |
    !        !  +--------X-------->
    !        
    !        !write(*,100) ' Element type: ', ElemType, 'whit', Nne, 'nodes per element'
    !        !coordinates of the nodes of the master element
    !        master_nodes = reshape([1, -1, -1, 1, 1, 1, -1, -1], [Nne,DimPr])
    !        ! dN(xi,eta)/dx = dN/dxi(dxi\dx) + dN/deta(deta/dx)
    !        ! dN(xi,eta)/dy = dN/dxi(dxi\dy) + dN/deta(deta/dy)
    !        ! Aqui se calculan as funciones de forma N y parte de las derivadas dN/dxi and dN_deta
    !        ! mas no las derivadas dN/dx and dN/dy completas
    !        !do loop: compute N, dN_dxi, dN_deta
    !        if (present(dN_dxi) .and. present(dN_deta))then
    !          !write(*,*) 'Using derivatives of shape functions', Nne
    !          allocate(dN_dxi(Nne,totGp) )
    !          allocate(dN_deta(Nne,totGp) )
    !          dN_dxi  = 0.0
    !          dN_deta = 0.0
    !          do j=1,totGp                              ! columns for point 1,2 ...
    !            xi=xi_vector(j);                      ! xi-coordinate of point j 
    !            eta=eta_vector(j);                    ! eta-coordinate of point j 
    !            do i=1,4                              ! rows for N1, N2, ...
    !              mn_xi = master_nodes(i,1)
    !              mn_eta= master_nodes(i,2)
    !              dN_dxi(i,j)= mn_xi*(1.0 + mn_eta*eta)/4.0             ! dNi/dxi(xi,eta)
    !              dN_deta(i,j)= mn_eta*(1.0 + mn_xi*xi )/4.0            ! dNi/deta(xi,eta
    !            end do
    !          end do
    !          Hesxieta(2,1)= 0.25
    !          Hesxieta(2,2)=-0.25
    !          Hesxieta(2,3)= 0.25
    !          Hesxieta(2,4)=-0.25
    !        else
    !          !write(*,*) 'Using only shape functions'
    !          continue
    !        endif
    !        
    !        do j=1,totGp                            ! columns for point 1,2 ...
    !          xi=xi_vector(j);                      ! xi-coordinate of point j 
    !          eta=eta_vector(j);                    ! eta-coordinate of point j 
    !          do i=1,4                              ! rows for N1, N2, ...
    !            mn_xi = master_nodes(i,1)
    !            mn_eta= master_nodes(i,2)
    !            N(i,j)=(1.0 + mn_xi*xi)*(1.0 + mn_eta*eta)/4.0        ! Ni(xi,eta)
    !          end do
    !        end do
    !        
    !      case DEFAULT
    !        write(*,*) 'Invalid number of nodes in the element.'
    !        stop
    !      end select
    !      
    !    CASE ("Triangle")      
    !      if(Nne .EQ. 3) then
    !        !  |
    !        !  |        o
    !        !  |       / \
    !        !  |      /   \
    !        !  Y     /     \
    !        !  |    /       \
    !        !  |   o---------o
    !        !  |
    !        !  +--------X-------->
    !        if (present(dN_dxi) .and. present(dN_deta))then
    !          !write(*,*) 'Using derivatives of shape functions', Nne
    !          allocate(dN_dxi(Nne,totGp) )
    !          allocate(dN_deta(Nne,totGp) )
    !          dN_dxi  = 0.0
    !          dN_deta = 0.0
    !          do j=1,totGp
    !            xi=xi_vector(j);  ! xi-coordinate of point j 
    !            eta=eta_vector(j); 
    !            N(1,j)      = 1.0-xi-eta
    !            N(2,j)      = xi
    !            N(3,j)      = eta
    !            dN_dxi(1,j) = -1.0
    !            dN_dxi(2,j) = 1.0
    !            dN_dxi(3,j) = 0.0
    !            dN_deta(1,j)= -1.0
    !            dN_deta(2,j)= 0.0
    !            dN_deta(3,j)= 1.0
    !          end do
    !        else
    !          !write(*,*) 'Using only shape functions'
    !          continue
    !        endif
    !        do j=1,totGp
    !          xi=xi_vector(j);    ! xi-coordinate of point j 
    !          eta=eta_vector(j); 
    !          N(1,j) = 1.0-xi-eta
    !          N(2,j) = xi
    !          N(3,j) = eta
    !        end do
    !      elseif(Nne .EQ. 6) then
    !        !  |
    !        !  |        o
    !        !  |       / \
    !        !  |      /   \
    !        !  Y     o     o
    !        !  |    /       \
    !        !  |   /         \
    !        !  |  o-----o-----o
    !        !  |
    !        !  +--------X-------->
    !        if (present(dN_dxi) .and. present(dN_deta))then
    !          allocate(dN_dxi(Nne,totGp) )
    !          allocate(dN_deta(Nne,totGp) )
    !          dN_dxi  = 0.0
    !          dN_deta = 0.0
    !          do j=1,totGp
    !            xi=xi_vector(j);  ! xi-coordinate of point j 
    !            eta=eta_vector(j); 
    !            N(1,j)=(2.0*(1-xi-eta)-1.0)*(1-xi-eta)
    !            N(2,j)=(2.0*xi-1.0)*xi
    !            N(3,j)=(2.0*eta-1.0)*eta
    !            N(4,j)= 4.0*(1-xi-eta)*xi
    !            N(5,j)= 4.0*xi*eta
    !            N(6,j)= 4.0*(1-xi-eta)*eta
    !            dN_dxi(1,j)=1.0-4.0*(1-xi-eta)
    !            dN_dxi(2,j)=4.0*xi-1.0
    !            dN_dxi(3,j)=0.0
    !            dN_dxi(4,j)=4.0*((1-xi-eta)-xi)
    !            dN_dxi(5,j)=4.0*eta
    !            dN_dxi(6,j)=-4.0*eta
    !            dN_deta(1,j)=1.0-4.0*(1-xi-eta)
    !            dN_deta(2,j)=0.0
    !            dN_deta(3,j)=4.0*eta-1.0
    !            dN_deta(4,j)=-4.0*xi
    !            dN_deta(5,j)=4.0*xi
    !            dN_deta(6,j)=4.0*((1-xi-eta)-eta)
    !            !dN_dxi(1,j) = -2*(1-xi-eta-0.5)-2*(1-xi-eta) 
    !            !dN_dxi(2,j) = 2*(xi-0.5) + 2*xi 
    !            !dN_dxi(3,j) = 0
    !            !dN_deta(1,j)= -2*(1-xi-eta-0.5)-2*(1-xi-eta) 
    !            !dN_deta(2,j)= 0
    !            !dN_deta(3,j)= 2*(eta-0.5) + 2*eta 
    !          end do
    !        else
    !          !write(*,*) 'Using only shape functions'
    !          continue
    !        endif
    !        N(1,j)=(2.0*(1-xi-eta)-1.0)*(1-xi-eta)
    !        N(2,j)=(2.0*xi-1.0)*xi
    !        N(3,j)=(2.0*eta-1.0)*eta
    !        N(4,j)= 4.0*(1-xi-eta)*xi
    !        N(5,j)= 4.0*xi*eta
    !        N(6,j)= 4.0*(1-xi-eta)*eta
    !      else
    !        print*, 'Invalid number of node in element' 
    !        stop
    !      end if
    !     
    !    case DEFAULT
    !      print*, 'En ShapeFucntions'
    !      write(*,*) 'Invalid type of element.'   
    !      stop
    !  end select
    ! 
    !end subroutine ShapeFunctions
    
  !end contains

end module biunit
