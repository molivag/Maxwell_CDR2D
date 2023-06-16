module sourceTerm
 use param 
 use geometry

  contains
    
    
    subroutine source_term(ielem, basis, xi_cor, yi_cor, EMsource)
      implicit none
      
      
      double precision, parameter :: pi = 4*atan(1.d0)
      
      double precision,dimension(nne), intent(in) :: basis, xi_cor, yi_cor
      integer                        , intent(in) :: ielem
      double precision :: dey_dydx,dex_dy2,dex_dx2,dey_dxdy,   dey_dx2,dex_dxdy,dex_dydx,dey_dy2
      double precision :: x, y, alpha, beta, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk
      double precision :: itan, senn, coss
      integer          :: ibase
      real             :: n
      double precision, dimension(ndofn), intent(out) :: EMsource
      !double precision :: xq, yq, Icurr
      
      
      EMsource = 0.0
      
      x = 0.0
      y = 0.0
      
      do ibase = 1, nne 
        x = x + basis(ibase)*xi_cor(ibase)
        y = y + basis(ibase)*yi_cor(ibase)
        !print"(3(1x,f10.7))",  basis(ibase), yi_cor(ibase), basis(ibase)*yi_cor(ibase)
      end do
      
      select case(simul)
       
        case(0)
          EMsource(1) = force(1)
        case(1)
          
          !***********************************************************!
          !The source term is given by:                               !
          !                                                           !
          !              u = grad(r^{2n/3}*sin((2n/3)*theta))         !
          ! where:                                                    !
          ! r = sqrt(x^2 + y^2)   ;   theta = atan(y/x)    
          !                                                           !
          ! and                                                       !
          !         f = Lu  ;   where L is the diferential operator   !
          !                                                           !
          !***********************************************************!
          
          beta  = Cu*lambda*(helem**2/ell**2)
          alpha = ell**2/lambda
          n = n_val
          
          !terms for derivatives
          aa   = (2.0/27.0)*n
          bb   = (x**2 + y**2) 
          cc   = (n/3.0 - 5.0/2.0)
          dd   = (8.0*n**2 - 24.0*n + 27.0) 
          ee   = (2.0/3.0)*n
          ff   = (y/x) 
          gg   = 4.0*(n-3.0)*n
          hh   = (x**2 - y**2)  
          ii   = (4.0*n**2 - 6.0*n + 9.0) 
          jj   = (2.0*n**2 - 9.0*n + 9.0) 
          kk   = (8.0*x*y)*(n-3.0)*n 
          itan = atan(ff) 
          senn = sin(ee*itan)
          coss = cos(ee*itan)
          
          !Derivatives in x-direction
          dey_dydx = aa * bb**cc *( x*y*dd *coss - gg*hh *senn )  
          dex_dy2  =-aa * bb**cc *( ((x**2)*ii - 2.0*(y**2)*jj)*senn - kk*coss )  
          dex_dx2  = aa * bb**cc *( (2.0*(x**2)*jj - (y**2)*ii)*senn - kk*coss )
          dey_dxdy = aa * bb**cc *( x*y*dd *coss - gg*hh *senn )  
          
          !Derivatives in y-direction
          dey_dx2  = aa * bb**cc *( (2.0*(x**2)*jj - (y**2)*ii)*coss + kk*senn )
          dex_dxdy = aa * bb**cc *( x*y*dd *senn + gg*hh *coss )
          dex_dydx = aa * bb**cc *( x*y*dd *senn + gg*hh *coss )  
          dey_dy2  =-aa * bb**cc *( ((x**2)*ii - 2.0*(y**2)*jj)*coss + kk*senn ) 
          
          !Source Term computation
          EMsource(1) = lambda*( dey_dydx - dex_dy2 + beta*( dex_dx2 + dey_dxdy) )
          EMsource(2) = lambda*(-dey_dx2 + dex_dxdy + beta*( dex_dydx+ dey_dy2 ) )
          EMsource(3) = force(3)
          
        case(2)
          !The transient source in a full-space geophysicall modelling is a Dirac Delta
          !EMsource(1) = force(1)*A_F((srcLoc-1)*ndofn,1)*u(t)
          !EMsource(2) = force(2)*A_F((srcLoc-1)*ndofn,1)*u(t)
          !EMsource(3) = force(3)*A_F((srcLoc-1)*ndofn,1)*u(t)
          
          !print*, 'The static case'
          EMsource(1) = force(1)
          EMsource(2) = force(2)
          EMsource(3) = force(3)
         
        case(3)
          !Source with NO DivDiv term
          
          !Derivatives in x-direction
          dey_dydx =-2.0*(6.0*x**2 - 6.0*x + 1.0)*(y-1.0)*y*(2.0*y - 1.0)
          dex_dy2  = 6*(x**2 - 2.0*x + 1.0)*(2.0*y-1.0)*x**2
          dex_dx2  = 2.0*(6.0*x**2 - 6.0*x +1.0)*y*(2.0*y**2 - 3.0*y + 1.0)
          dey_dxdy =-2.0*(6.0*x**2 - 6.0*x + 1.0)*(y-1.0)*y*(2.0*y - 1.0) 
          
          
          !Derivatives in y-direction
          dey_dx2  =-6.0*(2.0*x - 1.0)*(y**2 -2.0*y + 1.0)*y**2
          dex_dxdy = 2.0*(x-1.0)*x*(2.0*x - 1.0)*(6.0*y**2 -6.0*y + 1.0)
          dex_dydx = 2.0*(x-1.0)*x*(2.0*x - 1.0)*(6.0*y**2 -6.0*y + 1.0) 
          dey_dy2  =-2.0*x*(2*x**2 - 3.0*x +1.0)*(6.0*y**2 -6.0*y + 1.0)
          
          
          EMsource(1) = force(1)*lambda*( dey_dydx - dex_dy2)! - beta*(dex_dx2  - beta*dey_dxdy ) 
          EMsource(2) = force(2)*lambda*(-dey_dx2 + dex_dxdy)! - beta*(dex_dydx - beta*dey_dy2  )  
          EMsource(3) = force(3)
          
        case(4) !stokes
          
          
          !Derivatives in x-direction
          dex_dx2  = 2.0*(6.0*x**2 - 6.0*x +1.0)*y*(2.0*y**2 - 3.0*y + 1.0)
          dex_dy2  = 6*(x**2 - 2.0*x + 1.0)*(2.0*y-1.0)*x**2
          
          
          !Derivatives in y-direction
          dey_dx2  =-6.0*(2.0*x - 1.0)*(y**2 -2.0*y + 1.0)*y**2
          dey_dy2  =-2.0*x*(2*x**2 - 3.0*x +1.0)*(6.0*y**2 -6.0*y + 1.0)
          
          
          EMsource(1) = -lambda * ( dex_dx2 + dex_dy2 )
          EMsource(2) = -lambda * ( dey_dx2 + dey_dy2 )
          EMsource(3) = force(3)
          
        case(5)
          
          
          EMsource(1) = -lambda * force(1)
          EMsource(2) = -lambda * force(2)
          EMsource(3) =  force(3)
          
        case default
        
        ! - - - Coordinate of source location node 2501 at the center of the mesh 158x158
        !xq = 1.0
        !yq = 1.0     
        !Icurr = force(1)
        !
        !!x = 0.0
        !!Dos posibilidades de implementacion
        !! - - Primera opcion                                                             if(ProbType.ne.'TRAN')then
        !do inode = 1,nne
        !  if((xi_cor(inode).eq.xq).and.(yi_cor(inode).eq.yq))then
        !    EMsource = Icurr !Here is the source location
        !  else
        !    EMsource = 0.0
        !  endif
        !end do
       
        ! - - Segunda opcion
        !EMsource = Icurr *dirac(xi_cor - xq)*dirac(yi_cor - yq)
        !print*, EMsource
       
      end select                                                                            
      
    end subroutine source_term
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine currDensity(time,eTime,Jsource)
      
      implicit none
      
      !double precision,allocatable, dimension(:,:), intent(in) :: A_F
      double precision, intent(in)                        :: eTime
      integer :: tw, inode, time
      double precision, intent(out) :: Jsource(ntotv,1)
      
      !allocate(A_F(ntotv,1))
      
      Jsource = 0.0
      !locatx= (srcLoc(inode)-1)*ndofn+1
      !locaty= (srcLoc(inode)-1)*ndofn+2
      
      !print"(A,I0,A3,f5.3)", 'u(',time,')= ', etime 
      
      !selectcase(srcType)
       ! case(1) !source parallel to x-axis
          do inode=1,nodalSrc
            Jsource((srcLoc(inode)-1)*ndofn+1,1) = -Icurr(1)*eTime
            Jsource((srcLoc(inode)-1)*ndofn+2,1) = -Icurr(2)*eTime
          end do
          
        !case(2) !source parallel to y-axis
        !case(3) !source in a diagonal direction
        !case default
          !write(*,*) 'Source type not available'
      !end select
      
    end subroutine   
    !          
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =                              
    !          
    double precision function dirac(x)
      double precision, intent(in) :: x
      double precision :: absx
      absx = abs(x)
      dirac = merge(1.0,0.0,absx<epsilon(x))
    end function dirac
    !          
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =                              
    !          
    !FUNCTION DIRAC(X,XMIN,XMAX)
    !  ! Assumes XMIN<=XMAX...
    !  DIRAC=0.0D+00
    !  if(X.LT.XMIN.OR.X.GT.XMAX) DIRAC=1.0D+00
    !  RETURN
    !END
    !          
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =                              
    !          
    !program test_dirac
    !    implicit none
    !    
    !    integer :: i
    !    real(real64) :: x
    !    
    !    do i=-10,10
    !       x = real(i,real64)
    !       print *, 'x=',x,' dirac=',dirac(x)
    !    end do
    !    
    !    contains
    !    
    !    elemental real(real64) function dirac(x)
    !        real(real64), intent(in) :: x
    !        real(real64) :: absx
    !        absx = abs(x)
    !        dirac = merge(1.0_real64/absx,0.0_real64,absx<epsilon(x))
    !    end function dirac
    !end program test_dirac
    
    
    
end module sourceTerm

