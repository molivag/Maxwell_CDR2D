module sourceTerm
 use param 

  contains
    
    
    subroutine source_term(basis, xi_cor, yi_cor, EMsource)
      implicit none
      
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
      
      double precision, parameter :: pi = 4*atan(1.d0)
      
      double precision,dimension(nne), intent(in) :: basis, xi_cor, yi_cor
      integer :: ibase
      real    :: n
      double precision :: dey_dydx,dex_dy2,dex_dx2,dey_dxdy,   dey_dx2,dex_dxdy,dex_dydx,dey_dy2
      double precision :: dp_dx, dp_dy, d2p_dx2, d2p_dy2
      double precision :: x, y, alpha, beta, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk
      double precision :: itan, senn, coss
      double precision, dimension(ndofn), intent(out) :: EMsource
      
      
      beta  = Cu*(helem**2/ell**2)
      alpha = ell**2/mu
      n = n_val
      x = 0.0
      y = 0.0
      
      do ibase = 1, nne 
        x = x + basis(ibase)*xi_cor(ibase)
        y = y + basis(ibase)*yi_cor(ibase)
        !print"(3(1x,f10.7))",  basis(ibase), yi_cor(ibase), basis(ibase)*yi_cor(ibase)
      end do
      
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
      
      select case(simul)
      case(1)
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
        EMsource(1) = mu*( dey_dydx - dex_dy2 + beta*( dex_dx2 + dey_dxdy) )
        EMsource(2) = mu*(-dey_dx2 + dex_dxdy + beta*( dex_dydx+ dey_dy2 ) )
        EMsource(3) = force(3)
        
      case(2)
        !Derivatives in x-direction
        dey_dydx = exp(x) * ( sin(y) + y*cos(y) )
        dex_dy2  = exp(x) * (3.0*sin(y)+ y*cos(y) )
        dex_dx2  =-exp(x) * ( sin(y) + y*cos(y) )
        dey_dxdy = exp(x) * ( sin(y) + y*cos(y) )
        
        !Derivatives in y-direction
        dey_dx2  = exp(x) * y*sin(y)
        dex_dxdy = exp(x) * ( y*sin(y) - 2.0*cos(y) )
        dex_dydx = exp(x) * ( y*sin(y) - 2.0*cos(y) )
        dey_dy2  = exp(x) * ( 2.0*cos(y) - y*sin(y) )
        
        !Gradient of multiplier
        dp_dx = 0.5*pi * sin(0.5*pi*(y-1.0)) * cos(0.5*pi*(x-1))
        dp_dy = 0.5*pi * sin(0.5*pi*(x-1.0)) * cos(0.5*pi*(y-1))
        
        !Laplacian of multiplier
        d2p_dx2 = -0.25*pi**2 * sin(0.5*pi*(x-1.0)) * sin(0.5*pi*(y-1))
        d2p_dy2 = -0.25*pi**2 * sin(0.5*pi*(x-1.0)) * sin(0.5*pi*(y-1))
        
        !Source Term computation
        EMsource(1) = mu*( dey_dydx - dex_dy2 + beta*( dex_dx2 + dey_dxdy) ) + dp_dx
        EMsource(2) = mu*(-dey_dx2 + dex_dxdy + beta*( dex_dydx+ dey_dy2 ) ) + dp_dy
        EMsource(3) = alpha * (d2p_dx2 + d2p_dy2)
        
      case(3)
        !Derivatives in x-direction
        dey_dydx = 2.0-(4.0*y)
        dex_dy2  = 0.0 
        
        !Derivatives in y-direction
        dey_dx2  = 0.0
        dex_dxdy = (4.0*x)-2
        
        EMsource(1) = mu*(dey_dydx - dex_dy2)
        EMsource(2) = mu*(-dey_dx2  + dex_dxdy)
        EMsource(3) = 0.0!force(ndofn)
        
        !Source considering the DivDiv term
        
        ! !Derivatives in x-direction
        ! dey_dydx = 2.0-(4.0*y)
        ! dex_dy2  = 0.0 
        ! dex_dx2  = (4.0*y)-2.0
        ! dey_dxdy = 2.0-(4.0*y)
        ! 
        ! 
        ! !Derivatives in y-direction
        ! dey_dx2  = 0.0
        ! dex_dxdy = (4.0*x)-2
        ! dex_dydx = (4.0*x)-2
        ! dey_dy2  = 2.0-4.0*x
        ! 
        ! 
        ! EMsource(1) = mu*( dey_dydx - dex_dy2 + beta*( dex_dx2 + dey_dxdy) )
        ! EMsource(2) = mu*(-dey_dx2 + dex_dxdy + beta*( dex_dydx+ dey_dy2 ) )
        ! EMsource(3) = force(3)
        
      case(4)
        !Source considering the DivDiv term
        
        !Derivatives in x-direction
        dey_dydx = 
        dex_dy2  = 
        dex_dx2  = 
        dey_dxdy = 
        
        
        !Derivatives in y-direction
        dey_dx2  = 
        dex_dxdy = 
        dex_dydx = 
        dey_dy2  = 
        
        
        EMsource(1) = mu*( dey_dydx - dex_dy2 + beta*( dex_dx2 + dey_dxdy) )
        EMsource(2) = mu*(-dey_dx2 + dex_dxdy + beta*( dex_dydx+ dey_dy2 ) )
        EMsource(3) = force(3)
        
      end select
      
    end subroutine source_term

end module sourceTerm

