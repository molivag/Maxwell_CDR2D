module BoundVal
  use param
  use geometry

  contains
    
    subroutine SetBoundVal( nBVs, nBVscol )
      !========================================================================
      !Esta subroutina revisa todos los nodos de la malla y define el tipo de Dof de
      !los nodos en la frontera. Abre un archivo en donde comenzara a escribir, 
      !en la primer columna: el numero de nodo. 
      ! La segunda columna tendra el tipo de nodo
      ! 1 = Nodo con valor preescrito (nodo con valor definido)
      ! 0 = Nodo libre (su valor debera ser calculado)
      ! La tercera columna tendra el valor del nodo preescrito  
      !=========================================================================
      
      implicit none
      
      double precision, parameter :: pi = 4*atan(1.d0)
      character(len=*), parameter :: fileplace ="./"
      
      integer :: ierror, a ,b, c, d, e, f,i 
      double precision :: x, y, xmin, xmax, ymin, ymax, xmiddle, ymiddle, ux, uy, mult
      double precision :: aa, cc
      integer, intent(out) :: nBVs, nBVscol
      
      
      open(unit=200, file=fileplace//'ifpre.dat',Status= 'replace', action= 'write',iostat=ierror)
      open(unit=300, file=fileplace//'BoVal.dat',Status= 'replace', action= 'write',iostat=ierror)
      
      a = 0
      b = 0
      c = 0
      d = 0
      e = 0
      f = 0
      
      xmin = minval(coord(1,:)) !the smallest number in x column
      xmax = maxval(coord(1,:)) !the greatest number in x column
      ymin = minval(coord(2,:)) !the smallest number in y column
      ymax = maxval(coord(2,:)) !the greatest number in y column
      xmiddle =  0.0 !xmax/2.0
      ymiddle =  0.0 !ymax/2.0
      !print*, 'xmin = ', xmin
      !print*, 'xmax = ', xmax
      !print*, 'ymin = ', ymin
      !print*, 'ymax = ', ymax
      !print*, 'xmiddle = ', xmiddle
      !print*, 'ymiddle = ', ymiddle
      
      
      
      if(ndofn.eq.1) then
        do i =1, nnodes
          x=coord(1,i)
          y=coord(2,i)
          if(y.eq.ymax) then
            if(x.eq.xmax)then                   !right top corner 
              write(200,10) i, 1
              write(300,30) 1.0
            elseif(x.eq.xmin)then               !left top corner 
              write(200,10) i, 1
              write(300,30) 1.0
            else
              write(200,10) i, 1                !Top boundary
              write(300,30) 1.0
            end if
            a = a+1
          else if (y.eq.ymin)then
            if(x.eq.xmin)then                   !right bottom corner  
              write(200,10) i, 1
              write(300,30) 0.0
            elseif(x.eq.xmax)then               !left bottom corner
              write(200,10) i, 1
              write(300,30) 0.0
            else
              write(200,10) i, 1                !bottom boundary
              write(300,30) 0.0 
            end if
            b = b+1
          else if(x.eq.xmax)then                !right boundary
              write(200,10) i, 1
              write(300,30) 0.0
            c = c+1
          else if (x.eq.xmin)then               !left boundary
              write(200,10) i, 1
              write(300,30) 0.0
            d = d+1
           
          end if
          nBVs = a+b+c+d
         
        end do
        nBVscol = 3
        
        
      elseif(ndofn .eq. 2) then
        do i =1, nnodes
          x=coord(1,i)
          y=coord(2,i)
          if(y.eq.ymax) then 
            if(x.eq.xmax)then                         
              write(200,10) i, 1, 1
              write(300,40)   1.0, 0.0
            elseif(x.eq.xmin)then                     
              write(200,10) i, 1, 1
              write(300,40)   1.0, 0.0
            else
              write(200,10) i, 1, 1
              write(300,40)   1.0, 0.0
            end if
            a = a+2
          else if (y.eq.ymin)then
            if(x.eq.xmin)then                         
              write(200,10) i, 1, 1
              write(300,40)   0.0, 0.0
            elseif(x.eq.xmax)then                      
              write(200,10) i, 1, 1
              write(300,40)   0.0, 0.0
            else
              write(200,10) i, 1, 1
              write(300,40)   0.0, 0.0
            end if
            b = b+2
          else if(x.eq.xmax)then                      
              write(200,10) i, 1, 1
              write(300,40)   0.0, 0.0
            c = c+2
          else if (x.eq.xmin)then                     
              write(200,10) i, 1, 1
              write(300,40)   0.0, 0.0
            d = d+2
           
          end if
          
          nBVs = a+b+c+d+e
         
        end do
        nBVs = nBVs/2
        nBVscol = 5
        
      elseif(ndofn .eq. 3)then
        select case(simul)
          case(1)
            aa = (2.0/3.0)*n_val
            cc = (n_val/3.0) - (1.0/2.0)
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              uy = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then
                  ux = aa*(2.0)**cc *sin(n_val*(pi/6.0)) 
                  uy = aa*(2.0)**cc *cos(n_val*(pi/6.0))
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0                        !Right top Corner
                  
                elseif(x.eq.xmin)then
                  ux = aa*(2.0)**cc *sin(-n_val*(pi/6.0)) 
                  uy = aa*(2.0)**cc *cos(-n_val*(pi/6.0))
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0                        !Left top Corner
                  
                else
                  ux = aa*(x**2 + 1.0)**cc *sin(aa*datan(1.0/x))
                  write(200,10) i, 1, 0, 1
                  write(300,20)   ux, uy, 0.0                        !Top boundary
                end if
                
                !if(x.eq.xmiddle)then
                !  write(200,10) i, 1, 1, 1
                !  write(300,20)   1.0, 0.0, 0.0                     !middle top preasure boundary
                !end if
                a = a+3
                
              else if(y.eq.ymin)then
                if(x.eq.xmin)then
                  ux = aa*(2.0)**cc *sin(n_val*(pi/6.0))
                  uy = aa*(2.0)**cc *cos(n_val*(pi/6.0))
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0                        !Left bottom Corner
                 
                elseif(x.eq.xmiddle)then
                  ux = aa*(1.0)**cc * sin(n_val*(pi/3.0))
                  uy = aa*(1.0)**cc * cos(n_val*(pi/3.0))
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0                        !central bottom Corner (0,-1)
                  
                else
                  ux = aa*(x**2+1.0)**cc *sin(aa*datan(-1.0/x))
                  write(200,10) i, 1, 0, 1                           !Bottom boundary
                  write(300,20)   ux, uy, 0.0
                  
                end if
                b = b+3
                
              else if(x.eq.xmax)then
                if(y.gt.ymiddle)then
                  uy = aa*(1.0+y**2)**cc *cos(aa*datan(y))  
                  write(200,10) i, 0, 1, 1                           !Right boundary
                  write(300,20)   ux, uy, 0.0
                  
                elseif(y.eq.ymiddle)then
                  ux = 0.0 
                  uy = aa*(1.0 )**cc
                  write(200,10) i, 1, 1, 1                           !middle right corner
                  write(300,20)   ux, uy, 0.0
                  
                endif
                c = c+3
              else if (x.eq.xmin)then
                  uy = aa *(1.0+y**2)**cc *cos(aa*datan(-y))  
                  write(200,10) i, 0, 1, 1                           !Left boundary
                  write(300,20)   ux, uy, 0.0
                d = d+3
               
              else if(x.eq.xmiddle)then
                if(y .eq. ymiddle)then
                  ux = 0.0
                  uy = 0.0
                  write(200,10) i, 1, 1, 1                           !central corner  (0,0)
                  write(300,20)   ux, uy, 0.0
                  
                elseif((y .gt. ymin).and.(y.lt.0.0))then
                  uy = aa*(y**2)**cc *cos(n_val*(pi/3.0))
                  write(200,10) i, 0, 1, 1                     !central vertical boundary (x=0)
                  write(300,20)   ux, uy, 0.0
                end if
                e = e+3
               
              else if(y.eq.ymiddle)then
                if(x .gt. xmiddle)then
                  ux = 0.0
                  write(200,10) i, 1, 0, 1
                  write(300,20)   ux, uy, 0.0                 !horizontal boundary at the middle
                end if
               
                f = f+3
              end if
              
              nBVs = a+b+c+d+e+f
            end do
            nBVs = nBVs/3
            nBVscol = 7 
            
          case(2)               !Test 1 Houston L-domain
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux  = 0.0 
              uy  = 0.0 
              mult= 0.0
              if(y.eq.ymax) then
                if(x.eq.xmax)then
                  ux =-exp(x)*(y*cos(y) + sin(y))  !-3.756049227 
                  uy = exp(x)*y*sin(y)             !2.287355287
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, mult                       !Right top Corner
                  a = a+1
                elseif(x.eq.xmin)then
                  ux =-exp(x)*(y*cos(y) + sin(y)) ! 0.508325986 
                  uy = exp(x)*y*sin(y)            ! 2.287355287
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, mult                        !Left top Corner
                  a = a+1
                else
                  ux = -1.3817732*exp(x)
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 0, 1
                  write(300,20)   ux, uy, mult                        !Top boundary
                  a = a+1
                end if
                
              else if(y.eq.ymin)then
                if(x.eq.xmin)then
                  ux =-exp(x)*(y*cos(y) + sin(y))
                  uy = exp(x)*y*sin(y)          
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, mult                        !Left bottom Corner
                  b = b+1
                elseif(x.eq.xmiddle)then
                  ux =-exp(x)*(y*cos(y) + sin(y))
                  uy = exp(x)*y*sin(y)           
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, mult                        !central bottom Corner (0,-1)
                  b = b+1
                else
                  ux = 1.3817732*exp(x)
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 0, 1                           !Bottom boundary
                  write(300,20)   ux, uy, mult
                  b = b+1
                end if
                
              else if(x.eq.xmax)then
                if(y.gt.ymiddle)then
                  uy = 2.718281828*y*sin(y) 
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 0, 1, 1                           !Right boundary
                  write(300,20)   ux, uy, mult
                  c = c+1
                elseif(y.eq.ymiddle)then
                  ux =-exp(x)*(y*cos(y) + sin(y))
                  uy = exp(x)*y*sin(y)           
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 1, 1                           !middle right corner
                  write(300,20)   ux, uy, mult
                  c = c+1
                endif
                
              else if (x.eq.xmin)then
                  uy = 0.3678794412*y*sin(y) 
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 0, 1, 1                           !Left boundary
                  write(300,20)   ux, uy, mult
                  d = d+1
              else if(x.eq.xmiddle)then
                if(y .eq. ymiddle)then
                  ux =-exp(x)*(y*cos(y) + sin(y))
                  uy = exp(x)*y*sin(y)           
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 1, 1                           !central corner  (0,0)
                  write(300,20)   ux, uy, mult
                  d = d+1
                elseif(y .gt. ymin.and.(y.lt.0.0))then
                  uy = y*sin(y) 
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 0, 1, 1                    !central vertical boundary (x=0)
                  write(300,20)   ux, uy, mult
                  d = d+1
                end if
               
              else if(y.eq.ymiddle)then
                if(x .gt. xmiddle)then
                  ux =-exp(x)*(y*cos(y) + sin(y))
                  mult = sin(pi*(x-1.0)/2.0) * sin(pi*(y-1.0)/2.0)
                  write(200,10) i, 1, 0, 1
                  write(300,20)   ux, uy, mult                 !central horizontal boundary
                  e = e+1
                end if
               
              end if
              
              nBVs = a+b+c+d+e
            end do
            !nBVs = nBVs/3
            nBVscol = 7 
            
          case(3)               !Cuadratic function
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              uy = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                       !Upper Right Corner
                  ux = (x-x**2)*(1.0-2.0*y)
                  uy =-(1.0-2.0*x)*(y-y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                elseif(x.eq.xmin)then                  !Upper Left Corner
                  ux = (x-x**2)*(1.0-2.0*y)
                  uy =-(1.0-2.0*x)*(y-y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                else
                  ux = -(x - x**2)
                  uy = 0.0                              !Upper Boundary 
                  write(200,10) i, 1, 0, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                end if
                
              else if(y.eq.ymin)then
                if(x.eq.xmin)then
                  ux = (x-x**2)*(1.0-2.0*y)
                  uy =-(1.0-2.0*x)*(y-y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                elseif(x.eq.xmax)then
                  ux = (x-x**2)*(1.0-2.0*y)
                  uy =-(1.0-2.0*x)*(y-y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                else
                  ux = (x - x**2)
                  uy = 0.0                              !Down Border
                  write(200,10) i, 1,  0, 1
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                end if
                
              else if(x.eq.xmax)then
                ux = 0.0                                !Right Boundary Corner
                uy = (y - y**2)
                write(200,10) i, 0,  1, 1
                write(300,20)   ux, uy, 0.0
                c = c+1
              else if (x.eq.xmin)then
                ux = 0.0                                !Left Boundary Corner
                uy = -(y-y**2)
                write(200,10) i, 0,  1, 1
                write(300,20)   ux, uy, 0.0
                d = d+1
              end if
              
              nBVs = a+b+c+d
            end do
            nBVscol = 7 
            
          case(4)
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              uy = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                          !Upper Right Corner
                  ux = x**2(1.0-x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2*x**3 -3*x**2 +x)*y**2(1.0-y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                elseif(x.eq.xmin)then                      !Upper Left Corner
                  ux = x**2(1.0-x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2*x**3 -3*x**2 +x)*y**2(1.0-y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                else
                  ux = -(x - x**2)
                  uy = 0.0                                 !Upper Boundary 
                  write(200,10) i, 1, 0, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                end if
                
              else if(y.eq.ymin)then
                if(x.eq.xmin)then
                  ux = x**2(1.0-x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2*x**3 -3*x**2 +x)*y**2(1.0-y**2)
                  write(200,10) i, 1, 1, 1                 !Down Left Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                elseif(x.eq.xmax)then
                  ux = x**2(1.0-x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2*x**3 -3*x**2 +x)*y**2(1.0-y**2)
                  write(200,10) i, 1, 1, 1                 !Down Right Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                else
                  ux =
                  uy = 
                  write(200,10) i, 1,  0, 1                !Down Border
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                end if
                
              else if(x.eq.xmax)then
                ux = 
                uy =
                write(200,10) i, 0,  1, 1                  !Right Boundary
                write(300,20)   ux, uy, 0.0
                c = c+1
              else if (x.eq.xmin)then
                ux =
                uy = 
                write(200,10) i, 0,  1, 1                  !Left Boundary
                write(300,20)   ux, uy, 0.0
                d = d+1
              end if
              
              nBVs = a+b+c+d
            end do
            nBVscol = 7 
          case Default
            print*,'From BVS, case 3 and 4 are for 1 DoF'
            stop
        end select
       
      end if
      
      close(100)
      
      
      print*, 'No. of Bound Val', nBVS 
      
      
     
      10 format(I6,2x,3(1x,I2))
      20 format(3(e15.7,2x))
      30 format(f12.5)
      40 format(2(f12.5,2x))
      
    end subroutine SetBoundVal 
   
    
    
  !end contains 

end module BoundVal

