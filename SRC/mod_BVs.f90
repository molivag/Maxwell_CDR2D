module BoundVal
  use param
  use geometry

  contains
    
    subroutine SetBoundVal( boundaryPoints, dataColumns )
      !========================================================================
      !Esta subroutina revisa todos los nodos de la malla y define el tipo de Dof de
      !los nodos en la frontera. Abre un archivo en donde comenzara a escribir, 
      !en la primer columna: el numero de nodo. 
      ! La segunda columna tendra el tipo de nodo
      ! 1 = Nodo con valor preescrito (nodo con valor definido)
      ! 0 = Nodo libre (su valor debera ser calculado) eq to neumann integral = 0
      ! La tercera columna tendra el valor del nodo preescrito  
      !=========================================================================
      
      implicit none
      
      double precision, parameter :: pi = 4*atan(1.d0)
      character(len=*), parameter :: fileplace ="./"
      
      integer              :: ierror, a ,b, c, d, e, f,i 
      double precision     :: x, y, xmin, xmax, ymin, ymax, xmiddle, ymiddle
      double precision     :: ux, uy, exRe, eyRe, ezRe, exIm, eyIm, ezIm, pRe, pIm 
      double precision     :: aa, cc
      integer, intent(out) :: boundaryPoints, dataColumns
      
      
      open(unit=200, file=fileplace//'ifpre.dat',Status= 'replace', action= 'write',iostat=ierror)
      open(unit=300, file=fileplace//'BoVal.dat',Status= 'replace', action= 'write',iostat=ierror)
      
      a = 0; b = 0; c = 0;
      d = 0; e = 0; f = 0
      
      xmin = minval(coord(1,:)) !the smallest number in x column
      xmax = maxval(coord(1,:)) !the greatest number in x column
      ymin = minval(coord(2,:)) !the smallest number in y column
      ymax = maxval(coord(2,:)) !the greatest number in y column
      if(TwoHalf == 'Y')then
        xmiddle =  (abs(xmax)-abs(xmin))/2.0
        ymiddle =  ((ymax)-(ymin))/2.0
      else
        xmiddle =  (abs(xmax)-abs(xmin))/2.0
        ymiddle =  (abs(ymax)-abs(ymin))/2.0
      endif
      
      print*, ' '
      print*, '!================ Domain Dimensions ===========!'
      write(*,'(A)') ' -    xmin ,   xmax'
      write(*,'(A3,1x,I0,A,2x,I0)')' ', int(xmin),'  ,', int(xmax) 
      write(*,'(A9,I0)') ' - ymax = ', int(ymax)
      write(*,'(A9,I0)') ' - ymin = ', int(ymin)
      write(*,'(A9,f7.2)') ' - xhlf = ', real(xmiddle)
      write(*,'(A9,f7.2)') ' - yhlf = ', real(ymiddle)
      
      
      if(ndofn.eq.1) then
        print*, 'scalar Boundary Conditions ', ndofn 
        select case(BCsProb)
          case(6)
            print*, '!Resistivity'
              
            do i =1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                   !right top corner 
                  write(200,10) i, 0
                  write(300,30) ux 
                elseif(x.eq.xmin)then               !left top corner 
                  write(200,10) i, 0
                  write(300,30) ux 
                else
                  write(200,10) i, 0                !Top boundary
                  write(300,30) ux 
                end if
                a = a+1
              else if (y.eq.ymin)then
                if(x.eq.xmin)then                   !left bottom corner  
                  write(200,10) i, 1
                  write(300,30) ux
                elseif(x.eq.xmax)then               !right bottom corner
                  write(200,10) i, 1
                  write(300,30) ux
                else
                  write(200,10) i, 1                !bottom boundary
                  write(300,30) ux
                end if
                b = b+1
              else if(x.eq.xmax)then                !right boundary
                  write(200,10) i, 1
                  write(300,30) ux
                c = c+1
              else if (x.eq.xmin)then               !left boundary
                  write(200,10) i, 1
                  write(300,30) ux
                d = d+1
               
              end if
              boundaryPoints = a+b+c+d
             
            end do
            
          case(7)
            print*, ' !Double Line Source'
            do i =1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                   !right top corner 
                  write(200,10) i, 1
                  write(300,30) ux 
                elseif(x.eq.xmin)then               !left top corner 
                  write(200,10) i, 1
                  write(300,30) ux 
                else
                  write(200,10) i, 1                !Top boundary
                  write(300,30) ux 
                end if
                a = a+1
              else if (y.eq.ymin)then
                if(x.eq.xmin)then                   !left bottom corner  
                  write(200,10) i, 1
                  write(300,30) ux
                elseif(x.eq.xmax)then               !right bottom corner
                  write(200,10) i, 1
                  write(300,30) ux
                else
                  write(200,10) i, 1                !bottom boundary
                  write(300,30) ux
                end if
                b = b+1
              else if(x.eq.xmax)then                !right boundary
                  write(200,10) i, 1
                  write(300,30) ux
                c = c+1
              else if (x.eq.xmin)then               !left boundary
                  write(200,10) i, 1
                  write(300,30) ux
                d = d+1
               
              end if
              boundaryPoints = a+b+c+d
             
            end do
            
        endselect
        dataColumns = 3
        
        
      elseif(ndofn .eq. 2) then
        do i =1, nnodes
          x=coord(1,i)
          y=coord(2,i)
          if(y.eq.ymax) then 
            if(x.eq.xmax)then                    !right top corner         
              write(200,10) i, 1, 1                                    
              write(300,40)   1.0, 0.0                                 
            elseif(x.eq.xmin)then                !left top corner          
              write(200,10) i, 1, 1                                    
              write(300,40)   1.0, 0.0                                 
            else                                                       
              write(200,10) i, 1, 1              !Top boundary
              write(300,40)   1.0, 0.0                                 
            end if                                                     
            a = a+2                                                    
          else if (y.eq.ymin)then                                      
            if(x.eq.xmin)then                    !left bottom corner       
              write(200,10) i, 1, 1                                    
              write(300,40)   0.0, 0.0                                 
            elseif(x.eq.xmax)then                !right bottom corner       
              write(200,10) i, 1, 1                                    
              write(300,40)   0.0, 0.0                                 
            else                                                       
              write(200,10) i, 1, 1              !bottom boundary
              write(300,40)   0.0, 0.0                                 
            end if                                                     
            b = b+2                                                    
          else if(x.eq.xmax)then                 !right boundary           
              write(200,10) i, 1, 1                                    
              write(300,40)   0.0, 0.0                                 
            c = c+2                                                    
          else if (x.eq.xmin)then                !left boundary            
              write(200,10) i, 1, 1
              write(300,40)   0.0, 0.0
            d = d+2
           
          end if
          
          boundaryPoints = a+b+c+d+e
         
        end do
        boundaryPoints = nBVs/2
        dataColumns = 5
        
      elseif(ndofn .eq. 3)then
        select case(BCsProb)
          case(1)       !---------------------------------------------> Singular Solution
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
                e = e+2
               
              else if(y.eq.ymiddle)then
                if(x .gt. xmiddle)then
                  ux = 0.0
                  write(200,10) i, 1, 0, 1
                  write(300,20)   ux, uy, 0.0                 !horizontal boundary at the middle
                end if
                
                f = f+1
              end if
              
              boundaryPoints = a+b+c+d+e+f
            end do
            boundaryPoints = nBVs/3
            dataColumns = 7 
            
          case(2)               !Maxwell Problem nxE = 0 
            print*,'Maxwel BVs'
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              uy = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                      !Upper Right Corner
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                elseif(x.eq.xmin)then                  !Upper Left Corner
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                else                                    !Upper Border
                  write(200,10) i,  1, 0, 1  
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                end if
               
              else if(y.eq.ymin)then
                if(x.eq.xmin)then                     !left bottom corner  
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                elseif(x.eq.xmax)then                 !right bottom corner
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                else                                  !Down Border
                  write(200,10) i, 1,  0, 1
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                end if
                
              else if(x.eq.xmax)then                  !Right Boundary
                write(200,10) i, 0,  1, 1
                write(300,20)   ux, uy, 0.0
                c = c+1
              else if (x.eq.xmin)then                 !Left Boundary
                write(200,10) i, 0,  1, 1
                write(300,20)   ux, uy, 0.0
                d = d+1
              end if
              
              boundaryPoints = a+b+c+d
            end do
            dataColumns = 7 
            
          case(3)              !Maxwelll manufacture solution
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              uy = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                          !Upper Right Corner
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                elseif(x.eq.xmin)then                      !Upper Left Corner
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                else
                  ux = 0.0
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 0, 1                 !Upper Boundary 
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                end if
                
              else if(y.eq.ymin)then
                if(x.eq.xmin)then
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1                 !Down Left Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                elseif(x.eq.xmax)then
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1                 !Down Right Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                else                                             !Down Border
                  ux = 0.0
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1,  0, 1
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                end if
                
              else if(x.eq.xmax)then                              !Right Boundary
                ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                uy = 0.0 
                write(200,10) i, 0,  1, 1
                write(300,20)   ux, uy, 0.0
                c = c+1
              else if (x.eq.xmin)then                             !Left Boundary
                ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                uy = 0.0                                
                write(200,10) i, 0,  1, 1
                write(300,20)   ux, uy, 0.0
                d = d+1
              end if
              
              boundaryPoints = a+b+c+d
            end do
            dataColumns = 7 
            
          case(4)       !Stokes manufacture solution
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              uy = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                          !Upper Right Corner
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                elseif(x.eq.xmin)then                      !Upper Left Corner
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                else
                  ux = 0.0
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1                 !Upper Boundary 
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                end if
                
              else if(y.eq.ymin)then
                if(x.eq.xmin)then
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1                 !Down Left Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                elseif(x.eq.xmax)then
                  ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1                 !Down Right Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                else
                  ux = 0.0
                  uy =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(-2.0*y+1.0+y**2)
                  write(200,10) i, 1, 1, 1                !Down Border
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                end if
                
              else if(x.eq.xmax)then                       !Right Boundary
                ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                uy = 0.0 
                write(200,10) i, 1, 1, 1
                write(300,20)   ux, uy, 0.0
                c = c+1
              else if (x.eq.xmin)then                      !Left Boundary
                ux = x**2 * (-2.0*x+1.0+x**2)*(2*y**3 -3*y**2 +y)
                uy = 0.0                                
                write(200,10) i, 1,  1, 1
                write(300,20)   ux, uy, 0.0
                d = d+1
              end if
              
              boundaryPoints = a+b+c+d
            end do
            dataColumns = 7 
            
          case(5)       !Cavitty Driven Flow/Stokes
            print*,'Cavitty Driven Flow Boundary Conditions'
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              ux = 0.0 
              uy = 0.0 
              if(y.eq.ymax) then
                if(x.eq.xmax)then                    !Upper Right Corner
                  ux = 1.0
                  write(200,10) i, 1, 1, 0
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                elseif(x.eq.xmin)then                !Upper Left Corner
                  ux = 1.0
                  write(200,10) i, 1, 1, 0
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                elseif(x.eq.xmiddle)then             !Upper middle top preasure condition
                  print*,'centro', xmiddle
                  ux = 1.0
                  write(200,10) i, 1, 1, 1
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                else
                  ux = 1.0
                  write(200,10) i, 1, 1, 0           !Upper Boundary 
                  write(300,20)   ux, uy, 0.0
                  a = a+1
                end if
                
              else if(y.eq.ymin)then
                if(x.eq.xmin)then
                  ux = 0.0
                  write(200,10) i, 1, 1, 0            !Down Left Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                elseif(x.eq.xmax)then
                  ux = 0.0
                  write(200,10) i, 1, 1, 0                 !Down Right Corner
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                else
                  ux = 0.0
                  write(200,10) i, 1, 1, 0                !Down Boundary
                  write(300,20)   ux, uy, 0.0
                  b = b+1
                end if
                
              else if(x.eq.xmax)then
                ux = 0.0  
                write(200,10) i, 1, 1,  0                 !Right Boundary
                write(300,20)   ux, uy, 0.0
                c = c+1
              else if (x.eq.xmin)then
                ux = 0.0
                write(200,10) i, 1,  1, 0                  !Left Boundary
                write(300,20)   ux, uy, 0.0
                d = d+1
              end if
              
              boundaryPoints = a+b+c+d
            end do
            dataColumns = 7 
            
            
            
          case Default
            print*,'From BVS, case 3 and 4 are for 1 DoF'
            stop
        end select
      elseif(ndofn .eq. 8)then
            print*,'Maxwel BVs for 2.5D'
            exRe = 0.0 
            eyRe = 0.0 
            ezRe = 0.0 
            exIm = 0.0 
            eyIm = 0.0 
            ezIm = 0.0 
            pRe  = 0.0 
            pIm  = 0.0 
            
            do i = 1, nnodes
              x=coord(1,i)
              y=coord(2,i)
              if(y.eq.ymax) then
                if(x.eq.xmax)then                      !Upper Right Corner
                  write(200,11) i, 1,    1,    1,    1,   1,    1,    1,     1
                  write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm 
                  a = a+1
                elseif(x.eq.xmin)then                  !Upper Left Corner
                  write(200,11) i,  1,    1,    1,    1,  1,    1,    1,    1
                  write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm 
                  a = a+1
                else                                    !Upper Border
                  write(200,11) i, 1,    0,    0,    1,   1,    0,    0,    1
                  write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm  
                  a = a+1
                end if
               
              else if(y.eq.ymin)then
                if(x.eq.xmin)then                     !left bottom corner  
                  write(200,11) i,  1,    1,    1,    1,    1,    1,   1,   1
                  write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm 
                  b = b+1
                elseif(x.eq.xmax)then                 !right bottom corner
                  write(200,11) i,  1,    1,    1,    1,    1,    1,   1,   1
                  write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm
                  b = b+1
                else                                  !Down Border
                  write(200,11) i, 1,    0,    0,    1,   1,    0,    0,    1   
                  write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm 
                  b = b+1
                end if
                
              else if(x.eq.xmax)then                  !Right Boundary
                write(200,11) i, 0,    0,    1,    1,   0,    0,    1,    1   
                write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm 
                c = c+1
              else if (x.eq.xmin)then                 !Left Boundary
                write(200,11) i, 0,    0,    1,    1,   0,    0,    1,    1   
                write(300,21)   exRe, eyRe, ezRe, pRe, exIm, eyIm, ezIm,  pIm 
                d = d+1
              end if
              
              boundaryPoints = a+b+c+d
            end do
            dataColumns = 17 
       
      else
        print*,"- Degrees of Freedom exceeded while setting Boundary COnditions"
        print*, ' ' 
        print*, '>>>>> Program stopped' 
        print*, ' ' 
        stop
      end if
      close(200)
      close(300)
      
      
      !print*, 'Num. of Boundary nodes', nBVS 
      
      10 format(I6,2x,3(1x,I2))
      11 format(I6,2x,9(1x,I2))
      20 format(3(e15.7,2x))
      21 format(9(e15.7,2x))
      30 format(f12.5)
      40 format(2(f12.5,2x))
     
    end subroutine SetBoundVal 
    
    
  !end contains 

end module BoundVal

