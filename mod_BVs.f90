module BoundVal
  use param

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
      
      character(len=*), parameter :: fileplace ="./"
      integer :: ierror, a ,b, c, d, e, f,i 
      double precision  :: x, y, xmin, xmax, ymin, ymax, xhalf, yhalf, ux, uy
      integer, intent(out) :: nBVs, nBVscol
      
      
      open(unit=200, file=fileplace//'ifpre.dat',Status= 'replace', action= 'write',iostat=ierror)
      open(unit=300, file=fileplace//'BoVal.dat',Status= 'replace', action= 'write',iostat=ierror)
      
      a = 0
      b = 0
      c = 0
      d = 0
      e = 0
      f = 0
      
      !u(x) = u1/x1 - x0
      
      xmin = minval(coord(:,2)) !the smallest number in x column
      xmax = maxval(coord(:,2)) !the greatest number in x column
      ymin = minval(coord(:,3)) !the smallest number in y column
      ymax = maxval(coord(:,3)) !the greatest number in y column
      xhalf = xmax/2.0 !0.0
      yhalf = ymax/2.0 !0.0
      !print*, 'xmin = ', xmin
      !print*, 'xmax = ', xmax
      !print*, 'ymin = ', ymin
      !print*, 'ymax = ', ymax
      !print*, 'xhalf = ', xhalf
      !print*, 'yhalf = ', yhalf
      
      
      
      if(ndofn .eq. 1) then
        
        do i =1, nnodes
          x=coord(i,2)
          y=coord(i,3)
          if(y.eq.ymax) then 
            if(x.eq.xmax)then           
              write(100,70) i, 1, real(1)   !right top corner
            elseif(x.eq.xmin)then
              write(100,70) i, 1, real(0)   !left top corner
            else                
              write(100,70) i, 0, real(0)   !top edge
            end if
            a = a+1
          else if (y.eq.ymin)then
            if(x.eq.xmin)then
              write(100,70) i, 1, real(0)   !left bottom corner
            elseif(x.eq.xmax)then                                
              write(100,70) i, 1, real(1)   !right bottom corner
            else
              write(100,70) i, 0, real(0)   !botomm edge
            end if
            b = b+1
          else if(x.eq.xmax)then
              write(100,70) i, 1, real(1)   !right edge
            d = d+1
          else if (x.eq.xmin)then
              write(100,70) i, 1, real(0)   !left edge
            c = c+1
           
          end if
          nBVs = a+b+c+d
         
        end do
        nBVscol = 3
        
        
      elseif(ndofn .eq. 2) then
        do i =1, nnodes
          x=coord(i,2)
          y=coord(i,3)
          if(y.eq.ymax) then 
            if(x.eq.xmax)then                         
              write(100,60) i, 1,1, real(0), real(0)     !right top corner 
            elseif(x.eq.xmin)then                     
              write(100,60) i, 1,1, real(0), real(0)     !left top corner
            else
              write(100,60) i, 1,1, real(0), real(0)     !top edge 
            end if
            a = a+2
          else if (y.eq.ymin)then
            if(x.eq.xmin)then                         
              write(100,60) i, 1,1, real(0), real(0)    !left bottom corner 
            elseif(x.eq.xmax)then                      
              write(100,60) i, 1,1, real(0), real(0)    !right bottom corner
            else
              write(100,60) i, 1,1, real(0), real(0)    !botomm edge
            end if
            b = b+2
          else if(x.eq.xmax)then                      
              write(100,60) i, 1,1, real(0), real(0)    !right edge
            c = c+2
          else if (x.eq.xmin)then                     
              write(100,60) i, 1,1, real(0), real(0)    !left edge
            d = d+2
           
          else if(x.eq.xhalf)then
            !print*,'coordinate', coord(i,1)
            !print*,' ' 
            !print*, 'x  = ', x
            !
            !print*, 'y_ = ', y
            !print*, ' ' 
            
            if(y .gt. ymin)then
              !print*, 'y = ', y
              write(100,60) i, 1,1, real(0), real(0)     !half top edge
              e = e+2
            end if
            
          else if(y.eq.yhalf)then
            !print*,'coordinate', coord(i,1)
            !print*,' ' 
            !print*, 'x  = ', x
            !
            !print*, 'y_ = ', y
            !print*, ' ' 
            
            if(x .lt. xmax)then
              !print*, 'y = ', y
              write(100,60) i, 1,1, real(0), real(0)     !half top edge
              e = e+2
            end if
            
          end if
          
          
          !end if
          nBVs = a+b+c+d+e
         
        end do
        nBVs = nBVs/2
        nBVscol = 5
        
      elseif(ndofn .eq. 3)then
        do i = 1, nnodes
          x=coord(i,2)
          y=coord(i,3)
          
          if(y.eq.ymax) then
            if(x.eq.xmax)then
              !ux =2.0!1.0 + y
              !uy =1.0                !---->  Lineal malla 0-1
              ux = 0.0
              uy = (y - y**2)        !---->  malla 0-1
              !ux = -90.*(1.-2.*y)
              !uy =19.*(y-y**2)       !---->  malla 0-10
              write(200,10) i, 1, 1, 1
              write(300,20)   ux, uy, 0.0                                 !right top corner  
            elseif(x.eq.xmin)then                                                            
              !ux = 2.0!1.0 + y                                                              
              !uy = 2.!2.0             !---->  Lineal malla 0-1                              
              ux = 0.0                                                                       
              uy = -(y-y**2)          !---->  malla 0-1                                      
              !ux = 0.0                                                                      
              !uy = -(y-y**2)          !---->  malla 0-10                                    
              write(200,10) i, 1, 1, 1                                                       
              write(300,20)   ux, uy, 0.0                                 !left topt corner
            else
              !ux = 2.0
              !uy = 2.0 - x            !---->  Lineal malla 0-1
              ux = -(x - x**2)
              uy = 0.0                !---->  malla 0-1
              !ux =-19*(x-x**2)
              !uy = 90*(1-2*x)         !---->  malla 0-10
              write(200,10) i, 1, 0, 1
              write(300,20)   ux, uy, 0.0                                 !top edge
            end if
            
            !if(x.eq.xhalf)then
            !  write(100,50) i, 1, 1, 1, real(0), real(0), real(0)        !half top edge
            !end if
            a = a+3
            
          else if(y.eq.ymin)then
            if(x.eq.xmin)then
              !ux = 1.!1.0! + y
              !uy = 2.!2.0              !---->  Lineal malla 0-1 
              ux = 0.0
              uy = -(y-y**2)          !---->  malla 0-1
              !ux = 0.0
              !uy = -(y-y**2)           !---->  malla 0-10
              write(200,10) i, 1, 1, 1
              write(300,20)   ux, uy, 0.0                                 !left bottom corner
            elseif(x.eq.xmax)then
              !ux = 1.!1.0! + y
              !uy = 1.0             !---->  Lineal malla 0-1 
              ux = 0.0
              uy = (y - y**2)         !---->  malla 0-1
              !ux =-90.*(1.-2.*y)
              !uy = 19.*(y-y**2)        !---->  malla 0-10
              write(200,10) i, 1, 1, 1
              write(300,20)   ux, uy, 0.0                                 !right bottom corner
            else
              !ux = 1.0
              !uy = 2.0 - x            !---->  Lineal malla 0-1 
              ux = (x - x**2)
              uy = 0.0                !---->  malla 0-1
              !ux = (x-x**2)
              !uy = 0.0                 !---->  malla 0-10
              write(200,10) i, 1,  0, 1
              write(300,20)   ux, uy, 0.0                                 !bottom edge
            end if
            b = b+3
            
          else if(x.eq.xmax)then
              !ux = 1.0 + y            !---->  Lineal malla 0-1 
              !uy = 1.0
              ux = 0.0                !---->  malla 0-1
              uy = (y - y**2)
              !ux =-90.*(1.-2.*y)       !---->  malla 0-10
              !uy = 19.*(y-y**2)
              write(200,10) i, 0,  1, 1
              write(300,20)   ux, uy, 0.0                                 !right edge
            c = c+3
            
          else if (x.eq.xmin)then
              !ux = 1.0 + y            !---->  Lineal malla 0-1 
              !uy = 2!2.0
              ux = 0.0                !---->  malla 0-1
              uy = -(y-y**2)
              !ux = 0.0                 !---->  malla 0-10
              !uy = -(y-y**2)
              write(200,10) i, 0,  1, 1
              write(300,20)   ux, uy, 0.0                                 !left edge
            d = d+3
           
          end if
          nBVs = a+b+c+d
         
        end do
        nBVs = nBVs/3
        nBVscol = 7 
        
      end if
      
      close(100)
      
      
      10 format(I6,1x,3(1x,I2))
      20 format(3(f12.5,2x))
     ! 30 format(I6,2x,3(1x,I2),3(f12.5,2x))
      
      60 format(I6,2x,2I2,2x,2f12.5)
      70 format(I6,2x,1I2,2x,f12.5)
      
    end subroutine SetBoundVal 
   
    
    
  !end contains 

end module BoundVal

