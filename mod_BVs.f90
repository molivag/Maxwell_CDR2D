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
      integer :: ierror, a ,b, c, d, i
      real    :: x, y, xmin, xmax, ymin, ymax
      integer, intent(out) :: nBVs, nBVscol
      
      
      open(unit=100, file=fileplace//'BVs.dat',Status= 'replace', action= 'write',iostat=ierror)
      
      a = 0
      b = 0
      c = 0
      d = 0
      
      xmin = minval(coord(:,2)) !the smallest number in x column
      xmax = maxval(coord(:,2)) !the greatest number in x column
      ymin = minval(coord(:,3)) !the smallest number in y column
      ymax = maxval(coord(:,3)) !the greatest number in y column
      
      
      !print*, ' '
      !print*, 'xmin= ', xmin
      !print*, 'xmax= ', xmax
      !print*, 'ymin= ', ymin
      !print*, 'ymax= ', ymax
      !print*, ' '
      
      
      if(ndofn .eq. 3) then
        do i =1, nnodes
          x=coord(i,2)
          y=coord(i,3)
          if(y.eq.ymax) then !top edge 
            write(100,50) i, 1,1,1, real(1), real(0), real(1)
            a = a+3
          else if (y.eq.ymin)then !botomm edge
            write(100,50) i, 1,1,1, real(0), real(0), real(0)
            b = b+3
          else if (x.eq.xmin)then !left edge
            write(100,50) i, 1,1,1, real(0), real(0), real(0)
            c = c+3
          else if(x.eq.xmax)then !right edge
            write(100,50) i, 1,1,1, real(0), real(0), real(0)
            d = d+3
          end if
          nBVs = a+b+c+d
        end do
          nBVs = nBVs/3
          nBVscol = 7 
        
      elseif(ndofn .eq. 2) then
        do i =1, nnodes
          x=coord(i,2)
          y=coord(i,3)
          if(y.eq.ymax) then                          !top edge
            write(100,60) i, 1,1, real(0), real(0)
            a = a+2
          else if (y.eq.ymin)then                     !bottom edge
            write(100,60) i, 1,1, real(0), real(0)
            b = b+2
          else if ((x.eq.xmin) .and. (y.eq.ymin))then !left down corner 
            write(100,60) i, 1,1, real(1), real(1)
          else if ((x.eq.xmin) .and. (y.eq.ymax))then !left up corner 
            write(100,60) i, 1,1, real(1), real(1)
          else if (x.eq.xmin )then                    !left edge
            write(100,60) i, 1,1, real(1), real(1)
            c = c+2
          else if (x.eq.xmax)then                     !right edges
            write(100,60) i, 1,1, real(0), real(0)
            d = d+2
          end if
          nBVs = a+b+c+d
        end do
          nBVs = nBVs/2
          nBVscol = 5
        
      elseif(ndofn .eq. 1)then
        do i =1, nnodes
          x=coord(i,2)
          y=coord(i,3)
          if(y.eq.ymax) then                          !top edge
            write(100,70) i, 1, real(0)
            a = a+1
          else if (y.eq.ymin)then                     !botomm
            write(100,70) i, 1, real(0)
            b = b+1
          else if ((x.eq.xmin) .and. (y.eq.ymax))then !left up corner 
            write(100,70) i, 1, real(1)
          else if (x.eq.xmin)then                     !left edge
            write(100,70) i, 1, real(0)
            c = c+1
          else if ((x.eq.xmin) .and. (y.eq.ymin))then !left down corner 
            write(100,70) i, 1, real(1)
          else if (x.eq.xmax)then                     !right edge
            write(100,70) i, 1, real(0)
            d = d+1
          end if
          nBVs = a+b+c+d
        end do
          nBVscol = 3
      end if
      
      close(100)
      
      
      50 format(4I6,3f10.3)
      60 format(3I6,2f10.3)
      70 format(2I6,f10.3)
      
      
    end subroutine SetBoundVal 
   
    
    
  !end contains 

end module BoundVal

  
