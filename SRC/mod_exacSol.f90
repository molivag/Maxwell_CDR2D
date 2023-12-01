module E0field
  use param
  use geometry
  
  integer                      :: id_poin

  contains
    subroutine Efield_WholeSpace(time, tEz, Ez_r)
     
      ! Implements eq. (17) of paper by Zhang & Liu (2021)
      
      character(len=*), parameter  :: fileplace3 = "Exact_Sol_TEM/2D_DoubleLine_WholeSpace/"
      double precision, parameter  :: pi = 4*atan(1.d0)
      double precision, intent(in) :: tEz                       !correspond to nt 
      integer         , intent(in) :: time                      !correspond to time 
      double precision             :: sigma, S, x1,y1,x2,y2,x,y
      double precision             :: Curr_x,Curr_y, angle_loop
      double precision             :: facto, theta2, rho1, rho2, exp1, exp2
      double precision             :: Ez_r_x, Ez_r_y
      integer                      :: itest, inode, ii, ipoin
      double precision             :: Mr,Mx,My, x_profile, mu
      double precision, intent(out):: Ez_r(ntotv)
      double precision, dimension(1,nnodes) :: xcor, ycor
     
      
      xcor  = spread(coord(1,:),dim = 1, ncopies= 1)
      ycor  = spread(coord(2,:),dim = 1, ncopies= 1)
      
      Ez_r = 0.0 
      mu = 1.0/lambda
      angle_loop = 90.0
      sigma = 0.01
      x1 = coord(1,Srcloc(1))
      y1 = coord(2,Srcloc(1))
      x2 = coord(1,Srcloc(2))
      y2 = coord(2,Srcloc(2))
      S  = abs(x1*x2) ! not needed, cancels out
      
      angle_loop=angle_loop*pi/180. ! theta=-30 deg
      Mr=Icurr(1)*S
      Mx=Mr*sin(angle_loop)
      My=Mr*cos(angle_loop)
      Curr_x= Mx/S
      Curr_y= My/S
      
      ! loop coord.
      facto =4.0*pi*tEz
      theta2=mu*sigma/(4.*tEz) ! = theta^2
      
      do inode = 1, nnodes  !simple Function
        x = coord(1,inode)
        y = coord(2,inode)
        
        rho1 = sqrt( (x-x1)**2 + (y-y1)**2 )
        rho2 = sqrt( (x-x2)**2 + (y-y2)**2 )
        exp1 = exp(-rho1*rho1*theta2)
        exp2 = exp(-rho2*rho2*theta2)
        ! Source x
        Ez_r_x = -mu*(Curr_x+Curr_y)/facto*exp1
        Ez_r_y = +mu*(Curr_x+Curr_y)/facto*exp2
        !Electric Field
        if(ndofn.eq.1)then
          Ez_r(inode) = Ez_r_x + Ez_r_y
        else
          Ez_r(ndofn*inode-2) = Ez_r_x + Ez_r_y
        endif
       
      end do
      
      !Este archivo se creara solo en el t=0 para no repetir en cada tiempo
      if (time.eq.0) then
        ! Aqui se crea el perfil espacial, identificasndo los nodos sobre la y de interes
        open(unit=10, file= fileplace3//"Id_spatial_profile.dat", ACTION="write", STATUS="replace")
        id_poin = 0! 113 
        do ipoin =1,nnodes
          if(ycor(1,ipoin).eq.0.0)then
            id_poin = id_Poin+1
            write(10,906) ipoin, xcor(1,ipoin)
          else
            continue
          endif
        end do
        close(10)
        !Aqui se escribe el archivo para GNUPLOT
        if(tEz.eq.time_ini)then
          open(unit=5, file= fileplace3//"Id_spatial_profile.dat", status='old', action='read')
          open(unit=6, file=fileplace3//"ini_exact_Ex.dat", ACTION="write", STATUS="replace")
          write(6,'(A3,A,e10.3,A)') ' "t','=',tEz,'"'
          do ii = 1,id_poin
            read(5,*) ipoin, x_profile
            write(6,918) ipoin, x_profile, Ez_r(ipoin) 
          end do
          close(5)
          close(6)
        end if
        
      endif
      
      !open(unit=5, file= fileplace3//"E0exact.dat", status='old', action='read')
      !write(5,1)x,Ez_r_x,Ez_r_y,Ez_r
      918 format(I7,3x,E14.6,3x,E14.6) !format for res velocity
      906 format(I7,2(3x,f9.4)) !format for msh
    end subroutine Efield_WholeSpace
   
end module E0field
