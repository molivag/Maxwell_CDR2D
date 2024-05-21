module E0field
  use param
  use geometry
  use library
  
  integer                      :: id_poin

  contains
    subroutine Efield_WholeSpace(time, tEz, Ez_r)
     
      ! Implements eq. (17) of paper by Zhang & Liu (2021)
      
      character(len=*), parameter  :: fileplace3 = "Exact_Sol_TEM/2D_DoubleLine_WholeSpace/"
      double precision, parameter  :: pi = 4*atan(1.d0)
      double precision, intent(in) :: tEz             !correspond to nt 
      integer         , intent(in) :: time            !correspond to time 
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
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine Res_Matlab(solution)
      
      implicit none
      
      double precision, parameter :: pi = 4*atan(1.d0)
      character(len=*), parameter :: fileplace = "Res/Exact_Solutions/"
      character(len=*), parameter :: fileplace2 = "Res/Geometry/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      double precision, dimension(ntotv, 1)             :: E_field_exac
      double precision, dimension(t_steps+1)            :: Ex_field
     
      double precision, dimension(1, ntotv) :: solution_T
      !double precision, dimension(1,nnodes) :: xcoor, ycoor
      !double precision, dimension(nnodes)   :: x, y
      !double precision, dimension(t_steps) :: t
      double precision, dimension(nnodes)    :: exact_y, exact_x, exact_p, FEM_x, FEM_y, FEM_p
      double precision, dimension(t_steps)   :: Texact_y, Texact_x, Texact_z, t
      double precision     :: aa, bb, cc, dd, ee, ds, sigma, mu, SrcCurr, r_vec, spi
      double precision     :: x, y, z, x1, y1, x2, y2
      double precision     :: theta,ex,ey,ez, nt, arg
      double precision     :: sum_error, error_EM, error_p, errL2_x, errL2_y, errL2_p
      double precision     :: x_FEM, y_FEM, p_FEM, uxSol, uySol, multi
      double precision     :: fi, psi, der_fi, der_psi 
      double precision     :: rho1, rho2 
      !double precision     :: SrcCurr, z, sigma, ds
      character(len=4)     :: extension!, File_Solution
      integer              :: ipoin, ielem, inode, i, itime, ii
      
      
      extension =".dat"
      !File_Solution ="globsolution"
      
      solution_T = transpose(solution)
      
      errL2_x=0.0; errL2_y=0.0; errL2_p=0.0
      error_EM = 0.0
      error_p  = 0.0
      sum_error = 0.0
      exact_x = 0.0
      exact_y = 0.0
      exact_p = 0.0
      uxSol = 0.0 
      uySol = 0.0
      multi = 1E-10
      
      select case(exacSol)
        case(0)
          print*,'No analytic solution'
        case(1)
          aa = (2.0/3.0)*n_val
          bb = 0.0
          cc = (n_val/3.0) - 1.0
          !cc = (n_val/3.0) - (1.0/2.0)
          dd = 0.0
          ee = 0.0 
          !write(*,*) '       FEMx','            Ex_x', '            FEMy','           Ex_y'
          do inode = 1, nnodes  
            x = coord(1,inode)
            y = coord(2,inode)
            
           !exact solution
            bb    = ((x**2 + y**2))
            
            
            if(abs(x).ne.0.01.and.abs(x).le.1.0e-4)then
              dd    = y/x
              ee    = atan(dd)
            else
              ee    = pi/2.0
            end if
              
            uxSol = aa * bb**cc * ( x*sin(aa*ee) - y*cos(aa*ee) )
            uySol = aa * bb**cc * ( y*sin(aa*ee) + x*cos(aa*ee) )
            !uxSol = aa * bb**cc * sin(aa*ee)
            !uySol = aa * bb**cc * cos(aa*ee)
            
            !FEM solution
            x_FEM = solution_T(1,ndofn*inode-2)
            y_FEM = solution_T(1,ndofn*inode-1)
            
            !write(*,"(4(f15.5,1x))") x_FEM, uxSol, y_FEM, uySol
            !error = error + (x_FEM - uxSol)**2 + (y_FEM - uySol)**2 
            error_EM = error_EM + ( (uxSol - x_FEM)**2  + (uySol - y_FEM)**2 )
            !print*, 'error', error  
            
            !Write to plotting file 
            exact_x(inode) = uxSol 
            exact_y(inode) = uySol
            
          end do
        case(2)
          ! Implements eq. (2.50) of Nabighian 1988 EM Methods (EM Theory book) (p. 175)
          !Nota: Esta solucion analitica determina el campo electrico en un punto de la malla (x,y)
          !para usarse de comparacion se requiere ejecutar el codigo y luego en el post-proceso 
          !extraer en un punto determinado ux e uy para todos los tiempos simulados
          
          ! Define variables
          ! I*ds = dipole moment, is set to I*ds=1
          SrcCurr  =  Icurr(1)
          ds       =  1.0
          sigma    =  1.0
          mu       =  1.0/lambda 
          
          nt  = time_ini
          arg = 0.0
          aa  = 0.0
          ee  = 0.0
          ex  = 0.0; ey = 0.0; ez = 0.0
          E_field_exac  = 0.0 
          ii = 0.0
          
          !delta_t 1e-3!( time_fin - time_ini ) / (t_steps + 1.0)  
          
          do i=1,t_steps 
            t(i) = nt 
            nt   = nt + delta_t
          end do
          
          spi  = sqrt(pi)
          
          write(*,*) ' Writing exact solution'
          do i=1,t_steps
            !print*, i
            do inode = 1, nnodes  
              x = coord(1,inode)
              y = coord(2,inode)
              z = 0.0
              
              r_vec= sqrt(x*x+y*y+z*z)
              cc   = SrcCurr*ds/(4.0*pi*sigma*r_vec**3)
              
              theta = sqrt(mu*sigma/(4.0*t(i)))
              aa    = 4.0/spi*theta**3*r_vec**3 + 6.0/spi*theta*r_vec
              arg   = -theta**2*r_vec**2
              ee    = erfc(theta*r_vec)
              aa    = aa*exp(arg)+3.0*ee
              bb    = 4.0/spi*theta**3*r_vec**3 + 2.0/spi*theta*r_vec
              bb    = bb*exp(arg)+ee
              
              ! geometry term
              ex    = cc * (aa*x**2/r_vec**2 - bb)
              ey    = cc * aa*x*y/r_vec**2
              ez    = cc * aa*x*z/r_vec**2
              !write(*,'(2x,I5,1x,3(e15.6))') inode, ex, ey, ez 
              E_field_exac(ndofn*inode-2,1) = ex
              E_field_exac(ndofn*inode-1,1) = ey
              E_field_exac(ndofn*inode-0,1) = ez
            end do
              !write(*,'(I5, e15.5, 3(E14.6))') i-1, t(i), E_field_exac( 
              call GID_PostProcess(0,2,E_field_exac, 'res', i-1, t(i), time_fin, Ex_field) 
              ii = ii+1.0
          end do
           
          !call GID_PostProcess(2,E_field_exac, 'msh', 0, t(1), time_fin, Ex_field) 
          
        case(3) !polynomial solution fo Maxwell (MVAF) and Stokes problem
          
          do inode = 1, nnodes  !simple Function
            x = coord(1,inode)
            y = coord(2,inode)
            
            !exact solution
            uxSol   = x**2 * (1.0-2.0*x + x**2)*(2*y**3 -3*y**2 +y)
            uySol   =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(1.0-2.0*y + y**2)
            
            !FEM solution
            x_FEM = solution_T(1,ndofn*inode-2)
            y_FEM = solution_T(1,ndofn*inode-1)             
            p_FEM = solution_T(1,ndofn*inode)
           
            error_EM = error_EM + ( (uxSol - x_FEM)**2  + (uySol - y_FEM)**2 )
            error_p  = error_p  + (multi - p_FEM )**2
            
            !Write to plotting file 
            exact_x(inode) = uxSol
            exact_y(inode) = uySol
            exact_p(inode) = multi
            FEM_x(inode) = solution_T(1,ndofn*inode-2)
            FEM_y(inode) = solution_T(1,ndofn*inode-1)
            FEM_p(inode) = solution_T(1,ndofn*inode-0) 
          end do
          
        case(4)
          print*,'Double Line Source exact solution'
          E_field_exac  = 0.0 
          SrcCurr = Icurr(1)
          sigma = 1.0
          mu    = 1.0/lambda 
          ez    = 0.0
          ii    = 0.0
          nt    = time_ini
          
          do i=1,t_steps 
            t(i) = nt 
            nt   = nt + delta_t
          end do
          
          do i = 1, t_steps
            do inode = 1,nnodes
              x = coord(1,inode)
              y = coord(2,inode)
              theta= sqrt((mu*sigma/4.0*t(i)))
              x1 = coord(1,Srcloc(1)) 
              y1 = coord(2,Srcloc(1))
              x2 = coord(1,Srcloc(2))
              y2 = coord(2,Srcloc(2))
              rho1 = sqrt( (x - x1)**2 + (y-y1)**2 )
              rho2 = sqrt( (x - x2)**2 + (y-y2)**2 )
              aa = rho1**2
              bb = rho2**2
              
              ez = -(mu*SrcCurr/4.0*sigma*t(i)) * ( exp(aa*theta**2) + exp(bb*theta**2) )
              E_field_exac(ndofn*inode-2,1) = ez
            end do
              call GID_PostProcess(0,2,E_field_exac, 'res', i-1, t(i), time_fin, Ex_field) 
              ii = ii+1.0
          end do
         
      end select
     
      error_EM = sqrt(error_EM/nnodes)
      error_p = sqrt(error_p/nnodes)
     
      errL2_x=norm2(exact_x - FEM_x)/norm2(exact_x)
      errL2_y=norm2(exact_y - FEM_y)/norm2(exact_y)
      errL2_p=norm2(exact_p - FEM_p)/norm2(exact_p)
      
      ! = = = = = = Begin Error computations lines
      !== Lines to write the error computation in a file to be load it in matlab
      open(unit=777, file= fileplace//error_name//extension, ACTION="write", STATUS="replace")
      write(777,"(1x,E15.5,3x, A)") error_EM,'%error in electric field'
      write(777,"(1x,E15.5,3x, A)") error_p, '%error in multiplier'
      write(777,"(1x,E15.5,3x, A)") errL2_x,'%L2error in ex'
      write(777,"(1x,E15.5,3x, A)") errL2_y,'%L2error in ey'
      write(777,"(1x,E15.5,3x, A)") errL2_p,'%L2error in multiplier'
      close(777)
      print*, ' '
      print*, '!============== Error Estimation ==============!'
      !write(*,"(A10,f7.5,A25,E13.5)")' -For h = ', helem, 'the error estimation is '
      write(*,"(A14,E13.5)")' -For u     : ', error_EM
      write(*,"(A14,E13.5)")' -For p     : ', error_p
      write(*,"(A14,E13.5)")' -Norm L2 ux: ', errL2_x
      write(*,"(A14,E13.5)")' -Norm L2 uy: ', errL2_y
      write(*,"(A14,E13.5)")' -Norm L2 p : ', errL2_p
      print*, ' '
      ! = = = = = = End Error computations lines
      
      
      ! = = = = = = Begin writing FEM and Exact solutions
      open(unit=111, file= fileplace//File_Nodal_Vals//extension, ACTION="write", STATUS="replace")
      !write(*,*) '       FEMx','            Ex_x', '            FEMy','           Ex_y'
     
      if(ProbType.eq.'TIME')then
        print*, 'xxxxxxxxxxxxxxx'
        if(ndofn.eq.1)goto 10
        write(111,'(A)')'%  Time         step         ex             ey             ez'
        do itime = 1, t_steps
          write(111,908) itime-1, t(itime), Texact_x(itime), Texact_y(itime), Texact_z(itime)
        end do
      else
        !ESTO ESTA MAL Y HAY QUE ARREGLARLO
        if(exacSol.eq.2)goto 10 
        if(ndofn.eq.3)then
          do ipoin = 1, nnodes  !   uh_x    uh_y    uex_x   uex_y
            write(111,'(A)') '%      FEM_x             FEM_y             Exact_x           Exact_y'
            write(111,906)&
            &solution_T(1,ndofn*ipoin-2),solution_T(1,ndofn*ipoin-1), exact_x(ipoin), exact_y(ipoin)
          end do
          
        elseif(ndofn.eq.1)then 
          do ipoin = 1, nnodes  !   uh_x    uh_y    uex_x   uex_y
            write(111,'(A)') '%      FEM               Exact_x'
            write(111,906) solution_T(1, ipoin), exact_x(ipoin)
          end do
        else
          print*, 'In Res_Matlab, Problem type not defined'
          stop
        end if
        10 continue
      end if
      
      print*, ' '
      print*, '!====== Matlab file ======'
      write(*,"(A7,A16,A23,A)") ' -File ',File_Nodal_Vals//extension,'written succesfully in ',fileplace
      close(111)
      ! = = = = = = End writing FEM and Exact solutions
      
      ! = = = = = = Begin write geometry
      open(unit=444, file= fileplace2//coord_name//extension, ACTION="write", STATUS="replace")
      open(unit=333, file= fileplace2//conec_name//extension, ACTION="write", STATUS="replace")
      
      do ielem=1,nelem
        write(333,902) (lnods(ielem,inode),inode=1,nne)
      end do
      
      do ipoin = 1, nnodes
        write(444,904) ipoin, coord(1,ipoin), coord(2,ipoin)
      end do
      write(*,903)&
      ' -File ',coord_name//extension,' and ',conec_name//extension,'written succesfully in ',fileplace2
      print*, ' '
      close(333)
      close(444)
      ! = = = = = = End write geometry
      
     
      902 format(1x,i5,10(1x,i5))
      903 format(A7,A16,A5,A16,A23,A)
      904 format(I7,2(3x,f9.4) ) !format for msh
      906 format(6(E15.5, 3x))
      908 format(I5,1x, F15.5, 3(E15.6))
      
      115 continue
    end subroutine Res_Matlab
   
end module E0field
