program DebyeEXP

  implicit none
  
  character(len=*), parameter :: PLT_FILE = 'plot.plt' ! Gnuplot file.
  integer         , parameter :: nfrec=41

  real              ::sigma_inf, sigma_0, m, tau, c, gamaP, tauP
  !REAL :: gammaIM=0.0001 !Factor de regularización imaginarios
  !REAL :: gammaRE=0.0002 !Factor de regularización reales
  !REAL ::frecmin, frecmax

  double precision              :: ti(1:nfrec), w(1:nfrec)
  double precision              :: ColeCole_real(1:nfrec), ColeCole_imag(1:nfrec), DebyeRe(1:nfrec), DebyeIM(1:nfrec)
  double complex                :: ColeCole(1:nfrec), Debye(1:nfrec)
  integer                       :: i, p!, k, L
  !REAL, allocatable,dimension(:) :: frec
  !REAL,allocatable,dimension(:) :: w

  sigma_inf= 0.125
    sigma_0= 0.05
        tau= 0.125
          m= (sigma_inf - sigma_0)/ sigma_inf
          c= 0.8

     gamaP = 0.7
      tauP = 0.9



  print*, 'Parametros Cole-Cole'
  print*
  print*, 'Conductividad a frec. maxima',sigma_inf
  print*, 'Conductividad a frec. cero..',sigma_0
  print*, 'Cargabilidad................',m
  print*, 'Constante de tiempo.........',tau
  print*, 'Dependencia de la frecuencia',c
  print*


  !Muestras
  ti(1)=-1.0
  do i=2,nfrec
    ti(i)= ti(i-1) + 0.1
    w(i-1)= 1*10**ti(i-1)
  end do
  
  !Frecuencias de muestreo
  do i=1,nfrec
    w(i)= 1*10**ti(i)
  end do

  !Modelo Cole-Cole
  do i= 1,nfrec
    ColeCole(i)= (sigma_inf*( 1 - ( m /( 1+(cmplx(0,1)*w(i)*tau)**c  ) ) ) )
  end do

  print *,'Muestras                  Frecuencias de muestreo:'
  do i=1,nfrec
    write(*,'(1x,I3,2x,2(f15.5))') i, ti(i), w(i)
  end do

  do i=1,nfrec
    ColeCole_real(i)= real(ColeCole(i))
    ColeCole_imag(i)=aimag(ColeCole(i))
  end do


  !Modelo Debye
  do i =1,nfrec
    Debye(i) = sigma_inf + (sigma_0-sigma_inf) * (gamaP/(1+(cmplx(0,1)*w(i)*tauP) ) )
  end do

  do i=1,nfrec
    DebyeRe(i) = real(Debye(i))
    DebyeIM(i) = aimag(Debye(i))
  end do


   100 format (e15.6, 5x,7(e15.6))
  open(unit=1, file="ColeCole.dat")
  write(1,*)'freq CCre CCim DeRE DeIM'
  do i=1,nfrec
    write(unit=1, fmt=100) w(i), ColeCole_real(i),ColeCole_imag(i), DebyeRe(i), DebyeIM(i)
  end do
  close(unit=1)

  call execute_command_line('gnuplot ' // PLT_FILE)
  call execute_command_line('open ColeCole.pdf')






end program DebyeEXP
