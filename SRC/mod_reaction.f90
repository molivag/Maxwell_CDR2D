program DebyeEXP

  implicit none
  
  character(len=*), parameter :: PLT_FILE = 'plot.plt' ! Gnuplot file.

  real                                          :: sigma_inf, sigma_0, m, tau, c, gamaP, tauP, aa, bb,cc
  double precision, allocatable, dimension(:)   :: w,  ColeCole_real, ColeCole_imag, DebyeRe, DebyeIM, freq, fvec, x_sol
  double precision, allocatable, dimension(:,:) :: fjac

  double complex, allocatable, dimension(:)   :: ColeCole, Debye               
  double precision                            :: pio_dbuf(2), freqrange_ini(2), freqrange(2), x_sol(2)
  double precision x(2), fvec(36), fjac(36, 2)
  integer                                     :: i, p, nppdecades, ndat_max, ndat, nn

  sigma_inf= 0.125
    sigma_0= 0.05
        tau= 0.125
          m= (sigma_inf - sigma_0)/ sigma_inf
          c= 0.8


  print*, 'Parametros Cole-Cole'
  print*
  print*, 'Conductividad a frec. maxima',sigma_inf
  print*, 'Conductividad a frec. cero..',sigma_0
  print*, 'Cargabilidad................',m
  print*, 'Constante de tiempo.........',tau
  print*, 'Dependencia de la frecuencia',c
  print*

     nppdecades= 5
       ndat_max= 200
  pio_dbuf(1:2)= (/8.0982361D-05,4382.519D0/) ! from em1d settings
  freqrange_ini= pio_dbuf(1:2);
      freqrange= freqrange_ini ! init freq.-range for CCM data creation
  !call ColeColeIP_frqsample(freqrange,nppdecades,ndat_max,ndat)

  aa=log10(freqrange(1))
  bb=log10(freqrange(2))
  cc=bb-aa
  nn=max(1,int(cc)) ! no. of decades
  ndat=nn*nppdecades ! points/decade
  cc=cc/dble(ndat)

  print*,'aa', aa
  print*,'bb', bb
  print*,'cc', cc
  print*,'nn', nn
  print*, 'ndat', ndat
  print*, 'cc', cc
  
 allocate(w(ndat),ColeCole(1:ndat), Debye(1:ndat))
 allocate( ColeCole_real(ndat), ColeCole_imag(ndat), DebyeRe(ndat), DebyeIM(1:ndat)) 
 allocate( x_sol(n), fvec(ndat), fjac(m,n) )

  !Frecuencias de muestreo
  ! Aqui se debe cambiar por ndat y hacer dinamicos todos los arreglos
  ! en lugar de nfreq como dimension, debe ser ndat
  w=0
  do i=1,ndat
    w(i)= 10**aa
    aa=aa+cc
  end do
 !!Muestras
 !ti(1)=aa!-3.0
 !do i=2,ndat
 !  ti(i)= ti(i-1) + 0.1
 !  w(i-1)= 1*10**ti(i-1)
 !end do
 !!Frecuencias de muestreo
 !do i=1,ndat
 !  w(i)= 1*10**ti(i)
 !end do

  !Solucion inicial
  x_sol = (/ 0.5, 0.1/)!;  gamaP = x_sol(1);  tauP = x_sol(2)
  call fcn(ndat, 2, x_sol, fvec, fjac, m, 1)




  !Modelo Cole-Cole
  do i= 1,ndat
    ColeCole(i)= (sigma_inf*( 1 - ( m /( 1+(cmplx(0,1)*w(i)*tau)**c  ) ) ) )
    ColeCole_real(i)= real(ColeCole(i))
    ColeCole_imag(i)=aimag(ColeCole(i))
  end do


  !Modelo Debye
  do i =1,ndat
    Debye(i) = sigma_inf + (sigma_0-sigma_inf) * (gamaP/(1+(cmplx(0,1)*w(i)*tauP) ) )
  end do

  do i=1,ndat
    DebyeRe(i) = real(Debye(i))
    DebyeIM(i) = aimag(Debye(i))
  end do


  100 format (e15.6, 5x,7(e15.6))
  open(unit=1, file="ColeCole.dat")
  write(1,*)'freq CCre CCim DeRE DeIM'
  do i=1,ndat
    write(unit=1, fmt=100) w(i), ColeCole_real(i),ColeCole_imag(i), DebyeRe(i), DebyeIM(i)
  end do
  close(unit=1)

  call execute_command_line('gnuplot ' // PLT_FILE)
  call execute_command_line('open ColeCole.pdf')






end program DebyeEXP


! subroutine fcn(m, n, x, fvec, fjac, ldfjac, iflag)
!   integer m, n, ldfjac, iflag
!   double precision x(n), fvec(m), fjac(ldfjac, n)

!   ! Definir las derivadas parciales respecto a gamma_p y tau_p
!   double precision gamma_p, tau_p
!   double precision df_dgamma_p, df_dtau_p

!   ! Definir valores específicos
!   gamma_p = x(1) ! Asumiendo que gamma_p está en la primera posición de x
!   tau_p = x(2)   ! Asumiendo que tau_p está en la segunda posición de x
!   ! omega es constante

!   if (iflag == 1) then
!     ! Calcular las funciones en fvec
!     fvec(1) = gamma_p / (1.0d0 + i * omega * tau_p)
!     ! Puedes agregar más funciones aquí si es necesario
!   elseif (iflag == 2) then
!     ! Calcular la matriz Jacobiana en fjac
!     df_dgamma_p = 1.0d0 / (1.0d0 + i * omega * tau_p)**2
!     df_dtau_p = -i * gamma_p * omega / (1.0d0 + i * omega * tau_p)**2

!     fjac(1, 1) = df_dgamma_p
!     fjac(1, 2) = df_dtau_p
!     ! Puedes agregar más derivadas parciales aquí si es necesario
!   endif

!   return
! end subroutine fcn



subroutine ColeColeIP_frqsample(f12,ndec,ndat_max,ndat)
integer,intent(in):: ndec,ndat_max
real*8,intent(in):: f12(2)
real*8, allocatable, dimension(:) :: freq
integer,intent(out):: ndat
! Local
real*8:: a,b,c
a=log10(f12(1)); b=log10(f12(2)); c=b-a
! Number of samples depending on number of decades ndec
n=max(1,int(c)) ! no. of decades
ndat=n*ndec ! points/decade
if(ndat>ndat_max/2)ndat=ndat_max/2
if(allocated(freq))deallocate(freq)
allocate(freq(ndat+3))
c=c/dble(ndat)
print *,'Muestras                  Frecuencias de muestreo:',a
do n=1,ndat+1
   freq(n)=10.0**a
   a=a+c
  write(*,'(1x,I3,2x,2(f15.5))') n, a, freq(n)
enddo
ndat=ndat+1
! Number of data = number of frequencies times 2
ndat=ndat*2 ! Re + Im data



end subroutine ColeColeIP_frqsample



! subroutine fcn(m, n, x, fvec, fjac, ldfjac, iflag)
!     integer m, n, ldfjac, iflag
!     double precision x(n), fvec(m), fjac(ldfjac, n)

!     ! Parámetros a estimar (pueden ser globales o argumentos adicionales)
!     double precision, dimension(:), allocatable :: gamma, tau
!     integer :: p

!     ! Asigna valores a los parámetros a estimar (esto puede variar según tu implementación)
!     allocate(gamma(n), tau(n))
!     gamma = x(1:n)
!     tau = x(n+1:2*n)

!     ! Calcula la función y la matriz jacobiana según el valor de iflag
!     if (iflag == 1) then
!         do p = 1, m
!             fvec(p) = gamma(p) / (1.0d0 + cmplx(0.0d0, 1.0d0) * omega * tau(p))
!         end do
!     elseif (iflag == 2) then
!         do p = 1, m
!             fjac(p, 1:n) = 1.0d0 / (1.0d0 + cmplx(0.0d0, 1.0d0) * omega * tau(p))
!             fjac(p, n+1:2*n) = -cmplx(0.0d0, omega * gamma(p)) / ((1.0d0 + cmplx(0.0d0, 1.0d0) * omega * tau(p))**2)
!         end do
!     else
!         print *, 'Error: Valor de iflag no reconocido.'
!         return
!     end if

!     deallocate(gamma, tau)
! end subroutine fcn



subroutine fcn(m, n, x, fvec, fjac, ldfjac, iflag)
    integer m, n, ldfjac, iflag
    double precision x(n), fvec(m), fjac(ldfjac, n)

    ! Parámetros a estimar (pueden ser globales o argumentos adicionales)
    double precision, dimension(:), allocatable :: gamma, tau
    integer :: p

    ! Asigna valores a los parámetros a estimar (esto puede variar según tu implementación)
    allocate(gamma(n), tau(n))
    gamma = x(1:n)
    tau = x(n+1:2*n)

    ! Calcula la función y la matriz jacobiana según el valor de iflag
    if (iflag == 1) then
        do p = 1, m
            fvec(p) = gamma((p+1)/2) / (1.0d0 + cmplx(0.0d0, 1.0d0) * omega * tau((p+1)/2))
        end do
    elseif (iflag == 2) then
        do p = 1, m
            fjac(p, 1:n) = 1.0d0 / (1.0d0 + cmplx(0.0d0, 1.0d0) * omega * tau((p+1)/2))
            fjac(p, n+1:2*n) = -cmplx(0.0d0, omega * gamma((p+1)/2)) / ((1.0d0 + cmplx(0.0d0, 1.0d0) * omega * tau((p+1)/2))**2)
        end do
    else
        print *, 'Error: Valor de iflag no reconocido.'
        return
    end if

    deallocate(gamma, tau)
end subroutine fcn

