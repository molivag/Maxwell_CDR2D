program DebyeEXP

  implicit none
  
  integer, parameter::nfrec=41, nfrec_rel=12
  real              ::sigma_inf, sigma_0, m, tau, c
  integer           :: ND !No. de fun. debye debe ser 1 a 1.5 veces el numero de decadas de w
  integer           :: Nf !No. de frec. a las que se muestreara sigmaIM(w)
  !REAL :: gammaIM=0.0001 !Factor de regularización imaginarios
  !REAL :: gammaRE=0.0002 !Factor de regularización reales
  !REAL ::frecmin, frecmax

  double precision                             :: ti(1:nfrec), w(1:nfrec)
  double precision                             :: qr(1:nfrec_rel), taur(1:nfrec_rel), wr(1:nfrec_rel)
  double precision                             :: ColeCole_real(1:nfrec), ColeCole_imag(1:nfrec)
  double precision, allocatable, dimension(:)  :: ci, cr, c_tot, A_imag, A_real, A, DiTci
  double precision, allocatable, dimension(:,:):: Di, Dr, D, DiT, DtD_imag
  double complex                                    :: ColeCole(1:nfrec)
  integer :: i, p!, k, L
  !REAL, allocatable,dimension(:) :: frec
  !REAL,allocatable,dimension(:) :: w

  sigma_inf= 500.0
    sigma_0= 150.0
          m= (sigma_inf - sigma_0)/ sigma_inf
        tau= 1.0
          c= 0.75
         ND= 2
         Nf= 2*ND

  allocate(ci(Nf), cr(Nf), c_tot(Nf), A_imag(ND), A_real(ND), A(ND), DiTci(ND))

  allocate(Di(Nf,ND), Dr(Nf,ND), D(Nf,ND), DiT(ND,Nf),DtD_imag(ND,ND))


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
  end do

  !Frecuencias de muestreo
  do i=1,nfrec
   w(i)= 1*10**ti(i)
  end do

  print *,'Muestras                  Frecuencias de muestreo:'
  do i=1,nfrec
   write(*,'(1x,I3,2x,2(f15.5))') i, ti(i), w(i)
  end do

  ! * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * !
  print*, '~ ~ ~ Parametros para la expansión de Debye ~ ~ ~'
  
  !el espaciado en el muestreoo
  qr(1)=1.1
  do i=2,nfrec_rel
    qr(i)= qr(i-1) - 0.4
  end do

  !el tiempo de relajacion de muestreo (propuesto)
  do i=1,nfrec_rel
   taur(i)= 1*10**qr(i)
  end do
  !taur(1:nfrec_rel)=1*10**qr(1:nfrec_rel)

  !la frecuencia de muestreo (calculada)
  do i=1,nfrec_rel
    wr(i)= 1/taur(i)
  end do

  print*
  print *,'Frecuencias de relajación:'
  do i=1,nfrec_rel
    print *,wr(i)
  end do


  !Modelo Cole-Cole
  do i= 1,nfrec
   ColeCole(i)= sigma_inf*( 1 - ( m /( 1+(cmplx(0,1)*w(i)*tau)**c  ) ) )
  end do
  !tambien puedo realizarlo de esta forma:
  !ColeCole(1:nfrec) = sigma_inf*( 1 - ( m /( 1+(complex(0,1)*w(1:nfrec)*tau)**c  ) ) )
  !lo anterior es posible por haber definido parameter ademas de la forma en como defini
  !las dimensiones del arreglo.
  !print*
  !print*,'Modelo Cole-Cole'
  !print*,ColeCole
  !print*,realpart(ColeCole)
  !print*
  !print*,imagpart(ColeCole)
  
  do i=1,nfrec
    ColeCole_real(i)= real(ColeCole(i))
    ColeCole_imag(i)=aimag(ColeCole(i))
  end do

  !print*,'* * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
  !PRINT*,ColeCole_real
  !print*,'* * * * * * * * * * * * * * * * * * * * * * * * * * * * *'

  !pesos parte Imaginaria

  do i=1,Nf
    ci(i)=ColeCole_imag(i)/(sigma_0 - sigma_inf)
  end do

  do i=1,Nf
    do p=1,ND
      Di(i,p) = (w(i)/wr(p)) / (1+ (w(i)/wr(p))**2)
    end do
  end do

  print*,'* * * * * * * * * matriz Di * * ** * * * * * * *'
  do i= 1,Nf
    print*,(Di(i,p),p = 1,ND)
  enddo

  !DiT = transpose(Di)
  do p=1,ND
    do i=1,Nf
      DiT(p,i) = Di(i,p)
    end do
  end do


  print*,'* * * * * * * * * matriz Di transpose * * ** * * * * * * *'

  do p= 1,ND
    print*,(DiT(p,i),i = 1,Nf)
  end do




  !Calculo del producto DtD
  !do i=1,Nf
  !  do p =1,ND
  !    DtD_imag(i,p)=0.0
  !    do k=1,4
  !      DtD_imag(i,p) = DtD_imag(i,p) + DiT(p,k)*Di(k,p)
  !    end do
  !  end do
  !end do

  DtD_imag = matmul(DiT,Di)

  print*,'* * * * * * * * * matriz DtD* * ** * * * * * * *'
  do p= 1,ND
    print*,(DtD_imag(p,i),i=1,ND)
  end do
  !aqui falta calcular la inversa de DtD_imag =A_imag
  DiTci = matmul(DiT,ci)
  print*,'* * * * * * * * * matriz DtD* * ** * * * * * * *'
  do p= 1,ND
    print*,DiTci(p)
  end do
  !Con las multiplicaciones hechas, solo necesito la inversa de
  ! DtD_imag para obtener los pesos
  !A_imag = inv(DtD_imag)*DiTci

  100 format (f15.6, 5x, f15.6, 5x, f15.6)
  open(unit=1, file="ColeCole.dat")
  do i=1,nfrec
    write(unit=1, fmt=100) w(i), ColeCole_real(i),ColeCole_imag(i)
  end do
  close(unit=1)


end program DebyeEXP
