program VanDerPol
!Resuelucion la ecuacion de van der Pol X''(t)+beta*(1-x(t)**2)x'(t)+x(t)=F(t) con runge kutta de orden 2
implicit none

real(kind=8) t0,tmax,dt,x0,y0,beta
real(kind=8), allocatable, dimension (:) :: t,x,y,f
integer i,j,N

!**********************************************************************
t0   = 0.0d0
tmax = 80.0d0
N    = 10200
x0   = 0.19320d0
y0   = -1.8789d0
beta = 3.0d0

allocate(t(0:N),x(0:N),y(0:N),f(0:N))
!**********************************************************************
dt = (tmax - t0) / dble(N) !llenando vector temporal
do i=0,N
  t(i) = t0 + dt * dble(i)
end do
!**********************************************************************
 x(0) = x0 !valores iniciales
 y(0) = y0
 f    = 0
!**********************************************************************
do i=1,N !runge kutta
  do j=1,2
    if (j.eq.1) then
      !if(i.eq.200) then
      !  f(i-1) = 0.01d0 !F(t)= 0.01 cuando t=[(tmax-t0)/N]*200 y cero en todos los demas puntos
      !  y(i-1) = y(i-1) + f(i-1)
      !end if
      f(i)= 5.0d0*sin(t(i)*1.788d0)
      !y(i-1) = y(i-1) - f(i)
      x(i) = x(i-1) + y(i-1) * dt
      y(i) = y(i-1) + (beta*(1-x(i-1)**2)*y(i-1)-x(i-1)-f(i-1)) * dt
    else
      x(i) = 0.5d0 * (x(i-1) + y(i) * dt + x(i))
      y(i) = 0.5d0 * (y(i-1)+(beta*(1-x(i)**2)*y(i)-x(i)-f(i))*dt + y(i)) 
    end if
  end do
end do
!**********************************************************************
open(1,file='van.dat') !llenando archivo
do i=0,N,1
  write(1,*) t(i),x(i),y(i),f(i)
end do
close(1) 
!**********************************************************************
end program VanDerPol