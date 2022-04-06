program Anticipation
!Resuelucion la ecuacion de un oscilador forzado generalizado X'''(t)+alpha*X''(t)+beta*(1-X(t)**2)*X'(t)+X(t)=0 con runge kutta de orden 2
implicit none

real(kind=8) t0,tmax,dt,x0,y0,z0,alpha,beta
real(kind=8), allocatable, dimension (:) :: t,x,y,z
integer i,j,N

!**********************************************************************
t0    = 0.0d0
tmax  = 500.0d0
N     = 20500
x0    = 0.01d0
y0    = 0.0d0
z0    = 0.0d0
alpha = 0.5d0
beta  = -5.6d0

allocate(t(0:N),x(0:N),y(0:N),z(0:N))
!**********************************************************************
dt = (tmax - t0) / dble(N) !llenando vector temporal
do i=0,N
  t(i) = t0 + dt * dble(i)
end do
!**********************************************************************
 x(0) = x0 !valores iniciales
 y(0) = y0
 z(0) = z0
!**********************************************************************
do i=1,N !runge kutta
  do j=1,2
    if (j.eq.1) then
      x(i) = x(i-1) + y(i-1) * dt
      y(i) = y(i-1) + z(i-1) * dt
      z(i) = z(i-1) + (-alpha*z(i-1)-beta*(1-x(i-1)**2)*y(i-1)-x(i-1)) * dt
    else
      x(i) = 0.5d0 * (x(i-1) + y(i)*dt + x(i))
      y(i) = 0.5d0 * (y(i-1) + z(i)*dt + y(i))
      z(i) = 0.5d0 * (z(i-1) + (-alpha*z(i)-beta*(1-x(i)**2)*y(i)-x(i))*dt + z(i))
    end if
  end do
end do
!**********************************************************************
open(1,file='ant.dat') !llenando archivo
do i=0,N,1
  write(1,*) t(i),x(i),y(i),z(i)
end do
close(1) 
!**********************************************************************
end program Anticipation