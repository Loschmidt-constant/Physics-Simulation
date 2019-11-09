!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module constant_values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

integer,parameter		:: n_max = 1000
integer,parameter		:: m_max = 10

real(8)	,parameter		:: pi = 4d0*atan(1d0)

real(4),parameter			:: x_min = 0.0
real(4),parameter			:: x_max = 1.0
real(8),parameter 		:: dx   = (x_max -x_min)/dble(m_max)	!! 空間刻み幅
real(8),parameter 		:: dt   = 1.0/dble(n_max)						!! 時間刻み幅

real(4),parameter		:: kappa = 0.5	!! 拡散係数

end module constant_values


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module set_functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

end module set_functions


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use constant_values
use set_functions

implicit none

integer									:: i
integer									:: j
real,dimension(0:m_max)		:: u
real,dimension(0:m_max)		:: x
real											:: t
character(len=50) :: filename
character(len=50) :: tmp

do i = 0, m_max
 x(i) = dx * i
end do

u(0) = 0.0
u(m_max) = 10.0

call get_ICs(u)

do j = 0, n_max
	t = dt * j
	call calculation_unew(u,x,t)	
	if (mod(j,100)==0) then
		write(tmp,'(i4)') j 
		filename='diffusion-ana_'//trim(adjustl(tmp))//'.dat'
		open(10, file = filename, status="replace")
		
		do i = 0, m_max
			write(10,*) x(i), u(i)
		end do
		
		close(10)
		
	end if
end do

stop
end program main


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ICs(u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use constant_values
use set_functions

implicit none
integer								:: i
real										:: x
real,dimension(0:m_max)	:: u


do i = 1, m_max-1
	x = dx * i
	
	if (x >= 0.0 .and. x < 0.5) then
		u(i) = 0.0
	else if(x>=0.5 .and. x<=1) then
		u(i) = 10
	end if

end do


end subroutine get_ICs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculation_unew(u,x,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use constant_values
use set_functions

implicit none
integer			  					:: i
real,dimension(0:m_max)	:: u
real,dimension(0:m_max)	:: u_old
real,dimension(0:m_max)	:: x
real,dimension(0:m_max)	:: c
real										:: t

do i = 0, m_max
	u_old(i) = u(i)
end do

do i = 1, m_max
			c(i) = -(20/(pi*i))*(cos(pi*i*x(i))-cos(0.5*pi*i*x(i)))
			u(i) = u_old(i-1) + c(i)*exp(-0.5*t*(pi*i)**2)*sin(pi*i*x(i))
end do

!write(*,*) u(1)

end subroutine calculation_unew
