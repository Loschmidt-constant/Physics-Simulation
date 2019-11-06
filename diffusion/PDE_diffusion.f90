!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module constant_values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

integer,parameter		:: n_max = 1000
integer,parameter		:: m_max = 10

real(4),parameter		:: x_min = 0.0
real(4),parameter		:: x_max = 1.0
real(8),parameter 		:: dx   = 0.1	!! 空間刻み幅
real(8),parameter 		:: dt   = 0.00001						!! 時間刻み幅

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

integer										:: i
integer										:: j
real,dimension(0:m_max)		:: u
real,dimension(0:m_max)		:: x
real,dimension(0:n_max)			:: t
character(len=50) :: filename
character(len=50) :: tmp
real lambda

lambda = kappa * (dt/(dx**2))	!! 刻み幅定数
 
u(0) = 0.0
u(1) = 0.0

do i = 0, m_max
 x(i) = dx * i
end do

	call get_ICs(u)

	do j = 1, n_max
		t = dt * j
		call calculation_unew(u)
		if (mod(j,100)==0) then
			write(tmp,'(i4)') j 
			filename='diffusion_'//trim(adjustl(tmp))//'.dat'
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
integer									:: i
real										:: x
real,dimension(0:m_max)	:: u


do i = 1, m_max-1
	x = dx * i
	
	if (x >= 0.25 .and. x < 0.75) then
		u(i) = 1.0
	else 
		u(i) = 0.0
	end if

end do


end subroutine get_ICs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculation_unew(u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use constant_values
use set_functions

implicit none
integer			  						:: i
real,dimension(0:m_max)	:: u
real,dimension(0:m_max)	:: u_old
real lambda

lambda = kappa * (dt/(dx**2))	!! 刻み幅定数

do i = 0, m_max
	u_old(i) = u(i)
end do

do i = 1, m_max-1
	u(i) = u_old(i) + lambda * (u_old(i+1) - 2*u_old(i) + u_old(i-1))
end do


end subroutine calculation_unew
