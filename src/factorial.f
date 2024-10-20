	integer (kind=16) function factorial(p) 
	! length of block that would be added
	! to the map if we add a site to p-1 sites
	implicit none
	integer, intent(in):: p
	integer :: i
	factorial = 1;
	do i=2,p
		factorial = factorial*i
	end do
	return
	end function factorial
