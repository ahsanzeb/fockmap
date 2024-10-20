	integer (kind=16) function binomial(n,m) 
	! length of block that would be added
	! to the map if we add a site to p-1 sites
	implicit none
	integer, intent(in):: n,m
	integer :: i, smaller, greater
	integer(kind=16)  :: factorial
	external :: factorial

	if(n<m) then
		binomial = 0;
	elseif(n==m) then
		binomial = 1;		
	else
		!binomial = factorial(n)/(factorial(m)*factorial(n-m))	
		! remove the largest factorial
		if(m .ge. n-m) then
			greater = m
			smaller = n-m
		else
			greater = n-m
			smaller = m
		endif

		binomial = 1
		do i=greater+1,n
				binomial = binomial*i
		end do
		binomial = binomial/factorial(smaller)
			
	endif
	return
	end function binomial

