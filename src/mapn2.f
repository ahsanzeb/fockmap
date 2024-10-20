	!---------------------------------------------------
	! map for two modes/sites
	subroutine mapn2(m, map)
	implicit none
	integer, intent(in):: m
	integer, dimension(0:m,m+1), intent(out) :: map	
	integer :: k,i,j
	k = 0;
	do i=0,m ! vib
		do j=i,m ! basis states for single site
			map(i,j+1) = k ! +1 to col index to make it 1:ntot for basis
			map(j,i+1) = k
			k = k + 1;
		end do
	end do
	return
	end subroutine mapn2
