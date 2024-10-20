
!.................................................
! Arguments:
!input:
! 	n = number of subsystems
! mv = determined the number of levels, m+1.
!output:
! pntr = shifts in the index of the basis states for p-subsystems, p=0,1,...,n
! dims = sizes of the permutation symmetric spaces for p-subsystems, p=0,1,...,n
! ntot = sum(dims), combined size of perm sym subspaces for p-subsystems, p=0,1,...,n
!.................................................

	!===================================================================
	! calculate the sizes of permutation symmtric spaces for 0-n modes
	!	pntr(s)+1 is the start index for s-mode case
	subroutine sizes(n,mv, pntr, dims, ntot)
	implicit none
	integer, intent(in) :: n, mv
	integer, dimension(0:n), intent(out) :: pntr, dims
	integer, intent(out) :: ntot
	integer :: s
	integer(kind=16) :: binomial
	
	dims(0) = 1;
	pntr(0) = 0;

	do s=1,n
		dims(s) = binomial(s+mv,mv)
		pntr(s) = pntr(s-1) + dims(s-1);
	enddo
	ntot = pntr(n) + dims(n)

	pntr = pntr - 1;! shift by -1; use pntr as a shift
	return
	end subroutine sizes
	!===================================================================

