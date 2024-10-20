


C ********************************************************************
C calculates the ratio: sqrt(P_{p}(i)/P_{p+1}(j)), j=map(k,i)
C for all 0,1,2....,n sites/modes
C *********************** INPUT **************************************
C integer n              : number of subsystems
C integer mv             : determined the number of levels, m+1.
C integer pntr(n+1)      : shifts in the index of the basis states 
C                        : for p-subsystems, p=0,1,...,n
C integer ntot	           : number of permutation symmetric states
C integer freq(mv+1,ntot) : occupation frequencies for the basis states
C *********************** OUTPUT **************************************
C double ratio(mv+1,ntot) : squareroot of the ratio of the number of
C                            permutation of states with p subsystems 
C                            and states with p+1 subsystems that are 
C                            linked to them with the mapping.
C *********************************************************************

	!===================================================================
	! r = sqrt(P_{p}(i)/P_{p+1}(j)), j=map(k,i)
	subroutine RatiosAll(n,mv,ntot, pntr, freq, ratio)
	implicit none
	integer, intent(in) :: n, mv, ntot
	integer, dimension(0:n), intent(in) :: pntr
	integer, dimension(0:mv,0:ntot-1), intent(in) :: freq
	double precision, dimension(0:mv,0:ntot-1), intent(out) :: ratio
	integer :: i, k,p
	integer :: gk1
	do p=0,n-1
	 do i=pntr(p)+1,pntr(p+1)
	  do k=0,mv
	   gk1 = freq(k,i) + 1;
	   ratio(k,i) = dsqrt(gk1*1.0d0/(p+1));
	  end do ! k levels
	 end do ! i states
	end do! p sectors/modes/sites
	return
	end subroutine RatiosAll
