

C ********************************************************************
C calculates the ratio: sqrt(P_{p}(i)/P_{p+1}(j)), j=map(k,i)
C for a fixed number of sites/modes
C *********************** INPUT **************************************
C integer n              : number of subsystems
C integer mv             : determined the number of levels, m+1.
C integer pntr(n+1)      : shifts in the index of the basis states 
C                        : for p-subsystems, p fixed.
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
	subroutine PermRatios(p,mv,ntot, freq, ratio) ! for p sites
	implicit none
	integer, intent(in) :: p, mv, ntot
	integer, dimension(0:mv,ntot), intent(in) :: freq
	double precision, dimension(0:mv,ntot), intent(out) :: ratio
	integer :: i, k
	integer :: gk1

	if(p /= sum(freq(:,1))) then
	 write(*,'(a)') "Error(PermRatios2): wrong number of modes!"
	 stop
	endif

	do i= 1, ntot
	 !p = sum(freq(:,i)) ! dont need for every state if we know the sites/sec index from pntr,dims etc.
	 do k=0,mv
	  gk1 = freq(k,i) + 1;
	  ratio(k,i) = dsqrt(gk1*1.0d0/(p+1));
	 end do
	end do
	return
	end subroutine PermRatios
