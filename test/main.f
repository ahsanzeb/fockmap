
	program test
	implicit none

	integer :: n, mv, ntot
	integer, allocatable, dimension(:) :: pntr, dims
	integer, allocatable, dimension(:,:) :: freq
	double precision, allocatable, dimension(:):: Perm ! dble to avoid overflow; 
	double precision, allocatable, dimension(:,:):: ratio
	integer, allocatable, dimension(:,:) :: map	
	integer(kind=16) :: binomial

	integer:: i

	! uses libfockmap
	external :: sizes, Basis, BasisAll
	external :: Ratios, RatiosAll
	external :: mapping !, writemap
	external :: binomial

	
	write(*,'(a)') "Give N and M:"
	read(*,*) n, mv

	!------------------------------------------------------------------
	! 1. mapping:
	!------------------------------------------------------------------
	ntot = binomial(max(n,1)-1+mv,mv) ! n_min=1 (map_{0 to 1}).

	allocate(map(0:mv, ntot))
	call mapping(n,mv, ntot, map)
	!------------------------------------------------------------------
	do i=1,ntot
		write(*,'(1000i5)') map(0:mv,i)-1
	enddo



	!------------------------------------------------------------------
	! 2. basis states related variables
	!------------------------------------------------------------------
	! 	  2.1. frequencies and permutations
	!.....................................

	allocate(pntr(0:n))
	allocate(dims(0:n))
	! calculates the sizes of permut symm spaces
	! for s=0-n modes, each having 0:mv levels
	call sizes(n,mv, pntr, dims, ntot)

	write(*,'(a,1000i5)') "ntot: ", ntot
	write(*,'(a,1000i5)') "dim : ", dims
	write(*,'(a,1000i5)') "pntr: ", pntr

	allocate(perm(0:ntot-1))
	allocate(freq(0:mv,0:ntot-1))

	!call BasisAll(n,mv, pntr, dims, ntot, freq, perm) ! 0:n mode cases
	call Basis(n,mv, pntr, dims, ntot, freq, perm) ! only n-1,n mode case

	!write(*,'(a,10000f10.5)') "freq: all modes "
	!do i= 0, ntot-1
	!	write(*,'(1000i5)') freq(0:mv,i)
	!enddo

	! Basis has a bug. fixed now. but published version with comput paper has this issue. 
	if (1==0) then
	
	write(*,'(a,10000f10.5)') "freq: N-1 modes "
	do i= pntr(n-1)+1,pntr(n-1)+dims(n-1)
		write(*,'(1000i5)') freq(0:mv,i)
	enddo
	write(*,'(a,10000f10.5)') '***************************'
	write(*,'(a,10000f10.5)') "freq: N modes"
	do i= pntr(n)+1,pntr(n)+dims(n)
		write(*,'(1000i5)') freq(0:mv,i)
	enddo
	
	endif

	if (1==0) then
	
	call BasisAll(n,mv, pntr, dims, ntot, freq, perm) ! 0:n mode cases
	
	write(*,'(a,10000f10.5)') "ALL: freq: N-1 modes "
	do i= pntr(n-1)+1,pntr(n-1)+dims(n-1)
		write(*,'(1000i5)') freq(0:mv,i)
	enddo
	write(*,'(a,10000f10.5)') '***************************'
	write(*,'(a,10000f10.5)') "ALL: freq: N modes"
	do i= pntr(n)+1,pntr(n)+dims(n)
		write(*,'(1000i5)') freq(0:mv,i)
	enddo
	
	endif


	!.....................................
	!   2.2 ratios of permutations of states linked by the mapping
	!       ratio(k,i) = sqrt(P_{p}(i)/P_{p+1}(j)) where j=map(k,i)
	!.....................................

	!allocate(ratio(0:mv,0:pntr(n))) ! pntr(n) = ntot_{n-1}
	!call RatiosAll(n,mv, ntot, pntr, freq, ratio) ! 0:n-1 mode cases
	!write(*,'(a,10000f10.5)') "ratio: ", ratio

	!call Ratios(n,mv, ntot, freq, ratio) ! only n-1 to n case
	!------------------------------------------------------------------


	stop
	end program test
