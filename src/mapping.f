

C ********************************************************************
C Generates the mapping
C *********************** INPUT **************************************
C integer n              : number of subsystems
C integer m              : determined the number of levels, m+1.
C integer ntot	           : number of permutation symmetric states
C                          for n-1 subsystems.
C *********************** OUTPUT **************************************
C integer map(0:m,ntot)  : the mapping from n-1 to n subsystem.
C *********************************************************************

	subroutine mapping(n,m, ntot, map)
	implicit none
	integer, intent(in) :: n,m, ntot
	integer, dimension(0:m,ntot), intent(out) :: map	
	! local
	integer, dimension(m+1) :: colist
	integer, dimension(:), allocatable :: list1, list2
	integer :: lmax, iarg,ii,j1,l1,l2, karg,nn,nrows,iii,i,x
	integer :: lcol
	integer(kind=16)  :: binomial
	external :: mapn2, triangles, binomial

	map(:,:) = 0

	if(n==1) then
	 map(:,1) = (/ (ii, ii=0,m) /);
	 map = map + 1
	 return
	endif

	! map(0:m,1:m+1) set by mapn2(m)
	call mapn2(m, map(0:m,1:m+1)) ! sets map for two sites, the starting square block of m+1 x m+1.

	if (n < 3 ) then
		!ntotx = (m+1) ! number of states of n-1 sites in the map
		map = map + 1
		return
	endif

	! max size list1,2 can get:
	lmax = ntot !binomial(n+m,m) ! roughly speaking... 
	allocate(list1(lmax))
	allocate(list2(lmax))
	
	!karg = mapfull[-1,-1] +1;
	karg = map(m,m+1) + 1; ! argument for next block's triangles.
	list1(1) = m !(/ ( m ) /);
	l1 = 1; ! length of list1
	j1 = m+1; ! last baiss index in the map for two sites; 1 :=> 2 sites

	do nn=3,n ! no of sites
		!write(*,*) 'list1 = ',list1(1:l1)
		do ii=1,l1
			iarg = list1(ii)
			!py: colist = mapfull[-1,0:m-iarg+1]+1; py: x=[0,1,2]; x[0:1] =[0]
			colist(1:m-iarg+1) = map(0:m-iarg,j1) + 1 ! j1 last row so far
			nrows = (iarg+1)*iarg/2;
			lcol = m-iarg+1
			call triangles(nrows,lcol,colist(1:lcol),iarg,karg,m,j1+1,
     .                                      ntot, map) ! karg inout
			j1 = j1 + nrows
		end do 
		!call mkarglist(l1,list1(1:l1),l2,list2(1:l2))
		!..............................
		if(nn<n) then
			l2 = 0;
			do i=1,l1
				x = list1(i);
				list2(l2+1:l2+x) = (/ (ii,ii=x,1,-1) /)
				l2 = l2 + x;
			end do
			!..............................	
			!write(*,*) 'l1, l2= '	,l1,l2,'lst=',list2(1:l2)
			list1(1:l2) = list2(1:l2)!
			l1 = l2;
		endif
	end do

	map = map + 1 ! basis indices 1:*
	!ntotx = j1 ! col dim of map
	deallocate(list1,list2)
	
	return
	end subroutine mapping
