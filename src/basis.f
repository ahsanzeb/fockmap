
C ********************************************************************
C calculates the occupation frequencies and number of permutations
C for n-1 and n sites/modes.
C *********************** INPUT **************************************
C integer n              : number of subsystems
C integer mv             : determined the number of levels, m+1.
C integer pntr(n+1)      : shifts in the index of the basis states 
C                        : for p-subsystems, p=n-1,n
C integer dims(n+1)      : sizes of the permutation symmetric spaces
C integer ntot	           : number of permutation symmetric states
C                          ntot= dims(n-1)+dims(n)                         
C *********************** OUTPUT **************************************
C integer freq(mv+1,ntot) : occupation frequencies for the basis states
C double perm(ntot)      : number of permutations for the basis states
C *********************************************************************

	!===================================================================
	!..........................................
	! permutation symmetric vibrational basis: retruns data for n-1 and n site bases
	!..........................................
	subroutine Basis(n,mv, pntr, dims, ntot, freq, perm)
	!	calculates Perm, frequiencies of occup etc 
	implicit none

	integer, intent(in) :: n,mv, ntot
	integer, dimension(0:n), intent(in) :: pntr, dims
	integer, dimension(0:mv,0:ntot-1), intent(out)  :: freq
	! actually integers, but define dble to avoid overflow; 
	double precision, dimension(0:ntot-1), intent(out) :: Perm
	
	! local
	integer, dimension(:,:), allocatable :: seti, setf
	integer :: i,j,k,rk,s, p
	double precision :: fac
	! aux arrays
	double precision, allocatable, dimension(:) :: Px1,Px2
	integer, allocatable, dimension(:,:) :: fx1,fx2

	
	! with n sites, mv+1 states per mode
	! aux arrays to hold actual sets of occup for all
	! perm sym basis for a given number of sites.
	allocate(seti(n,dims(n)))
	allocate(setf(n,dims(n)))

	!..........................................
	! manulally set values for 0,1 sites
	Perm(0:0) = (/ 1.0 /)
	freq(:,0) = (/(0,i=0,mv)/)
	Perm(1:mv+1) = 1.0; !(/(1.0,i=0,mv)/)
	freq(:,1:mv+1) = 0
	do i=1,mv+1!pntr(1)+1,pntr(1)+dims(1)
		freq(i-1,i) = 1; ! every occup in every site has freq 1
	enddo

	!========================================================
	!proceed only if n > 1 
	if (n > 1) then 
	!========================================================
	! n > 1: aux arrays
	allocate(fx1(0:mv,dims(n)))
	allocate(fx2(0:mv,dims(n)))
	allocate(Px1(dims(n)))
	allocate(Px2(dims(n)))
	fx1 = 0; fx2 = 0;
	Px1 = 0.0d0; Px2 = 0.0d0;

	!..........................................
	! one site: n=1; mv+1 states= 0,1,2,...,mv
	! prepend sites one by one, so fist site is at nth position.
	seti(n,1:mv+1) = (/(i,i=0,mv)/) ! set nth site occupations, add sites in reverse direction.
	fx1(:,1:mv+1) = freq(:,1:mv+1);
	Px1(1:mv+1) = Perm(1:mv+1);
	!..........................................
	! use data for 1 site to generate data for 2, and so on....	
	!..........................................
	do s=1,n-1 ! take a given s sector, make s+1 sector.
		j = 0;
		do i= 1, dims(s)         !pntr(s),pntr(s)+dims(s) !1,ntoti
			!write(*,*) 'seti(n-s+1,i), j ,i= ',seti(n-s+1,i), j,i
			do k= seti(n-s+1,i), mv ! prepend only >= occupations 
				j = j + 1 ! advance index of final sets
				setf(n-s+1:n,j) = seti(n-s+1:n,i) ! copy all occup in s case 
				setf(n-s,j) = k ! set occup of additional site, prepended. 

				!Nv(j2) = Nv(i1) + k ! number of vib quanta
				fx2(:,j) = fx1(:,i) ! freq of parent set
				fx2(k,j) = fx2(k,j) + 1 ! freq of k increased by 1
				! { P(j) = P(i) * (s+1)*rk!/(rk+1)! }  ! no need to repeat comput.
				!	use table of factorials: fact(i) = i!; i in [0,n] ; 
				!													table computed by a function once for all.
				rk = fx1(k,i);
				fac = (s+1)*1.0d0/(rk+1) !(s+1)*fact(rk)/fact(rk+1);
				Px2(j) = Px1(i)*fac
			end do! k
		end do ! i
		seti(n-s:n,1:j) = setf(n-s:n,1:j)
		
		if(s < n-1) then ! dont do for the last iteration: n-1 to n
			fx1(:,1:dims(s+1)) = fx2(:,1:dims(s+1))
			Px1(1:dims(s+1)) = Px2(1:dims(s+1))
		endif
		
	enddo ! s
	!..........................................
	deallocate(seti,setf)

	! set output variables:
	 Perm(pntr(n-1)+1:pntr(n-1)+dims(n-1) ) = Px1(1:dims(n-1))
	 Perm(pntr(n)+1:pntr(n)+dims(n) ) = Px2(1:dims(n))

	 freq(:,pntr(n-1)+1:pntr(n-1)+dims(n-1) ) = fx1(:,1:dims(n-1))
	 freq(:,pntr(n)+1:pntr(n)+dims(n) ) = fx2(:,1:dims(n))

	endif ! n > 1
	!========================================================

	!write(*,'(a,10000f10.5)') "freq: N-1 modes "
	!do i= pntr(n-1)+1,pntr(n-1)+dims(n-1)
	!	write(*,'(1000i5)') freq(0:mv,i)
	!enddo

!		open(1,file='basis-new.dat', form="formatted", action="write")
!		do i=0,ntot-1
!		 write(1,'(1000000i10)') freq(:,i)
!		enddo
!		write(1,'(1000000f10.2)') Perm
!		close(1)
		
	return
	end 	subroutine Basis
