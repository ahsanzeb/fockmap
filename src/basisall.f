


C ********************************************************************
C calculates the occupation frequencies and number of permutations
C for all 0,1,....n sites/modes.
C *********************** INPUT **************************************
C integer n              : number of subsystems
C integer mv             : determined the number of levels, m+1.
C integer pntr(n+1)      : shifts in the index of the basis states 
C                        : for p-subsystems, p=0,1,...,n
C integer dims(n+1)      : sizes of the permutation symmetric spaces
C integer ntot	           : number of permutation symmetric states=sum(dims)                         
C *********************** OUTPUT **************************************
C integer freq(mv+1,ntot) : occupation frequencies for the basis states
C double perm(ntot)      : number of permutations for the basis states
C *********************************************************************

	!===================================================================
	!..........................................
	! permutation symmetric vibrational basis
	!..........................................
	subroutine BasisAll(n,mv, pntr, dims, ntot, freq, perm)
	!	calculates Perm, frequiencies of occup etc 
	implicit none

	integer, intent(in) :: n,mv, ntot
	integer, dimension(0:n), intent(in) :: pntr, dims
	integer, dimension(0:mv,0:ntot-1), intent(out)  :: freq
	! actually integers, but define dble to avoid overflow; 
	double precision, dimension(0:ntot-1), intent(out) :: Perm
	
	! local
	integer, dimension(:,:), allocatable :: seti, setf
	integer :: i,j,k,rk,s, p, i1,j2
	double precision :: fac
	
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
	do i=pntr(1)+1,pntr(1)+dims(1)
		freq(i-1,i) = 1; ! every occup in every site has freq 1
	enddo
	!..........................................
	! one site: n=1; mv+1 states= 0,1,2,...,mv
	! prepend sites one by one, so fist site is at nth position.
	seti(n,1:mv+1) = (/(i,i=0,mv)/) ! set nth site occupations, add sites in reverse direction.
	!..........................................
	! use data for 1 site to generate data for 2, and so on....	
	!..........................................
	do s=1,n-1 ! take a given s sector, make s+1 sector.
		j = 0;
		do i= 1, dims(s)         !pntr(s),pntr(s)+dims(s) !1,ntoti
			!write(*,*) 'seti(n-s+1,i), j ,i= ',seti(n-s+1,i), j,i
			do k= seti(n-s+1,i), mv ! prepend only >= occupations 
				j = j + 1 ! advance index of final sets
				!write(*,*) 'k, j = ', k, j
				setf(n-s+1:n,j) = seti(n-s+1:n,i) ! copy all occup in s case 
				setf(n-s,j) = k ! set occup of additional site, prepended. 

	      ! shift index
	      i1= pntr(s) + i ! s sites
	      j2= pntr(s+1) + j ! s+1 sites

				freq(:,j2) = freq(:,i1) ! freq of parent set
				freq(k,j2) = freq(k,j2) + 1 ! freq of k increased by 1
				! { P(j) = P(i) * (s+1)*rk!/(rk+1)! }  ! no need to repeat comput.
				!	use table of factorials: fact(i) = i!; i in [0,n] ; 
				!													table computed by a function once for all.
				rk = freq(k,i1);
				fac = (s+1)*1.0d0/(rk+1) !(s+1)*fact(rk)/fact(rk+1);
				Perm(j2) = Perm(i1)*fac
			end do! k
		end do ! i
		seti(n-s:n,1:j) = setf(n-s:n,1:j)
	enddo ! s
	!..........................................
	deallocate(seti,setf)



!		open(1,file='basis-new.dat', form="formatted", action="write")
!		do i=0,ntot-1
!		 write(1,'(1000000i10)') freq(:,i)
!		enddo
!		write(1,'(1000000f10.2)') Perm	
!		close(1)
		
	return
	end 	subroutine BasisAll
