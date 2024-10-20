	!---------------------------------------------------
	! makes blocks of map for more than two sites using some info of prev blocks
	! used the map variable of the main routine.
	subroutine triangles(nrows,lcol,colist,iin,k,m, j1,ntot,map)
	implicit none
	integer, intent(in) ::nrows,lcol, iin, m, j1 , ntot
	integer, dimension(lcol), intent(in) :: colist
	integer, intent(inout) :: k
	integer, dimension(0:m,ntot), intent(inout) :: map
	! local
	integer :: i,j,jj,ii,i0,yi,yf,mjj,mj,j2,km,iii,j11

	!nrows = (iin+1)*iin//2;

	! fill all left side columns for this iin
	j2 = j1+nrows
	i = 0;
	do j=1,lcol
		jj = colist(j)
		map(i,j1:j2-1) = (/(ii,ii=jj,jj+nrows-1)/) ! py: range(0:1) =[0]
		i = i + 1;
	end do

	!write(*,*)'--------------------'
	
	! fill the triangles
	yi = 0;
	do i = 0, iin-1 ! iin \in {1,2,...,M}
		i0 = i+m-iin+1;
		yf = (i+1)*iin - i*(i+1)/2;
		do jj = yi, yf-1
			j = jj - yi;
			do mjj = i0+j, m
				mj = mjj - i0;
				map(mjj,j1+jj) = k
				map(i0+j,j1+yi+mj) = k
				if (j == 0 .and. mjj==m) km = k
				k = k + 1;
			end do ! mjj
		end do ! jj
		!write(*,*)'l:183 maps: i0 = ', i0
		!do iii=1,21
		!write(*,'(6i4)') map(:,iii)
		!end do
		!write(*,*) 'i0, j1+yf, j2-1, km+1, k-1 =',i0,j1+yf,j2-1,km+1,k-1
		!py: Map[yf:,i0] = range(km+1,k)
		map(i0,j1+yf:j2-1) = (/(ii,ii=km+1,k-1)/) ! left column below this triangle
		yi = yf;
	end do !i
	! k has advanced.... goes in the output
	return
	end subroutine triangles
	!---------------------------------------------------

