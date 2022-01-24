
! From:
!
!     http://icl.cs.utk.edu/lapack-for-windows/lapack/
!

program main

	implicit none

	double precision, parameter :: pi = 4.d0 * atan(1.d0)

	double complex, allocatable :: a(:,:), x(:), b(:), s(:)

	integer :: i, n, ns, nnonzero, nmat, lda, nrhs, ldb, info

	integer, allocatable :: inz(:,:), ipiv(:)

	!*******************************************************************

	! Dan Piponi's example:
	!
	!     https://twitter.com/sigfpe/status/1484960640284258309?s=20
	!

	! Size of matrices
	n = 8

	! Size of generator sample set
	ns = 3

	allocate(a(n,n))
	allocate(x(n))
	allocate(b(n))

	! Generator sample set
	allocate(s(ns))

	! Cube roots of 1
	s(1) = complex(1.d0, 0.d0)
	s(2) = complex(cos(2.d0 * pi / 3.d0), sin(2.d0 * pi / 3.d0))
	s(3) = complex(cos(4.d0 * pi / 3.d0), sin(4.d0 * pi / 3.d0))

	!print *, 's(2) = ', s(2)
	print *, 's = ', s

	! Number of non-zeros

	! Upper triangle, main diagonal, and first band below diagonal
	nnonzero = n * (n + 1) / 2 + n - 1

	print *, 'nnonzero = ', nnonzero

	! Indices of non-zeros
	allocate(inz(2, nnonzero))

	! TODO: set inz in a non-insane way

	inz(:, 1) = [1, 1]
	inz(:, 2) = [1, 2]
	inz(:, 3) = [1, 3]
	inz(:, 4) = [1, 4]
	inz(:, 5) = [1, 5]
	inz(:, 6) = [1, 6]
	inz(:, 7) = [1, 7]
	inz(:, 8) = [1, 8]

	inz(:, 9) = [2, 1]
	inz(:,10) = [2, 2]
	inz(:,11) = [2, 3]
	inz(:,12) = [2, 4]
	inz(:,13) = [2, 5]
	inz(:,14) = [2, 6]
	inz(:,15) = [2, 7]
	inz(:,16) = [2, 8]

	inz(:,17) = [3, 2]
	inz(:,18) = [3, 3]
	inz(:,19) = [3, 4]
	inz(:,20) = [3, 5]
	inz(:,21) = [3, 6]
	inz(:,22) = [3, 7]
	inz(:,23) = [3, 8]

	inz(:,24) = [4, 3]
	inz(:,25) = [4, 4]
	inz(:,26) = [4, 5]
	inz(:,27) = [4, 6]
	inz(:,28) = [4, 7]
	inz(:,29) = [4, 8]

	inz(:,30) = [5, 4]
	inz(:,31) = [5, 5]
	inz(:,32) = [5, 6]
	inz(:,33) = [5, 7]
	inz(:,34) = [5, 8]

	inz(:,35) = [6, 5]
	inz(:,36) = [6, 6]
	inz(:,37) = [6, 7]
	inz(:,38) = [6, 8]

	inz(:,39) = [7, 6]
	inz(:,40) = [7, 7]
	inz(:,41) = [7, 8]

	inz(:,42) = [8, 7]
	inz(:,43) = [8, 8]

	a = complex(0.d0, 0.d0)
	do i = 1, nnonzero
		a(inz(1,i), inz(2,i)) = s(mod(i,ns)+1)
	end do

	! There is no Fortran edit descriptor for complex numbers :(
	print *, 'a = '
	print '(16es14.4)', a
	!print *, 'a = ', a

	! For top-left element only of inverse, we don't need to invert.  Just solve
	! a linear system with b = [1, 0, 0, 0, ...]

	b = 0.d0
	b(1) = complex(1.d0, 0.d0)

	! Solve the system

	lda = n
	nrhs = 1
	ldb = n

	allocate(ipiv(n))

	x = b
	!print *, 'x = ', x
	call zgesv(n, nrhs, a, lda, ipiv, x, ldb, info)

	print *, 'info = ', info
	!print *, 'ipiv = ', ipiv
	print *, 'x = ', x
	!print *, 'a = ', a

	! TODO: iterate through all permutations of a.  Save results for plotting

end program main

