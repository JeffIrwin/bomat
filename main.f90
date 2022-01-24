
! From:
!
!     http://icl.cs.utk.edu/lapack-for-windows/lapack/
!

program main

	integer, parameter :: n=2, lda=2, nrhs=1, ldb=2

	double precision :: a(lda,n), b(ldb,nrhs)

	integer :: ipiv(n), info, i, j

	a(1,1) = 1
	a(1,2) = 2
	a(2,1) = 3
	a(2,2) = 4

	b(1,1) = 5
	b(2,1) = 6

	print *, "Hello world"

	! Solve linear system a*x = b, for x
	call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

	print *, "info = ", info
	print *, "a (factored) = ", a

	! Expected output: x = -4.0  4.5
	print *, "x = ", b

	print *, "done main"

end program main

