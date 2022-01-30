
!===============================================================================

module mbomat

	implicit none

	! My name a bomat.  Very nice!
	character(len = *), parameter :: me = "bomat"

	character, parameter :: nullc = char(0), newline = char(10)

	double complex, parameter :: ic = cmplx(0.d0, 1.d0, kind = 8)

	double precision, parameter :: pi = 4.d0 * atan(1.d0)

	! C++ routines from colormapper submodule
	integer, external :: load_colormap, writepng

	integer, parameter :: nrgb = 3

	integer, parameter :: &
			ERR_WRITEPNG = -1, &
			ERR_COLORMAP = -2, &
			ERR_CMD_ARGS = -3

	!********

	type bomat_settings

		! User-facing input settings

		character(len = :), allocatable :: f, fcolormap, colormap, fjson, &
				fdata, fmeta

		integer :: n, np, nx, ny
		integer(kind = 8) :: nsample
		integer, allocatable :: inz(:,:)

		double complex, allocatable :: p(:)

		double precision :: xmin, xmax, ymin, ymax

		! TODO: change arg_all default?  i.e. require -a flag?
		!
		!   - change "-e" to on by default.  it's painful to forget it.  add an
		!     other arg to turn it off during a 1-pass run.  rename to "export"
		!     instead of "eig"
		!
		logical :: arg_all = .true., arg_eig = .false., arg_plot = .false., &
				arg_help = .false.

	end type bomat_settings

	!********

	type bomat_data

		! Internal data passed between subroutines

		! 2D eigenvalue density histogram
		integer, allocatable :: hist(:,:)

		! Pixel values for plot
		character, allocatable :: img(:)

	end type bomat_data

	!********

contains

!===============================================================================

subroutine bomat(io)

	! High-level driver

	integer, intent(out) :: io

	!********

	type(bomat_data) :: d

	type(bomat_settings) :: s

	write(*,*)
	write(*,*) "Starting "//me
	write(*,*)

	call load_args(s, io)
	if (io /= 0 .or. s%arg_help) return

	call load_settings(s, io)
	if (io /= 0) return

	write(*,*) 'Loading colormap ...'
	write(*,*)

	! This sets static C variables which are later used by map()
	io = load_colormap(s%fcolormap//nullc, s%colormap//nullc)
	if (io /= 0) then
		io = ERR_COLORMAP
		return
	end if

	! Load eigenvalues from previous run or calculate from scratch?
	if (s%arg_plot) then

		!print *, 'calling load_eigenvalues ...'
		call load_eigenvalues(s, d, io)
		if (io /= 0) return

	else

		call get_eigenvalues(s, d, io)
		if (io /= 0 .or. s%arg_eig) return

	end if

	call draw_plot(s, d, io)
	if (io /= 0) return

	! TODO: set output based on fjson
	io = writepng(d%img, s%nx, s%ny, "scratch/bomat.png"//nullc)
	if (io /= 0) then
		io = ERR_WRITEPNG
		return
	end if
	write(*,*)

	write(*,*) "Done "//me
	write(*,*)

end subroutine bomat

!===============================================================================

subroutine get_eigenvalues(s, d, io)

	! Calculate eigenvalues.  If running with "-e", export them to a .data file

	type(bomat_settings), intent(in) :: s

	type(bomat_data), intent(out) :: d

	integer, intent(out) :: io

	!********

	character :: jobvl, jobvr

	double complex, allocatable :: a(:,:), w(:), vl(:,:), vr(:,:), work(:)

	double precision :: r, dx, dy, wx, wy, wxmin, wxmax, wymin, wymax, &
			wmargin, dwx, dwy
	double precision, allocatable :: rwork(:)

	integer :: i, lda, info, nseed, ix, iy, ip, ip0, np, ldvl, ldvr, &
			lwork, id, im
	integer(kind = 8) :: is, ist, neigen

	integer, allocatable :: seed(:)

	logical :: maxinit = .false.

	! Not actually thrown from here
	io = 0

	! Matrix
	allocate(a(s%n, s%n))

	! LAPACK solver args

	! zgeev (eigenvectors/eigenvalues)

	lda = s%n
	jobvl = 'N'
	jobvr = 'N'
	ldvl = 1
	ldvr = 1
	lwork = 2 * s%n

	allocate(w(s%n))
	allocate(vl(ldvl, s%n))
	allocate(vr(ldvr, s%n))
	allocate(work(lwork))
	allocate(rwork(2 * s%n))

	call random_seed(size = nseed)
	allocate(seed(nseed))

	! Consistently seeded
	seed = 0
	call random_seed(put = seed)

	!call random_seed(get = seed)  ! different every time, print for repeatability
	!call print_seed(seed)

	! Console progress bar size (in chars) and trackers
	np = 70
	ip0 = 0
	ist = 0

	dx = s%xmax - s%xmin
	dy = s%ymax - s%ymin

	! Eigenvalue density histogram
	allocate(d%hist(s%nx, s%ny))
	d%hist = 0

	neigen = 0

	open(newunit = id, file = s%fdata, access = "stream")

	write(*,*) 'Sampling eigenvalues ...'
	write(*,*)

	! Progress bar header
	write(*,*) '|'//repeat('-', np)//'|'
	write(*, '(a)', advance = 'no') ' |'

	!$OMP parallel do default(shared) &
	!$OMP&    private(a, i, r, info, ix, iy, w, vl, vr, work, rwork, wx, wy)
	do is = 1, s%nsample

		! Progress bar
		!$OMP critical
		ist = ist + 1
		ip = int(dble(np) * ist / s%nsample)
		if (ip > ip0) then
			ip0 = ip
			write(*, '(a)', advance = 'no') '='
		end if
		!$OMP end critical

		! Fill non-zeros of matrix a randomly from population
		a = 0.d0
		do i = 1, size(s%inz, 2)
			call random_number(r)
			a(s%inz(1,i), s%inz(2,i)) = s%p(floor(r * s%np) + 1)
		end do

		!! There is no Fortran edit descriptor for complex numbers :(
		!print *, "a = "
		!print "(16es14.4)", a
		!print *, "a = ", a

		! Call LAPACK to get eigenvalues (but not eigenvectors)
		call zgeev(jobvl, jobvr, s%n, a, lda, w, vl, ldvl, vr, ldvr, &
				work, lwork, rwork, info)

		if (info == 0) then

			do i = 1, s%n

				!! Remove real line.  TODO: make this optional.  Handle in
				!! draw_plot instead?  Could set tolerance in pixels there
				!! instead of complex units

				!if (abs(imagpart(w(i))) > 1.d-6) then

					wx =  real(w(i))
					wy = aimag(w(i))

					!$OMP critical

					! TODO: these could all be done with OMP reduction
					! operators

					neigen = neigen + 1

					if (.not. maxinit) then
						maxinit = .true.
						wxmin = wx
						wxmax = wx
						wymin = wy
						wymax = wy
					else
						wxmin = min(wxmin, wx)
						wxmax = max(wxmax, wx)
						wymin = min(wymin, wy)
						wymax = max(wymax, wy)
					end if

					!$OMP end critical

					if (s%arg_eig) then

						!$OMP critical
						write(id) wx, wy
						!$OMP end critical

					else

						! Pixel coordinates.  Invert so y is up
						ix = int(1    + (wx - s%xmin) / dx * s%nx)
						iy = int(s%ny - (wy - s%ymin) / dy * s%ny)

						if (ix >= 1 .and. ix <= s%nx .and. &
						    iy >= 1 .and. iy <= s%ny) then

							!$OMP critical
							d%hist(ix, iy) = d%hist(ix, iy) + 1
							!$OMP end critical

						end if
					end if
				!end if
			end do
		end if
	end do
	!$OMP end parallel do

	write(*, '(a)') '|'
	write(*,*)

	close(id)

	! Metadata
	open(newunit = im, file = s%fmeta, access = "stream")
	write(im) neigen, wxmin, wxmax, wymin, wymax
	close(im)

	! TODO:  warn if eigenvalues are outside of plot bounds? Not an issue with
	! 2-pass run

	write(*,*) 'Eigenvalue bounds:'
	write(*,*) 'Re in [', wxmin, ', ', wxmax, ']'
	write(*,*) 'Im in [', wymin, ', ', wymax, ']'
	write(*,*)

	dwx = wxmax - wxmin
	dwy = wymax - wymin

	wmargin = 0.05

	!write(*,*) 'With 5% margin:'
	!write(*,*) 'Re in [', wxmin - wmargin*dx, ', ', wxmax + wmargin*dx, ']'
	!write(*,*) 'Im in [', wymin - wmargin*dy, ', ', wymax + wmargin*dy, ']'
	!write(*,*)

end subroutine get_eigenvalues

!=======================================================================

subroutine draw_plot(s, d, io)

	! Draw a Bohemian matrix eigenvalue plot in memory in data variable d

	type(bomat_settings), intent(in) :: s

	type(bomat_data), intent(inout) :: d

	integer, intent(out) :: io

	!********

	character :: rgb(3)

	double precision :: xi, xi0

	integer :: i, j, nnzh, nxy
	integer, allocatable :: idx(:), hist1(:)

	! Not actually thrown from here
	io = 0

	write(*,*) 'Coloring pixels ...'
	write(*,*)

	! TODO: reshape to img(nrgb, s%nx(,*)s%ny)?  C++ call shouldn't care

	!print *, 'd%hist = ', d%hist(1:100,1)

	! TODO: can C allocate idx inside sortidx() function?

	! Note that sortidx treats hist as a rank-1 array and returns 0-based rank-1
	! indices in idx
	nxy = s%nx * s%ny
	allocate(idx(nxy))
	call sortidx(d%hist, nxy, idx)

	!print *, 'idx = ', idx(1: 100)

	allocate(d%img(nxy * nrgb))

	xi0 = 0.d0

	! Number of non-zeros in histogram
	nnzh = count(d%hist /= 0)

	write(*,*) 'Non-zero pixel fraction = ', real(nnzh) / nxy
	write(*,*)

	! TODO:  avoid reshaping.  Do 1d/2d index math instead?  Or can we use
	! move_alloc for different ranks?
	allocate(hist1(nxy))
	hist1 = reshape(d%hist, [nxy])

	j = 0
	do i = 1, nxy

		if (hist1(idx(i) + 1) == 0) then
			xi = 0.d0
		else if (i > 0 .and. hist1(idx(i) + 1) == hist1(idx(i-1) + 1)) then
			j = j + 1
			xi = xi0
		else
			j = j + 1
			xi = dble(j) / nnzh
		end if
		xi0 = xi

		call map(xi, rgb)
		d%img(idx(i) * nrgb + 1: idx(i) * nrgb + 3) = rgb

	end do

	! TODO:  deallocate arrays as soon as they're done being used

	!deallocate(d%hist)

end subroutine draw_plot

!===============================================================================

subroutine load_eigenvalues(s, d, io)

	! Load eigenvalues from a previous run saved in a .data file

	use iso_fortran_env

	type(bomat_settings), intent(inout) :: s

	type(bomat_data), intent(out) :: d

	integer, intent(out) :: io

	!********

	integer, parameter :: nchunk = 1024 ** 2  ! optimal value?

	double precision :: dx, dy, wxmin, wxmax, wymin, wymax, &
			wmargin, dwx, dwy

	!! Parameterized array sizes crash gfortran for large nchunk.  Allocatable
	!! arrays work better, may be stack/heap difference
	!double precision :: ws(2,nchunk)
	double precision, allocatable :: ws(:,:)

	integer :: ix, iy, id, im, ip, ip0, np
	integer(kind = 8) :: i, j, i0, neigen, ist

	!print *, 'starting load_eigenvalues ...'
	!print *, ''

	! Not actually thrown from here
	io = 0

	write(*,*) 'Meta file = "'//s%fmeta//'"'
	write(*,*) 'Data file = "'//s%fdata//'"'
	write(*,*)

	! Metadata
	open(newunit = im, file = s%fmeta, access = "stream")
	read(im) neigen, wxmin, wxmax, wymin, wymax
	close(im)

	! TODO: error for empty file(s)?  Should happen automatically from compiler

	! TODO: include optional bounds override e.g. for zooming
	s%xmin = wxmin
	s%xmax = wxmax
	s%ymin = wymin
	s%ymax = wymax
	dx = s%xmax - s%xmin
	dy = s%ymax - s%ymin
	if (dy > dx) then
		s%ny = int(dy / dx * s%nx)
	else
		s%nx = int(dx / dy * s%ny)
	end if

	! Eigenvalue density histogram
	allocate(d%hist(s%nx, s%ny))
	d%hist = 0

	write(*, '(a,i0)') ' Number of eigenvalues = ', neigen
	write(*, '(a,i0)') ' Eigenvalue chunk size = ', nchunk
	write(*,*)

	write(*,*) 'Loading eigenvalues ...'
	write(*,*)

	open(newunit = id, file = s%fdata, access = "stream")

	!allocate(wxs(nchunk), wys(nchunk), ws(2,nchunk))
	allocate(ws(2,nchunk))

	! Console progress bar size (in chars) and trackers
	np = 70
	ip0 = 0
	ist = 0

	! Progress bar header
	write(*,*) '|'//repeat('-', np)//'|'
	write(*, '(a)', advance = 'no') ' |'

	! TODO: parallelize outer loop.  Will threaded seek/read work?
	i = 0
	do while (i < neigen)
		i0 = i
		i = i + nchunk

		! Progress bar
		ist = i
		ip = int(dble(np) * ist / neigen)
		if (ip > ip0) then
			! TODO: ip > ip0 + 1.  Don't skip characters
			ip0 = ip
			write(*, '(a)', advance = 'no') '='
		end if

		i = min(i, neigen)

		!read(id) (wxs(j), wys(j), j = 1, i - i0)
		!read(id) (ws(:,j), j = 1, i - i0)
		read(id) ws(:, 1: i - i0)

		!$OMP parallel do default(shared) private(j, ix, iy)
		do j = 1, i - i0

			!! Progress bar
			!ist = ist + 1
			!ip = int(dble(np) * ist / neigen)
			!if (ip > ip0) then
			!	ip0 = ip
			!	write(*, '(a)', advance = 'no') '='
			!end if

			!! Remove real line.  TODO: make this optional
			!if (abs(imagpart(w(i))) > 1.d-6) then

				! Pixel coordinates.  Invert so y is up
				ix = int(1    + (ws(1,j) - s%xmin) / dx * s%nx)
				iy = int(s%ny - (ws(2,j) - s%ymin) / dy * s%ny)
				!ix = int(1    + (wxs(j) - s%xmin) / dx * s%nx)
				!iy = int(s%ny - (wys(j) - s%ymin) / dy * s%ny)

				if (ix >= 1 .and. ix <= s%nx .and. &
				    iy >= 1 .and. iy <= s%ny) then

					d%hist(ix, iy) = d%hist(ix, iy) + 1

				end if
			!end if
		end do
		!$OMP end parallel do

	end do

	write(*, '(a)') '|'
	write(*,*)

	close(id)

	write(*,*) 'Eigenvalue bounds:'
	write(*,*) 'Re in [', wxmin, ', ', wxmax, ']'
	write(*,*) 'Im in [', wymin, ', ', wymax, ']'
	write(*,*)

	!! TODO: benchmarking only
	!stop

	dwx = wxmax - wxmin
	dwy = wymax - wymin

	wmargin = 0.05

	!write(*,*) 'With 5% margin:'
	!write(*,*) 'Re in [', wxmin - wmargin*dx, ', ', wxmax + wmargin*dx, ']'
	!write(*,*) 'Im in [', wymin - wmargin*dy, ', ', wymax + wmargin*dy, ']'
	!write(*,*)

end subroutine load_eigenvalues

!===============================================================================

subroutine load_args(s, io)

	! Parse command line arguments

	type(bomat_settings), intent(out) :: s

	integer, intent(out) :: io

	!********

	character :: argv*256
	character(len = :), allocatable :: help
	character(len = *), parameter :: &
			id_h    = "-h"    , &
			id_help = "--help", &
			id_all  = "-a"    , &
			id_plot = "-p"    , &
			id_eig  = "-e"

	integer :: argc, i, ipos

	io = 0

	help = "" &
			//"Usage: "//me//" ["//id_h//"] ["//id_plot//"] ["//id_eig//"] " &
			//"FILE.JSON"//newline &
			//newline &
			//"Calculate Bohemian matrix eigenvalues and export a plot to a " &
			//"PNG file"//newline &
			//newline &
			//"Positional arguments:"//newline &
			//"FILE.JSON   Configuration filename for setting inputs"//newline &
			//newline &
			//"Optional arguments:"//newline &
			//id_h//", "//id_help//"  Show this help message and exit"//newline &
			//id_eig //"          Calculate and export eigenvalues without plotting"//newline &
			//id_plot//"          Plot eigenvalues from previous job"//newline &
			//newline

	argc = command_argument_count()

	!! TODO:  re-enable this check to require input file.  print help too
	!if (argc < 1) then
	!	io = ERR_CMD_ARGS
	!	return
	!end if

	ipos = 0
	i = 0
	do while (i < argc)
		i = i + 1
		call get_command_argument(i, argv)

		! Optional arguments
		if (argv == id_all) then
			s%arg_all = .true.

		else if (argv == id_plot) then
			s%arg_plot = .true.

		else if (argv == id_eig) then
			s%arg_eig = .true.

		else if (argv == id_h .or. argv == id_help) then
			write(*, '(a)') help
			s%arg_help = .true.
			return

		else

			! Positional arguments
			if (ipos == 0) then
				ipos = ipos + 1
				s%fjson = trim(argv)

			else
				write(*,*) 'Warning:  unknown cmd argument '//trim(argv)
				write(*,*)

			end if
		end if
	end do

	!print *, 'arg_all  = ', s%arg_all
	!print *, 'arg_eig  = ', s%arg_eig
	!print *, 'arg_plot = ', s%arg_plot
	!print *, 'fjson    = ', s%fjson
	!write(*, '(a)') help

	write(*,*) "Configuration file = """//s%fjson//""""

	i = len_trim(s%fjson)
	if (i > 0) then
		do while (s%fjson(i: i) /= '.')
			i = i - 1
		end do
	end if
	s%f = s%fjson(1: i - 1)

	s%fdata = s%f//".data"
	s%fmeta = s%f//".meta"

	write(*,*) "Base filename      = """//s%f//""""
	write(*,*)

end subroutine load_args

!===============================================================================

subroutine load_settings(s, io)

	! Load hard-coded settings.  TODO: json config file TBD, this API is not
	! stable

	type(bomat_settings), intent(inout) :: s

	integer, intent(out) :: io

	!********

	integer :: i, j, k, nnonzero
	integer, allocatable :: template(:), t2(:,:)

	! Not actually thrown from here
	io = 0

	! Size of matrices
	s%n = 8

	! Size of population
	s%np = 6

	! Generator sample set.  This is called the "population"
	allocate(s%p(s%np))

	s%p(1) = cmplx( 1.d0,  0.d0 , kind = 8)
	s%p(2) = cmplx(-1.d0,  0.d0 , kind = 8)
	s%p(3) = cmplx( 0.d0,  0.5d0, kind = 8)

	!s%p(1) = cmplx(1.d0, 0.d0, kind = 8)
	!s%p(1) = cmplx(0.5d0, 0.d0, kind = 8)
	!s%p(2) = cmplx(cos(2.d0 * pi / 3.d0), sin(2.d0 * pi / 3.d0), kind = 8)
	!s%p(3) = cmplx(cos(4.d0 * pi / 3.d0), sin(4.d0 * pi / 3.d0), kind = 8)
	!s%p(2) = cmplx(0.d0,  1.d0, kind = 8)
	!s%p(3) = cmplx(0.d0, -1.d0, kind = 8)

	s%p(1) = cmplx(1.d0, 0.d0, kind = 8)
	s%p(2) = cmplx(cos(2.d0 * pi / 3.d0), sin(2.d0 * pi / 3.d0), kind = 8)
	s%p(3) = cmplx(cos(4.d0 * pi / 3.d0), sin(4.d0 * pi / 3.d0), kind = 8)
	s%p(4:6) = -0.1 * s%p(1:3)
	s%p = s%p * exp(ic * 0.5d0 * pi)

	!s%p(1) = cmplx( 1.d0,  0.d0 , kind = 8)
	!s%p(2) = cmplx( 1.d0,  0.1d0, kind = 8)
	!s%p(3) = cmplx(-1.d0,  0.d0 , kind = 8)
	!s%p(4) = cmplx(-1.d0, -0.1d0, kind = 8)

	!s%p(1) = cmplx( 1.d0, -0.3d0, kind = 8)
	!s%p(2) = cmplx( 1.d0,  0.1d0, kind = 8)
	!s%p(3) = cmplx(-1.d0, -0.3d0, kind = 8)
	!s%p(4) = cmplx(-1.d0,  0.1d0, kind = 8)

	!s%p(1) = cmplx( 1.d0,  0.d0, kind = 8)
	!s%p(2) = cmplx( 0.d0,  1.d0, kind = 8)
	!s%p(3) = cmplx(-1.d0,  0.d0, kind = 8)
	!s%p(4) = cmplx( 0.d0, -1.d0, kind = 8)

	!s%p(1) = cmplx( 1.d0,  0.d0 , kind = 8)
	!s%p(2) = cmplx( 0.d0,  0.1d0, kind = 8)
	!s%p(3) = cmplx(-1.d0,  0.d0 , kind = 8)
	!s%p(4) = cmplx( 0.d0, -0.1d0, kind = 8)

	!!s%p(1) = cmplx(1.d0, 0.d0, kind = 8)
	!s%p(1) = cmplx(0.1d0, 0.d0, kind = 8)
	!s%p(2) = cmplx(cos(2.d0 * pi / 5.d0), sin(2.d0 * pi / 5.d0), kind = 8)
	!s%p(3) = cmplx(cos(4.d0 * pi / 5.d0), sin(4.d0 * pi / 5.d0), kind = 8)
	!s%p(4) = cmplx(cos(6.d0 * pi / 5.d0), sin(6.d0 * pi / 5.d0), kind = 8)
	!s%p(5) = cmplx(cos(8.d0 * pi / 5.d0), sin(8.d0 * pi / 5.d0), kind = 8)

	!s%p(1) = cmplx(1.d0, 0.d0, kind = 8)
	!s%p(2) = cmplx(cos(2.d0 * pi / 6.d0), sin(2.d0 * pi / 6.d0), kind = 8)
	!s%p(3) = cmplx(cos(4.d0 * pi / 6.d0), sin(4.d0 * pi / 6.d0), kind = 8)
	!s%p(4) = cmplx(cos(6.d0 * pi / 6.d0), sin(6.d0 * pi / 6.d0), kind = 8)
	!s%p(5) = cmplx(cos(8.d0 * pi / 6.d0), sin(8.d0 * pi / 6.d0), kind = 8)
	!s%p(6) = cmplx(cos(1.d1 * pi / 6.d0), sin(1.d1 * pi / 6.d0), kind = 8)

	!s%p(1) = cmplx(1.d0, 0.d0, kind = 8)
	!s%p(2) = cmplx(cos( 2.d0 * pi / 7.d0), sin( 2.d0 * pi / 7.d0), kind = 8)
	!s%p(3) = cmplx(cos( 4.d0 * pi / 7.d0), sin( 4.d0 * pi / 7.d0), kind = 8)
	!s%p(4) = cmplx(cos( 6.d0 * pi / 7.d0), sin( 6.d0 * pi / 7.d0), kind = 8)
	!s%p(5) = cmplx(cos( 8.d0 * pi / 7.d0), sin( 8.d0 * pi / 7.d0), kind = 8)
	!s%p(6) = cmplx(cos( 1.d1 * pi / 7.d0), sin( 1.d1 * pi / 7.d0), kind = 8)
	!s%p(7) = cmplx(cos(1.2d1 * pi / 7.d0), sin(1.2d1 * pi / 7.d0), kind = 8)

	!print *, "s%p(2) = ", s%p(2)
	!print *, "s%p = ", s%p

	! Non-zero pattern template

	allocate(template(s%n * s%n))

	!! Upper triangle, main diagonal, and first band below diagonal
	!template = [                &
	!	1, 1, 1, 1, 1, 1, 1, 1, &
	!	1, 1, 1, 1, 1, 1, 1, 1, &
	!	0, 1, 1, 1, 1, 1, 1, 1, &
	!	0, 0, 1, 1, 1, 1, 1, 1, &
	!	0, 0, 0, 1, 1, 1, 1, 1, &
	!	0, 0, 0, 0, 1, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 0, 1, 1  &
	!	]

	! Tridiagonal (TODO: add enum options for things like this, Toeplitz,
	! Hermitian, symmetric, skew-symmetric, fully-dense, etc.)

	! 8x8
	template = [                &
		1, 1, 0, 0, 0, 0, 0, 0, &
		1, 1, 1, 0, 0, 0, 0, 0, &
		0, 1, 1, 1, 0, 0, 0, 0, &
		0, 0, 1, 1, 1, 0, 0, 0, &
		0, 0, 0, 1, 1, 1, 0, 0, &
		0, 0, 0, 0, 1, 1, 1, 0, &
		0, 0, 0, 0, 0, 1, 1, 1, &
		0, 0, 0, 0, 0, 0, 1, 1  &
		]

	!! 7x7
	!template = [                &
	!	1, 1, 0, 0, 0, 0, 0, &
	!	1, 1, 1, 0, 0, 0, 0, &
	!	0, 1, 1, 1, 0, 0, 0, &
	!	0, 0, 1, 1, 1, 0, 0, &
	!	0, 0, 0, 1, 1, 1, 0, &
	!	0, 0, 0, 0, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 1, 1  &
	!	]

	!! 3.5-ish diagonal
	!template = [                &
	!	1, 1, 0, 0, 0, 0, 0, 0, &
	!	1, 1, 1, 1, 0, 0, 0, 0, &
	!	0, 1, 1, 1, 0, 0, 0, 0, &
	!	0, 0, 1, 1, 1, 1, 0, 0, &
	!	0, 0, 0, 1, 1, 1, 0, 0, &
	!	0, 0, 0, 0, 1, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 0, 1, 1  &
	!	]

	!! Quad-diagonal.  More non-zeros seem to just fill in the sparse areas and
	!! make the plot more washed-out and less interesting
	!template = [                &
	!	1, 1, 1, 0, 0, 0, 0, 0, &
	!	1, 1, 1, 1, 0, 0, 0, 0, &
	!	0, 1, 1, 1, 1, 0, 0, 0, &
	!	0, 0, 1, 1, 1, 1, 0, 0, &
	!	0, 0, 0, 1, 1, 1, 1, 0, &
	!	0, 0, 0, 0, 1, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 0, 1, 1  &
	!	]

	!! Pentadiagonal.  More non-zeros seem to just fill in the sparse areas and
	!! make the plot more washed-out and less interesting
	!template = [                &
	!	1, 1, 1, 0, 0, 0, 0, 0, &
	!	1, 1, 1, 1, 0, 0, 0, 0, &
	!	1, 1, 1, 1, 1, 0, 0, 0, &
	!	0, 1, 1, 1, 1, 1, 0, 0, &
	!	0, 0, 1, 1, 1, 1, 1, 0, &
	!	0, 0, 0, 1, 1, 1, 1, 1, &
	!	0, 0, 0, 0, 1, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 1, 1, 1  &
	!	]

	!! Fully-dense
	!template = 1

	! Number of random samples to take
	s%nsample = 1000000

	! Colormap file and name
	s%fcolormap = "submodules/colormapper/submodules/colormaps/ColorMaps5.6.0.json"

	!s%colormap = "Inferno (matplotlib)"
	!s%colormap = "Viridis (matplotlib)"
	s%colormap = "Magma (matplotlib)"
	!s%colormap = "Plasma (matplotlib)"

	!s%colormap = "erdc_blue2green_BW"
	!s%colormap = "Black, Blue and White"
	!s%colormap = "Blue Orange (divergent)"

	!s%colormap = "Asymmtrical Earth Tones (6_21b)"

	!! Similar to erdc_blue2green_BW.  Not bad
	!s%colormap = "gist_earth"

	!! Meh

	!s%colormap = "Blue - Green - Orange"
	!s%colormap = "Blue to Red Rainbow"

	!! This is just plain inferior to inferno.  They're very similar, but
	!! inferno has some purple.
	!s%colormap = "Black-Body Radiation"

	!! Custom
	!s%fcolormap = "custom-maps.json"
	!s%colormap = "bisexual"
	!s%colormap = "bisexual-light"

	!! Pretty ugly
	!s%colormap = "cyan-magenta"

	! Image size

	!s%nx = 3840
	!s%nx = 1920
	s%nx = 1080
	!s%nx = 720

	s%ny = s%nx
	!s%ny = 4352
	!s%ny = 4292

	! Image bounds in complex plane.  Careful with aspect ratio

	s%xmin = -2.5
	s%xmax =  2.5
	s%ymin = -2.3
	s%ymax =  2.7

	!s%xmin = -2.73
	!s%xmax =  3.23
	!s%ymin = -2.88
	!s%ymax =  2.88

	!s%xmin = -3.0
	!s%xmax =  3.0
	!s%ymin = -3.0
	!s%ymax =  3.0

	!s%xmin = -3.9
	!s%xmax =  3.9
	!s%ymin = -3.9
	!s%ymax =  3.9

	!********

	! Do some quick basic processing on some of the inputs

	! Mark non-zero locations from template 1/0's matrix

	! Number of non-zeros in matrix
	nnonzero = count(template /= 0)

	!print *, "nnonzero = ", nnonzero

	allocate(t2(s%n, s%n))
	t2 = reshape(template, [s%n, s%n])
	deallocate(template)

	! Indices of non-zeros
	allocate(s%inz(2, nnonzero))

	k = 0
	do i = 1, s%n
	do j = 1, s%n

		! Note the transpose from row-major template to Fortran default
		! column-major
		if (t2(j, i) /= 0) then
			k = k + 1
			s%inz(:, k) = [i, j]
		end if

	end do
	end do
	deallocate(t2)

end subroutine load_settings

!===============================================================================

subroutine print_seed(seed)
  
	! Print RNG seed in a format that can be pasted back into Fortran code (e.g.
	! for repeatability or debugging)

	integer, intent(in), allocatable :: seed(:)

	!********

	character :: com

	integer :: i

	write(*,*)   "   seed = [        &"
	com = ","
	do i = 1, size(seed)
		if (i == size(seed)) com = " "
		write(*,*) "   ", seed(i), com, "  &"
	end do
	write(*,*)   "          ]"

end subroutine print_seed

!===============================================================================

end module mbomat

!===============================================================================

program main

	use mbomat
	implicit none

	integer :: io

	call bomat(io)
	call exit(io)

end program main

!===============================================================================

