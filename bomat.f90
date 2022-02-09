
!===============================================================================

module mbomat

	implicit none

	! My name a bomat.  Very nice!
	character(len = *), parameter :: me = "bomat"

	character, parameter :: nullc = char(0), t = char(9), newline = char(10)

	! JSON keys
	character(len = *), parameter :: &
		population_id = 'Population'      , &
		struct_id     = 'Matrix structure', &
		mat_size_id   = 'Matrix size'     , &
		template_id   = 'Template matrix' , &
		samples_id    = 'Samples'         , &
		img_size_id   = 'Image size'      , &
		margin_id     = 'Margin'          , &
		fcolormap_id  = 'Colormap file'   , &
		colormap_id   = 'Colormap name'

	! Other identifiers
	character(len = *), parameter :: &
		toeplitz_id    = "Toeplitz"      , &
		tridiagonal_id = "Tridiagonal"   , &
		hessenberg_id  = "Hessenberg"    , &
		symmetric_id   = "Symmetric"     , &
		skew_sym_id    = "Skew symmetric", &
		hermitian_id   = "Hermitian"     , &
		dense_id       = "Dense"

	!double complex, parameter :: ic = cmplx(0.d0, 1.d0, kind = 8)
	!double precision, parameter :: pi = 4.d0 * atan(1.d0)

	! C++ routines from colormapper submodule
	integer, external :: load_colormap, writepng

	integer, parameter :: debug = 0

	integer, parameter :: nrgb = 3

	integer, parameter :: &
			ERR_JSON_SYNTAX = -5, &
			ERR_LOAD_JSON   = -4, &
			ERR_CMD_ARGS    = -3, &
			ERR_COLORMAP    = -2, &
			ERR_WRITEPNG    = -1

	!********

	type bomat_settings

		! User-facing input settings

		character(len = :), allocatable :: f, fcolormap, colormap, fjson, &
				fdata, fmeta, fpng

		integer :: n, np, nx = 0, ny
		integer(kind = 8) :: nsample = 0
		integer, allocatable :: inz(:,:)

		double complex, allocatable :: p(:)

		double precision :: xmin, xmax, ymin, ymax, margin = 0.d0

		logical :: size_defined = .false., template_defined = .false.

		! Command-line arguments
		logical :: arg_in_core = .false., arg_eig = .false., &
				arg_plot = .false., arg_help = .false.

		! TODO:  add n-diagonal option, with n a variable int, e.g. tridiagonal,
		! pentadiagonal, etc. (generalization of tridiagonal)

		! Matrix structure options
		logical :: toeplitz = .false., tridiagonal = .false., &
				hessenberg = .false., symmetric = .false., &
				skew_sym = .false., hermitian = .false., &
				dense = .false.

		contains
			procedure :: print => bomat_settings_print

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

	! Calculate eigenvalues from scratch?
	if (.not. s%arg_plot) then

		call calc_eigenvalues(s, d, io)
		if (io /= 0 .or. s%arg_eig) return

	end if

	! Load eigenvalues from disk from a previous calc_eigenvalues() call?
	if (s%arg_plot .or. .not. s%arg_in_core) then

		!print *, 'calling load_eigenvalues ...'
		call load_eigenvalues(s, d, io)
		if (io /= 0) return

	end if

	call draw_plot(s, d, io)
	if (io /= 0) return

	io = writepng(d%img, s%nx, s%ny, s%fpng//nullc)
	if (io /= 0) then
		io = ERR_WRITEPNG
		return
	end if
	write(*,*)

	write(*,*) "Done "//me
	write(*,*)

end subroutine bomat

!===============================================================================

subroutine calc_eigenvalues(s, d, io)

	! Calculate eigenvalues.  Unless running with in-core "-i" flag, export them
	! to a .data file

	!use iso_fortran_env

	type(bomat_settings), intent(in) :: s

	type(bomat_data), intent(out) :: d

	integer, intent(out) :: io

	!********

	character :: jobvl, jobvr

	double complex, allocatable :: a(:,:), w(:), vl(:,:), vr(:,:), work(:)

	double precision :: dx, dy, wx, wy, wxmin, wxmax, wymin, wymax
	double precision, allocatable :: rwork(:)

	integer :: i, lda, info, nseed, ix, iy, ip, ip0, np, ldvl, ldvr, &
			lwork, id, im
	integer(kind = 8) :: is, ist, neigen

	integer, allocatable :: seed(:)

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

	if (s%arg_in_core) then
		! Eigenvalue density histogram
		allocate(d%hist(s%nx, s%ny))
		d%hist = 0
	end if

	neigen = 0

	open(newunit = id, file = s%fdata, access = "stream")

	write(*,*) 'Sampling eigenvalues ...'
	write(*,*)

	! Progress bar header
	write(*,*) '|'//repeat('-', np)//'|'
	write(*, '(a)', advance = 'no') ' |'

	!$OMP parallel do default(shared) schedule(static) &
	!$OMP&    reduction(min: wxmin, wymin) reduction(+: neigen) &
	!$OMP&    reduction(max: wxmax, wymax) &
	!$OMP&    private(a, i, info, ix, iy, w, vl, vr, work, rwork, wx, wy)
	do is = 1, s%nsample

		! Progress bar
		!$OMP critical
		ist = ist + 1
		ip = int(dble(np) * ist / s%nsample)
		if (ip > ip0) then
			ip0 = ip
			write(*, '(a)', advance = 'no') '='
			!flush(output_unit)
		end if
		!$OMP end critical

		!$OMP critical
		! TODO: use abstract interfaces.  c.f. blog
		if (s%toeplitz) then
			a = random_toeplitz(s)
		else
			a = random_matrix(s)
		end if
		!$OMP end critical

		!! There is no Fortran edit descriptor for complex numbers :(
		!print *, "a = "
		!print "(16es14.4)", a
		!!print *, "a = ", a

		! Call LAPACK to get eigenvalues (but not eigenvectors)
		call zgeev(jobvl, jobvr, s%n, a, lda, w, vl, ldvl, vr, ldvr, &
				work, lwork, rwork, info)

		if (info == 0) then

			do i = 1, s%n

				wx =  real(w(i))
				wy = aimag(w(i))

				neigen = neigen + 1
				wxmin = min(wxmin, wx)
				wxmax = max(wxmax, wx)
				wymin = min(wymin, wy)
				wymax = max(wymax, wy)

				if (.not. s%arg_in_core) then

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
			end do
		end if
	end do
	!$OMP end parallel do

	write(*, '(a)') repeat('=', np - ip0)//'|'
	write(*,*)

	close(id)

	! Metadata
	open(newunit = im, file = s%fmeta, access = "stream")
	write(im) neigen, wxmin, wxmax, wymin, wymax
	close(im)

	write(*,*) 'Eigenvalue bounds:'
	write(*,*) 'Re in [', wxmin, ', ', wxmax, ']'
	write(*,*) 'Im in [', wymin, ', ', wymax, ']'
	write(*,*)

end subroutine calc_eigenvalues

!=======================================================================

function random_matrix(s) result(a)

	! Make a random matrix with no specific structure with entries sampled from
	! the population

	type(bomat_settings), intent(in) :: s

	double complex :: a(s%n, s%n)

	!********

	double precision :: r

	integer :: i

	! Fill non-zeros of matrix a randomly from population
	a = 0.d0
	do i = 1, size(s%inz, 2)
		call random_number(r)
		a(s%inz(1,i), s%inz(2,i)) = s%p(floor(r * s%np) + 1)
	end do

end function random_matrix

!=======================================================================

function random_toeplitz(s) result(a)

	! Make a random Toeplitz matrix.  If you don't know what Toeplitz is, a 4x4
	! is like this:
	!
	! [
	!     e f g h
	!     d e f g
	!     c d e f
	!     b c d e
	! ]

	type(bomat_settings), intent(in) :: s

	double complex :: a(s%n, s%n)

	!********

	double complex :: diags(-(s%n - 1): s%n - 1)

	double precision :: r

	integer :: i, j, k

	do i = -(s%n - 1), s%n - 1
		call random_number(r)
		diags(i) = s%p(floor(r * s%np) + 1)
	end do

	! Template-based sparse Toeplitz
	a = 0.d0
	do k = 1, size(s%inz, 2)
		i = s%inz(1,k)
		j = s%inz(2,k)
		a(i, j) = diags(i - j)
	end do

end function random_toeplitz

!=======================================================================

subroutine draw_plot(s, d, io)

	! Draw a Bohemian matrix eigenvalue plot in memory in data variable d

	type(bomat_settings), intent(in) :: s

	type(bomat_data), intent(inout) :: d

	integer, intent(out) :: io

	!********

	character :: rgb(3)

	double precision :: xi, xi0

	integer :: i, j, ix, iy, nnzh, nxy, h, h0
	integer, allocatable :: idx(:)!, hist1(:)

	! Not actually thrown from here
	io = 0

	write(*,*) 'Coloring pixels ...'
	write(*,*)

	!print *, 'd%hist = ', d%hist(1:100,1)

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

	j = 0
	h0 = 0
	do i = 1, nxy

		! Convert 1D 0-based index to 2D 1-based indices.  This is more complex
		! than reshaping the hist array, but it doesn't require copying
		! a potentially large amount of data
		iy = idx(i) / s%nx + 1
		ix = mod(idx(i), s%nx) + 1

		h = d%hist(ix,iy)

		if (h == 0) then
			xi = 0.d0
		else if (i > 0 .and. h == h0) then
			j = j + 1
			xi = xi0
		else
			j = j + 1
			xi = dble(j) / nnzh
		end if
		xi0 = xi

		call map(xi, rgb)
		d%img(idx(i) * nrgb + 1: idx(i) * nrgb + 3) = rgb

		h0 = h
	end do

	! Deallocate arrays as soon as they're done being used?  This is almost the
	! last subroutine anyway, other than writepng.
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
			dwx, dwy

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

	! TODO: include optional bounds override e.g. for zooming or margin

	dwx = wxmax - wxmin
	dwy = wymax - wymin

	s%xmin = wxmin - s%margin * dwx
	s%xmax = wxmax + s%margin * dwx
	s%ymin = wymin - s%margin * dwy
	s%ymax = wymax + s%margin * dwy

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
		if (ip > ip0 .and. ip < np) then
			write(*, '(a)', advance = 'no') repeat('=', ip - ip0)
			ip0 = ip
		end if

		i = min(i, neigen)

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

			! Pixel coordinates.  Invert so y is up
			ix = int(1    + (ws(1,j) - s%xmin) / dx * s%nx)
			iy = int(s%ny - (ws(2,j) - s%ymin) / dy * s%ny)

			if (ix >= 1 .and. ix <= s%nx .and. &
			    iy >= 1 .and. iy <= s%ny) then

				!$OMP critical
				d%hist(ix, iy) = d%hist(ix, iy) + 1
				!$OMP end critical

			end if
		end do
		!$OMP end parallel do

	end do


	write(*, '(a)') repeat('=', np - ip0)//'|'
	write(*,*)

	close(id)

	write(*,*) 'Eigenvalue bounds:'
	write(*,*) 'Re in [', wxmin, ', ', wxmax, ']'
	write(*,*) 'Im in [', wymin, ', ', wymax, ']'
	write(*,*)

	!! TODO: benchmarking only
	!stop

end subroutine load_eigenvalues

!===============================================================================

subroutine load_args(s, io)

	! Parse command line arguments

	type(bomat_settings), intent(out) :: s

	integer, intent(out) :: io

	!********

	character :: argv*256
	character(len = :), allocatable :: help, help_short
	character(len = *), parameter :: &
			id_file = "FILE.JSON" , &
			id_h    = "-h"        , &
			id_help = "--help"    , &
			id_inco = "-i"        , &
			id_plot = "-p"        , &
			id_eig  = "-e"

	integer :: argc, i, ipos

	io = 0

	help_short = "" &
			//"Usage: "//me//" ["//id_h//"] ["//id_plot//"] ["//id_eig//"] " &
			//id_file//newline &
			//newline &
			//"Calculate Bohemian matrix eigenvalues and export a plot to a " &
			//"PNG file"//newline &
			//newline &
			//"Positional arguments:"//newline &
			//id_file//"   Configuration filename for setting inputs"//newline &
			//newline &
			//"Optional arguments:"//newline &
			//id_h//", "//id_help//"  Show this help message and exit"//newline &
			//id_eig //"          Calculate and export eigenvalues without plotting"//newline &
			//id_plot//"          Plot eigenvalues from previous job"//newline &
			//newline

	! TODO: add structures to help text, e.g. Hessenberg, Toeplitz, ...
	help = help_short &
			//"Sample "//id_file//" contents are like this:"//newline &
			//'{'//newline &
			//newline &
			//'	# This is a comment (non-standard JSON extension)'//newline &
			//newline &
			//'	# Complex numbers'//newline &
			//'	"'//population_id//'":'//newline &
			//'	['//newline &
			//'		# Re, Im'//newline &
			//'		 0        ,  1  ,'//newline &
			//'		-0.8660254, -0.5,'//newline &
			//'		 0.8660254, -0.5'//newline &
			//'	],'//newline &
			//newline &
			//'	# Integers-only.  Zeros will remain 0, non-zeros will be sampled randomly'//newline &
			//'	# from population.  This JSON array is rank-1, but it is reshaped in'//newline &
			//'	# a row-major sense into a rank-2 matrix'//newline &
			//'	"'//template_id//'":'//newline &
			//'	['//newline &
			//'		1, 1, 0, 0, 0, 0, 0, 0,'//newline &
			//'		1, 1, 1, 0, 0, 0, 0, 0,'//newline &
			//'		0, 1, 1, 1, 0, 0, 0, 0,'//newline &
			//'		0, 0, 1, 1, 1, 0, 0, 0,'//newline &
			//'		0, 0, 0, 1, 1, 1, 0, 0,'//newline &
			//'		0, 0, 0, 0, 1, 1, 1, 0,'//newline &
			//'		0, 0, 0, 0, 0, 1, 1, 1,'//newline &
			//'		0, 0, 0, 0, 0, 0, 1, 1'//newline &
			//'	],'//newline &
			//newline &
			//'	"'//samples_id//'": 1000000,'//newline &
			//'	"'//img_size_id//'": 1920,'//newline &
			//newline &
			//'	"'//fcolormap_id//'": "submodules/colormapper/submodules/colormaps/ColorMaps5.6.0.json",'//newline &
			//'	"'//colormap_id//'": "Magma (matplotlib)"'//newline &
			//newline &
			//'}'//newline &
			//newline &
			//"For more, see:"//newline &
			//newline &
			//"	https://github.com/JeffIrwin/bomat"//newline &
			//newline

	argc = command_argument_count()

	if (argc < 1) then
		write(*,*) 'Error: required positional argument '//id_file &
				//' is not defined'
		write(*,*)
		write(*, '(a)') help_short
		io = ERR_CMD_ARGS
		return
	end if

	ipos = 0
	i = 0
	do while (i < argc)
		i = i + 1
		call get_command_argument(i, argv)

		! Optional arguments
		if (argv == id_inco) then

			! TODO: not documented and not recommended.  Plot size/bounds cannot
			! be controlled (yet)
			s%arg_in_core = .true.

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
	s%fpng  = s%f//".png"

	write(*,*) "Base filename      = """//s%f//""""
	write(*,*)

end subroutine load_args

!===============================================================================

subroutine load_settings(s, io)

	! Load input settings from JSON config file

	use json_module

	type(bomat_settings), intent(inout) :: s

	integer, intent(out) :: io

	!********

	integer :: i, j, k, nnonzero

	type(json_file) :: json

	! Not actually thrown from here
	io = 0

	write(*,*) 'Loading json file "'//s%fjson//'" ...'
	write(*,*)

	! Set defaults. TODO pick default values and error out if required JSON keys
	! are not defined.
	s%n = 0
	s%fcolormap = ""
	s%colormap  = ""

	! Image bounds in complex plane.  Careful with aspect ratio.  These are
	! reset automatically later in a 2-pass run.  Check thesis.  Is ~sqrt(n)
	! better?
	s%xmin = -3.0
	s%xmax =  3.0
	s%ymin = -3.0
	s%ymax =  3.0

	call json%initialize()
	call json%load(filename = s%fjson)
	if (json%failed()) then
		write(*,*) 'Error:'
		write(*,*) 'Could not load file "'//s%fjson//'"'
		write(*,*)
		call json%print_error_message()
		write(*,*)
		io = ERR_LOAD_JSON
		return
	end if

	!! One number per line in arrays.  Not great for template matrix
	!call json%print()

	call json%traverse(traverse_bomat_json)
	if (io /= 0) return

	! Mark non-zeros for string-specified structures
	!
	! TODO: implement other options besides Hessenberg
	if (s%hessenberg) then

		! Upper Hessenberg.  Lower Hessenberg is not implemented, but would it
		! make a difference?
		nnonzero = s%n * (s%n + 1) / 2 + s%n - 1

		! Indices of non-zeros
		allocate(s%inz(2, nnonzero))

		! Mark non-zero locations from template 1/0's matrix
		k = 0
		do j = 1, s%n
		do i = 1, s%n

			if (i <= j + 1) then
				k = k + 1
				s%inz(:, k) = [i, j]
			end if

		end do
		end do

	else if (s%tridiagonal) then

		nnonzero = s%n + 2 * (s%n - 1)

		! Indices of non-zeros
		allocate(s%inz(2, nnonzero))

		! Mark non-zero locations from template 1/0's matrix
		k = 0
		do j = 1, s%n
		do i = 1, s%n

			if (i <= j + 1 .and. i >= j - 1) then
				k = k + 1
				s%inz(:, k) = [i, j]
			end if

		end do
		end do

	end if

	call s%print()

	!print *, 'size(s%inz) = ', size(s%inz)
	!print *, 'size(s%inz)  = ', size(s%inz)
	!print *, 's%fcolormap = ', s%fcolormap

contains

!===============================================================================

subroutine traverse_bomat_json(json, p, finished)

	! Parse the bomat json config file, called iteratively from load_settings()
	! via the traverse callback.  This is a nested subroutine because it needs
	! to access settings s from caller, but json traverse signature cannot take
	! extra arguments.

	! Based on the traverser here:
	!
	!     https://github.com/jacobwilliams/json-fortran/issues/204
	!

	use json_module

	class(json_core), intent(inout)       :: json
	type(json_value), pointer, intent(in) :: p
	logical(json_LK), intent(out)         :: finished

	!********

	character(kind=json_CK, len=:), allocatable :: key, sval

	integer(json_IK) :: ival, ncount, ij
	integer, allocatable :: template(:), t2(:,:)

	logical(json_LK) :: found

	real(json_RK) :: rvalx, rvaly

	type(json_value), pointer :: pc

	! Get the name of the key and the type of its value
	call json%info(p, name = key)

	! TODO: parse bounds, image y, etc.  Not required for 2-pass run

	!print *, 'key = "'//key//'"'

	if (key == fcolormap_id) then

		! Colormap file
		call json%get(p, '@', s%fcolormap)

	else if (key == colormap_id) then

		! Colormap name
		call json%get(p, '@', s%colormap)

	else if (key == samples_id) then

		! Number of random samples to take.  json-fortran doesn't seem to support
		! integer*8.  Cast from real instead to handle large values >~ 2 billion
		call json%get(p, '@', rvalx)
		!call json%get(p, '@', s%nsample)
		s%nsample = int(rvalx, kind = 8)

	else if (key == margin_id) then

		! Margin (as a fraction of 1) for plot boundary
		call json%get(p, '@', s%margin)

	else if (key == img_size_id) then

		! Image size.  In a 2-pass run, one of these is automatically
		! resized later for an appropriate aspect ratio
		call json%get(p, '@', s%nx)
		s%ny = s%nx

	else if (key == population_id) then

		! Generator sample set.  This is called the "population"
		ncount = json%count(p)

		! Size of population (real + imaginary pairs)
		s%np = ncount / 2

		allocate(s%p(s%np))

		do ij = 0, s%np - 1

			call json%get_child(p, ij*2 + 1, pc, found)
			call json%get(pc, '@', rvalx)
			call json%get_child(p, ij*2 + 2, pc, found)
			call json%get(pc, '@', rvaly)

			s%p(ij + 1) = cmplx(rvalx, rvaly, kind = 8)

		end do

		!print *, 's%p = ', s%p

	else if (key == mat_size_id) then

		s%size_defined = .true.
		call json%get(p, '@', s%n)

	else if (key == struct_id) then

		! Structure string array
		ncount = json%count(p)

		do ij = 1, ncount

			call json%get_child(p, ij, pc, found)
			call json%get(pc, '@', sval)

			!print *, 'sval = ', sval

			! TODO:  implement other structures.  Check for conflicts, e.g.
			! tridiagonal and Hessenberg.  Explicit template conflicts with
			! tridiagonal, Hessenberg, and dense.

			if (sval == toeplitz_id) then
				s%toeplitz = .true.

			else if (sval == hessenberg_id) then
				s%hessenberg = .true.

			else if (sval == tridiagonal_id) then
				s%tridiagonal = .true.

			else if (sval == symmetric_id) then
				s%symmetric = .true.

			else if (sval == skew_sym_id) then
				s%skew_sym = .true.

			else if (sval == hermitian_id) then
				s%hermitian = .true.

			else if (sval == dense_id) then
				s%dense = .true.

			else
				write(*,*) 'Error: unknown '//struct_id//' string'
				write(*,*) 'String: '//sval
				write(*,*)
				finished = .true.
				io = ERR_JSON_SYNTAX

			end if

		end do

	else if (key == template_id) then

		! Non-zero pattern template matrix
		s%template_defined = .true.

		ncount = json%count(p)

		! Size of matrices
		if (s%size_defined) then

			if (s%n * s%n /= ncount) then
				write(*,*) 'Error: '//template_id//' size does not match ' &
						//mat_size_id
				write(*,*)
				finished = .true.
				io = ERR_JSON_SYNTAX
			end if

		else

			! This is dangerous.  There's no Fortran integer sqrt
			s%n = int(sqrt(dble(ncount)))

			!print *, 's%n = ', s%n

			if (s%n * s%n /= ncount) then
				write(*,*) 'Error: '//template_id//' is not square'
				write(*,*)
				finished = .true.
				io = ERR_JSON_SYNTAX
			end if

		end if

		allocate(template(s%n * s%n))

		do ij = 1, ncount
			call json%get_child(p, ij, pc, found)
			call json%get(pc, '@', ival)
			!print *, 'ival = ', ival
			template(ij) = ival
		end do

		! Number of non-zeros in matrix
		nnonzero = count(template /= 0)

		allocate(t2(s%n, s%n))
		t2 = reshape(template, [s%n, s%n])
		deallocate(template)

		! Indices of non-zeros
		allocate(s%inz(2, nnonzero))

		! Mark non-zero locations from template 1/0's matrix
		k = 0
		do j = 1, s%n
		do i = 1, s%n

			! Note the transpose from row-major template to Fortran default
			! column-major
			if (t2(j, i) /= 0) then
				k = k + 1
				s%inz(:, k) = [i, j]
			end if

		end do
		end do
		deallocate(t2)

		! TODO: if possible, update pointer p to tail of array, so this callback
		! doesn't reiterate through each element individually
		!
		! json_get_tail() ?

	else if (key /= "" .and. key /= s%fjson) then

		! TODO: consider making this an error, unless running with a "loose
		! syntax" cmd arg.  Same idea for unknown cmd args.  Allowing unknown
		! keys is good for future compatibility but bad for users who might make
		! typos.

		write(*,*) 'Warning:  unknown JSON key'
		write(*,*) 'Key    :"'//key//'"'
		write(*,*)

	end if

	if (json%failed()) then

		write(*,*) 'Error:'
		write(*,*) 'Could not load file "'//s%fjson//'"'
		write(*,*)
		call json%print_error_message()
		write(*,*)
		finished = .true.
		io = ERR_JSON_SYNTAX

	end if

	! always false, since we want to traverse all nodes:
	finished = .false.

end subroutine traverse_bomat_json

!===============================================================================

end subroutine load_settings

!===============================================================================

subroutine bomat_settings_print(s)

	use iso_fortran_env

	class(bomat_settings), intent(in) :: s

	character(len = *), parameter :: dlm = ': ', vdm = ', '

	!********

	integer :: i, j, iu
	integer, allocatable :: t2(:,:)

	! This could be an optional argument for printing to file
	iu = output_unit

	! Settings like fjson, arg_plot, etc. are not printed.  Only JSON file
	! settings

	write(iu, *) fcolormap_id, dlm, s%fcolormap
	write(iu, *) colormap_id , dlm, s%colormap
	write(iu, *) img_size_id , dlm, s%nx
	write(iu, *) margin_id   , dlm, s%margin
	write(iu, *) samples_id  , dlm, s%nsample

	if (s%template_defined .or. debug > 0) then

		! Recreate template matrix just to print it
		allocate(t2(s%n, s%n))
		t2 = 0
		do i = 1, size(s%inz, 2)
			t2(s%inz(1,i), s%inz(2,i)) = 1
		end do

		write(iu, *) template_id, dlm
		write(iu, *) '['
		do i = 1, s%n
			write(iu, '(a)', advance = 'no') t
			do j = 1, s%n
				write(iu, '(i0,a)', advance = 'no') t2(i,j), vdm
			end do
		write(iu, *)
		end do
		write(iu, *) ']'

	end if

	if (s%size_defined) then
		write(iu, *) mat_size_id, dlm, s%n
	end if

	write(iu, *) struct_id, dlm
	write(iu, *) '['

	if (s%toeplitz)    write(iu, *) t, toeplitz_id   , vdm
	if (s%tridiagonal) write(iu, *) t, tridiagonal_id, vdm
	if (s%hessenberg)  write(iu, *) t, hessenberg_id , vdm
	if (s%symmetric)   write(iu, *) t, symmetric_id  , vdm
	if (s%skew_sym)    write(iu, *) t, skew_sym_id   , vdm
	if (s%hermitian)   write(iu, *) t, hermitian_id  , vdm
	if (s%dense)       write(iu, *) t, dense_id      , vdm

	write(iu, *) ']'

	write(iu, *) population_id, dlm
	write(iu, *) '['
	do i = 1, s%np
		write(iu, '(a,es14.4,a,es14.4,a)') t, real(s%p(i)), vdm, &
		                                     aimag(s%p(i)), vdm
	end do
	write(iu, *) ']'

	write(iu, *)

end subroutine bomat_settings_print

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

