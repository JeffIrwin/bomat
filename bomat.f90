
!===============================================================================

module mbomat

	implicit none

	! My name a bomat.  Very nice!
	character(len = *), parameter :: me = "bomat"

	character, parameter :: nullc = char(0), t = char(9), newline = char(10)

	! JSON keys
	character(len = *), parameter :: &
		population_id = 'Population'     , &
		template_id   = 'Template matrix', &
		samples_id    = 'Samples'        , &
		img_size_id   = 'Image size'     , &
		fcolormap_id  = 'Colormap file'  , &
		colormap_id   = 'Colormap name'

	double complex, parameter :: ic = cmplx(0.d0, 1.d0, kind = 8)

	double precision, parameter :: pi = 4.d0 * atan(1.d0)

	! C++ routines from colormapper submodule
	integer, external :: load_colormap, writepng

	integer, parameter :: nrgb = 3

	integer, parameter :: &
			ERR_LOAD_JSON = -4, &
			ERR_CMD_ARGS  = -3, &
			ERR_COLORMAP  = -2, &
			ERR_WRITEPNG  = -1

	!********

	type bomat_settings

		! User-facing input settings

		character(len = :), allocatable :: f, fcolormap, colormap, fjson, &
				fdata, fmeta, fpng

		integer :: n, np, nx, ny
		integer(kind = 8) :: nsample
		integer, allocatable :: inz(:,:)

		double complex, allocatable :: p(:)

		double precision :: xmin, xmax, ymin, ymax

		logical :: arg_in_core = .false., arg_eig = .false., &
				arg_plot = .false., arg_help = .false.

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

subroutine bomat_settings_print(s)

	use iso_fortran_env

	class(bomat_settings), intent(in) :: s

	character(len = *), parameter :: dlm = ': ', vdm = ', '

	!********

	integer :: i, j, iu
	integer, allocatable :: t2(:,:)

	! This could be an optional argument for printing to file
	iu = output_unit

	write(iu, *) fcolormap_id, dlm, s%fcolormap
	write(iu, *) colormap_id , dlm, s%colormap
	write(iu, *) img_size_id , dlm, s%nx
	write(iu, *) samples_id  , dlm, s%nsample

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
	! off-core run

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

end subroutine calc_eigenvalues

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
		if (ip > ip0 .and. ip < np) then
			write(*, '(a)', advance = 'no') repeat('=', ip - ip0)
			ip0 = ip
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

	! TODO: finish any remainder progress

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

	! TODO: use json key parameters in help text

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

	help = help_short &
			//"Sample "//id_file//" contents are like this:"//newline &
			//'{'//newline &
			//newline &
			//'	# This is a comment (non-standard JSON extension)'//newline &
			//newline &
			//'	# Complex numbers'//newline &
			//'	"Population":'//newline &
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
			//'	"Template matrix":'//newline &
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
			//'	"Samples": 1000000,'//newline &
			//'	"Image size": 1920,'//newline &
			//newline &
			//'	"Colormap file": "submodules/colormapper/submodules/colormaps/ColorMaps5.6.0.json",'//newline &
			//'	"Colormap name": "Magma (matplotlib)"'//newline &
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

	call json%initialize()
	call json%load(filename = s%fjson)
	if (json%failed()) then
		write(*,*) 'Error:'
		write(*,*) 'Could not load file "'//s%fjson//'"'
		write(*,*)
		call json%print_error_message()
		io = ERR_LOAD_JSON
		return
	end if

	!call json%print()

	call json%traverse(traverse_bomat_json)

	call s%print()

	!print *, 'size(s%inz) = ', size(s%inz)
	!print *, 'size(s%inz)  = ', size(s%inz)
	!print *, 's%fcolormap = ', s%fcolormap

	! Tridiagonal (TODO: add enum options for things like this, Toeplitz,
	! Hermitian, symmetric, skew-symmetric, fully-dense, etc.)

	! Image bounds in complex plane.  Careful with aspect ratio.  These are
	! reset automatically later in a 2-pass run.  Check thesis.  Is ~sqrt(n)
	! better?
	s%xmin = -3.0
	s%xmax =  3.0
	s%ymin = -3.0
	s%ymax =  3.0





	!s%p(1) = cmplx( 1.d0,  0.d0 , kind = 8)
	!s%p(2) = cmplx(-1.d0,  0.d0 , kind = 8)
	!s%p(3) = cmplx( 0.d0,  0.5d0, kind = 8)

	!s%p(1) = cmplx(1.d0, 0.d0, kind = 8)
	!s%p(1) = cmplx(0.5d0, 0.d0, kind = 8)
	!s%p(2) = cmplx(cos(2.d0 * pi / 3.d0), sin(2.d0 * pi / 3.d0), kind = 8)
	!s%p(3) = cmplx(cos(4.d0 * pi / 3.d0), sin(4.d0 * pi / 3.d0), kind = 8)
	!s%p(2) = cmplx(0.d0,  1.d0, kind = 8)
	!s%p(3) = cmplx(0.d0, -1.d0, kind = 8)

	!s%p(1) = cmplx(1.d0, -1.d0, kind = 8)
	!s%p(2) = cmplx(1.d0,  0.d0, kind = 8)
	!s%p(3) = cmplx(1.d0,  1.d0, kind = 8)

	!s%p(1) = cmplx(1.d0, 0.d0, kind = 8)
	!s%p(2) = cmplx(cos(2.d0 * pi / 3.d0), sin(2.d0 * pi / 3.d0), kind = 8)
	!s%p(3) = cmplx(cos(4.d0 * pi / 3.d0), sin(4.d0 * pi / 3.d0), kind = 8)
	!s%p(4:6) = -0.1 * s%p(1:3)
	!s%p = s%p * exp(ic * 0.5d0 * pi)

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

	!! 8x8
	!template = [                &
	!	1, 1, 0, 0, 0, 0, 0, 0, &
	!	1, 1, 1, 0, 0, 0, 0, 0, &
	!	0, 1, 1, 1, 0, 0, 0, 0, &
	!	0, 0, 1, 1, 1, 0, 0, 0, &
	!	0, 0, 0, 1, 1, 1, 0, 0, &
	!	0, 0, 0, 0, 1, 1, 1, 0, &
	!	0, 0, 0, 0, 0, 1, 1, 1, &
	!	0, 0, 0, 0, 0, 0, 1, 1  &
	!	]

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

	!! Colormap file and name
	!s%fcolormap = "submodules/colormapper/submodules/colormaps/ColorMaps5.6.0.json"

	!s%colormap = "Inferno (matplotlib)"
	!s%colormap = "Viridis (matplotlib)"
	!s%colormap = "Magma (matplotlib)"
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

	integer(json_IK) :: var_type, ival, ncount, ij
	integer :: i, j, k, nnonzero
	integer, allocatable :: template(:), t2(:,:)

	logical(json_LK) :: found

	real(json_RK) :: rvalx, rvaly

	type(json_value), pointer :: pc

	! Get the name of the key and the type of its value
	call json%info(p, name = key, var_type = var_type)

	! TODO: parse bounds, image y, etc.  Not required for 2-pass run

	! TODO: remove var_type conditions, they aren't necessary.  Simply use
	! %get() inside each key condition instead

	! String values
	if (var_type == json_string) then

		call json%get(p, '@', sval)

		!print *, 'key = "'//key//'"'
		!print *, 'val = "'//sval//'"'

		if (key == fcolormap_id) then

			! Colormap file
			s%fcolormap = sval

		else if (key == colormap_id) then

			! Colormap name
			s%colormap = sval

		else if (key /= "") then

			! TODO: can this be done generically for any type?
			write(*,*) 'Warning:  unknown string JSON key'
			write(*,*) 'Key    : "'//key//'"'
			write(*,*) 'Value  : "'//sval//'"'
			write(*,*)

		end if

	! Integer values
	else if (var_type == json_integer) then

		call json%get(p, '@', ival)

		!print *, 'key = "'//key//'"'
		!print *, 'val = ', ival

		if (key == samples_id) then

			! Number of random samples to take
			s%nsample = ival

		else if (key == img_size_id) then

			! Image size.  In a 2-pass run, one of these is automatically
			! resized later for an appropriate aspect ratio
			s%nx = ival
			s%ny = ival

		else if (key /= "") then
			write(*,*) 'Warning:  unknown integer JSON key'
			write(*,*) 'Key    : "'//key//'"'
			write(*,*) 'Value  : ', ival
			write(*,*)

		end if

	else if (var_type == json_array) then

		ncount = json%count(p)

		!print *, 'array key = "'//key//'"'
		!print *, 'ncount = ', ncount
		!print *, ''

		if (key == population_id) then

			! Generator sample set.  This is called the "population"

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

		else if (key == template_id) then

			! Non-zero pattern template matrix

			! Size of matrices
			!
			! TODO: this is dangerous.  There's no Fortran integer sqrt
			s%n = int(sqrt(dble(ncount)))

			!print *, 's%n = ', s%n

			! TODO:  add a "Matrix size" input to fix this.  It will be required
			! for e.g. "tridiagonal" input
			if (s%n * s%n /= ncount) then
				write(*,*) 'Error: matrix is not square'
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

		else if (key /= "") then
			write(*,*) 'Warning:  unknown array JSON key'
			write(*,*) 'Key    : "'//key//'"'
			!write(*,*) 'Value  : "'//sval//'"'
			write(*,*)

		end if

		! TODO: if possible, update pointer p to tail of array, so this callback
		! doesn't reiterate through each element individually
		!
		! json_get_tail() ?

	else if (key /= "" .and. key /= s%fjson) then

		write(*,*) 'Warning:  unknown JSON key with unexpected type'
		write(*,*) 'Key    :"'//key//'"'
		write(*,*) 'Type   :', var_type
		write(*,*)

	end if

	! always false, since we want to traverse all nodes:
	finished = .false.

end subroutine traverse_bomat_json

!===============================================================================

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

