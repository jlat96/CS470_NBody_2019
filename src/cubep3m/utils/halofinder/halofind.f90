! -------------------------------------------------------------------------------------------------------
! Post-processing halofinder utility written by JD Emberson September 2013.
! Determines mass of halos with overdensity halo_odc (parameter below) by
! searching over the local particle distribution around each fine mesh density
! peak. Use flag -DHVIR if aslo want to determine the mass of the halo evaulated
! at its virial radius (a function of redshift). 
! -------------------------------------------------------------------------------------------------------
 
program halofind 

    implicit none

    include 'mpif.h'
    include '../../parameters'

    !! Halofind parameters 
    real(4), parameter    :: halo_odc = 200.
    integer(4), parameter :: max_halo_np = nc**2 
    real(4), parameter    :: den_peak_cutoff = 100.
    integer(4), parameter :: min_halo_particles = 100
    integer(4), parameter :: ngrid_max = 420 
    integer(4), parameter :: nc_halo_max = 64
    integer(4), parameter :: max_maxima= 5*nc_halo_max**3
    integer(4), parameter :: nlist = 5*(nc_halo_max+1)**3
    logical, parameter    :: complete_shell = .true.
    logical, parameter    :: halo_write = .true.

    !! Physical constants
    real(4), parameter :: pi = 3.141592654
    real(4), parameter :: G = 1.0 / 6.0 / pi
    real(4), parameter :: threeover4pi = 3. / 4. / pi

    !! Parallelization variables
    integer(4), dimension(6) :: cart_neighbor
    integer(4), dimension(3) :: cart_coords
    integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

    !! Internal parametes 
    integer(4), parameter :: hc = nc / 2
    integer(4), parameter :: np = hc
    integer(4), parameter :: nc_node_dim = nc / nodes_dim
    integer(4), parameter :: np_node_dim = np / nodes_dim
    integer(4), parameter :: np_buffer = np_node_dim**3/2
    integer(4), parameter :: max_np = np_node_dim**3 + np_buffer
    integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
    integer(4), parameter :: tiles_node = tiles_node_dim * tiles_node_dim * tiles_node_dim
    integer(4), parameter :: nf_physical_tile_dim = nf_tile - 2 * nf_buf
    integer(4), parameter :: nf_physical_dim = nf_physical_tile_dim * tiles_node_dim * nodes_dim
    integer(4), parameter :: nf_physical_node_dim = nf_physical_tile_dim * tiles_node_dim
    integer(4), parameter :: mesh_scale = 4
    integer(4), parameter :: nc_buf = nf_buf / mesh_scale
    integer(4), parameter :: nc_tile_dim = ( nf_tile - 2 * nf_buf ) / mesh_scale
    integer(4), parameter :: hoc_nc_l = 1 - nc_buf
    integer(4), parameter :: hoc_nc_h = nc_tile_dim*tiles_node_dim + nc_buf
    integer(4), parameter :: hoc_pass_depth = 2*nc_buf
    real(4), parameter    :: rnf_buf = real(nf_buf)
    integer(4) :: np_local
    real(4) :: mass_p

    !! Checkpoint parameters
    character(len=*), parameter :: checkpoints = cubepm_root//'/input/halofinds'
    integer(4), parameter :: max_checkpoints = 100
    real(4), dimension(max_checkpoints) :: z_checkpoint
    integer(4) num_checkpoints, cur_checkpoint

    !! Particle arrays
    real(4), dimension(6, max_np) :: xvp
    integer(4), dimension(max_np) :: ll
    real(4), dimension(nf_tile+2, nf_tile, nf_tile) :: rho_f
    real(4), dimension(6, np_buffer) :: xp_buf
    real(4), dimension(6*np_buffer) :: send_buf, recv_buf
    integer(4) :: hoc(hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h, hoc_nc_l:hoc_nc_h) 

    !! Halofind arrays
    integer(4), dimension(max_maxima) :: isortpeak
    integer(4), dimension(max_halo_np) :: isortpos
    integer(4), dimension(nlist) :: isortdist
    real(4), dimension(3, max_maxima) :: ipeak
    integer(4), dimension(3, nlist) :: idist
    real(4), dimension(max_maxima) :: den_peak
    real(4), dimension(nlist) :: rdist
    integer(1) :: hpart_odc(max_np)
    real(4), dimension(max_maxima) :: halo_mesh_mass
    real(4), dimension(ngrid_max, ngrid_max, ngrid_max) :: finegrid
    integer(4) :: ilist_odc(max_halo_np)
    real(4), dimension(max_halo_np, 4) :: pos
#ifdef HVIR
    integer(1) :: hpart_vir(max_np)
    integer(4) :: ilist_vir(max_halo_np)
#endif

    !! Halofind parameters
    integer(4) :: irtot, nhalo
    real(4) :: z_write, a
    integer(4) :: pp, i, j, k, fstat
    integer(4), dimension(3) :: tile
    character (len=200) :: ofile
    character (len=7) :: z_s
    character (len=3) :: t_s
    character (len=5) :: r_s
    integer(4) :: num_candidates
    integer(4) :: jj, ii, iloc
    integer(8) :: search_fail, imass_odc, i_odc 
    real(4) :: mass_proxy, mass_odc, r_odc 
    real(4) :: v_disp
    real(4), dimension(3) :: x_mean, x2_mean, var_x, v_mean, v2_mean, l, offset, dx, l_CM
    real(4), dimension(3) :: r_wrt_halo, v_wrt_halo, v2_wrt_halo, hpos
    integer, parameter :: N_p = 50
    real(4), dimension(N_p) :: E
#ifdef PID_FLAG
    integer(8), dimension(N_p) :: pid_halo
    real(4), dimension(6,N_p) ::xv_halo
#endif
    real(4) :: dist, speed, E_tmp
    real(4), dimension(6) :: I_ij
    real(8) :: st1, st2
    integer(8) :: np_halo_local_odc, np_halo_odc, nhalo_tot, num_candidates_tot
#ifdef HVIR
    real(4) :: xflat, halo_vir, mass_vir, r_vir
    integer(8) :: imass_vir, i_vir
    integer(8) :: np_halo_local_vir, np_halo_vir
#endif

    real(8) :: wt1, wt2

    equivalence(isortpos, isortpeak)

#ifndef HVIR
    common xvp, xp_buf, send_buf, recv_buf, ll, hoc, rho_f, halo_mesh_mass, finegrid, ilist_odc, hpart_odc
#else
    common xvp, xp_buf, send_buf, recv_buf, ll, hoc, rho_f, halo_mesh_mass, finegrid, ilist_odc, hpart_odc, ilist_vir, hpart_vir
#endif

! -------------------------------------------------------------------------------------------------------
! MAIN
! -------------------------------------------------------------------------------------------------------

    call mpi_initialize
    call initialize_halofind
    call read_checkpoint_list

    do cur_checkpoint = 1, num_checkpoints

        st1 = mpi_wtime(ierr)

        call initvar

        a = 1. / (1. + z_checkpoint(cur_checkpoint))

        if (rank == 0) write(*,*) "Finding halos for z = ", z_checkpoint(cur_checkpoint)

        !
        ! Read particles and pass amongst nodes 
        !

        call read_particles
        call pass_particles

        !
        ! Construct linked list of particles and then add buffer particles from
        ! neighbouring nodes (so that we don't miss out on halos near edges). 
        !

        call link_list
        call buffer_particles

        !
        ! Find halo candidates based on local overdensities for each tile
        !

        wt1 = mpi_wtime(ierr)
        num_candidates = 0

        do i = 1, tiles_node
            tile(3) = (i-1) / (tiles_node_dim * tiles_node_dim)
            j = i - tile(3) * tiles_node_dim * tiles_node_dim
            tile(2) = (j-1) /  tiles_node_dim
            j = j - tile(2) * tiles_node_dim
            tile(1) = j - 1
            call find_halo_candidates(tile, num_candidates)
        enddo

        wt2 = mpi_wtime(ierr)
        if (rank == 0) write(*,*) "Finished find_halo_candidates ... elapsed time = ", wt2-wt1 
        
        !
        ! Sort density maxima 
        !

        isortpeak(:num_candidates) = (/ (i, i=1, num_candidates) /)
        call indexedsort(num_candidates, den_peak(:), isortpeak(:))
        ipeak(:, :num_candidates)  = ipeak(:, isortpeak(:num_candidates))
        halo_mesh_mass(:num_candidates) = halo_mesh_mass(isortpeak(:num_candidates))

        !
        ! Determine total number of halo candidates
        !

        call mpi_reduce(int(num_candidates, kind=8), num_candidates_tot, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
       
#ifdef HVIR
        !
        ! Determine Delta_vir using equation (6) of Bryan et al. (1997) for a flat universe
        !

        xflat    = (omega_m / a**3) / (omega_m / a**3 + omega_l) - 1.
        halo_vir = 18.*pi**2 + 82.*xflat - 39.*xflat**2
#endif

        !
        ! Open halo file
        !

        z_write = z_checkpoint(cur_checkpoint)
        call mpi_bcast(z_write, 1, mpi_real, 0, mpi_comm_world, ierr)
        write(z_s, "(f7.3)") z_write
        z_s = adjustl(z_s)

        write(r_s, "(i5)") rank
        r_s = adjustl(r_s)

        ofile = output_path//z_s(1:len_trim(z_s))//"halo"//r_s(1:len_trim(r_s))//".dat"

#ifdef BINARY
        open(unit=12, file=ofile, status="replace", iostat=fstat, form="binary")
#else
        open(unit=12, file=ofile, status="replace", iostat=fstat, form="unformatted")
#endif

        if (fstat /= 0) then
            write(*,*) "Error opening halo catalog for write"
            write(*,*) "rank", rank, "file:", ofile
            stop
        endif

        !
        ! Determine which candidates are to be considered halos and write their properties to file.
        !

        !! Halo ticker
        nhalo = 0

        !! Ticker recording how many particle searches end up empty handed
        search_fail = 0

        !! Node offsets
        offset(1) = cart_coords(3)*nf_physical_node_dim
        offset(2) = cart_coords(2)*nf_physical_node_dim
        offset(3) = cart_coords(1)*nf_physical_node_dim

        !! Initialize so that no particles are yet part of a halo
        hpart_odc = 0
#ifdef HVIR
        hpart_vir = 0
#endif

        !! Write file header (will rewrite nhalo later)
        write(12) nhalo
        write(12) halo_odc
#ifdef HVIR
        write(12) halo_vir
#endif       

        if (rank == 0) then

            write(*,*) "Searching over ", num_candidates_tot, " halo candidates..."
#ifndef HVIR
            write(*,*) "halo_odc = ", halo_odc
#else
            write(*,*) "halo_odc, halo_vir = ", halo_odc, halo_vir
#endif

        endif

        do iloc = num_candidates, 1, -1

            !! Fine mesh mass and position of this halo
            mass_proxy = halo_mesh_mass(iloc)
            hpos(:)    = ipeak(:, iloc)

            !! Search for particles by looking at local particle distribution
            call find_halo_particles(halo_odc, mass_proxy, hpos(:), r_odc, i_odc)
#ifdef HVIR
            call find_halo_particles(halo_vir, mass_proxy, hpos(:), r_vir, i_vir, DOVIR=1) 
#endif
    
            !! The following conditions must pass to be considered a halo
#ifndef HVIR
            if (i_odc >= min_halo_particles) then
#else
            if (i_odc >= min_halo_particles .and. i_vir >= min_halo_particles) then
#endif

                !! Halo ticker
                nhalo = nhalo + 1

                !! Initialize halo properties
                imass_odc = 0
                x_mean = 0.
                x2_mean = 0.
                v_mean = 0.
                v2_mean = 0.
                v_wrt_halo = 0.
                v2_wrt_halo = 0.
                l = 0.
                l_CM = 0.

                !! Loop over particles in this halo to get halo properties 
                do ii = 1, i_odc
                    pp = ilist_odc(ii)
                    imass_odc = imass_odc + 1
                    x_mean = x_mean + xvp(:3, pp)
                    x2_mean = x2_mean + xvp(:3, pp)**2
                    v_mean = v_mean + xvp(4:, pp)
                    v2_mean = v2_mean + xvp(4:, pp)**2
                    dx = hpos(:) - xvp(:3, pp)
                    l(1) = l(1) + (dx(3)*xvp(5,pp) - dx(2)*xvp(6,pp))
                    l(2) = l(2) + (dx(1)*xvp(6,pp) - dx(3)*xvp(4,pp))
                    l(3) = l(3) + (dx(2)*xvp(4,pp) - dx(1)*xvp(5,pp))

                    !! Remove this particle from future halo candidates
                    hpart_odc(pp) = 1
                enddo

#ifdef HVIR
                imass_vir = 0
                do ii = 1, i_vir
                    pp = ilist_vir(ii)
                    imass_vir = imass_vir + 1

                    !! Remove this particle from future halo candidates
                    hpart_vir(pp) = 1 
                enddo
                mass_vir = mass_p * imass_vir
#endif

                mass_odc = mass_p * imass_odc
                hpos(:) = hpos(:) + offset
                x_mean = x_mean/real(imass_odc) + offset
                x2_mean = x2_mean/real(imass_odc)
                v_mean = v_mean/real(imass_odc)
                v2_mean = v2_mean/real(imass_odc)
                l = l/real(imass_odc)
                l_CM(1) = l(1) - (x_mean(3)*v_mean(2) - x_mean(2)*v_mean(3))
                l_CM(2) = l(2) - (x_mean(1)*v_mean(3) - x_mean(3)*v_mean(1))
                l_CM(3) = l(3) - (x_mean(2)*v_mean(1) - x_mean(1)*v_mean(2))
                v_disp = sqrt(v2_mean(1) + v2_mean(2) + v2_mean(3))
                var_x = real(imass_odc)/(real(imass_odc-1)) * (x2_mean - (x_mean-offset)**2)

#ifdef PID_FLAG
                pid_halo = 0
#endif
                E = 0.
                E_tmp = 0.
                I_ij = 0.

                ii = 1
                do ii = 1, i_odc
                    pp = ilist_odc(ii)

                    r_wrt_halo = xvp(:3,pp) - (x_mean - offset)
                    v_wrt_halo = xvp(4:,pp) - v_mean
                    v2_wrt_halo = v2_wrt_halo + v_wrt_halo(:)**2
                    dist = sqrt(r_wrt_halo(1)**2 +r_wrt_halo(2)**2 + r_wrt_halo(3)**2)
                    I_ij(1) = I_ij(1)+r_wrt_halo(1)*r_wrt_halo(1)
                    I_ij(2) = I_ij(2)+r_wrt_halo(1)*r_wrt_halo(2)
                    I_ij(3) = I_ij(3)+r_wrt_halo(1)*r_wrt_halo(3)
                    I_ij(4) = I_ij(4)+r_wrt_halo(2)*r_wrt_halo(2)
                    I_ij(5) = I_ij(5)+r_wrt_halo(2)*r_wrt_halo(3)
                    I_ij(6) = I_ij(6)+r_wrt_halo(3)*r_wrt_halo(3)

                    speed = sqrt(v_wrt_halo(1)**2 +v_wrt_halo(2)**2 + v_wrt_halo(3)**2)
                    E_tmp = 0.5*(speed)**2 - imass_odc*mass_p*G/dist

                    !! Find the most bound particles
                    do jj = 1, N_p
                        if (E_tmp < E(jj)) then
                            E(jj+1:N_p) = E(jj:N_p-1)
                            E(jj) = E_tmp
#ifdef PID_FLAG
                            pid_halo(jj+1:N_p) = pid_halo(jj:N_p-1)
                            pid_halo(jj) = PID(pp)
                            xv_halo(:,jj) = xvp(:,pp)
#endif
                            exit
                        endif
                    enddo
                enddo

                v2_wrt_halo(:) = v2_wrt_halo(:)/real(imass_odc)

#ifdef PID_FLAG
#ifdef HVIR
                if (halo_write) write(12) hpos(:), mass_odc, mass_vir, r_odc, r_vir, x_mean, v_mean, l_CM, v2_wrt_halo, var_x, pid_halo, xv_halo
#else
                if (halo_write) write(12) hpos(:), mass_odc, r_odc, x_mean, v_mean, l_CM, v2_wrt_halo, var_x, pid_halo, xv_halo
#endif
#else
#ifdef HVIR
                if (halo_write) write(12) hpos(:), mass_odc, mass_vir, r_odc, r_vir, x_mean, v_mean, l_CM, v2_wrt_halo, var_x, I_ij
#else
                if (halo_write) write(12) hpos(:), mass_odc, r_odc, x_mean, v_mean, l_CM, v2_wrt_halo, var_x, I_ij
#endif
#endif

            endif !! i_odc test 

        enddo !! iloc loop 

        !
        ! Rewrite nhalo in the header
        !
        
        close(12)

#ifdef BINARY
        open (unit=12, file=ofile, status="old", iostat=fstat, form="binary", access="direct", recl=4)
        write(12, rec=1) nhalo
#else
        open (unit=12, file=ofile, status="old", iostat=fstat, form="unformatted", access="direct", recl=4)
        write(12, rec=2) nhalo
#endif
        close(12)

        !! Stop timer
        st2 = mpi_wtime(ierr)

        call mpi_barrier(mpi_comm_world, ierr)
        
        !
        ! Count particles that were found within halos
        !

        np_halo_local_odc = 0
        do ii = 1, np_local
            if (hpart_odc(ii) == 1) np_halo_local_odc = np_halo_local_odc + 1
        enddo
        call mpi_reduce(np_halo_local_odc, np_halo_odc, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

#ifdef HVIR
        np_halo_local_vir = 0
        do ii = 1, np_local
            if (hpart_vir(ii) == 1) np_halo_local_vir = np_halo_local_vir + 1
        enddo
        call mpi_reduce(np_halo_local_vir, np_halo_vir, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)
#endif

        call mpi_reduce(int(nhalo, kind=8), nhalo_tot, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

        !
        ! Print some info to screen
        !

        if (rank == 0) then

            write(*,*) "Done ... found ", nhalo_tot, " halos"
#ifndef HVIR
            write(*,*) "Total Particles found in halos = ", np_halo_odc
#else
            write(*,*) "Total Particles found in halos = ", np_halo_odc, np_halo_vir
#endif
            write(*,*) "ELAPSED TIME: ", st2-st1            

            if (search_fail > 0) then
                write(*,*) "CONSIDER INCREASING search_ratio ... search_fail = ", search_fail
            endif

        endif

        call mpi_barrier(mpi_comm_world, ierr)

    enddo !! cur_checkpoint

    call mpi_finalize(ierr)

! -------------------------------------------------------------------------------------------------------
! SUBROUTINES
! -------------------------------------------------------------------------------------------------------

contains

! -------------------------------------------------------------------------------------------------------

subroutine mpi_initialize

    implicit none

    integer(4) :: i, j, nodes_returned
    integer(4) :: dims(3), ndim
    logical :: periodic(3), reorder

    !! Set up global mpi communicator

    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)

    call mpi_comm_size(mpi_comm_world, nodes_returned, ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)

    if (nodes_returned /= nodes ) then
      write(*,*) 'cic_pow compiled for a different number of nodes'
      write(*,*) 'mpirun nodes = ', nodes_returned, ' cic_pow nodes = ',nodes
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    if (mod(nc, nodes) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into slabs'
      write(*,*) 'nc = ', nc, ' nodes = ', nodes, ' mod(nc, nodes) != 0'
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ierr)

    if (rank == 0) then
      write(*,*) 'cic_pow running on ', nodes, ' nodes'
      write(*,*) 'using cubic distribution: ', nodes_dim, ' nodes per dimension'
      write(*,*) nc, ' cells in mesh'
    endif

    !! Create cartesian communicator based on cubic decomposition

    dims(:) = nodes_dim
    periodic(:) = .true.
    reorder = .false.
    ndim = 3

    call mpi_cart_create(mpi_comm_world, ndim, dims, periodic, &
                       reorder, mpi_comm_cart, ierr)
    call mpi_comm_rank(mpi_comm_cart, cart_rank, ierr)
    call mpi_cart_coords(mpi_comm_cart, cart_rank, ndim,  &
                         cart_coords, ierr)

    ! cart_neighbor(1) -> down (negative z)
    ! cart_neighbor(2) -> up (positive z)
    ! cart_neighbor(3) -> back (negative y)
    ! cart_neighbor(4) -> front (positive y)
    ! cart_neighbor(5) -> left (negative x)
    ! cart_neighbor(6) -> right (positive x)

    do i = 0, ndim-1
      call mpi_cart_shift(mpi_comm_cart, i, 1, cart_neighbor(2*(i+1)-1), &
                          cart_neighbor(2*(i+1)), ierr)
    enddo

end subroutine mpi_initialize

! -------------------------------------------------------------------------------------------------------

subroutine initialize_halofind
    !
    ! Initializes halo cell search arrays 
    !    

    implicit none

    integer(4) :: ii, i, j, k
    real(4) :: r
    integer(4), allocatable , dimension(:,:) :: idist_tmp
    allocate(idist_tmp(3,nlist))

    ! Loop through a box of length 2*nc_halo_max
    ! if cell is within sphere of radius = box length / 2
    ! include distince in rdist at entry ii
    ! ordered bottom left to top right

    ii = 0
    do i = -nc_halo_max, nc_halo_max
        do j = -nc_halo_max, nc_halo_max
            do k = -nc_halo_max, nc_halo_max
                r = sqrt(real(i)**2 + real(j)**2 + real(k)**2)
                if (r > nc_halo_max) cycle
                ii = ii + 1
                if (ii > nlist) then
                    write(*,*) 'ii exceeded ', nlist
                    pause
                endif
                idist(:, ii) = (/i, j, k/)
                rdist(ii)=r
            enddo
        enddo
    enddo
    irtot = ii

    ! sorts the rdist array from lowest to highest radial position
    ! from center of sphere saves rdist array position values in idist

    isortdist(:ii) = (/ (i, i = 1, ii) /)
    call indexedsort(ii, rdist, isortdist)
    idist_tmp(:,:ii) = idist(:, isortdist(:ii))
    idist(:,:ii) = idist_tmp(:, :ii)
    deallocate(idist_tmp)

end subroutine initialize_halofind

! -------------------------------------------------------------------------------------------------------

subroutine read_checkpoint_list
    !
    ! Read in the list of redshift checkpoints for which to calculate spectra
    ! for 
    !

    implicit none

    integer :: i, fstat

    if (rank == 0) then

        open(11, file=checkpoints, status='old', iostat=fstat)

        !! Check for opening error 
        if (fstat /= 0) then
            print *,'ERROR: Cannot open checkpoint list file'
            print *,'rank ', rank, ' file: ', checkpoints
            call mpi_abort(mpi_comm_world, ierr, ierr)
        endif

        !! Read the redshifts
        do num_checkpoints = 1, max_checkpoints
            read(unit=11, err=51, end=41, fmt='(f20.10)') z_checkpoint(num_checkpoints)
        enddo

        !! Tabulate total number of checkpoints
   41  num_checkpoints = num_checkpoints - 1
   51  close(11)

        !! Print to screen
        print *, 'Checkpoints to recompose:'
        do i = 1, num_checkpoints
            write(*, '(f5.1)') z_checkpoint(i)
        enddo

    endif

    call mpi_bcast(num_checkpoints, 1, mpi_integer, 0, mpi_comm_world, ierr)

end subroutine read_checkpoint_list

! -------------------------------------------------------------------------------------------------------

subroutine initvar
    !
    ! Initialize data arrays
    !

    implicit none

    integer :: k

    do k = 1, max_np
       xvp(:, k) = 0.
       ll(k) = 0
    enddo

    do k = 1, np_buffer
        xp_buf(:, k) = 0.
    enddo

    do k = 1, 6*np_buffer
        send_buf(k) = 0.
        recv_buf(k) = 0.
    enddo

    hoc(:, :, :) = 0
    finegrid(:, :, :) = 0.

    return

end subroutine initvar

! -------------------------------------------------------------------------------------------------------

subroutine read_particles
    !
    ! Read x, y, z positions and velocities and store in xvp
    !

    implicit none

    real z_write, np_total
    integer j, fstat
    character(len=7) :: z_string
    character(len=4) :: rank_string
    character(len=100) :: check_name

    real(8) :: sec1, sec2

    !! These are unnecessary headers from the checkpoint
    real(4) :: a, t, tau, dt_f_acc, dt_c_acc, dt_pp_acc
    integer(4) :: nts, sim_checkpoint, sim_projection, sim_halofind

    sec1 = mpi_wtime(ierr)

    !! Generate checkpoint names on each node
    if (rank==0) then
      z_write = z_checkpoint(cur_checkpoint)
    endif

    call mpi_bcast(z_write, 1, mpi_real, 0, mpi_comm_world, ierr)

    !! Determine the file name
    write(z_string,'(f7.3)') z_write
    z_string=adjustl(z_string)

    write(rank_string,'(i4)') rank
    rank_string=adjustl(rank_string)

    check_name = output_path//z_string(1:len_trim(z_string))//'xv'// &
               rank_string(1:len_trim(rank_string))//'.dat'

    !! Open the file    
#ifdef BINARY
    open(unit=21, file=check_name, status='old', iostat=fstat, form='binary')
#else
    open(unit=21, file=check_name, status='old', iostat=fstat, form='unformatted')
#endif

    !! Check for opening error
    if (fstat /= 0) then
      write(*,*) 'ERROR: Cannot open checkpoint position file'
      write(*,*) 'rank', rank, ' file: ',check_name
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    !! Read in checkpoint header data
#ifdef PPINT
    read(21) np_local, a, t, tau, nts, dt_f_acc, dt_pp_acc, dt_c_acc, sim_checkpoint, &
               sim_projection, sim_halofind, mass_p
#else
    read(21) np_local, a, t, tau, nts, dt_f_acc, dt_c_acc, sim_checkpoint, &
               sim_projection, sim_halofind, mass_p
#endif

    !! Check for memory problems
    if (np_local > max_np) then
      write(*,*) 'ERROR: Too many particles to store in memory!'
      write(*,*) 'rank', rank, 'np_local', np_local, 'max_np', max_np
      call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    !! Tally up total number of particles
    call mpi_reduce(real(np_local, kind=4), np_total, 1, mpi_real, &
                         mpi_sum, 0, mpi_comm_world, ierr)

    if (rank == 0) write(*,*) 'Total number of particles = ', int(np_total, 8)

    do i = 1, np_local
        read(21) xvp(:, i)
    enddo

    close(21)

#ifdef KAISER

    !! Red Shift Distortion: x_z -> x_z +  v_z/H(Z)   
    !! Converting seconds into simulation time units
    !! cancels the H0...

    !! xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i)) +
    !xv(6,ip+1:ip+nploc(i))*1.5*sqrt(omegam/cubepm_a)
    !! xv(3,ip+1:ip+nploc(i))=xv(3,ip+1:ip+nploc(i)) +
    !xv(6,ip+1:ip+nploc(i))*1.5/sqrt(cubepm_a*(1+cubepm_a*omegak/omegam +
    !omegav/omegam*cubepm_a**3))

    xvp(3, :) = xvp(3, :) + xvp(6, :) * 1.5 / sqrt(a * (1 + a * (1 - omega_m - omega_l) / omega_m + omega_l / omega_m * a**3))

    call pass_particles

    if(rank == 0) then
       write(*,*) '**********************'
       write(*,*) 'Included Kaiser Effect'
       write(*,*) 'Omega_m = ', omega_m, ' a = ', a
       !write(*,*) '1/H(z) =', 1.5*sqrt(omegam/cubepm_a)
       write(*,*) '1 / H(z) = ', 1.5 / sqrt(a * (1 + a * (1 - omega_m - omega_l) / omega_m + omega_l / omega_m * a**3))
       write(*,*) '**********************'
    endif
#endif

    sec2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "Finished read_particles ... elapsed time = ", sec2 - sec1

end subroutine read_particles

! -------------------------------------------------------------------------------------------------------

subroutine pass_particles
    !
    ! Pass particles inside buffer space to their appropriate nodes.
    !

    implicit none

    integer i,pp,np_buf,np_exit,npo,npi
    real x(6),lb,ub,np_finalr
    integer(8) :: np_final
    integer, dimension(mpi_status_size) :: status,sstatus,rstatus
    integer :: tag,srequest,rrequest,sierr,rierr
    real(4), parameter :: eps = 1.0e-03 
    real(8) :: sec1, sec2

    sec1 = mpi_wtime(ierr)
 
    !
    ! Identify particles within the buffer
    !

    lb = 0.
    ub = real(nc_node_dim)

    np_buf = 0
    pp = 1

    do

        if (pp > np_local) exit

        !! Read its position  
        x = xvp(:, pp)

        !! See if it lies within the buffer
        if (x(1) < lb .or. x(1) >= ub .or. x(2) < lb .or. x(2) >= ub .or. &
            x(3) < lb .or. x(3) >= ub ) then

            !write (*,*) 'PARTICLE OUT', xvp(:, pp)

            !! Make sure we aren't exceeding the maximum
            np_buf = np_buf + 1

            if (np_buf > np_buffer) then
                print *, rank, 'np_buffer =', np_buffer, 'exceeded - np_buf =', np_buf
                call mpi_abort(mpi_comm_world, ierr, ierr)
            endif

            xp_buf(:, np_buf) = xvp(:, pp)
            xvp(:, pp)        = xvp(:, np_local)
            np_local          = np_local - 1

            cycle

        endif

        pp = pp + 1

    enddo

    call mpi_reduce(np_buf, np_exit, 1, mpi_integer, mpi_sum, 0, &
                    mpi_comm_world, ierr)

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'np_exit = ', np_buf
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    if (rank == 0) print *,'Total exiting particles = ',np_exit

    !
    ! Pass +x
    !

    !! Find particles that need to be passed

    tag = 11
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(1, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'np_out=', npo
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(6), &
                              tag, cart_neighbor(5), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf + pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(1, np_buf + pp) = max(xp_buf(1, np_buf+pp) - ub, lb)
    enddo

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'x+ np_local=', np_local
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    pp = 1

    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf + pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi - 1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'x+ np_exit=', np_buf, np_local
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif 

    !
    ! Pass -x
    !

    tag = 12
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(1, pp) < lb) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(5), &
                              tag, cart_neighbor(6), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(1, np_buf+pp) = min(xp_buf(1,np_buf+pp) + ub, ub-eps)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi - 1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

    !
    ! Pass +y
    !

    tag = 13
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(2, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(4), &
                              tag, cart_neighbor(3), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(2, np_buf+pp) = max(xp_buf(2, np_buf+pp)-ub, lb)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi-1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

    !
    ! Pass -y
    !

    tag = 14
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(2,pp) < lb) then
        npo = npo+1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(3), &
                              tag, cart_neighbor(4), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(2, np_buf+pp) = min(xp_buf(2, np_buf+pp)+ub, ub-eps)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local+1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

    !
    ! Pass +z
    !

    tag = 15
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(3, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(2), &
                              tag,cart_neighbor(1),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*6+1:pp*6)
      xp_buf(3,np_buf+pp)=max(xp_buf(3,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp=pp+1
    enddo

    np_buf=np_buf+npi

    !
    ! Pass -z
    !

    tag=16
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(3,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*6+1:npo*6)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(1), &
                              tag,cart_neighbor(2),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*6+1:pp*6)
      xp_buf(3,np_buf+pp)=min(xp_buf(3,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp=pp+1
    enddo

    np_buf=np_buf+npi

#ifdef DEBUG
    do i=0,nodes-1
      if (rank==i) print *,rank,'particles left in buffer=',np_buf
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif 

    call mpi_reduce(np_buf,np_exit,1,mpi_integer,mpi_sum,0, &
                    mpi_comm_world,ierr)

    if (rank == 0) print *,'total buffered particles =',np_exit

    !! Tally up total number of particles
    call mpi_reduce(real(np_local, kind=4), np_finalr, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
    np_final = int(np_finalr, kind=8)

    if (rank == 0) then
      print *, "np_final = ", np_final
      if (np_final /= np**3) then
        print *,'ERROR: total number of particles incorrect after passing'
      endif
    endif

!!  Check for particles out of bounds

    do i=1,np_local
      if (xvp(1,i) < 0 .or. xvp(1,i) >= nc_node_dim .or. &
          xvp(2,i) < 0 .or. xvp(2,i) >= nc_node_dim .or. &
          xvp(3,i) < 0 .or. xvp(3,i) >= nc_node_dim) then
        print *,'particle out of bounds',rank,i,xvp(:3,i),nc_node_dim
      endif
    enddo

    sec2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "Finished pass_particles ... elapsed time = ", sec2 - sec1

end subroutine pass_particles

! -------------------------------------------------------------------------------------------------------

subroutine buffer_particles
    !
    ! Exchange coarse mesh buffers for finding halos near node edges. Add these
    ! particles to the linked list.
    !

    implicit none

    integer :: pp, i, j, k
    integer :: np_buf, nppx, npmx, nppy, npmy, nppz, npmz 
    integer, dimension(mpi_status_size) :: status, sstatus, rstatus
    integer :: tag, srequest, rrequest, sierr, rierr
    real(4), parameter :: eps = 1.0e-03
    integer(4) :: np_local_start, np_buffer_sent_local
    integer(8) :: np_buffer_sent
    real(8) :: sec1, sec2

    sec1 = mpi_wtime(ierr)

    np_local_start = np_local

    !
    ! Pass +x
    ! 

    tag = 11
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
            do i = hoc_nc_h-hoc_pass_depth, hoc_nc_h
                pp = hoc(i, j, k)    
                do
                    if (pp == 0) exit
                    if (xvp(1, pp) >= nf_physical_node_dim-rnf_buf) then
                        np_buf = np_buf + 1
                        send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp) 
                    endif
                    pp = ll(pp)
                enddo
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    nppx = np_buf

    call mpi_sendrecv_replace(nppx, 1, mpi_integer, cart_neighbor(6), tag, &
                              cart_neighbor(5), tag, mpi_comm_world, status, ierr)

    if (np_local+nppx > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer particles (nppx): ", rank, np_local, nppx, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, nppx*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do i = 1, nppx
        xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
        xvp(1, np_local+i) = max(xvp(1,np_local+i)-nf_physical_node_dim, -rnf_buf)
    enddo

    np_local = np_local + nppx

    !
    ! Pass -x
    !

    tag = 12
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_h
            do i = hoc_nc_l, hoc_nc_l+hoc_pass_depth
                pp = hoc(i, j, k)
                do
                    if (pp == 0) exit
                    if (xvp(1, pp) < rnf_buf) then
                        np_buf = np_buf + 1
                        send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                    endif
                    pp = ll(pp)
                enddo
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    npmx = np_buf

    call mpi_sendrecv_replace(npmx, 1, mpi_integer, cart_neighbor(5), tag, &
                              cart_neighbor(6), tag, mpi_comm_world, status, ierr)

    if (np_local+npmx > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer particles (npmx): ", rank, np_local, npmx, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npmx*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do i = 1, npmx
        xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
        if (abs(xvp(1, np_local+i)) .lt. eps) then
            if(xvp(1, np_local+i) < 0.) then
                xvp(1, np_local+i) = -eps
            else
                xvp(1, np_local+i) = eps
            endif
        endif
        xvp(1, np_local+i) = min(xvp(1, np_local+i)+real(nf_physical_node_dim,4), &
                                 nf_physical_node_dim+rnf_buf-eps)
    enddo

    np_local = np_local + npmx

    !
    ! Add additional particles to the linked list
    !

    pp = np_local-npmx-nppx + 1
    do
        if (pp > np_local) exit
        i = floor(xvp(1, pp)/mesh_scale) + 1
        j = floor(xvp(2, pp)/mesh_scale) + 1
        k = floor(xvp(3, pp)/mesh_scale) + 1
#ifdef DIAG
        if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
            j < hoc_nc_l .or. j > hoc_nc_h .or. &
            k < hoc_nc_l .or. k > hoc_nc_h) then
            write (*, *) 'BUFFER PARTICLE DELETED', xvp(:, pp)      
            xvp(:,pp) = xvp(:,np_local)
            np_local = np_local - 1
            cycle
        endif
#endif
        ll(pp) = hoc(i, j, k)
        hoc(i, j, k) = pp
        pp = pp + 1
    enddo
   
    !
    ! Pass +y
    ! 

    tag = 13
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_h-hoc_pass_depth, hoc_nc_h 
            do i = hoc_nc_l, hoc_nc_h 
                pp = hoc(i, j, k)
                do
                    if (pp == 0) exit
                    if (xvp(2, pp) >= nf_physical_node_dim-rnf_buf) then
                        np_buf = np_buf + 1
                        send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                    endif
                    pp = ll(pp)
                enddo
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    nppy = np_buf

    call mpi_sendrecv_replace(nppy, 1, mpi_integer, cart_neighbor(4), tag, &
                              cart_neighbor(3), tag, mpi_comm_world, status, ierr)

    if (np_local+nppy > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer particles (nppy): ", rank, np_local, nppy, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, nppy*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do i = 1, nppy
        xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
        xvp(2, np_local+i) = max(xvp(2,np_local+i)-nf_physical_node_dim, -rnf_buf)
    enddo

    np_local = np_local + nppy

    !
    ! Pass -y
    !

    tag = 14
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_h
        do j = hoc_nc_l, hoc_nc_l+hoc_pass_depth 
            do i = hoc_nc_l, hoc_nc_h 
                pp = hoc(i, j, k)
                do
                    if (pp == 0) exit
                    if (xvp(2, pp) < rnf_buf) then
                        np_buf = np_buf + 1
                        send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                    endif
                    pp = ll(pp)
                enddo
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    npmy = np_buf

    call mpi_sendrecv_replace(npmy, 1, mpi_integer, cart_neighbor(3), tag, &
                              cart_neighbor(4), tag, mpi_comm_world, status, ierr)

    if (np_local+npmy > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer particles (npmy): ", rank, np_local, npmy, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npmy*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do i = 1, npmy
        xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
        if (abs(xvp(2, np_local+i)) .lt. eps) then
            if(xvp(2, np_local+i) < 0.) then
                xvp(2, np_local+i) = -eps
            else
                xvp(2, np_local+i) = eps
            endif
        endif
        xvp(2, np_local+i) = min(xvp(2,np_local+i)+real(nf_physical_node_dim,4), &
                                 nf_physical_node_dim+rnf_buf-eps)
    enddo

    np_local = np_local + npmy

    !
    ! Add additional particles to the linked list 
    !

    pp = np_local-npmy-nppy + 1
    do
        if (pp > np_local) exit
        i = floor(xvp(1, pp)/mesh_scale) + 1
        j = floor(xvp(2, pp)/mesh_scale) + 1
        k = floor(xvp(3, pp)/mesh_scale) + 1
#ifdef DIAG
        if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
            j < hoc_nc_l .or. j > hoc_nc_h .or. &
            k < hoc_nc_l .or. k > hoc_nc_h) then
            write (*, *) 'BUFFER PARTICLE DELETED', xvp(:, pp)
            xvp(:,pp) = xvp(:,np_local)
            np_local = np_local - 1
            cycle
        endif
#endif
        ll(pp) = hoc(i, j, k)
        hoc(i, j, k) = pp
        pp = pp + 1
    enddo

    !
    ! Pass +z
    ! 

    tag = 15
    np_buf = 0

    do k = hoc_nc_h-hoc_pass_depth, hoc_nc_h 
        do j = hoc_nc_l, hoc_nc_h
            do i = hoc_nc_l, hoc_nc_h 
                pp = hoc(i, j, k)
                do
                    if (pp == 0) exit
                    if (xvp(3, pp) >= nf_physical_node_dim-rnf_buf) then
                        np_buf = np_buf + 1
                        send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                    endif
                    pp = ll(pp)
                enddo
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    nppz = np_buf

    call mpi_sendrecv_replace(nppz, 1, mpi_integer, cart_neighbor(2), tag, &
                              cart_neighbor(1), tag, mpi_comm_world, status, ierr)

    if (np_local+nppz > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer particles (nppz): ", rank, np_local, nppz, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, nppz*6, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do i = 1, nppz
        xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
        xvp(3, np_local+i) = max(xvp(3,np_local+i)-nf_physical_node_dim, -rnf_buf)
    enddo

    np_local = np_local + nppz

    !
    ! Pass -z
    !

    tag = 16
    np_buf = 0

    do k = hoc_nc_l, hoc_nc_l+hoc_pass_depth 
        do j = hoc_nc_l, hoc_nc_h 
            do i = hoc_nc_l, hoc_nc_h
                pp = hoc(i, j, k)
                do
                    if (pp == 0) exit
                    if (xvp(3, pp) < rnf_buf) then
                        np_buf = np_buf + 1
                        send_buf((np_buf-1)*6+1:np_buf*6) = xvp(:, pp)
                    endif
                    pp = ll(pp)
                enddo
            enddo
        enddo
    enddo

    if (6*np_buf > np_buffer) then
        write(*,*) "ERROR: Not enough space to send buffer particles: ", rank, np_buf, np_buffer
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    npmz = np_buf

    call mpi_sendrecv_replace(npmz, 1, mpi_integer, cart_neighbor(1), tag, &
                              cart_neighbor(2), tag, mpi_comm_world, status, ierr)

    if (np_local+npmz > max_np) then
        write(*,*) "ERROR: Not enough space to receive buffer particles (npmz): ", rank, np_local, npmz, max_np
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    call mpi_isend(send_buf, np_buf*6, mpi_real, cart_neighbor(1), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npmz*6, mpi_real, cart_neighbor(2), &
                   tag, mpi_comm_world, rrequest, rierr)

    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do i = 1, npmz
        xvp(:, np_local+i) = recv_buf((i-1)*6+1:(i-1)*6+6)
        if (abs(xvp(3, np_local+i)) .lt. eps) then
            if(xvp(3, np_local+i) < 0.) then
                xvp(3, np_local+i) = -eps
            else
                xvp(3, np_local+i) = eps
            endif
        endif
        xvp(3,np_local+i) = min(xvp(3,np_local+i)+real(nf_physical_node_dim,4), &
                                 nf_physical_node_dim+rnf_buf-eps)
    enddo

    np_local = np_local + npmz

    !
    ! Add additional particles to the linked list 
    !

    pp = np_local-npmz-nppz + 1
    do
        if (pp > np_local) exit
        i = floor(xvp(1, pp)/mesh_scale) + 1
        j = floor(xvp(2, pp)/mesh_scale) + 1
        k = floor(xvp(3, pp)/mesh_scale) + 1
#ifdef DIAG
        if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
            j < hoc_nc_l .or. j > hoc_nc_h .or. &
            k < hoc_nc_l .or. k > hoc_nc_h) then
            write (*, *) 'BUFFER PARTICLE DELETED', xvp(:, pp)
            xvp(:,pp) = xvp(:,np_local)
            np_local = np_local - 1
            cycle
        endif
#endif
        ll(pp) = hoc(i, j, k)
        hoc(i, j, k) = pp
        pp = pp + 1
    enddo

    !
    ! Collect some statistics
    !

    np_buffer_sent_local = np_local - np_local_start
    call mpi_reduce(int(np_buffer_sent_local, kind=8), np_buffer_sent, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

    if (rank == 0) write(*,*) "Total Buffer Particles in Halo Search: ", np_buffer_sent

    sec2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "Finished buffer_particles ... elapsed time = ", sec2 - sec1

end subroutine buffer_particles

! -------------------------------------------------------------------------------------------------------

subroutine link_list
    !
    ! Makes linked list of particles in each coarse mesh cell.
    !

    implicit none

    integer(4) :: i, j, k, pp
    real(8) :: sec1, sec2

    sec1 = mpi_wtime(ierr)

    pp = 1
    do
        if (pp > np_local) exit
        i = floor(xvp(1, pp)/mesh_scale) + 1
        j = floor(xvp(2, pp)/mesh_scale) + 1
        k = floor(xvp(3, pp)/mesh_scale) + 1
        if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
            j < hoc_nc_l .or. j > hoc_nc_h .or. &
            k < hoc_nc_l .or. k > hoc_nc_h) then
            write (*, *) 'PARTICLE DELETED', xvp(:,pp)
            xvp(:,pp) = xvp(:,np_local)
            np_local = np_local - 1
            cycle
        else
            ll(pp) = hoc(i, j, k)
            hoc(i, j, k) = pp
        endif
        pp = pp + 1
    enddo

    sec2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "Finished link_list ... elapsed time = ", sec2 - sec1

end subroutine link_list

! -------------------------------------------------------------------------------------------------------

subroutine find_halo_particles(HODC, HMASS, HPOS, RODC, ITOT, DOVIR)
    !
    ! Finds refined density maxima and then sorts particles based on their
    ! distance from here, moving one particle at a time until the target radius
    ! and overdensity are reached.
    !

    implicit none

    real(4), intent(in) :: HODC, HMASS
    real(4), dimension(3), intent(inout) :: HPOS
    real(4), intent(out) :: RODC
    integer(8), intent(out) :: ITOT
    integer(4), intent(in), optional :: DOVIR
    logical :: HALOVIR

    real :: r, dgrid, maxfinegrid, rrefine, rsearch
    integer :: pp, np_search, ii, jj, kk, i, j, k
    integer :: crbox(3, 2), csbox(3, 2), frbox(3, 2)
    integer :: ngrid(3)
    real :: dr(3), p(3)

    real :: odci, odcj, r1, r2, d1, d2, w1, w2

    integer, parameter :: search_ratio = 2 
    integer, parameter :: refine_ratio = 5

    !
    ! Determine which overdensity we are trying to reach (only matters ifdef HVIR)
    !

    if (present(DOVIR)) then

        HALOVIR = .true.

    else

        HALOVIR = .false.

    endif

    !
    ! Initialize parameters
    !

    !! Use the mass proxy to guess at the size of the halo
    !! This will be the radius within which the refined density peak will be found 
    rrefine = (3. * HMASS / 4. / pi / HODC)**(1./3.)

    !! This (larger) radius will be used to store particle positions to determine which are part of the halo
    rsearch = search_ratio * rrefine

    !! Coarse mesh cells within the refined region 
    crbox(:, 1) = int((HPOS(:)-rrefine)/mesh_scale) + 1
    crbox(:, 2) = int((HPOS(:)+rrefine)/mesh_scale) + 1

    !! Coarse mesh cells within the search region
    csbox(:, 1) = int((HPOS(:)-rsearch)/mesh_scale) + 1
    csbox(:, 2) = int((HPOS(:)+rsearch)/mesh_scale) + 1

    !! Boundary of the refined region in fine mesh cell units
    frbox(:, 1) = mesh_scale*(crbox(:, 1) - 1)
    frbox(:, 2) = mesh_scale*crbox(:, 2)

    !! Number of extra refined cells in this region and their spacing
    ngrid(:) = refine_ratio*(frbox(:, 2) - frbox(:, 1))
    dgrid    = 1./refine_ratio

    if (ngrid(1) > ngrid_max .or. ngrid(2) > ngrid_max .or. ngrid(3) > ngrid_max) write(*,*) "ERROR: ngrid = ", ngrid

    !! Exclude coarse mesh cells in the search region that lie within the tile buffer
    csbox(1, 1) = max(csbox(1, 1), hoc_nc_l)
    csbox(1, 2) = min(csbox(1, 2), hoc_nc_h)
    csbox(2, 1) = max(csbox(2, 1), hoc_nc_l)
    csbox(2, 2) = min(csbox(2, 2), hoc_nc_h)
    csbox(3, 1) = max(csbox(3, 1), hoc_nc_l)
    csbox(3, 2) = min(csbox(3, 2), hoc_nc_h)

    !
    ! Store particle positions within the search radius at the same time as NGP
    ! interpolating particles within the refined radius in order to find the density peak.
    !

!    finegrid(:, :, :) = 0.
    np_search = 0

#ifdef HVIR
    if (.not. HALOVIR) then
#endif

        do k = csbox(3,1), csbox(3,2)
            do j = csbox(2,1), csbox(2,2)
                do i = csbox(1,1), csbox(1,2)
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (hpart_odc(pp) == 0) then !! particle is not yet part of a halo
                            p(:) = xvp(:3, pp)
                            dr   = HPOS(:) - p(:)
                            r    = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                            if (r < rsearch) then 
                                np_search = np_search + 1
                                pos(np_search, 1:3) = p(:)
                                ilist_odc(np_search) = pp
                                if (r < rrefine) then
                                    ii = int((p(1)-frbox(1,1))/dgrid) + 1
                                    jj = int((p(2)-frbox(2,1))/dgrid) + 1
                                    kk = int((p(3)-frbox(3,1))/dgrid) + 1
                                    finegrid(ii, jj, kk) = finegrid(ii, jj, kk) + 1
                                endif
                            endif
                        endif
                        pp = ll(pp)
                    enddo !! pp loop
                enddo !! i loop
            enddo !! j loop
        enddo !! k loop

#ifdef HVIR
    else

        do k = csbox(3,1), csbox(3,2)
            do j = csbox(2,1), csbox(2,2)
                do i = csbox(1,1), csbox(1,2)
                    pp = hoc(i, j, k)
                    do
                        if (pp == 0) exit
                        if (hpart_vir(pp) == 0) then !! particle is not yet part of a halo
                            p(:) = xvp(:3, pp)
                            dr   = HPOS(:) - p(:)
                            r    = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                            if (r < rsearch) then
                                np_search = np_search + 1
                                pos(np_search, 1:3) = p(:)
                                ilist_vir(np_search) = pp
                                if (r < rrefine) then
                                    ii = int((p(1)-frbox(1,1))/dgrid) + 1
                                    jj = int((p(2)-frbox(2,1))/dgrid) + 1
                                    kk = int((p(3)-frbox(3,1))/dgrid) + 1
                                    finegrid(ii, jj, kk) = finegrid(ii, jj, kk) + 1
                                endif
                            endif
                        endif
                        pp = ll(pp)
                    enddo !! pp loop
                enddo !! i loop
            enddo !! j loop
        enddo !! k loop

    endif
#endif

    !
    ! Fine refined mesh density maximum
    !

    maxfinegrid = 0.

    do k = 1, ngrid(3)
       do j = 1, ngrid(2)
          do i = 1, ngrid(1)
             if (finegrid(i, j, k) > maxfinegrid) then
                maxfinegrid = finegrid(i, j, k)
                HPOS(1) = frbox(1,1) + (i-0.5)*dgrid
                HPOS(2) = frbox(2,1) + (j-0.5)*dgrid
                HPOS(3) = frbox(3,1) + (k-0.5)*dgrid
             endif
             finegrid(i, j, k) = 0. !! Set to zero for next candidate
          enddo
       enddo
    enddo

    !
    ! Sort the particles within the search region based on their distance from the centre.
    !

    if (np_search > max_halo_np) write(*,*) "ERROR: np_search, max_halo_np = ", np_search, max_halo_np

    do i = 1, np_search
        dr = HPOS(:) - pos(i, 1:3)
        pos(i, 4) = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
    enddo

    isortpos(:np_search) = (/ (i, i=1, np_search) /)
    call indexedsort(np_search, pos(:, 4), isortpos(:))
#ifdef HVIR
    if (.not. HALOVIR) then
#endif
        ilist_odc(:np_search) = ilist_odc(isortpos(:np_search)) 
#ifdef HVIR
    else
        ilist_vir(:np_search) = ilist_vir(isortpos(:np_search))
    endif
#endif    

    !
    ! Move one particle at a time until we reach the radius for which the desired overdensity is obtained.
    !

    RODC = 0.
    ITOT = 0

    do j = 2, np_search !! Assume the first particle is part of the halo
    
        odcj = threeover4pi * j * mass_p / pos(j, 4)**3 

        if (odcj <= HODC) then

            odci = threeover4pi * (j-1) * mass_p / pos(j-1, 4)**3

            !! Interpolate halo radius
            r2 = log10(pos(j, 4))
            r1 = log10(pos(j-1, 4))
            d2 = log10(odcj) 
            d1 = log10(odci) 
            w1 = log10(HODC) - d2
            w2 = d1 - log10(HODC)

            RODC = 10**((w1 * r1 + w2 * r2) / (d1 - d2))

            if (pos(j, 4) <= RODC) then
                ITOT = j
            else
                ITOT = j-1        
            endif

            exit

        endif

    enddo

    if (ITOT == 0) search_fail = search_fail + 1 

end subroutine find_halo_particles

! -------------------------------------------------------------------------------------------------------

subroutine find_halo_candidates(tile, ic)
    !
    ! Finds density maxima and returns a sorted array containing potential halo  candidates.
    !

    implicit none

    integer(4), dimension(3) :: offset
    integer(4), dimension(3) :: cic_l,cic_h,tile
    integer(4), intent(inout) :: ic

    integer(4) :: i, j, k, pp, thread
    integer(4) :: ii, ix, iy, iz, ix0, iy0, iz0 
    real(4) :: amtot, denmax

    !! Tile offset in local coordinates

    offset = tile * nf_physical_tile_dim - nf_buf

    !! Initialize density
    thread = 1
    rho_f(:, :, :) = 0.

    !! Limits for mass assignment. Ignore outermost buffer (4 fine cells)
    cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
    cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1

    !! Calculate fine mesh density for tile

    do k = cic_l(3), cic_h(3)
        do j = cic_l(2), cic_h(2)
            do i = cic_l(1), cic_h(1)
                pp = hoc(i, j, k)
#ifdef NGPH
                call fine_ngp_mass(pp, tile)
#else
                call fine_cic_mass(pp, tile)
#endif
            enddo
        enddo
    enddo

    !
    ! Find density maxima 
    !

    do k = 1 + nf_buf, nf_buf + nf_physical_tile_dim     
        do j = 1 + nf_buf, nf_buf + nf_physical_tile_dim
            do i = 1 + nf_buf, nf_buf + nf_physical_tile_dim

                !! Maximum density around nearest cells
                denmax   = maxval(rho_f(i-1:i+1, j-1:j+1, k-1:k+1))

                !! Continue if this cell is a large enough local maxima 
                if (denmax == rho_f(i, j, k) .and. denmax > den_peak_cutoff) then

                    if (ic > max_maxima -2) then
                        write(*,*) "Too many halos: ic, max_maxima", ic, max_maxima
                        exit
                    endif

                    !! Find the fine mesh mass of this peak
                    amtot = 0.
                    ix0 = i
                    iy0 = j
                    iz0 = k
                    do ii = 1, irtot
                        ix = ix0 + idist(1, ii)
                        if (ix < 5 .or. ix > nf_tile-4) cycle
                        iy = iy0 + idist(2, ii)
                        if (iy < 5 .or. iy > nf_tile-4) cycle
                        iz = iz0 + idist(3, ii)
                        if (iz < 5 .or. iz > nf_tile-4) cycle
                        amtot = amtot + rho_f(ix, iy, iz)
                        if (complete_shell .and. rdist(ii) == rdist(ii+1)) cycle
                        if (ii > 18 .and. amtot/(real(ii)) < halo_odc) exit 
                    enddo

                    !! Consider this a halo candidate if the fine mesh mass is large enough
                    if (amtot > mass_p*min_halo_particles/2.) then
                
                        ic = ic + 1

                        !! Store integer coordinates of the peak as well as its mesh mass
                        ipeak(:, ic) = (/real(i), real(j), real(k)/) - 0.5 + offset
                        den_peak(ic) = denmax
                        halo_mesh_mass(ic) = amtot

                    endif

                endif !! denmax

            enddo !! i loop
        enddo !! j loop
    enddo !! k loop

end subroutine find_halo_candidates 

! -------------------------------------------------------------------------------------------------------

subroutine fine_cic_mass(pp, tile)
    !
    ! CIC interpolation of partilces onto PM mesh
    !

    implicit none

    integer(4) :: pp
    integer(4), dimension(3) :: tile

    integer(4), dimension(3) :: i1, i2
    real(4), dimension(3) :: x, offset, dx1, dx2

    offset(:) = - tile(:) * nf_physical_tile_dim + nf_buf !- 0.5 

    do
      if (pp == 0) exit
      x(:) = xvp(1:3, pp) + offset(:)
      i1(:) = floor(x(:)) + 1
      i2(:) = i1(:) + 1
      dx1(:) = i1(:) - x(:)
      dx2(:) = 1 - dx1(:)

      dx1(1) = mass_p * dx1(1)
      dx2(1) = mass_p * dx2(1)

      rho_f(i1(1),i1(2),i1(3)) = rho_f(i1(1),i1(2),i1(3)) &
                                       + dx1(1) * dx1(2) * dx1(3)
      rho_f(i2(1),i1(2),i1(3)) = rho_f(i2(1),i1(2),i1(3)) &
                                       + dx2(1) * dx1(2) * dx1(3)
      rho_f(i1(1),i2(2),i1(3)) = rho_f(i1(1),i2(2),i1(3)) &
                                       + dx1(1) * dx2(2) * dx1(3)
      rho_f(i2(1),i2(2),i1(3)) = rho_f(i2(1),i2(2),i1(3)) &
                                       + dx2(1) * dx2(2) * dx1(3)
      rho_f(i1(1),i1(2),i2(3)) = rho_f(i1(1),i1(2),i2(3)) &
                                       + dx1(1) * dx1(2) * dx2(3)
      rho_f(i2(1),i1(2),i2(3)) = rho_f(i2(1),i1(2),i2(3)) &
                                       + dx2(1) * dx1(2) * dx2(3)
      rho_f(i1(1),i2(2),i2(3)) = rho_f(i1(1),i2(2),i2(3)) &
                                       + dx1(1) * dx2(2) * dx2(3)
      rho_f(i2(1),i2(2),i2(3)) = rho_f(i2(1),i2(2),i2(3)) &
                                       + dx2(1) * dx2(2) * dx2(3)
      pp = ll(pp)
    enddo

end subroutine fine_cic_mass

! -------------------------------------------------------------------------------------------------------

subroutine fine_ngp_mass(pp, tile)
    !
    ! NGP interpolation of partilces onto PM mesh
    !

    implicit none

    integer(4) :: pp
    integer(4), dimension(3) :: tile, i1
    real(4), dimension(3) :: x, offset

    offset(:) = - tile(:) * nf_physical_tile_dim + nf_buf
! removed the half-cell offset so that fine mesh cells will line up with coarse
! mesh cells
!    offset(:)= - tile(:) * nf_physical_tile_dim + nf_buf - 0.5

    do
        if (pp == 0) exit
        x(:) = xvp(1:3,pp) + offset(:)
        i1(:) = floor(x(:)) + 1
        rho_f(i1(1),i1(2),i1(3)) = rho_f(i1(1),i1(2),i1(3)) + mass_p
        pp = ll(pp)
    enddo

end subroutine fine_ngp_mass

! -------------------------------------------------------------------------------------------------------

end program halofind

