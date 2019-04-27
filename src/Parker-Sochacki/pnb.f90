!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOTAL_ENERGY %%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE ENERGY (n, G, M, X, V, KE, PE, Total_Energy)
!
!-----------------------------------------------------------------------
!
!  Computes the total energy of an N-body system as the sum of 
!  kinetic (KE) and potential energies (PE)
!
!  INPUTS:  n = integer number of bodies
!           G = real constant, universal gravitational constant
!           M = n vector containing masses
!           X = 3 by n array containing position coordinates for bodies
!           V = 3 by n array containing velocity components  for bodies
!
!  OUTPUT: Total_Energy = KE + PE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), INTENT(IN) :: G
      REAL(KIND=8), DIMENSION(n),   INTENT(IN) :: M
      REAL(KIND=8), DIMENSION(3,n), INTENT(IN) :: X, V
      REAL(KIND=8), INTENT(OUT) :: KE, PE, Total_Energy
!
!  local variables
!
      INTEGER :: i, j
      REAL(KIND=8) :: Vsq, radius
      REAL(KIND=8) :: dX1, dX2, dX3
      
      KE = 0.0d0  ! sum the kinetic energies
      DO j = 1,n
         Vsq = V(1,j)*V(1,j) + V(2,j)*V(2,j) + V(3,j)*V(3,j)
         KE = KE + M(j) * Vsq 
      END DO
      KE = 0.5d0 * KE

      PE = 0.0d0  ! sum the potential energies
      DO i = 1,n
         DO j = 1,i-1
            dX1 = X(1,i) - X(1,j)
            dX2 = X(2,i) - X(2,j)
            dX3 = X(3,i) - X(3,j)
            radius = dX1*dX1 + dX2*dX2 + dX3*dX3
            radius = SQRT(radius) ! separation distance
            PE = PE - G*M(i)*M(j) / radius
         END DO
      END DO

      Total_Energy = KE + PE

      RETURN
   END SUBROUTINE ENERGY
!
!%%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR MOMENTA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE LINEAR_MOMENTA (n, M, V, LM)
!
!-----------------------------------------------------------------------
!
!  Computes the three components of linear momentum
!
!  INPUTS:  n = integer number of bodies
!           M = n vector containing masses
!           V = 3 by n array containing velocity components for bodies
!
!  OUTPUT: Three components of linear momentum (LM1,LM2,LM3)
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), DIMENSION(n),   INTENT(IN) :: M
      REAL(KIND=8), DIMENSION(3,n), INTENT(IN) :: V
      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: LM

      INTEGER :: i, j

      DO i = 1,3
         LM(i) = 0.0D0
         DO j = 1,n
            LM(i) = LM(i) + M(j)*V(i,j) 
         END DO
      END DO

      RETURN
   END SUBROUTINE LINEAR_MOMENTA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%% ANGULAR MOMENTA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE ANGULAR_MOMENTA (n, M, R, V, AM)
!
!-----------------------------------------------------------------------
!
!  Computes the three components of angular momentum
!
!  INPUTS:  n = integer number of bodies
!           M = n vector containing masses
!           R = 3 by array contain radius vectors for n bodies
!           V = 3 by n array containing velocity components for bodies
!
!  OUTPUT: Three components of angular momentum (AM1,AM2,AM3)
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), DIMENSION(n),   INTENT(IN) :: M
      REAL(KIND=8), DIMENSION(3,n), INTENT(IN) :: R, V
      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: AM

      INTEGER :: i, j

      AM(1:3) = 0.0D0
      DO j = 1,n  ! loop over number of bodies
         AM(1) = AM(1) + M(j)*(R(2,j)*V(3,j)-R(3,j)*V(2,j))
         AM(2) = AM(2) + M(j)*(R(3,j)*V(1,j)-R(1,j)*V(3,j))
         AM(3) = AM(3) + M(j)*(R(1,j)*V(2,j)-R(2,j)*V(1,j))
      END DO

      RETURN
   END SUBROUTINE ANGULAR_MOMENTA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELTA_T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Computes allowable time-step increment (dt) to satisfy global error
!  tolerance on position, given order (m) and velocity vector (V)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE DELTA_T (n, m, u, V, dt)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n         ! system size
      INTEGER, INTENT(IN) :: m         ! maximum order
      REAL(KIND=8), INTENT(IN) :: u    ! error tolerance per unit time 
      REAL(KIND=8), DIMENSION(1:3,1:n), INTENT(IN) :: V ! velocity vector
      REAL(KIND=8), INTENT(OUT) :: dt  ! time increment

      INTEGER :: i, j
      REAL(KIND=8) :: Vtemp, Vmax

      Vmax = 0.0
      EUCLIDEAN_NORM: DO j = 1, n
         Vtemp = V(1,j)*V(1,j)
         DO i = 2,3
            Vtemp = Vtemp + V(i,j)*V(i,j)
         END DO
         IF (Vtemp > Vmax) Vmax = Vtemp
      END DO EUCLIDEAN_NORM
      dt = (u / SQRT(Vmax))**(1./m)

      RETURN

   END SUBROUTINE DELTA_T
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATIONAL_COST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Estimates the computational cost (flops/unit time) of a given time step
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE COMPUTATIONAL_COST (n, m, dt, cost, ComputeFlops)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n, m
      REAL(KIND=8), INTENT(IN) :: dt
      REAL(KIND=8), INTENT(OUT) :: cost 
      LOGICAL, INTENT(IN) :: ComputeFlops
!
! local variables
!
      REAL(KIND=8) :: wt_a=1.0D0, wt_m=1.0D0, wt_s=1.0D0, wt_d=1.0D0
      REAL(KIND=8) :: weighted_flops
      INTEGER :: add, sub, mul, div
      INTEGER :: n2, mm1, nm, mp1, mp2, n2m, mmp1, nnm1, nnm1m, &
                 mp1mp2, n2mmp1, nnm1mmp1, nnm1mp1mp2
!
! compute number of operations
!
      IF (ComputeFlops) THEN
         n2 = n*n
         nnm1 = n * (n-1)
         mm1 = m - 1
         nm  = n*mm1
         mp1 = mm1 + 1
         mp2 = mm1 + 2
         n2m = n2*mm1
         mmp1 = mm1 * mp1
         nnm1m = nnm1 * mm1
         mp1mp2 = mp1 * mp2
         n2mmp1 = n2 * mmp1
         nnm1mmp1 = nnm1*mmp1
         nnm1mp1mp2 = nnm1*mp1mp2
         add = 4.0*n2+6.0*nm+0.5*nnm1m+1.5*n2mmp1+1.25*nnm1mp1mp2
         sub = 5.0*n2+nm+3.0*n2mmp1+nnm1mmp1+2.*nnm1mp1mp2
         mul = n2+9.0*nm+nnm1+6.0*n2m+nnm1m+1.5*n2mmp1+0.25*nnm1mmp1 &
             + 1.25*nnm1mp1mp2
         div = nm
         weighted_flops = add*wt_a + mul*wt_m + sub*wt_s + div*wt_d
      END IF

      cost = weighted_flops / dt

      RETURN

   END SUBROUTINE COMPUTATIONAL_COST
 
PROGRAM Parallel_N_body
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Parallel ODE solver for the N-body problem based on Picard iteration
!  ala Parker + Sochacki + Rudmin
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Serial-Code Authors:  Dave Pruett, Joe Rudmin, and Justin Lacy
!            JMU Departments of Mathematics & Statistics and Physics
!            pruettcd@jmu.edu / rudminjw@jmu.edu
!
!  Parallel-Code Authors: Dave Pruett and Bill Ingham
!            JMU Departments of Math & Stats and Physics & Astronomy
!            pruettcd@jmu.edu / inghamwh@jmu.edu
!            Refinements: 10/2009
!
!  Notes: 1) This version requires the input of the maximum order (mo)
!         2) Parameter-free version
!         3) Receive counts for GATHER and ALLGATHER are *per process*
!         4) Because Fortran arrays have column preference, GATHER
!            and ALLGATHER fill up arrays by columns.
!         5) Modified 06/2008 to back off of unit round-off error by
!            one order of magnitude.
!         6) Instrumented 09/2009 for computation of total energy
!            on head node.
!         7) Modified 10/15/09 for cases in which n is not an integer
!            multiple of number of processes np.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   IMPLICIT NONE
   INCLUDE 'mpif.h'
   INTEGER :: ierr, status(MPI_STATUS_SIZE)
   INTEGER :: np                  ! number of processes
   INTEGER :: npmax               ! maximum no. of processes allowed
   INTEGER :: id                  ! process id

!  REAL(KIND=8), PARAMETER :: G = 6.673D-11 ! grav. const. (SI)
   REAL(KIND=8), PARAMETER :: G = 1.0D0     ! grav. const.

   REAL(KIND=8) :: tmin,  &       ! initial value of t
                   tmax,  &       ! terminal value of t
                   t,     &       ! time-like variable of integration
                   tout,  &       ! next output time
                   dt,    &       ! time-like increment
                   dtout, &       ! time increment for output only 
                   t_temp,&       ! trial value of time
                   dX,    &       ! separation in a component of X
                   s,     &       ! separation distance btw. centers
                   EV,    &       ! maximum absolute error on velocity
                   u,     &       ! absolute error per unit of time
                   Vmax,  &       ! per-step Euclidean err. in velocity  
                   Vtemp, &       ! scalar temporary for velocity
                   cost,  &       ! computational cost for current order
                   prev_cost, &   ! computational cost for previous order
                   sum_cost       ! cumulative cost

   REAL(KIND=8), DIMENSION(3) :: LM ! linear momenta
   REAL(KIND=8), DIMENSION(3) :: AM ! angular momenta
   REAL(KIND=8) :: KE,    &         ! kinetic energy diagnostic
                   PE,    &         ! potential energy diagnostic
                   Total_Energy     ! total energy diagnostic

   INTEGER :: n,        &         ! number of bodies
              nout,     &         ! number of bodies output
              mo,       &         ! max. order of Maclaurin series 
              mo1,      &         ! order of series plus one
              i,        &         ! index for 3-vectors X and V
              j,        &         ! global index of body
              k,        &         ! also global index of body
              m,        &         ! index of order
              outcount, &         ! count for output calls
              steps,    &         ! time-step counter
              sv                  ! allocate status variable

   INTEGER :: my_n,      &        ! storage length per process
              my_nb,     &        ! number of bodies per process
              my_j,      &        ! local index of body for process
              m3,        &        ! my_n*3 (send/receive counts for gather)
              j0,        &        ! global starting index for process
              j1,        &        ! global stopping index for process
              b                   ! bodies per process

   LOGICAL :: TimeToQuit = .FALSE.  ! termination control variable
   LOGICAL :: TimeToPrint           ! print control variable
   LOGICAL :: ComputeFlops = .TRUE. !  
   LOGICAL :: Terminate             ! exit condition for order loop  
   LOGICAL :: Diagnostics           ! turns on diagnostics (energy, etc)

   REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: Mass ! mass(j)
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: X    ! position
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: V    ! velocity
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: IS1, &
                                                  IS2, &
                                                  IS3  ! separations
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: A    ! action term
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: all_X, all_V
   REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: delta_X, delta_V
!
!----------------------------- TEMPORARIES -----------------------------
!
   INTEGER :: mm1, mmq, q, mm1mq                      ! counters
   REAL(KIND=8) :: mi                                 ! reciprocal of m
   REAL(KIND=8) :: sum1, sum2                         ! accumulators
   REAL(KIND=8) :: start_time, stop_time, system_time ! for timing
   CHARACTER(12) :: filename1                         ! for output
!
!-----------------------------INITIALIZATION----------------------------
!
!  define system size
!
   CALL MPI_init(ierr)
   CALL MPI_comm_size(MPI_COMM_WORLD, np, ierr) ! everyone should know
   CALL MPI_comm_rank(MPI_COMM_WORLD, id, ierr) ! who am I?

   IF (id == 0) THEN  ! root process
      READ (5,*) n, nout
      WRITE (6,*) ' number of bodies, n = ', n
      IF (nout < 0) THEN
         nout = n
      END IF
      READ (5,*) mo
      WRITE (6,*) ' number of processes = ', np
      WRITE (6,*) ' maximum order of Maclaurin expansion, mo = ', mo
      mo1 = mo + 1
      IF (np > n) THEN
            WRITE (6,*) ' Number of processes > number of bodies: STOP '
            TimeToQuit = .TRUE.    
      ELSE IF (MOD(n,np) == 0) THEN ! perfectly load-balanced
         my_n = n / np
      ELSE   ! impose efficiency constraint if not load-balanced
         npmax = CEILING(SQRT(REAL(n))) + 1
         IF (np > npmax) THEN ! violates efficiency/distribution constaint
            WRITE (6,*) ' '
            WRITE (6,*) ' EITHER: '
            WRITE (6,*) '    (1) specify number of processes to '
            WRITE (6,*) '        be an integer divisor of ', n
            WRITE (6,*) ' OR: '
            WRITE (6,*) '    (2) reduce number of processes to '
            WRITE (6,*) '        at most ', npmax
            WRITE (6,*) ' '
            WRITE (6,*) ' Aborting due to computational inefficiency '
            TimeToQuit = .TRUE.    
         ELSE
            my_n = FLOOR( REAL(n + (np - 1)) / REAL(np))
         END IF
      END IF
   END IF

   CALL MPI_BCAST (TimeToQuit,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   IF (TimeToQuit) STOP

   CALL MPI_BCAST (n   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) ! from root
   CALL MPI_BCAST (my_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) ! from root
   CALL MPI_BCAST (mo1, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) ! from root
   mo = mo1 - 1                ! everyone computes
   m3 = my_n*3                 ! everyone computes

   IF (id == np-1) THEN        
      my_nb = n - (np-1)*my_n  ! very last process gets the leftovers
      IF (my_nb < 0) THEN ! should never be encountered w. eff. constraint
         WRITE (6,*) ' Last process assigned negative number of bodies. '  
         WRITE (6,*) ' Reduce number of processes: STOP '
         TimeToQuit = .TRUE.
      END IF
   ELSE
      my_nb = my_n             ! everyone but last gets the same
   END IF
   j0 = id*my_n                ! everyone computes his starting index
   j1 = MIN((id+1)*my_n,n)     ! everyone computes his stopping index 
   WRITE (6,*) ' Process ',id, ' gets ', my_nb, ' bodies: from ', j0+1, j1
!
!  abnormal stop
!
   CALL MPI_BCAST (TimeToQuit,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   IF (TimeToQuit) STOP
!
!  allocate same storage to everyone (even if not fully unused)
!
   ALLOCATE (Mass(n), STAT=sv)
   IF (sv /= 0) STOP ' allocation failure for Mass'
   ALLOCATE (X(3,my_n,0:mo1), V(3,my_n,0:mo1), STAT=sv)
   IF (sv /= 0) STOP ' allocation failure for X and/or V'
   ALLOCATE (A(n,my_n,0:mo1), STAT=sv)
   IF (sv /= 0) STOP ' allocation failure for A'
   ALLOCATE (IS1(n,my_n,0:mo1), IS2(n,my_n,0:mo1), IS3(n,my_n,0:mo1), STAT=sv)
   IF (sv /= 0) STOP ' allocation failure for IS1, IS2, and/or IS3'
   ALLOCATE (delta_X(n,my_n,3,0:mo1), delta_V(n,my_n,3,0:mo1), STAT=sv)
   IF (sv /= 0) STOP ' allocation failure for delta_X, delta_V'
   ALLOCATE (all_X(3,my_n*np), all_V(3,my_n*np), STAT=sv)
   IF (sv /= 0) STOP ' allocation failure for all_X, all_V'
!
!  read data and define control parameters
!
   IF (id == 0) THEN
      READ (5,*) tmin, tmax, dtout
      WRITE (6,*) ' time-like interval [tmin,tmax] = ', tmin, tmax
      WRITE (6,*) ' print interval, dtout = ', dtout
      READ (5,*) EV, Diagnostics 
      IF (EV < 0.0D0) EV = EPSILON(EV)*10.0D0 ! back off one order
      WRITE (6,*) ' absolute error limit on velocities, EV = ', EV
      WRITE (6,*) ' diagnostics trigger, Diagnostics = ', Diagnostics
      u = EV / (tmax - tmin)
      DO j = 1,n
         READ (5,*) Mass(j), (all_X(i,j),i=1,3), (all_V(i,j),i=1,3)
      END DO
!
!  open output files (host node only)
!
      OPEN (2,file='comp_cost_est',status='UNKNOWN')
      WRITE (2,*) ' step, order, cost, cum_cost, dt, t '
      OPEN (7,file='sys_summary',status='UNKNOWN')
   END IF
!
!------------------ OUTPUT THE INITIAL CONDITION -------------------
!
   steps = 0
   outcount = 1
   t = tmin
   IF (id == 0) THEN
      DO j = 1,nout
         WRITE (filename1,'(8hparticle,i4.4)') j
         OPEN (8,FILE=filename1,STATUS='UNKNOWN')
         WRITE (8,'(    4(1x,g19.13))') t, (all_X(i,j),i=1,3)
         WRITE (8,'(20x,3(1x,g19.13))')    (all_V(i,j),i=1,3)
         CLOSE (8)
      END DO
      IF (Diagnostics) OPEN ( 9,FILE='energies',STATUS='UNKNOWN')
      IF (Diagnostics) OPEN (10,FILE='momenta1',STATUS='UNKNOWN')
      IF (Diagnostics) OPEN (11,FILE='momenta2',STATUS='UNKNOWN')
      start_time = MPI_WTIME() ! start the clock
   END IF
!
!---------------------- TIME INTEGRATION ---------------------------
!
!  Distribute the initial data
!
   CALL MPI_BCAST( Mass,    n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(all_X,np*m3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(all_V,np*m3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   TIME_ADVANCE: DO 

      steps = steps + 1
      tout = tmin + dtout*outcount
      IF (id == 0) prev_cost = 1.0D28 ! set outrageously large to initialize

      ORDER_ZERO_POSITIONS_AND_VEL: DO my_j = 1, my_nb
         j = j0 + my_j
         X(1:3,my_j,0) = all_X(1:3,j)
         V(1:3,my_j,0) = all_V(1:3,j)
      END DO ORDER_ZERO_POSITIONS_AND_VEL

      del_XV0: DO my_j = 1, my_nb
         j = j0 + my_j
         DO k = 1, n
            delta_X(k,my_j,1:3,0) = all_X(1:3,k) - all_X(1:3,j)
            delta_V(k,my_j,1:3,0) = all_V(1:3,k) - all_V(1:3,j)
         END DO  
         delta_X(j,my_j,1:3,0) = 0.0D0
         delta_V(j,my_j,1:3,0) = 0.0D0
      END DO del_XV0

      INITIAL_SEPARATIONS: DO my_j = 1, my_nb  ! IS1 is symmetric
         j = j0 + my_j
         DO k = 1, j-1
            s = 0.0D0
            DO i = 1, 3
               dX = delta_X(k,my_j,i,0)
               dX = dX * dX
               s = s + dX
            END DO   
            s = SQRT(s)
            IS1(k,my_j,0) = 1.0D0 / s
         END DO
         IS1(j,my_j,0)         = 0.0D0
         DO k = j+1, n
            s = 0.0D0
            DO i = 1, 3
               dX = delta_X(k,my_j,i,0)
               dX = dX * dX
               s = s + dX
            END DO   
            s = SQRT(s)
            IS1(k,my_j,0) = 1.0D0 / s
         END DO  
      END DO INITIAL_SEPARATIONS

      ORDER_ZERO_SEPARATIONS: DO my_j = 1, my_nb  
         j = j0 + my_j
         DO k = 1, n
            IS2(k,my_j,0) = IS1(k,my_j,0) * IS1(k,my_j,0) ! symmetric
            IS3(k,my_j,0) = IS1(k,my_j,0) * IS2(k,my_j,0) ! symmetric
         END DO
         IS2(j,my_j,0) = 0.0D0
         IS3(j,my_j,0) = 0.0D0
      END DO ORDER_ZERO_SEPARATIONS

      ORDER_ZERO_ACTION: DO my_j = 1, my_nb
         j = j0 + my_j
         DO k = 1, n
            sum1 = 0.0D0
            DO i = 1, 3
               sum1 = sum1 + delta_X(k,my_j,i,0) * delta_V(k,my_j,i,0)
            END DO
            A(k,my_j,0) = sum1
         END DO
         A(j,my_j,0) = 0.0D0
      END DO ORDER_ZERO_ACTION

      HIGHER_ORDER: DO m = 1, mo1   ! want this loop on outside
         mi = 1.D0 / REAL(m)
         mm1 = m - 1

         POSITIONS: DO my_j = 1, my_nb
            X(1:3,my_j,m) = mi * V(1:3,my_j,mm1) ! V lags X in order
         END DO POSITIONS

         CALL MPI_ALLGATHER (X(1:3,1:my_n,m),m3,MPI_DOUBLE_PRECISION, &
              all_X,m3,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

         DEL_X: DO my_j = 1, my_nb  ! del_X is antisymmetric
            j = j0 + my_j
            DO k = 1, n
               delta_X(k,my_j,1:3,m) = all_X(1:3,k) - all_X(1:3,j)
            END DO
            delta_X(j,my_j,1:3,m) = 0.0D0
         END DO DEL_X

         VELOCITIES: DO my_j = 1,my_nb
            DO i = 1, 3
               sum1 = 0.0D0
               DO k = 1, n
                  sum2 = 0.0D0
                  DO q = 0, mm1
                     mm1mq = mm1 - q
                     sum2 = sum2 + delta_X(k,my_j,i,q)*IS3(k,my_j,mm1mq)
                  END DO
                  sum1 = sum1 + mi * G * Mass(k) * sum2 
               END DO
               V(i,my_j,m) = sum1
            END DO
         END DO VELOCITIES

         CALL MPI_ALLGATHER (V(1:3,1:my_n,m),m3,MPI_DOUBLE_PRECISION, &
              all_V,m3,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

         DEL_V: DO my_j = 1,my_nb  ! del_V is antisymmetric
            j = j0 + my_j 
            DO k = 1, n
               delta_V(k,my_j,1:3,m) = all_V(1:3,k) - all_V(1:3,j)
            END DO
            delta_V(j,my_j,1:3,m) = 0.0D0
         END DO DEL_V

         SEPARATIONS: DO my_j = 1, my_nb
            j = j0 + my_j
            DO k = 1, n
               sum1 = 0.0D0      ! compute IS1 
               DO q = 0, mm1
                  mm1mq = mm1 - q
                  sum1 = sum1 - IS3(k,my_j,q)*A(k,my_j,mm1mq) 
               END DO
               IS1(k,my_j,m) = mi * sum1
               sum1 = 0.0D0      ! compute IS2
               DO q = 0, m
                 mmq = m - q
                 sum1 = sum1 + IS1(k,my_j,q)*IS1(k,my_j,mmq)
               END DO
               IS2(k,my_j,m) = sum1
               sum1 = 0.0D0      ! compute actions A and IS3
               sum2 = 0.0D0
               DO q = 0, m
                  mmq = m - q
                  sum2 = sum2 + IS2(k,my_j,q)*IS1(k,my_j,mmq)
                  DO i = 1, 3
                     sum1 = sum1 + delta_X(k,my_j,i,q) * delta_V(k,my_j,i,mmq)
                  END DO
               END DO
               A(k,my_j,m)   = sum1
               IS3(k,my_j,m) = sum2
            END DO
            A(j,my_j,m)   = 0.0D0 ! overwrite diagonal values to zero
            IS1(j,my_j,m) = 0.0D0 ! overwrite diagonal values to zero
            IS2(j,my_j,m) = 0.0D0 ! overwrite diagonal values to zero
            IS3(j,my_j,m) = 0.0D0 ! overwrite diagonal values to zero
         END DO SEPARATIONS
!
! determine whether or not to stop at order m based on relative cost
!
         Terminate = .FALSE.  ! everyone assumes not time to quit
         IF (id == 0) THEN    ! host node assesses the situation
            IF (m > 2) THEN   ! always at least second order
               CALL DELTA_T (n, m-1, u, all_V, dt) ! everyone has all_V 
               CALL COMPUTATIONAL_COST (n, m, dt, cost, ComputeFlops) 
               IF (cost < prev_cost) THEN
                  ! WRITE (6,*) ' m, cost = ', m, cost
                  prev_cost = cost
               ELSE
                  mo = m-1
                  ! WRITE (6,*) ' final mo, cost = ', mo, cost
                  Terminate = .TRUE.
               END IF
            END IF
         END IF

         CALL MPI_BCAST (Terminate,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
         IF (Terminate) THEN  ! everyone tests exit condition
            CALL MPI_BCAST (mo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            EXIT  ! exit before maximum order attained
         END IF

      END DO HIGHER_ORDER
!
!-------------------- UPDATE VIA HORNER'S ALGORITHM --------------------
!
! compute the maximum time increment (dt) allowed by the error criterion
!
      IF (id == 0) THEN
!
! adjust the time increment (dt) if necessary to hit desired output time
!
         t_temp = t + dt
         IF (t_temp > tmax) THEN        ! final time exceeded ?
            dt = tmax - t
            t = tmax
            TimeToPrint = .TRUE.
            TimeToQuit  = .TRUE.
         ELSE IF (dtout > 0.0D0) THEN   ! desired output time ?
            IF (t_temp > tout) THEN
               dt = tout - t
               t = tout
               TimeToPrint = .TRUE.
            ELSE
               t = t_temp
               TimeToPrint = .FALSE.
            END IF
         ELSE                           ! print every time step (dtout < 0)
            t = t_temp
            TimeToPrint = .TRUE.
         END IF
         sum_cost = sum_cost + cost*dt
         WRITE (2,'(1x,i6,1x,i6,4(1x,g12.6))') steps, mo, &
                                               cost, sum_cost, dt, t
      END IF

      CALL MPI_BCAST (dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST (TimeToQuit,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!
! Horner update with adjusted time increment (dt)
!
      HORNER: DO m = mo-1, 0, -1 ! accumulate in highest order register   
         X(1:3,1:my_n,mo) = dt*X(1:3,1:my_n,mo) + X(1:3,1:my_n,m)
         V(1:3,1:my_n,mo) = dt*V(1:3,1:my_n,mo) + V(1:3,1:my_n,m)
      END DO HORNER
!
! Reset initial conditions; every process ready to start over
!
      CALL MPI_ALLGATHER (X(1:3,1:my_n,mo),m3,MPI_DOUBLE_PRECISION, &
              all_X,m3,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLGATHER (V(1:3,1:my_n,mo),m3,MPI_DOUBLE_PRECISION, &
              all_V,m3,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!
!------------------- PERIODICALLY OUTPUT POSITIONS ---------------------
!
      IF (id == 0) THEN
         WRITE (6,'(2(1x,i5),2(1x,g20.14))') steps, mo, dt, t
         IF (diagnostics) THEN
            CALL ENERGY (n,G,Mass,all_X,all_V, KE,PE,Total_Energy)
            WRITE ( 9,'(4(1x,G19.13))') t, KE, PE, Total_Energy
            CALL  LINEAR_MOMENTA (n, Mass, all_V, LM)
            WRITE (10,'(4(1x,G19.13))') t, LM(1:3)
            CALL ANGULAR_MOMENTA (n, Mass, all_X, all_V, AM)
            WRITE (11,'(4(1x,G19.13))') t, AM(1:3)
         END IF
         IF (TimeToPrint) THEN
            DO j = 1,nout
               WRITE (filename1,'(8hparticle,i4.4)') j
               OPEN (8,FILE=filename1,STATUS='OLD',POSITION='APPEND')
               WRITE (8,'(    4(1x,g19.13))') t, (all_X(i,j),i=1,3)
               WRITE (8,'(20x,3(1x,g19.13))')    (all_V(i,j),i=1,3)
               CLOSE (8)
            END DO
            outcount = outcount + 1
         END IF
      END IF
      IF (TimeToQuit) EXIT 

   END DO TIME_ADVANCE

   IF (id == 0) THEN
      stop_time = MPI_WTIME()
      system_time = stop_time - start_time
      WRITE (6,*) ' system time = ', system_time
      WRITE (7,*) n, np, mo1, system_time, sum_cost
   END IF
!
!---------------------------- CLEAN UP ------------------------------
!
   DEALLOCATE (Mass,X,V, A,IS1,IS2,IS3, delta_X,delta_V, all_X,all_V)
   IF (id == 0) THEN
      CLOSE (2)
      CLOSE (7)
      IF (Diagnostics) CLOSE ( 9)
      IF (Diagnostics) CLOSE (10)
      IF (Diagnostics) CLOSE (11)
   END IF
 
   CALL MPI_FINALIZE(ierr)

   STOP

END PROGRAM Parallel_N_body
