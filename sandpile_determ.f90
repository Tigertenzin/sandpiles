program sandpile                                                             ! Program for running NN BTW sandpile in 2D 
    integer, parameter :: dp = selected_real_kind(15, 500)
    
    ! Initializing system-specific parameters  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    real :: lambda                                                          ! Dissipation probability
    integer(kind=8) :: lambdRate
    integer :: epsilon                                                      ! Driving amount (sand added per drive)
    integer(kind=8) :: L                                                    ! Linear system size (system is L by L)
    integer(kind=4), allocatable :: lattice(:,:)                            ! The system configuration (Dimension L x L for 2D system)
    integer(kind=4), allocatable :: latticeT(:,:)                           ! The system configuration (Dimension L x L for 2D system)
    real(kind=4), allocatable :: tAlattice(:,:)                             ! Record of system's time average (Dimension L x L for 2D system)
    integer(kind=8) :: timer                                                ! Total number of avalanches to run the system for. 
    integer(kind=8) :: lambdaTimer                                          ! 1/lambda: number of topplings before a dissipation. 
    character(len=1024) :: filename1                                        ! File name for the lattice after transience. 
    character(len=1024) :: filename2                                        ! File name for the avalanche sizes. 
    character(len=1024) :: filename3                                        ! File name for the avalanche durations. 
    character(len=1024) :: filename4                                        ! File name for the lattice snapshots at 100,000 timesteps. 


    INTEGER, PARAMETER :: LONG=SELECTED_INT_KIND(18)                        ! Necessray variables for pseudo-random number generator
    INTEGER(LONG) :: seed
    REAL rnd48, r
    
    integer(kind=8) :: k=1,j=1, m=1, x=1,y=1, t=1                           ! variables for loop indices

    
    read *,timer, L, lambda, epsilon, seed                                  ! Read in system parameters from console
    lambdRate = int(1/lambda, kind=8)
    
    allocate(lattice(0:L-1, 0:L-1))                                         ! Allocate appropriate indices for lattice (linear size L: 0 -> L-1
    allocate(latticeT(0:L-1, 0:L-1))                                         ! Allocate appropriate indices for lattice (linear size L: 0 -> L-1
    allocate(tAlattice(0:L-1, 0:L-1))                                         ! Allocate appropriate indices for time-averaged lattice (linear size L: 0 -> L-1
    
    ! Prepare files to be written to  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (filename1, fmt='(a,I0,a)') "sp",seed,".txt"                      ! Test file of lattice configuration
    open(unit=4,file=filename1)
    write(4,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(4,*) seed, L, lambda, epsilon, timer
    
    write (filename2, fmt='(a,I0,a)') "sp",seed,"_s.txt"                    ! Test file of avalanche sizes
    open(unit=8,file=filename2)
    write(8,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(8,*) seed, L, lambda, epsilon, timer

    write (filename3, fmt='(a,I0,a)') "sp",seed,"_d.txt"                    ! Test file of avalanche duration
    open(unit=9,file=filename3)
    write(9,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(9,*) seed, L, lambda, epsilon, timer    
    
    write (filename4, fmt='(a,I0,a)') "sp",seed,"_tA.txt"                    ! Test file of avalanche duration
    open(unit=10,file=filename4)
    write(10,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(10,*) seed, L, lambda, epsilon, timer    
    
    ! Set up the 2D lattice !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j=0, L-1                                                             ! Initialize lattice as all 1’s (except the left and right) 
        do k=0,L-1
            lattice(j,k) = 1
            latticeT(j,k) = 1
            tAlattice(j,k) = 1
        End do
    end do 
    
    ! Check lattice was initialized correctly
    

    ! Run dynamics for set number of time steps !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lambdaTimer = 0
    do t = 1,1000000                                                       ! Run system for transience, 1,000,000 avalanches
        call addGrainTransient(seed, L, lambdRate, epsilon, lattice, lambdaTimer)
        
        latticeT = latticeT + lattice
        tAlattice = latticeT / t
        if (mod(t,100000) .eq. 0) then 
            write(10,*) tAlattice(j,:)
        end if
    end do

    do j=0,L+1
        write(4,*) lattice(j,:)
    end do
    
    do t = 1,timer
        call addGrain(seed, L, lambdRate, epsilon, lattice, lambdaTimer)
        
        latticeT = latticeT + lattice
        tAlattice = latticeT / t
        if (mod(t,100000) .eq. 0) then 
            write(10,*) tAlattice(j,:)
        end if
    end do

    close(5) 
    
end program sandpile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extra Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to add energy grain randomly

Subroutine addGrain(seed, L, lambdRate, epsilon, lattice, lambdaTimer)

	integer(kind=8), intent(in) :: lambdRate
	integer(kind=8) :: lambdaTimer
	integer:: epsilon
	integer(kind=8), intent(in) :: L, seed
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer :: x, y
    integer :: sizeS=0, durationD=0 

    x =  nint(rnd48(seed)*(L-1))                                    ! Pick random x-coord to drive
    y =  nint(rnd48(seed)*(L-1))                                    ! Pick random y-coord to drive
    
    lattice(x,y) = lattice(x,y) + epsilon                           ! Bring site (x,y) to failure
    if (lattice(x,y) .gt. 4)  then                                  ! Check if driven site is above threshold, z_(x,y) > 4
    	call relaxLattice(seed, L, lambdRate, sizeS, durationD, lattice, lambdaTimer)          
                                                                    ! Call subroutine to topple active site
    	write(8,'(I0)') sizeS
    	write(9,'(I0)') durationD
    	sizeS=0
    	durationD = 0
	end if

end Subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutin to topple all critical sites 

Subroutine relaxLattice(seed, L, lambdRate, sizeS, durationD, lattice, lambdaTimer)

	logical :: active
	integer(kind=8), intent(in) :: lambdRate
	integer(kind=8) :: lambdaTimer
	integer(kind=8), intent(in) :: L, seed
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer :: j, k
    integer(kind=4) :: toppleThese(0:L-1, 0:L-1)
    integer :: bulkDissipationLayer(0:3)
    integer, intent(inout) :: sizeS,durationD 
    
    ! do j = 0,3                                                      ! Initialize dissipation matrix (a 1x4 where each NN either dissipates or doesn't)
    !     if (rnd48(seed) .ge. lambda) then                           ! If random number above lambda probability, 
    !         bulkDissipationLayer(j) = 1                             ! then set to not dissipate
    !     else		                                                  ! otherwise
    !         bulkDissipationLayer(j) = 0                             ! set to dissipate on a given site....
    !     end if                                                      !
    ! end do                                                          !
    
    active = .true.
    do while (active .eqv. .true.)										! Continue relaxaing until active = .false., i.e. lattice is no longer active
	    do j = 0, L-1                                                   ! Initialize the toppling matrix
	    do k = 0, L-1                                                   ! this will be added to the lattice after all currently active sites are toppled
	        toppleThese(j,k) = 0                                        ! (hence in parallel)
	    end do                                                          !
	    end do                                                          !
	    
	    do j=0, L-1                                                     ! Run over lattice, check for sites to topple
	    do k=0, L-1
	        if (lattice(j,k) .gt. 4) then                               ! If a given site is above threshold...            
	            toppleThese(j,k) = toppleThese(j,k)-4                   ! active site loses 4 grains
	            sizeS = sizeS + 1                                       ! Add to total size, s, of avalanche (because site topples)
	            
	            lambdaTimer = lambdaTimer + 1                           ! Incriment lambda timer. 
	            if (mod(lambdaTimer,lambdRate) .ne. 0) then             ! If lambda timer is modulus of rate, then dont dissipate... 
    	            
    	            if (j .ne. 0) then                                      ! Left NN gains 1 (if not over boundary)
    	               toppleThese(j-1,k) = toppleThese(j-1,k) + 1
    	            end if
    	            
    	            if (j .ne. L-1) then                                    ! Right NN gains 1 (if not over boundary)
    	                toppleThese(j+1,k) = toppleThese(j+1,k) + 1
    	            end if
    	            
    	            if (k .ne. 0) then                                      ! Top NN gains 1 (if not over boundary)
    	                toppleThese(j,k-1) = toppleThese(j,k-1) + 1
    	            end if
    	            
    	            if (k .ne. L-1) then                                    ! Bottom NN gains 1 (if not over boundary)
    	                toppleThese(j,k+1) = toppleThese(j,k+1) + 1
    	            end if
	            end if 
	        end if
	    end do
	    end do
	    
	    lattice = lattice + toppleThese                                 ! Add total Toppling Matrix to current lattice
	    durationD = durationD + 1										! Add to total duration, D, of avalanche 
	    
	    active = .false.
	    do j=0,L-1                                                      ! Check each site if over 4
	    do k=0,L-1
	    	if (lattice(j,k) .gt. 4) then
	    		active = .true.
	    		exit
	    	end if
	    end do 
	    	if (active .eqv. .true.) then
	     		exit
	    	end if 
	    end do

	end do


end Subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to add energy grain randomly

Subroutine addGrainTransient(seed, L, lambdRate, epsilon, lattice, lambdaTimer)

   	integer(kind=8), intent(in) :: lambdRate
	integer(kind=8) :: lambdaTimer
	integer:: epsilon
	integer(kind=8), intent(in) :: L, seed
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer :: x, y
    integer :: sizeS=0, durationD=0 

    x =  nint(rnd48(seed)*(L-1))                                    ! Pick random x-coord to drive
    y =  nint(rnd48(seed)*(L-1))                                    ! Pick random y-coord to drive
    
    lattice(x,y) = lattice(x,y) + epsilon                           ! Bring site (x,y) to failure
    if (lattice(x,y) .gt. 4)  then                                  ! Check if driven site is above threshold, z_(x,y) > 4
    	call relaxLattice(seed, L, lambdRate, sizeS, durationD, lattice, lambdaTimer)          
                                                                    ! Call subroutine to topple active site
    	sizeS=0
    	durationD = 0
	end if
	
end Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extra Functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION rnd48(seed)

  IMPLICIT NONE
  INTEGER, PARAMETER :: LONG=SELECTED_INT_KIND(18)
  INTEGER, PARAMETER :: REAL8=SELECTED_REAL_KIND(15,300)
  INTEGER(LONG), PARAMETER :: SEED_M=25214903917_LONG, SEED_A=11_LONG
  REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
  REAL rnd48
  INTEGER(LONG) :: seed

  seed=IAND(SEED_M*seed+SEED_A,281474976710655_LONG)
  rnd48=TWONEG48*seed

END FUNCTION rnd48
