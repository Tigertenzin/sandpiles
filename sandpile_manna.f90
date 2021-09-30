program sandpile                                                             ! Program for running NN BTW sandpile in 2D 
    integer, parameter :: dp = selected_real_kind(15, 500)
    
    ! Initializing system-specific parameters  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    real :: lambda                                                          ! Dissipation probability
    real(kind=4) :: tmMetric                                                ! TM Metric
    integer(kind=4) :: lambda01
    integer(kind=8) :: lambdRate
    integer :: epsilon                                                      ! Driving amount (sand added per drive)
    integer(kind=8) :: L                                                    ! Linear system size (system is L by L)
    integer(kind=8) :: pU                                                   ! Total number of plate updates
    integer(kind=4), allocatable :: lattice(:,:)                            ! The system configuration (Dimension L x L for 2D system)
    integer(kind=4), allocatable :: latticeT(:,:)                           ! The system configuration (Dimension L x L for 2D system)
    real(kind=4), allocatable :: tAlattice(:,:)                             ! Record of system's time average (Dimension L x L for 2D system)
    integer(kind=8) :: timer                                                ! Total number of avalanches to run the system for. 
    integer(kind=8) :: lambdaTimer                                          ! 1/lambda: number of topplings before a dissipation. 
    character(len=1024) :: filename1                                        ! File name for the lattice after transience. 
    character(len=1024) :: filename2                                        ! File name for the avalanche sizes. 
    character(len=1024) :: filename3                                        ! File name for the avalanche durations. 
    character(len=1024) :: filename4                                        ! File name for the lattice snapshots at 100,000 timesteps. 
    character(len=1024) :: filename5                                        ! File name for height of several selected sites
    character(len=1024) :: filename6                                        ! File name for TM metric


    INTEGER, PARAMETER :: LONG=SELECTED_INT_KIND(18)                        ! Necessray variables for pseudo-random number generator
    INTEGER(LONG) :: seed
    REAL rnd48, r
    
    integer(kind=8) :: k=1,j=1, m=1, x=1,y=1, t=1                           ! variables for loop indices

    
    read *,timer, L, lambda, epsilon, seed                                  ! Read in system parameters from console
    
    allocate(lattice(0:L-1, 0:L-1))                                         ! Allocate appropriate indices for lattice (linear size L: 0 -> L-1
    allocate(latticeT(0:L-1, 0:L-1))                                         ! Allocate appropriate indices for lattice (linear size L: 0 -> L-1
    allocate(tAlattice(0:L-1, 0:L-1))                                         ! Allocate appropriate indices for time-averaged lattice (linear size L: 0 -> L-1
    
    ! Prepare files to be written to  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(6,*) "opening the file"

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
    
    write (filename4, fmt='(a,I0,a)') "sp",seed,"_tA.txt"                    ! Test file of rolling time Averaged lattice configus
    open(unit=10,file=filename4)
    write(10,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(10,*) seed, L, lambda, epsilon, timer    
    
    write (filename5, fmt='(a,I0,a)') "sp",seed,"_tS.txt"                    ! Test file of height of chosen sites vs plate updates
    open(unit=11,file=filename5)
    write(11,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(11,*) seed, L, lambda, epsilon, timer
    
    write (filename6, fmt='(a,I0,a)') "sp",seed,"_tm.txt"                ! Test file of height of chosen sites vs plate updates
    open(unit=12,file=filename6)
    write(12,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(12,*) seed, L, lambda, epsilon, timer
    
    ! Set up the 2D lattice !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(6,*) "Setting up lattices"
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
    pU = 0
    ! write(6,*) "Start running the transience"

    do t = 1,1000000                                                       ! Run system for transience, 1,000,000 avalanches
        call addGrainTransient(seed, L, lambda, epsilon, lattice, lambdaTimer, pU) ! add a grain routin
        
        latticeT = latticeT + lattice
        tAlattice = real(latticeT, kind=4) / t
        
        call calc_tmMetric(L, t, tAlattice)

        if (mod(t,50000) .eq. 0) then 
            write(10,*) tAlattice(:,:)
        end if
    end do

    ! write(6,*) "Record the latice post-transience"
    do j=0,L+1                                                              ! Record lattice after exiting transience
        write(4,*) lattice(j,:)
    end do
    
    do t = 1,timer
        call addGrain(seed, L, lambda, epsilon, lattice, lambdaTimer, pU)
        
        latticeT = latticeT + lattice
        tAlattice = real(latticeT, kind=4) / t
        if (mod(t,50000) .eq. 0) then 
            write(10,*) tAlattice(:,:)
        end if
    end do

    close(4) 
    close(8)
    close(9)
    close(10)
    close(11)
    close(12)
    
end program sandpile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extra Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to add energy grain randomly

Subroutine addGrain(seed, L, lambda, epsilon, lattice, lambdaTimer, pU)

    real :: lambda
	integer(kind=8) :: lambdaTimer
	integer:: epsilon
	integer(kind=8), intent(in) :: L, seed
	integer(kind=8) :: pU
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer :: x, y
    integer :: sizeS=0, durationD=0 

    x =  nint(rnd48(seed) * L)                                                        ! Pick random x-coord to drive
    y =  nint(rnd48(seed) * L)                                                        ! Pick random y-coord to drive
    
    lattice(x,y) = lattice(x,y) + epsilon                           ! Bring site (x,y) to failure
    if (lattice(x,y) .gt. 1)  then                                  ! Check if driven site is above threshold, z_(x,y) >= 2
    	call relaxLattice_manna2(seed, L, lambda, sizeS, durationD, lattice, lambdaTimer, pU)          
                                                                    ! Call subroutine to topple active site
        pU = pU + 1 
    	write(8,'(I0)') sizeS
    	write(9,'(I0)') durationD
    	write(11,*) pU, lattice(L/2-2,L/2), lattice(L/2-3,L/2), lattice(L/2-4,L/2), lattice(L/2-5,L/2), lattice(L/2-6,L/2)
    	sizeS=0
    	durationD = 0
	end if

end Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to add energy grain randomly

Subroutine addGrainTransient(seed, L, lambda, epsilon, lattice, lambdaTimer, pU)

    real :: lambda
	integer(kind=8) :: lambdaTimer
	integer:: epsilon
	integer(kind=8), intent(in) :: L, seed
	integer(kind=8) :: pU
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer :: x, y
    integer :: sizeS=0, durationD=0 
    ! write(6,*) "adding a grain"

    x =  nint(rnd48(seed) * L)                                    ! Pick random x-coord to drive
    y =  nint(rnd48(seed) * L)                                    ! Pick random y-coord to drive
    
    lattice(x,y) = lattice(x,y) + epsilon                           ! Bring site (x,y) to failure
    if (lattice(x,y) .gt. 1)  then                                  ! Check if driven site is above threshold, z_(x,y) >= 2
        ! write(6,*) "Relaxing site"
    	call relaxLattice_manna2(seed, L, lambda, sizeS, durationD, lattice, lambdaTimer, pU)  ! Call subroutine to topple active site
    	sizeS=0
    	durationD = 0
	end if
	
end Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutin to topple all critical sites through manna toppling rules
! Input lambda is treated as a probability for a given site to lose grain rather than receieve it. 
! The toppling site has fixed residual height (h_r = 0)

Subroutine relaxLattice_manna1(seed, L, lambda, sizeS, durationD, lattice, lambdaTimer, pU)

	logical :: active
	real :: lambda
	integer(kind=8) :: lambdRate
	integer(kind=8) :: lambdaTimer
	integer(kind=8) :: enerLost
	
	integer(kind=8), intent(in) :: L, seed
	integer(kind=8) :: pU
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer(kind=4) :: lambda01
    integer :: j, k, chosenSite
    integer(kind=4) :: toppleThese(0:L-1, 0:L-1)
    integer :: bulkDissipationLayer(0:3)
    integer, intent(inout) :: sizeS,durationD 
        
    active = .true.
    do while (active .eqv. .true.)										! Continue relaxaing until active = .false., i.e. lattice is no longer active
	    do j = 0, L-1                                                   ! Initialize the toppling matrix
	    do k = 0, L-1                                                   ! this will be added to the lattice after all currently active sites are toppled
	        toppleThese(j,k) = 0                                        ! (hence in parallel)
	    end do                                                          !
	    end do                                                          !
	    
	    do j=0, L-1                                                     ! Run over lattice, check for sites to topple
	    do k=0, L-1
	        if (lattice(j,k) .gt. 1) then                               ! If a given site is above threshold...
	            sizeS = sizeS + 1                                       ! Add to total size, s, of avalanche (because a site is toppling)
	            
	            do while (lattice(j,k) .gt. -1*toppleThese(j,k))        ! Continuously subtract from critical site until hits h_r = 0
	                chosenSite = nint(4*rnd48(seed)-0.5)                ! Randomly choose site to topple to (randomly chosen from 0-3 for each of the four neighbors)
	                if (chosenSite .eq. 0) then
	                    if (j .ne. 0) then
                            toppleThese(j-1,k) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 1) then
                        if (j .ne. L-q) then
    	                    toppleThese(j+1,k) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 2) then
                        if (k .ne. 0) then
    	                    toppleThese(j,k-1) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 3) then
                        if (k .ne. L-1) then
    	                    toppleThese(j,k+1) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1
	                end if
	            end do
	        end if
	    end do
	    end do
	    
	    lattice = lattice + toppleThese                                 ! Add total Toppling Matrix to current lattice
	    durationD = durationD + 1										! Add to total duration, D, of avalanche 
	    
	                                                                    ! Record height of 6 site to the left of the driven site, 
	    
	    active = .false.
	    do j=0,L-1                                                      ! Check each site if over 4
	    do k=0,L-1
	    	if (lattice(j,k) .gt. 1) then
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
! Subroutin to topple all critical sites through manna toppling rules
! Input lambda is treated as a probability for a given site to lose grain rather than receieve it. 
! The toppling site has fixed change in height (Delta h = 2)

Subroutine relaxLattice_manna2(seed, L, lambda, sizeS, durationD, lattice, lambdaTimer, pU)

	logical :: active
	real :: lambda
	integer(kind=8) :: lambdRate
	integer(kind=8) :: lambdaTimer
	integer(kind=8) :: enerLost
	
	integer(kind=8), intent(in) :: L, seed
	integer(kind=8) :: pU
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer(kind=4) :: lambda01
    integer :: j, k, h, chosenSite
    integer(kind=4) :: toppleThese(0:L-1, 0:L-1)
    integer :: bulkDissipationLayer(0:3)
    integer, intent(inout) :: sizeS,durationD 
        
    active = .true.
    do while (active .eqv. .true.)										! Continue relaxaing until active = .false., i.e. lattice is no longer active
	    do j = 0, L-1                                                   ! Initialize the toppling matrix
	    do k = 0, L-1                                                   ! this will be added to the lattice after all currently active sites are toppled
	        toppleThese(j,k) = 0                                        ! (hence in parallel)
	    end do                                                          !
	    end do                                                          !
	    ! write(6,*) "starting to topple"

	    do j=0, L-1                                                     ! Run over lattice, check for sites to topple
	    do k=0, L-1
	        if (lattice(j,k) .gt. 1) then                               ! If a given site is above threshold...
	            sizeS = sizeS + 1                                       ! Add to total size, s, of avalanche (because a site is toppling)
	            ! write(6,*) "added size of avalanche"

	            do h = 1,2                                              ! Topple two grains of sand (I.e. fixed Delta h)
	                chosenSite = nint(4*rnd48(seed)-0.5)                ! Randomly choose site to topple to (randomly chosen from 0-3 for each of the four neighbors)
	                if (chosenSite .eq. 0) then
	                    if (j .ne. 0) then
                            toppleThese(j-1,k) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 1) then
                        if (j .ne. L-q) then
    	                    toppleThese(j+1,k) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 2) then
                        if (k .ne. 0) then
    	                    toppleThese(j,k-1) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 3) then
                        if (k .ne. L-1) then
    	                    toppleThese(j,k+1) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1
	                end if
	            end do
	        end if
	    end do
	    end do
	    
	    lattice = lattice + toppleThese                                 ! Add total Toppling Matrix to current lattice
	    durationD = durationD + 1										! Add to total duration, D, of avalanche 
	    
	                                                                    ! Record height of 6 site to the left of the driven site, 
	    
	    active = .false.
	    do j=0,L-1                                                      ! Check each site if over 4
	    do k=0,L-1
	    	if (lattice(j,k) .gt. 1) then
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
! Function to compute TM at time t for height:
! Height (h_j) is the observable. O_j = tAlattice_j (at site j)
 
subroutine calc_tmMetric(L, t, tAlattice)

    implicit none
    integer(kind=8), intent(in) :: L, t
    real(kind=4), intent(in) :: tAlattice(0:L-1, 0:L-1)
    real(kind=4) :: tmMetric, tAstress
    integer :: j, k
    tmMetric = 0                                                    ! tmMetric = the resulting TM metric at time t
    tAstress = 0                                                    ! tAstress = spatial averagre of the (rolling time-averaged) height, < O_j >
    
    do j=0,L-1
        do k=0,L-1
            tAstress = tAstress + tAlattice(j,k) / L / L
        end do
    end do 
    
    do j=0, L-1
        do k=0,L-1
            tmMetric = tmMetric + (tAlattice(j,k) - tAstress)**2    ! tmMetric = sum over each lattice site of (O_j + < O_j >)^2
        end do
    end do 
    
    write(12,*) t, tmMetric / L / L

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function to decide on bulk dissipation (0-no, 1-yes)

Function lambda01(lambda,seed)
    
    implicit none
    real, intent(in) :: lambda
    integer(kind=8)  :: seed
    real :: rnd48
    integer(kind=4) :: lambda01
    
    if (rnd48(seed) .le. lambda) then
        lambda01 = 0
    else
        lambda01 = 1
    end if 
    
end Function lambda01


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



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
