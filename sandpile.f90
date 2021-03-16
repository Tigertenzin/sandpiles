program sandpile                                                             ! Program for running NN BTW sandpile in 2D 
    integer, parameter :: dp = selected_real_kind(15, 500)
    
    ! Initializing system-specific parameters  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    real :: lambda                                                          ! Dissipation probability
    integer(kind=4) :: lambda01
    integer :: epsilon                                                      ! Height threhsold (hc)
    integer(kind=8) :: L                                                    ! Linear system size (system is L by L)
    integer(kind=4), allocatable :: lattice(:,:)                            ! The system configuration (Dimension L x L for 2D system)
    integer(kind=4), allocatable :: latticeTriggers(:,:)                    ! Record the no. of times a site on the lattice triggers
    integer(kind=8) :: timer                                                ! Total number of avalanches to run the system for. 
    
    character(len=1024) :: filename1
    character(len=1024) :: filename2
    character(len=1024) :: filename3
    character(len=1024) :: filename4 
    character(len=1024) :: filename5


    INTEGER, PARAMETER :: LONG=SELECTED_INT_KIND(18)                        ! Necessray variables for pseudo-random number generator
    INTEGER(LONG) :: seed
    REAL rnd48, r
    
    integer(kind=8) :: k=1,j=1, m=1, x=1,y=1, t=1                           ! variables for loop indices

    
    read *,timer, L, lambda, epsilon, seed                                  ! Read in system parameters from console
    
    allocate(lattice(0:L-1, 0:L-1))                                         ! Allocate appropriate indices for lattice (linear size L: 0 -> L-1
    allocate(latticeTriggers(0:L-1, 0:L-1))                                 ! Allocate appropriate indices for time-averaged lattice (linear size L: 0 -> L-1
    
    ! Prepare files to be written to  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    write (filename1, fmt='(a,I0,a)') "sp",seed,"_s.txt"                    ! file of avalanche sizes
    open(unit=8,file=filename1)
    write(8,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(8,*) seed, L, lambda, epsilon, timer

    write (filename2, fmt='(a,I0,a)') "sp",seed,"_d.txt"                    ! file of avalanche duration
    open(unit=9,file=filename2)
    write(9,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(9,*) seed, L, lambda, epsilon, timer    
    
    write (filename3, fmt='(a,I0,a)') "sp",seed,"_a.txt"                   ! file of avalanche size w/o double counting
    open(unit=14,file=filename3)
    write(14,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(14,*) seed, L, lambda, epsilon, timer
	
    write (filename4, fmt='(a,I0,a)') "sp",seed,"_rga.txt"                    ! file of radius of gyration (rg_a)
    open(unit=10,file=filename4)
    write(10,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(10,*) seed, L, lambda, epsilon, timer    

    write (filename5, fmt='(a,I0,a)') "sp",seed,"_rgs.txt"                   ! file of radius of gyration (rg_s)
    open(unit=13,file=filename5)
    write(13,*) "seed", "SystemSize", "dissipation", "driving", "timer"
    write(13,*) seed, L, lambda, epsilon, timer
  
    
    ! Set up the 2D lattice !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j=0, L-1                                                             ! Initialize lattice as all 1’s (except the left and right) 
        do k=0,L-1
            lattice(j,k) = 0
            latticeTriggers(j,k) = 0
        End do
    end do 
    
    ! Check lattice was initialized correctly
    
	write(6,*) 0
    ! Run dynamics for set number of time steps !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do t = 1,1000000                                                      ! Run system for transience, 1,000,000 avalanches
        call addGrainTransient(seed, L, lambda, epsilon, lattice, latticeTriggers)

    end do

    do t = 1,timer
        call addGrain(seed, L, lambda, epsilon, lattice, latticeTriggers)
    end do
    
    close(8)
    close(9)
    close(10)
    close(13)
    close(14)
    
end program sandpile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extra Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to add energy grain randomly

Subroutine addGrain(seed, L, lambda, epsilon, lattice, latticeTriggers)

	logical :: unstable
    real :: lambda
	integer:: epsilon
	integer(kind=8), intent(in) :: L, seed
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer(kind=4), intent(inout) :: latticeTriggers(0:L-1, 0:L-1)
    integer :: x, y
    integer :: sizeS, durationD, sizeA
	
    sizeS = 0
    durationD = 0
    sizeA = 0
    do j=0,L-1                                                  ! Calculate and record the cluster size w/o counting (from Triggers array)
    do k=0,L-1
        latticeTriggers(j,k) = 0
    end do
    end do 
	
	unstable = .false.
	do while (unstable .eqv. .false.)
		
	    x =  nint(rnd48(seed)*(L-1))                                    ! Pick random x-coord to drive
	    y =  nint(rnd48(seed)*(L-1))                                    ! Pick random y-coord to drive
    
	    lattice(x,y) = lattice(x,y) + 1                           		! Bring site (x,y) to failure
	    if (lattice(x,y) .ge. epsiilon)  then                           ! Check if driven site is above threshold, z_(x,y) > 4
	    	call relaxLattice_stoch_fD(seed, L, lambda, epsilon, sizeS, durationD, lattice, latticeTriggers)          
	                                                                    ! Call subroutine to topple active site
        
	        do j=0,L-1                                                  ! Calculate and record the cluster size w/o counting (from Triggers array)
	        do k=0,L-1
	            if (latticeTriggers(j,k) .ne. 0) then                   ! Record a site failing if it didn't trigger
	                sizeA = sizeA + 1
	            end if
	        end do
	        end do 
			
			call calc_rg(L, sizeS, latticeTriggers)                     ! Cacluate and record the radius of gyration. 
			call calc_rga(L, sizeA, latticeTriggers)                     ! Cacluate and record the radius of gyration. 
        
	    	write(8,'(I0)') sizeS
	    	write(9,'(I0)') durationD
	    	write(14,'(I0)') sizeA
    	
	    	sizeS=0
	    	durationD = 0
	    	sizeA=0
			
			unstable = .true.
		end if
	
	end do

end Subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to add energy grain randomly

Subroutine addGrainTransient(seed, L, lambda, epsilon, lattice, latticeTriggers)

	logical :: unstable
    real :: lambda
	integer:: epsilon
	integer(kind=8), intent(in) :: L, seed
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer(kind=4), intent(inout) :: latticeTriggers(0:L-1, 0:L-1)
    integer :: x, y
    integer :: sizeS, durationD, sizeA
	
    sizeS = 0
    durationD = 0
    sizeA = 0
	
	unstable = .false.
	do while (unstable .eqv. .false.)
		
	    x =  nint(rnd48(seed)*(L-1))                                    ! Pick random x-coord to drive
	    y =  nint(rnd48(seed)*(L-1))                                    ! Pick random y-coord to drive
    
	    lattice(x,y) = lattice(x,y) + 1                           		! Bring site (x,y) to failure
	    if (lattice(x,y) .ge. epsilon)  then                            ! Check if driven site is above threshold, z_(x,y) > 4
	    	call relaxLattice_stoch_fD(seed, L, lambda, epsilon, sizeS, durationD, lattice, latticeTriggers)          
	                                                                    ! Call subroutine to topple active site
			unstable = .true.
		end if
	
	end do
	
end Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutin to topple all critical sites stochastically
! Input lambda is treated as a probability for a given site to lose grain rather than receieve it. 
! Toppling is fixed delta h, i.e. fixed change in height for tpopled site

Subroutine relaxLattice_stoch_fD(seed, L, lambda, epsilon, sizeS, durationD, lattice, latticeTriggers)

	logical :: active
	real :: lambda
	integer:: epsilon
	integer(kind=8), intent(in) :: L, seed
    integer(kind=4), intent(inout) :: lattice(0:L-1, 0:L-1)
    integer(kind=4), intent(inout) :: latticeTriggers(0:L-1, 0:L-1)
    integer(kind=4) :: lambda01
    integer :: j, k, h
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
            if (lattice(j,k) .ge. epsilon) then                               ! If a given site is above threshold...            
                latticeTriggers(j,k) = latticeTriggers(j,k) + 1         ! Record that site (j,k) as triggered once
	            sizeS = sizeS + 1                                       ! Add to total size, s, of avalanche (because site topples)
	            
				do h = 1,epsilon
					chosenSite = nint(4*rnd48(seed)-0.5)				! Randomly choose site to topple to (randomly chosen from 0-3 for each of the four neighbors)
					
	                if (chosenSite .eq. 0) then
	                    if (j .ne. 0) then
                            toppleThese(j-1,k) = toppleThese(j-1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 1) then
                        if (j .ne. L-q) then
    	                    toppleThese(j+1,k) = toppleThese(j+1,k) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 2) then
                        if (k .ne. 0) then
    	                    toppleThese(j,k-1) = toppleThese(j,k-1) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1

	                else if (chosenSite .eq. 3) then
                        if (k .ne. L-1) then
    	                    toppleThese(j,k+1) = toppleThese(j,k+1) + 1 * lambda01(lambda,seed)
                        end if
	                    toppleThese(j,k)   = toppleThese(j,k)   - 1
	                end if
					
				end do 
	            
	        end if
	    end do
	    end do
	    
	    lattice = lattice + toppleThese                                 ! Add total Toppling Matrix to current lattice
	    durationD = durationD + 1										! Add to total duration, D, of avalanche 
	    
	    active = .false.
	    do j=0,L-1                                                      ! Check each site if over 4
	    do k=0,L-1
	    	if (lattice(j,k) .ge. epsilon) then
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
 
subroutine calc_rg(L, sizeS, latticeTriggers)

    implicit none
    integer(kind=8), intent(in) :: L
    integer(kind=4), intent(inout) :: latticeTriggers(0:L-1, 0:L-1)
    real(kind=4) :: comX, comY, rg
    integer :: j, k, sizeS
    comX = 0                                                         ! x-coordcom
    comY = 0                                                         ! y-coord com
    rg = 0                                                          ! rg

    
    do j=0,L-1
        do k=0,L-1
            comX = comX + j * real(latticeTriggers(j,k), kind=4) / real(sizeS, kind=4)
            comY = comY + k * real(latticeTriggers(j,k), kind=4) / real(sizeS, kind=4)
        end do
    end do 
    
    do j=0, L-1
        do k=0,L-1
            rg = rg + real(latticeTriggers(j,k), kind=4) * (real(j, kind=4) - comX)**2
            rg = rg + real(latticeTriggers(j,k), kind=4) * (real(k, kind=4) - comY)**2
        end do
    end do 
    
    write(13,*) rg/real(sizeS, kind=4)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
subroutine calc_rga(L, sizeA, latticeTriggers)

    implicit none
    integer(kind=8), intent(in) :: L
    integer(kind=4), intent(inout) :: latticeTriggers(0:L-1, 0:L-1)
	integer(kind=4) :: latticeTriggeredOnce(0:L-1, 0:L-1)
    real(kind=4) :: comX, comY, rg
    integer :: j, k, sizeA
    comX = 0                                                         ! x-coordcom
    comY = 0                                                         ! y-coord com
    rg = 0                                                          ! rg
	
    do j=0,L-1
        do k=0,L-1
            if (latticeTriggers(j,k) .ne. 0) then                   ! Record a site failing if it didn't trigger
                latticeTriggeredOnce(j,k) = 1
            end if
            if (latticeTriggers(j,k) .eq. 0) then                   ! Record a site failing if it didn't trigger
                latticeTriggeredOnce(j,k) = 0 
            end if
        end do
    end do 
    
    do j=0,L-1
        do k=0,L-1
            comX = comX + j * real(latticeTriggeredOnce(j,k), kind=4) / real(sizeA, kind=4)
            comY = comY + k * real(latticeTriggeredOnce(j,k), kind=4) / real(sizeA, kind=4)
        end do
    end do 
    
    do j=0, L-1
        do k=0,L-1
            rg = rg + real(latticeTriggeredOnce(j,k), kind=4) * (real(j, kind=4) - comX)**2
            rg = rg + real(latticeTriggeredOnce(j,k), kind=4) * (real(k, kind=4) - comY)**2
        end do
    end do 
    
    write(10,*) rg/real(sizeA, kind=4)

end subroutine

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function to decide on bulk dissipation (0-no, 1-yes)

Function lambda01(lambda,seed)
    
    implicit none
    real, intent(in) :: lambda
    integer(kind=8)  :: seed
    real :: rnd48
    integer(kind=4) :: lambda01

    lambda01 = 0                                                            ! Assume dissipation    
    if (rnd48(seed) .ge. lambda) then
        lambda01 = 1                                                        ! If random number exceeds lambda, then don't dissipate (when lambda=0, then never dissipating)
    end if 
    
end Function lambda01