program ising_2d
    implicit none

    ! Parameters
    integer :: N, L, steps                 ! Size of lattice in x and y direction
    integer :: i, j, k, j1, k1, kp, km, jp, jm, t, m
    real :: ran, temp, energy, magnetization, neighbors, w, dE
    character(len=50) :: arg1, arg2  ! Adjust length according to expected input size
    integer, dimension(:,:), allocatable :: spins

    ! Check if there are exactly two command-line arguments
  if (command_argument_count() /= 2) then
      print*, "Usage: ./ising-2d (Length of Matrix) (MCS)"
      stop
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  
  read(arg1, *) L
  read(arg2, *) steps
  
  N = L * L
  allocate(spins(L, L))

!print*, "Size of the matrix is:", L
!print*, "Number of MCS is:", steps
!print*, "Total number of lattice sites is:", N

 !   open(1, file='ising.dat')


    ! Initialize lattice with random spins
    do j = 1, L
        do k = 1, L
            ran = rand()
            if (ran < 0.5) then
                spins(j, k) = 1     ! Spin up
            else
                spins(j, k) = 1    ! Spin down
            end if
        end do
    end do
    
   ! Initial Magnetization
        magnetization = 0.
        do  j = 1,L
            do  k = 1,L
            magnetization = magnetization + spins(j,k)
            enddo
        enddo
        magnetization = magnetization
        
  !     print*, "magnetization when all spins are up", magnetization
        
        ! Initial Energy
        energy = 0.
	   do  j = 1,L
            do  k = 1,L
                    jp = mod(j , L) + 1
                    jm = mod(j - 2 + L, L) + 1
                    kp = mod(k , L) + 1
                    km = mod(k - 2 + L, L) + 1
		energy = energy - spins(j,k)*(spins(jp, k) + spins(jm, k) + &
                                spins(j, kp) + spins(j, km))
                end do
                end do
                energy = energy/2
    !      print*, "energy when all spins are up", energy
 ! Print the spin matrix
do t = 0, 70
        temp = 0.50 + t * 0.05
do m =  1, steps/10
        ! Equilibration phase
        do i = 1, N
            ran = rand()
	    j = 1 + int(ran * L)
            ran = rand()
            k = 1 + int(ran * L)
                    jp = mod(j , L) + 1
                    jm = mod(j - 2 + L, L) + 1
                    kp = mod(k , L) + 1
                    km = mod(k - 2 + L, L) + 1

                    neighbors = spins(jp, k) + spins(jm, k) + &
                                spins(j, kp) + spins(j, km)

                    dE = 2. * neighbors * spins(j, k)

                    if (dE < 0) then
                        spins(j, k) = -spins(j, k)
                        magnetization =  magnetization + 2*spins(j, k)
                        energy = energy + dE
                    else
                        ran = rand()
                        w = exp(-(dE/temp))
                        if (ran <= w) then
                            spins(j, k) = -spins(j, k)
                            magnetization =  magnetization + 2*spins(j, k)
                            energy = energy + dE
                        endif
                    endif
                  !  print*, energy, magnetization
                   ! energy = energy
                end do
          end do
          
          if (temp == 0.177) then   !can vary parameter of temp to see when changes occur
    ! Write spins matrix to file
    open(unit=10, file='spins_initial_matrix.txt', status='replace', action='write')

    do j = 1, L
        write(10, '(10I4)') (spins(j, k), k = 1, L)
    end do

    close(10)
endif


        ! Main Monte Carlo simulation

do m = 1, steps - (steps/10)
        do i = 1, N
            ran = rand()
	    j = 1 + int(ran * L)
            ran = rand()
            k = 1 + int(ran * L)
            
                    jp = mod(j , L) + 1
                    jm = mod(j - 2 + L, L) + 1
                    kp = mod(k , L) + 1
                    km = mod(k - 2 + L, L) + 1

                    neighbors = spins(jp, k) + spins(jm, k) + &
                                spins(j, kp) + spins(j, km)

                    dE = 2. * neighbors * spins(j, k)

                    if (dE < 0) then
                        spins(j, k) = -spins(j, k)
                        magnetization =  magnetization + 2*spins(j, k)
                        energy = energy + dE
                    else
                        ran = rand()
                        w = exp(-(dE/temp))
                        if (ran <= w) then
                            spins(j, k) = -spins(j, k)
                            magnetization =  magnetization + 2*spins(j, k)
                            energy = energy + dE
                        endif
                    endif

                    
        end do
    end do
enddo

 open(unit=11, file='spins_final_matrix.txt', status='replace', action='write')

    do j = 1, L
        write(11, '(10I4)') (spins(j, k), k = 1, L)
    end do

    close(11)

 !   close(1)
     deallocate(spins)

end program ising_2d

