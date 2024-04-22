program ising_2d
    implicit none

    ! Parameters
    integer :: N, L               ! Size of lattice in x and y direction
    integer, parameter :: steps = 5000  
    integer :: i, j, k, j1, k1, kp, km, jp, jm, t, m
    real :: ran, temp, energy, energy_sq, magnetization, magnetization_sq, neighbors, w, dE
    real ::  energy_ave, energy_sq_ave, magnetization_ave, &
    magnetization_sq_ave, magnetization_abs_ave
    real :: specific_heat, magnetic_susceptibility, magnetic_susceptibility_abs
    real :: max_susceptibility, max_heat, max_heat_temp, max_susceptibility_temp
    real, dimension(:,:), allocatable :: spins


do L = 0, 100, 10
N = L*L
allocate(spins(L, L))

! Initialize lattice with all spins pointed up like for low temperature system
    do j = 1, L
        do k = 1, L
                spins(j, k) = 1     ! Spin up
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
               
max_susceptibility = 0.

do t = 0, 1000
        temp = 2 + t * 0.001
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

        ! Main Monte Carlo simulation


energy_ave = 0.
energy_sq_ave = 0.
magnetization_ave = 0.
magnetization_abs_ave = 0.
magnetization_sq_ave = 0.

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
        energy_ave = energy_ave + energy
        energy_sq_ave = energy_sq_ave + energy*energy
        magnetization_ave = magnetization_ave + magnetization
        magnetization_abs_ave = magnetization_abs_ave + abs(magnetization)
        magnetization_sq_ave = magnetization_sq_ave + magnetization * magnetization

    end do
    specific_heat = (energy_sq_ave/(steps - (steps/10)) - &
    (energy_ave/(steps - (steps/10)))**2)/(N*temp**2)
    magnetic_susceptibility = (magnetization_sq_ave/(steps - (steps/10)) - &
    (magnetization_ave/(steps - (steps/10)))**2) / (temp*N)
     magnetic_susceptibility_abs = (magnetization_sq_ave/(steps - (steps/10)) - &
    (magnetization_abs_ave/(steps - (steps/10)))**2) / (temp*N)
     
    
  !  print *, L, T, magnetic_susceptibility, max_susceptibility
    if (magnetic_susceptibility > max_susceptibility) then
                max_susceptibility = magnetic_susceptibility
                max_susceptibility_temp = temp
    endif
            
             if (specific_heat > max_heat) then
                max_heat = specific_heat
                max_heat_temp = temp
            endif
                
enddo ! T loop

            print*, max_heat_temp, max_susceptibility_temp, L

 deallocate(spins)
 
enddo ! lattice size loop

end program ising_2d


