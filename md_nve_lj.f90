! md_nve_lj.f90
! Molecular dynamics, NVE ensemble
PROGRAM md_nve_lj

  ! Takes in a configuration of atoms (positions, velocities)
  ! Cubic periodic boundary conditions
  ! Conducts molecular dynamics using Velocity Verlet algorithm
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! However, input configuration, output configuration, most calculations, and all results 
  ! are given in simulation/Natural units defined by the model
  ! Interaction between the particles are under Lennard Jones Potential.
  ! Note : Values of sigma, epislon , kb and mass are all set to 1 in the simulation.

  ! Despite the program name, there is nothing here specific to Lennard-Jones
  ! The model is defined in md_module

  USE, INTRINSIC :: iso_fortran_env,    ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  USE               config_io_module,   ONLY : read_randm_particles, write_cnf_atoms, bcc_latt, PBC, check_pbc, &
       &                                     fcc_latt
  USE               averages_module,    ONLY : run_begin, run_end, blk_begin, blk_end, blk_add
  USE               md_module,          ONLY : introduction, conclusion, allocate_arrays, deallocate_arrays, &
       &                                     force, r, v, f, vir,lap,pot!potential_type
  !USE               verlet_list_module, ONLY : 

  IMPLICIT NONE

  ! Most important variables
  REAL*8 :: Density ! Reduced density of the material rho' = rho*(sigma**3/m)
  REAL*8 :: Box     ! Box length
  REAL*8 :: dt      ! Time step
  REAL*8 :: r_cut   ! Potential cutoff distance
  REAL*8 :: amu, mass, sigma, epsilon

  REAL*8 :: dtsq, temp, KE_init, Aheat, scale
  REAL*8 :: Epot, DEnr, DEk
  INTEGER            :: dimnsn, n, seed, blk, stp, nstep, nblock, ioerr, k
  LOGICAL            :: yes_vel, yes_pos
  REAL, DIMENSION(3) :: vcm

  CHARACTER(len=4), PARAMETER :: cnf_prefix = 'cnf_'
  CHARACTER(len=3), PARAMETER :: inp_tag    = 'inp'
  CHARACTER(len=3), PARAMETER :: out_tag    = 'out'
  CHARACTER(len=3)            :: sav_tag    = 'sav' ! May be overwritten with block number
  CHARACTER(len=3),PARAMETER  :: err_tag    = 'err' ! Tag for overlapping configuration
  CHARACTER(len=4), PARAMETER :: filetype   = '.dat' ! File type

  NAMELIST /nml/ nblock, nstep, r_cut, dt

  WRITE ( unit=output_unit, fmt='(a)' ) 'md_nve_lj'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE ensemble'
  WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass=1 throughout'
  CALL introduction

  ! Set sensible default run parameters for testing
  nblock = 10
  nstep  = 600000
  r_cut  = 2.25d0
  dt     = 0.0005d0
  ioerr  = 0
  dimnsn = 3
  dtsq = dt**2

  !!!! Not required for the simulation but may be used in certain cases if we want in Real units.
  amu = 1.67E-24
  sigma = 3.4E-8
  mass = 39.95
  !kb = 1.38E-23
  !!!!! Write out run parameters
  PRINT*, 'Number of blocks         ', nblock
  PRINT*, 'Number of steps per block', nstep
  PRINT*, 'Potential cutoff distance', r_cut
  PRINT*, 'Time step                ', dt
  PRINT*, 'Degrees of freedom       ',dimnsn
  
  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
 ! READ ( unit=input_unit, nml=nml, iostat=ioerr )
 ! IF ( ioerr /= 0 ) THEN
 !    WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
 !    IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
 !    IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
 !    STOP 'Error in md_nve_lj'
 ! END IF

  ! Number of particles (this can be input by the user if the following lines are uncommented)
  !PRINT*, "Input number of atoms : "
  !Read*, N
  N = 256

  !!!!! Flags for random velocities and position. Seed can be changed.
  yes_pos = .FALSE.
  yes_vel = .TRUE.
  seed = 4848462

  Density = 0.8094d0
  PRINT*, "Temperature seed for velocities :"
  READ*, Temp

  !!!!! Box dimension
  Box = (DBLE(n)/density)**(1.0d0/3.0d0)
  
  !!!!! Read in initial configuration and allocate necessary arrays

  WRITE ( unit=output_unit, fmt='(a,t50,i15)'   ) 'Number of particles :',   n
  WRITE ( unit=output_unit, fmt='(a,t50,f15.10)' ) 'Actual Density of the material ( in g/cm**3):', (Density*mass*amu)/sigma**3
  WRITE ( unit=output_unit, fmt='(a,t50,f15.10)' ) 'Intial temperature of the system :', Temp
  CALL allocate_arrays ( dimnsn, N, Box, r_cut )
  CALL read_randm_particles(N, box, dimnsn, R, V, seed, temp, yes_pos, yes_vel)

  CALL fcc_latt(R,N,Dimnsn,Box)                             ! Generate a fcc lattice
  WRITE( unit=output_unit, fmt='(a,t50,f15.10)'   ) 'Simulation Box length :' , Box
  CALL Check_PBC(R,N,Dimnsn,box)                        ! Periodic boundaries
  vcm(:) = SUM ( v(:,:), dim=2 ) / REAL(n)                  ! Centre-of mass velocity
  v(:,:) = v(:,:) - SPREAD ( vcm(:), dim = 2, ncopies = n ) ! Set COM velocity to zero

  !!!! Initial forces, potential, etc plus overlap check
  CALL force ( Box, n, dimnsn, 1.0d0, 1.0d0, r_cut, vir, lap, pot)

  PRINT*, "Initial Potential Energy :", Pot/Dble(N)
  PRINT*, "Initial Kinetic Energy   :", 0.5d0*SUM(v**2)/DBLE(N)

  !!!! Initialize arrays for averaging and write column headings
  CALL run_begin ( calc_variables() )
  CALL write_cnf_atoms ( cnf_prefix//inp_tag//filetype, dimnsn, n, Box, r, v ) ! Write out initial configuration

  OPEN(UNIT = 12, FILE = "Run_data.dat")
  DO blk = 1, nblock ! Begin loop over blocks

     CALL blk_begin

     DO stp = 1, nstep ! Begin loop over steps
        !!!! Velocity Verlet algorithm
        v(:,:) = v(:,:) + 0.5d0 * dt * f(:,:) ! Kick half-step

        r(:,:) = r(:,:) + dt * v(:,:) ! Drift step

        CALL Check_PBC(R,N,dimnsn,box) ! Check periodic boundary conditions

        CALL force ( Box, n, dimnsn, 1.0d0, 1.0d0, r_cut, vir, lap, pot) ! Force evaluation

        v(:,:) = v(:,:) + 0.5d0 * dt * f(:,:) ! Kick half-step

        ! Calculate and accumulate variables for this step
        CALL blk_add ( calc_variables() )

     END DO ! End loop over steps

     CALL blk_end ( blk )                                           ! Output block averages
     
     IF ( nblock < 1000 ) WRITE(sav_tag,'(i3.3)') blk               ! Number configuration by block
     CALL write_cnf_atoms ( cnf_prefix//sav_tag//filetype, dimnsn, n, box, r, v ) ! Save configuration

  END DO ! End loop over blocks

  CALL run_end ( calc_variables() ) ! Output run averages

  CALL force ( Box, n, dimnsn, 1.0d0, 1.0d0, r_cut,vir,lap,pot)

  CALL write_cnf_atoms ( cnf_prefix//out_tag//filetype, dimnsn, n, Box, r, v ) ! Write out final configuration
   
  CALL deallocate_arrays
  CALL conclusion

CONTAINS

  FUNCTION calc_variables ( ) RESULT ( variables )
    USE lrc_module,      ONLY : potential_lrc, pressure_lrc
    USE md_module,       ONLY : hessian
    USE averages_module, ONLY : variable_type, msd, cke
    IMPLICIT NONE
    TYPE(variable_type), DIMENSION(8) :: variables ! The 8 variables listed below

    ! This routine calculates all variables of interest and (optionally) writes them out
    ! They are collected together in the variables array, for use in the main program

    TYPE(variable_type) :: e_s, p_s, e_f, p_f, t_k, t_c, c_s, conserved_msd
    REAL*8                :: vol, rho, fsq, kin, eng, hes, tmp

    ! Preliminary calculations
    vol = Box**3                  ! Volume
    rho = DBLE(N)/vol          ! Density
    kin = 0.5d0*SUM(v**2)           ! Kinetic energy
    tmp = 2.0d0*kin/DBLE(3*n-3) ! Remove three degrees of freedom for momentum conservation
    fsq = SUM( f**2 )            ! Total squared force
    hes = hessian(n,box,r_cut)      ! Total Hessian
    eng = kin + pot         ! Total energy for simulated system

    ! Variables of interest, of type variable_type, containing three components:
    !   %val: the instantaneous value
    !   %nam: used for headings
    !   %method: indicating averaging method
    ! If not set below, %method adopts its default value of avg
    ! The %nam and some other components need only be defined once, at the start of the program,
    ! but for clarity and readability we assign all the values together below

    ! Internal energy (cut-and-shifted) per atom
    ! Total KE plus total cut-and-shifted PE divided by N
    e_s = variable_type ( nam = 'E/N cut&shifted', val = eng/REAL(n) )

    ! Internal energy (full, including LRC) per atom
    ! LRC plus total KE plus total cut (but not shifted) PE divided by N
    e_f = variable_type ( nam = 'E/N full', val = potential_lrc(rho,r_cut) + (kin+pot)/REAL(n) )

    ! Pressure (cut-and-shifted)
    ! Ideal gas contribution plus total virial divided by V
    p_s = variable_type ( nam = 'P cut&shifted', val = rho*tmp + vir/vol )

    ! Pressure (full, including LRC)
    ! LRC plus ideal gas contribution plus total virial divided by V
    p_f = variable_type ( nam = 'P full', val = pressure_lrc(rho,r_cut) + rho*tmp + vir/vol )

    ! Kinetic temperature
    t_k = variable_type ( nam = 'T kinetic', val = tmp )

    ! Configurational temperature
    ! Total squared force divided by total Laplacian with small Hessian correction
    t_c = variable_type ( nam = 'T config', val = fsq/(lap-(2.0*hes/fsq)) )

    ! MSD of kinetic energy, intensive
    ! Use special method to convert to Cv/N
    c_s = variable_type ( nam = 'Cv/Ncut&s hifted', val = kin/SQRT(REAL(n)), method = cke, instant = .FALSE. )

    ! Mean-squared deviation of conserved energy per atom
    conserved_msd = variable_type ( nam = 'Conserved MSD', val = eng/REAL(n), &
         &                          method = msd, e_format = .TRUE., instant = .FALSE. )

    ! Collect together for averaging
    variables = [ e_s, p_s, e_f, p_f, t_k, t_c, c_s, conserved_msd ]

  END FUNCTION calc_variables

END PROGRAM md_nve_lj

