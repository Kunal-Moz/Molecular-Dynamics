! md_lj_vl_module.f90
! Force routine for MD, LJ atoms, using Verlet neighbour list
MODULE md_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: introduction, conclusion, allocate_arrays, deallocate_arrays
  PUBLIC :: force, hessian

  ! Public data
  INTEGER,                              PUBLIC :: n ! Number of atoms
  REAL*8,                               PUBLIC :: vir, lap, pot
  REAL*8,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: r ! Positions (3,n)
  REAL*8,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: v ! Velocities (3,n)
  REAL*8,    DIMENSION(:,:), ALLOCATABLE, PUBLIC :: f ! Forces (3,n)

  ! Public derived type
  !TYPE, PUBLIC :: potential_type ! A composite variable for interactions comprising
  !   REAL*8    :: cut ! the potential energy cut (but not shifted) at r_cut and
  !   REAL*8    :: pot ! the potential energy cut-and-shifted at r_cut and
  !   REAL*8    :: vir ! the virial and
  !   REAL*8    :: lap ! the Laplacian and
  !   LOGICAL :: ovr ! a flag indicating overlap (i.e. pot too high to use)
  ! CONTAINS
  !   PROCEDURE :: add_potential_type
  !   GENERIC   :: OPERATOR(+) => add_potential_type
  !END TYPE potential_type

CONTAINS

  !FUNCTION add_potential_type ( a, b ) RESULT (c)
  !  IMPLICIT NONE
  !  TYPE(potential_type)              :: c    ! Result is the sum of
  !  CLASS(potential_type), INTENT(in) :: a, b ! the two inputs
  !  c%cut = a%cut  +   b%cut
  !  c%pot = a%pot  +   b%pot
  !  c%vir = a%vir  +   b%vir
  !  c%lap = a%lap  +   b%lap
  !  c%ovr = a%ovr .OR. b%ovr
  !END FUNCTION add_potential_type

  SUBROUTINE introduction
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)' ) 'Lennard-Jones potential'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut-and-shifted version for dynamics'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Cut (but not shifted) version also calculated'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Diameter, sigma = 1'
    WRITE ( unit=output_unit, fmt='(a)' ) 'Well depth, epsilon = 1'

  END SUBROUTINE introduction

  SUBROUTINE conclusion
    IMPLICIT NONE

    WRITE ( unit=output_unit, fmt='(a)') 'Program ends'

  END SUBROUTINE conclusion

  SUBROUTINE allocate_arrays (dimnsn, n, box, r_cut )
    USE verlet_list_module, ONLY : initialize_list
    IMPLICIT NONE
    INTEGER, INTENT(in) :: dimnsn !Dimension of phase space
    INTEGER, INTENT(in) :: N      ! # of particles
    REAL*8, INTENT(in)  :: box   ! Simulation box length
    REAL*8, INTENT(in)  :: r_cut ! Potential cutoff distance

    !REAL*8 :: r_cut_box

    ALLOCATE ( r(dimnsn,n), v(dimnsn,n), f(dimnsn,n) )
    !r_cut_box = r_cut / box
    !IF ( r_cut_box > 0.5 ) THEN
    !   WRITE ( unit=error_unit, fmt='(a,f15.6)' ) 'r_cut/box too large ', r_cut_box
    !   STOP 'Error in allocate_arrays'
    !END IF

    CALL initialize_list ( n, r_cut, box )

  END SUBROUTINE allocate_arrays

  SUBROUTINE deallocate_arrays
    USE verlet_list_module, ONLY : finalize_list
    IMPLICIT NONE

    DEALLOCATE ( r, v, f )
    CALL finalize_list

  END SUBROUTINE deallocate_arrays
  
  SUBROUTINE force ( Box, n, dimnsn, sigma, epsi, r_cut,vir,lap,pot)
    USE verlet_list_module, ONLY : point, list, make_list
    USE config_io_module, ONLY : pbc, Check_PBC
    IMPLICIT NONE
    INTEGER                             :: N, dimnsn     ! # of particles
    REAL*8,                 INTENT(in)  :: box   ! Simulation box length
    REAL*8,                 INTENT(in)  :: r_cut ! Potential cutoff distance
    REAL*8,                 INTENT(in)  :: sigma, epsi  !LJ parameters
    REAL*8,                 INTENT(out)  :: vir,lap,pot
    !TYPE(potential_type), INTENT(out) :: total ! Composite of pot, vir, lap etc

    ! pot is the nonbonded potential energy for whole system
    ! vir is the corresponding virial
    ! lap is the corresponding Laplacian
    ! This routine also calculates forces and stores them in the array f
    ! Forces are derived from pot

    ! Forces are calculated in units where sigma = 1 and epsilon = 1
    ! Uses a Verlet list

    INTEGER              :: i, j, k, ii
    REAL*8                 :: r_cut_sq, box_sq, rij2,src2,src6,src12,boxh
    REAL*8                 :: sr2, sr6, sr12, pot_cut,sg2, x,y,z, rij,fij, vij,virij
    !REAL*8, DIMENSION(dimnsn)   :: rij, fij
    REAL*8, PARAMETER      :: sr2_ovr = 1.0d0 ! Overlap threshold (pot > 100)
    !TYPE(potential_type) :: pair
    boxh = 0.5d0*box
    r_cut_sq = r_cut ** 2
    sg2 = sigma**2
    ! Calculate potential at cutoff
    src2     = sg2 / r_cut_sq ! in sigma=1 units
    src6     = src2 ** 3
    src12    = src6 **2
    pot_cut = src12 - src6 ! Without numerical factor 4

    CALL make_list ( n, dimnsn, r, box )

    ! Initialize
    f = 0.0d0
    pot=0.0d0
    vir=0.0d0
    lap=0.0d0

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO k = point(i), point(i+1) - 1 ! Begin inner loop over neighbour atoms (if any)

          j = list(k) ! Neighbour atom index


          x = R(1,i) - R(1,j) ; y = R(2,i) - R(2,j) ; z = R(3,i) - R(3,j)

          CALL PBC(x,-boxh,boxh,box)
          CALL PBC(y,-boxh,boxh,box)
          CALL PBC(z,-boxh,boxh,box)

          rij2 = x**2 + y**2 + z**2
          rij = DSQRT(rij2)
          sr2 = sg2/rij2

          IF(rij .LT. r_cut) THEN
            sr6  = sr2**3
            sr12 = sr2**6
            vij = sr12 - sr6
            pot = pot + vij               ! LJ pair potential (cut but not shifted)
            virij = vij + sr12
            vir = vir + virij                ! LJ pair virial
            lap = lap + ( 22.0*sr12 - 5.0*sr6 ) * sr2   ! LJ pair Laplacian
            fij = virij/rij
           
            f(1,i) = f(1,i) + fij*(x/rij)
            f(2,i) = f(2,i) + fij*(y/rij)
            f(3,i) = f(3,i) + fij*(z/rij)
            f(1,j) = f(1,j) - fij*(x/rij)
            f(2,j) = f(2,j) - fij*(y/rij)
            f(3,j) = f(3,j) - fij*(z/rij)

          END IF

       END DO ! End inner loop over neighbour atoms (if any)

    END DO ! End outer loop over atoms

    ! Multiply results by numerical factors
    f   = f   * 24.0d0*epsi            ! 24*epsilon
    pot = pot * 4.0d0*epsi             ! 4*epsilon
    vir = vir * 24.0d0 / 3.0d0         ! 24*epsilon and divide virial by 3
    lap = lap * 24.0d0 * 2.0d0         ! 24*epsilon and factor 2 for ij and ji


  END SUBROUTINE force



  FUNCTION hessian ( n, box, r_cut ) RESULT ( hes )
    USE verlet_list_module, ONLY : point, list, make_list
    USE config_io_module, ONLY : pbc, Check_PBC
    IMPLICIT NONE
    INTEGER, INTENT(in) :: N     ! # of particles
    REAL*8              :: hes   ! Returns total Hessian
    REAL*8, INTENT(in)  :: box   ! Simulation box length
    REAL*8, INTENT(in)  :: r_cut ! Potential cutoff distance

    ! Calculates Hessian function (for 1/N correction to config temp)
    ! This routine is only needed in a constant-energy ensemble
    ! It is assumed that positions are in units where box = 1
    ! but the result is given in units where sigma = 1 and epsilon = 1
    ! It is assumed that forces have already been calculated 
    ! Uses a Verlet list

    INTEGER            :: i, j, k,ii
    REAL*8               :: r_cut_box, r_cut_sq, box_sq, rij_sq
    REAL*8               :: sr2, sr6, sr8, sr10, rf, ff, v1, v2
    REAL*8, DIMENSION(3) :: rij, fij

    !r_cut_box    = r_cut / box
    r_cut_sq = r_cut ** 2
    box_sq       = box ** 2

    CALL make_list ( n,3, r, box )
    
    hes = 0.0

    DO i = 1, n - 1 ! Begin outer loop over atoms

       DO k =  point(i), point(i+1) - 1 ! Begin inner loop over neighbour atoms (if any)

          j = list(k) ! Neighbour atom index

          rij(:) = r(:,i) - r(:,j)           ! Separation vector
          DO ii=1,3,1
            CALL PBC(rij(ii),0.0d0,box,box)! Periodic boundary conditions
          END DO 
          rij_sq = SUM ( rij**2 )            ! Squared separation

          IF ( rij_sq < r_cut_sq ) THEN

             !rij_sq = rij_sq * box_sq ! Now in sigma=1 units
             !rij(:) = rij(:) * box    ! Now in sigma=1 units
             fij(:) = f(:,i) - f(:,j) ! Difference in forces

             ff   = DOT_PRODUCT(fij,fij)
             rf   = DOT_PRODUCT(rij,fij)
             sr2  = 1.0 / rij_sq
             sr6  = sr2 ** 3
             sr8  = sr6 * sr2
             sr10 = sr8 * sr2
             v1   = 24.0 * ( 1.0 - 2.0 * sr6 ) * sr8
             v2   = 96.0 * ( 7.0 * sr6 - 2.0 ) * sr10
             hes  = hes + v1 * ff + v2 * rf**2

          END IF

       END DO ! End inner loop over neighbour atoms (if any)

    END DO ! End outer loop over atoms

  END FUNCTION hessian

!  SUBROUTINE vel_scale(dtsq,DEnr,scale)
!  IMPLICIT NONE
!  REAL*8 :: dtsq, DEnr, scale, Epot, DEk, scale

!    Epot = Pot/N
!    DEk = DEnr - Epot
!    scale = DSQRT(2.0d0*N*Dek*dtsq/SUM(v**2))
!  END SUBROUTINE
END MODULE md_module
