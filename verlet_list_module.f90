! verlet_list_module.f90
! Verlet list handling routines for MD simulation
MODULE verlet_list_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  
  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: initialize_list, finalize_list, make_list

  ! The initialize_list routine sets the value of r_list_factor. If you wish to read it
  ! from standard input using a NAMELIST nml_list, just uncomment the indicated statements
  ! It is assumed that all positions and displacements are divided by box
  ! r_list_box is set to r_cut_box*r_list_factor

  ! Public data
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: point ! index to neighbour list (n)
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: list  ! Verlet neighbour list (nl)

  ! Private data
  INTEGER                           :: nl         ! Size of list
  REAL*8                              :: r_list   ! List range parameter
  REAL*8                              :: r_skin   ! List skin parameter
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: r_save     ! Saved positions for list (3,n)
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: dr         ! Displacements (3,n)

CONTAINS

  SUBROUTINE initialize_list ( n, r_cut, box )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n         ! number of particles
    REAL*8,    INTENT(in) :: r_cut, box ! r_cut, box

    REAL*8    :: r_list_factor

    REAL*8, PARAMETER :: pi = 4.0d0*DATAN(1.0d0)

    ! Uncomment the following statements if you wish to read in the value of r_list_factor
!!$    INTEGER :: ioerr
!!$    NAMELIST /nml_list/ r_list_factor

    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Verlet list based on r_cut =', r_cut

    ! Sensible default for r_list_factor
    r_list_factor = 1.1d0

    ! Uncomment the following statements if you wish to read in the value of r_list_factor
!!$    READ ( unit=input_unit, nml=nml_list, iostat=ioerr ) ! namelist input
!!$    IF ( ioerr /= 0 ) THEN
!!$       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error reading namelist nml_list from standard input', ioerr
!!$       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
!!$       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
!!$       STOP 'Error in initialize_list'
!!$    END IF

    !WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Verlet list factor = ', r_list_factor
    !IF ( r_list_factor <= 1.0 ) THEN
    !   WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_list_factor must be > 1', r_list_factor
    !   STOP 'Error in initialize_list'
    !END IF
    r_list = r_cut*r_list_factor
    IF ( r_list .GT. 0.5d0*box ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_list too large', r_list
       !STOP 'Error in initialize_list'
       !r_list = 0.5d0*box
    END IF
    r_skin = r_list - r_cut
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Verlet list range = ', r_list
    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Verlet list skin  = ', r_skin

    ! Estimate list size based on density + 10 per cent
    nl = CEILING ( 1.1d0*(4.0d0*pi/3.0d0)*(r_list**3)*DBLE(n*(n-1)) / 2.0d0 )
    WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'Verlet list size = ', nl

    ALLOCATE ( r_save(3,n), dr(3,n), point(n), list(nl) )

  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list
    IMPLICIT NONE

    DEALLOCATE ( r_save, dr, point, list )

  END SUBROUTINE finalize_list

  SUBROUTINE resize_list ! reallocates list array, somewhat larger
    IMPLICIT NONE
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER                            :: nl_new

    nl_new = CEILING ( 1.25d0*DBLE(nl) )
    WRITE( unit=output_unit, fmt='(a)', advance='no' ) 'Warning: new Verlet list array size = '

    ALLOCATE ( tmp(nl_new) ) ! new size for list
    tmp(1:nl) = list(:)      ! copy elements across

    CALL MOVE_ALLOC ( tmp, list )
    nl = SIZE(list)
    WRITE( unit=error_unit, fmt='(t60,i15)') nl

  END SUBROUTINE resize_list

  SUBROUTINE make_list ( n,dimnsn, r, box )
    USE config_io_module, ONLY : PBC, Check_PBC      
    IMPLICIT NONE
    REAL*8,                  INTENT(in) :: box
    INTEGER,                 INTENT(in) :: n, dimnsn
    REAL*8,    DIMENSION(dimnsn,n), INTENT(in) :: r

    INTEGER            :: i, j, k, ii
    REAL*8               :: r_list_sq, rij_sq, dr_sq_max,boxh
    REAL*8, DIMENSION(dimnsn) :: rij

    LOGICAL, SAVE :: first_call = .TRUE.

    IF ( .NOT. first_call ) THEN
       !PRINT*, Size(dr), Size(r_save), size(r)
       dr = r - r_save                             ! Displacement since last list update
       !CALL Check_PBC(dr,N,dimnsn,box)
       dr = dr - ANINT( dr/box )*box                  ! Periodic boundaries in box=1 units
       dr_sq_max = MAXVAL ( SUM(dr**2,dim=1) )     ! Squared maximum displacement
       IF ( 4.0d0*dr_sq_max .LT. r_skin**2 ) RETURN ! No need to make list
    END IF

    first_call = .FALSE.

    k = 0
    boxh = 0.5d0*box
    r_list_sq = r_list ** 2

    DO i = 1, n - 1 ! Begin outer loop over atoms

       point(i) = k + 1

       DO j = i + 1, n ! Begin inner loop over partner atoms

          rij(:) = r(:,i) - r(:,j)
          !DO ii =1,dimnsn
          !  CALL PBC(rij(ii),-boxh,boxh,box) ;  ! Periodic boundary conditions
          !END DO
          rij(:) = rij(:) - box * ANINT ( rij(:) / box )
          rij_sq = SUM ( rij**2 )

          IF ( rij_sq .LT. r_list_sq ) THEN

             k = k + 1
             IF ( k > nl ) CALL resize_list
             list(k) = j

          END IF

       END DO ! End inner loop over partner atoms

    END DO ! End outer loop over atoms

    point(n) = k + 1

    r_save = r

  END SUBROUTINE make_list
  
END MODULE verlet_list_module
