! config_io_module.f90
! Routines for atomic/molecular configuration I/O
MODULE config_io_module

  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit, iostat_end, iostat_eor, output_unit

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: read_cnf_atoms, write_cnf_atoms, read_cnf_mols, write_cnf_mols, read_randm_particles, bcc_latt, pbc, Check_pbc, FCC_Latt

CONTAINS

   SUBROUTINE read_cnf_atoms ( filename, n, r, v ) ! Read in atomic configuration
     IMPLICIT NONE
     CHARACTER(len=*),               INTENT(in)    :: filename ! Supplied filename
     INTEGER,                        INTENT(in) :: n        ! Number of atoms
     !REAL,                           INTENT(out)   :: box      ! Simulation box length (assumed cubic)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r        ! Atomic positions (3,n)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: v        ! Atomic velocities (3,n)
 
     ! The first call of this routine is just to get n and box, used to allocate arrays
     ! The second call attempts to read in the atomic positions and optionally velocities
 
     INTEGER :: cnf_unit, ioerr, i
 
     ! Open given filename, will terminate on any errors
 
     OPEN ( newunit=cnf_unit, file=filename, status='old', action='read', iostat=ioerr )
     IF ( ioerr /= 0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
        STOP 'Error in read_cnf_atoms'
     END IF
 
     ! Read number of atoms from first record
     
     !READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) n
     !IF ( ioerr /= 0 ) THEN
     !   WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading n from ', filename, ioerr
     !   IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     !   IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     !   STOP 'Error in read_cnf_atoms'
     !END IF
 
     ! Read box length from second record
     
     !READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) box
     !IF ( ioerr /= 0 ) THEN
     !   WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading box from ', filename, ioerr
     !   IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
     !   IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
     !   STOP 'Error in read_cnf_atoms'
     !END IF
 
     ! Expected format is one record per atom containing either r(:,i), v(:,i) or just r(:,i)
 
     IF ( PRESENT ( r ) ) THEN
 
        IF ( n > SIZE ( r, dim=2 ) ) THEN
           WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
           STOP 'Error in read_cnf_atoms'
        END IF
 
        IF ( PRESENT ( v ) ) THEN
 
           IF ( n > SIZE ( v, dim=2 ) ) THEN
              WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
              STOP 'Error in read_cnf_atoms'
           END IF
 
           ! Read positions, velocities
           DO i = 1, n
              READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), v(:,i)
              IF ( ioerr /= 0 ) THEN
                 WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, v from ', filename, ioerr
                 IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                 IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                 STOP 'Error in read_cnf_atoms'
              END IF
           END DO
 
        ELSE
 
           ! Read positions
           DO i = 1, n
              READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i)
              IF ( ioerr /= 0 ) THEN
                 WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r from ', filename, ioerr
                 IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                 IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                 STOP 'Error in read_cnf_atoms'
              END IF
           END DO
 
        END IF
 
     END IF
 
     CLOSE ( unit=cnf_unit )
 
   END SUBROUTINE read_cnf_atoms
 
   SUBROUTINE write_cnf_atoms ( filename, dimnsn, n, box, r, v ) ! Write out atomic configuration
     IMPLICIT NONE
     CHARACTER(len=*),                 INTENT(in) :: filename ! Supplied filename
     INTEGER,                          INTENT(in) :: n        ! Number of atoms
     INTEGER,                          INTENT(in) :: dimnsn   ! Dimensions
     REAL*8,                           INTENT(in) :: box      ! Simulation box length (assumed cubic)
     REAL*8, DIMENSION(:,:),           INTENT(in) :: r        ! Atomic positions (3,n)
     REAL*8, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v        ! Atomic velocities
 
     INTEGER :: cnf_unit, ioerr, i
 
     ! Open given filename, replacing it if it already exists, will terminate on any errors
     OPEN ( newunit=cnf_unit, file=filename, status='replace', iostat=ioerr )
     IF ( ioerr /= 0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
        STOP 'Error in write_cnf_atoms'
     END IF
 
     ! Write number of atoms to first record, box length to second record
     !WRITE ( unit=cnf_unit, fmt='(i15)'  ) n
     !WRITE ( unit=cnf_unit, fmt='(f15.8)') box
 
     IF ( n > SIZE ( r, dim=2 ) ) THEN
        WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
        STOP 'Error in write_cnf_atoms'
     END IF
 
     IF ( PRESENT ( v ) ) THEN
 
        IF ( n > SIZE ( v, dim=2 ) ) THEN
           WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
           STOP 'Error in write_cnf_atoms'
        END IF
 
        ! Write positions, velocities
        DO i = 1, n
           WRITE ( unit=cnf_unit, fmt='(*(f15.10))' ) r(:,i), v(:,i)
        END DO
 
     ELSE
 
        ! Write positions
        DO i = 1, n
           WRITE ( unit=cnf_unit, fmt='(*(f15.10))' ) r(:,i)
        END DO
 
     END IF
 
     CLOSE ( unit=cnf_unit )
 
   END SUBROUTINE write_cnf_atoms
 
   SUBROUTINE read_cnf_mols ( filename, n, box, r, e, v, w ) ! Read in molecular configuration
     IMPLICIT NONE
     CHARACTER(len=*),               INTENT(in)    :: filename ! Supplied filename
     INTEGER,                        INTENT(inout) :: n        ! Number of molecules
     REAL,                           INTENT(out)   :: box      ! Simulation box length (assumed cubic)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r        ! Molecular positions (3,n)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: e        ! Orientation vectors (3,n) or quaternions (4,n)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: v        ! Molecular velocities (3,n)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: w        ! Angular velocities/momenta (3,n)
 
     ! The first call of this routine is just to get n and box, used to allocate arrays
     ! The second call attempts to read in the molecular positions, orientations
     ! and optionally velocities, angular velocities/momenta
 
     INTEGER :: cnf_unit, ioerr, i
 
     ! Open given filename, will terminate on any errors
 
     OPEN ( newunit=cnf_unit, file=filename, status='old', action='read', iostat=ioerr )
     IF ( ioerr /= 0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
        STOP 'Error opening file in read_cnf_mols'
     END IF
 
     ! Read number of molecules from first record
     
     READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) n
     IF ( ioerr /= 0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading n from ', filename, ioerr
        IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
        IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
        STOP 'Error in read_cnf_mols'
     END IF
 
     ! Read box length from second record
     
     READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) box
     IF ( ioerr /= 0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading box from ', filename, ioerr
        IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
        IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
        STOP 'Error in read_cnf_mols'
     END IF
 
     ! Expected format is one line per atom containing either r(:,i), e(:,i), v(:,i), w(:,i)  or just r(:,i), e(:,i)
     ! The first dimension of the e array can be 3 elements (vector) or 4 elements (quaternion)
 
     IF ( PRESENT ( r ) ) THEN
 
        IF ( .NOT. PRESENT ( e )    ) THEN
           WRITE ( unit=error_unit, fmt='(a,a,i15)') 'r and e arguments must be present together'
           STOP 'Error in read_cnf_mols'
        END IF
        IF ( n > SIZE ( r, dim=2 ) ) THEN
           WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
           STOP 'Error in read_cnf_mols'
        END IF
        IF ( n > SIZE ( e, dim=2 ) ) THEN
           WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of e ', n, SIZE ( e, dim=2 )
           STOP 'Error in read_cnf_mols'
        END IF
 
        IF ( PRESENT ( v ) ) THEN
 
           IF ( .NOT. PRESENT ( w )    ) THEN
              WRITE ( unit=error_unit, fmt='(a,a,i15)') 'v and w arguments must be present together'
              STOP 'Error in read_cnf_mols'
           END IF
           IF ( n > SIZE ( v, dim=2 ) ) THEN
              WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
              STOP 'Error in read_cnf_mols'
           END IF
           IF ( n > SIZE ( w, dim=2 ) ) THEN
              WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of w ', n, SIZE ( w, dim=2 )
              STOP 'Error in read_cnf_mols'
           END IF
 
           ! Read positions, orientation vectors or quaternions, velocities, angular velocities/momenta
           DO i = 1, n
              READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), e(:,i), v(:,i), w(:,i)
              IF ( ioerr /= 0 ) THEN
                 WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, e, v, w from ', filename, ioerr
                 IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                 IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                 STOP 'Error in read_cnf_mols'
              END IF
           END DO
 
        ELSE
 
           ! Read positions, orientation vectors or quaternions
           DO i = 1, n
              READ ( unit=cnf_unit, fmt=*, iostat=ioerr ) r(:,i), e(:,i)
              IF ( ioerr /= 0 ) THEN
                 WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error reading r, e from ', filename, ioerr
                 IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
                 IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
                 STOP 'Error in read_cnf_mols'
              END IF
           END DO
 
        END IF
 
     END IF
 
     CLOSE ( unit=cnf_unit )
 
   END SUBROUTINE read_cnf_mols
 
   SUBROUTINE write_cnf_mols ( filename, n, box, r, e, v, w ) ! Write out molecular configuration
     IMPLICIT NONE
     CHARACTER(len=*),               INTENT(in) :: filename ! Supplied filename
     INTEGER,                        INTENT(in) :: n        ! Number of molecules
     REAL,                           INTENT(in) :: box      ! Simulation box length (assumed cubic)
     REAL, DIMENSION(:,:),           INTENT(in) :: r        ! Molecular positions (3,n)
     REAL, DIMENSION(:,:),           INTENT(in) :: e        ! Orientation vectors (3,n) or quaternions (0:3,n)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: v        ! Molecular velocities (3,n)
     REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) :: w        ! Angular velocities/momenta (3,n)
 
     INTEGER :: cnf_unit, ioerr, i
 
     ! Open given filename, replacing it if it already exists, will terminate on any errors
     OPEN ( newunit=cnf_unit, file=filename, status='replace', iostat=ioerr )
     IF ( ioerr /= 0 ) THEN
        WRITE ( unit=error_unit, fmt='(a,a,i15)') 'Error opening ', filename, ioerr
        STOP 'Error in write_cnf_mols'
     END IF
 
     ! Write number of molecules to first record, box length to second record
     !WRITE(cnf_unit,'(i15)'  ) n
     !WRITE(cnf_unit,'(f15.8)') box
 
     IF ( n > SIZE ( r, dim=2 ) ) THEN
        WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of r ', n, SIZE ( r, dim=2 )
        STOP 'Error in write_cnf_mols'
     END IF
     IF ( n > SIZE ( e, dim=2 ) ) THEN
        WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of e ', n, SIZE ( e, dim=2 )
        STOP 'Error in write_cnf_mols'
     END IF
 
     IF ( PRESENT ( v ) ) THEN
        IF ( .NOT. PRESENT ( w )    ) THEN
           WRITE ( unit=error_unit, fmt='(a,a,i15)') 'v and w arguments must be present together'
           STOP 'Error in write_cnf_mols'
        END IF
        IF ( n > SIZE ( v, dim=2 ) ) THEN
           WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of v ', n, SIZE ( v, dim=2 )
           STOP 'Error in write_cnf_mols'
        END IF
        IF ( n > SIZE ( w, dim=2 ) ) THEN
           WRITE ( unit=error_unit, fmt='(a,2i15)') 'Error in size of w ', n, SIZE ( w, dim=2 )
           STOP 'Error in write_cnf_mols'
        END IF
 
        ! Write positions, orientation vectors or quaternions, velocities, angular velocities/momenta
        DO i = 1, n
           WRITE ( unit=cnf_unit, fmt='(*(f15.10))') r(:,i), e(:,i), v(:,i), w(:,i)
        END DO
 
     ELSE
 
        ! Write positions, orientation vectors or quaternions
        DO i = 1, n
           WRITE ( unit=cnf_unit, fmt='(*(f15.10))') r(:,i), e(:,i)
        END DO
 
     END IF
 
     CLOSE ( unit=cnf_unit )
 
   END SUBROUTINE write_cnf_mols
   
   SUBROUTINE read_randm_particles(N,box,dimnsn,R,V,seed,temp,yes_pos,yes_vel)
     IMPLICIT NONE
     INTEGER,                          INTENT(in)    :: n        ! Number of atoms
     INTEGER,                          INTENT(in)    :: dimnsn   ! Dimensions
     INTEGER,                          INTENT(in)    :: seed     ! Random seed for the assigning positions and velocity
     REAL*8,                           INTENT(out)   :: box      ! Simulation box length (assumed cubic)
     REAL*8, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: r        ! Atomic positions (3,n)
     REAL*8, DIMENSION(:,:), OPTIONAL, INTENT(out)   :: v        ! Atomic velocities (3,n)
     REAL*8,                           INTENT(in)    :: temp     ! Temperature of the system
     LOGICAL,                          INTENT(in)    :: yes_vel  ! If we have random velocities or not
     LOGICAL,                          INTENT(in)    :: yes_pos  ! If we have ranodm positions or not
     
     INTEGER :: i, seed1, seed2
     REAL*8  :: boxh, cc
     !cc = 8.0d0/3.141592653589793d0
     cc = 3.0d0
     boxh = 0.5d0*box
     seed1 = seed+1
     seed2 = seed-1 
     CALL SRAND(seed1)
     ! FORTRAN Random number generator only generates numbers between 0 .le. x .lt. 1
     ! We generate the numbers such that the positions lie between (-boxh,boxh).  
     DO i = 1,N,1
        IF (yes_pos) THEN
          R(1,i) = 2.0d0*box*RAND() - box
          R(2,i) = 2.0d0*box*RAND() - box
          R(3,i) = 2.0d0*box*RAND() - box
        ELSE IF (yes_vel) THEN
          V(1,i) = 2.0d0*DSQRT(cc*temp)*RAND() - DSQRT(cc*temp)
          V(2,i) = 2.0d0*DSQRT(cc*temp)*RAND() - DSQRT(cc*temp)
          V(3,i) = 2.0d0*DSQRT(cc*temp)*RAND() - DSQRT(cc*temp)
        ELSE
          V(:,i) = 0.0d0
          R(:,i) = 0.0d0
        END IF 

     END DO

   END SUBROUTINE read_randm_particles   
   
   SUBROUTINE BCC_LATT(R,N,Dimnsn,B)
     IMPLICIT NONE
     REAL*8,DIMENSION(Dimnsn,N) :: R
     REAL*8 :: A, B       !Lattice Spacing, Box length
     INTEGER :: N      !Number of particles
     INTEGER :: Dimnsn !Dimensions
     INTEGER :: Nc     !Number of unit cell lattices
     INTEGER :: Ncheck, M, i, j, k, ij, kct

     Nc = INT((DBLE(N)/2.0d0)**(1.0d0/3.0d0) + 0.1d0)
     Ncheck = 2*Nc**3
     IF( Ncheck .LT. N) THEN
       PRINT*, "Something is wrong!!!"
       Nc = Nc + 1
     END IF

     A = B/DBLE(Nc)

     R(1,1) = 0.0d0   ; R(2,1) = 0.0d0   ; R(3,1) = 0.0d0
     R(1,2) = A/2.0d0 ; R(2,2) = A/2.0d0 ; R(3,2) = A/2.0d0

     M = 0
     kct = 0
     DO i=1,Nc,1
     DO j=1,Nc,1
     DO k=1,Nc,1
       DO ij=1,2
         IF (kct .LT. N) THEN
           R(1,ij+m) = R(1,ij) + A*DBLE(k-1) ;
           R(2,ij+m) = R(2,ij) + A*DBLE(j-1) ;
           R(3,ij+m) = R(3,ij) + A*DBLE(i-1)

         END IF
         kct = kct + 1
      END DO
      M = M + 2
    END DO
    END DO
    END DO

  END SUBROUTINE BCC_LATT

   SUBROUTINE FCC_Latt(R,N,Dimnsn,B)
     IMPLICIT NONE
     REAL*8,DIMENSION(Dimnsn,N) :: R
     REAL*8 :: A, B
     INTEGER :: Dimnsn, N, Nc
     INTEGER :: Ncheck, i, j, k, ij, kct, m
     
     Nc = INT((DBLE(N)/4.0d0)**(1.0d0/3.0d0) + 0.1d0)
     Ncheck = 4*Nc**3
     IF( Ncheck .LT. N) THEN
       PRINT*, "Something is wrong!!!"
       Nc = Nc + 1
     END IF

     A = 0.5d0*B/DBLE(Nc)
     WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Lattice spacing =', 2.0d0*A

     R(1,1) = 0.0d0   ; R(2,1) = 0.0d0   ; R(3,1) = 0.0d0
     R(1,2) = 0.0d0   ; R(2,2) = A       ; R(3,2) = A
     R(1,3) = A       ; R(2,3) = 0.0d0   ; R(3,3) = A
     R(1,4) = A       ; R(2,4) = A       ; R(3,4) = 0.0d0


     M = 0
     kct = 0
     DO i=1,Nc,1
     DO j=1,Nc,1
     DO k=1,Nc,1
       DO ij=1,4
         IF (kct .LT. N) THEN
           R(1,ij+m) = R(1,ij) + 2.0d0*A*DBLE(k-1)
           R(2,ij+m) = R(2,ij) + 2.0d0*A*DBLE(j-1)
           R(3,ij+m) = R(3,ij) + 2.0d0*A*DBLE(i-1)
         END IF
         kct = kct + 1
       END DO
       M = M + 4
     END DO
     END DO
     END DO

  END SUBROUTINE FCC_Latt

   SUBROUTINE PBC(pos,blow,bhi,edge)
     IMPLICIT NONE
     REAL*8 ::pos,blow,bhi,edge
     IF(pos .LT. blow) THEN
       pos = pos + edge
     ELSE IF(pos .gt. bhi) THEN
       pos = pos - edge
     END IF
   END SUBROUTINE PBC

   SUBROUTINE Check_pbc(R,N,dimnsn,box)
     IMPLICIT NONE
     REAL*8, INTENT(in)  :: Box
     INTEGER, INTENT(in) :: Dimnsn,N
     REAL*8,DIMENSION(dimnsn,N), INTENT(inout) :: R
     INTEGER :: i, j

     DO i = 1,N
       DO j = 1, Dimnsn
         CALL PBC(r(j,i),0.0d0,box,box)
       END DO
     END DO
   END SUBROUTINE Check_pbc

END MODULE config_io_module
