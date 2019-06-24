!! Molecular Dynamics for Soft Spherese !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Mol_Dyn
      IMPLICIT REAL*8 (A-H, O-Z)
      LOGICAL LG(4), LENER, LSCALE
      COMMON /INTGRS/ IN(7)
      COMMON /REALS / RL(26)

      EQUIVALENCE (IN(1), KSAMPL), (IN(4),NATOM), (IN(6),NSTEP)
      EQUIVALENCE (LG(1), LENER), (LG(2), LSCALE), (IN(5),Natom1)
      EQUIVALENCE (RL(16),STPSQH)

      !!!! Initialization !!!!
      CALL SSFLAG(KRDF,KWRITE,MAXSTP,MAXEQB,NDEAD,DENER)
      CALL SSPARM
      CALL SSFCC
      CALL SSPRNT
      CALL SSIVEL
      CALL SSEVAL
      CALL SSX2SC(STPSQH,NATOM)
      !!!! Equilibriation !!!!
      DO 5 NSTEP = 1, MAXEQB
        CALL SSPDCT
        CALL SSEVAL
        CALL SSCORR
        CALL SSPRTY
        IF (MOD(NSTEP,KWRITE) .EQ. 0) CALL SSMNTR (1)
        IF (NSTEP .LT. (MAXEQB-NDEAD)) THEN 
                IF (LSCALE .OR. LENER) CALL SSCALE(DENER)
        END IF
    5 CONTINUE
      

      CALL SSRSET

      DO 7 NSTEP = 1, MAXSTP
        !PRINT*, NSTEP
        CALL SSPDCT
        CALL SSEVAL
        CALL SSCORR
        
        IF (MOD(NSTEP,KSAMPL) .EQ. 0) CALL SSPRTY
        IF (MOD(NSTEP,KWRITE) .EQ. 0) CALL SSMNTR(1)
        IF (MOD(NSTEP,KRDF) .EQ. 0) CALL SSGOFR
    7 CONTINUE
      
      STOP

END PROGRAM

!!!!!! Run parameters and flags typically set by the user !!!!!!!
SUBROUTINE SSFLAG(KRDF, KWRITE, MAXSTP, MAXEQB, NDEAD, DENER)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LG(4),LENER,LSCALE,LSHIFT
      COMMON / INTGRS / IN(7)
      COMMON / LGCLS / LG
      COMMON / REALS / RL(26)

      EQUIVALENCE (IN(1),KSAMPL), (IN(2),KSQRT), (IN(4), NATOM)
      EQUIVALENCE (IN(6), NSTEP)
      EQUIVALENCE (LG(1),LENER), (LG(2),LSCALE), (LG(3), LSHIFT)
      EQUIVALENCE (RL(4),DENSITY), (RL(11),RCUT), (RL(14),STEP)
      EQUIVALENCE (RL(20),TEMP)

      !!!! Run Parameters !!!!
      NATOM = 500
      DENSITY = 0.7d0
      TEMP = 1.2d0
      STEP = 0.0001d0
      RCUT = 2.5d0

      !!!! Run Flags !!!!
      MAXEQB =  100000
      MAXSTP = 500000
      KWRITE = 200
      KRDF = 25000
      KSQRT = 10
      Ksampl = 10
      Nstep = 0
      Ndead = 1000

      !!!! To scale velocities for temperarture during equilibrium, set
      !LSCALE = .TRUE.
      Lscale = .TRUE.

      !!!! To scale velocitites for set point total energy, set LENER =
      !.TRUE.
      IF( .NOT. LScale) THEN
              LENER = .FALSE.
              ! if scaling energy load set point value
              IF (LENER) DENER = -2.0d0
      END IF
      ! if using shifted force potential set lshift = .TRUE.
      ! for potentials that are truncated set lshift = .FALSE. 
      LSHIFT = .FALSE. 

      RETURN
END SUBROUTINE

!!!!!!! Calculate parameters from run parameters given in SSFLAG
SUBROUTINE SSPARM
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON /ALPHA / ALFA0, ALFA1, ALFA3, ALFA4, ALFA5
       COMMON /INTGRS / IN(7)
       COMMON /REALS / RL(26)
       EQUIVALENCE (IN(4),NATOM), (IN(5),NATOM1), (IN(7),NVDELS)
       EQUIVALENCE (RL(1),AHEAT) ,(RL(2),CUBE),(RL(3),CUBEH),(RL(4),DENSITY), (RL(7), Eshft), (RL(8), Etail), (RL(9), Fnatom) 
       EQUIVALENCE (RL(4), Fshft), (Rl(11),Rcut), (RL(12),Rdel), (RL(13),Rlist), (RL(14),Step), (RL(15),Stepsq), (RL(16),Stepsqh)
       EQUIVALENCE (RL(20), TEMP), (RL(22),Vdel), (RL(25),Vmax), (RL(26),Vtail)

    !functions of particle number and time-step
    Fnatom = Natom
    !PRINT*, Natom, Fnatom
    Natom1 = Natom - 1
    !PRINT*, Natom1
    Stepsq = Step*Step
    Stepsqh = 0.5d0*Stepsq
    
    !size of cube, cut-off and list distances, scale factor for velocities, and smapling increments for g(r)
    Vol = Fnatom/Density
    Cube = Vol**(1.0d0/3.0d0)
    CubeH = 0.5d0*Cube
    Rlist = Rcut + 0.3d0
    IF(Rlist .GT. CubeH) Rlist = CubeH
    IF(Rcut .GT. CubeH) Rcut = CubeH - 0.35d0
    Aheat = 3.0d0*Fnatom*Stepsq*temp
    Rdel = 0.025d0

    !parameters in the predictor-corrector method
    Alfa0 = 3.0d0/16.0d0
    Alfa1 = 251.0d0/360.0d0
    Alfa3 = 11.0d0/18.0d0
    Alfa4 = 1.0d0/6.0d0
    Alfa5 = 1.0d0/60.0d0
    
    !Shifted Force Constants 
    RCinv = 1.0d0/Rcut
    RC6inv = RCinv**6
    Eshft = RC6inv*(28.0d0-52.0d0*RC6inv)
    Fshft = 48.0d0*RCinv*RC6inv*(RC6inv-0.5d0)

    !long range corrections to pressure and internal energy and parameters for velocity distribution 
    PI = 3.14159265d0
    RC3 = Rcut**3
    RC9 = RC3**3
    Etail = 8.0d0*PI*Density*(1.0d0/(9.0d0*RC9) - 1.0d0/(3.0d0*RC3))
    Vtail = 96.0d0*PI*Density*(0.5d0/(3.0d0*RC3) - 1.0d0/(9.0d0*RC9))

    Vmax = 5.0
    Vdel = 0.05
    NVdels = INT(2.0*Vmax/Vdel + 1.01)

    RETURN
END SUBROUTINE 


!!!!!! Assign intial positons based on face-centred cubic lattice !!!!!!
SUBROUTINE SSFCC 
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON /INTGRS / IN(7)
       COMMON /REALS / RL(26)
       COMMON /POSIT / X(10000), Y(10000), Z(10000)
       EQUIVALENCE (IN(4),NATOM)
       EQUIVALENCE (RL(2),CUBE), (RL(5),Dist), (RL(9), Fnatom)

       !number of basic lattice units, then check whether user asked for non-fcc number of atoms 
       Nunit = INT((Fnatom/4.0)**(1.0/3.0) + 0.1)
   5   Ncheck = 4*(Nunit**3)
       IF (Ncheck .LT. Natom) THEN
               Nunit = Nunit + 1 
               GOTO 5
       END IF

       !Set lattice distance based on cube of side= cube
       DIST = 0.5d0*Cube/DFLOAT(Nunit)

       !Assign positions of first four atoms
       X(1) = 0.0d0
       Y(1) = 0.0d0
       Z(1) = 0.0d0
       X(2) = 0.0d0
       Y(2) = DIST
       Z(2) = DIST
       X(3) = DIST
       Y(3) = 0.0d0
       Z(3) = DIST
       X(4) = DIST
       Y(4) = DIST
       Z(4) = 0.0d0

       !Replicate first four postions over NUnits
       M=0
       KCT = 0 
       DO 12 i=1,Nunit
       DO 12 j=1,Nunit
       DO 12 k=1,Nunit
         DO 10 ij=1,4
         !If number of positions assingned exceeds the number of atoms, then get out
           IF (KCT .LT. Natom) THEN
                X(ij+M) = X(ij) + 2.0d0*dist*(k-1)
                Y(ij+M) = Y(ij) + 2.0d0*dist*(j-1)
                Z(ij+M) = Z(ij) + 2.0d0*dist*(i-1)
           END IF
           KCT = KCT + 1
   10    CONTINUE
         M = M+4
   12  CONTINUE
       !DO i=1,Natom       
       !  PRINT*, X(i), Y(i), Z(i) 
       !END DO 
       RETURN

END SUBROUTINE

!!!!!! Assign Random Intial Velocities !!!!!!!
SUBROUTINE SSIVEL
        IMPLICIT REAL*8 (A-H,O-Z)
        !REAL*4 ROULET
        COMMON /INTGRS / IN(7)
        COMMON /REALS / RL(26)
        COMMON /VEL / X1(10000), Y1(10000), Z1(10000)

        EQUIVALENCE (IN(4),Natom)
        EQUIVALENCE (RL(9),Fnatom), (RL(23),Velsq)

        DATA Sumx, Sumy, Sumz, Velsq/4*0.0d0/
        MSEED = -30509
        CALL SRAND(Mseed)
        ! Assign Random Velocity Components on (-1,+1) : Need ROULET for this
        DO 100 i=1,Natom
            X1(i) = 2.0d0*RAND() - 1.0d0 !ROULET(MSEED)!DBLE(ROULET(MSEED))!
            Y1(i) = 2.0d0*RAND() - 1.0d0 !ROULET(MSEED)!DBLE(ROULET(MSEED+1))
            Z1(i) = 2.0d0*RAND() - 1.0d0 !ROULET(MSEED)!DBLE(ROULET(MSEED+2))!
            Sumx = Sumx + X1(i)
            Sumy = Sumy + Y1(i)
            Sumz = Sumz + Z1(i)
   100  CONTINUE     
                  
        !Scale velocities so that total linear momentum is zero
        DO 120 i=1,Natom
            X1(i) = X1(i) - Sumx/Fnatom
            Y1(i) = Y1(i) - Sumy/Fnatom
            Z1(i) = Z1(i) - Sumz/Fnatom
            Velsq = Velsq + X1(i)**2 + Y1(i)**2 + Z1(i)**2
   120  CONTINUE

        !Scale the velocities to set point temperature or to set-point total energy
        CALL SSCALE(0.0d0)
        !DO i=1,Natom        
        !  PRINT*, X1(i), Y1(i), Z1(i)
        !END DO
        RETURN
END SUBROUTINE

! Use fifth order Taylor series to predict positions and their derivatives at the next time-step
SUBROUTINE SSPDCT
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON /ACCEL / X2(10000), Y2(10000), Z2(10000)
        COMMON /DERIV3 / X3(10000), Y3(10000), Z3(10000)
        COMMON /DERIV4 / X4(10000), Y4(10000), Z4(10000)
        COMMON /DERIV5 / X5(10000), Y5(10000), Z5(10000)
        COMMON /FORCE / Fx(10000), Fy(10000), Fz(10000)
        COMMON /INTGRS / IN(7)
        COMMON / POSIT /  X(10000), Y(10000), Z(10000)
        COMMON /VEL /  X1(10000), Y1(10000), Z1(10000)
        EQUIVALENCE (IN(4),Natom)

        !velocities and higher derivative are implicitly scaled by appropriate factors of the time-step (Subroutine ODINIT)
        !Note that the numerical coefficient in the Taylor series form a Pascal Triangle
        DO 200 i=1,Natom
            X(i) = X(i) + X1(i)+X2(i)+X3(i)+X4(i)+X5(i)
            Y(i) = Y(i) + Y1(i)+Y2(i)+Y3(i)+Y4(i)+Y5(i)
            Z(i) = Z(i) + Z1(i)+Z2(i)+Z3(i)+Z4(i)+Z5(i)
           
            X1(i) =  X1(i)+2.0d0*X2(i)+3.0d0*X3(i)+4.0d0*X4(i)+5.0d0*X5(i)
            Y1(i) =  Y1(i)+2.0d0*Y2(i)+3.0d0*Y3(i)+4.0d0*Y4(i)+5.0d0*Y5(i)
            Z1(i) =  Z1(i)+2.0d0*Z2(i)+3.0d0*Z3(i)+4.0d0*Z4(i)+5.0d0*Z5(i)

            X2(i) =  X2(i)+3.0d0*X3(i)+6.0d0*X4(i)+10.0d0*X5(i)
            Y2(i) =  Y2(i)+3.0d0*Y3(i)+6.0d0*Y4(i)+10.0d0*Y5(i)
            Z2(i) =  Z2(i)+3.0d0*Z3(i)+6.0d0*Z4(i)+10.0d0*Z5(i)

            X3(i) =  X3(i)+4.0d0*X4(i)+10.0d0*X5(i)
            Y3(i) =  Y3(i)+4.0d0*Y4(i)+10.0d0*Y5(i)
            Z3(i) =  Z3(i)+4.0d0*Z4(i)+10.0d0*Z5(i)

            X4(i) =  X4(i)+5.0d0*X5(i)
            Y4(i) =  Y4(i)+5.0d0*Y5(i)
            Z4(i) =  Z4(i)+5.0d0*Z5(i)
           
           ! Initialize Forces arrays
           Fx(i) = 0.0d0
           Fy(i) = 0.0d0
           Fz(i) = 0.0d0
    200 CONTINUE 
        RETURN

END SUBROUTINE 

!Evaluate Forces on each atom using predicted positions
SUBROUTINE SSEVAL
        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL LG(4), Lupdat
        COMMON / INTGRS/ IN(7)
        COMMON / LGCLS / LG
        COMMON / NABLST / LIST(15000000) , Npoint(10000)
        COMMON / POSIT / X(10000), Y(10000), Z(10000)
        COMMON / REALS / RL(26)
        COMMON / Zloop / XI, YI, ZI, RIJ, I, J, Jbegin, Jend

        EQUIVALENCE (IN(2),KSQRT), (IN(3),Nabors), (IN(4),Natom)
        EQUIVALENCE (IN(5), NATOM1), (IN(6), Nstep)
        EQUIVALENCE (RL(6),Energy), (RL(24),Virial), (LG(4), Lupdat)
       !Initialize instantaneos potential energy and virial
        Energy = 0.0d0
        Virial = 0.0d0
        
       !set counter and flag for neighbour list update at intervals
        Nabors = 0
        Lupdat = .FALSE.
        IF(MOD(Nstep,Ksqrt) .EQ. 0) Lupdat = .TRUE.
        
        !Begin loop over atoms
        DO 300 i=1,Natom1 
       !Neighbor list to find neighbours of atom 1
          Jbegin = Npoint(i)
          Jend = Npoint(i+1) - 1

       !at Ksort intervals, use all (Ntaom - 1) atoms as neighbours
          IF (Lupdat) THEN
               Npoint(i) = NAbors + 1
               Jbegin = i+1
               Jend = Natom
          END IF
       !store position of atom i
          Xi = X(i)
          Yi = Y(i)
          Zi = Z(i)
          CALL SSINLP
    300 CONTINUE
         
        IF (Lupdat) Npoint(Natom) = Nabors + 1
        RETURN
END SUBROUTINE

!!!! Perform inner loop over the atoms to get force on i
SUBROUTINE SSINLP
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER k
        LOGICAL LG(4),Lshift, Lupdat
        COMMON / INTGRS / IN(7)
        COMMON /FORCE / Fx(10000), Fy(10000), Fz(10000)
        COMMON /LGCLS / LG
        COMMON /NABLIST / List(15000000), Npoint(10000)
        COMMON /POSIT / X(10000), Y(10000), Z(10000)
        COMMON /REALS / RL(26)
        COMMON /ZLOOP / Xi, Yi, Zi, Rij, i,j, Jbegin,Jend
        EQUIVALENCE (IN(4), Natom)   
        EQUIVALENCE (RL(2),Cube), (RL(3),CubeH), (RL(6),Energy)                                                              
        EQUIVALENCE (RL(7),Eshft), (RL(10),Fshft), (RL(11),Rcut)                     
        EQUIVALENCE (RL(24),Virial), (LG(3),Lshift), (LG(4),Lupdat)
                
        
        DO 400 jx = Jbegin,JEnd
          j = Jx
          IF (.NOT. Lupdat) J = List(jx) 
          
        !components of vector between atoms i and j
          Xij = Xi-X(j)                              
          Yij = Yi-Y(j)                                     
          Zij = Zi-Z(j)                                                                                                                                         
                                            
        !find image of j closest to i
             IF (Xij .LT. -CubeH) Xij = Xij + Cube                                         
             IF (Xij .GT.  CubeH) Xij = Xij - Cube                                                     
             IF (Yij .LT. -CubeH) Yij = Yij + Cube                                         
             IF (Yij .GT.  CubeH) Yij = Yij - Cube
             IF (Zij .LT. -CubeH) Zij = Zij + Cube                                         
             IF (Zij .GT.  CubeH) Zij = Zij - Cube
          Rsq = Xij*Xij + Yij*Yij + Zij*Zij
          Rij = DSQRT(Rsq)
          IF (Lupdat .AND. RIJ .LE. CubeH) CALL SSUPDT
         
          IF (Rij .LE. Rcut) THEN
             Rsqinv = 1.0d0/Rsq
             R6inv = Rsqinv*Rsqinv*Rsqinv
            
        !Enr is U(r), for is (-1/r)(du/dr)
             Enr = 4.0d0*R6inv*(R6inv - 1.0d0)
             For = Rsqinv*48.0d0*R6inv*(R6inv - 0.5d0)
             IF (Lshift) Enr = Enr + Eshft + Rij*Fshft                                         
             IF (Lshift) For = For - Rsqinv*Rij*Fshft              
               Fx(i) = Fx(i) + For*Xij                                 
               Fx(j) = Fx(j) - For*Xij
               Fy(i) = Fy(i) + For*Yij                                 
               Fy(j) = Fy(j) - For*Yij
               Fz(i) = Fz(i) + For*Zij                                 
               Fz(j) = Fz(j) - For*Zij
             Energy = Energy + Enr
             Virial = Virial - For*Rsq
          END IF
   400  CONTINUE
        !DO k=1,Natom
        !  PRINT*, k, Fx(k), Fy(k), Fz(k)
        !END DO        
        RETURN
END SUBROUTINE


!!!!!! Correct predicted postions and their derivatives !!!!!!!
SUBROUTINE SSCORR
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON /ACCEL / X2(10000), Y2(10000), Z2(10000)
        COMMON /ALPHA / Alfa0, Alfa1, Alfa3, Alfa4, Alfa5
        COMMON /DERIV3 / X3(10000), Y3(10000), Z3(10000)
        COMMON /DERIV4 / X4(10000), Y4(10000), Z4(10000)
        COMMON /DERIV5 / X5(10000), Y5(10000), Z5(10000)
        COMMON /FORCE / Fx(10000), Fy(10000), Fz(10000)
        COMMON /INTGRS / IN(7)
        COMMON / POSIT /  X(10000), Y(10000), Z(10000)
        COMMON / REALS / RL(26)
        COMMON /VEL /  X1(10000), Y1(10000), Z1(10000)
        
        EQUIVALENCE (IN(4), NATOM), (RL(2),CUBE), (RL(16),Stpsqh)
         

        DO 500 i=1,Natom
          Xerror = Stpsqh*Fx(i) - X2(i)
          Yerror = Stpsqh*Fy(i) - Y2(i)
          Zerror = Stpsqh*Fz(i) - Z2(i) 
          
          X(i)  = X(i)  + Xerror*Alfa0
          X1(i) = X1(i) + Xerror*Alfa1
          X2(i) = X2(i) + Xerror
          X3(i) = X3(i) + Xerror*Alfa3
          X4(i) = X4(i) + Xerror*Alfa4
          X5(i) = X5(i) + Xerror*Alfa5

          Y(i)  = Y(i)  + Yerror*Alfa0
          Y1(i) = Y1(i) + Yerror*Alfa1
          Y2(i) = Y2(i) + Yerror
          Y3(i) = Y3(i) + Yerror*Alfa3
          Y4(i) = Y4(i) + Yerror*Alfa4
          Y5(i) = Y5(i) + Yerror*Alfa5

          Z(i)  = Z(i)  + Zerror*Alfa0
          Z1(i) = Z1(i) + Zerror*Alfa1
          Z2(i) = Z2(i) + Zerror
          Z3(i) = Z3(i) + Zerror*Alfa3
          Z4(i) = Z4(i) + Zerror*Alfa4
          Z5(i) = Z5(i) + Zerror*Alfa5

         !Apply periodic boundary conditions
          CALL SSPBC (X(i), Y(i), Z(i), 0.0d0, Cube, Cube) 
  500   CONTINUE            
        RETURN
END SUBROUTINE 


!!!!! Apply periodic boundary conditions to point (x,y,z) from Bhi to Blow !!!!!
SUBROUTINE SSPBC(X, Y, Z, Blow, Bhi, Edge) 
        IMPLICIT REAL*8 (A-H,O-Z)
        
        IF (X .LT. Blow) THEN
          X = X + EDGE
        ELSEIF (X .GT. Bhi) THEN
          X = X - EDGE
        END IF

        IF (Y .LT. Blow) THEN
          Y = Y + EDGE
        ELSEIF (Y .GT. Bhi) THEN
          Y = Y - EDGE
        END IF

        IF (Z .LT. Blow) THEN
          Z = Z + EDGE
        ELSEIF (Z .GT. Bhi) THEN
          Z = Z - EDGE
        END IF

        RETURN
END SUBROUTINE


!!!!!! Accumulate property values every Ksampl time-steps !!!!!!
SUBROUTINE SSPRTY
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON / INTGRS / IN(7)
        COMMON /REALS / RL(26)
        COMMON / Vel / X1(10000), Y1(10000), Z1(10000) 
        EQUIVALENCE (IN(4),Natom)
        EQUIVALENCE (RL(6),Energy), (RL(17),Sumenr), (RL(18),Sumvir)
        EQUIVALENCE (RL(19),Sumvsq), (RL(23),Velsq), (RL(24),Virial) 
        
       !Sum square of velocities
        Velsq = 0.0d0
        DO 700 i=1,Natom
          Velsq = velsq + X1(i)**2 + Y1(i)**2 + Z1(i)**2
    700 CONTINUE
        
       !Accumulate sums for running averages
        Sumenr = sumenr + Energy
        Sumvir = Sumvir + Virial
        Sumvsq = Sumvsq + Velsq
        RETURN 
END SUBROUTINE

!Print Properties at intervals 
SUBROUTINE SSMNTR(Ksampl)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON / INTGRS / IN(7)
        COMMON / REALS / RL(26)
        EQUIVALENCE (IN(3),Nabors) , (IN(6),Nstep)
        EQUIVALENCE (RL(4), Density), (RL(6), Energy), (RL(8),Etail)
        EQUIVALENCE (RL(9),Fnatom), (RL(15),Stepsq), (RL(17),Sumenr)
        EQUIVALENCE (RL(18),Sumvir), (RL(19),Sumvsq), (RL(21),Tmpave)
        EQUIVALENCE (RL(23),Velsq), (RL(24),Virial), (RL(26),Vtail)

        !Instantaneous property values at current time-step
        Ekin = Velsq/(2.0d0*Fnatom*Stepsq)
        Epot = Energy/Fnatom + Etail
        Etotal = Ekin + Epot
        Tmp = 2.0d0*Ekin/3.0d0
        Pres = Density*(Tmp - Virial/(3.0d0*Fnatom) - Vtail/3.0d0)

        !Average property values over duration of run
        Denom = Fnatom*DFLOAT(Nstep/Ksampl)
        Ekave = Sumvsq/(2.0d0*Denom*Stepsq)
        Epave = Sumenr/Denom + Etail
        Tmpave = 2.0d0*Ekave/3.0d0
        Preave = Density*(Tmpave - Sumvir/(3.0d0*Denom) - Vtail/3.0d0)
        !print*, "Nstep : ", Nstep
        !PRINT*, "Sampling time : ", Ksampl
        !PRINT*, "Denominator :", Denom

        !Get the order parameters and H-function
        CALL SSORDR(Alamda,Alam1)
        !CALL SSVELO(Hinst)
        WRITE(6,940) Nstep, Tmp, Tmpave, Epot, Epave, Pres, Preave, Alam1, Hinst, Nabors, Etotal
   940  FORMAT(2X,I6,6F9.4,2(1X,F7.3),1X,I7,1X,F8.5)
        RETURN
END SUBROUTINE

!!!!! At Ksqrt intervals, sample for g(r) and update the neighbour list !!!!!
SUBROUTINE SSUPDT
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON / INTGRS / IN(7)
        COMMON / NABLIST / LIST(15000000), Npoint(10000)
        COMMON / RDF / Ngofr(60000), Gr(60000)
        COMMON / REALS / RL(26)
        COMMON / ZLOOP / Xi, Yi, Zi, Rij, i, j, Jbegin, Jend
        EQUIVALENCE (IN(3), Nabors), (RL(12),Rdel), (RL(13),Rlist)

        !Sampling for g(r)
        Nshell = INT(Rij/Rdel + 0.5d0)
        Ngofr(Nshell) = Ngofr(Nshell) + 1

        IF(Rij .LE. Rlist) THEN
          Nabors = Nabors + 1
          LIST(Nabors) = J
        END IF
        RETURN
END SUBROUTINE

!!!!!!! Compute Verlet's translational order parameter !!!!!!!!
SUBROUTINE SSORDR(Alamda, Alam1)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON / INTGRS / IN(7)
        COMMON / POSIT / X(10000), Y(10000), Z(10000)
        COMMON / REALS / RL(26)
        EQUIVALENCE (IN(4),Natom), (RL(5),Dist)
          Pi4D = 4.0d0*3.14159265d0/Dist
          Alamda = DCOS(Pi4D*X(1)) + DCOS(Pi4D*Y(1)) + DCOS(Pi4D*Z(1))
          Alam1 = 0.0d0
          Fatom = DFLOAT(Natom)

        DO 1000 i=2,Natom
          Alamda = Alamda + DCOS(Pi4D*X(i)) + DCOS(Pi4D*Y(i)) + DCOS(Pi4D*Z(i))
          Alam1 = Alam1 + DCOS(Pi4D*(X(i) - X(1))) + DCOS(Pi4D*(Y(i) - Y(1))) + DCOS(Pi4D*(Z(i) - Z(1)))
 1000   CONTINUE
          Alamda = Alamda/(3.0d0*Fatom)
          Alam1 = Alam1/(3.0d0*(Fatom - 1.0d0))
        RETURN
END SUBROUTINE

!!!!!!! Accumalate H-function !!!!!!!!
!SUBROUTINE SSVELO(Hinst)
!        IMPLICIT REAL*8 (A-H,O-Z)
!        COMMON / INTGRS / IN(7)
!        COMMON / REALS / RL(26)
!        COMMON / Vel / X1(10000), Y1(10000), Z1(10000)
!        DIMENSION NVinx(60000), NViny(60000), NVinz(60000)
!        EQUIVALENCE (IN(4),Natom), (IN(7),Nvdels), (RL(9),Fnatom)
!        EQUIVALENCE (RL(14),STEP), (RL(22),Vdel), (RL(25),Vmax)
!
!        DO 870 i=1,Nvdels
!          Nvinx(i) = 0
!          Nviny(i) = 0
!          Nvinz(i) = 0
!  870   CONTINUE
!
!        !Sample for instantaneous velocity distribution
!        Vfact = Vdel*STEP
!        Offset = 1.5 + Vmax/Vdel
!        DO 880 i=1,Natom
!          IF (DABS(X1(i)) .LT. Vmax) THEN
!            Nshell = INT(X1(i)/Vfact + Offset)
!            NVinx(Nshell) = NVinx(Nshell) + 1
!              Nshell = INT(Y1(i)/Vfact + Offset)
!              NViny(Nshell) = NViny(Nshell) + 1
!            Nshell = INT(Z1(i)/Vfact + Offset)
!            NVinz(Nshell) = NVinz(Nshell) + 1
!          END IF
!   880  CONTINUE
!
!        !compute instantaneous H-function
!        HH = 0.0d0
!        DO 890 k=1,Nvdels
!          IF (NVinx(k) .GT. 0) THEN
!            FOFV = DFLOAT(NVinx(k))/Fnatom
!            HH = HH + FOFV*DLOG(FOFV)
!          END IF
!          IF (NViny(k) .GT. 0) THEN
!            FOFV = DFLOAT(NViny(k))/Fnatom
!            HH = HH + FOFV*DLOG(FOFV)
!          END IF
!          IF (NVinz(k) .GT. 0) THEN
!            FOFV = DFLOAT(NVinz(k))/Fnatom
!            HH = HH + FOFV*DLOG(FOFV)
!          END IF
!   890  CONTINUE
!        Hinst = HH*Vdel/3.0
!        RETURN
!END SUBROUTINE

!!!!!! Normalize counters for radial distribution function !!!!!!
SUBROUTINE SSGOFR
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON / INTGRS / IN(7)
        COMMON / REALS / RL(26)
        COMMON / RDF / Ngofr(60000), Gr(60000)
        EQUIVALENCE (IN(7),Ksqrt), (IN(4),Natom), (IN(6),Nstep)
        EQUIVALENCE (RL(3),CubeH), (RL(4),Density), (RL(12),Rdel)
        EQUIVALENCE (RL(21),Tmpave)

        PI = 3.14159265d0

       !print heading
        WRITE(6,963) Nstep
        WRITE(6,965) Natom,Density,Tmpave
        WRITE(6,967)

        Nrdels = INT(CubeH/Rdel) - 1
        Origns = (Nstep/Ksqrt)*(Natom/2)

        !Loop over radial shells
        DO 800 j=1,Nrdels
          Gr(j) = 0.0d0
          IF(Ngofr(j) .GT. 0) THEN
            Radius = Rdel*DFLOAT(j)
            Volshl = 4.0d0*PI*Rdel*Radius*Radius + (PI*Rdel**3)/3.0d0
            Gr(j) = DFLOAT(Ngofr(j))/(Density*Origns*Volshl)
            WRITE(6,969) j,Radius,Ngofr(j),Gr(j)
          END IF
   800  CONTINUE

        WRITE(6,971)
   963    FORMAT('1'////T20,'Radial Distribution Function at time-','Step',I6/)
   965    FORMAT(T18,'Natom = ',I4,3X,'Density = ',F6.3,3X,'Ave Temp = ',F6.3//)
   967    FORMAT(29X,'I'6X,'R',6X,'Numeratr',3X,'G(r)'/)
   969    FORMAT(28X,I3,3X,F6.3,3X,I6,F8.3)
   971    FORMAT('1'///)
        RETURN
END SUBROUTINE

!!!!!! Scale Velocities during equilibriation for set-point !!!!!!!
SUBROUTINE SSCALE(DENER)
        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL LG(4), LENER,Lscale
        COMMON / INTGRS / IN(7)
        COMMON / LGCLS / LG
        COMMON / REALS / RL(26)
        COMMON / VEL / X1(10000), Y1(10000), Z1(10000)

        EQUIVALENCE (IN(4),Natom), (LG(1),LENER), (LG(2),Lscale)
        EQUIVALENCE (RL(1),Aheat), (RL(6),Energy), (RL(8),Etail)
        EQUIVALENCE (RL(15),Stepsq), (RL(23),Velsq)
        !factor = 1.0d0
        !scale factor for temperature
        IF (Lscale) factor = DSQRT(Aheat/Velsq)

        !scale factor for total energy
        IF (LENER) THEN
          Epot = Energy/Natom + Etail
          Dek = DENER - Epot
          factor = DSQRT(2.0*Natom*Dek*Stepsq/Velsq)
        END IF
        !Print*, factor
        !apply the scaling
        DO 600 i=1,Natom
          X1(i) = X1(i)*factor
          Y1(i) = Y1(i)*factor
          Z1(i) = Z1(i)*factor
          !PRINT*, i, X1(i), Y1(i), Z1(i)
  600   CONTINUE

        RETURN
END SUBROUTINE

!!!!!!!!! At time zero, scale intial forces to units of acceleration !!!!!!!!
SUBROUTINE SSX2SC(Stpsqh,Natom)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON /ACCEL / X2(10000), Y2(10000), Z2(10000)
        COMMON /FORCE / Fx(10000), Fy(10000), Fz(10000)
          DO 350 i=1,Natom
            X2(i) = Fx(i)*Stpsqh
            Y2(i) = Fy(i)*Stpsqh
            Z2(i) = Fz(i)*Stpsqh
  350     CONTINUE
        RETURN
END SUBROUTINE

!!!!!! Print initial values of run parameters
SUBROUTINE SSPRNT
        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL LG(4), Lshift
        COMMON / INTGRS / IN(7)
        COMMON / LGCLS / LG
        COMMON / REALS / RL(26)
        EQUIVALENCE (IN(4),Natom), (RL(2),Cube), (RL(3),CubeH)
        EQUIVALENCE (RL(4),Density), (RL(11),Rcut), (RL(13),RList)
        EQUIVALENCE (RL(14),Step), (RL(20),TEMP), (LG(3),Lshift)

        WRITE(6,900)
        WRITE(6,902)
        WRITE(6,904)
        WRITE(6,906) Natom
        WRITE(6,904)
          IF (Lshift)   WRITE(6,907)
          IF (.NOT. Lshift) WRITE(6,909)
        WRITE(6,904)
        WRITE(6,902)
        WRITE(6,904)
        WRITE(6,908) Density, Temp
        WRITE(6,910) Step
        WRITE(6,912) Cube,CubeH
        WRITE(6,914) Rcut,RList
        WRITE(6,904)
        WRITE(6,902)
        WRITE(6,916)
        WRITE(6,918)

        !check that order parameter is unity
        CALL SSORDR(Alamda,Alam1)
        WRITE(6,921) Alamda

  900   FORMAT(1H1///)
  902   FORMAT(T20,49('*'))
  904   FORMAT(T20,'*',T68,'*')
  906   FORMAT(T20,'*',6X,'Molecular Dynamics for',I4,'LJ Atoms',T68,'*')
  907   FORMAT(T20,'*',3X,'Fixed Density W/shift-Force Potential',T68,'*')
  908   FORMAT(T20,'*',2X,'Density = ',F7.3,T45,'Temp =',F7.3,T68,'*')
  909   FORMAT(T20,'*',5X,'Fixed Density W/Truncated Potential',T68,'*')
  910   FORMAT(T20,'*',2X,'Time Step = ',F6.3,T68,'*')
  912   FORMAT(T20,'*',2X,'Cube = ',F7.3,T45,'Cube/2 = ',F6.3,T68,'*')
  914   FORMAT(T20,'*',2X,'Rcut = ',F7.3,T45,'Rlist = ',F6.3,T68,'*')
  916   FORMAT(////)
  918   FORMAT(3X,'Tstep',4X,'Temp',4X,'<Temp>',2X,'P Enrgy',2X,'<P enrgy>',4X,'Pres', &
        &      4X,'<Pres>',3X,'Order',2X,'H-func',3X,'Nabors',1X,'TTL Enr')
  921   FORMAT(7X,'0',T65,F6.3)

        RETURN
END SUBROUTINE

!!!!!!! At end of equilibriation, Reinitialize property sums  !!!!!!!
SUBROUTINE SSRSET
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON / RDF / Ngofr(60000), Gr(60000)
        COMMON / REALS / RL(26)
        EQUIVALENCE (RL(3),CubeH), (RL(12),Rdel)
        EQUIVALENCE (RL(17),Sumenr), (RL(18),Sumvir), (RL(19),Sumvsq)

        Sumenr = 0.0d0
        Sumvir = 0.0d0
        Sumvsq = 0.0d0

        Nrdels = INT(CubeH/Rdel)
        DO 650 i=1,Nrdels
          Ngofr(i) = 0
  650   CONTINUE

        WRITE(6,999)
  999   FORMAT(/10X,'**** Equilibriation Complete ***'//)

        RETURN
END SUBROUTINE

BLOCK DATA
IMPLICIT REAL*8 (A-H,O-Z)
COMMON /DERIV3 / X3(10000), Y3(10000), Z3(10000)
COMMON /DERIV4 / X4(10000), Y4(10000), Z4(10000)
COMMON /DERIV5 / X5(10000), Y5(10000), Z5(10000)
COMMON /FORCE / Fx(10000), Fy(10000), Fz(10000)
COMMON /Nablist / LIST(15000000), Npoint(10000)
COMMON /RDF / Ngofr(60000), Gr(60000)
COMMON /REALS / RL(26)

EQUIVALENCE (RL(17),Sumenr), (RL(18),Sumvir), (RL(19),Sumvsq)

DATA X3/ 10000*0.d0/,Y3/10000*0.0d0/,Z3/10000*0.0d0/
DATA X4/ 10000*0.d0/,Y4/10000*0.0d0/,Z4/10000*0.0d0/
DATA X5/ 10000*0.d0/,Y5/10000*0.0d0/,Z5/10000*0.0d0/
DATA FX/ 10000*0.d0/,FY/10000*0.0d0/,FZ/10000*0.0d0/
DATA Ngofr/60000*0/
DATA List/15000000*0/ , Npoint/10000*0/
DATA Sumenr,Sumvir,Sumvsq/3*0.0d0/

END


!!!!!!! Psuedo-Random Number generator !!!!!!!
FUNCTION Roulet(Iseed)
        REAL RR(97)
        COMMON / RAND / RR
        INTEGER Ix1, Ix2, Ix3
        !Parameter value overflows at 2**24
        DATA M1, IA1, IC1/ 31104, 625, 6571/        
        DATA M2, IA2, IC2/ 12960, 1741, 2731/         
        DATA M3, IA3, IC3/ 14000, 1541, 2957/           
        Rm1 = 1.0/M1
        Rm2 = 1.0/M2
        Ix1 = 0 
        Ix2 = 0
        Ix3 = 0
        !Initialize the shuffling vector RR to hold the random numbers
        IF (Iseed .LT. 0) THEN 
          Ix1 = MOD(IC1 - Iseed, M1)
          Ix1 = MOD(IA1*Ix1+ IC1, M1)
          Ix2 = MOD(Ix1, M2)
          Ix1 = MOD(IA1*Ix1+IC1,M1)
          Ix3 = MOD(Ix1,M3) 

        !Load vector RR 
          DO 11 j=1,97        
            Ix1 = MOD(IA1*Ix1+IC1,M1)
            Ix2 = MOD(IA2*Ix2+IC2,M2)
            RR(j) = (FLOAT(Ix1) + FLOAT(Ix2)*RM2)*RM1
    11    CONTINUE
          Iseed = 1
        END IF

        !randomly sample vector RR
        Ix3 = MOD(IA3*Ix3+Ic3,M3)
        j = 1 + INT((97*Ix3)/M3)

        IF(j .GT. 97 .OR. j .LT. 1) WRITE(6,99)
   99     FORMAT(//5X,'Array size for RR violated in Roulet'/)

        !change interval from [0,1] to [-1,1]
          ROULET = 2.0*RR(j) - 1.0

        !replace RR(j) with the next value in the sequence
        Ix1 = MOD(IA1*Ix1 + IC1,M1)
        Ix2 = MOD(IA2*Ix2 + IC2,M2)
        RR(j) = (FLOAT(Ix1) + FLOAT(Ix2)*RM2)*RM1

        RETURN
END 

