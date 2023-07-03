      PROGRAM ISING
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60, N=3600, MED=3600)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            DIMENSION EXP_RESULTS(-4:4), IPOS(0:L+1)
            CHARACTER*1 SYMBOLS(-1:1)
            COMMON /LATTICE/ LATTICE
            COMMON /EXP_RESULTS/ EXP_RESULTS
            COMMON /PARAM/ ENERGY, RMAG
            COMMON /IPOS/ IPOS

            SYMBOLS(-1) = '-'
            SYMBOLS(1) = '+'

            ITOTAL = L*L

            OPEN(1,FILE='./C/c1mag2602.dat',STATUS='UNKNOWN')
            OPEN(2,FILE='./C/c1ene2602.dat',STATUS='UNKNOWN')

            CALL SETUPRANDOM()
            CALL POS()
            CALL EXPONENTIAL(0.d0)

            rmag = RFMAG(ITOTAL)
            energy = REnergy()

            STEP = 0.001
            maxCicle = 1750

            WRITE(1,*) 0, rmag/ITOTAL
            WRITE(2,*) 0, ',', energy/ITOTAL

            DO 10 ICICLE = 1, maxCicle
                  WRITE(*,*) ICICLE
                  BETA = STEP*ICICLE
                  CALL EXPONENTIAL(BETA)

                  AcumEnergy = 0
                  DO 20 I=1,N
                        CALL FLIP(ITOTAL)
                        AcumEnergy = AcumEnergy + ENERGY/N
 20               END DO

                  WRITE(1,*) ICICLE*10, rmag
                  WRITE(2,*) ICICLE*10, ',', AcumEnergy/N

 10         END DO 

            DO 30 ICICLE = (maxCicle-1), 0, -1
                  WRITE(*,*) ICICLE
                  BETA = STEP*ICICLE
                  CALL EXPONENTIAL(BETA)

                  AcumEnergy = 0
                  DO 40 I=1,N
                        CALL FLIP(ITOTAL)
                        AcumEnergy = AcumEnergy + ENERGY/N
 40               END DO

                  WRITE(1,*) ICICLE*10, rmag
                  WRITE(2,*) ICICLE*10,',', AcumEnergy/N
 30         END DO 
            
            CLOSE(1) 
            CLOSE(2)

      END PROGRAM ISING

      SUBROUTINE EXPONENTIAL(BETA)
            IMPLICIT REAL*8 (A-H,O-Z)
            REAL*8 EXP_RESULTS(-4:4)
            COMMON  /EXP_RESULTS/ EXP_RESULTS

            DO I = -4, 4
            EXP_RESULTS(I)=DEXP(-BETA*I)/(DEXP(-BETA*I)+DEXP(BETA*I))
            END DO

      END SUBROUTINE EXPONENTIAL

      SUBROUTINE POS()
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60)
            INTEGER IPOS(0:L+1)
            COMMON /POSITION/ IPOS

            DO 10 I=1,L
                  IPOS(I) = I
 10         END DO
                  
            IPOS(0) = L
            IPOS(L+1) = 1

      END SUBROUTINE POS

      SUBROUTINE SETUPRANDOM()
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            COMMON /LATTICE/ LATTICE

            DO 10 I=1,L
                  DO 20 K=1,L
                        if (rand() < 0.5) then
                              LATTICE(I,K) = 1
                        else
                              LATTICE(I,K) = -1
                        end if
 20               END DO 
 10         END DO
      END SUBROUTINE SETUPRANDOM

      SUBROUTINE FLIP(ITOTAL)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60, J=1.0)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            DIMENSION EXP_RESULTS(-4:4), IPOS(0:L+1)
            COMMON /LATTICE/ LATTICE
            COMMON /EXP_RESULTS/ EXP_RESULTS
            COMMON /POSITION/ IPOS
            
            I = 1 + FLOOR(RAND()*L)
            K = 1 + FLOOR(RAND()*L)

            I_S_DELTA_M =  J*(LATTICE(IPOS(I-1),K) +
     &                        LATTICE(IPOS(I+1),K) +
     &                        LATTICE(I,IPOS(K-1)) +
     &                        LATTICE(I,IPOS(K+1)))
     &                        *LATTICE(I,K)

            IF (RAND() < EXP_RESULTS(I_S_DELTA_M)) THEN
                  LATTICE(I,K) = -LATTICE(I,K)
                  CALL DeltaMag(LATTICE(I,K), ITOTAL)
                  CALL DeltaEnergy(I,K)
            END IF

      END SUBROUTINE FLIP

      FUNCTION RFMAG(ITOTAL)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            COMMON /LATTICE/ LATTICE

            RFMAG = 0
            DO 10 I=1,L
                  DO 20 J=1,L
                        RFMAG = RFMAG + LATTICE(I,J)
 20               END DO
 10         END DO
            RFMAG = RFMAG/ITOTAL
      END FUNCTION RFMAG

      SUBROUTINE DeltaMag(pm, ITOTAL)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60)
            BYTE PM
            COMMON /PARAM/ ENERGY, RMAG

            RMAG = RMAG + (PM/ITOTAL)*2
      
      END SUBROUTINE DeltaMag

      FUNCTION REnergy()
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60, J=1.0)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            DIMENSION IPOS(0:L+1)
            COMMON /POSITION/ IPOS
            COMMON /LATTICE/ LATTICE

            REnergy = 0
            DO 10 I=1,L
                  DO 20 K=1,L
                        REnergy = REnergy
     &                 + LATTICE(I,K)*
     &                  (LATTICE(IPOS(I-1),K) +
     &                   LATTICE(IPOS(I+1),K) +
     &                   LATTICE(I,IPOS(K-1)) +
     &                   LATTICE(I,IPOS(K+1)))
 20               END DO
 10         END DO
            REnergy = -J*REnergy/2
      END FUNCTION REnergy

      SUBROUTINE DeltaEnergy(i, k)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60, J=1.0)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            DIMENSION IPOS(0:L+1)
            COMMON /LATTICE/ LATTICE
            COMMON /PARAM/ ENERGY, RMAG
            COMMON /POSITION/ IPOS

            Delta = 2*J*LATTICE(i,k)*(LATTICE(IPOS(i-1),k) +
     &                                LATTICE(IPOS(i+1),k) +
     &                                LATTICE(i,IPOS(k-1)) +
     &                                LATTICE(i,IPOS(k+1)))

            ENERGY = ENERGY - Delta

      END SUBROUTINE DeltaEnergy
