      PROGRAM ISING
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=60, N=10000, MED = 10000)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            DIMENSION EXP_RESULTS(-4:4), IPOS(0:L+1)
            CHARACTER*1 SYMBOLS(-1:1)
            COMMON /POSITION/ IPOS
            COMMON /LATTICE/ LATTICE
            COMMON /EXP_RESULTS/ EXP_RESULTS
            COMMON /PARAM/ ENERGY, RMAG

            SYMBOLS(-1) = '-'
            SYMBOLS(1) = '+'

            ITOTAL = L*L
            STEP = 0.001

            OPEN(1,FILE='./B/b1mag.DAT',STATUS='UNKNOWN')
            OPEN(2,FILE='./B/b1ene.DAT',STATUS='UNKNOWN')
            OPEN(3,FILE='./B/b1lattice.DAT',STATUS='UNKNOWN')

            CALL POS()
            CALL SETUPRANDOM()

            rmag = RFMAG(ITOTAL)
            energy = REnergy()

            WRITE(1,*) 0,rmag/Itotal
            WRITE(2,*) 0,energy/Itotal

            DO 10 ICICLE = 1, 3000
                  WRITE(*,*) ICICLE
                  BETA = STEP*ICICLE
                  CALL EXPONENTIAL(BETA)

                  DO 30 I=1,N
                        CALL FLIP(ITOTAL)
 30               END DO

                  WRITE(1,*) ICICLE, RMAG/ITOTAL
                  WRITE(2,*) ICICLE, ENERGY/ITOTAL

 10         END DO 

            DO 50 K=1,L
                  WRITE(3,'(60A2)') (SYMBOLS(LATTICE(I,K)),I=1,L)
50          END DO
            
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
                  DO 20 K=1,L
                        RFMAG = RFMAG + LATTICE(I,K)
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
            COMMON /LATTICE/ LATTICE
            COMMON /POSITION/ IPOS

            REnergy = 0
            DO 10 I=1,L
                  DO 20 K=1,L
                        REnergy = REnergy
     &                        + LATTICE(I,K)*(LATTICE(IPOS(I-1),K) +
     &                        LATTICE(IPOS(I+1),K) +
     &                        LATTICE(I,IPOS(K-1)) +
     &                        LATTICE(I,IPOS(K+1)))
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
     &                        LATTICE(IPOS(i+1),k) +
     &                        LATTICE(i,IPOS(k-1)) +
     &                        LATTICE(i,IPOS(k+1)))

            ENERGY = ENERGY - Delta

      END SUBROUTINE DeltaEnergy
