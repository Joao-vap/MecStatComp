      PROGRAM ISING
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (Lm = 10, N=10000, MED = 10000)
            BYTE, DIMENSION (1:Lm,1:Lm) :: LATTICE
            DIMENSION EXP_RESULTS(-4:4), IPOS(0:Lm+1)
            CHARACTER*1 SYMBOLS(-1:1)
            COMMON /POSITION/ IPOS
            COMMON /LATTICE/ LATTICE
            COMMON /EXP_RESULTS/ EXP_RESULTS
            COMMON /PARAM/ ENERGY, RMAG

            SYMBOLS(-1) = '-'
            SYMBOLS(1) = '+'

            STEP = 0.01

            OPEN(1,FILE='dout.DAT',STATUS='UNKNOWN')

            DO 100 L = 4, 10
                  CALL POS(L)

                  BETA = 0.5
                  ITOTAL = L*L
                  CALL EXPONENTIAL(BETA)
                  CALL SETUPRANDOM(L)

                  rmag = RFMAG(ITOTAL, L)
                  energy = REnergy(L)

                  flag = 0
                  icontador = 0
                  iquantas = 0
                  Medtime = 0

                  DO 10 ICICLE = 0, 3000

                        rmag_old = rmag

                        DO 30 I=1,N
                              CALL FLIP(ITOTAL, L)
 30                     END DO

                        IF(rmag*rmag_old < 0) THEN
                              WRITE(*,*) Medtime
                              Medtime = Medtime + icontador
                              flag = 0
                              icontador = 0
                              iquantas = iquantas + 1
                        END IF

                        icontador = icontador + 1

 10               END DO 

                  rMedTime = Medtime/iquantas
                  WRITE(1,*) L, rMedTime
      
 100        END DO
            
            CLOSE(1)

      END PROGRAM ISING

      SUBROUTINE EXPONENTIAL(BETA)
            IMPLICIT REAL*8 (A-H,O-Z)
            REAL*8 EXP_RESULTS(-4:4)
            COMMON  /EXP_RESULTS/ EXP_RESULTS

            DO I = -4, 4
            EXP_RESULTS(I)=DEXP(-BETA*I)/(DEXP(-BETA*I)+DEXP(BETA*I))
            END DO

      END SUBROUTINE EXPONENTIAL

      SUBROUTINE POS(L)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (Lm=10)
            INTEGER IPOS(0:Lm+1)
            COMMON /POSITION/ IPOS

            DO 10 I=1,L
                  IPOS(I) = I
 10         END DO
                  
            IPOS(0) = L
            IPOS(L+1) = 1

      END SUBROUTINE POS


      SUBROUTINE SETUPRANDOM(L)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (Lm=10)
            BYTE, DIMENSION (1:Lm,1:Lm) :: LATTICE
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

      SUBROUTINE FLIP(ITOTAL, L)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (Lm=10, J=1.0)
            BYTE, DIMENSION (1:Lm,1:Lm) :: LATTICE
            DIMENSION EXP_RESULTS(-4:4), IPOS(0:Lm+1)
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
                  CALL DeltaEnergy(I, K)
            END IF

      END SUBROUTINE FLIP

      FUNCTION RFMAG(ITOTAL, L)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (Lm=10)
            BYTE, DIMENSION (1:Lm,1:Lm) :: LATTICE
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
            PARAMETER (Lm=10)
            BYTE PM
            COMMON /PARAM/ ENERGY, RMAG

            RMAG = RMAG + PM*2.d0/ITOTAL
      
      END SUBROUTINE DeltaMag

      FUNCTION REnergy(L)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (Lm=10, J=1.0)
            BYTE, DIMENSION (1:Lm,1:Lm) :: LATTICE
            DIMENSION IPOS(0:Lm+1)
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
            PARAMETER (Lm=10, J=1.0)
            BYTE, DIMENSION (1:Lm,1:Lm) :: LATTICE
            DIMENSION IPOS(0:Lm+1)
            COMMON /LATTICE/ LATTICE
            COMMON /PARAM/ ENERGY, RMAG
            COMMON /POSITION/ IPOS

            Delta = 2*J*LATTICE(i,k)*(LATTICE(IPOS(i-1),k) +
     &                                LATTICE(IPOS(i+1),k) +
     &                                LATTICE(i,IPOS(k-1)) +
     &                                LATTICE(i,IPOS(k+1)))

            ENERGY = ENERGY + Delta

      END SUBROUTINE DeltaEnergy
