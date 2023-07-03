      PROGRAM ISING
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=100, N=10000, MED = 10000)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            DIMENSION EXP_RESULTS(-4:4) ,IPOS(0:L+1)
            CHARACTER*1 SYMBOLS(-1:1)
            COMMON /POSITION/ IPOS
            COMMON /LATTICE/ LATTICE
            COMMON /EXP_RESULTS/ EXP_RESULTS
            COMMON /PARAM/ RMAG

            SYMBOLS(-1) = '-'
            SYMBOLS(1) = '+'

            ITOTAL = L*L

            OPEN(1,FILE='./A/ISING.dat',STATUS='UNKNOWN')
            OPEN(2,FILE='./A/ISINGD.dat',STATUS='UNKNOWN')

            CALL EXPONENTIAL()
            CALL SETUPORDERED()
            CALL POS()

            rmag = RFMAG(ITOTAL)

            ite = 0
            WRITE(2,*) ite, rmag

            DO 20 WHILE (flag < 10) 
      
                  rmag_old = rmag

                  DO 10 I=1,N
                        CALL FLIP(ITOTAL)
 10               END DO

                  WRITE(2,*) ite, rmag

                  IF(abs(rmag-rmag_old)<0.01)then
                        flag = flag + 1
                  ELSE
                        flag = 0
                  END IF

                  ite = ite + 1
                        
 20         END DO

            Acum_mag = 0
            DO 40 I=1,MED
                  CALL FLIP(ITOTAL)
                  Acum_mag = Acum_mag + rmag
 40         END DO
            Acum_mag = Acum_mag/MED
            WRITE(*,*) Acum_mag

            DO 50 K=1,L
                  WRITE(1,'(100A2)') (SYMBOLS(LATTICE(I,K)),I=1,L)
 50         END DO
            
            CLOSE(1) 

      END PROGRAM ISING

      SUBROUTINE EXPONENTIAL()
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (BETA=0.1)
            REAL*8 EXP_RESULTS(-4:4)
            COMMON  /EXP_RESULTS/ EXP_RESULTS

            DO I = -4, 4
            EXP_RESULTS(I)=DEXP(-BETA*I)/(DEXP(-BETA*I)+DEXP(BETA*I))
            END DO

      END SUBROUTINE EXPONENTIAL

      SUBROUTINE POS()
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=100)
            INTEGER IPOS(0:L+1)
            COMMON /POSITION/ IPOS

            DO 10 I=1,L
                  IPOS(I) = I
 10         END DO
                  
            IPOS(0) = L
            IPOS(L+1) = 1

      END SUBROUTINE POS

      SUBROUTINE SETUPORDERED()
            PARAMETER (L=100)
            BYTE, DIMENSION (1:L,1:L) :: LATTICE
            COMMON /LATTICE/ LATTICE

            DO 10 I=1,L
                  DO 20 K=1,L
                        LATTICE(I,K) = 1
 20               END DO 
 10         END DO
      END SUBROUTINE SETUPORDERED

      SUBROUTINE FLIP(ITOTAL)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=100, J=1.0)
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
            END IF

      END SUBROUTINE FLIP

      FUNCTION RFMAG(ITOTAL)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=100)
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

      SUBROUTINE DeltaMag(PM, ITOTAL)
            IMPLICIT REAL*8 (A-H,O-Z)
            PARAMETER (L=100)
            BYTE PM
            COMMON /PARAM/ RMAG

            RMAG = RMAG + (PM/ITOTAL)*2
      
      END SUBROUTINE DeltaMag