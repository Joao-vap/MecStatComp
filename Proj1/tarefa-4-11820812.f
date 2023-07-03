      PROGRAM TRANSINVERSE
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16) , DIMENSION (1000) :: data
            COMMON /DAT/ data

c     Tomar dados do arquivo do espectro
            OPEN(newunit=io, file="saidas/saida-1.1-11820812.out")

            DO 10 i=1,1000
                  READ(io,*,end=1) freq, data(i)
  10        END DO
  1         CONTINUE

            CLOSE(io)

            M = i - 1
            dt  = 0.04d0

c     Chamar a rotina de fourrier que escrevera em data.out
            CALL INVERSE(M, dt)

      END PROGRAM TRANSINVERSE
      
      SUBROUTINE INVERSE (M, dt)
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16) :: yj, yjf

            pi = acos(-1.d0)

            OPEN(newunit=io, file="saidas/saida-4-11820812.out")
            DO 10 j=0, (2*M)-1
                  yj = yjf(M, j, dt)
                  time = j * dt
                  WRITE(io, *) time, yj
  10        END DO
            CLOSE(io)

      END SUBROUTINE INVERSE

      COMPLEX (KIND=16) FUNCTION yjf (M, j, dt)
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16) , DIMENSION (1000) :: Y
            COMMON /DAT/ Y
            COMPLEX (KIND=16) :: zi

            pi = ACOS(-1.d0)
            zi = (0.d0, 1.d0)

            yjf = 0.d0
            DO 10 k=0,M-1
                  arg = k * j / (2.d0 * M)
                  yjf = yjf + Y(k+1) * EXP(-zi * 2.d0 * pi * arg)
  10        END DO
      
            RETURN
      END FUNCTION yjf