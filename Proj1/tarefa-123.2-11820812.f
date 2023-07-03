      PROGRAM TRANSFOURRIER
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16), DIMENSION (1000) :: data
            COMMON /DAT/ data

c     Tomar dados do arquivo do sinal
            OPEN(newunit=io, file="entradas/entrada-3.2-11820812.in")

            DO 10 i=1,1000
                  READ(io,*,end=1) time, data(i)
  10        END DO
  1         CONTINUE

            CLOSE(io)

            N = i - 1
            dt = time/(N-1)

c     Chamar a rotina de fourrier que escrevera em data.out
            CALL FOURRIER(N, dt)

      END PROGRAM TRANSFOURRIER

      SUBROUTINE FOURRIER (N, dt)
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16) :: Yk, Ykf
            
            pi = acos(-1.d0)

            OPEN(newunit=io, file="saidas/saida-3.2-11820812.out")
            DO 10 k=0,N/2-1
                  Yk = Ykf(N, k) / (N/2)
                  freq = k / (N*dt)
                  WRITE(io, *) freq, Yk
  10        END DO
            CLOSE(io)

      END SUBROUTINE FOURRIER

      COMPLEX (KIND=16) FUNCTION Ykf (N, k)
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16), DIMENSION (1000) :: y
            COMMON /DAT/ y
            COMPLEX (KIND=16) :: zi

            pi = ACOS(-1.d0)
            zi = (0.d0, 1.d0)

            Ykf = 0.d0
            DO 10 j=0,N-1
                  Ykf  = Ykf +  y(j+1) * EXP(2.d0*pi*zi * k * j/N)
  10        END DO

            RETURN
      END FUNCTION Ykf