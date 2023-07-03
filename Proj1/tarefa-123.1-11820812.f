      PROGRAM GENERATE
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            pi = acos(-1.d0)

            N = 200
            dt = 0.04d0
            a1 = 2.d0
            a2 = 4.d0
            w1 = 4*pi
            w2 = 2.5*pi

            CALL ESCREVER(N, dt, a1, a2, w1, w2)

      END PROGRAM GENERATE

      SUBROUTINE ESCREVER (N, dt, a1, a2, w1, w2)
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16) sig, SINAL

            open(newunit=io, file="entradas/entrada-1.1-11820812.in")

            DO 10 i=0,N-1
                  sig = SINAL(a1, a2, w1, w2, dt*i)
                  time = dt*i
                  WRITE (io, *) time, sig
  10        END DO

            close(io)
      END SUBROUTINE ESCREVER

      COMPLEX (KIND=16) FUNCTION SINAL (a1, a2, w1, w2, t)
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            sinal = CMPLX(a1*cos(w1*t) + a2*sin(w2*t), 0.d0, 16)
            RETURN
      END FUNCTION SINAL