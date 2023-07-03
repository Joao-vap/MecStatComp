      module Frequencias

      contains

            SUBROUTINE FOURRIER (y, N, dt)
                  IMPLICIT REAL (KIND=8) (a-h,o-z)
                  COMPLEX (KIND=16) :: Yk, y(:)
                  
                  pi = acos(-1.d0)

                  OPEN(newunit=io, file="out/saida-6-11820812.out")
                  DO 10 k=0,N/2-1
                        Yk = Ykf(y, N, k) / (N/2)
                        freq = k / (N*dt)
                        WRITE(io, *) freq, (REAL(Yk)**2 + AIMAG(Yk)**2)
 10               END DO
                  CLOSE(io)

            END SUBROUTINE FOURRIER

            COMPLEX (KIND=16) FUNCTION Ykf (y, N, k)
                  IMPLICIT REAL (KIND=8) (a-h,o-z)
                  COMPLEX (KIND=16) :: y(:), zi
                  pi = ACOS(-1.d0)
                  zi = (0.d0, 1.d0)

                  Ykf = 0.d0
                  DO 10 j=0,N-1
                        Ykf  = Ykf +  y(j+1) * EXP(2.d0*pi*zi * k * j/N)
 10               END DO

                  RETURN
            END FUNCTION Ykf
            
      end module Frequencias