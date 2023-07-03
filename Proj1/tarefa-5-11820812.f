      PROGRAM TIMEFOURRIER
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16), DIMENSION (1000) :: data
            COMMON /DAT/ data
            CHARACTER*200 :: filename, Nin(4)
            INTEGER :: N(4)
            REAL t(2)

            N =  [ 50, 100, 200, 400 ]
            Nin = ['050-11820812.in', '100-11820812.in',
     +            '200-11820812.in', '400-11820812.in' ]


            open(newunit=ioT, file="saidas/saida-5.1-11820812.out")
            open(newunit=ioTL, file="saidas/saida-5.2-11820812.out")
            write(ioT,*) 'N ', 'tempo'
            write(ioTL,*) 'N ', 'SqrtTempo'

            e = dtime( t )
            eaux = e

            DO 10 i=1,4

                  filename = "entradas/entrada-5-"//Nin(i)

                  OPEN(newunit=io, file=filename)
                  DO 20 j=1,1000
                        READ(io,*,end=1) time, data(j)
  20              END DO
  1               CONTINUE

                  CLOSE(io)

                  M = j
                  dt = time/(M-1)

                  DO 100 k=1,500
                        CALL FOURRIER(N(i), dt)
  100             END DO                   

                  e = dtime( t )
                  write(ioT,*) N(i), (e-eaux)/100
                  write(ioTL,*) N(i), SQRT((e-eaux)/100)
                  eaux = e

  10        END DO

      END PROGRAM TIMEFOURRIER

      SUBROUTINE FOURRIER (N, dt)
            IMPLICIT REAL (KIND=8) (a-h,o-z)
            COMPLEX (KIND=16) :: Yk, Ykf
            
            pi = acos(-1.d0)

            DO 10 k=0,N/2-1
                  Yk = Ykf(N, k) / (N/2)
                  freq = k / (N*dt)
  10        END DO

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