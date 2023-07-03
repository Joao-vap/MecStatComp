      PROGRAM PROJ3

c           Este programa utiliza dois modulos
c           O modulo Ondas contem as rotinas para a resolucao da equacao
c           de onda em um meio homogeneo e isotropo
c           O modulo Frequencias contem a rotina para a transformada de
c           fourrier

c           Para compilar este programa, deve-se utilizar o comando
c           gfortran NomeArquivo modFourrier.f modOndas.f exec.out

            use Ondas
            use Frequencias

            IMPLICIT REAL (KIND=8) (a-h,l,o-z)
            DIMENSION grade(100, 3)
            COMPLEX (KIND=16) y(2000)
            PARAMETER (ixs=100, its=2000, c=300.d0, L=1.d0, r=1.d0)

            dx = L / ixs

            open(1, file='out/tarefa-1-aux-11820812.out')

c           condicoes iniciais
            DO 10 n = 1,2
                  DO 20 i=1,ixs
                        grade(i,n) = GAUSSIAN(i*dx, L)
  20              END DO
                  grade(1,n) = 0.d0
                  grade(ixs,n) = 0.d0
                  write(1,fmt='(100(F12.4))') (grade(i,n), i=1,ixs)
                  y(n) = grade(25,n)
  10        END DO

c           definir a evolucao temporal (bordas fixas)
            DO 40 n=3,its
                  CALL UPDATEFIXOSOLTO(grade, ixs, r)
                  write(1,fmt='(100(F12.4))') (grade(i,3), i=1,ixs)
                  y(n) = grade(25,3)
  40        END DO

            close(1)

            N = its
            dt = r*dx/c

            write(*,*) 'N = ', N, 'dt = ', dt

c           Chamar a rotina de fourrier que escrevera em data.out

            CALL FOURRIER(y, N, dt)

      END PROGRAM PROJ3