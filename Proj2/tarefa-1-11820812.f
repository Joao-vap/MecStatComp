      PROGRAM GENERATE
            IMPLICIT REAL (KIND=8) (a-h,l,o-z)
            DIMENSION grade(100, 3)
            COMMON /ONDA/ grade

c           grade
            ixs = 100
            its = 15

c           parametros
            c = 300.d0
            L = 1.d0
            dx = L / ixs
            r = 2.d0

            open(1, file='out/tarefa-1-11820812.out', status='unknown')

c           condicoes iniciais (bordas devem ser zero)
            DO 10 n = 1,2
                  DO 20 i=1,ixs
                        grade(i,n) = CONDINICIAIS(i*dx, L)
  20              END DO
                  write(1,fmt='(100(F20.4))') (grade(i,n), i=1,ixs)
  10        END DO

c           definir a evolucao temporal
            DO 40 n=3,its
                  CALL UPDATE(ixs, r)
                  write(1,fmt='(100(F20.8))') (grade(i,3), i=1,ixs)
  40        END DO

            close(1)

      END PROGRAM GENERATE

      SUBROUTINE UPDATE(ixs, r)
            IMPLICIT REAL (KIND=8) (a-h,l,o-z)
            DIMENSION grade(100, 3)
            COMMON /ONDA/ grade

c           paredes fixas
            grade(1, 3) = grade(1, 2)
            grade(ixs, 3) = grade(ixs, 2)

            r2 = r*r

c           nova linha (sem paredes)
            DO 10 i=2,ixs-1
                  grade(i, 3) = 2.d0*(1.d0-r2)*grade(i,2) +
     +                  r2*(grade(i+1,2)+grade(i-1,2))-grade(i,1)
 10         END DO

            grade(:,1) = grade(:,2)
            grade(:,2) = grade(:,3)

      END SUBROUTINE UPDATE

      FUNCTION CONDINICIAIS(x, L)
            IMPLICIT REAL (KIND=8) (a-h,l,o-z)
c           F0(x) =  exp[−(x − x0)^2 / σ^2]
c           σ = L/30
c           x0 = L/3
            CONDINICIAIS = EXP(-((x-L/3)**2)/(L/30)**2)
      END FUNCTION CONDINICIAIS
