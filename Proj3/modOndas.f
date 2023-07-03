      module Ondas

      contains

            SUBROUTINE UPDATEFIXOSOLTO(grade, ixs, r)
                  IMPLICIT REAL (KIND=8) (a-h,l,o-z)
                  DIMENSION grade(:,:)

c                 parede fixas à esquerda
                  grade(1, 3) = grade(1, 2)

                  r2 = r*r

c                 nova linha (parede esquerda fixa)
                  DO 10 i=2,ixs-1
                  grade(i, 3) = 2.d0*(1.d0-r2)*grade(i,2) +
     +                  r2*(grade(i+1,2)+grade(i-1,2))-grade(i,1)
 10               END DO

c                 parede solta à direita
                  grade(ixs, 3) = grade(ixs-1, 3)

                  grade(:,1) = grade(:,2)
                  grade(:,2) = grade(:,3)

            END SUBROUTINE UPDATEFIXOSOLTO

            SUBROUTINE UPDATEFIXOFIXO(grade, ixs, r)
                  IMPLICIT REAL (KIND=8) (a-h,l,o-z)
                  DIMENSION grade(:,:)

c                 paredes fixas
                  grade(1, 3) = grade(1, 2)
                  grade(ixs, 3) = grade(ixs, 2)

                  r2 = r*r

c                 nova linha (sem paredes)
                  DO 10 i=2,ixs-1
                  grade(i, 3) = 2.d0*(1.d0-r2)*grade(i,2) +
     +                  r2*(grade(i+1,2)+grade(i-1,2))-grade(i,1)
 10               END DO

                  grade(:,1) = grade(:,2)
                  grade(:,2) = grade(:,3)

            END SUBROUTINE UPDATEFIXOFIXO

            FUNCTION GAUSSIAN(x, L)
            IMPLICIT REAL (KIND=8) (a-h,l,o-z)
c               F0(x) =  exp[−(x − x0)^2 / s^2]
c               s = L/30
c               x0 = L/2
            GAUSSIAN = EXP(-((x-2*L/3)**2)/(L/30)**2)
            END FUNCTION GAUSSIAN

      end module Ondas