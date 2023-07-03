      FUNCTION VIOLAO(x, L)
            IMPLICIT REAL (KIND=8) (a-h,l,o-z)
c           FROM 0 TO L/4 -> F(X) = X
c           FROM L/4 TO L -> F(X) = L/3 - x/3
            IF (x .LE. L/4) THEN
                  VIOLAO = x
            ELSE
                  VIOLAO = L/3 - x/3
            END IF
      END FUNCTION VIOLAO