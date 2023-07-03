      PROGRAM RAIO
C           this program computes the 2d dla - diffusion limited aggregation
            PARAMETER (N=100, iparticles=5000)
            INTEGER COORD(0:N,0:N), FREE(2), heigth
            character(5) :: lum(2)
            COMMON /COORD/ COORD
            LOGICAL GONE, TOUCH

            lum = (/'.','$'/)
            heigth = 0

C           initialize the seed on all the line y=0
            DO 10 I=0,N
                  COORD(0,I) = 1
 10         END DO

C           initialize the random number generator
            CALL SRAND(998784433)

C           loop over all particles
            DO 100 I=1, iparticles

C                 initialize a new free particle in a random position
                  CALL RANDOM_POSITION(FREE, N, heigth)

                  IF (heigth >= N - 1) GOTO 100

C                 loop its random walk while not gone or not touched
                  DO 30 WHILE (.NOT. GONE(FREE))
                        CALL WALK(FREE)
                        IF (TOUCH(FREE)) THEN
                              IF (heigth .LT. FREE(1)) heigth = FREE(1)
                              COORD(FREE(1),FREE(2)) = 1
                              GOTO 20
                        END IF
 30               END DO
 20               CONTINUE

 100        END DO

C           print the cluster
            open(1, file='dlaRaioexp3.dat')
            DO 40 l=0,N
                  write(1,'(101A1)') (lum(COORD(N-l,N-k)+1),k=0,N)
 40         END DO

            close(1)

      END PROGRAM RAIO

C     this subroutine generates a random position for a new particle
C     the position is a random point on a square of side = radius*5
C     the center of the square is the origin
      SUBROUTINE RANDOM_POSITION(FREE, N, heigth)
            INTEGER FREE(2), heigth
  
            FREE(1) = heigth + 15
            FREE(2) = CEILING(RAND()*N)

            RETURN
      END SUBROUTINE RANDOM_POSITION

C     this function checks if a particle has gone out of the square
      LOGICAL FUNCTION GONE(FREE)
            INTEGER FREE(2)
            PARAMETER (N=100)

            GONE = .FALSE.
            IF (ABS(FREE(1)) .GT. N*5/4) GONE = .TRUE.
            IF (ABS(FREE(2) - N/2) .GT. (N*5/4)) GONE = .TRUE.

            RETURN
      END FUNCTION GONE

C     this subroutine makes a random walk of a particle
C     the particle moves in a random direction by a random distance
      SUBROUTINE WALK(FREE)
            INTEGER FREE(2)

c           walk in a random direction -1, 0, 1 in x and -1, 0 in y
            FREE(1) = FREE(1) + FLOOR(RAND()*3)-1
            FREE(2) = FREE(2) + FLOOR(RAND()*3)-1

            RETURN
      END SUBROUTINE WALK

C     this function checks if a particle has touched the cluster
      LOGICAL FUNCTION TOUCH(FREE)
            PARAMETER (N=100)
            INTEGER FREE(2), COORD(0:N,0:N)
            COMMON /COORD/ COORD

            TOUCH = .FALSE.
C           CHECK BORDERS
            IF (FREE(1) + 1 > N .OR. FREE(1) - 1 < 0) GOTO 30
            IF (FREE(2) + 1 > N .OR. FREE(2) - 1 < 0) GOTO 30
            DO 10 I=-1,1
                  DO 20 J=-1,1
                        IF(COORD(FREE(1)+I,FREE(2)+J).EQ.1) THEN
                              TOUCH=.TRUE.
                              GOTO 30
                        END IF
 20               END DO     
 10         END DO
 30         CONTINUE
            RETURN
      
      END FUNCTION TOUCH