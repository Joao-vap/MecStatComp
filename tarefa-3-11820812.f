      PROGRAM DLA3D
C           this program computes the 2d dla - diffusion limited aggregation
            PARAMETER (N=100, iparticles=100000)
            INTEGER COORD(-N:N,-N:N,-N:N), FREE(3), radius, r0, r2, part
            character(5) :: lum(2)
            COMMON /COORD/ COORD
            LOGICAL GONE, TOUCH
      
            lum = (/'.','$'/)
            open(1, file='dla3d/radius03D.dat')
      
            radius = 1
            part = 0
      
C           initialize the seed
            COORD(0,0,0) = 1
      
C           initialize the random number generator
            CALL SRAND(987991650)
      
C           loop over all particles
            DO 100 I=1, iparticles
      
                  IF (radius >= 99) GOTO 100
      
C                 initialize a new free particle in a random position
                  CALL RANDOM_POSITION(FREE, radius, r0)
      
C                 loop its random walk
                  DO 10 WHILE (.NOT. GONE(FREE, r0))
                        CALL WALK(FREE)
                        IF (TOUCH(FREE)) THEN
                              COORD(FREE(1),FREE(2),FREE(3)) = 1
                              part = part + 1
                              r2 = FREE(1)**2+FREE(2)**2+FREE(3)**2
                              IF(r2.GT.radius**2) radius=ceiling(r**0.5)
                              GOTO 20
                        END IF
 10               END DO
 20               CONTINUE

                  write(1,*) radius, part
      
 100        END DO
      
            close(1)
      
      END PROGRAM DLA3D
      
C     this subroutine generates a random position for a new particle
C     the position is a random point on a square of side = radius*5
C     the center of the square is the origin
      SUBROUTINE RANDOM_POSITION(FREE, radius, r0)
            INTEGER FREE(3), radius, r0, sgnx, sgny, sgnz
      
c           isg is the sign of the random number  
            sgnx= sign(1,FlOOR(RAND()-0.5))
            sgny = sign(1,FlOOR(RAND()-0.5))
            sgnz = sign(1,FlOOR(RAND()-0.5))

c           random points on a square of side radius*5 at Origin
c           except for the center square of side radius+5
            FREE(1)=sgnx*((radius+5)+CEILING(RAND()*radius*3))
            FREE(2)=sgny*((radius+5)+CEILING(RAND()*radius*3))
            FREE(3)=sgnz*((radius+5)+CEILING(RAND()*radius*3))
c           r0 is the distance of the particle from the origin
            r0 = max(abs(FREE(1)), abs(FREE(2)), abs(FREE(3)))
      
            RETURN
      END SUBROUTINE RANDOM_POSITION
      
C     this subroutine makes a random walk of a particle
C     the particle moves in a random direction by a random distance
      SUBROUTINE WALK(FREE)
            INTEGER FREE(3)
c           walk in a random direction -1, 0, 1 in x and -1, 0, 1 in y
            FREE(1) = FREE(1) + FLOOR(RAND()*3)-1
            FREE(2) = FREE(2) + FLOOR(RAND()*3)-1
            FREE(3) = FREE(3) + FLOOR(RAND()*3)-1
            RETURN
      END SUBROUTINE WALK
      
C     this function checks if a particle has gone out of the square
      LOGICAL FUNCTION GONE(FREE, r0)
            INTEGER FREE(3), r0
            GONE = .FALSE.
            IF (ABS(FREE(1)) .GT. r0*4) GONE = .TRUE.
            IF (ABS(FREE(2)) .GT. r0*4) GONE = .TRUE.
            IF (ABS(FREE(3)) .GT. r0*4) GONE = .TRUE.
            RETURN
      END FUNCTION GONE
      
C     this function checks if a particle has touched the cluster
      LOGICAL FUNCTION TOUCH(FREE)
            PARAMETER (N=100)
            INTEGER FREE(3), COORD(-N:N,-N:N,-N:N)
            COMMON /COORD/ COORD
      
            TOUCH = .FALSE.
C           CHECK BORDERS
      IF(ABS(FREE(1))>=N.OR.ABS(FREE(2))>=N.OR.ABS(FREE(3))>=N) GOTO 40
            DO 10 I=-1,1
            DO 20 J=-1,1
            DO 30 K=-1,1
                  IF(COORD(FREE(1)+I,FREE(2)+J,FREE(3)+K).EQ.1) THEN
                        TOUCH=.TRUE.
                        GOTO 40
                  END IF
 30         END DO
 20         END DO     
 10         END DO
 40         CONTINUE
            RETURN
            
      END FUNCTION TOUCH