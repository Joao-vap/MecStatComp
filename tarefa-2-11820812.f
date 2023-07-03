      PROGRAM DLA
C           this program computes the 2d dla - diffusion limited aggregation
            PARAMETER (N=100, iparticles=5000)
            INTEGER COORD(-N:N,-N:N), FREE(2), radius, r0, r2, part
            character(5) :: lum(2)
            COMMON /COORD/ COORD
            LOGICAL GONE, TOUCH

            lum = (/'.','$'/)
            open(2, file='dla2d/radius1.dat')

            radius = 1
            part = 0

C           initialize the seed
            COORD(0,0) = 1

C           initialize the random number generator
            CALL SRAND(987339365)

C           loop over all particles
            DO 100 I=1, iparticles

                  IF (radius >= 99) GOTO 100

C                 initialize a new free particle in a random position
                  CALL RANDOM_POSITION(FREE, radius, r0)

C                 loop its random walk
                  DO 10 WHILE (.NOT. GONE(FREE, r0))
                        CALL WALK(FREE)
                        IF (TOUCH(FREE)) THEN

c                             add the particle to the cluster
                              COORD(FREE(1),FREE(2)) = 1

c                             update the radius and count
                              part = part + 1
                              r2 = FREE(1)**2 + FREE(2)**2
                              if (r2 > radius**2) then 
                                    radius=Ceiling(r2**0.5)
                                    write(2,*) radius, part
                              end if
                              GOTO 20

                        END IF
 10               END DO
 20               CONTINUE

 100        END DO

C           print the cluster
            open(1, file='dla2d/dla1.dat')
            DO 30 l=-N,N
                  write(1,'(201A1)') (lum(COORD(k,l)+1),k=-N,N)
 30         END DO

            close(1)
            close(2)

      END PROGRAM DLA

C     this subroutine generates a random position for a new particle
C     the position is a random point on a square of side = radius*5
C     the center of the square is the origin
      SUBROUTINE RANDOM_POSITION(FREE, radius, r0)
            INTEGER FREE(2), radius, r0, sgnx, sgny

c           isg is the sign of the random number  
            sgnx= sign(1,FlOOR(RAND()-0.5))
            sgny = sign(1,FlOOR(RAND()-0.5))
c           random points on a square of side radius*5 at Origin
c           except for the center square of side radius+5
            FREE(1)=sgnx*((radius+5)+CEILING(RAND()*radius*3))
            FREE(2)=sgny*((radius+5)+CEILING(RAND()*radius*3))
c           r0 is the distance of the particle from the origin
            r0 = max(abs(FREE(1)), abs(FREE(2)))

            RETURN
      END SUBROUTINE RANDOM_POSITION

C     this subroutine makes a random walk of a particle
C     the particle moves in a random direction by a random distance
      SUBROUTINE WALK(FREE)
            INTEGER FREE(2)
c           walk in a random direction -1, 0, 1 in x and -1, 0, 1 in y
            FREE(1) = FREE(1) + FLOOR(RAND()*3)-1
            FREE(2) = FREE(2) + FLOOR(RAND()*3)-1
            RETURN
      END SUBROUTINE WALK

C     this function checks if a particle has gone out of the square
      LOGICAL FUNCTION GONE(FREE, r0)
            INTEGER FREE(2), r0
            GONE = .FALSE.
            IF (ABS(FREE(1)) .GT. r0*4) GONE = .TRUE.
            IF (ABS(FREE(2)) .GT. r0*4) GONE = .TRUE.
            RETURN
      END FUNCTION GONE

C     this function checks if a particle has touched the cluster
      LOGICAL FUNCTION TOUCH(FREE)
            PARAMETER (N=100)
            INTEGER FREE(2), COORD(-N:N,-N:N)
            COMMON /COORD/ COORD

            TOUCH = .FALSE.
C           CHECK BORDERS
            IF (ABS(FREE(1))+1 >= N .OR. ABS(FREE(2)) + 1 >= N) GOTO 30
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