      PROGRAM REVOLUCAO

C     SIMULAREMOS O AGREGAMENTO DE PARTICULAS EM UMA SEMENTE LIVRE
C     SE MOVENDO NO PLANO E SE AGREGANDO A PARTICULAS AO CONTATO
      
      PARAMETER(N=200)
      LOGICAL INFIELD, FULL
      INTEGER PARTICLES((2*N)**2+1,2),COORD(-N-1:N+1,-N-1:N+1),
     + REF(2),IMAX, ESCRITA(-N:N,-N:N), Raio, Raux
      COMMON /C/ COORD, NFREEPARTICLES
      COMMON /S/ PARTICLES, REF, NLOCKEDPARTICLES
      character(5) :: lum(2)

C     initialize the random number generator
      CALL SRAND(998784433)

      NLOCKEDPARTICLES=0
      NFREEPARTICLES=0
      IMAX = 0
      Raio = 0
      Raux = 0
      Lum = (/'.','$'/)

C     ARQUIVO PARA CALCULO DE DIMENSÃO FRACTAL
      OPEN(UNIT=2,FILE='revolucao01D.dat',STATUS='UNKNOWN')

C     INCIALIZAREMOS NOSSAS CONDIÇÕES INICIAIS DE PARTICULAS NO PLANO 
C     DE COOR DE FORMA SIMILAR AO CAMPO MINADO, OU SEJA:
C                 0  0  0  0  0  0  0  0  0  0
C                 0  1  1  1  0  0  1  1  1  0
C                 0  1  10 1  0  1  2  11 1  0
C                 0  1  1  1  0  1  11 2  1  0
C                 0  0  0  0  0  1  1  1  0  0 
C                 0  0  0  0  0  0  0  0  0  0
C     ONDE 10 É SOMADO AO PONTO DA PARTICULA E 1 AOS ARREDORES

      CALL SETUPFIELD()

C     INICIALIZAMOS NOSSO VETOR DE PARTICULAS E SEMENTE INICIAL
C     NOTE QUE USAREMOS REF PARA MOVIMENTAR O BLOCO SEM ATUALIZAR
C     O VETOR DE PARTICULAS, APENAS SOMANDO VETORES
      CALL SETUPPARTICLES()

C     COMEÇAMOS A AGREGAR E MOVIMENTAR AS PARTICULAS
      CALL AGGREGATE(Raio, Raux)

      IF (NLOCKEDPARTICLES<NFREEPARTICLES) FULL=.FALSE.

      DO WHILE (INFIELD().AND..NOT.FULL.AND.NLOCKEDPARTICLES<10000)

C     QUANDO NÃO HOUVER MAIS PARTICULAS ENTORNO, MOVIMENTAMOS O BLOCO
            CALL MOVEBLOCK()
            
C     AGREGATE VAI CHECAR POR PARTICULAS NO ARREDOR, REMOVE-LAS DA GRADE
C     E ADICIONA-LAS AO BLOCO DE PARTICULAS
      
            CALL AGGREGATE(Raio, Raux)

            IF (NLOCKEDPARTICLES<NFREEPARTICLES) FULL=.FALSE.
            IMAX=IMAX+1

      END DO

      OPEN(UNIT=1,FILE='revolucao0.dat',STATUS='UNKNOWN')
      M = N
      DO I =-M,M
            DO J=-M,M
                  ESCRITA(I,J)=0
            END DO
      END DO
      DO I=1,NLOCKEDPARTICLES
            ESCRITA(PARTICLES(I,1),PARTICLES(I,2))=1
      END DO
      write(*,*) "ok"
      DO I=-M,M
            write(1,'(402A1)') (lum(ESCRITA(I,J)+1),J=-M,M)
      END DO
      close(1)

      END PROGRAM REVOLUCAO

      SUBROUTINE SETUPFIELD()
            PARAMETER(N=200, p=0.0)
            INTEGER COORD(-N-1:N+1,-N-1:N+1)
            COMMON /C/ COORD, NFREEPARTICLES

            DO I=-N-1,N+1
                  DO J=-N-1,N+1
                        COORD(I,J)=0
                  END DO
            END DO
          
            DO I=-N,N
                  DO J=-N,N
                        IF(RAND() < p) THEN
                                  CALL ADDPARTICLE(I, J)
                                  NFREEPARTICLES=NFREEPARTICLES+1
                        END IF
                  END DO
            END DO
            
      END SUBROUTINE SETUPFIELD

      SUBROUTINE SETUPPARTICLES()
            PARAMETER(N=200)
            INTEGER PARTICLES ((2*N)**2+1,2), REF(2)
            COMMON /S/ PARTICLES, REF, NLOCKEDPARTICLES

            DO I=1,N
                  DO J=1,2
                        PARTICLES(I,J)=0
                  END DO
            END DO

            PARTICLES(1,1)=0
            PARTICLES(1,2)=0
            PARTICLES(2,1)=10101010
            PARTICLES(2,2)=10101010

            REF(1)=0
            REF(2)=0
            NLOCKEDPARTICLES=1

      END SUBROUTINE SETUPPARTICLES

      SUBROUTINE ADDPARTICLE(I, J)
            PARAMETER(N=200)
            INTEGER I, J
            INTEGER COORD(-N-1:N+1,-N-1:N+1)
            COMMON /C/ COORD, NFREEPARTICLES

            COORD(I,J)=COORD(I,J)+9
            DO K=-1,1
                  DO L=-1,1
                        COORD(I+K,J+L)=COORD(I+K,J+L)+1
                  END DO           
            END DO

      END SUBROUTINE ADDPARTICLE

      SUBROUTINE REMOVEPARTICLE(I, J)
            PARAMETER(N=200)
            INTEGER I, J
            INTEGER COORD(-N-1:N+1,-N-1:N+1)
            COMMON /C/ COORD, NFREEPARTICLES

            COORD(I,J)=COORD(I,J)-9
            DO K=-1,1
                  DO L=-1,1
                        COORD(I+K,J+L)=COORD(I+K,J+L)-1
                  END DO
            END DO

      END SUBROUTINE REMOVEPARTICLE

      SUBROUTINE LOCKPARTICLE(I, J)
            PARAMETER(N=200)
            INTEGER PARTICLES((2*N)**2+1,2), REF(2), NLOCKEDPARTICLES
            COMMON /S/ PARTICLES, REF, NLOCKEDPARTICLES

            PARTICLES(NLOCKEDPARTICLES+1,1)=I
            PARTICLES(NLOCKEDPARTICLES+1,2)=J
            PARTICLES(NLOCKEDPARTICLES+2,1)=10101010
            PARTICLES(NLOCKEDPARTICLES+2,2)=10101010

            NLOCKEDPARTICLES=NLOCKEDPARTICLES+1

      END SUBROUTINE LOCKPARTICLE

      LOGICAL FUNCTION INFIELD()
            PARAMETER(N=200)
            INTEGER PARTICLES ((2*N)**2+1,2), REF(2)
            COMMON /S/ PARTICLES, REF, NLOCKEDPARTICLES

            INFIELD=.TRUE.
            IF (ABS(REF(1)) > N) INFIELD=.FALSE.
            IF (ABS(REF(2)) > N) INFIELD=.FALSE.

            RETURN
      END FUNCTION INFIELD

      SUBROUTINE MOVEBLOCK()
            PARAMETER(N=200)
            INTEGER PARTICLES ((2*N)**2+1,2), REF(2)
            COMMON /S/ PARTICLES, REF, NLOCKEDPARTICLES

            REF(1)=REF(1)+floor(RAND()*3)-1
            REF(2)=REF(2)+floor(RAND()*3)-1

      END SUBROUTINE MOVEBLOCK

      SUBROUTINE AGGREGATE(Raio, Raux)
            PARAMETER(N=200)
            INTEGER PARTICLES((2*N)**2+1, 2), REF(2),
     +              COORD(-N-1:N+1, -N-1:N+1), X, Y, C1, C2, Raio, Raux
            COMMON /S/ PARTICLES, REF, NLOCKEDPARTICLES
            COMMON /C/ COORD, NFREEPARTICLES

            I = 1
            DO 10 WHILE (PARTICLES(I,1) /= 10101010)

                  X = PARTICLES(I,1) + REF(1)
                  Y = PARTICLES(I,2) + REF(2)

                  IF (ABS(X) > N .OR. ABS(Y) > N) GOTO 100
                  IF (COORD(X,Y) == 0) THEN
                        GOTO 100
                  ELSE
                        DO 20 K=-1,1
                              DO 30 L=-1,1
                              IF (K == 0 .AND. L == 0) GOTO 200

                              C1 = X+K
                              C2 = Y+L

                              IF (COORD(C1,C2) >= 10) THEN

                              IF(x**2+y**2>Raio**2)Raio=(x**2+y**2)**0.5
                        
                              CALL REMOVEPARTICLE(C1,C2)
                              CALL LOCKPARTICLE(C1-REF(1),C2-REF(2))

                              if (Raio - Raux > 10) then
                                    write(2,*) Raio, NLOCKEDPARTICLES
                                    Raux = Raio
                              end if

                              END IF
 200                          CONTINUE
 30                           END DO
 20                     END DO
                  END IF

 100              CONTINUE
                  I = I + 1
 10         END DO

      END SUBROUTINE AGGREGATE




            