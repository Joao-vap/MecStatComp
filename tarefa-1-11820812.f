      PROGRAM AUTOMATA
            INTEGER r(8), celulata(102,2), rule(3), nrule
            character(50) :: filename, lum(2)
            COMMON /RULE/ r
            COMMON /CELULATA/ celulata

            lum = (/'.','$'/)

            rule = (/51, 232, 254/)

            DO I=1,3

                  nrule = rule(I)
                  
                  WRITE(filename, '(A, I0, A)') 'validate/', nrule
                  OPEN(UNIT=1,FILE=filename,STATUS='UNKNOWN')

                  call setcelulata()
                  call setrule(nrule)

                  write(1,'(102A1)') (lum(celulata(k,1)+1), k=1,102)

                  DO J=1,50
                        CALL UPDATECELULATA()
                        write(1,'(102A1)')(lum(celulata(k,1)+1),k=1,102)
                  END DO

                  CLOSE(UNIT=1)                                

            END DO

      END PROGRAM AUTOMATA

      SUBROUTINE SETCELULATA()
            INTEGER celulata(102,2)
            COMMON /CELULATA/ celulata

            DO N = 1,102
                  r = rand()
                  if (r.gt.0.5) then
                        celulata(N,1)=1
                  else
                        celulata(N,1)=0
                  end if
            END DO

      END SUBROUTINE SETCELULATA

      SUBROUTINE UPDATECELULATA()
            INTEGER celulata(102,2)
            COMMON /CELULATA/ celulata

            DO I=2,101
                  CALL UPDATECELL(I)
            END DO

            celulata(1,2)=celulata(101,2)
            celulata(102,2)=celulata(2,2)

            DO I=1,102
                  celulata(I,1)=celulata(I,2)
            END DO

      END SUBROUTINE UPDATECELULATA

      SUBROUTINE UPDATECELL(id)
            INTEGER celulata(102,2), r(8), case(3), BINARYTOINT
            COMMON /CELULATA/ celulata
            COMMON /RULE/ r

            case = (/celulata(id-1,1),celulata(id,1),celulata(id+1,1)/)
            icase = BINARYTOINT(case)

            SELECT CASE (icase)
                  CASE (0)
                        celulata(id,2)=r(8)
                  CASE (1)
                        celulata(id,2)=r(7)
                  CASE (2)
                        celulata(id,2)=r(6)
                  CASE (3)
                        celulata(id,2)=r(5)
                  CASE (4)    
                        celulata(id,2)=r(4)
                  CASE (5)
                        celulata(id,2)=r(3)
                  CASE (6)
                        celulata(id,2)=r(2)
                  CASE (7)
                        celulata(id,2)=r(1)
                  CASE DEFAULT
                        write(*,*) 'error'
            END SELECT

      END SUBROUTINE UPDATECELL

      SUBROUTINE SETRULE(N)
            INTEGER r(8), B(8)
            COMMON /RULE/ r

            Naux = N
            DO I=1,8
                  B(I)=MOD(Naux,2)
                  Naux = Naux/2
            END DO

            DO I=1,8
                  R(9-I)=B(I)
            END DO
      END SUBROUTINE SETRULE

      INTEGER FUNCTION BINARYTOINT(B)
            INTEGER B(3)
            N=0
            DO I=1,3
                  N=N+B(I)*2**(3-I)
            END DO
            BINARYTOINT=N
      END FUNCTION BINARYTOINT
      
