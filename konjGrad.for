C Bestimmung der Lösung eines linearen Gleichungssystems
C mit SYMETRISCHER POSITIV DEFINITER MATRIX
	PROGRAM konjGrad
C Vereinbarungen
C NMAX - Beschränkung der Dimension
C n - RANG
C A - geg. Matrix
C b - rechte Seite
C x - Startvektor
	IMPLICIT NONE
	EXTERNAL LGSIN,KONGRAD
	INTEGER NMAX,n,i,control,kmax
	PARAMETER (NMAX=50)
	REAL A(NMAX,NMAX),b(NMAX),x(NMAX),eps
C
C Anweisung
C Einlesen der Startwerte
	CALL LGSIN(A,b,x,n,kmax,eps)
C Ausgabe des Startwertes
	DO 10 i=1,n
		PRINT'(/A,I2,A,F15.10/)','vor: x(',i,')=',x(i)
 10	CONTINUE
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Aufruf der SUBROUTINE, x enthält hoffentlich die Lösung
	CALL KONGRAD(A,x,b,n,kmax,eps)
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Kontrollausgabe der Endwerte
	DO 20 i=1,n
		PRINT'(/A,I2,A,F15.10/)','nach: x(',i,')=',x(i)
 20	CONTINUE
 	PRINT'(/A,I4/)','benoetigte Iterationsschritte: ',kmax
C	Dient nur der Vorführung (system("PAUSE");-Ersatz)
 	READ*,control
	END



C ********************************************
C Unterprogramm Verfahren der KONjungierten GRADienten
	SUBROUTINE KONGRAD( A , x , b , n , kmax , eps )
C Vereinbarung
	IMPLICIT NONE
	EXTERNAL MAP,SKP,NORM
	INTEGER n,NMAX,i,j,kmax
	PARAMETER (NMAX=50)
	REAL A(NMAX,NMAX),x(NMAX),b(NMAX),P(NMAX),S(NMAX),R(NMAX),
     &alpha,betha,SKP,eps,NORM
C Anweisung
C Bestimmung von r0, 
C R enthält Startresiduum, P die Startrichtung
	j=0
	CALL MAP(A,x,R,n)
	DO 30 i=1,n
		R(i)=R(i)-b(i)
		P(i)=-R(i)
 30	CONTINUE
C Hauptschleife
40	CONTINUE
		  CALL MAP(A,P,S,n)
		  alpha=SKP(R,R,n)/SKP(P,S,n)
		  DO 50 i=1,n
			  x(i)=x(i)+alpha*P(i)
			  R(i)=R(i)+alpha*S(i)
 50	CONTINUE
		  betha=SKP(R,S,n)/SKP(P,S,n)
		  DO 60 i=1,n
			  P(i)=-R(i)+betha*P(i)
 60	CONTINUE
	DO 61 i=1,n
		PRINT'(A,I2,A,F15.10)','innen: x(',i,')=',x(i)
 61	CONTINUE
 	j=j+1
 	IF ( NORM(R,n) .LE. eps .OR. j .GE. kmax )GOTO 62
	GOTO 40
 62	CONTINUE
 	kmax=j
 	END



C ********************************************
C Unterprogramm MAtrixProdukt
	SUBROUTINE MAP( A , P , S , n )
C Vereinbarung
	IMPLICIT NONE
	INTEGER n,NMAX,i,j
	PARAMETER (NMAX=50)
	REAL A(NMAX,NMAX),P(NMAX),S(NMAX)
C Anweisung
	DO 70 i=1,n
	S(i)=0
 70	CONTINUE
	DO 80 i=1,n
		DO 90 j=1,n
			S(i)=S(i)+A(i,j)*P(j)
 90	CONTINUE
 80	CONTINUE
	END



C ********************************************
C Funktion SKalarProdukt
	REAL FUNCTION SKP( Vek1 , Vek2 , n )
C Vereinbarung
	IMPLICIT NONE
	INTEGER n,NMAX,i
	PARAMETER (NMAX=50)
	REAL Vek1(NMAX),Vek2(NMAX),a
C Anweisung
	a=0
	DO 110 i=1,n
		a=a+Vek1(i)*Vek2(i)
 110	CONTINUE
	SKP=a
	END


C ********************************************
C Funktion Vektornorm (euklidische Norm)
	REAL FUNCTION NORM( Vek , n )
C Vereinbarung
	IMPLICIT NONE
	INTEGER n,NMAX,i
	PARAMETER (NMAX=50)
	REAL Vek(NMAX),a
C Anweisung
	a=0
	DO 210 i=1,n
		a=a+Vek(i)*Vek(i)
 210	CONTINUE
	NORM=SQRT(a)
	END



C *******************************************
C Einlesen des LGS und restlicher Daten
      SUBROUTINE LGSIN( A , b , x , n , kmax , eps )
C Vereinbarungsteil
      IMPLICIT NONE
      INTEGER NMAX,n,i,j,UNITNR,kmax
      PARAMETER(NMAX=50,UNITNR=20)
      REAL b(NMAX),x(NMAX),A(NMAX,NMAX),eps
C Anweisungsteil
C Datenfile oeffnen
C
	OPEN(UNIT=UNITNR,FILE='LGS.dat',STATUS='OLD')
C
C Dimension der Matrix bzw. des Vektors einlesen
C
      READ(UNITNR,*) n
      WRITE(*,'(/A,I2/)') 'eingelesene Dimension des Problems: ', n
C
C Maximale Anzahl der Interationen
C
      READ(UNITNR,*) kmax
      WRITE(*,'(/A,I4/)') 'Maximale Anzahl der Iterationen: ', kmax
C
C Genauigkeit
C
      READ(UNITNR,*) eps
      WRITE(*,'(/A,F15.10/)') 'gewuenschte Genauigkeit: ', eps
C
C Werte der Matrixkomponenten einlesen 
C
      WRITE(*,'(/A/)') 'Eingelesene Matrix:'
      DO 120 i=1,n
          READ(UNITNR,*,ERR=99,END=99) (A(i,j), j=1,n)      
          WRITE(*,*) (A(i,j), j=1,n)      
 120  CONTINUE
C
C Werte der Vektorkomponenten einlesen
C
      WRITE(*,'(//A/)') 'Eingelesener b-Vektor:'
      DO 130 i=1,n
          READ(UNITNR,*,ERR=99,END=99) b(i)      
          WRITE(*,'((F15.10))') b(i)
 130  CONTINUE
C
	WRITE(*,'(//A/)') 'Eingelesener x0-Vektor:'
      DO 140 i=1,n
          READ(UNITNR,*,ERR=99,END=99) x(i)      
          WRITE(*,'((F15.10))') x(i)
 140  CONTINUE
C
C Datenfile schliessen
C
      CLOSE(UNITNR)
C
	GOTO 199
C
C Fehler-Behandlung
C     
 99   WRITE(*,*) 'Fehler auf dem Datenfile! Programmabbruch'
C
 199  CONTINUE
 	END  
