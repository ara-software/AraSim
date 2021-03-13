      DOUBLE PRECISION FUNCTION DGAUSS(F,A,B,EPS)                               
      DOUBLE PRECISION F,A,B,EPS                                                
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST                   
      save
      LOGICAL MFLAG,RFLAG                                                       
      EXTERNAL F                                                                
C                                                                               
C     ******************************************************************        
C                                                                               
C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                            
C                                                                               
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF           
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER           
C     EPS.                                                                      
C                                                                               
C     ******************************************************************        
C                                                                               
      DATA W / 0.10122 85362 90376 259D0,                                       
     1         0.22238 10344 53374 471D0,                                       
     2         0.31370 66458 77887 287D0,                                       
     3         0.36268 37833 78361 983D0,                                       
     4         0.27152 45941 17540 949D-1,                                      
     5         0.62253 52393 86478 929D-1,                                      
     6         0.95158 51168 24927 848D-1,                                      
     7         0.12462 89712 55533 872D0,                                       
     8         0.14959 59888 16576 732D0,                                       
     9         0.16915 65193 95002 538D0,                                       
     A         0.18260 34150 44923 589D0,                                       
     B         0.18945 06104 55068 496D0/                                       
                                                                                
      DATA X / 0.96028 98564 97536 232D0,                                       
     1         0.79666 64774 13626 740D0,                                       
     2         0.52553 24099 16328 986D0,                                       
     3         0.18343 46424 95649 805D0,                                       
     4         0.98940 09349 91649 933D0,                                       
     5         0.94457 50230 73232 576D0,                                       
     6         0.86563 12023 87831 744D0,                                       
     7         0.75540 44083 55003 034D0,                                       
     8         0.61787 62444 02643 748D0,                                       
     9         0.45801 67776 57227 386D0,                                       
     A         0.28160 35507 79258 913D0,                                       
     B         0.95012 50983 76374 402D-1/                                      
C                                                                               
C     ******************************************************************        
C                                                                               
C  START.                                                                       
      DGAUSS=0.0D0                                                              
      IF(B.EQ.A) RETURN                                                         
      CONST=0.005D0/(B-A)                                                       
      BB=A                                                                      
C                                                                               
C  COMPUTATIONAL LOOP.                                                          
    1 AA=BB                                                                     
      BB=B                                                                      
    2    C1=0.5D0*(BB+AA)                                                       
         C2=0.5D0*(BB-AA)                                                       
         S8=0.0D0                                                               
         DO 3 I=1,4                                                             
            U=C2*X(I)                                                           
            S8=S8+W(I)*(F(C1+U)+F(C1-U))                                        
    3    CONTINUE                                                               
         S8=C2*S8                                                               
         S16=0.0D0                                                              
         DO 4 I=5,12                                                            
            U=C2*X(I)                                                           
            S16=S16+W(I)*(F(C1+U)+F(C1-U))                                      
    4    CONTINUE                                                               
         S16=C2*S16                                                             
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5                       
         BB=C1                                                                  
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2                              
      DGAUSS=0.0D0                                                              
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG)                                  
      IF(MFLAG) THEN                                                            
         IF(LGFILE.EQ.0) THEN                                                   
            WRITE(*,6)                                                          
         ELSE                                                                   
            WRITE(LGFILE,6)                                                     
         ENDIF                                                                  
      ENDIF                                                                     
      IF(.NOT. RFLAG) CALL ABEND                                                
      RETURN                                                                    
    5 DGAUSS=DGAUSS+S16                                                         
      IF(BB.NE.B) GO TO 1                                                       
      RETURN                                                                    
C                                                                               
    6 FORMAT( 4X, 'FUNCTION DGAUSS ... TOO HIGH ACCURACY REQUIRED')             
      END                                                                       

c -------------------------------------------------------------------

          SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)                        
          save
          PARAMETER(KOUNTE  =  27)                                    
          CHARACTER*6         ERCODE,   CODE(KOUNTE)                            
          LOGICAL             MFLAG,    RFLAG                                   
          INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)                  
          DATA      LOGF      /  0  /                                           
          DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 100, 100 /           
          DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 100, 100 /           
          DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 100, 100 /           
          DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 100, 100 /           
          DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 100, 100 /           
          DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 100, 100 /           
          DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 100, 100 /           
          DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 100, 100 /           
          DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 100, 100 /           
          DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 100, 100 /           
          DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 100, 100 /           
          DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 100, 100 /           
          DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 100, 100 /           
          DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 100, 100 /           
          DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 100, 100 /           
          DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 100, 100 /           
          DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 100, 100 /           
          DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 100, 100 /           
          DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 100, 100 /           
          DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 100, 100 /           
          DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 100, 100 /           
          DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 100,   0 /           
          DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 100,   0 /           
          DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 100,   0 /           
          DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 100,   0 /           
          DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 100, 100 /           
          DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 100, 100 /           
          LOGF  =  LGFILE                                                       
          IF(ERCODE .EQ. ' ')  THEN                                             
             L  =  0                                                            
          ELSE                                                                  
             DO 10  L = 1, 6                                                    
                IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12                            
  10            CONTINUE                                                        
  12         CONTINUE                                                           
          ENDIF                                                                 
          DO 14     I  =  1, KOUNTE                                             
             IF(L .EQ. 0)  GOTO 13                                              
             IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14                         
  13         KNTM(I)  =  LIMITM                                                 
             KNTR(I)  =  LIMITR                                                 
  14         CONTINUE                                                           
            RETURN                                                                
          ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)                                  
          LOG  =  LOGF                                                          
          DO 20     I  =  1, KOUNTE                                             
             IF(ERCODE .EQ. CODE(I))  GOTO 21                                   
  20         CONTINUE                                                           
          WRITE(*,1000)  ERCODE                                                 
          CALL ABEND                                                            
          RETURN                                                                
  21      RFLAG  =  KNTR(I) .GE. 1                                              
          IF(RFLAG  .AND.  (KNTR(I) .LT. 100))  KNTR(I)  =  KNTR(I) - 1         
          MFLAG  =  KNTM(I) .GE. 1                                              
          IF(MFLAG  .AND.  (KNTM(I) .LT. 100))  KNTM(I)  =  KNTM(I) - 1         
          IF(.NOT. RFLAG)  THEN                                                 
             IF(LOGF .LT. 1)  THEN                                              
                WRITE(*,1001)  CODE(I)                                          
             ELSE                                                               
                WRITE(LOGF,1001)  CODE(I)                                       
             ENDIF                                                              
          ENDIF                                                                 
          IF(MFLAG .AND. RFLAG)  THEN                                           
             IF(LOGF .LT. 1)  THEN                                              
                WRITE(*,1002)  CODE(I)                                          
             ELSE                                                               
                WRITE(LOGF,1002)  CODE(I)                                       
             ENDIF                                                              
          ENDIF                                                                 
          RETURN                                                                
1000      FORMAT(' KERNLIB LIBRARY ERROR. ' /                                   
     +           ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',                 
     +           ' ERROR MONITOR. RUN ABORTED.')                                
1001      FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',               
     +           'CONDITION ',A6)                                               
1002      FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)                    
          END                                                                   
      SUBROUTINE ABEND                                                          
      save
C                                                                               
C CERN PROGLIB# Z035    ABEND           .VERSION KERNVAX  1.10  811126          
                                                                                
      STOP '*** ABEND ***'                                                      
      END                                                                       
