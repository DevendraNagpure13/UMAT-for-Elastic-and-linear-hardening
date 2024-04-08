      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4)
      
	  
C Elastic constants from the .inp file
      E=PROPS(1)
      PNU= PROPS(2)
      ALPHA=PROPS(3)
      ETREF=PROPS(4)
      
C Correction of Young's modulus for temperature
      E=E*(1-ALPHA*(TEMP-ETREF))

C Definining three terms for DDSDDE matrix elements
      TERM1= E/(1.0+PNU)/(1.0 - 2.0*PNU)
      TERM2= 1.0 - PNU
      TERM3= (1.0 -2.0*PNU)/2.0

C Initialise DDSDDE to make sure that DDSDDE elements are zero initially.

      DO I=1,NTENS
      DO J=1,NTENS 
	     DDSDDE(I,J)=0.0
      END DO 
      END DO

C Based on the matrix on right define DDSDDE elements
      DDSDDE(1,1)= TERM1*TERM2 
      DDSDDE(2,2)= TERM1* TERM2
      DDSDDE(3,3)= TERM1*TERM2 
      DDSDDE(4,4)= TERM1*TERM3
      DDSDDE(5,5)= TERM1*TERM3 
      DDSDDE(6,6)= TERM1*TERM3
      DDSDDE(1,2)= TERM1*PNU
      DDSDDE(1,3)= TERM1*PNU
      DDSDDE(2,3)= TERM1*PNU
      DDSDDE(2,1)= TERM1*PNU
      DDSDDE(3,1)= TERM1*PNU
      DDSDDE(3,2)= TERM1*PNU

C UPDATE STRESS BASED ON STRAIN INCREMENTS AND JACOBIAN

      DO I=1,NTENS
      DO J=1,NTENS
         STRESS(I) = STRESS(I) + DDSDDE(I,J)*DSTRAN(J)
      END DO
      END DO

      RETURN
      END