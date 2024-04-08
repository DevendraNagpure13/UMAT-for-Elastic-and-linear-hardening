      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      REAL STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4)
      REAL E, PNU
      REAL C(6,6)
      REAL SIGMA(6), SIGMA_D(6), N_VECTOR(6)
      REAL d_eps(6)
      REAL h, d_lambda
      REAL YIELD_STRENGTH, EFF_STRESS, F
      REAL d_stress(6),d_epsp(6), epsp(6)
      real DDSDDE(6,6)
      integer i
      E= PROPS(1)
      PNU= PROPS(2)
      h=PROPS(3)
      SIGMA= STRESS
      d_eps= DSTRAN
      epsp= STATEV
      YIELD_STRENGTH = PROPS(4)
      ! for abaqus sub routine need to remove this do loop because this main program will be run as subroutine in abaqus
      ! the loop will be decided by time step

      call stiffness_matrix(E,PNU,C)
      CALL DEVIATORIC_STRESS(SIGMA, SIGMA_D)
      CALL EFFECTIVE_STRESS(SIGMA,EFF_STRESS)
      CALL F_FACTOR(YIELD_STRENGTH, EFF_STRESS,F)
      CALL NORMAL_VECTOR(SIGMA_D,EFF_STRESS,N_VECTOR)
      
      if (F .LT. 0) then
          d_lambda= 0
          DDSDDE= C
      else
          CALL PLASTIC_MUL(N_VECTOR,d_eps,C,h,d_lambda)
          CALL jacob_explicit(C,N_VECTOR,h,DDSDDE)
      end if
      
      CALL increments(C,d_eps,d_lambda,N_VECTOR,d_stress,d_epsp)
      call jacob_explicit(C,N_VECTOR,h,DDSDDE)
      
      ! Update values 
      SIGMA = SIGMA+ d_stress
      STRESS= SIGMA
      epsp= epsp + d_epsp
      STATEV= epsp
      !YIELD_STRENGTH= YIELD_STRENGTH + (h*d_lambda)
      
      RETURN
      END
      
!C 100 determine stiffness matrix      
      subroutine stiffness_matrix(E,PNU,C)
      real E,PNU
      real C(6,6)
      integer i,j
      do i= 1,6
          do j= 1,6
              C(i,j)= 0
          end do
      end do
      T1= E/((1+PNU)*(1-2*PNU))
      T2= 1 - PNU
      T3= (1-2*PNU)/2
      do i= 1,3
          C(i,i)= T2
          C(i+3,i+3)= T3
      end do
      C(1,2)=PNU      ! assigning rest of the term
      C(1,3)=PNU
      C(2,1)=PNU
      C(2,3)=PNU
      C(3,1)=PNU
      C(3,2)=PNU
      C= C*T1
      
      return 
      end subroutine stiffness_matrix
!C 100 determine stiffness matrix 

C effective stresss subroutine copied from file 200_cal_effective_stress      
      SUBROUTINE EFFECTIVE_STRESS(S,EFF_STRESS)
      REAL S(6)  ! S= state of stress given by abaqus as a variable name STRESS
      T1= ((S(1)-S(2))**2 + (S(2)-S(3))**2 + (S(3)-S(1))**2)/2.0   ! equation is writen in parts for better readability 
      T2= 3.0*(S(4)**2 + S(5)**2 + S(6)**2)
      EFF_STRESS = SQRT(T1+T2)
      RETURN
      END

      SUBROUTINE F_FACTOR(YIELD_STRENGTH, EFFECTIVE_STRESS,F) 
      REAL YIELD_STRENGTH, EFFECTIVE_STRESS
      F= EFFECTIVE_STRESS - YIELD_STRENGTH 
      RETURN
      END
      ! f factor is a value of a yield funtion at a perticular state of stress
      ! f factor should always be less then zero 'f<0' 
      ! f factor is used to check whether curve is entered in plastic zone  
C 200 subroutine to calculate 'effective stress' and 'f factor' of given stress of state
      
C 210 subroutine to cal deviatoric stress       
      SUBROUTINE DEVIATORIC_STRESS(SIGMA, SIGMA_D)
      IMPLICIT NONE
      REAL SIGMA(6), SIGMA_D(6)
      REAL SIGMA_MEAN
      INTEGER I
      SIGMA_MEAN = (SIGMA(1)+SIGMA(2)+SIGMA(3))/3
      
      DO 5 I= 1,3
          SIGMA_D(I)= SIGMA(I)-SIGMA_MEAN
5     END DO
      
      DO 10 I= 4,6
          SIGMA_D(I)= SIGMA(I)
10    END DO 
      RETURN 
      END SUBROUTINE DEVIATORIC_STRESS
C 210 subroutine to cal deviatoric stress          
      
C 220 subroutine to cal norma vector 
      SUBROUTINE NORMAL_VECTOR(SIGMA_D,EFF_STRESS,N_VECTOR)
      REAL SIGMA_D(6)
      REAL EFF_STRESS
      REAL N_VECTOR(6)
      N_VECTOR= (3/2)*(SIGMA_D/EFF_STRESS)
      RETURN
      END SUBROUTINE NORMAL_VECTOR

C 230 subroutine to cal plastic multipier
      SUBROUTINE PLASTIC_MUL(N_VECTOR,d_eps,C,h,d_lambda)
      REAL N_VECTOR(6), d_eps(6), NVec_mat(6,1), d_eps_mat(6,1)
      REAL C(6,6)
      REAL h, h_mat(1,1)
      REAL d_lambda, d_lambda_mat(1,1)
      real T1(6,1),T2(6,1),N_VEC_TRAN(1,6)
      integer i
      do i=1,6
          NVec_mat(i,1)= N_VECTOR(i)
          d_eps_mat(i,1)= d_eps(i)
      end do
      h_mat(1,1)= h
      T1= MATMUL(C,d_eps_mat)
      T2= MATMUL(C,NVec_mat)
      N_VEC_TRAN= TRANSPOSE(NVec_mat)
      
      d_lambda_mat=MATMUL(N_VEC_TRAN,T1)/(MATMUL(N_VEC_TRAN,T2) + h_mat)
      d_lambda= d_lambda_mat(1,1)
      RETURN
      END SUBROUTINE PLASTIC_MUL
C 230 subroutine to calculate plastic multipier 
      
C 250 subroutine to calculate stress increament, plastic strain increment
      subroutine increments(C,d_eps,d_lambda,N_VECTOR,d_stress,d_epsp)
      real C(6,6),d_eps(6),N_VECTOR(6), d_stress(6),d_epsp(6)
      real d_lambda
      real d_epse(6), d_epse_mat(6,1),d_stress_mat(6,1)
      integer i
      d_epsp= d_lambda*N_VECTOR
      d_epse= d_eps - d_epsp
      do i= 1,6
          d_epse_mat(i,1)=d_epse(i)
      end do
      d_stress_mat= matmul(C,d_epse_mat)
      do i= 1,6
          d_stress(i)= d_stress_mat(i,1)
      end do
      return
      end subroutine increments
C 250 subroutine to calculate stress increament, plastic strain increment

C 260 subriutine to calculate jacobian matrix for explicite
      subroutine jacob_explicit(C,N_VECTOR,h,J)
      real C(6,6), N_VECTOR(6), h, J(6,6)
      real NVec_mat(6,1), Cn_mat(6,1), NVec_Trans(1,6), h_mat(1,1)
      real Cn_mat_trans(1,6),J1_mat(6,6),J2_mat(1,1),J2
      integer i
      h_mat(1,1)= h
      do i= 1,6
          NVec_mat(i,1)= N_VECTOR(i)
      end do
      NVec_Trans= transpose(NVec_mat)
      Cn_mat= matmul(C,NVec_mat)
      Cn_mat_trans= transpose(Cn_mat)
      J1_mat= matmul(Cn_mat,Cn_mat_trans)
      J2_mat= matmul(NVec_Trans,Cn_mat)+ h_mat
      J2= J2_mat(1,1)
      J= C- J1_mat/J2
      return
      end subroutine 
C 260 subriutine to calculate jacobian matrix for explicite  
      
          