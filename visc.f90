MODULE EMPIRICAL_COEFFICIENTS

!****** Description: 
!****** This model is using to 
!****** calculate the viscosity of the 
!****** different species in function of  
!****** the temperature. 
!****** Author: Matheus Silva
!****** Update: 23/11/2017   
!Sutheland's formula
!mi(j)   = (1.458d-6)*(Tb(j)*20.0d0 + 293.15d0)**(3.0d0/2.0d0)/((Tb(j)*20.0d0+293.15d0)+110.4d0)/mi0
!paper:Transport Coefficients for the NASA Çewis Chemical Equilibrium Program
!author: Roger Svehla 1995

USE GLOBAL

CONTAINS

SUBROUTINE MOLECULAR_MI()

IMPLICIT NONE
INTEGER :: I,J

MI0   = DEXP(A_MI*DLOG(TEMP_0) + B_MI/TEMP_0 + C_MI/TEMP_0/TEMP_0 + D_MI)

DO I=I_S,D_S
    DO J=S_S,E_S
        MI2(I,J)   = DEXP(A_MI*DLOG(TEMP_0*TEMP(I,J)) + B_MI/TEMP_0/TEMP(I,J) + C_MI/TEMP_0/TEMP(I,J)/TEMP_0/TEMP(I,J) + D_MI)
        MI(I,J)    = MI2(I,J)/MI0
    END DO 
END DO  

END SUBROUTINE MOLECULAR_MI

SUBROUTINE INIT_MI()
IMPLICIT NONE
INTEGER :: I,J
CALL  VISC_N()
mib0   = dexp(A_mi*dlog(1.0D0*TEMP_0) + B_mi/1.0D0/TEMP_0 + C_mi/1.0D0/TEMP_0/1.0D0/TEMP_0 + D_mi) 

!tem que calcular o mi1 e mi2
!mi1 com a temperatura dimensional um, sabendo qual é T_0
!mib1   = dexp(A_mi*dlog(1.0D0*TEMP_0*T1) + B_mi/1.0D0/TEMP_0*T1 + C_mi/1.0D0/TEMP_0*T1/1.0D0/TEMP_0*T1 + D_mi) 

!mi2 com a temperatura dimensional dois, sabendo qual é T_0
!mib1   = dexp(A_mi*dlog(1.0D0*TEMP_0*T2) + B_mi/1.0D0/TEMP_0*T2 + C_mi/1.0D0/TEMP_0*T2/1.0D0/TEMP_0*T2 + D_mi) 

!Só basta calcular uma vez o perfil, para que dois Mib2 e mib?

!Agora tem que adimensionalizar essa visocidades com mi1

!mi_1=mi1/mi1=1
!mi_2=mi2/mi1

!agora com o perfil de temperatura q a gente tem, tb(j), calcula todas as outras viscoidades e 
!as adimensionaliza dividindo por mi_1

!DO J=S_S,E_S
! mib(j)   = dexp(A_mi*dlog(Tb(j)*TEMP_0) + B_mi/Tb(j)/TEMP_0 + C_mi/Tb(j)/TEMP_0/Tb(j)/TEMP_0 + D_mi)
!END DO

!e pronto vc tem o perfil de viscoidade adimensional 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mathues já gastamos tempo demais aqui, temos que andar pq se nao 
! vc se atrasa, ainda tem muitas coisas por fazer. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO J=S_S,E_S
    mib2(j)   = dexp(A_mi*dlog(Tb(j)*TEMP_0) + B_mi/Tb(j)/TEMP_0 + C_mi/Tb(j)/TEMP_0/Tb(j)/TEMP_0 + D_mi)
    mib(j)    = mib2(j)/mib0   
END DO 

DO I=I_S,D_S
    DO J=S_S,E_S
        MI(I,J) = MIB(J)
    END DO 
END DO 
END SUBROUTINE INIT_MI 

!###### NITROGEN ######
SUBROUTINE VISC_N()
IMPLICIT NONE
A_MI    = 0.82926975D0 
B_MI    =  405.82833D0 
C_MI    = -159002.42D0 
D_MI    = 0.17740763D0 

A_K = 0.82928303d0
b_k = 405.77643d0
c_k = -158950.37d0
d_k = 0.97751362d0
END SUBROUTINE VISC_N

END MODULE

