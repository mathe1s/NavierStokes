MODULE EMPIRICAL_COEFFICIENTS

!****** Description: 
!****** This model is using to 
!****** calculate the viscosity of the 
!****** different species in function of  
!****** the temperature. 
!****** Author: Matheus Silva
!****** Update: 19/02/2018   
!Sutheland's formula
!mi(j)   = (1.458d-6)*(Tb(j)*20.0d0 + 293.15d0)**(3.0d0/2.0d0)/((Tb(j)*20.0d0+293.15d0)+110.4d0)/mi0
!paper:Transport Coefficients for the NASA Lewis Chemical Equilibrium Program
!author: Roger Svehla 1995

USE GLOBAL

CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!RETURNS THE NONDIMENSIONAL VALUE FOR THE VISCOSITY COEFFICIENT MI(I,J)
SUBROUTINE MOLECULAR_MI(T_inst)
IMPLICIT NONE
INTEGER :: I,J
REAL(KIND=IP),DIMENSION(1,1:imax,1:jmax)::T_inst


!nao pode ser calculado com a temperatura pertubadada vai dar zero!  

DO I=1,imax
    DO J=1,jmax
      MI(I,J)   = 1.0D0/MI0*(DEXP(A_MI*DLOG(Tref*T_inst(1,I,J)) + B_MI/Tref/T_inst(1,I,J) +& 
                  C_MI/Tref/T_inst(1,I,J)/Tref/T_inst(1,I,J) + D_MI))
      !!!viscosidade perturbada 
      !MI(I,J)   = MI(I,J)-MIb(I,J)
    END DO 
END DO  
END SUBROUTINE MOLECULAR_MI

!RETURNS THE NONDIMENSIONAL VALUE FOR THE THERMAL CONDUCTION K1(I,J) COEFFICIENT
subroutine k_coefficient(T_inst)
implicit none
INTEGER :: I,J
REAL(KIND=IP),DIMENSION(1,1:imax,1:jmax)::T_inst
DO I=1,imax
    DO J=1,jmax
            !instantaneus 
            K1(I,J)     = 1.0D0/K0*(DEXP(A_K*DLOG(Tref*T_inst(1,I,J)) + B_K/Tref/T_inst(1,I,J) + & 
                        C_K/Tref/Tref/T_inst(1,I,J)/T_inst(1,I,J) + D_K))
            !perturbate
            !K1(I,J)   = K1(I,J)-K1b(I,J)
    END DO 
END DO  
end subroutine k_coefficient


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!INSIDE THIS SUBROUTINE, THE USER CHOOSES THE CHEMICAL SPECIE AND GETS
!THE REFERENCE VALUE FOR THE VISCOSITY AND THERMAL CONDUCTION COEFFICIENTS
SUBROUTINE GET_COEFFICIENTS()
IMPLICIT NONE
INTEGER :: I,J
!CHOOSE HERE THE CHEMICAL SPECIE
CALL  SPECIE_N()

!DEFINING THE REFERENCE VALUES
k0     = DEXP(A_K *DLOG(1.0d0*Tref) + B_K /Tref  + C_K /Tref/Tref + D_K)
MI0    = DEXP(A_MI*DLOG(1.0D0*Tref) + B_MI/Tref  + C_MI/Tref/Tref + D_MI)

!RETURNS VALUES FOR THE BASE FLOW VERTICAL PROFILE OF VISCOSITY COEFFICIENT
   ! DO J=1,jmax
   !     mibase(j) = 1.0D0/MI0 * (dexp(A_mi*dlog(Tb(j)*Tref) + B_mi/Tb(j)/Tref + C_mi/Tb(j)/Tref/Tb(j)/Tref + D_mi))  
   ! END DO 

!REPLIES THE BASE FLOW VERTICAL PROFILE FOR THE INTIRE FLOW
!DO I=1,imax
!    DO J=1,jmax
!        MIB(I,J) = MiBase(J)
!    END DO 
!END DO 

!DO J=1,jmax    
!    DO I=1,imax
!K1B(I,J)=1.0D0/K0*(DEXP(A_K*DLOG(Tref*TEMPB(I,J))+B_K/Tref/TEMPB(I,J)+C_K/Tref/Tref/TEMPB(I,J)/TEMPB(I,J)+D_K))
!    END DO 
!END DO  

END SUBROUTINE GET_COEFFICIENTS 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BELOW, SOME CHEMICAL SPECIES TO BE USED 
!INSIDE THE SUBROUTINE "GET COEFFICIENTS"
!###### NITROGEN (N) ######
   !paper:Transport Coefficients for the NASA Çewis Chemical Equilibrium Program
   !author: Roger Svehla 1995
SUBROUTINE SPECIE_N()
IMPLICIT NONE
A_MI    = 0.62526577D0 
B_MI    =  -31.779652D0 
C_MI    = -1640.7913D0 
D_MI    = 1.7454992D0 
A_K 	= 0.85372829d0
B_K 	= 105.18665d0
C_K 	= -12229.753d0
D_K 	= 0.48299104d0
END SUBROUTINE SPECIE_N

END MODULE

