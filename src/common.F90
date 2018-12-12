#include "fabm_driver.h"

module PMbench_common

   use fabm_types
   use fabm_standard_variables

   implicit none
   
   public
   
   ! Aggregate diagnostics for e.g., carbon budgets.
   type (type_bulk_standard_variable),parameter :: total_PPR = type_bulk_standard_variable(name='total_PPR',units='mmolC/m^3/d',aggregate_variable=.true.)
   
   contains
   
!-----------------------------------------------------------------------
!
! !IROUTINE: Light limitation for the FlexPFT model 
!
! !INTERFACE:
   REALTYPE function SIT(aI,mu0hat,I,ThH,Tfac)
!
! !DESCRIPTION:
! Here, the light limitation term (Pahlow and Oschlies MEPS 2013) is calculated. 
! This term also depends on T, because the growth rate depends on T.  
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                 :: aI,mu0hat,I, ThH, Tfac 
!
! !REVISION HISTORY:
!  Original author(s): S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !dependence of growth rate on light (eq. 7 in Smith et al 2016)
   SIT = 1.0 - exp( - aI * ThH * I / (Tfac * mu0hat) )
   return
 end function SIT
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! !IROUTINE: Nutrient uptake rate based on Optimal Uptake (OU) kinetics 
!
! !INTERFACE:
   REALTYPE function vOU(N,Apot,Vpot)
!
! !DESCRIPTION:
!            The nutrient uptake rate calculated by Optimal Uptake (OU) kinetics assuming instantaneous 
!            physiological acclimation of the allocation of internal resources for nutrient uptake 
!            (Pahlow. MEPS, 2005; Smith et al. MEPS, 2009). In the FlexP model, based on the model of 
!            Pahlow and Oschlies (MEPS, 2013), VNhat is the potential maximum rate of nutrient uptake 
!            at the current ambient concentration, N. 
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: N, Apot, Vpot
!
! !REVISION HISTORY:
!  Original author(s): S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! eq. 20
   !Original:
   !vOU = Vpot * Apot * N / ( Vpot + 2*sqrt(Vpot*Apot*n) + Apot*N )
   vOU = Vpot * N / ( Vpot/Apot  + 2*sqrt(Vpot*N/Apot) + N )
   return
 end function vOU
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
! !IROUTINE: Nutrient uptake rate based on Affinity-based kinetics with a trade-off 
!
! !INTERFACE:
   REALTYPE function vAff(N,f,Apot,Vpot)
!
! !DESCRIPTION:
!            The nutrient uptake rate calculated by the affinity-based equation including the trade-off
!            postulated by Pahlow (MEPS, 2005) between affinity (initial slope) and maximum uptake rate. 
!            This function accepts a value for the allocation factor, f, (f_A as defined by Pahlow), allowing the 
!            uptake rate to be calculated for any allocation factor, whether it be optimal or not. 
!            Optimal Uptake (OU) kinetics assuming instantaneous 
!            In the FlexP model, based on the model of 
!            Pahlow and Oschlies (MEPS, 2013), VNhat is the potential maximum rate of nutrient uptake 
!            at the current ambient concentration, N. 
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: N, f, Apot, Vpot
!
! !REVISION HISTORY:
!  Original author(s): S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! eq. 16 in Smith et al. 2016
   vAff = (1-f)*Vpot * f*Apot * N / ( (1-f)*Vpot + f*Apot*N )
   return
 end function vAff
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
! !IROUTINE: Temperature dependence of rates for the FlexPFT model 
!
! !INTERFACE:
   REALTYPE function FofT(tC)
!
! !DESCRIPTION:
! Here, Arrhenius type temperature dependence is calcuated, for a reference temperature, Tr, 
! and activation energy Ea [ J / mol ]. 
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: tC            ! temperature [ degrees C ]
   REALTYPE, parameter                 :: R = 8.31446   ! gas constant [ J /mol /K ]
   REALTYPE, parameter                 :: Ea=4.82e4        ! [ J / mol ]
   REALTYPE, parameter                 :: Tr=20.0          ! [ degrees C ]
!
! !REVISION HISTORY:
!  Original author(s):  S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   FofT = exp( - (Ea/R)*( 1/(273.15+tC) - 1/(273.15+Tr) ) )
   return
   end function FofT
!EOC
!-----------------------------------------------------------------------   

end module

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
