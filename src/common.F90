#include "fabm_driver.h"

module nflexpd_common

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
   real(rk) function SIT(aI,mu0hat,I,ThH)
!
! !DESCRIPTION:
! Here, the light limitation term (Pahlow and Oschlies MEPS 2013) is calculated. 
! This term also depends on T, because the growth rate depends on T.  
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                 :: aI,mu0hat,I, ThH 
!
! !REVISION HISTORY:
!  Original author(s): S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   SIT = 1.0 - exp( - aI * ThH * I / (mu0hat) ) !Eq.22 in K20
   return
 end function SIT
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! !IROUTINE: Nutrient uptake rate based on Optimal Uptake (OU) kinetics 
!
! !INTERFACE:
   real(rk) function vOU(N,Apot,Vpot)
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
   real(rk), intent(in)                :: N, Apot, Vpot
!
! !REVISION HISTORY:
!  Original author(s): S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! eq. 20 in Smith et al 2016
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
   real(rk) function vAff(N,f,Apot,Vpot)
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
   real(rk), intent(in)                :: N, f, Apot, Vpot
!
! !REVISION HISTORY:
!  Original author(s): S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! eq. 16 in Smith et al. 2016
   vAff = (1-f)*Vpot * f*Apot * N / ( (1-f)*Vpot + f*Apot*N ) !Eq.17 in K20
   return
 end function vAff
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
! !IROUTINE: Temperature dependence of rates for the FlexPFT model 
!
! !INTERFACE:
   real(rk) function FofT(tC)
!
! !DESCRIPTION:
! Here, Arrhenius type temperature dependence is calcuated, for a reference temperature, Tr, 
! and activation energy Ea [ J / mol ]. 
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: tC            ! temperature [ degrees C ]
   real(rk), parameter                 :: R = 8.31446   ! gas constant [ J /mol /K ]
   real(rk), parameter                 :: Ea=4.82e4        ! [ J / mol ]
   real(rk), parameter                 :: Tr=20.0          ! [ degrees C ]
!
! !REVISION HISTORY:
!  Original author(s):  S. Lan Smith, 20141213
!
!EOP
!-----------------------------------------------------------------------
!BOC
   FofT = exp( - (Ea/R)*( 1/(273.15+tC) - 1/(273.15+Tr) ) ) !Eq.28 in K20
   return
   end function FofT
!EOC
!-----------------------------------------------------------------------   


!-----------------------------------------------------------------------
!
! !IROUTINE: Fractional Day Length
!
! !INTERFACE:
   real(rk) function FDL(L,doy)
!
! !DESCRIPTION:
! Here, the sunrise and sunset are calculated based on latitude and day of year (doy)
! based on astronomical formulations
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in)    :: L,doy
   real(rk)                :: P,D,nom,denom
   INTEGER                 :: J
! !CONSTANTS   
   real(rk),parameter      :: pi=3.14159
   
!
! !REVISION HISTORY:
!  Original author(s):  O. Kerimoglu 27.11.2018
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !http://mathforum.org/library/drmath/view/56478.html
   !...
   !_Ecological Modeling_, volume 80 (1995) pp. 87-95, "A Model  Comparison for Daylength as a Function of Latitude and Day of the  Year." This article presented a model that apparently does a very good job of estimating the daylight - the error is less than one minute within 40 degrees of the equator, and less than seven minutes within  60 degrees and usually within two minutes for these latitudes.

   !D = daylength
   !L = latitude
   !J = day of the year
   
   J=floor(doy)

   P = asin(.39795*cos(.2163108 + 2*atan(.9671396*tan(.00860*(J-186)))))
   nom= sin(0.8333*pi/180) + sin(L*pi/180)*sin(P)
   denom= cos(L*pi/180)*cos(P)  
   D = 24 - (24/pi)*acos(nom/denom)
   
   !write(*,*)'L,J,D:',L,J,D
   
   FDL=D/24.
   
   return
   end function FDL
!EOC

end module
!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
