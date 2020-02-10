#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NflexPD - phytoplankton component
! Original Authors: S. Lan Smith, 2014-12-09
! FABM implementation: O. Kerimoglu 20181122
!
! !INTERFACE:
   module NflexPD_phy
!
! !DESCRIPTION:
!
! !USES:
   use lambert
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_NflexPD_phy
!     Variable identifiers
      type (type_state_variable_id)        :: id_phyN
      type (type_state_variable_id)        :: id_din,id_don,id_detn
      type (type_dependency_id)            :: id_parW,id_temp,id_par_dmean
      type (type_horizontal_dependency_id) :: id_FDL
      type (type_diagnostic_variable_id)   :: id_Q,id_Chl2C,id_mu,id_fV,id_fA,id_ThetaHat
      type (type_diagnostic_variable_id)   :: id_PPR
      
      
!     Model parameters
      real(rk) :: kc,w_phy
      real(rk) :: zetaN,zetaChl,kexc,M0p,Mpart,RMChl
      real(rk) :: mu0hat,aI
      real(rk) :: A0hat,V0hat,Q0
      real(rk) :: fA_fixed,fV_fixed,TheHat_fixed
      logical  :: fV_opt,fA_opt,Theta_opt
      real(rk) :: dic_per_n

      contains

      procedure :: initialize
      procedure :: do
   end type
   
   type (type_bulk_standard_variable),parameter :: total_PPR = type_bulk_standard_variable(name='total_PPR',units='mmolC/m^3/d',aggregate_variable=.true.)
   
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the NflexPD_phy namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_NflexPD_phy), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk)            :: w_phy
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day and are converted here to values per second.
   ! General:
   call self%get_parameter(self%kc,   'kc',   'm2 mmolN-1','specific light extinction',               default=0.03_rk)
   call self%get_parameter(self%w_phy,       'w_phy',  'm d-1',    'vertical velocity (<0 for sinking)',      default=-1.0_rk, scale_factor=d_per_s)
   !optimality switches
   call self%get_parameter(self%Theta_opt, 'Theta_opt','-', 'whether to optimize theta', default=.false.)
   call self%get_parameter(self%fA_opt, 'fA_opt','-', 'whether to optimize fA', default=.false.)
   call self%get_parameter(self%fV_opt, 'fV_opt','-', 'whether to optimize fV', default=.false.)
   !light-related
   call self%get_parameter(self%TheHat_fixed, 'TheHat_fixed','gChl molC-1', 'Theta_Hat to use when Theta_opt=false', default=0.6_rk)
   call self%get_parameter(self%RMchl, 'RMchl','d-1', 'loss rate of chlorophyll', default=0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%mu0hat, 'mu0hat','d-1', 'max. potential growth rate', default=5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%aI, 'aI','(m^2 E-1 molC gChl-1)', 'Chl-specific slope of the PI curve', default=1.0_rk) ! really /mol or /micromol?
   !nutrient-related
   call self%get_parameter(self%fA_fixed, 'fA_fixed','-', 'fA to use when fa_opt=false', default=0.5_rk)
   call self%get_parameter(self%fV_fixed, 'fV_fixed','-', 'fV to use when fv_opt=false', default=0.25_rk)
   call self%get_parameter(self%Q0, 'Q0','molN molC-1', 'Subsistence cell quota', default=0.039_rk)
   call self%get_parameter(self%V0hat, 'V0hat','molN molC-1 d-1', 'Potential maximum uptake rate', default=5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%A0hat, 'A0hat','m3 mmolC-1 d-1', 'Potential maximum nutrient affinity', default=0.15_rk,scale_factor=d_per_s)
   !mortality/loss/respiration
   call self%get_parameter(self%zetaN, 'zetaN','molC molN-1', 'C-cost of N uptake', default=0.6_rk)
   call self%get_parameter(self%zetaChl, 'zetaChl','molC gChl-1', 'C-cost of Chlorophyll synthesis', default=0.8_rk)
   call self%get_parameter(self%kexc,  'kexc',  '-',    'excreted fraction of primary production',                          default=0.01_rk)
   call self%get_parameter(self%M0p, 'M0p', 'm3 molN-1 d-1', 'sp. quad. mortality rate',              default=0.1_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Mpart, 'Mpart', '-',   'part of the mortality that goes to detritus',default=0.5_rk)
   

   ! Register state variables
   call self%register_state_variable(self%id_phyN,'N','mmolN/m^3','bound-N concentration',0.0_rk,minimum=0.0_rk,vertical_movement=w_phy, specific_light_extinction=self%kc)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phyN)

   ! Register dependencies on external state variables
   call self%register_state_dependency(self%id_din, 'din',   'mmolN/m^3','dissolved inorganic nitrogen')
   call self%register_state_dependency(self%id_don, 'don','mmolN/m^3','dissolved organic nitrogen')
   call self%register_state_dependency(self%id_detN,'detN','mmolN/m^3','detrital nitrogen')

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_Q, 'Q','molN/molC',    'cellular nitrogen Quota',           &
                                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Chl2C, 'Chl2C','gChl/molC',    'cellular chlorophyll content',           &
                                     output=output_instantaneous)                                     
   call self%register_diagnostic_variable(self%id_mu, 'mu','/d',    'net sp. growth rate',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fV, 'fV','-',    'fV',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fA, 'fA','-',    'fA',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_ThetaHat, 'ThetaHat','-', 'ThetaHat',           &
                                     output=output_time_step_averaged)
                                     
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmolC/m^3/d','Primary production rate',      &
                                     output=output_time_step_averaged)
   call self%add_to_aggregate_variable(total_PPR,self%id_PPR)

   ! Register environmental dependencies
   call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_par_dmean, 'PAR_dmean','E/m^2/d','photosynthetically active radiation, daily averaged')
   call self%register_horizontal_dependency(self%id_FDL, 'FDL','-',       'fractional day length')
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   
   
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of the NflexPD model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_NflexPD_phy), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: din,phyN,parW,par,par_dm,Ld
   real(rk)                   :: ThetaHat,vNhat,muIhat
   real(rk)                   :: Q,Theta,fV,fA,Rchl,I_zero,ZINT
   real                       :: larg !argument to WAPR(real(4),0,0) in lambert.f90
   real(rk)                   :: tC,Tfac
   real(rk)                   :: mu,exc,mort,Pprod
   real(rk)                   :: f_din_phy,f_phy_don,f_phy_detn
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_phyN,phyN)         ! phytoplankton
   _GET_(self%id_din,din) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_(self%id_parW,parW)             ! local photosynthetically active radiation
   par=parW* 4.6 * 1e-6   !molE/m2/s
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   _GET_(self%id_par_dmean,par_dm) !in molE/m2/d
   par_dm=par_dm/secs_pr_day !convert to molE/m2/s
   
   if ( par_dm .lt. 0.0 ) then
     par_dm=par
   end if
   
   !get Ld (fractional day length)
   _GET_HORIZONTAL_(self%id_FDL,Ld)
   
   _GET_(self%id_temp,tC) ! temperature in Celcius
   
   !Calculate intermediate terms:
   !Temperature factor 
   Tfac = FofT(tC)
   
   ! Primary production
   ! Optimization of ThetaHat (optimal Chl content in the chloroplasts)
   I_zero = self%zetaChl * self%RMchl * Tfac / (Ld*self%aI)   ! Threshold irradiance
   
   !why is this at the end in the original code?
   if( self%theta_opt ) then
     if( par .gt. I_zero ) then
       !argument for the Lambert's W function
       larg = (1.0 + self%RMchl * Tfac/(Ld*self%mu0hat*Tfac)) * exp(1.0 + self%aI*par_dm/(self%mu0hat*Tfac*self%zetaChl))
       ! eq. 8 in Smith et al 2016
       ThetaHat = 1.0/self%zetaChl + ( 1.0 -  WAPR(larg, 0, 0) ) * self%mu0hat*Tfac/(self%aI*par_dm)
       ThetaHat=max(0.01,ThetaHat)
     else
       ThetaHat = 0.01  !  a small positive value 
     end if
   else
     ThetaHat = self%TheHat_fixed
   end if
   
   ! Light limited growth rate (eq. 6 in Smith et al 2016)
   muIhat = self%mu0hat * Tfac * SIT(self%aI,self%mu0hat,par,ThetaHat,Tfac)
   
   !Optimal allocation for affinity vs max. uptake
   if( self%fA_opt ) then
     ! eq. 17 in Smith et al 2016) 
     fA = 1.0 / ( 1.0 + sqrt(self%A0hat * din /(Tfac * self%V0hat)) ) 
   else
     fA =  self%fA_fixed
   end if
   
   ! T-dependence only for V0, not for A0 (as suggested by M. Pahlow)
   vNhat = vAff( din, fA, self%A0hat, self%V0hat * Tfac )
   
   ! Alternative way of calculating vNhat
   !vNhat2=vOU(din,self%A0hat,self%V0hat * Tfac)
   
   !Optimization of fV (synthesis vs nut. uptake)
   !Intermediate term  in brackets that appears in Smith et al 2016, eqs. 13 & 14
   ZINT = (self%zetaN + muIhat/vNhat) * self%Q0 / 2.0
   !write(*,'(A,4F12.5)')'  (phy) ZINT, muIhat/vNhat:',ZINT,muIhat/vNhat
   
   if( self%fV_opt .and.  par .gt. I_zero ) then
     ! eq. 13  in Smith et al 2016
     fV = (-1.0 + sqrt(1.0 + 1.0 / ZINT) ) * (self%Q0 / 2.0) * muIhat / vNhat
   else
     fV = self%fV_fixed 
   end if
   
   !!$ ***  Calculating the optimal cell quota, based on the term ZINT, as calculated above
   if( self%fV_opt ) then
     ! eq. 14 in Smith et al 2016
     Q = ( 1.0 + sqrt(1.0 + 1.0/ZINT) )*(self%Q0/2.0)
   else
     Q = 1.0/6.67 !Almost Redfield? (106/16=6.625)
   end if
!!$ *** To correctly apply the Balanced Growth Assumption with this Non-adaptive model *** 
!!$ * * * Calculate QN, the N cell quota [mol N / mol C],
!!$ * * * based on the Balanced Growth Assumption ( V = mu*Q  <=>  Q = V/mu )
!!$ *** Note:  This does NOT assume optimal resource allocation ***  
!!$ * * * IF, fA and fV have both been optimized, this should give the same result as the caculation
!!$ * * * based on the term ZINT. 
!!$      Q(ci) = ( (Q0/2.0)*muIhat + fV*VNhat ) / ( (1.0-fV)*muIhat - fV*zetaN*VNhat )    

!!$      dQdNbyQ(ci) = 0.0

   ! Losses due to Chlorophyll
   ! eq. 26 in Smith et al 2016
   Rchl = (muIhat + self%RMchl*Tfac) * ( 1 - fV - self%Q0/(2.0*Q) ) * self%zetaChl * ThetaHat
   
   !  Net specific growth rate, assuming instantantaneous optimal resource allocation. 
      !  Either equation gives the same result, provided fA, fV and QN have been optimized.  
      !  The term with ZINT accounts for the cost of N assimilation, but not for chl maintenance. 
      !  mu = muIhat * ( 1 + 2*( ZINT - sqrt(ZINT*(1 + ZINT)) ) )           - Rchl 
!!$      mu = muIhat*(1.0 + 2.0*(ZINT - sqrt(ZINT*(1.0+ZINT))) ) - Rchl
   ! eq. 5 in Pahlow and Oschlies 2013 (-Rchl)
   mu = muIhat * ( 1 - fV - self%Q0/(2.0*Q) ) - self%zetaN*fV*vNhat - Rchl ![/s]
   
   !Just for the diagnostics:
   !Primary production rate:
   PProd = (1.0-self%kexc) * max(0.0, mu*phyN / Q )  ! PP [ mmolC / m3 / s ]

   !Total Chl content per C in Cell (eq. 10 in Smith et al 2016)
   Theta= (1 - self%Q0 / 2 / Q - fV)* Q
   
   !Excretion:
   exc = self%kexc * mu * phyN
   ! Mortality
   mort=self%M0p * Tfac * PhyN**2
   
   !Calculate fluxes between pools
   f_din_phy = mu * phyN
   f_phy_detn =       self%Mpart  * mort 
   f_phy_don = (1.0 - self%Mpart) * mort + exc
   
   
   ! Set temporal derivatives
   _SET_ODE_(self%id_phyN, f_din_phy - f_phy_don - f_phy_detn)

   ! If externally maintained dim,dom und det pools are coupled:
   _SET_ODE_(self%id_din, -f_din_phy)
   _SET_ODE_(self%id_don,  f_phy_don)
   _SET_ODE_(self%id_detN, f_phy_detn)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_Q, Q)
   _SET_DIAGNOSTIC_(self%id_Chl2C, Theta)
   _SET_DIAGNOSTIC_(self%id_fV, fV)
   _SET_DIAGNOSTIC_(self%id_fA, fA)
   _SET_DIAGNOSTIC_(self%id_mu, mu * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_ThetaHat, ThetaHat) 
   _SET_DIAGNOSTIC_(self%id_PPR, PProd*secs_pr_day) !*s_p_d such that output is in d-1

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC


!-----------------------------------------------------------------------
!
! !IROUTINE: Light limitation for the FlexPFT model 
!
! !INTERFACE:
   real(rk) function SIT(aI,mu0hat,I,ThH,Tfac)
!
! !DESCRIPTION:
! Here, the light limitation term (Pahlow and Oschlies MEPS 2013) is calculated. 
! This term also depends on T, because the growth rate depends on T.  
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                 :: aI,mu0hat,I, ThH, Tfac 
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
   ! eq. 20
   !Original:
   !vOU = Vpot * Apot * N / ( Vpot + 2*sqrt(Vpot*Apot*n) + Apot*N )
   vOU = Vpot * N / ( Vpot/Apot  + 2*sqrt(Vpot*N/Apot) + N )
   return
 end function vOU
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
   vAff = (1-f)*Vpot * f*Apot * N / ( (1-f)*Vpot + f*Apot*N )
   return
 end function vAff
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
   FofT = exp( - (Ea/R)*( 1/(273.15+tC) - 1/(273.15+Tr) ) )
   return
   end function FofT
!EOC

   end module NflexPD_phy

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
