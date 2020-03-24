#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: nflexpd - DOQ (Dynamically Optimized Quota) phytoplankton component
! Original Authors: O. Kerimoglu 20181204
!
! !INTERFACE:
   module nflexpd_phy_DOQ
!
! !DESCRIPTION:
!
! !USES:
   use lambert
   use nflexpd_common
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_nflexpd_phy_DOQ
!     Variable identifiers
      type (type_state_variable_id)        :: id_phyC,id_phyN
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
!  Here, the nflexpd_phy_DOQ namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_nflexpd_phy_DOQ), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk), parameter :: N2C_RF = 16._rk/106._rk !Redfield N:C ratio
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
   call self%register_state_variable(self%id_phyC,'C','mmolC/m^3','bound-C concentration',0.0_rk,minimum=0.0_rk,vertical_movement=w_phy)
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
! !IROUTINE: Right hand sides of the nflexpd model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_nflexpd_phy_DOQ), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: din,phyC,phyN,parW,par,par_dm,Ld
   real(rk)                   :: ThetaHat,vNhat,muIhat
   real(rk)                   :: Q,Theta,fV,fA,Rchl,I_zero,ZINT,valSIT
   real(rk)                   :: vN,Vhat_fNT
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
   _GET_(self%id_phyC,phyC)  ! phytoplankton-C
   _GET_(self%id_phyN,phyN)  ! phytoplankton-N
   _GET_(self%id_din,din)    ! nutrients
   
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
   ! in cmo, mu0=phy%V0*ft/(phy%rdl + phy%daylen)
   I_zero = self%zetaChl * self%RMchl * Tfac / (Ld*self%aI)   ! Threshold irradiance
   if( self%theta_opt ) then
     if( par_dm .gt. I_zero ) then !in cmo: .and. (mu0>0.0)
       !argument for the Lambert's W function
       larg = (1.0 + self%RMchl * Tfac/(Ld*self%mu0hat*Tfac)) * exp(1.0 + self%aI*par_dm/(self%mu0hat*Tfac*self%zetaChl))
       !larg=min(1e38,larg) !larg can explode if aI is too large compared to mu0hat*zetaChl
       ! eq. 8 in Smith et al 2016
       ThetaHat = 1.0/self%zetaChl + ( 1.0 -  WAPR(larg, 0, 0) ) * self%mu0hat*Tfac/(self%aI*par_dm)
       ThetaHat=max(0.1,ThetaHat) !  a small positive value 
       !if (ThetaHat .lt. 0.09)then
       !  write(*,*)'larg, self%aI*par_dm, self%mu0hat*Tfac*self%zetaChl',larg, self%aI*par_dm, self%mu0hat*Tfac*self%zetaChl
         !write(*,*)'ThetaHat,larg,WAPR',ThetaHat,larg,WAPR(larg,0,0)
       !end if
     else
       !write(*,*)'par_dm,I_0',par_dm,I_zero
       ThetaHat = 0.1  !  a small positive value
       !in cmo: ThetaHat=0.0 
     end if
   else
     ThetaHat = self%TheHat_fixed
   end if
   
   ! Light limited growth rate (eq. 6 in Smith et al 2016)
   if (par_dm .gt. I_zero) then !in cmo: .and. (mu0>0.0)
     valSIT=SIT(self%aI,self%mu0hat,par,ThetaHat,Tfac)
     muIhat = self%mu0hat * Tfac * valSIT 
   else
     valSIT=0.0
     muIhat=0.0
   end if
   
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
   
   !Dynamically calculated quota is needed for calculating some rates below
   Q= phyN/phyC
   
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
   PProd = (1.0-self%kexc) * max(0.0, mu*phyC )  ! PP [ mmolC / m3 / s ]

   !Total Chl content per C in Cell (eq. 10 in Smith et al 2016)
   Theta= (1 - self%Q0 / 2 / Q - fV) * ThetaHat
   
   !Uptake rate
   vN = fV * vNhat
   
   !Excretion:
   exc = self%kexc * mu * phyN
   ! Mortality
   mort=self%M0p * Tfac * PhyN**2
   
   !Calculate fluxes between pools
   f_din_phy = vN * phyC
   f_phy_detn =       self%Mpart  * mort 
   f_phy_don = (1.0 - self%Mpart) * mort + exc
   
   
   ! Set temporal derivatives
   _SET_ODE_(self%id_phyC, mu*phyC - f_phy_don/Q - f_phy_detn/Q)
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

   end module nflexpd_phy_DOQ

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
