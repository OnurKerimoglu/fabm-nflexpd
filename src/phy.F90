#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: nflexpd_phy: amorphous phytoplankton class that can be instanced as:
! 1) N-based models (dynQN=false), with:
!!! 1a) (FS) constant-stoichiometry without any flexibility, mimics Monod model
!!! 1b) (IA)(*) 'Instantaneously Acclimating' N:C based on flexible fA,fV (&ThetaHat)
!!! 1c) (IAfix)(**): 'Instantaneously Adjusting' N:C based on fixed fA,fV,ThetaHat 
! 2) N&C-based models (dynQN=true), with:
!!! 2a) (DA) 'Dynamically Acclimating' N:C, based on flexible fA,fV (&ThetaHat)
!!! 2b) (DAfix) 'Dynamically Adjusting' N:C based on flex fA,fV (&ThetaHat), mimics Droop model
!
!(*,**): 'Flex' and 'Control' models in  Smith et al., 2016, J. Plankton Res. 38, 977–992
! 
! Original Authors: O. Kerimoglu (v0:December 2018; v1: May 2020)
!
! !INTERFACE:
   module nflexpd_phy
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
   type,extends(type_base_model),public :: type_nflexpd_phy
!     Variable identifiers
      type (type_state_variable_id)        :: id_phyC,id_phyN
      type (type_state_variable_id)        :: id_din,id_don,id_doc,id_detn,id_detc
      type (type_dependency_id)            :: id_parW,id_temp,id_par_dmean,id_depth
      type (type_horizontal_dependency_id) :: id_FDL
      type (type_diagnostic_variable_id)   :: id_Q,id_d_phyC,id_Chl,id_Chl2C,id_fV,id_fA,id_ThetaHat
      type (type_diagnostic_variable_id)   :: id_PPR,id_fdinphy_sp,id_mu,id_respN,id_respChl
      type (type_diagnostic_variable_id)   :: id_fQ,id_fNmonod,id_fN,id_fL,id_Tfac
      type(type_diagnostic_variable_id)    :: id_fphydoc,id_fphydon,id_fphydetc,id_fphydetn
      
!     Model parameters
      real(rk) :: kc,w_phy,mindin
      real(rk) :: zetaN,zetaChl,M0p,Mpart,RMChl
      real(rk) :: mu0hat,aI
      real(rk) :: A0hat,V0hat,Q0,Qmax,KN_monod
      real(rk) :: fA_fixed,fV_fixed,TheHat_fixed,Q_fixed
      logical  :: dynQN,fV_opt,fA_opt,Theta_opt,mimic_Monod
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
!  Here, the nflexpd_phy namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_nflexpd_phy), intent(inout), target :: self
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
   call self%get_parameter(self%mindin,       'min_din',  'mmolN m-3',    'when provided, minimum din concentration that allows growth and uptake',      default=0.0_rk)
   !general switches
   call self%get_parameter(self%dynQN, 'dynQN','-', 'whether dynamically resolve QN', default=.true.)
   !optimality switches
   call self%get_parameter(self%Theta_opt, 'Theta_opt','-', 'whether to optimize theta', default=.false.)
   call self%get_parameter(self%fA_opt, 'fA_opt','-', 'whether to optimize fA', default=.false.)
   call self%get_parameter(self%fV_opt, 'fV_opt','-', 'whether to optimize fV', default=.false.)
   call self%get_parameter(self%mimic_Monod, 'mimic_Monod','-', 'whether to mimic Monod model', default=.false.)
   !light-related
   call self%get_parameter(self%TheHat_fixed, 'TheHat_fixed','gChl molC-1', 'Theta_Hat to use when Theta_opt=false', default=0.6_rk)
   call self%get_parameter(self%RMchl, 'RMchl','d-1', 'loss rate of chlorophyll', default=0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%mu0hat, 'mu0hat','d-1', 'max. potential growth rate', default=5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%aI, 'aI','(m^2 E-1 molC gChl-1)', 'Chl-specific slope of the PI curve', default=1.0_rk) ! really /mol or /micromol?
   !nutrient-related
   call self%get_parameter(self%fA_fixed, 'fA_fixed','-', 'fA to use when fa_opt=false', default=0.5_rk)
   call self%get_parameter(self%fV_fixed, 'fV_fixed','-', 'fV to use when fv_opt=false', default=0.25_rk)
   call self%get_parameter(self%Q_fixed, 'Q_fixed','-', 'Q to use when provided, dynQN=false and fV_opt=false', default=-1.0_rk)
   call self%get_parameter(self%Qmax, 'Qmax','molN molC-1', 'Maximum cell quota', default=0.3_rk)
   call self%get_parameter(self%Q0, 'Q0','molN molC-1', 'Subsistence cell quota', default=0.039_rk)
   call self%get_parameter(self%V0hat, 'V0hat','molN molC-1 d-1', 'Potential maximum uptake rate', default=5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%A0hat, 'A0hat','m3 mmolC-1 d-1', 'Potential maximum nutrient affinity', default=0.15_rk,scale_factor=d_per_s)
   call self%get_parameter(self%KN_monod, 'KN_monod','mmolN m-3', 'Half saturation constant for growth [when Monod model is mimicked]', default=-1.0_rk)
   
   !consistency checks
   if (self%Q_fixed .gt. 0.0 ) then
     !assume that user wants to mimic Monod model, check if other options are consistent:
     self%mimic_Monod=.true.
     if (self%dynQN) then
       call self%fatal_error('phy.F90/initialize:','for '//trim(self%name)// ' a valid Q_fixed was provided that signals intention to mimic Monod model, but dynQN needs to bet set to .false.')
     end if
     if (self%fV_opt) then
       call self%fatal_error('phy.F90/initialize:','for '//trim(self%name)// ' a valid Q_fixed was provided that signals intention to mimic Monod model, but fV_opt needs to be set to .false.')
     end if
   end if
     
   if (self%mimic_Monod) then
     if (self%dynQN) then
       call self%fatal_error('phy.F90/initialize:','for '//trim(self%name)// ' mimic_Monod and dynQN cannot be simultaneously true')
     end if
     if (self%fV_opt) then
       call self%fatal_error('phy.F90/initialize:','for '//trim(self%name)// ' mimic_Monod and fV_opt cannot be simultaneously true')
     end if
     if (self%Q_fixed .lt. 0.0_rk) then
       call self%fatal_error('phy.F90/initialize:','for '//trim(self%name)// ' for mimicking Monod model, a valid Q_fixed (<1.0) is required')
     end if
   end if
   
   !mortality/loss/respiration
   call self%get_parameter(self%zetaN, 'zetaN','molC molN-1', 'C-cost of N uptake', default=0.6_rk)
   call self%get_parameter(self%zetaChl, 'zetaChl','molC gChl-1', 'C-cost of Chlorophyll synthesis', default=0.8_rk)
   call self%get_parameter(self%M0p, 'M0p', 'm3 molN-1 d-1', 'sp. quad. mortality rate',              default=0.1_rk, scale_factor=d_per_s)
   call self%get_parameter(self%Mpart, 'Mpart', '-',   'part of the mortality that goes to detritus',default=0.5_rk)
   

   ! Register state variables
   if ( self%dynQN ) then
     call self%register_state_variable(self%id_phyC,'C','mmolC/m^3','bound-C concentration',0.0_rk,minimum=0.0_rk,vertical_movement=w_phy, specific_light_extinction=0.0_rk)
   else
     call self%register_diagnostic_variable(self%id_d_phyC, 'C','mmolC/m^3', 'bound-C concentration (diagnostic)', &
                                     output=output_instantaneous)
   end if
   call self%register_state_variable(self%id_phyN,'N','mmolN/m^3','bound-N concentration',0.0_rk,minimum=0.0_rk,vertical_movement=w_phy, specific_light_extinction=self%kc)
   
   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phyN)

   ! Register dependencies on external state variables
   call self%register_state_dependency(self%id_din, 'din',   'mmolN/m^3','dissolved inorganic nitrogen')
   call self%register_state_dependency(self%id_don, 'don','mmolN/m^3','dissolved organic nitrogen')
   call self%register_state_dependency(self%id_detN,'detN','mmolN/m^3','detrital nitrogen')
   call self%register_state_dependency(self%id_doc, 'doc','mmolC/m^3','dissolved organic carbon')
   call self%register_state_dependency(self%id_detC,'detC','mmolC/m^3','detrital carbon')

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_Q, 'Q','molN/molC',    'cellular nitrogen Quota',           &
                                     output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Chl, 'Chl','mgChl/m^3',    'Chlorophyll concentration',           &
                                     output=output_instantaneous)                                  
   call self%register_diagnostic_variable(self%id_Chl2C, 'Chl2C','gChl/molC',    'cellular chlorophyll content',           &
                                     output=output_instantaneous)
                                     
   call self%register_diagnostic_variable(self%id_mu, 'mu','/d',    'net sp. growth rate',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fdinphy_sp, 'V_N','molN/molC/d',    'net sp. N uptake',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_respN, 'R_N','/d',    'Respiration cost of N uptake',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_respChl, 'R_Chl','/d',    'Respiration cost of Chl uptake',     &
                                     output=output_time_step_averaged)
                                     
   call self%register_diagnostic_variable(self%id_fV, 'fV','-',    'fV',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fA, 'fA','-',    'fA',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_ThetaHat, 'ThetaHat','-', 'ThetaHat',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fN, 'fN','-',    'N limitation function',&
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fL, 'fL','-',    'Light limitation function',&
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_Tfac, 'Tfac','-',    'Temperature factor',&
                                     output=output_time_step_averaged)
                                     
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmolC/m^3/d','Primary production rate',      &
                                     output=output_time_step_averaged)
   call self%add_to_aggregate_variable(total_PPR,self%id_PPR)
   
   !optional diagnostics
   if ( self%dynQN ) then
     call self%register_diagnostic_variable(self%id_fQ, 'fQ', '-',    'Down-regulation term (only for dynQN)',   &
                                     output=output_instantaneous) 
   end if
   if ( self%mimic_Monod ) then
     call self%register_diagnostic_variable(self%id_fNmonod, 'fN_monod','-',    'Monod function of DIN',   &
                                     output=output_instantaneous) 
   end if
   
   call self%register_diagnostic_variable(self%id_fphydon, 'f_phy_don','molN/m^3/d',    'bulk phy-N loss to detritus',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydoc, 'f_phy_doc','molC/m^3/d',    'bulk phy-C loss to detritus',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydetn, 'f_phy_detn','molN/m^3/d',    'bulk phy-N loss to DOM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydetc, 'f_phy_detc','molN/m^3/d',    'bulk phy-C loss to DOM',           &
                                     output=output_time_step_averaged)
                                     
   ! Register environmental dependencies
   call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_par_dmean, 'PAR_dmean','E/m^2/d','photosynthetically active radiation, daily averaged')
   call self%register_horizontal_dependency(self%id_FDL, 'FDL','-',       'fractional day length')
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   
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
   class (type_nflexpd_phy), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: din,phyC,phyN,parW,par,par_dm,Ld
   real(rk)                   :: ThetaHat,vNhat,muIhat
   real(rk)                   :: Q,Theta,fV,fQ,fA,Rchl,I_zero,ZINT,valSIT
   real(rk)                   :: vN,Vhat_fNT
   real                       :: larg !argument to WAPR(real(4),0,0) in lambert.f90
   real(rk)                   :: tC,Tfac,depth
   real(rk)                   :: mu,respN,mort,Pprod,muIN,KN_monod
   real(rk)                   :: fL,fN,fN_monod
   real(rk)                   :: f_din_phy,f_phy_don,f_phy_detn,f_phy_doc,f_phy_detc
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
   _GET_(self%id_depth,depth)     ! depth
   _GET_(self%id_phyN,phyN)  ! phytoplankton-N
   ! Retrieve current (local) state variable values.
   if ( self%dynQN ) then
     _GET_(self%id_phyC,phyC)  ! phytoplankton-C
     !Dynamically calculated quota is needed for calculating some rates below
     Q= phyN/phyC
     !calculate a down-regulation term for nutrient uptake to avoid Q>Qmax
     fQ=min(1.0,max(0.01, (self%Qmax-Q)/(self%Qmax-self%Q0/2.0)))
     !fQ=1.0
   end if
     
   _GET_(self%id_din,din)    ! nutrients
   
   ! Retrieve current environmental conditions.
   _GET_(self%id_parW,parW)             ! local photosynthetically active radiation
   par=parW* 4.6 * 1e-6   !molE/m2/s
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   _GET_(self%id_par_dmean,par_dm) !in molE/m2/d
   !write(*,*)'P.L275:depth,par_dm',depth,par_dm
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
       ThetaHat=max(0.0,ThetaHat) !  a small positive value 
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
   fL=SIT(self%aI,self%mu0hat,par_dm,ThetaHat,Tfac)
   !write(*,*)'depth,par_dm,SIT',depth,par_dm,1.0-exp(-self%aI*ThetaHat*par_dm/(Tfac*self%mu0hat))
   muIhat = self%mu0hat * Tfac * fL
   
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
   
   if( self%fV_opt ) then
     ! eq. 13  in Smith et al 2016
     fV = (-1.0 + sqrt(1.0 + 1.0 / ZINT) ) * (self%Q0 / 2.0) * muIhat / vNhat
   else
     fV = self%fV_fixed
   end if
   
   if (.not. self%dynQN) then
     !!$ ***  Calculating the optimal cell quota, based on the term ZINT, as calculated above
     if ( self%mimic_Monod ) then
       Q = self%Q_fixed !6.67 !Almost Redfield? (106/16=6.625)
       if (self%KN_monod .le. 0.0_rk) then
          KN_monod = self%V0hat*Tfac/self%A0hat
       else
          KN_monod =self%KN_monod
       end if
       fN_monod = din / ( KN_monod + din)
     else
       ! eq. 14 in Smith et al 2016
       Q = ( 1.0 + sqrt(1.0 + 1.0/ZINT) )*(self%Q0/2.0)       
     end if
     phyC=phyN/Q
   end if
   
   ! Losses due to Chlorophyll
   ! eq. 26 in Smith et al 2016
   Rchl = (muIhat + self%RMchl*Tfac) * ( 1 - fV - self%Q0/(2.0*Q) ) * self%zetaChl * ThetaHat
   
   !  Net specific growth rate, assuming instantantaneous optimal resource allocation. 
      !  Either equation gives the same result, provided fA, fV and QN have been optimized.  
      !  The term with ZINT accounts for the cost of N assimilation, but not for chl maintenance. 
      !  mu = muIhat * ( 1 + 2*( ZINT - sqrt(ZINT*(1 + ZINT)) ) )           - Rchl 
   !!$      mu = muIhat*(1.0 + 2.0*(ZINT - sqrt(ZINT*(1.0+ZINT))) ) - Rchl
   ! eq. 5 in Pahlow and Oschlies 2013 (-Rchl)
   !to prevent model crashing:
   if (din .gt. self%mindin) then !can be interpreted as 'din detection limit' for phytoplankton
     if ( self%mimic_Monod ) then
       fN = fN_monod
     else
       fN = 1.0_rk - fV - self%Q0/(2.0*Q) 
     end if
   else
     fN = 0.0_rk
   end if
   
   muIN =  muIhat * fN
   
   !Total Chl content per C in Cell (eq. 10 in Smith et al 2016)
   Theta= (1 - self%Q0 / 2 / Q - fV) * ThetaHat
   
   if ( self%dynQN ) then !Explicit uptake rate
     !to prevent model crashing:
     if (din .gt. self%mindin) then !can be interpreted as 'din detection limit' for phytoplankton
       vN = fV*vNhat*fQ !molN/molC/d
     else
       vN = 0.0_rk
     end if
   else
       vN = muIN*Q !/d * molN/molC 
   end if
   
   respN=self%zetaN*vN !molC/molN *molN/molC/d = /d
   mu = muIN - (respN+Rchl)
   
   !Calculate fluxes between pools
   if ( self%dynQN ) then 
     f_din_phy = vN * phyC
   else
     f_din_phy = mu*phyN
   end if
   
   ! Mortality
   mort=self%M0p * Tfac * PhyN**2
   f_phy_detn =       self%Mpart  * mort 
   f_phy_don = (1.0 - self%Mpart) * mort
   f_phy_detc = f_phy_detn/Q
   f_phy_doc = f_phy_don/Q
   
   ! Set temporal derivatives
   if ( self%dynQN ) then
     _SET_ODE_(self%id_phyC, mu*phyC - f_phy_doc - f_phy_detc)
   end if
   _SET_ODE_(self%id_phyN, f_din_phy - f_phy_don - f_phy_detn)
   
   ! If externally maintained dim,dom und det pools are coupled:
   _SET_ODE_(self%id_din, -f_din_phy)
   _SET_ODE_(self%id_don,  f_phy_don)
   _SET_ODE_(self%id_detN, f_phy_detn)
   _SET_ODE_(self%id_doc,  f_phy_doc)
   _SET_ODE_(self%id_detC, f_phy_detc)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_Q, Q)
   if ( self%dynQN ) then
     _SET_DIAGNOSTIC_(self%id_fQ, fQ)
   else
      _SET_DIAGNOSTIC_(self%id_d_phyC, phyC)
   end if
   if ( self%mimic_Monod ) then
     _SET_DIAGNOSTIC_(self%id_fNmonod,fN_monod)
   end if
   _SET_DIAGNOSTIC_(self%id_Chl, Theta*phyC)
   _SET_DIAGNOSTIC_(self%id_Chl2C, Theta)
   _SET_DIAGNOSTIC_(self%id_fV, fV)
   _SET_DIAGNOSTIC_(self%id_fA, fA)
   _SET_DIAGNOSTIC_(self%id_fN, fN)
   _SET_DIAGNOSTIC_(self%id_fL, fL)
   _SET_DIAGNOSTIC_(self%id_Tfac, Tfac)
   _SET_DIAGNOSTIC_(self%id_fdinphy_sp, f_din_phy/phyC * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_mu, mu * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_respN, respN * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_respChl, Rchl * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_ThetaHat, ThetaHat) 
   _SET_DIAGNOSTIC_(self%id_PPR, mu*phyC*secs_pr_day) !*s_p_d such that output is in d-1
   
   !Export diagnostic bulk fluxes
   _SET_DIAGNOSTIC_(self%id_fphydoc, f_phy_doc * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_fphydon, f_phy_don * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_fphydetc, f_phy_detc * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_fphydetn, f_phy_detn * secs_pr_day) !*s_p_d such that output is in d-1
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC
!-----------------------------------------------------------------------

   end module nflexpd_phy

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
