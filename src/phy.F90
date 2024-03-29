#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: nflexpd_phy: amorphous phytoplankton class that can be instanced as:
! 1) N-based models (dynQN=false), with:
!!! 1a) (FS) constant-stoichiometry without any flexibility, mimics Monod model
!!! 1b) (IA)(*) 'Instantaneously Acclimating' N:C based on flexible fA,fV (&ThetaHat)
! 2) N&C-based models (dynQN=true), with:
!!! 2a) (DA) 'Dynamically Acclimating' N:C, based on flexible fA,fV (&ThetaHat)
!!! 2b) (DQ) 'Dyamical Quota': variable N:C, based on fixed fV (&fA,ThetaHat)
!
!(*): 'Flex' model in  Smith et al., 2016, J. Plankton Res. 38, 977–992; FABM implementation described by Kerimoglu O., Anugerahanti, P., Smith, S.L., 2020, GMDD (hereafter K20)
! 
! Original Author(s): 
! S. Lan Smith 2014-12-09 (gotm-bio), 
! O. Kerimoglu 2018-11-27 (first fabm version), 2020-11-19 (version in 'K20')
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
      type (type_diagnostic_variable_id)   :: id_muhatNET,id_ZINT,id_Vhat,id_Ahat,id_KN
      type (type_diagnostic_variable_id)   :: id_Q,id_d_phyC,id_Chl,id_Chl2C,id_fV,id_fA,id_ThetaHat
      type (type_diagnostic_variable_id)   :: id_PPR,id_fdinphy_sp,id_mu,id_muNET,id_muG,id_vN,id_respN,id_respChl
      type (type_diagnostic_variable_id)   :: id_muhatG,id_vNhat,id_respHatChl
      type (type_diagnostic_variable_id)   :: id_fQ,id_limfunc_Nmonod,id_fC,id_limfunc_L,id_Tfac
      type(type_diagnostic_variable_id)    :: id_fphydoc,id_fphydon,id_fphydetc,id_fphydetn
      
!     Model parameters
      real(rk) :: kc,w_phy,mindin
      real(rk) :: zetaN,zetaChl,M0p,Mpart,RMChl
      real(rk) :: mu0hat,aI
      real(rk) :: A0hat,V0hat,Q0,Qmax,KN_monod
      real(rk) :: fA_fixed,fV_fixed,TheHat_fixed,Q_fixed,ThetaHat_min
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
   call self%get_parameter(self%ThetaHat_min, 'ThetaHat_min','gChl molC-1', 'Minimum allowed value, which is also used when par<I_0', default=0.1_rk)
   call self%get_parameter(self%RMchl, 'RMchl','d-1', 'loss rate of chlorophyll', default=0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%mu0hat, 'mu0hat','d-1', 'max. potential growth rate', default=5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%aI, 'aI','(m^2 E-1 molC gChl-1)', 'Chl-specific slope of the PI curve', default=1.0_rk) ! really /mol or /micromol?
   !nutrient-related
   call self%get_parameter(self%fA_fixed, 'fA_fixed','-', 'fA to use when fa_opt=false', default=-9.9_rk)
   call self%get_parameter(self%fV_fixed, 'fV_fixed','-', 'fV to use when fv_opt=false', default=-9.9_rk)
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
   call self%register_diagnostic_variable(self%id_Chl2C, 'Chl2C','gChl/gC',    'cellular chlorophyll content',           &
                                     output=output_instantaneous)
                                     
   call self%register_diagnostic_variable(self%id_muNET, 'muNET','/d',    'net cell-specific growth rate before respN',           &
                                     output=output_time_step_averaged)                                  
   call self%register_diagnostic_variable(self%id_muG, 'muG','/d',    'gross cell-specific growth rate ',           &
                                     output=output_time_step_averaged) 
   call self%register_diagnostic_variable(self%id_mu, 'mu','/d',    'net sp. growth rate',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_vN, 'vN','molN/molC/d',    'Specific N uptake rate',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_muhatG, 'muhatG','/d',    'Gross growth rate within chloroplast',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_muhatNET, 'muhatNET','/d',   'Net growth rate within chloroplast',           &
                                     output=output_time_step_averaged)                                  
   call self%register_diagnostic_variable(self%id_vNhat, 'vNhat','molN/molC/d',    'Potential specific N uptake rate',           &
                                     output=output_time_step_averaged)                                     
                                     
   call self%register_diagnostic_variable(self%id_fdinphy_sp, 'f_dinphy','molN/molC/d',    'effective net sp. N uptake',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_respN, 'R_N','/d',    'Respiration cost of N uptake',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_respChl, 'R_Chl','/d',    'Cell-specific respiration cost of Chl uptake',     &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_respHatChl, 'Rhat_Chl','/d',    'Chloroplast-specific respiration cost of Chl uptake',     &
                                     output=output_time_step_averaged)                                     
                                     
   call self%register_diagnostic_variable(self%id_ZINT, 'ZINT','-',    'Qs*(muhatNET/vNhat+zetaN)',           &
                                     output=output_time_step_averaged)                                     
   call self%register_diagnostic_variable(self%id_fV, 'fV','-',    'fV',           &
                                     output=output_time_step_averaged)                                 
                                     
   call self%register_diagnostic_variable(self%id_fA, 'fA','-',    'fA',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_ThetaHat, 'ThetaHat','gChl/gC', 'ThetaHat',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fC, 'fC','-',    'fractional allocation to Carbon fixation (=1-fV-Qs/Q)',&
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_limfunc_L, 'limfunc_L','-',    'Light limitation function',&
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_Tfac, 'Tfac','-',    'Temperature factor',&
                                     output=output_time_step_averaged)
                                     
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmolC/m^3/d','Net primary production rate',      &
                                     output=output_time_step_averaged)
   call self%add_to_aggregate_variable(total_PPR,self%id_PPR)
   
   !optional diagnostics
   if ( self%dynQN ) then
     call self%register_diagnostic_variable(self%id_fQ, 'fQ', '-',    'Down-regulation term (only for dynQN)',   &
                                     output=output_instantaneous) 
   end if
   if ( self%mimic_Monod ) then
     call self%register_diagnostic_variable(self%id_limfunc_Nmonod, 'limfunc_Nmonod','-',    'Monod function of DIN',   &
                                     output=output_instantaneous)
   else
     call self%register_diagnostic_variable(self%id_Vhat, 'Vhat','molN molC-1 d-1',    '(1-fA)*V0hat*fT',   &
                                     output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_Ahat, 'Ahat','m3 mmolC-1 d-1',    'fA*A0hat',   &
                                     output=output_instantaneous)
     call self%register_diagnostic_variable(self%id_KN, 'K_N_equivalent','mmolN/m^3',    '(1-fA)*V0hat*fT/(fA*A0hat)',   &
                                     output=output_instantaneous)                                
   end if
   
   call self%register_diagnostic_variable(self%id_fphydon, 'f_phy_don','mmolN/m^3/d',    'bulk phy-N loss to detritus',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydoc, 'f_phy_doc','mmolC/m^3/d',    'bulk phy-C loss to detritus',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydetn, 'f_phy_detn','mmolN/m^3/d',    'bulk phy-N loss to DOM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydetc, 'f_phy_detc','mmolN/m^3/d',    'bulk phy-C loss to DOM',           &
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
   real(rk)                   :: ThetaHat,vNhat,muhatG,RhatChl,muhatNET
   real(rk)                   :: Q,Theta,fV,fQ,fA,Rchl,I_zero,ZINT,valSIT
   real(rk)                   :: vN,Vhat_fNT,RMchl,zetaChl
   real                       :: larg !argument to WAPR(real(4),0,0) in lambert.f90
   real(rk)                   :: tC,Tfac,depth
   real(rk)                   :: mu,respN,mort,Pprod,muNET,muG,KN_monod
   real(rk)                   :: limfunc_L,fC,limfunc_Nmonod
   real(rk)                   :: f_din_phy,f_phy_don,f_phy_detn,f_phy_doc,f_phy_detc
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
   _GET_(self%id_depth,depth)     ! depth
   !get Ld (fractional day length)
   _GET_HORIZONTAL_(self%id_FDL,Ld)
   _GET_(self%id_temp,tC) ! temperature in Celcius
   _GET_(self%id_din,din)    ! nutrients
   
   _GET_(self%id_phyN,phyN)  ! phytoplankton-N
   
   ! Retrieve current (local) state variable values.
   if ( self%dynQN ) then
     _GET_(self%id_phyC,phyC)  ! phytoplankton-C
     !Dynamically calculated quota is needed for calculating some rates below
     Q= phyN/phyC !Eq.9 in K20
     !relative Quota (used as a down-regulation term in classical droop model)
     fQ=min(1.0,max(0.0, (self%Qmax-Q)/(self%Qmax-self%Q0/2.0))) !not used in K20
   end if
   
   ! Retrieve current environmental conditions.
   _GET_(self%id_parW,parW)             ! local photosynthetically active radiation
   par=parW* 4.6 * 1e-6   !molE/m2/s
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   _GET_(self%id_par_dmean,par_dm) !par is assumed to be daytime average PAR, in molE/m2/d 
   !write(*,*)'P.L275:depth,par_dm',depth,par_dm
   par_dm=(par_dm/secs_pr_day) !in molE/m2/s 
   
   if ( par_dm .lt. 0.0 ) then
     par_dm=par
   end if
   
   !when par_dm=0 (eg., during the first day, where par_dm is not yet calculated, muhatNET and fV becomes both 0, which makes Q=NaN
   !so we set it to a very small value
   if (par_dm .eq. 0.0 ) then
     par_dm=1e-6
   end if
   
   !Calculate intermediate terms:
   !Temperature factor 
   Tfac = FofT(tC)
   
   ! Primary production
   ! Optimization of ThetaHat (optimal Chl content in the chloroplasts)
   I_zero = self%zetaChl * self%RMchl * Tfac / (Ld*self%aI)   ! Eq.27, Threshold irradiance
   zetaChl=self%zetaChl
   RMchl=self%RMchl
   if( self%theta_opt ) then
     if( par_dm .gt. I_zero ) then !in cmo: .and. (mu0>0.0)
       !argument for the Lambert's W function
       larg = (1.0 + self%RMchl * Tfac/(Ld*self%mu0hat*Tfac)) * exp(1.0 + self%aI*par_dm/(self%mu0hat*Tfac*self%zetaChl)) !Eq.26, term in brackets
       ThetaHat = 1.0/self%zetaChl + ( 1.0 -  WAPR(larg, 0, 0) ) * self%mu0hat*Tfac/(self%aI*par_dm) !Eq.26
     else
       !write(*,*)'par_dm,I_0',par_dm,I_zero
       ThetaHat = self%ThetaHat_min  !Eq.26, if I<=IC (=0 in K20)
       !Setting zetaChl=RMchl=0 was necessary for using non-zero ThetaHat_min. 
       !I think we shouldn't do that anymore, since ThatHat_min is perfectly consistent, and 
       !Relatedly, ThataHat_min should be a constant, and not parameter
       !zetaChl=0.0 
       !RMchl=0.0
     end if
   else
     ThetaHat = self%TheHat_fixed !Used for FS in K20
   end if
   
   ! Light limited growth rate
   limfunc_L=SIT(self%aI,self%mu0hat,par_dm,ThetaHat,Tfac) !Eq.22 in K20
   !write(*,*)'depth,par_dm,SIT',depth,par_dm,1.0-exp(-self%aI*ThetaHat*par_dm/(Tfac*self%mu0hat))
   muhatG = Ld * self%mu0hat * Tfac * limfunc_L !Eq.21 in K20
   
   !'Net' light limited growth rate, muhatNET
   !muhatNET=muhatG*(1.0-zetaChl*ThetaHat)-Tfac*RMchl*zetaChl*ThetaHat
   !Chloroplast-specific respiration rate:
   RhatChl=muhatG*zetaChl*ThetaHat + Tfac*RMchl*zetaChl*ThetaHat !Eq.23 in K20
   !Chloroplast-specific net growth rate
   muhatNET=muhatG-RhatChl !Eq. 23 in K20 (= A-cursive in Pahlow etal 2013, Appendix 1)
   if (.not. self%mimic_Monod .and. self%fV_opt .and. par_dm .gt. I_zero .and. muhatNET .lt. 0.0) then
     write(*,'(A,F10.8,A,F10.8,A,F5.2,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8)')'Ld:',Ld,'  fT:',Tfac,'  depth:',depth,'  I_C:',I_zero*86400,'  Idm:',par_dm*86400,'  WAPR:',WAPR(larg, 0, 0),'  ThetaHat:',ThetaHat,'  SI:',limfunc_L,'  muhatNET:',muhatNET*86400
   end if
   
   !Optimal allocation for affinity vs max. uptake
   if( self%fA_opt ) then
     fA = 1.0 / ( 1.0 + sqrt(self%A0hat * din /(Tfac * self%V0hat)) ) !Eq.18 in K20
   else
     fA =  self%fA_fixed !Used for FS in K20
   end if
   
   ! Optimal VNhat
   vNhat = vAff( din, fA, self%A0hat, self%V0hat * Tfac ) !Eq.16,17 in K20
   ! Equivalent way of calculating optimal vNhat
   !vNhat2=vOU(din,self%A0hat,self%V0hat * Tfac)

   if (.not. self%dynQN) then
     if ( self%mimic_Monod ) then
       Q = self%Q_fixed !6.67 !Almost Redfield? (106/16=6.625)
       !Nutrient limitation function in Monod-form
       if (self%KN_monod .le. 0.0_rk) then
          KN_monod = (1.0-self%fA_fixed)*self%V0hat*Tfac/(self%A0hat*self%fA_fixed) !Similar to Eq.19, although KN is explicitly provided there.
       else
          KN_monod =self%KN_monod !Used for FS in K20
       end if
       limfunc_Nmonod = din / ( KN_monod + din)
       !Nutrient limitation function in affinity form:
       !limfunc_Nmonod = fA*self%A0hat*din / ((1.0-fA)*self%V0hat*Tfac+fA*self%A0hat*din)
     else
       !Optimal Q: 
       ZINT = (self%zetaN + muhatNET/vNhat) * self%Q0 / 2.0 !Eq.10, denominator term in sqrt
       !write(*,'(A,4F12.5)')'  (phy) ZINT, muhatG/vNhat:',ZINT,muhatG/vNhat
       Q = ( 1.0 + sqrt(1.0 + 1.0/ZINT) )*(self%Q0/2.0) !Eq.10 in K20 (=Eq10. in PO13)
     end if
     phyC=phyN/Q
   end if
   
   if( self%fV_opt ) then
     fV=(self%Q0/2.0)/Q - self%zetaN*(Q - self%Q0) !Eq.14 in K20 (=Eq.9 in PO13)
     if (fV .lt. 0.0) then
       write(*,*)'nudging negative fV:',fV,' to 0.0'
       fV=0.0
     end if
     !write(*,'(A,F5.2,A,F7.5,A,F7.5,A,F7.5,A,F7.5)')'depth:',depth,'  Q:',Q,'  fV:',fV,'  Q0/2/Q:',(self%Q0/2.0)/Q,'  zN*(Q-Q0):',self%zetaN*(Q - self%Q0)
   else
     fV = self%fV_fixed !Used for FS in K20
   end if
   
   !fC
   fC = 1.0_rk - fV - self%Q0/(2.0*Q) !Eq.11 in K20
   
   !For FS (mimic_Monod):s pecification of a positive fV_fixed implies that a constant Chl:C is desired. For this case, we use fC as calculated by constant fV and Q, and use this to scale RhatChl and ThetaHat
   if ( self%mimic_Monod .and. self%fV_fixed .le. 0.0) then
     !If however no fV is provided,  FS becomes more consistent (with regard to the calcualtion of RChl), although it's not a typical classical model constant C:Chl ratio anymore
     fC = limfunc_Nmonod
   end if 
   
   !can help avoiding model crashing:
   if (din .lt. self%mindin) then ! 'din detection limit' for phytoplankton
     if ( self%mimic_Monod ) then
       limfunc_Nmonod=0.0_rk ! used for FS (as implied by Eq.15)
     else
       fC = 0.0_rk
     end if
   end if
   
   ! Losses due to Chlorophyll
   ! eq. 26 in Smith et al 2016
   !For FS, fV>0 implies constant Chl:C, which necessitates scaling of RhatChl with (1 - fV - self%Q0/(2.0*Q))
   if ( self%mimic_Monod) then
     muG = muhatG * limfunc_Nmonod
   else
     !Default behavior is to scale with fC
     muG = muhatG * fC
     !Rchl = (muhatG + Tfac*RMchl) * fC * zetaChl * ThetaHat
   end if
   
   RChl = RhatChl * fC
   
   !Cellular light-limited growth potential
   muNET = muG - RChl ! Eq.13 in K20; (fC*muhatG=muG, RChl=fC*RhatChl, Eq.25)
   
   !Total Chl content per C in Cell (eq. 10 in Smith et al 2016)
   Theta= fC * ThetaHat
   
   !Calculate respN:
   if ( self%dynQN ) then !Explicit uptake rate
     if ( .not. self%fV_opt) then
       !!For the non-acclimative DA variant: downgregulation based on Q
       respN=self%zetaN*fV*vNhat*fQ !molC/molN *molN/molC/d = /d
     else
       respN=self%zetaN*fV*vNhat !molC/molN *molN/molC/d = /d
     end if
   else
     if ( self%mimic_Monod ) then
       !solve  respN=zetaN*vN=zetaN*muQ from mu=muhatNET*fC-mu*Q*zetaN; mu=muhatNET*fC/(1+Q*zetaN)
       if (self%fV_fixed .gt. 0.0) then
         respN=self%zetaN*muhatNET*limfunc_Nmonod*Q/(1.0+Q*self%zetaN) ! : This can become negative again for muNET>Rchl
       else
         respN=self%zetaN*muhatNET*fC*Q/(1.0+Q*self%zetaN) ! : This can become negative again for muNET>Rchl
       end if
     else
       !like DA: RN based on V = fV*vNhat
       respN=self%zetaN*fV*vNhat ![molC/molN *molN/molC/d = /d]
     end if  
   end if
   
   !Net growth rate
   mu = muNET - respN !Eq.7 Note that muNET already contains -Rchl
   !for FS, this becomes (see Appendix A1 in K20)
   !mu= muNET - zetaN*muNET*Q/(1.0+Q*zetaN)= (muNET (1+QzetaN) - muNET Q ZetaN ) / (1+QzetaN) = muNET/(1+QzetaN)
   !mu+muQzetaN = muNET -> mu=muNET-muQ*zetaN
   
   !if (depth .lt. 1.0) write(*,'(A,6F5.2)')'Ld,mu,muNET,respN,muG,RChl',Ld,mu*secs_pr_day,muNET*secs_pr_day,respN*secs_pr_day,muG*secs_pr_day,RChl*secs_pr_day
   
   if ( self%dynQN ) then !Explicit uptake rate
     !to prevent model crashing:
     if (din .gt. self%mindin) then !can be interpreted as 'din detection limit' for phytoplankton
       if ( .not. self%fV_opt) then
         !For the non-acclimative DA variant: downgregulation based on Q
         vN = fV*vNhat*fQ !molN/molC/d
       else
         vN = fV*vNhat    !Eq.12 in K20 [molN/molC/d] 
       end if  
     else
       vN = 0.0_rk
     end if
   else
       !Balanced growth:
       vN=mu*Q            !Eq.6 in K20 [molN/molC/d]
   end if
   
   !Calculate fluxes between pools
   f_din_phy = vN * phyC  !Eq.5 in K20
   
   ! Mortality
   mort=self%M0p * Tfac * PhyN**2
   f_phy_detn =       self%Mpart  * mort  !Table 1
   f_phy_detc = f_phy_detn/Q              !Table 1
   f_phy_don = (1.0 - self%Mpart) * mort  !Doesn't appear in K20, since Mpart=1 -> f_phy_don=0 
   f_phy_doc = f_phy_don/Q                !Doesn't appear in K20, since Mpart=1 -> f_phy_don=0 
   
   ! Set temporal derivatives
   _SET_ODE_(self%id_phyN, f_din_phy - f_phy_don - f_phy_detn) !Eq.1a (f_phy_don doesn't appear in K20, since it's=0 (see above))
   if ( self%dynQN ) then
     _SET_ODE_(self%id_phyC, mu*phyC - f_phy_doc - f_phy_detc) !Eq.1b, Eq8 (for mu*phyC) (f_phy_doc doesn't appear in K20, since it's=0 (see above))
   end if
   
   ! If externally maintained dim,dom und det pools are coupled:
   _SET_ODE_(self%id_detN, f_phy_detn)  !Eq.2a in K20
   _SET_ODE_(self%id_detC, f_phy_detc)  !Eq.2b in K20
   _SET_ODE_(self%id_don,  f_phy_don)   !Doesn't appear in K20 since f_phy_don = 0 (see above)
   _SET_ODE_(self%id_doc,  f_phy_doc)   !Doesn't appear in K20 since f_phy_doc = 0 (see above)
   _SET_ODE_(self%id_din, -f_din_phy)   !eq.4 in K20
   
   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_Q, Q)
   if ( self%dynQN ) then
     _SET_DIAGNOSTIC_(self%id_fQ, fQ)
   else
      _SET_DIAGNOSTIC_(self%id_d_phyC, phyC)
   end if
   if ( self%mimic_Monod ) then
     _SET_DIAGNOSTIC_(self%id_limfunc_Nmonod,limfunc_Nmonod)
   else
     _SET_DIAGNOSTIC_(self%id_Vhat,(1.0-fA)*self%V0hat*Tfac*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_Ahat,fA*self%A0hat*secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_KN,(1.0-fA)*self%V0hat*Tfac/(fA*self%A0hat)) 
   end if
   _SET_DIAGNOSTIC_(self%id_Chl, Theta*phyC) 
   _SET_DIAGNOSTIC_(self%id_Chl2C, Theta/12.0) !gChl/molC*1molC/12gC =gChl/gC
   if (.not. self%dynQN .and. .not. self%mimic_Monod) then
     _SET_DIAGNOSTIC_(self%id_ZINT, ZINT)
   end if
   _SET_DIAGNOSTIC_(self%id_fV, fV)
   _SET_DIAGNOSTIC_(self%id_fA, fA)
   _SET_DIAGNOSTIC_(self%id_fC, fC)
   _SET_DIAGNOSTIC_(self%id_limfunc_L, limfunc_L)
   _SET_DIAGNOSTIC_(self%id_Tfac, Tfac)
   _SET_DIAGNOSTIC_(self%id_fdinphy_sp, f_din_phy/phyC * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_mu, mu * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_muNET, muNET * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_muG, muG * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_vN, vN * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_muhatG, muhatG * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_muhatNET, muhatNET * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_vNhat, vNhat * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_respN, respN * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_respChl, Rchl * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_respHatChl, RhatChl * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_ThetaHat, ThetaHat/12.0) !gChl/molC*1molC/12gC =gChl/gC
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
