#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NflexPD - phytoplankton component, Carbon-based version
! Original Authors: S. Lan Smith, 2014-12-09
! FABM implementation: O. Kerimoglu 20181122
!
! !INTERFACE:
   module NflexPD_phyabio_Cbased
!
! !DESCRIPTION:
!
! !USES:
   use lambert
   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_NflexPD_phyabio_Cbased
!     Variable identifiers
      !phy
      type (type_state_variable_id)        :: id_phyC,id_phyN
      type (type_dependency_id)            :: id_parW,id_temp,id_par_dmean,id_depth,id_depFDL
      !type (type_horizontal_dependency_id) :: id_depFDL
      type (type_diagnostic_variable_id)   :: id_muhatNET,id_ZINT,id_Vhat,id_Ahat,id_KN
      type (type_diagnostic_variable_id)   :: id_Q,id_d_phyN,id_Chl,id_Chl2C,id_fV,id_fA,id_ThetaHat
      type (type_diagnostic_variable_id)   :: id_PPR,id_fdinphy_sp,id_mu,id_muNET,id_muhatG,id_vNhat,id_vN,id_respN,id_respChl
      type (type_diagnostic_variable_id)   :: id_fQ,id_limfunc_Nmonod,id_fC,id_limfunc_L,id_Tfac
      type(type_diagnostic_variable_id)    :: id_fphydoc,id_fphydon,id_fphydetc,id_fphydetn
      type (type_dependency_id)            :: id_dep_delta_t,id_dep_delta_din,id_dep_delta_par
      
      !abio
      type (type_state_variable_id)     :: id_din,id_don,id_doc,id_detn,id_detc
      !type (type_dependency_id)         :: id_depth!,id_temp,id_parW,
      type (type_dependency_id)         :: id_parW_dmean
      type (type_horizontal_dependency_id)  :: id_lat
      type (type_global_dependency_id)  :: id_doy
      type (type_diagnostic_variable_id):: id_dPAR,id_dPAR_dmean,id_dFDL
      !type (type_surface_diagnostic_variable_id):: id_dFDL
      type (type_diagnostic_variable_id):: id_fdetdon,id_fdetdoc,id_fdondin,id_fdocdic
      type (type_bottom_diagnostic_variable_id):: id_detn_sed,id_detc_sed
      
      !for saving and accessing the doy,din and par of the previous integration time step
      !for doy and delta_t, global diagnostic and dependency would be better but they don't exist
      type (type_diagnostic_variable_id)   :: id_ddoy,id_delta_t
      type (type_dependency_id)            :: id_ddoy_dep
      !for delta_par and delta_din we need 3-D diagnostics anyway
      type (type_diagnostic_variable_id)   :: id_ddin,id_delta_din,id_delta_par
      type (type_dependency_id)            :: id_ddin_dep,id_dpardm_dep
      
!     Model parameters
      !phy
      real(rk) :: kc,w_phy,mindin
      real(rk) :: zetaN,zetaChl,M0p,Mpart,RMChl
      real(rk) :: mu0hat,aI
      real(rk) :: A0hat,V0hat,Q0,Qmax,KN_monod
      real(rk) :: fA_fixed,fV_fixed,TheHat_fixed,Q_fixed,ThetaHat_min
      logical  :: dynQN,fV_opt,fA_opt,Theta_opt,mimic_Monod
      real(rk) :: dic_per_n
      
      !abio
      real(rk) :: w_det,kdet,kdon,par0_dt0,kc_dt0
      

      contains

      procedure :: initialize
      !procedure :: do_surface
      procedure :: do_bottom  
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
!  Here, the NflexPD_phyabio_Cbased namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_NflexPD_phyabio_Cbased), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk), parameter :: N2C_RF = 16._rk/106._rk !Redfield N:C ratio
   real(rk)            :: w_phy
   real(rk)            :: w_det, kc
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
   
   !abio
   call self%get_parameter(self%w_det,     'w_det','m d-1',    'vertical velocity (<0 for sinking)',default=-5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(kc,      'kc', 'm2 mmol-1','specific light extinction',         default=0.03_rk)
   call self%get_parameter(self%kdet,'kdet','d-1',      'sp. rate for f_det_don',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%kdon,'kdon','d-1',      'sp. rate for f_don_din',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%par0_dt0,'par0_dt0','W m-2', 'daily average par at the surface on the first time step',  default=4.5_rk)
   call self%get_parameter(self%kc_dt0,'kc_dt0','m-1', 'attenuaton coefficient on the first time step',  default=0.2_rk)

   ! Register state and diagnostic variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
   !ABIO:
   call self%register_state_variable(self%id_din,'din','mmolN/m^3','DIN concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true., &
                                specific_light_extinction=0.0_rk)
   call self%register_state_variable(self%id_don,'don','mmolN/m^3','DON concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true., &
                                specific_light_extinction=0.0_rk)
   call self%register_state_variable(self%id_doc,'doc','mmolC/m^3','DOC concentration',     &
                                6.625_rk,minimum=0.0_rk,no_river_dilution=.true., &
                                specific_light_extinction=0.0_rk)
   call self%register_state_variable(self%id_detn,'detn','mmolN/m^3','Det-N concentration',    &
                                1.0_rk,minimum=0.0_rk,vertical_movement=self%w_det, &
                                specific_light_extinction=kc)
   call self%register_state_variable(self%id_detc,'detc','mmolC/m^3','Det-C concentration',    &
                                6.625_rk,minimum=0.0_rk,vertical_movement=self%w_det, &
                                specific_light_extinction=0.0_rk)
                                
   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_din)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_don)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_detn)
   
   ! Register diagnostic variables
   !ABIO
   call self%register_bottom_diagnostic_variable(self%id_detn_sed,'detn_sed','mmolN/m^2/d','sedimentation rate of detN')
   call self%register_bottom_diagnostic_variable(self%id_detc_sed,'detc_sed','mmolC/m^2/d','sedimentation rate of detC')
   
   !call self%register_surface_diagnostic_variable(self%id_dFDL,'FDL','-',       'fractional day length')
   call self%register_diagnostic_variable(self%id_dFDL,'FDL','-',       'fractional day length')
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','E/m^2/d',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_dPAR_dmean, 'PAR_dmean','E/m^2/s','photosynthetically active radiation, daily averaged')
   
   call self%register_diagnostic_variable(self%id_fdetdon, 'f_det_don','mmolN/m^3/d',    'bulk N flux from detritus to DOM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fdetdoc, 'f_det_doc','mmolC/m^3/d',    'bulk C flux from detritus to DOM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fdondin, 'f_don_din','mmolN/m^3/d',    'bulk N flux from DOM to DIM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fdocdic, 'f_doc_dic','mmolC/m^3/d',    'bulk C flux from DOM to DIM',           &
                                     output=output_time_step_averaged) 
                                     
   ! Register environmental dependencies
   !moved to the end
   
   !for saving and accessing the doy,din and par of the previous integration time step
   call self%register_diagnostic_variable(self%id_ddoy,'ddoy','d', 'diagn_number_of_days_since_start_of_the_year',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_ddoy_dep,'ddoy','d','diagn_number_of_days_since_start_of_the_year')
   call self%register_diagnostic_variable(self%id_delta_t,'delta_t','s','diff betw current and prev time step',&
                     output=output_instantaneous)
   
   
   call self%register_diagnostic_variable(self%id_ddin,'ddin','mmolN/m^3', 'diagn. din conc',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_ddin_dep,'ddin','mmolN/m^3', 'diagn din conc')
   call self%register_diagnostic_variable(self%id_delta_din,'delta_din','mmolN/m^3','diff betw current and prev time step',&
                     output=output_instantaneous)
   
   !dependencies moved to the end
   !call self%register_dependency(self%id_dpardm_dep,'PAR_dmean','E/m^2/s',       'photosynthetically active radiation, daily averaged')
   !call self%register_diagnostic_variable(self%id_delta_par,'delta_par','E/m^2/s','diff betw current and prev time step',&
   !                  output=output_instantaneous)
                     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
   !PHY
   call self%register_state_variable(self%id_phyC,'C','mmolC/m^3','bound-C concentration',0.0_rk,minimum=0.0_rk,vertical_movement=w_phy, specific_light_extinction=0.0_rk)
   if ( self%dynQN ) then
     call self%register_state_variable(self%id_phyN,'N','mmolN/m^3','bound-N concentration',0.0_rk,minimum=0.0_rk,vertical_movement=w_phy, specific_light_extinction=0.0_rk)
     ! Register contribution of diagnostic to global aggregate variables.
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phyN)
   else
     call self%register_diagnostic_variable(self%id_d_phyN, 'N','mmolN/m^3',    'bound-N concentration (diag)',           &
                                     output=output_instantaneous)
     ! Register contribution of diagnostic to global aggregate variables.
     call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_d_phyN)                                     
   end if
   
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_Q, 'Q','molN/molC',    'cellular nitrogen Quota',           &
                                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Chl, 'Chl','mgChl/m^3',    'Chlorophyll concentration',           &
                                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Chl2C, 'Chl2C','gChl/gC',    'cellular chlorophyll content',           &
                                     output=output_instantaneous)
                                     
   call self%register_diagnostic_variable(self%id_muNET, 'muNET','/d',    'muNET - Rchl',           &
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
   call self%register_diagnostic_variable(self%id_respChl, 'R_Chl','/d',    'Respiration cost of Chl uptake',     &
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
                                     
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmolC/m^3/d','Primary production rate',      &
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
   
   call self%register_diagnostic_variable(self%id_fphydon, 'f_phy_don','molN/m^3/d',    'bulk phy-N loss to detritus',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydoc, 'f_phy_doc','molC/m^3/d',    'bulk phy-C loss to detritus',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydetn, 'f_phy_detn','molN/m^3/d',    'bulk phy-N loss to DOM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fphydetc, 'f_phy_detc','molN/m^3/d',    'bulk phy-C loss to DOM',           &
                                     output=output_time_step_averaged)
                                     
   ! Register environmental dependencies
   !ABIO:
   
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_parW_dmean,temporal_mean(self%id_parW,period=1._rk*86400._rk,resolution=1._rk))
   
   call self%register_dependency(self%id_dpardm_dep,'PAR_dmean','E/m^2/d',       'photosynthetically active radiation, daily averaged')
   call self%register_diagnostic_variable(self%id_delta_par,'delta_par','E/m^2/d','diff betw current and prev time step',&
                     output=output_instantaneous)
   
   !PHY
   !call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_par_dmean, 'PAR_dmean','E/m^2/d','photosynthetically active radiation, daily averaged')
   !call self%register_horizontal_dependency(self%id_depFDL, 'FDL','-',       'fractional day length dependency')
   call self%register_dependency(self%id_depFDL, 'FDL','-',       'fractional day length dependency')
   !call self%register_dependency(self%id_temp,standard_variables%temperature)
   !call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   
   call self%register_dependency(self%id_dep_delta_t, 'delta_t','s','diff betw current and prev time step')
   call self%register_dependency(self%id_dep_delta_din, 'delta_din','mmolN/m^3','diff in DIN betw current and prev time step')
   call self%register_dependency(self%id_dep_delta_par, 'delta_par','E/m^2/d','diff in PAR betw current and prev time step')
   
   !call self%register_dependency(self%id_dep_ThetaHat, 'ThetaHat','-', 'ThetaHat')
   
   end subroutine initialize
!EOC


! !-----------------------------------------------------------------------
! !BOP
! !
! ! !IROUTINE: Right hand sides of Abiotic model
! !
! ! !INTERFACE:
!    subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
! !
! ! !INPUT PARAMETERS:
!    class (type_NflexPD_phyabio_Cbased), intent(in)     :: self
!    _DECLARE_ARGUMENTS_DO_SURFACE_
! !
! ! !LOCAL VARIABLES:
!    real(rk)                   :: lat,doy,Ld
! !EOP
! !-----------------------------------------------------------------------
! !BOC
!    ! Enter spatial loops (if any)
!    _HORIZONTAL_LOOP_BEGIN_
!    
!    _GET_HORIZONTAL_(self%id_lat,lat)
!    _GET_GLOBAL_(self%id_doy,doy)
!    Ld=FDL(lat,doy)
!    _SET_SURFACE_DIAGNOSTIC_(self%id_dFDL,Ld) !Fractional day length
!    ! Leave spatial loops (if any)
!    _HORIZONTAL_LOOP_END_
! 
!    end subroutine do_surface
! !EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Abiotic model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_NflexPD_phyabio_Cbased), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: spsedrate
   real(rk)                   :: detn,detc
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_
   
   _GET_(self%id_detn,detn)
   _GET_(self%id_detc,detc)
   spsedrate=-1.0_rk*self%w_det
   !write(*,*)'detn,det_sed',detn,detn*spsedrate*secs_pr_day
   _SET_BOTTOM_DIAGNOSTIC_(self%id_detn_sed,detn*spsedrate*secs_pr_day)
   _SET_BOTTOM_DIAGNOSTIC_(self%id_detc_sed,detc*spsedrate*secs_pr_day)
   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
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
   class (type_NflexPD_phyabio_Cbased), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   !phy
   real(rk)                   :: din,phyC,phyN,parW,parE,parE_dm,Ld
   real(rk)                   :: ThetaHat,vNhat,muhatG,RhatChl,muhatNET
   real(rk)                   :: Q,Theta,fV,fQ,fA,Rchl,I_zero,ZINT,valSIT
   real(rk)                   :: vN,Vhat_fNT,RMchl,zetaChl
   real                       :: larg !argument to WAPR(real(4),0,0) in lambert.f90
   real(rk)                   :: tC,Tfac,depth
   real(rk)                   :: mu,respN,mort,Pprod,muNET,KN_monod
   real(rk)                   :: limfunc_L,fC,limfunc_Nmonod
   real(rk)                   :: f_din_phy,f_phy_don,f_phy_detn,f_phy_doc,f_phy_detc
   real(rk)                   :: delQ_delt,delQ_delI,delQ_delN,dI_dt,dN_dt
   real(rk)                   :: delQ_delZ,delZ_delI,delZ_delN
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
   
   !abio
   real(rk)                   :: detn, detc, don, doc
   real(rk)                   :: f_det_don, f_det_doc, f_don_din, f_doc_dic
   !real(rk)                   :: parE,parE_dm,depth,tC,parW !already declared for phy
   real(rk)                   :: lat,doy,doy_prev,delta_t
   real(rk)                   :: din_prev,delta_din !din
   real(rk)                   :: parW_dm,parEdm_prev,delta_parE
   !real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
   !ABIO:
   _GET_(self%id_depth,depth)     ! depth
   ! Retrieve current (local) state variable values.
   _GET_(self%id_detc,detc) ! detrital carbon
   _GET_(self%id_detn,detn) ! detrital nitrogen
   _GET_(self%id_doc,doc) ! dissolved organic carbon
   _GET_(self%id_don,don) ! dissolved organic nitrogen
   ! Retrieve environmental dependencies
   _GET_(self%id_temp,tC) ! temperature in Celcius
   
   !Calculate Fractional day length
   _GET_HORIZONTAL_(self%id_lat,lat)
   _GET_GLOBAL_(self%id_doy,doy)
   Ld=FDL(lat,doy)
   
   _GET_(self%id_parW,parW) ! local photosynthetically active radiation (PAR)
   _GET_(self%id_parW_dmean,parW_dm) !current daily average PAR
   if ( parW_dm .lt. 0.0 ) then
     parW_dm=parW
   end if
   parE = parW * 4.6 * 1e-6* secs_pr_day ![mol/m2/d]
   parE_dm= parW_dm * 4.6 * 1e-6* secs_pr_day/Ld  ![mol/m2/d]
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   
   !For providing the delta_t,delta_din and delta_par between the current and previous time step
   _GET_(self%id_din,din) ! din
   _GET_(self%id_ddoy_dep,doy_prev)  ! day of year at the previous time step
   !write(*,*)' (abio.1) doy_prev(s),doy(s),Ld',doy_prev*secs_pr_day,doy*secs_pr_day,Ld
   !Access the par and din at the previous time step and set the diagnostic only if the time step has really advanced
   if (doy .gt. doy_prev) then
     
     !Access the values at the prev. time step as recorded by the diagnostic variables
     _GET_(self%id_ddin_dep,din_prev)
     _GET_(self%id_dPARdm_dep,parEdm_prev) !mol/m2/d
     
     !in the first time step, strange things may happen, as the diagnostics are not available yet
     if (doy_prev .lt. 0.0) then
       doy_prev = -1.0 ! just an arbitrary finite number, as the delta_din&par will be 0      
       din_prev=din ! such that delta_din=0
       parEdm_prev=parE_dm ! such that delta_par=0
     end if
     
     !calculate the deltas
     delta_t=(doy-doy_prev)*secs_pr_day !days to secs
     delta_din=din-din_prev      
     delta_parE= parE_dm-parEdm_prev !mol/m2/d
     
     !assume no change in par within the first day
     if (doy .lt. 1.0+3*delta_t/86400) then
       delta_parE=0.0
     end if
     
     !write(*,*)' (abio.2) pardm_prev,pardm,delta_par',parEdm_prev,parE_dm,delta_parE,'  din_prev,din',din_prev,din
     
     !set the diagnostics
     
     ! Export diagnostic variables
     _SET_DIAGNOSTIC_(self%id_dFDL,Ld) !Fractional day length
     _SET_DIAGNOSTIC_(self%id_ddoy,doy)
     _SET_DIAGNOSTIC_(self%id_ddin, din)
     _SET_DIAGNOSTIC_(self%id_dPAR, parE) ! mol/m2/d
     _SET_DIAGNOSTIC_(self%id_dPAR_dmean, parE_dm) !mol/m2/d
     
     !Export diagnostic bulk fluxes
     _SET_DIAGNOSTIC_(self%id_fdetdon, f_det_don * secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_fdetdoc, f_det_doc * secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_fdondin, f_don_din * secs_pr_day)
     _SET_DIAGNOSTIC_(self%id_fdocdic, f_doc_dic * secs_pr_day)
     
     !Diagnostics required to calculate balance fluxes
     _SET_DIAGNOSTIC_(self%id_delta_t,delta_t)
     _SET_DIAGNOSTIC_(self%id_delta_din,delta_din)
     _SET_DIAGNOSTIC_(self%id_delta_par,delta_parE) !mol/m2/d
     
   else
     _GET_(self%id_par_dmean,parE_dm) !in molE/m2/d
     _GET_(self%id_dep_delta_t,delta_t)
     _GET_(self%id_dep_delta_din,delta_din)
     _GET_(self%id_dep_delta_par,delta_parE) !mol/m2/d
   end if
   
   !Calculate intermediate terms:
   !Temperature factor 
   Tfac = FofT(tC)
   
   !Calculate fluxes between pools
   f_det_don = self%kdet * Tfac * detn !Table 1 in K20
   f_det_doc = self%kdet * Tfac * detc !Table 1 in K20
   f_don_din = self%kdon * Tfac * don  !Table 1 in K20
   f_doc_dic = self%kdon * Tfac * doc  !Table 1 in K20
   
   ! Set temporal derivatives
   ! moved below to combine with the terms resulting from phyto growth
   !_SET_ODE_(self%id_detn, -f_det_don)             !Sink term in Eq.2a
   !_SET_ODE_(self%id_detc, -f_det_doc)             !Sink term in Eq.2b
   !_SET_ODE_(self%id_don,   f_det_don - f_don_din) !Eq.3a
   !_SET_ODE_(self%id_doc,   f_det_doc - f_doc_dic) !Eq.3b
   !_SET_ODE_(self%id_din,   f_don_din)             !Source term in Eq.4
   
   !PHY
   ! Retrieve current (local) state variable values.
   !_GET_(self%id_depth,depth)     ! depth
   !get Ld (fractional day length)
   !_GET_HORIZONTAL_(self%id_FDL,Ld)
   !_GET_(self%id_temp,tC) ! temperature in Celcius
   !_GET_(self%id_din,din)    ! nutrients
   
   _GET_(self%id_phyC,phyC)  ! phytoplankton-C
   if ( self%dynQN ) then
     _GET_(self%id_phyN,phyN)  ! phytoplankton-N
     !Dynamically calculated quota is needed for calculating some rates below
     Q= phyN/phyC !Eq.9 in K20
     !relative Quota (used as a down-regulation term in classical droop model)
     fQ=min(1.0,max(0.0, (self%Qmax-Q)/(self%Qmax-self%Q0/2.0))) !not used in K20
   end if
   
   ! Retrieve current environmental conditions.
   !_GET_(self%id_parW,parW)             ! local photosynthetically active radiation
   !parE=parW* 4.6 * 1e-6   !molE/m2/s
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   !_GET_(self%id_par_dmean,parE_dm) !par is assumed to be daytime average PAR, in molE/m2/d 
   !write(*,*)'P.L275:depth,parE_dm',depth,par_dm
   parE_dm=(parE_dm/secs_pr_day) !in molE/m2/s 
   
   if ( parE_dm .lt. 0.0 ) then
     parE_dm=parE
   end if
   
   !when parE_dm=0 (eg., during the first day, where parE_dm is not yet calculated, muhatNET and fV becomes both 0, which makes Q=NaN
   !so we set it to a very small value
   if (parE_dm .eq. 0.0 ) then
     parE_dm=1e-6
   end if
   
   ! delta_t,delta_parE,delta_din
   !_GET_GLOBAL_(self%id_doy,doy)  ! day of year
   !_GET_HORIZONTAL_(self%id_ddoy_dep,doy_prev)  ! day of year at the previous time step
   !_GET_(self%id_dep_delta_t,delta_t)
   !_GET_(self%id_dep_delta_din,delta_din)
   !_GET_(self%id_dep_delta_parE,delta_parE) !mol/m2/d
   delta_parE=(delta_parE/secs_pr_day) !mol/m2/s
   !write(*,'(A,F10.1,5F12.5)')'  (phy.1) doy(s),delta_I,delta_N,parE,din,phyC:',doy*secs_pr_day,delta_parE,delta_din,parE_dm,din,phyC
   
   !Calculate intermediate terms:
   !Temperature factor 
   Tfac = FofT(tC)
   
   ! Primary production
   ! Optimization of ThetaHat (optimal Chl content in the chloroplasts)
   I_zero = self%zetaChl * self%RMchl * Tfac / (Ld*self%aI)   ! Threshold irradiance
   zetaChl=self%zetaChl
   RMchl=self%RMchl
   if( self%theta_opt ) then
      if( parE_dm .gt. I_zero ) then
       !argument for the Lambert's W function
       larg = (1.0 + self%RMchl * Tfac/(Ld*self%mu0hat*Tfac)) * exp(1.0 + self%aI*parE_dm/(self%mu0hat*Tfac*self%zetaChl)) !Eq.26, term in brackets
       ThetaHat = 1.0/self%zetaChl + ( 1.0 -  WAPR(larg, 0, 0) ) * self%mu0hat*Tfac/(self%aI*parE_dm) !Eq.26
    else
       !write(*,*)'parE_dm,I_0',parE_dm,I_zero
       ThetaHat = self%ThetaHat_min  !Eq.26, if I<=IC (=0 in K20)
       !Setting zetaChl=RMchl=0 was necessary for using non-zero ThetaHat_min. 
       !I think we shouldn't do that anymore, since ThatHat_min is perfectly consistent, and 
       !Relatedly, ThataHat_min should be a constant, and not parameter
       !zetaChl=0.0 
       !RMchl=0.0 
    end if
 else
     ThetaHat = self%TheHat_fixed
  end if


   ! Light limited growth rate
   limfunc_L=SIT(self%aI,self%mu0hat,parE_dm,ThetaHat,Tfac) !Eq.22 in K20
   !write(*,*)'depth,parE_dm,SIT',depth,parE_dm,1.0-exp(-self%aI*ThetaHat*parE_dm/(Tfac*self%mu0hat))
   muhatG = Ld * self%mu0hat * Tfac * limfunc_L !Eq.21 in K20
   
   !'Net' light limited growth rate, muhatNET
   !muhatNET=muhatG*(1.0-zetaChl*ThetaHat)-Tfac*RMchl*zetaChl*ThetaHat
   !Chloroplast-specific respiration rate:
   RhatChl=muhatG*zetaChl*ThetaHat + Tfac*RMchl*zetaChl*ThetaHat !Eq.23 in K20
   !Chloroplast-specific net growth rate
   muhatNET=muhatG-RhatChl !Eq. 23 in K20 (= A-cursive in Pahlow etal 2013, Appendix 1)
   if (.not. self%mimic_Monod .and. self%fV_opt .and. parE_dm .gt. I_zero .and. muhatNET .lt. 0.0) then
     write(*,'(A,F10.8,A,F10.8,A,F5.2,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8)')'Ld:',Ld,'  fT:',Tfac,'  depth:',depth,'  I_C:',I_zero*86400,'  Idm:',parE_dm*86400,'  WAPR:',WAPR(larg, 0, 0),'  ThetaHat:',ThetaHat,'  SI:',limfunc_L,'  muhatNET:',muhatNET*86400
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
     phyN=phyC*Q
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
   !can help avoiding model crashing:
   if (din .gt. self%mindin) then ! 'din detection limit' for phytoplankton
     if ( self%mimic_Monod ) then
       fC = limfunc_Nmonod ! used for FS (as implied by Eq.15)
     else
       fC = 1.0_rk - fV - self%Q0/(2.0*Q) !Eq.11 in K20
     end if
   else
     fC = 0.0_rk
   end if
   
   ! Losses due to Chlorophyll
   ! eq. 26 in Smith et al 2016
   !For FS, fV>0 implies constant Chl:C, which necessitates scaling of RhatChl with (1 - fV - self%Q0/(2.0*Q))
   if ( self%mimic_Monod .and. self%fV_fixed .gt. 0.0 ) then
     RChl = RhatChl * ( 1 - fV - self%Q0/(2.0*Q) )
     !Rchl = (muhatG + Tfac*RMchl) * ( 1 - fV - self%Q0/(2.0*Q) ) * zetaChl * ThetaHat
   else
     !Default behavior is to scale with fC
     RChl = RhatChl * fC
     !Rchl = (muhatG + Tfac*RMchl) * fC * zetaChl * ThetaHat
   end if
   
   !write(*,*)'depth,DIN,fC,fV,Q,Q0/(2.0*Q):',depth,din,fC,fV,Q,self%Q0/(2.0*Q)   
   !muNET =  muhatNET * fC
   muNET =fC*muhatG - RChl ! Eq.13 in K20; (fC*muhatG=muG, RChl=fC*RhatChl, Eq.25)
   
   !Total Chl content per C in Cell (eq. 10 in Smith et al 2016)
   if ( self%mimic_Monod .and. self%fV_fixed .gt. 0.0 ) then
     !Specification of a positive fV_fixed implies that a constant Chl:C is desired. This can considered to be inconsistent with regard to the calculation of RChl (with fC)
     Theta= (1 - self%Q0 / 2 / Q - fV) * ThetaHat !Used for FS in K20
   else
     !With this, FS becomes more consistent (with regard to the calcualtion of RChl), although it's not a typical classical model constant C:Chl ratio anymore
     Theta= fC * ThetaHat !Eq.24 in K20
   end if
   
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
       !Attempt0: just like the others: this gives entirely unrealistic results (too high at the bottom layers):
       !respN=self%zetaN*fV*vNhat !molC/molN *molN/molC/d = /d
       !Attempt1: vN=mu*Q; #this is wrong, as it becomes negative for Rtot>mu
       !Attempt2: solve  respN=zetaN*vN=zetaN*muQ from mu=muhatNET*fC-mu*Q*zetaN; mu=muhatNET*fC/(1+Q*zetaN);
       respN=self%zetaN*muhatNET*fC*Q/(1.0+Q*self%zetaN) ! : This can become negative again for muNET>Rchl
       !Attempt3: take up proportional to light limited gross growth rate
       !respN=self%zetaN*muhatG*fC*Q  !where, muhatG * fC *Q=vN
     else
       !like DA: RN based on V = fV*vNhat
       respN=self%zetaN*fV*vNhat ![molC/molN *molN/molC/d = /d]
       !like for the FS: RN calculated based on V = mu*Q.
       !respN=self%zetaN*muhatNET*fC*Q/(1.0+Q*self%zetaN)
       !Problem: as muhatNET=0 in deep layers,RN becomes 0, which makes vN>0, which in turn makes mu>0
     end if  
   end if
   
   !Net growth rate
   mu = muNET - respN !Eq.7 Note that muNET already contains -Rchl
   !for FS, this becomes (see Appendix A1 in K20)
   !mu= muNET - zetaN*muNET*Q/(1.0+Q*zetaN)= (muNET (1+QzetaN) - muNET Q ZetaN ) / (1+QzetaN) = muNET/(1+QzetaN)
   !mu+muQzetaN = muNET -> mu=muNET-muQ*zetaN
   
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
     if ( self%mimic_Monod ) then
       vN=mu*Q            !Eq.6 in K20 [molN/molC/d]
     else
       ! Calculate the balance-flux
       Vhat_fNT= self%V0hat*din/(self%V0hat/self%A0hat + 2.0 * sqrt((self%V0hat*din/self%A0hat)) + din) !eq.20 in S16
       !write(*,'(A,5F12.5)')'  (phy) Vhat_fNT, self%V0hat*din, self%V0hat/self%A0hat, 2.0*sqrt((self%V0hat*din/self%A0hat)), din', Vhat_fNT, self%V0hat*din, self%V0hat/self%A0hat, 2.0*sqrt((self%V0hat*din/self%A0hat)), din
       !N-uptake:
       !through changes in Q (re-location of N)
       !ZINT=0.01170 obtained in the first time step is too small
       delQ_delZ= -self%Q0/(4*ZINT*sqrt(ZINT*(1.+ZINT))) !delQ/delZ, eq. A-2&3 (Z=ZINT)
       !write(*,'(A,4F12.5)')'  (phy) delQ_delZ,self%Q0,ZINT,4*ZINT*sqrt(ZINT*(1.+ZINT)):',delQ_delZ,self%Q0,ZINT,4*ZINT*sqrt(ZINT*(1.+ZINT))
       !(exp(-ai*ThetaHat*parE_dm/(mu0hat*Tfac))=1-valSIT

       !!$   delZ_delI= self%Q0*self%aI*ThetaHat/(2*din*Vhat_fNT)*(1-valSIT) !delZ/delI, eq.A-4 in S16
       !!$ SLS (20181205): factor of 'din' in denominator not found in FlexPFT code, omitting here
       delZ_delI= self%Q0*self%aI*ThetaHat/(2*Vhat_fNT)*(1-valSIT) !delZ/delI, eq.A-4 in S16

       !write(*,'(A,4F15.5)')'  (phy.2) delZ_delI, ThetaHat, Vhat_fNT, (1-valSIT)',delZ_delI, ThetaHat, Vhat_fNT, (1-valSIT)
       delZ_delN= -self%Q0*muhatNET/(2*din*Vhat_fNT)*(1-(Vhat_fNT/self%V0hat)-(Vhat_fNT/sqrt(self%V0hat*self%A0hat*din))) !delZ/delN, eq.A-5 in S16
       delQ_delI=delQ_delZ*delZ_delI !delQ/delI, eq. A-2 in S16
       !write(*,'(A,3F15.5)')'  (phy.3) delQ_delI,delQ_delZ,delZ_delI:',delQ_delI,delQ_delZ,delZ_delI
       delQ_delN=delQ_delZ*delZ_delN !delQ/delN, eq. A-3 in S16
       dI_dt = delta_parE / delta_t   !dI/dt
       dN_dt = delta_din / delta_t  !dN/dt
       !!$   dI_dt = 0.0  ! with this, the model runs: very small changes in I, during darkness, cause crashes.

       delQ_delt=delQ_delI*dI_dt + delQ_delN*dN_dt !delQ/delt, eq. A-6 in S16
       !write(*,'(A,5F20.10)')'  (phy.4) delQ_delI,dI_dt,delQ_delI*dI_dt,delQ_delN,dN_dt:',delQ_delI,dI_dt,delQ_delI*dI_dt,delQ_delN,dN_dt
       vN = mu*Q + delQ_delt !eq. A-6 in S16
       !write(*,'(A,3F15.10)')'  (phy.5) mu*Q,delQ_delI*dI_dt,delQ_delN*dN_dt:',mu*Q,delQ_delI*dI_dt,delQ_delN*dN_dt
     end if  
   end if
   
   !Calculate fluxes between pools
   f_din_phy = vN * phyC  !Eq.5 in K20
   
   ! Mortality
   mort=self%M0p * Tfac * PhyN**2
   f_phy_detn = self%Mpart  * mort
   f_phy_detc = f_phy_detn/Q              !Table 1
   f_phy_don = (1.0 - self%Mpart) * mort
   f_phy_doc = f_phy_don/Q                !Doesn't appear in K20, since Mpart=1 -> f_phy_don=0 
   
   !write(*,'(A,5F12.5)')'  (phy) dphyC*dt,vN, f_din_phy/Q, -f_phy_don/Q, -f_phy_detn/Q: ', (f_din_phy/Q - f_phy_don/Q - f_phy_detn/Q)*12,vN, f_din_phy/Q, -f_phy_don/Q, -f_phy_detn/Q
   ! Set temporal derivatives
   _SET_ODE_(self%id_phyC, mu*phyC - f_phy_don/Q - f_phy_detn/Q)  !  f_din_phy/Q - f_phy_don/Q - f_phy_detn/Q)
   !write(*,'(A,2F15.10)')'  (phy.6) phyC,delta_phyC',phyC,(mu*phyC - f_phy_don/Q - f_phy_detn/Q)*delta_t
   
   ! If externally maintained dim,dom und det pools are coupled:
   !dN/dt: based on explicit delQ/delN*deltaN/deltat (eq.44 in OptScale.v2)
   !_SET_ODE_(self%id_din, (f_don_din -f_din_phy))
   !dN/dt: after rearranging eq. 44 and isolating dN/dt (eq.45)
   _SET_ODE_(self%id_din, (f_don_din - (mu*Q+delQ_delI*dI_dt)*phyC)/(1+phyC*delQ_delN))
   _SET_ODE_(self%id_don,  f_phy_don + f_det_don - f_don_din)
   _SET_ODE_(self%id_doc,   f_phy_doc + f_det_doc - f_doc_dic)  !Eq.3b
   _SET_ODE_(self%id_detN, f_phy_detn - f_det_don)
   _SET_ODE_(self%id_detc, f_phy_detc-f_det_doc) !Sink term in Eq.2b

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_Q, Q)
   if ( self%dynQN ) then
     _SET_DIAGNOSTIC_(self%id_fQ, fQ)
   else
      _SET_DIAGNOSTIC_(self%id_d_phyN, phyN)
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
   _SET_DIAGNOSTIC_(self%id_vN, vN * secs_pr_day) !*s_p_d such that output is in d-1
   _SET_DIAGNOSTIC_(self%id_muhatG, muhatG * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_muhatNET, muhatNET * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_vNhat, vNhat * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_respN, respN * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_respChl, Rchl * secs_pr_day)
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

   end module NflexPD_phyabio_Cbased

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
