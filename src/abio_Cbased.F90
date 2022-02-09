#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NflexPD_abio_Cbased - abiotic component
! Original Author(s): S. Lan Smith 2014-12-09, O. Kerimoglu 2018-11-27
!
! !INTERFACE:
   module NflexPD_abio_Cbased
!
! !DESCRIPTION:
! This model describes the abiotic processes
!
! !USES:
   use nflexpd_common
   use fabm_types
   use fabm_expressions

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_NflexPD_abio_Cbased
!     Variable identifiers
      type (type_state_variable_id)     :: id_dic,id_din,id_don,id_doc,id_detn,id_detc
      type (type_dependency_id)         :: id_temp,id_depth,id_parW,id_parW_dmean
      type (type_horizontal_dependency_id)  :: id_lat
      type (type_global_dependency_id)  :: id_doy
      type (type_diagnostic_variable_id):: id_dPAR,id_dPAR_dmean,id_dFDL,id_dFDLdt
      !type (type_horizontal_diagnostic_variable_id):: id_dFDL
      type (type_diagnostic_variable_id):: id_fdetdon,id_fdetdoc,id_fdondin,id_fdocdic
      type (type_bottom_diagnostic_variable_id):: id_detn_sed,id_detc_sed
      
      !for saving and accessing the doy,din and par of the previous integration time step
      !for doy and delta_t, global diagnostic and dependency would be better but they don't exist
      type (type_diagnostic_variable_id)   :: id_ddoy,id_delta_t
      type (type_dependency_id)            :: id_ddoy_dep
      !for delta_par and delta_din we need 3-D diagnostics anyway
      type (type_dependency_id)            :: id_dep_delta_t,id_dep_delta_din,id_dep_delta_par,id_dep_delta_temp
      type (type_diagnostic_variable_id)   :: id_dtemp,id_ddin,id_delta_din,id_delta_par,id_delta_temp
      type (type_dependency_id)            :: id_ddin_dep,id_dtemp_dep,id_dpardm_dep
      !dependencies from the phyto components
      !type (type_diagnostic_variable_id)   :: id_total_del_phyn_din
      !type (type_dependency_id)            :: id_dep_total_del_phyn_din 
      !type (type_dependency_id)            :: id_dep_del_phyn_din 
      
!     Model parameters
      logical :: PAR_dmean_FDL,PAR_ext_inE
      real(rk) :: w_det,kdet,kdon,par0_dt0,kc_dt0

      contains

      procedure :: initialize
      !procedure :: do_surface
      procedure :: do_bottom
      procedure :: do
      

   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the Detritus model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the NflexPD_abio_Cbased namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_NflexPD_abio_Cbased), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk)            :: w_det, kc
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day and are converted here to values per second.
   call self%get_parameter(self%PAR_ext_inE, 'PAR_ext_inE','-', 'PAR provided externally are in Einsten/Quanta [mol/m2/d]', default=.false.)
   call self%get_parameter(self%PAR_dmean_FDL, 'PAR_dmean_FDL','-', 'PAR as day time average (PAR_dm/FDL, FDL: fractional day length)', default=.true.)
   call self%get_parameter(self%w_det,     'w_det','m d-1',    'vertical velocity (<0 for sinking)',default=-5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(kc,      'kc', 'm2 mmol-1','specific light extinction',         default=0.03_rk)
   call self%get_parameter(self%kdet,'kdet','d-1',      'sp. rate for f_det_don',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%kdon,'kdon','d-1',      'sp. rate for f_don_din',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%par0_dt0,'par0_dt0','W m-2', 'daily average par at the surface on the first time step',  default=4.5_rk)
   call self%get_parameter(self%kc_dt0,'kc_dt0','m-1', 'attenuaton coefficient on the first time step',  default=0.2_rk)
   
   ! Register state variables
   call self%register_state_variable(self%id_dic,'dic','mmolC/m^3','DIC concentration',     &
                                1000.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_din,'din','mmolN/m^3','DIN concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_don,'don','mmolN/m^3','DON concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_doc,'doc','mmolC/m^3','DOC concentration',     &
                                6.625_rk,minimum=0.0_rk,no_river_dilution=.true., &
                                specific_light_extinction=0.0_rk)                             
   call self%register_state_variable(self%id_detn,'detn','mmolN/m^3','Det-N concentration',    &
                                4.5_rk,minimum=0.0_rk,vertical_movement=self%w_det, &
                                specific_light_extinction=kc)
   call self%register_state_variable(self%id_detc,'detc','mmolC/m^3','Det-C concentration',    &
                                6.625_rk,minimum=0.0_rk,vertical_movement=self%w_det, &
                                specific_light_extinction=0.0_rk)                             
   
   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_dic)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_doc)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_detc)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_din)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_don)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_detn)
   
   ! Register diagnostic variables
   call self%register_bottom_diagnostic_variable(self%id_detn_sed,'detn_sed','mmolN/m^2/d','sedimentation rate of detN')
   call self%register_bottom_diagnostic_variable(self%id_detc_sed,'detc_sed','mmolC/m^2/d','sedimentation rate of detC')
   
   !call self%register_diagnostic_variable(self%id_dFDL,'FDL','-',       'fractional day length',source=source_do_surface) !,domain=domain_surface)
   call self%register_diagnostic_variable(self%id_dFDL,'FDL','-',       'fractional day length')
   call self%register_diagnostic_variable(self%id_dFDLdt,'dFDL_dt','/s',       'time derivative of fractional day length')
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','E/m^2/d',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_dPAR_dmean, 'PAR_dmean','E/m^2/d','photosynthetically active radiation averaged during day light')
   
   call self%register_diagnostic_variable(self%id_fdetdon, 'f_det_don','mmolN/m^3/d',    'bulk N flux from detritus to DOM',           &
                                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fdetdoc, 'f_det_doc','mmolC/m^3/d',    'bulk C flux from detritus to DOM',           &
                                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fdondin, 'f_don_din','mmolN/m^3/d',    'bulk N flux from DOM to DIM',           &
                                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fdocdic, 'f_doc_dic','mmolC/m^3/d',    'bulk C flux from DOM to DIM',           &
                                     output=output_instantaneous)
                                     
   ! Register environmental dependencies
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_parW_dmean,temporal_mean(self%id_parW,period=1._rk*86400._rk,resolution=1._rk))
   
   ! import individual phyN/DIN sensitivities, export total phyN/DIN sensitivity
   !call self%register_dependency(self%id_dep_del_phyn_din,'del_phyn_din','molN/molN','Rate of change in phyto1-N per change in DIN')
   !including the line above results in seg fault. Maybe an importing module cannot export at the same time?
   !call self%register_diagnostic_variable(self%id_total_del_phyn_din, 'total_del_phyn_din','molN/molN','Rate of change in total phyto-N per change in DIN')
   !call self%register_dependency(self%id_dep_total_del_phyn_din,'total_del_phyn_din', 'molN/molN','Previous value of the rate of change in total phyto-N per change in DIN')
   
   !for saving and accessing the doy,din and par of the previous integration time step
   call self%register_diagnostic_variable(self%id_ddoy,'ddoy','d', 'diagn_number_of_days_since_start_of_the_year',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_ddoy_dep,'ddoy','d','prev. val of diagn_number_of_days_since_start_of_the_year')
   
   call self%register_diagnostic_variable(self%id_delta_t,'delta_t','s','diff betw current and prev time step',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_dep_delta_t, 'delta_t','s','prev. val of diff betw current and prev time step dependency')
   
   call self%register_diagnostic_variable(self%id_dtemp,'dtemp','Degree_Celsius', 'diagn. water temperature',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_dtemp_dep,'dtemp','Degree_Celsius', 'prev. val of diagn. water temperature')
   call self%register_diagnostic_variable(self%id_delta_temp,'delta_temp','Degree_Celsius','diff betw current and prev time step',&
                     output=output_instantaneous)
                     
   call self%register_diagnostic_variable(self%id_ddin,'ddin','mmolN/m^3', 'diagn. din conc',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_ddin_dep,'ddin','mmolN/m^3', 'prev. val of diagn. din conc')
   call self%register_diagnostic_variable(self%id_delta_din,'delta_din','mmolN/m^3','diff betw current and prev time step',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_dep_delta_din, 'delta_din','mmolN/m^3','prev. val of diff in DIN betw current and prev time step')
   
   call self%register_dependency(self%id_dpardm_dep,'PAR_dmean','E/m^2/d',       'prev. val of photosynthetically active radiation, daily averaged')
   call self%register_diagnostic_variable(self%id_delta_par,'delta_par','E/m^2/d','diff betw current and prev time step',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_dep_delta_par, 'delta_par','E/m^2/d','prev. val of diff in PAR betw current and prev time step')
   call self%register_dependency(self%id_dep_delta_temp, 'delta_temp','Degree_Celsius','diff in temp betw current and prev time step')
   
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Abiotic model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_NflexPD_abio_Cbased), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: dic, din, detn, detc, don, doc
   real(rk)                   :: f_det_don, f_det_doc, f_don_din, f_doc_dic
   real(rk)                   :: Ld,dLd_dt,Tfac,parE,parE_dm
   real(rk)                   :: lat,depth,doy,tC,parW,parW_dm
   real(rk)                   :: doy_prev,delta_t
   real(rk)                   :: din_prev,delta_din
   real(rk)                   :: tC_prev,delta_temp
   real(rk)                   :: parEdm_prev,delta_parE
   real(rk)                   :: total_del_phyn_din,del_phyn_din
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_depth,depth)     ! depth
   ! Retrieve current (local) state variable values.
   _GET_(self%id_detc,detc) ! detrital carbon
   _GET_(self%id_detn,detn) ! detrital nitrogen
   _GET_(self%id_doc,doc) ! dissolved organic carbon
   _GET_(self%id_don,don) ! dissolved organic nitrogen
   ! Retrieve environmental dependencies
   _GET_(self%id_temp,tC) ! temperature in Celcius
   _GET_GLOBAL_(self%id_doy,doy) ! day of year 
   
   _GET_(self%id_parW,parW) ! local photosynthetically active radiation (PAR)
   _GET_(self%id_parW_dmean,parW_dm) !current daily average PAR
   if ( parW_dm .lt. 0.0 ) then
     parW_dm=parW
   end if
   if (self%PAR_ext_inE) then ![mol/m2/d]
     parE = parW ![mol/m2/d]
     parE_dm = parW_dm ![mol/m2/d]
   else !assume to be in !W/m2
     parE = parW * 4.6 * 1e-6* secs_pr_day ![mol/m2/d]
     parE_dm= parW_dm * 4.6 * 1e-6* secs_pr_day ![mol/m2/d]
     ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   end if
   
   !Calculate daytime average light based on fractional day length 
   if (self%PAR_dmean_FDL) then
    _GET_HORIZONTAL_(self%id_lat,lat)
    call calc_daylength(lat,doy,Ld,dLd_dt) 
   else
    Ld=1.0
    dLd_dt=0.0
   end if
   !Convert average irradiance throughout the day to average irradiance during day light, as needed by the phy module
   parE_dm=parE_dm/Ld ![mol/m2/d]
   
   !For providing the delta_t,delta_din and delta_par between the current and previous time step
   _GET_(self%id_ddoy_dep,doy_prev)  ! day of year at the previous time step
   !write(*,*)' (abio.1) doy_prev(s),doy(s),Ld',doy_prev*secs_pr_day,doy*secs_pr_day,Ld
   !Access the par and din at the previous time step and set the diagnostic only if the time step has really advanced
   
   !Access the values at the prev. time step as recorded by the diagnostic variables
   if (doy .ne. doy_prev) then !i.e., if it's a real time step (and, e.g., not a 'fake' step of a RK4 ode method)
     
     !Access the values at the prev. time step as recorded by the diagnostic variables
     _GET_(self%id_dic,dic) ! dic
     _GET_(self%id_din,din) ! din
     _GET_(self%id_dtemp_dep,tC_prev)
     _GET_(self%id_ddin_dep,din_prev)
     _GET_(self%id_dPARdm_dep,parEdm_prev) !mol/m2/d
     
     !write(*,*)'doy_prev,din_prev,parEdm_prev',doy_prev,din_prev,parEdm_prev 

     !in the first time step, strange things may happen, as the diagnostics are not available yet
     if (doy_prev .lt. 0.0) then
       doy_prev = -1.0 ! just an arbitrary finite number, as the delta_din&par will be 0      
       tC_prev=tC ! such that delta_tC=0
       din_prev=din ! such that delta_din=0
       parEdm_prev=parE_dm ! such that delta_par=0
     end if
     
     !calculate the deltas
     delta_t=(doy-doy_prev)*secs_pr_day !days to secs
     delta_temp=tC-tC_prev
     !write(*,*)'abio:tC,tC_prev',tC,tC-delta_temp
     delta_din=din-din_prev      
     delta_parE= parE_dm-parEdm_prev !mol/m2/d
     
     !write(*,*)' (abio.L273) pardm_prev,pardm,delta_par,delta_par,dI/dt',parEdm_prev,parE_dm,delta_parE,delta_parE/(delta_t/secs_pr_day)
     
     !set the diagnostics
     !_GET_(self%id_dep_del_phyn_din,del_phyn_din)
     !total_del_phyn_din=del_phyn_din !temporary shortcut for 1-species case for experimental purposes
     !_SET_DIAGNOSTIC_(self%id_total_del_phyn_din,total_del_phyn_din)
     !write(*,*)' (abio.3a) del_phy1n_din, tot_del_phyn_din',del_phyn_din,total_del_phyn_din
     
     ! Export diagnostic variables
     _SET_DIAGNOSTIC_(self%id_dFDL,Ld) !Fractional day length
     _SET_DIAGNOSTIC_(self%id_dFDLdt,dLd_dt) !/s Time derivative of the fractional day length
     _SET_DIAGNOSTIC_(self%id_ddoy,doy)
     _SET_DIAGNOSTIC_(self%id_ddin, din)
     _SET_DIAGNOSTIC_(self%id_dtemp, tC)
     _SET_DIAGNOSTIC_(self%id_dPAR, parE) ! mol/m2/d
     _SET_DIAGNOSTIC_(self%id_dPAR_dmean, parE_dm) !mol/m2/d
     
     _SET_DIAGNOSTIC_(self%id_delta_t,delta_t)
     _SET_DIAGNOSTIC_(self%id_delta_temp,delta_temp)
     _SET_DIAGNOSTIC_(self%id_delta_din,delta_din)
     _SET_DIAGNOSTIC_(self%id_delta_par,delta_parE) !mol/m2/d
   else
     _GET_(self%id_dep_delta_t,delta_t) !secs
     _GET_(self%id_dep_delta_din,delta_din)
     _GET_(self%id_dep_delta_par,delta_parE) !mol/m2/d
     _GET_(self%id_dep_delta_temp,delta_temp)
     
     !_GET_(self%id_dep_del_phyn_din,del_phyn_din)
     !total_del_phyn_din=del_phyn_din !temporary shortcut for 1-species case for experimental purposes
     !write(*,*)' (abio.3b) prev_tot_del_phyn_din',total_del_phyn_din
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
   _SET_ODE_(self%id_detn, -f_det_don)
   _SET_ODE_(self%id_detc, -f_det_doc) !Sink term in Eq.2b
   _SET_ODE_(self%id_don,   f_det_don - f_don_din)
   _SET_ODE_(self%id_doc,  f_det_doc - f_doc_dic)  !Eq.3b
   _SET_ODE_(self%id_dic, f_doc_dic)
   
   !Standard: i.e., DA or FS approaches:
   !_SET_ODE_(self%id_din,   f_don_din)
   !IA approach:
   !_SET_ODE_(self%id_din,   f_don_din/(1.0+total_del_phyn_din)) ! (eq A-8)
   
   !Export diagnostic bulk fluxes
   _SET_DIAGNOSTIC_(self%id_fdondin, f_don_din * secs_pr_day) !mmolN m-3 d-1
   _SET_DIAGNOSTIC_(self%id_fdocdic, f_doc_dic * secs_pr_day) !mmolC m-3 d-1
   _SET_DIAGNOSTIC_(self%id_fdetdon, f_det_don * secs_pr_day) !mmolN m-3 d-1
   _SET_DIAGNOSTIC_(self%id_fdetdoc, f_det_doc * secs_pr_day) !mmolC m-3 d-1
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
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
!    class (type_NflexPD_abio_Cbased), intent(in)     :: self
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
!    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dFDL,Ld) !Fractional day length
!    
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
   class (type_NflexPD_abio_Cbased), intent(in)     :: self
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



 
   end module NflexPD_abio_Cbased

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
