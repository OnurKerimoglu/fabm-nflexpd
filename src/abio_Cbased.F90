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
   use fabm_types
   use fabm_expressions

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_NflexPD_abio_Cbased
!     Variable identifiers
      type (type_state_variable_id)     :: id_din,id_don,id_detn
      type (type_dependency_id)         :: id_temp,id_depth,id_parW,id_parW_dmean
      type (type_horizontal_dependency_id)  :: id_lat
      type (type_global_dependency_id)  :: id_doy
      type (type_diagnostic_variable_id):: id_dPAR,id_dPAR_dmean
      type (type_horizontal_diagnostic_variable_id):: id_dFDL
      
      !for saving and accessing the doy,din and par of the previous integration time step
      !for doy and delta_t, global diagnostic and dependency would be better but they don't exist
      type (type_diagnostic_variable_id)   :: id_ddoy,id_delta_t
      type (type_dependency_id)            :: id_ddoy_dep
      !for delta_par and delta_din we need 3-D diagnostics anyway
      type (type_dependency_id)            :: id_dep_delta_t,id_dep_delta_din,id_dep_delta_par
      type (type_diagnostic_variable_id)   :: id_ddin,id_delta_din,id_delta_par
      type (type_dependency_id)            :: id_ddin_dep,id_dpardm_dep
      !dependencies from the phyto components
      type (type_dependency_id)            :: id_del_phyn_din_dep !experimental 1-phytoplanton case
      type (type_diagnostic_variable_id)   :: id_total_del_phyn_din
      type (type_dependency_id)            :: id_total_del_phyn_din_dep !for being able to access the previous value
      
!     Model parameters
      real(rk) :: kdet,kdon,par0_dt0,kc_dt0

      contains

      procedure :: initialize
      procedure :: do_surface
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
   call self%get_parameter(w_det,     'w_det','m d-1',    'vertical velocity (<0 for sinking)',default=-5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(kc,      'kc', 'm2 mmol-1','specific light extinction',         default=0.03_rk)
   call self%get_parameter(self%kdet,'kdet','d-1',      'sp. rate for f_det_don',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%kdon,'kdon','d-1',      'sp. rate for f_don_din',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%par0_dt0,'par0_dt0','W m-2', 'daily average par at the surface on the first time step',  default=4.5_rk)
   call self%get_parameter(self%kc_dt0,'kc_dt0','m-1', 'attenuaton coefficient on the first time step',  default=0.2_rk)
   
   ! Register state variables
   call self%register_state_variable(self%id_din,'din','mmolN/m^3','DIN concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_don,'don','mmolN/m^3','DON concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_detn,'detn','mmolN/m^3','Det-N concentration',    &
                                4.5_rk,minimum=0.0_rk,vertical_movement=w_det, &
                                specific_light_extinction=kc)
   
   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_din)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_don)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_detn)
   
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dFDL,'FDL','-',       'fractional day length',source=source_do_surface) !,domain=domain_surface)
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','E/m^2/d',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_dPAR_dmean, 'PAR_dmean','E/m^2/s','photosynthetically active radiation, daily averaged')
                                     
   ! Register environmental dependencies
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_parW_dmean,temporal_mean(self%id_parW,period=1._rk*86400._rk,resolution=1._rk))
   
   ! import individual phyN/DIN sensitivities, export total phyN/DIN sensitivity
   !call self%register_dependency(self%id_del_phyn_din_dep,'del_phyn_din','molN/molN','Rate of change in phyto1-N per change in DIN')
   !including the line above results in seg fault. Maybe an importing module cannot export at the same time?
   !call self%register_diagnostic_variable(self%id_total_del_phyn_din, 'total_del_phyn_din','molN/molN','Rate of change in total phyto-N per change in DIN')
   !call self%register_dependency(self%id_total_del_phyn_din_dep,'total_del_phyn_din', 'molN/molN','Previous value of the rate of change in total phyto-N per change in DIN')
   
   !for saving and accessing the doy,din and par of the previous integration time step
   call self%register_diagnostic_variable(self%id_ddoy,'ddoy','d', 'diagn_number_of_days_since_start_of_the_year',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_ddoy_dep,'ddoy','d','prev. val of diagn_number_of_days_since_start_of_the_year')
   
   call self%register_diagnostic_variable(self%id_delta_t,'delta_t','s','diff betw current and prev time step',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_dep_delta_t, 'delta_t','s','prev. val of diff betw current and prev time step dependency')
   
   
   call self%register_diagnostic_variable(self%id_ddin,'ddin','mmolN/m^3', 'diagn. din conc',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_ddin_dep,'ddin','mmolN/m^3', 'prev. val of diagn din conc')
   
   call self%register_diagnostic_variable(self%id_delta_din,'delta_din','mmolN/m^3','diff betw current and prev time step',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_dep_delta_din, 'delta_din','mmolN/m^3','prev. val of diff in DIN betw current and prev time step')
   
   call self%register_dependency(self%id_dpardm_dep,'PAR_dmean','E/m^2/s',       'prev. val of photosynthetically active radiation, daily averaged')
   call self%register_diagnostic_variable(self%id_delta_par,'delta_par','E/m^2/s','diff betw current and prev time step',&
                     output=output_instantaneous)
   call self%register_dependency(self%id_dep_delta_par, 'delta_par','E/m^2/s','prev. val of diff in PAR betw current and prev time step')
   
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
   real(rk)                   :: din, detn, don
   real(rk)                   :: f_det_don, f_don_din
   real(rk)                   :: Ld,Tfac,parE,parE_dm
   real(rk)                   :: lat,depth,doy,tC,parW,parW_dm
   real(rk)                   :: doy_prev,delta_t
   real(rk)                   :: din_prev,delta_din
   real(rk)                   :: parEdm_prev,delta_par
   real(rk)                   :: total_del_phyn_din,del_phyn_din
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_detn,detn) ! detrital nitrogen
   _GET_(self%id_don,don) ! dissolved organic nitrogen
   ! Retrieve environmental dependencies
   _GET_(self%id_temp,tC) ! temperature in Celcius
   
   _GET_(self%id_parW,parW) ! local photosynthetically active radiation (PAR)
   _GET_(self%id_parW_dmean,parW_dm) !current daily average PAR
   parE = parW * 4.6 * 1e-6 ![mol/m2/s]
   parE_dm= parW_dm * 4.6 * 1e-6  ![mol/m2/s]
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   if (parE_dm .lt. 0.0) then
     parE_dm=0.0
   end if
   !For providing the delta_t,delta_din and delta_par between the current and previous time step
   
   _GET_(self%id_ddoy_dep,doy_prev)  ! day of year at the previous time step
   _GET_GLOBAL_(self%id_doy,doy)  ! day of year
   write(*,*)' (abio.1) doy_prev(s),doy(s)',doy_prev*secs_pr_day,doy*secs_pr_day
   !Access the par and din at the previous time step and set the diagnostic only if the time step has really advanced
   if (doy .gt. doy_prev) then
     
     !Access the values at the prev. time step as recorded by the diagnostic variables
     _GET_(self%id_din,din) ! din
     _GET_(self%id_ddin_dep,din_prev)
     _GET_(self%id_dPARdm_dep,parEdm_prev) !mol/m2/s
     
     !in the first time step, strange things may happen, as the diagnostics are not available yet
     !if (parEdm_prev .lt.  0.01/secs_pr_day) then
     if (doy_prev .lt. 0.0) then
       doy_prev = -1.0 ! just an arbitrary finite number, as the delta_din&par will be 0      
       din_prev=din ! such that delta_din=0
       parEdm_prev=parE_dm ! such that delta_par=0
     end if
     
     !calculate the deltas
     delta_t=(doy-doy_prev)*secs_pr_day !days to secs
     delta_din=din-din_prev      
     delta_par= parE_dm-parEdm_prev !mol/m2/s
     
     !assume no change in par within the first day
     if (doy .lt. 1.0+3*delta_t/86400) then
       delta_par=0.0
     end if
     
!     write(*,'(A,2F12.5,A,2F12.5)')' (abio.2) parE_prev,parE',parE_prev,parE,'  din_prev,din',din_prev,din
     write(*,*)' (abio.2) pardm_prev,pardm,delta_par',parEdm_prev,parE_dm,delta_par,'  din_prev,din',din_prev,din
     
     !set the diagnostics
     del_phyn_din=0.0_rk
     !_GET_(self%id_del_phyn_din_dep,del_phyn_din)
     total_del_phyn_din=del_phyn_din !temporary shortcut for 1-species case for experimental purposes
     !_SET_DIAGNOSTIC_(self%id_total_del_phyn_din,total_del_phyn_din)
     write(*,*)' (abio.3a) del_phy1n_din, tot_del_phyn_din',del_phyn_din,total_del_phyn_din
     
     ! Export diagnostic variables
     _SET_DIAGNOSTIC_(self%id_ddoy,doy)
     _SET_DIAGNOSTIC_(self%id_ddin, din)
     _SET_DIAGNOSTIC_(self%id_dPAR, parE) ! mol/m2/s
     _SET_DIAGNOSTIC_(self%id_dPAR_dmean, parE_dm) !mol/m2/s
     
     _SET_DIAGNOSTIC_(self%id_delta_t,delta_t)
     _SET_DIAGNOSTIC_(self%id_delta_din,delta_din)
     _SET_DIAGNOSTIC_(self%id_delta_par,delta_par) !mol/m2/s
   else
     !_GET_(self%id_par_dmean,par_dm) !in molE/m2/s
     _GET_(self%id_parW_dmean,parW_dm) !current daily average PAR
     parE = parW * 4.6 * 1e-6 ![mol/m2/s]
     parE_dm= parW_dm * 4.6 * 1e-6  ![mol/m2/s]
     _GET_(self%id_dep_delta_t,delta_t)
     _GET_(self%id_dep_delta_din,delta_din)
     _GET_(self%id_dep_delta_par,delta_par) !mol/m2/s
     !_GET_(self%id_total_del_phyn_din_dep,total_del_phyn_din)
     total_del_phyn_din=0.0
     write(*,*)' (abio.3b) prev_tot_del_phyn_din',total_del_phyn_din
   end if
   
   !Calculate intermediate terms:
   !Temperature factor 
   Tfac = FofT(tC)
   
   !Calculate fluxes between pools
   f_det_don = self%kdet * Tfac * detn 
   f_don_din = self%kdon * Tfac * don
   
   ! Set temporal derivatives
   _SET_ODE_(self%id_detn, -f_det_don)
   _SET_ODE_(self%id_don,   f_det_don - f_don_din)
   !_SET_ODE_(self%id_din,   f_don_din)
   _SET_ODE_(self%id_din,   f_don_din/(1.0+total_del_phyn_din)) ! (eq A-8)
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Abiotic model
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !INPUT PARAMETERS:
   class (type_NflexPD_abio_Cbased), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk)                   :: lat,doy,Ld
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_
   
   _GET_HORIZONTAL_(self%id_lat,lat)
   _GET_GLOBAL_(self%id_doy,doy)
   Ld=FDL(lat,doy)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dFDL,Ld) !Fractional day length
   
   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
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
 
 
 !-----------------------------------------------------------------------
 
   end module NflexPD_abio_Cbased

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
