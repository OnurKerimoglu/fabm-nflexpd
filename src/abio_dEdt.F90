#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: PMbench_abio_dEdt - abiotic component that provides the delta_N,delta_I,delta_t required by the IOQ model
! Original Author(s): O. Kerimoglu 2018-11-27
!
! !INTERFACE:
   module PMbench_abio_dEdt
!
! !DESCRIPTION:
! This model describes the abiotic processes
!
! !USES:
   use fabm_types
   use PMbench_common
   use fabm_expressions

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_PMbench_abio_dEdt
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
      type (type_diagnostic_variable_id)   :: id_ddin,id_delta_din,id_delta_par
      type (type_dependency_id)            :: id_ddin_dep,id_dpar_dep
      
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
!  Here, the PMbench_abio_dEdt namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_PMbench_abio_dEdt), intent(inout), target :: self
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
   call self%register_horizontal_diagnostic_variable(self%id_dFDL,'FDL','-',       'fractional day length')
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','E/m^2/d',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_dPAR_dmean, 'PAR_dmean','E/m^2/d','photosynthetically active radiation, daily averaged')
                                     
   ! Register environmental dependencies
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_parW_dmean,temporal_mean(self%id_parW,period=1._rk*86400._rk,resolution=1._rk))
   
   
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
   
   call self%register_dependency(self%id_dpar_dep,'PAR','E/m^2/d',       'diagn. photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_delta_par,'delta_par','E/m^2/d','diff betw current and prev time step',&
                     output=output_instantaneous)
   
   
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
   class (type_PMbench_abio_dEdt), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: din, detn, don
   real(rk)                   :: f_det_don, f_don_din
   real(rk)                   :: Ld,Tfac,parE,parE_dm
   real(rk)                   :: lat,depth,doy,tC,parW,parW_dm
   real(rk)                   :: doy_prev,delta_t
   real(rk)                   :: din_prev,delta_din
   real(rk)                   :: parE_prev,delta_par
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
   
   _GET_(self%id_parW,parW) ! local photosynthetically active radiation
   _GET_(self%id_parW_dmean,parW_dm)
   ! for the first day, the daily mean value doesn't yet exist (with values in the order of -1e19), 
   ! so just restore it with a value obtained with an exp decay function to account for depth
   if ( parW_dm .lt. 0.0 ) then 
     _GET_(self%id_depth,depth)
     parW_dm=self%par0_dt0*exp(-depth*self%kc_dt0)
     !write(*,*)'*** depth,parW_dm',depth,parW_dm
   !else
   !write(*,*)'          parW_dm',parW_dm
   end if
   parE = parW * 4.6 * 1e-6 * secs_pr_day
   parE_dm= parW_dm * 4.6 * 1e-6 * secs_pr_day
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   
   !For providing the delta_t,delta_din and delta_par between the current and previous time step
   
   _GET_(self%id_ddoy_dep,doy_prev)  ! day of year at the previous time step
   _GET_GLOBAL_(self%id_doy,doy)  ! day of year
   
   !write(*,*)'.'
   write(*,'(A,2F10.1)')'doy_prev(s),doy(s)',doy_prev*secs_pr_day,doy*secs_pr_day
   !Access the par and din at the previous time step and set the diagnostic only if the time step has really advanced
   if (doy .gt. doy_prev) then
     delta_t=(doy-doy_prev)*secs_pr_day !days to secs
     
     _GET_(self%id_din,din) ! din
     _GET_(self%id_ddin_dep,din_prev)
     delta_din= din-din_prev
     _GET_(self%id_dPAR_dep,parE_prev)
     delta_par= parE-parE_prev
     
     !write(*,'(A,2F12.5,A,2F12.5)')'  parE_prev,parE',parE_prev,parE,'  din_prev,din',din_prev,din
     
     !in the first time step, strange things may happen, as the diagnostics are not available yet
     !if (delta_t .gt. 1e10) then
     ! delta_t=360
     !end if
     
     !set the diagnostics
     
     ! Export diagnostic variables
     _SET_DIAGNOSTIC_(self%id_ddoy,doy)
     _SET_DIAGNOSTIC_(self%id_ddin, din)
     _SET_DIAGNOSTIC_(self%id_dPAR, parE) ! mol/m2/d-1
     _SET_DIAGNOSTIC_(self%id_dPAR_dmean, parE_dm) !mol/m2/d
     
     _SET_DIAGNOSTIC_(self%id_delta_t,delta_t)
     _SET_DIAGNOSTIC_(self%id_delta_din,delta_din)
     _SET_DIAGNOSTIC_(self%id_delta_par,delta_par)
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
   _SET_ODE_(self%id_din,   f_don_din)
   
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
   class (type_PMbench_abio_dEdt), intent(in)     :: self
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
 
   end module PMbench_abio_dEdt

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
