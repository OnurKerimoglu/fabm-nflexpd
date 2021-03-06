#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: nflexpd_abio - abiotic component
! Original Author(s): 
! S. Lan Smith 2014-12-09 (gotm-bio)
! O. Kerimoglu 2018-11-27 (first fabm version), 2020-11-19 (version in 'K20')
!
! !INTERFACE:
   module nflexpd_abio
!
! !DESCRIPTION:
! This model describes the abiotic processes
!
! !USES:
   use fabm_types
   use nflexpd_common
   use fabm_expressions

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_nflexpd_abio
!     Variable identifiers
      type (type_state_variable_id)     :: id_din,id_don,id_doc,id_detn,id_detc
      type (type_dependency_id)         :: id_temp,id_depth,id_parW,id_parW_dmean
      type (type_horizontal_dependency_id)  :: id_lat,id_depFDL
      type (type_global_dependency_id)  :: id_doy
      type (type_diagnostic_variable_id):: id_dPAR,id_dPAR_dmean
      type (type_diagnostic_variable_id):: id_fdetdon,id_fdetdoc,id_fdondin,id_fdocdic
      type (type_surface_diagnostic_variable_id):: id_dFDL
      type (type_bottom_diagnostic_variable_id):: id_detn_sed,id_detc_sed
      
!     Model parameters
      real(rk) :: w_det,kdet,kdon,par0_dt0,kc_dt0

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
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
!  Here, the nflexpd_abio namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_nflexpd_abio), intent(inout), target :: self
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
   call self%get_parameter(self%w_det,     'w_det','m d-1',    'vertical velocity (<0 for sinking)',default=-5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(kc,      'kc', 'm2 mmol-1','specific light extinction',         default=0.03_rk)
   call self%get_parameter(self%kdet,'kdet','d-1',      'sp. rate for f_det_don',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%kdon,'kdon','d-1',      'sp. rate for f_don_din',             default=0.003_rk,scale_factor=d_per_s)
   call self%get_parameter(self%par0_dt0,'par0_dt0','W m-2', 'daily average par at the surface on the first time step',  default=4.5_rk)
   call self%get_parameter(self%kc_dt0,'kc_dt0','m-1', 'attenuaton coefficient on the first time step',  default=0.2_rk)
   
   
   ! Register state variables
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
   call self%register_surface_diagnostic_variable(self%id_dFDL,'FDL','-',       'fractional day length')
   call self%register_bottom_diagnostic_variable(self%id_detn_sed,'detn_sed','mmolN/m^2/d','sedimentation rate of detN')
   call self%register_bottom_diagnostic_variable(self%id_detc_sed,'detc_sed','mmolC/m^2/d','sedimentation rate of detC')
   
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','E/m^2/d',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_dPAR_dmean, 'PAR_dmean','E/m^2/d','photosynthetically active radiation, daytime average')
   
   call self%register_diagnostic_variable(self%id_fdetdon, 'f_det_don','mmolN/m^3/d',    'bulk N flux from detritus to DOM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fdetdoc, 'f_det_doc','mmolC/m^3/d',    'bulk C flux from detritus to DOM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fdondin, 'f_don_din','mmolN/m^3/d',    'bulk N flux from DOM to DIM',           &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_fdocdic, 'f_doc_dic','mmolC/m^3/d',    'bulk C flux from DOM to DIM',           &
                                     output=output_time_step_averaged)                                  
                                     
   ! Register environmental dependencies
   call self%register_global_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
    call self%register_horizontal_dependency(self%id_depFDL, 'FDL','-',       'fractional day length')
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_parW, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_parW_dmean,temporal_mean(self%id_parW,period=1._rk*86400._rk,resolution=1._rk))
   
   
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
   class (type_nflexpd_abio), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: detn,detc,don,doc
   real(rk)                   :: f_det_don,f_det_doc,f_don_din,f_doc_dic
   real(rk)                   :: Ld,Tfac,parE,parE_dm
   real(rk)                   :: lat,depth,doy,tC,parW,parW_dm
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_detc,detc) ! detrital carbon
   _GET_(self%id_detn,detn) ! detrital nitrogen
   _GET_(self%id_doc,doc) ! dissolved organic carbon
   _GET_(self%id_don,don) ! dissolved organic nitrogen
   ! Retrieve environmental dependencies
   _GET_(self%id_depth,depth)     ! depth
   _GET_(self%id_temp,tC) ! temperature in Celcius
   
   _GET_(self%id_parW,parW) ! local photosynthetically active radiation
   _GET_(self%id_parW_dmean,parW_dm)
   
   _GET_HORIZONTAL_(self%id_depFDL,Ld)
   
   ! for the first day, the daily mean value doesn't yet exist (resulting in values in the order of -1e19), 
   ! so just restore it with a value obtained with an exp decay function to account for depth
   if ( parW_dm .lt. 0.0_rk ) then
    parW_dm=0.0_rk
    !parW_dm=self%par0_dt0*exp(-depth*self%kc_dt0)
   end if
   parE = parW * 4.6 * 1e-6 * secs_pr_day
   parE_dm= parW_dm * 4.6 * 1e-6 * secs_pr_day/Ld !division by Ld converts 24h average to daytime average
   ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   !write(*,*)'A.L173:depth,par_dm',depth,parE_dm

   !Calculate intermediate terms:
   !Temperature factor 
   Tfac = FofT(tC)
   
   !Calculate fluxes between pools
   f_det_don = self%kdet * Tfac * detn !Table 1 in K20
   f_det_doc = self%kdet * Tfac * detc !Table 1 in K20
   f_don_din = self%kdon * Tfac * don  !Table 1 in K20
   f_doc_dic = self%kdon * Tfac * doc  !Table 1 in K20
   
   ! Set temporal derivatives
   _SET_ODE_(self%id_detn, -f_det_don) !Sink term in Eq.2a
   _SET_ODE_(self%id_detc, -f_det_doc) !Sink term in Eq.2b
   _SET_ODE_(self%id_don,   f_det_don - f_don_din)  !Eq.3a
   _SET_ODE_(self%id_doc,   f_det_doc - f_doc_dic)  !Eq.3b
   _SET_ODE_(self%id_din,   f_don_din)              !Source term in Eq.4
   
   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR, parE) ! mol/m2/d-1
   _SET_DIAGNOSTIC_(self%id_dPAR_dmean, parE_dm) ! mol/m2/d-1
    ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
    
    !Export diagnostic bulk fluxes
   _SET_DIAGNOSTIC_(self%id_fdetdon, f_det_don * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_fdetdoc, f_det_doc * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_fdondin, f_don_din * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_fdocdic, f_doc_dic * secs_pr_day)
   
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
   class (type_nflexpd_abio), intent(in)     :: self
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
   _SET_SURFACE_DIAGNOSTIC_(self%id_dFDL,Ld) !Fractional day length
   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Abiotic model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_nflexpd_abio), intent(in)     :: self
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

   end module nflexpd_abio

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
