#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NflexPD_abio - abiotic component
!
! !INTERFACE:
   module NflexPD_abio
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
   type,extends(type_base_model),public :: type_NflexPD_abio
!     Variable identifiers
      type (type_state_variable_id)     :: id_din,id_don,id_detn
      type (type_dependency_id)         :: id_temp,id_parW,id_parW_dmean
      type (type_diagnostic_variable_id)   :: id_dPAR,id_dPAR_dmean

!     Model parameters
      real(rk) :: kdet,kdon

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
! !IROUTINE: Initialise the Detritus model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the NflexPD_abio namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_NflexPD_abio), intent(inout), target :: self
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
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','E/m^2/d',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_dPAR_dmean, 'PAR_dmean','E/m^2/d','photosynthetically active radiation, daily averaged')
                                     
   ! Register environmental dependencies
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
   class (type_NflexPD_abio), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: detn, don
   real(rk)                   :: f_det_don, f_don_din
   real(rk)                   :: tC,Tfac,parW,parW_dm
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
   ! so just restore the instantaneous value
   if ( parW_dm .lt. 0.0 ) then
     parW_dm=parW
     !write(*,*)'restored parW_dm'
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
   
   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR, parW * 4.6 * 1e-6 * secs_pr_day) ! mol/m2/d-1
   _SET_DIAGNOSTIC_(self%id_dPAR_dmean, parW_dm * 4.6 * 1e-6 * secs_pr_day) !mol/m2/d
    ! 1 W/m2 ≈ 4.6 μmole/m2/s: Plant Growth Chamber Handbook (chapter 1, radiation; https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC


!-----------------------------------------------------------------------
!
! !IROUTINE: Temperature dependence of rates for the FlexPFT model 
!
! !INTERFACE:
   REALTYPE function FofT(tC)
!
! !DESCRIPTION:
! Here, Arrhenius type temperature dependence is calcuated, for a reference temperature, Tr, 
! and activation energy Ea [ J / mol ]. 
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: tC            ! temperature [ degrees C ]
   REALTYPE, parameter                 :: R = 8.31446   ! gas constant [ J /mol /K ]
   REALTYPE, parameter                 :: Ea=4.82e4        ! [ J / mol ]
   REALTYPE, parameter                 :: Tr=20.0          ! [ degrees C ]
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
 
   end module NflexPD_abio

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
