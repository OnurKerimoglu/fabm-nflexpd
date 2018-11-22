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

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_NflexPD_abio
!     Variable identifiers
      type (type_state_variable_id)     :: id_din
      type (type_state_variable_id)     :: id_detn

!     Model parameters
      real(rk) :: rdn

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd

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
   call self%get_parameter(self%rdn,'rdn','d-1',      'remineralization rate',             default=0.003_rk,scale_factor=d_per_s)

   ! Register state variables
   call self%register_state_variable(self%id_din,'din','mmol m-3','concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_detn,'detn','mmol m-3','concentration',    &
                                4.5_rk,minimum=0.0_rk,vertical_movement=w_det, &
                                specific_light_extinction=kc)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_din)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_detn)

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
   real(rk)                   :: detn
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_detn,detn) ! detrital nitrogen
   
   ! Set temporal derivatives
   !Mineralization (det->dim)
   _SET_ODE_(self%id_detn,-self%rdn*detn)
   _SET_ODE_(self%id_din, self%rdn*detn)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Abiotic model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!
! !INPUT PARAMETERS:
   class (type_NflexPD_abio), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
! !LOCAL VARIABLES:
   real(rk)                   :: detn
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_detn,detn) ! detritus

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
   ! are automatically assigned to pp(j,i) as well.
   _SET_DD_SYM_(self%id_detn,self%id_din,self%rdn*detn)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd
!EOC

!-----------------------------------------------------------------------

   end module NflexPD_abio

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
