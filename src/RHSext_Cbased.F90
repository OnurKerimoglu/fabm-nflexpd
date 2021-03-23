#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NflexPD_RHSext_Cbased - RHSext component
! Original Author(s): O. Kerimoglu 2021-03-23
!
! !INTERFACE:
   module NflexPD_RHSext_Cbased
!
! !DESCRIPTION:
! This model collects and sets Right-Hand-Side of Dissolved Inorganic Nutrients
!
! !USES:
   use fabm_types
   use fabm_expressions

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_NflexPD_RHSext_Cbased
!     Variable identifiers
      type (type_state_variable_id)        :: id_din
      type (type_dependency_id)            :: id_dep_fdondin
      type (type_dependency_id)            :: id_dep_fdinphy
      type (type_dependency_id)            :: id_dep_del_phyn_din,id_dep_vN_dQdt_I
      
!     Model parameters
      logical  :: IA

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
!  Here, the NflexPD_RHSext_Cbased namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_NflexPD_RHSext_Cbased), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   call self%get_parameter(self%IA, 'IA','-', 'whether IA approach is used', default=.true.)
   
   ! Register state variables
   call self%register_state_variable(self%id_din,'din','mmolN/m^3','DIN concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   
   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_din)
   
   call self%register_dependency(self%id_dep_fdondin, 'f_don_din','mmolN/m^3/d',    'Bulk N flux from DOM to DIM')
   
   if ( self%IA ) then !when using IA approach
     write(*,*)'  IA mode: ON'
     call self%register_dependency(self%id_dep_del_phyn_din,'del_phyn_din', '-','Change in phyto-N per change in DIN')
     call self%register_dependency(self%id_dep_vN_dQdt_I,'vN_dQdt_I', 'molN/m^3/d','phyC*(vN+delQ/delI*dI/dt)')
   else
     write(*,*)'  IA mode: OFF'
     call self%register_dependency(self%id_dep_fdinphy,'f_din_phy', 'mmolN/m^3/d','f_din_phy')
   end if  
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of RHSext model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_NflexPD_RHSext_Cbased), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: f_don_din
   real(rk)                   :: total_del_phyn_din,vN_dQdt_I
   real(rk)                   :: f_din_phy
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   _GET_(self%id_dep_fdondin,f_don_din)
   
   if ( self%IA ) then !when using IA approach
     _GET_(self%id_dep_del_phyn_din,total_del_phyn_din)
     _GET_(self%id_dep_vN_dQdt_I,vN_dQdt_I) !vN_dQdt = (vN+delQ_delI*dI_dt)*phyC
     !_SET_ODE_(self%id_din, - (vN+delQ_delI*dI_dt)*phyC/(1+total_del_phyn_din))
     !_SET_ODE_(self%id_din,   (f_don_din )/(1.0+total_del_phyn_din)) ! (eq A-8)
     _SET_ODE_(self%id_din, (f_don_din/secs_pr_day - vN_dQdt_I/secs_pr_day)/(1+total_del_phyn_din))
   else
     _GET_(self%id_dep_fdinphy,f_din_phy)
     _SET_ODE_(self%id_din, f_don_din/secs_pr_day-f_din_phy/secs_pr_day)
   end if
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

   end module NflexPD_RHSext_Cbased

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
