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

   public
   
   integer, parameter :: max_phynum = 10 !maximum number of phy. Adjust and recompile if necessary
   
   !structure to collect prey parameters (variable identifiers, stoichiometric ratios, etc)
   type,public                          :: phy_pars
     !type (type_state_variable_id)      :: id_?
     type (type_dependency_id)          :: id_dep_vN_dQdt_I,id_dep_del_phyn_din,id_dep_fdinphy
     !type (type_diagnostic_variable_id) :: id_?
     !real(rk)                           :: ?
     !logical                            :: ?
   end type

   !structure to collect phydata
   type,public                         :: phy_data
     real(rk),dimension(max_phynum)   :: f_din_phy
     real(rk),dimension(max_phynum)   :: vN_dQdt_I
     real(rk),dimension(max_phynum)   :: del_phyn_din
   end type
      
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_NflexPD_RHSext_Cbased
!     Variable identifiers
      type (type_state_variable_id)        :: id_din
      type (type_dependency_id)            :: id_dep_fabiodin
      type (type_dependency_id)            :: id_dep_fdinphy
      type (type_dependency_id)            :: id_dep_del_phyn_din,id_dep_vN_dQdt_I
      type(phy_pars),dimension(:),allocatable :: phypar
   
!     Model parameters
      logical  :: IA, DA
      integer  :: num_phy
      real(rk) :: kc

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
! !LOCAL VARIABLES
!parameters of the self structure
   integer                   :: i
   character(len=2)          :: istr !i converted to character string
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   call self%get_parameter(self%IA, 'IA','-', 'whether IA approach is used', default=.true.)
   call self%get_parameter(self%num_phy, 'numphy','-', 'number of phytoplankton instances to be coupled', default=2)
   
   ! Register state variables
   call self%register_state_variable(self%id_din,'din','mmolN/m^3','DIN concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)
   
   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_din)
   
   call self%register_dependency(self%id_dep_fabiodin, 'f_abio_din','mmolN/m^3/d',    'Bulk N flux from DOM to DIM')
   
   allocate(self%phypar(self%num_phy))
   if ( self%IA ) then !when using IA approach
     write(*,*)'  IA mode: ON'
     !call self%register_dependency(self%id_dep_del_phyn_din,'del_phyn_din', '-','Change in phyto-N per change in DIN')
     !call self%register_dependency(self%id_dep_vN_dQdt_I,'vN_dQdt_I', 'molN/m^3/d','phyC*(vN+delQ/delI*dI/dt)')
     DO i=1,self%num_phy
      write (istr,'(i0)') i
      call self%register_dependency(self%phypar(i)%id_dep_del_phyn_din,'del_phyn_din'//trim(istr), 'mmolN/m^3/s','Change in phyto-N per change in DIN for phy-'//trim(istr))
       call self%register_dependency(self%phypar(i)%id_dep_vN_dQdt_I,'vN_dQdt_I'//trim(istr), 'mmolN/m^3/s','phyC*(vN+delQ/delI*dI/dt)_for_phy'//trim(istr))
     END DO
   else
     write(*,*)'  IA mode: OFF'
     !call self%register_dependency(self%id_dep_fdinphy,'f_din_phy', 'mmolN/m^3/d','f_din_phy')
     DO i=1,self%num_phy
      write (istr,'(i0)') i
      call self%register_dependency(self%phypar(i)%id_dep_fdinphy,'f_din_phy'//trim(istr), 'mmolN/m^3/s','flux_between_din_and_phy'//trim(istr))
     END DO
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
   real(rk)                   :: f_abio_din
   real(rk)                   :: total_del_phyn_din,vN_dQdt_I
   real(rk)                   :: f_din_phy_sum,vN_dQdt_I_sum,del_phyn_din_sum
   integer                    :: i
   character(len=2)          :: istr !i converted to character string
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
   type(phy_data)             :: phydat
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   
   phydat%vN_dQdt_I=0.0
   phydat%del_phyn_din=0.0
   phydat%f_din_phy=0.0
   
   _GET_(self%id_dep_fabiodin,f_abio_din)
   
   !Collect the data from coupling targets
   DO i=1,self%num_phy
     write (istr,'(i0)') i
     if ( self%IA ) then !when using IA approach
       _GET_(self%phypar(i)%id_dep_vN_dQdt_I,phydat%vN_dQdt_I(i))
       _GET_(self%phypar(i)%id_dep_del_phyn_din,phydat%del_phyn_din(i))
       !write(*,*)'i,vN_dQdt_I,del_phyn_din',i,phydat%vN_dQdt_I(i),phydat%del_phyn_din(i)
     else
       _GET_(self%phypar(i)%id_dep_fdinphy,phydat%f_din_phy(i))
       !write(*,*)'i,f_din_phy',phydat%f_din_phy(i)
     end if
   END DO
       
   !Set the ODEs    
   if ( self%IA ) then !when using IA approach
     vN_dQdt_I_sum=sum(phydat%vN_dQdt_I)
     del_phyn_din_sum=sum(phydat%del_phyn_din)
     _SET_ODE_(self%id_din, (f_abio_din/secs_pr_day - vN_dQdt_I_sum/secs_pr_day)/(1+del_phyn_din_sum))
   else
     f_din_phy_sum=sum(phydat%f_din_phy)  
     _SET_ODE_(self%id_din, f_abio_din/secs_pr_day-f_din_phy_sum/secs_pr_day)
   end if
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

   end module NflexPD_RHSext_Cbased

!-----------------------------------------------------------------------
! Copyright Onur Kerimoglu (kerimoglu.o@gmail.com) - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
