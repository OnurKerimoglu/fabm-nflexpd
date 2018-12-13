module pmbench_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use pmbench_phy_classic
   use pmbench_phy_DOQ
   use pmbench_phy_IOQ
   use pmbench_phy_IOQ_Nbased
   use pmbench_abio
   use pmbench_abio_dEdt
   ! Add new PMbench modules here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: pmbench_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('phy_classic'); allocate(type_pmbench_phy_classic::model)
         case ('phy_DOQ'); allocate(type_pmbench_phy_DOQ::model)
         case ('phy_IOQ'); allocate(type_pmbench_phy_IOQ::model)
         case ('phy_IOQ_Nbased'); allocate(type_pmbench_phy_IOQ_Nbased::model)
         case ('abio'); allocate(type_pmbench_abio::model)
         case ('abio_dEdt'); allocate(type_pmbench_abio_dEdt::model)
         ! Add new PMbench models here
      end select
   end subroutine create

end module pmbench_model_library
