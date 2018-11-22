module nflexpd_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use nflexpd_det
   use nflexpd_nut
   use nflexpd_phy
   ! Add new NflexPD modules here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: nflexpd_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('nut'); allocate(type_nflexpd_nut::model)
         case ('phy'); allocate(type_nflexpd_phy::model)
         case ('det'); allocate(type_nflexpd_det::model)
         ! Add new NflexPD models here
      end select
   end subroutine create

end module nflexpd_model_library
