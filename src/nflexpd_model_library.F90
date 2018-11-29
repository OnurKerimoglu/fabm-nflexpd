module nflexpd_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use nflexpd_phy
   use nflexpd_phy_Cbased
   use nflexpd_abio
   use nflexpd_abio_Cbased
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
         case ('phy'); allocate(type_nflexpd_phy::model)
         case ('phy_Cbased'); allocate(type_nflexpd_phy_Cbased::model)
         case ('abio'); allocate(type_nflexpd_abio::model)
         case ('abio_Cbased'); allocate(type_nflexpd_abio_Cbased::model)
         ! Add new NflexPD models here
      end select
   end subroutine create

end module nflexpd_model_library
