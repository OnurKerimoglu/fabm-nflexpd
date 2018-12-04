module nflexpd_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use nflexpd_phy_DOQ
   use nflexpd_phy_IOQ
   use nflexpd_phy_IOQ_Nbased
   use nflexpd_abio
   use nflexpd_abio_dEdt
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
         case ('phy_DOQ'); allocate(type_nflexpd_phy_DOQ::model)
         case ('phy_IOQ'); allocate(type_nflexpd_phy_IOQ::model)
         case ('phy_IOQ_Nbased'); allocate(type_nflexpd_phy_IOQ_Nbased::model)
         case ('abio'); allocate(type_nflexpd_abio::model)
         case ('abio_dEdt'); allocate(type_nflexpd_abio_dEdt::model)
         ! Add new NflexPD models here
      end select
   end subroutine create

end module nflexpd_model_library
