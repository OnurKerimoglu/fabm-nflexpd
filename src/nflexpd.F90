module NflexPD_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use NflexPD_phy
   use NflexPD_nut
   use NflexPD_det
   use NflexPD_version
   ! Add new components modules here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: initialize
      procedure :: create
   end type

   type (type_factory),save,target,public :: NflexPD_model_factory

contains
   
   subroutine initialize(self)
      class (type_factory), intent(inout) :: self
      call self%register_version('NflexPD',git_commit_id//' ('//git_branch_name//' branch)')
   end subroutine initialize

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('phy'); allocate(type_NflexPD_phy::model)
         case ('nut'); allocate(type_NflexPD_nut::model)
         case ('det'); allocate(type_NflexPD_det::model)
         ! Add new components models here
      end select
   end subroutine create

end module NflexPD_model_library
