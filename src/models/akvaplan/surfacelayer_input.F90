#include "fabm_driver.h"

module akvaplan_surfacelayer_input

   ! A module for inserting a 2d field (flux_in) across the surface layer to a  specified depth (d) and making it available to other fabm modules (flux_out)
   ! If the water column is shallower than d then it will NOT correct for this, so the total flux_out will be less than flux_in

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_surfacelayer_input
      type (type_surface_dependency_id)    :: id_flux_in
      type (type_dependency_id)            :: id_centre_depth, id_layer_thickness 
      type (type_diagnostic_variable_id)   :: id_flux_out

      ! Parameters
      real(rk) :: d ! Target depth over which to split the flux

   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self,configunit)
      class (type_surfacelayer_input),intent(inout),target :: self
      integer,                     intent(in)           :: configunit

      ! Obtain parameter values
      call self%get_parameter(self%d, 'd', 'm', 'depth (below surface) over which to split the flux')

      ! Register dependencies
      call self%register_dependency(self%id_flux_in,'flux_in','quantity m-2 s-1','flux into model (2d)')
      call self%register_dependency(self%id_centre_depth,standard_variables%depth)
      call self%register_dependency(self%id_layer_thickness,standard_variables%cell_thickness)

      ! And output variable
      call self%register_diagnostic_variable(self%id_flux_out,'flux_out','quantity m-3 s-1','depth-explicit flux output')

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_surfacelayer_input),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: flux_in,centre_depth,layer_thickness,adj_thickness

      _LOOP_BEGIN_
         _GET_(self%id_centre_depth,centre_depth)
         _GET_(self%id_layer_thickness,layer_thickness)

         if (self%d>(centre_depth-0.5*layer_thickness) .and. self%d<=(centre_depth+0.5*layer_thickness)) then
            _GET_SURFACE_(self%id_flux_in,flux_in)
            adj_thickness = self%d - (centre_depth-0.5*layer_thickness)
            _SET_DIAGNOSTIC_(self%id_flux_out,(flux_in*(adj_thickness/self%d))/layer_thickness)
         else if (self%d>(centre_depth-0.5*layer_thickness)) then
            _GET_SURFACE_(self%id_flux_in,flux_in)
            _SET_DIAGNOSTIC_(self%id_flux_out,(flux_in*(layer_thickness/self%d))/layer_thickness) 
         else
            _SET_DIAGNOSTIC_(self%id_flux_out,0)
         end if
         
      _LOOP_END_

   end subroutine do

end module
