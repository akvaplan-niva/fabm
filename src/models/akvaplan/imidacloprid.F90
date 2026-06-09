#include "fabm_driver.h"

module akvaplan_imidacloprid

   ! A simple tracer model with support for sinking/floating and temperature-dependent decay.
   ! The temperature dependence of decay can be described with a Q10 formulation or an Arrhenius formulation.
   ! Copyright (C) 2016 - Akvaplan-niva

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_imidacloprid
      ! The model's own variables
      type (type_state_variable_id)            :: id_c ! tracer concentration
      type (type_bottom_state_variable_id)     :: id_c_bot !

      ! Environmental dependencies
      type (type_dependency_id) :: id_T ! temperature
      type (type_horizontal_dependency_id) :: id_shear !Bottom shear
      type (type_dependency_id) :: id_dswr
      type (type_dependency_id) :: id_h ! layer thickness
      type (type_global_dependency_id) :: id_max_dt ! model timestep

      ! Parameters
      real(rk) :: w
      logical :: do_sed
      integer :: resusp_meth
      real(rk) :: crt_shear,erate,bed_por
      real(rk) :: DT50P, DT50H, DT50M
      real(rk) :: Ia
      real(rk) :: Iz_min

      ! Diagnostics
      type (type_diagnostic_variable_id) :: id_light_out
 
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
   end type

   type (type_bulk_standard_variable),parameter :: non_water_volume_fraction = type_bulk_standard_variable(name='non_water_volume_fraction',units='-',      aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: non_water_density         = type_bulk_standard_variable(name='non_water_density',        units='kg m-3', aggregate_variable=.true.)
   
   !  private data members
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.0000000001_rk
   real(rk), parameter :: R = 8.3144598_rk    ! universal gas constant (J mol-1 K-1)
   real(rk), parameter :: Kelvin = 273.15_rk  ! offset between degrees Celsius and Kelvin


contains

   subroutine initialize(self,configunit)
      class (type_imidacloprid),intent(inout),target :: self
      integer,            intent(in)           :: configunit

      real(rk)                        :: sp_vol, sp_dens
      logical                         :: conserved
      character(len=attribute_length) :: standard_name


     ! Parameter values for degredation. Defaults are from Mould et al. 2020 
     call self%get_parameter(self%DT50P,'DT50P','','Photolysis half life in hours', default=20.38_rk)  ! Lab Irradiance Wm-2
     call self%get_parameter(self%Ia,'Ia','','Lab irradiance for photolysis half life', default=27.0_rk)  ! Lab Irradiance Wm-2 
     call self%get_parameter(self%Iz_min,'Iz_min', '', 'Lower cut off for photolysis', default=1.503_rk) !Iz min. corresponding to t50max = ln(2)/0.002536 (Dark cond. half life)
     call self%get_parameter(self%DT50H,'DT50H', '', 'Hydrolysis half life in hours', default=366.0_rk) !half life when dark, i.e. Hydrolysis (Mould 2020)
     call self%get_parameter(self%DT50M,'DT50M','', 'Microbial half life in hours', default=3840.0_rk) !Microbial half life in hours according to (Roberts, T. R. and Hutson, D. H., 1999), 720 best case 3840 worst case

      ! Parameters for sedimentation
      call self%get_parameter(self%w, 'w', 'm d-1', 'sinking velocity',                default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%do_sed,'do_sed','','sedimentation switch (.true.: include sedimentation, .false.: disregard sedimentation)',default=.false.)
      ! Include resuspension, default = 0 (no resusp): build out gradually
      call self%get_parameter(self%resusp_meth,'resusp_meth','','resuspension method (0: no resuspension)', default=0, minimum=0,maximum=1)
      !critical shear stress (bottom) for resuspension
      call self%get_parameter(self%crt_shear,    'crt_shear',    'N m-2',             'critical shear stress',                                     default=0.005_rk)
      call self%get_parameter(self%erate,        'erate',        'kg m-2 s-1',        'bed erodibility constant',                                  default=0.0_rk)
      call self%get_parameter(self%bed_por,      'bed_por',      'm3 m-3',            'porosity [vol of voids/total vol]',                         default=0.0_rk)
     
      ! Register the model's own variables.
      call self%register_state_variable(self%id_c,'c','quantity m-3','concentration',initial_value=1.0_rk,vertical_movement=-self%w,minimum=0.0_rk)   
      call self%register_state_variable(self%id_c_bot,'c_bot','quantity m-2','concentration at bottom',minimum=0.0_rk)
  
      ! Register environmental dependencies.
      call self%register_dependency(self%id_T,standard_variables%temperature)
      call self%register_dependency(self%id_shear,  standard_variables%bottom_stress)
      call self%register_dependency(self%id_dswr,standard_variables%downwelling_shortwave_flux)

      if (self%do_sed) then  
        call self%register_dependency(self%id_h,standard_variables%cell_thickness)
        call self%register_dependency(self%id_max_dt,standard_variables%maximum_time_step)
      end if

      ! Density hook based on specific volume and density of the tracer.
      call self%get_parameter(sp_vol,  'specific_volume', 'm3 quantity-1', 'specific volume', default=0.0_rk)
      call self%get_parameter(sp_dens, 'density','kg m-3', 'density (tracer mass/tracer volume)', default=0.0_rk)

      call self%add_to_aggregate_variable(non_water_volume_fraction,self%id_c,scale_factor=sp_vol)
      call self%add_to_aggregate_variable(non_water_volume_fraction,self%id_c_bot,scale_factor=sp_vol)
      call self%add_to_aggregate_variable(non_water_density,        self%id_c,scale_factor=sp_vol*sp_dens) ! converting from kg tracer/tracer_volume to kg tracer/total_volume
      call self%add_to_aggregate_variable(non_water_density,        self%id_c_bot,scale_factor=sp_vol*sp_dens)

      ! Optionally register the tracer as a "conserved quantity". This will prompt the hydrodynamic model to compute budgets across the entire domain.
      ! The name of the conserved quantity will be the short name of the model itself, with "_total" appended.
      call self%get_parameter(conserved,'conserved','','whether this is a conserved quantity (activates budget tracking in host)',default=.false.)
      if (conserved) then
         standard_name = get_safe_name(trim(self%get_path())//'_total')
         call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:),units='quantity m-3',aggregate_variable=.true.,conserved=.true.),self%id_c)
      end if

      ! Add a diagnostic for the light being used
      call self%register_diagnostic_variable(self%id_light_out, 'light_out', 'W', 'Light field used for photolysis')

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_imidacloprid),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: c, T, k, DT50R, Iz, crate

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Obtain model state and environment from FABM
         _GET_(self%id_c,c) ! tracer concentration
         _GET_(self%id_T,T) ! temperature (degree Celsius)
         _GET_(self%id_dswr,Iz) !downwelling shortwave radiation

         if (Iz>self%Iz_min) then ! Cut off 
            DT50R = self%DT50P*(self%Ia/Iz)
            k = LOG(2.0_rk)/DT50R
         else
            k = LOG(2.0_rk)/self%DT50H ! Hydrolysis applied when dark
         endif
         crate = -k/3600_rk ! adjust to per second

         ! Send instantaneous rate of change to FABM
         _SET_ODE_(self%id_c,c*crate) ! k already negative

         ! Set diagnostics
         _SET_DIAGNOSTIC_(self%id_light_out, Iz)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do
   
   
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      
      class (type_imidacloprid),intent(in) :: self
         _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c
      real(rk) :: c_bot
      real(rk) :: shear,erosion,h,max_dt
      real(rk) :: kmicrob,crate_bot
      
      ! Enter spatial loops (if any)
         _HORIZONTAL_LOOP_BEGIN_
  
          ! Environment
          _GET_(self%id_c,c)
          _GET_HORIZONTAL_(self%id_c_bot,c_bot)
          _GET_HORIZONTAL_(self%id_shear,shear)
          _GET_(self%id_h,h)
          _GET_GLOBAL_(self%id_max_dt,max_dt)

          kmicrob = LOG(2.0_rk)/self%DT50M
          crate_bot = -kmicrob/3600_rk ! adjust to per second

          ! Sedimentation
          if (self%do_sed) then
                ! Resuspension/Erosion
                select case(self%resusp_meth)

                  case(0)
                  !Do nothing
                    erosion = 0.0_rk
                  case(1)
                  !resuspension as described in Ariathurai and Arulanandan (1978). 
                  !Also used in fvcom and ROMS (Warner et al. 2008)
                    erosion = min(c_bot/max_dt, self%erate*(1.0_rk - self%bed_por)*max(shear/self%crt_shear - 1.0_rk, 0.0_rk))
                end select
                
                _SET_BOTTOM_EXCHANGE_(self%id_c,-min(h/max_dt,self%w)*c+erosion)
                _SET_BOTTOM_ODE_(self%id_c_bot, +min(h/max_dt,self%w)*c-erosion+c_bot*crate_bot)

          end if

         _HORIZONTAL_LOOP_END_
  
        end subroutine do_bottom
   
end module
