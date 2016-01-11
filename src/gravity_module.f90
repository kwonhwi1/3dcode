module gravity_module
  use config_module
  use grid_module
  implicit none
  private
  public :: t_gravity

  type t_gravity
    private
    integer :: stencil,npv,ndv,ngrd
    real(8) :: gravity(3)
    real(8), pointer :: dv(:),grd(:),pv(:,:)
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setpv
      procedure :: setdv
      procedure :: setgrd
      procedure :: getgravitysource
  end type t_gravity

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(gravity,config,grid)
      implicit none
      class(t_gravity), intent(out) :: gravity
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid

      gravity%stencil = config%getstencil()
      gravity%npv = config%getnpv()
      gravity%ndv = config%getndv()

      select case(config%getgravity())
      case(-1)
        gravity%gravity =(/-9.80665d0,0.d0,0.d0/)
      case(1)
        gravity%gravity =(/9.80665d0,0.d0,0.d0/)
      case(-2)
        gravity%gravity =(/0.d0,-9.80665d0,0.d0/)
      case(2)
        gravity%gravity =(/0.d0,9.80665d0,0.d0/)
      case(-3)
        gravity%gravity =(/0.d0,0.d0,-9.80665d0/)
      case(3)
        gravity%gravity =(/0.d0,0.d0,9.80665d0/)
      case default
        gravity%gravity =(/0.d0,0.d0,0.d0/)
      end select

      gravity%ngrd = grid%getngrd()

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(gravity)
      implicit none
      class(t_gravity), intent(inout) :: gravity

      if(associated(gravity%pv))  nullify(gravity%pv)
      if(associated(gravity%dv))  nullify(gravity%dv)
      if(associated(gravity%grd)) nullify(gravity%grd)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(gravity,pv)
      implicit none
      class(t_gravity), intent(inout) :: gravity
      real(8), intent(in), target :: pv(gravity%stencil,gravity%npv)
      
      gravity%pv => pv
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(gravity,dv)
      implicit none
      class(t_gravity), intent(inout) :: gravity
      real(8), intent(in), target :: dv(gravity%ndv)

      gravity%dv => dv

    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(gravity,grd)
      implicit none
      class(t_gravity), intent(inout) :: gravity
      real(8), intent(in), target :: grd(gravity%ngrd)

      gravity%grd => grd

    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getgravitysource(gravity)
      implicit none
      class(t_gravity), intent(in) :: gravity
      real(8) :: getgravitysource(gravity%npv)

      getgravitysource = 0.d0
      getgravitysource(2:4) = gravity%dv(1)*gravity%grd(1)*gravity%gravity(1:3)
      getgravitysource(5) = gravity%dv(1)*gravity%grd(1)*(gravity%gravity(1)*gravity%pv(2,2) &
                                                        + gravity%gravity(2)*gravity%pv(2,3) &
                                                        + gravity%gravity(3)*gravity%pv(2,4))

    end function getgravitysource
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module gravity_module