module gravity_module
  use config_module
  use variable_module
  use grid_module
  implicit none
  private
  public :: t_gravity

  type t_gravity
    private
    integer :: npv,ndv,ngrd
    real(8) :: gravity(3)
    real(8), pointer :: dv(:),grd(:)
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setdv
      procedure :: setgrd
      procedure :: getgravitysource
  end type t_gravity

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(gravity,config,grid,variable)
      implicit none
      class(t_gravity), intent(out) :: gravity
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable


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

      gravity%npv = variable%getnpv()
      gravity%ndv = variable%getndv()
      gravity%ngrd = grid%getngrd()

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(gravity)
      implicit none
      class(t_gravity), intent(inout) :: gravity

      if(associated(gravity%dv))  nullify(gravity%dv)
      if(associated(gravity%grd)) nullify(gravity%grd)

    end subroutine destruct
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

    end function getgravitysource
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module gravity_module