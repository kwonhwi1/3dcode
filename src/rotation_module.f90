module rotation_module
  use config_module
  use grid_module
  implicit none
  private
  public :: t_rotation

  type t_rotation
    private
    integer :: stencil,npv,ndv,ngrd
    real(8) :: omega(3)
    real(8), pointer :: pv(:,:),dv(:),grd(:)
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setpv
      procedure :: setdv
      procedure :: setgrd
      procedure :: getrotationsource
  end type t_rotation

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(rotation,config,grid)
      implicit none
      class(t_rotation), intent(out) :: rotation
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid

      rotation%stencil = config%getstencil()
      rotation%npv = config%getnpv()
      rotation%ndv = config%getndv()
      rotation%omega = config%getomega()
      rotation%ngrd = grid%getngrd()

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(rotation)
      implicit none
      class(t_rotation), intent(inout) :: rotation

      if(associated(rotation%pv)) nullify(rotation%pv)
      if(associated(rotation%dv))  nullify(rotation%dv)
      if(associated(rotation%grd)) nullify(rotation%grd)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(rotation,pv)
      implicit none
      class(t_rotation), intent(inout) :: rotation
      real(8), intent(in), target :: pv(rotation%stencil,rotation%npv)

      rotation%pv => pv

    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(rotation,dv)
      implicit none
      class(t_rotation), intent(inout) :: rotation
      real(8), intent(in), target :: dv(rotation%ndv)

      rotation%dv => dv

    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(rotation,grd)
      implicit none
      class(t_rotation), intent(inout) :: rotation
      real(8), intent(in), target :: grd(rotation%ngrd)

      rotation%grd => grd

    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getrotationsource(rotation)
      implicit none
      class(t_rotation), intent(in) :: rotation
      real(8) :: getrotationsource(rotation%npv)
      real(8) :: corioli(3),centrifugal(3)
      getrotationsource = 0.d0

      corioli(1) = 2.d0*(rotation%omega(2)*rotation%pv(2,4)-rotation%omega(3)*rotation%pv(2,3))
      corioli(2) = 2.d0*(rotation%omega(3)*rotation%pv(2,2)-rotation%omega(1)*rotation%pv(2,4))
      corioli(3) = 2.d0*(rotation%omega(1)*rotation%pv(2,3)-rotation%omega(2)*rotation%pv(2,2))

      centrifugal(1) = rotation%omega(1)*(rotation%omega(2)*rotation%grd(3)+rotation%omega(3)*rotation%grd(4)) &
                     - rotation%grd(2)*(rotation%omega(2)**2+rotation%omega(3)**2)
      centrifugal(2) = rotation%omega(2)*(rotation%omega(1)*rotation%grd(2)+rotation%omega(3)*rotation%grd(4)) &
                     - rotation%grd(3)*(rotation%omega(1)**2+rotation%omega(3)**2)
      centrifugal(3) = rotation%omega(3)*(rotation%omega(1)*rotation%grd(2)+rotation%omega(2)*rotation%grd(3)) &
                     - rotation%grd(4)*(rotation%omega(1)**2+rotation%omega(2)**2)

      getrotationsource(2) = rotation%dv(1)*rotation%grd(1)*(corioli(1)+centrifugal(1))
      getrotationsource(3) = rotation%dv(1)*rotation%grd(1)*(corioli(2)+centrifugal(2))
      getrotationsource(4) = rotation%dv(1)*rotation%grd(1)*(corioli(3)+centrifugal(3))

    end function getrotationsource
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module rotation_module
