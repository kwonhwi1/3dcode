module gravity_module
  use config_module
  use grid_module
  implicit none
  private
  public :: t_gravity

  type t_gravity
    private
    integer :: stencil,npv,ndv,ngrd,ndata
    real(8), dimension(:), allocatable :: tdata
    real(8), dimension(:,:), allocatable :: gravity
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
      integer :: m,io

      gravity%stencil = config%getstencil()
      gravity%npv = config%getnpv()
      gravity%ndv = config%getndv()

      if(config%getgravity().eq.4) then ! user-defined acceleration
        open(newunit=io,file='gravity.dat')
          read(io,*) gravity%ndata
          allocate(gravity%gravity(3,gravity%ndata), &
                   gravity%tdata(gravity%ndata))
          do m=1,gravity%ndata
            read(io,*) gravity%tdata(m),gravity%gravity(:,m)
          end do
        close(io)
      else ! constant gravity
        gravity%ndata = 1
        allocate(gravity%gravity(3,1),gravity%tdata(1))
        gravity%tdata=0.d0
        select case(config%getgravity())
        case(-1)
          gravity%gravity = reshape((/-9.80665d0,0.d0,0.d0/),(/3,1/))
        case(1)
          gravity%gravity = reshape((/9.80665d0,0.d0,0.d0/),(/3,1/))
        case(-2)
          gravity%gravity = reshape((/0.d0,-9.80665d0,0.d0/),(/3,1/))
        case(2)
          gravity%gravity = reshape((/0.d0,9.80665d0,0.d0/),(/3,1/))
        case(-3)
          gravity%gravity = reshape((/0.d0,0.d0,-9.80665d0/),(/3,1/))
        case(3)
          gravity%gravity = reshape((/0.d0,0.d0,9.80665d0/),(/3,1/))
        case default
          gravity%gravity = reshape((/0.d0,0.d0,0.d0/),(/3,1/))
        end select
      end if

      gravity%ngrd = grid%getngrd()

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(gravity)
      implicit none
      class(t_gravity), intent(inout) :: gravity

      if(allocated(gravity%tdata))   deallocate(gravity%tdata)
      if(allocated(gravity%gravity)) deallocate(gravity%gravity)
      if(associated(gravity%pv))     nullify(gravity%pv)
      if(associated(gravity%dv))     nullify(gravity%dv)
      if(associated(gravity%grd))    nullify(gravity%grd)

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
    function getgravitysource(gravity,time)
      implicit none
      class(t_gravity), intent(in) :: gravity
      real(8), intent(in) :: time
      real(8) :: getgravitysource(gravity%npv)
      real(8) :: accel(3)
      integer :: i

      do i=1,gravity%ndata-1
        if((time.ge.gravity%tdata(i)).and.(time.lt.gravity%tdata(i+1))) then
          accel(:) = gravity%gravity(:,i)+(gravity%gravity(:,i+1)-gravity%gravity(:,i))* &
                    (time-gravity%tdata(i))/(gravity%tdata(i+1)-gravity%tdata(i))
          exit
        end if
      end do
      if(time.ge.gravity%tdata(gravity%ndata)) accel(:) = gravity%gravity(:,gravity%ndata)

      getgravitysource = 0.d0
      getgravitysource(2:4) = gravity%dv(1)*gravity%grd(1)*accel(1:3)
      getgravitysource(5) = gravity%dv(1)*gravity%grd(1)*(accel(1)*gravity%pv(2,2) &
                                                        + accel(2)*gravity%pv(2,3) &
                                                        + accel(3)*gravity%pv(2,4))

    end function getgravitysource
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module gravity_module
