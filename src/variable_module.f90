module variable_module
  use config_module
  use grid_module
  implicit none
  private
  public :: t_variable
  
  type t_variable
    private
    integer :: npv,ndv,ntv,nqq
    real(8), dimension(:,:,:,:), allocatable :: pv   ! p,u,v,w,t,y1,y2,k,o
    real(8), dimension(:,:,:,:), allocatable :: dv   ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
    real(8), dimension(:,:,:,:), allocatable :: tv   ! vis,cond,emut
    real(8), dimension(:,:,:,:,:), allocatable :: qq ! unsteady n-1,n-2 conservative variable
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: getnpv
      procedure :: getndv
      procedure :: getntv
      procedure :: getnqq
      procedure :: getpv
      procedure :: getdv
      procedure :: gettv
      procedure :: getqq
      procedure :: setpv
      procedure :: setdv
      procedure :: settv
      procedure :: setqq
  end type t_variable
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(variable,config,grid)
      implicit none
      class(t_variable), intent(out) :: variable
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
          
      select case(config%getiturb())
      case(-3)
        variable%npv = 7
        variable%ntv = 0  
      case(-2)
        variable%npv = 7
        variable%ntv = 2      
      case(-1,0)
        variable%npv = 9
        variable%ntv = 3      
      end select
      variable%ndv = 18
      select case(config%getnsteady())
      case(0)
        variable%nqq = 0
      case(1)
        variable%nqq = 2
      end select      
      
      allocate(variable%pv(variable%npv,-1:grid%getimax()+3,-1:grid%getjmax()+3,-1:grid%getkmax()+3))
      allocate(variable%dv(variable%ndv,-1:grid%getimax()+3,-1:grid%getjmax()+3,-1:grid%getkmax()+3))
      allocate(variable%tv(variable%ntv,-1:grid%getimax()+3,-1:grid%getjmax()+3,-1:grid%getkmax()+3))
      allocate(variable%qq(variable%npv,variable%nqq,2:grid%getimax(),2:grid%getjmax(),2:grid%getkmax()))

      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(variable)
      implicit none
      class(t_variable), intent(inout) :: variable
          
      if(allocated(variable%pv)) deallocate(variable%pv)
      if(allocated(variable%dv)) deallocate(variable%dv)
      if(allocated(variable%tv)) deallocate(variable%tv)
      if(allocated(variable%qq)) deallocate(variable%qq)
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getnpv(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getnpv
      
      getnpv = variable%npv
      
    end function getnpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getndv(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getndv
      
      getndv = variable%ndv
      
    end function getndv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getntv(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getntv
      
      getntv = variable%ntv
      
    end function getntv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getnqq(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getnqq
      
      getnqq = variable%nqq
      
    end function getnqq
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getpv(variable,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k
      real(8) :: getpv(variable%npv)
      
      getpv = variable%pv(:,i,j,k)
      
    end function getpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getdv(variable,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k
      real(8) :: getdv(variable%ndv)
      
      getdv = variable%dv(:,i,j,k)
      
    end function getdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function gettv(variable,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k
      real(8) :: gettv(variable%ntv)
      
      gettv = variable%tv(:,i,j,k)
      
    end function gettv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getqq(variable,n,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k,n
      real(8) :: getqq(variable%npv)
      
      getqq = variable%qq(:,n,i,j,k)
      
    end function getqq
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var
      
      variable%pv(n,i,j,k) = var
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var
      
      variable%dv(n,i,j,k) = var
      
    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine settv(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var
      
      variable%tv(n,i,j,k) = var
      
    end subroutine settv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setqq(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var(variable%npv)
      
      variable%qq(:,n,i,j,k) = var
      
    end subroutine setqq
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module variable_module
