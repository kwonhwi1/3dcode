module unsteady_module
  use config_module
  use variable_module
  implicit none
  private
  public :: t_unsteady
 
  type t_unsteady
    private
    integer :: npv
    real(8) :: pref,dt_phy
    real(8), pointer, public :: pv(:),dv(:),grd(:),qq1(:),qq2(:) 
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: unsteadysource
  end type t_unsteady

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(unsteady,config,variable)
      implicit none
      class(t_unsteady), intent(out) :: unsteady
      type(t_config), intent(in) :: config
      type(t_variable), intent(in) :: variable
 
      unsteady%pref = config%getpref()
      unsteady%dt_phy = config%getdt_phy()

      unsteady%npv = variable%getnpv()

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(unsteady)
      implicit none
      class(t_unsteady), intent(inout) :: unsteady
      
      if(associated(unsteady%grd)) nullify(unsteady%grd)      
      if(associated(unsteady%pv))  nullify(unsteady%pv) 
      if(associated(unsteady%dv))  nullify(unsteady%dv) 
      if(associated(unsteady%qq1)) nullify(unsteady%qq1)
      if(associated(unsteady%qq2)) nullify(unsteady%qq2)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteadysource(unsteady) result(fx)
      implicit none
      class(t_unsteady), intent(in) :: unsteady
      real(8) :: fx(unsteady%npv)
      integer :: n
      
      fx(1) = unsteady%dv(1)
      fx(2) = unsteady%dv(1)*unsteady%pv(2)
      fx(3) = unsteady%dv(1)*unsteady%pv(3)
      fx(4) = unsteady%dv(1)*unsteady%pv(4)
      fx(5) = unsteady%dv(1)*(unsteady%dv(2)+0.5d0*(unsteady%pv(2)**2+unsteady%pv(3)**2+unsteady%pv(4)**2))-(unsteady%pv(1)+unsteady%pref)
      do n = 6,unsteady%npv
        fx(n) = unsteady%dv(1)*unsteady%pv(n)
      end do
      
      fx = (1.5d0*fx - 2.d0*unsteady%qq1 + 0.5d0*unsteady%qq2)*unsteady%grd(1)/unsteady%dt_phy       
      
    end function unsteadysource
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module unsteady_module
