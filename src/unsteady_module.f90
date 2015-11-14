module unsteady_module
  use config_module
  use grid_module
  implicit none
  private
  public :: t_unsteady
 
  type t_unsteady
    private
    integer :: stencil,ngrd,npv,ndv
    real(8) :: pref,dt_phy
    real(8), pointer :: pv(:,:),dv(:),grd(:),qq1(:),qq2(:) 
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setgrd       ! (grd) vol,xcen,ycen,zcen,ydns
      procedure :: setpv        ! (pv) p,u,v,w,t,y1,y2,k,o
      procedure :: setdv        ! (dv) rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      procedure :: setqq        ! (qq1,qq2)
      procedure :: unsteadysource
  end type t_unsteady

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(unsteady,config,grid)
      implicit none
      class(t_unsteady), intent(out) :: unsteady
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
 
      unsteady%stencil = config%getstencil()
      unsteady%npv = config%getnpv()
      unsteady%ndv = config%getndv()
      unsteady%pref = config%getpref()
      unsteady%dt_phy = config%getdt_phy()

      unsteady%ngrd = grid%getngrd()

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
    subroutine setgrd(unsteady,grd)
      implicit none
      class(t_unsteady), intent(inout) :: unsteady
      real(8), intent(in), target :: grd(unsteady%ngrd)
      
      unsteady%grd => grd
      
    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(unsteady,pv)
      implicit none
      class(t_unsteady), intent(inout) :: unsteady
      real(8), intent(in), target :: pv(unsteady%stencil,unsteady%npv)
      
      unsteady%pv => pv
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(unsteady,dv)
      implicit none
      class(t_unsteady), intent(inout) :: unsteady
      real(8), intent(in), target :: dv(unsteady%ndv)
      
      unsteady%dv => dv
      
    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setqq(unsteady,qq1,qq2)
      implicit none
      class(t_unsteady), intent(inout) :: unsteady
      real(8), intent(in), target :: qq1(unsteady%npv),qq2(unsteady%npv)
      
      unsteady%qq1 => qq1
      unsteady%qq2 => qq2
      
    end subroutine setqq
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteadysource(unsteady) result(fx)
      implicit none
      class(t_unsteady), intent(in) :: unsteady
      real(8) :: fx(unsteady%npv)
      integer :: n
      
      fx(1) = unsteady%dv(1)
      fx(2) = unsteady%dv(1)*unsteady%pv(2,2)
      fx(3) = unsteady%dv(1)*unsteady%pv(2,3)
      fx(4) = unsteady%dv(1)*unsteady%pv(2,4)
      fx(5) = unsteady%dv(1)*(unsteady%dv(2)+0.5d0*(unsteady%pv(2,2)**2+unsteady%pv(2,3)**2+unsteady%pv(2,4)**2))-(unsteady%pv(2,1)+unsteady%pref)
      do n = 6,unsteady%npv
        fx(n) = unsteady%dv(1)*unsteady%pv(2,n)
      end do
      
      fx = (1.5d0*fx - 2.d0*unsteady%qq1 + 0.5d0*unsteady%qq2)*unsteady%grd(1)/unsteady%dt_phy       
      
    end function unsteadysource
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module unsteady_module
