module lhs_module
  use config_module
  use grid_module
  implicit none
  private
  public :: t_lhs,t_lhs_flowonly,t_lhs_flowturball,t_lhs_flowonly_ex,t_lhs_flowturball_ex
  
  type, abstract :: t_lhs
    private
    integer :: npv,ndv,ntv,ngrd,nsteady
    real(8) :: uref,str,dt_phy
    real(8) :: omega(3)
    real(8), pointer :: cx1(:),cx2(:),ex1(:),ex2(:),tx1(:),tx2(:)
    real(8), pointer :: pv(:),tv(:),dv(:),grd(:),dt
    real(8), pointer :: c(:),t(:)
    real(8), pointer :: c0(:),t0(:)
    procedure(p_getsndp2), pointer :: getsndp2
    procedure(p_geteigenvis), pointer :: geteigenvis
    procedure(p_getenthalpy), pointer :: getenthalpy
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setnorm ! (cx1,cx2,ex1,ex2,tx1,tx2)
      procedure :: setgrd  ! vol,xcen,ycen,zcen,ydns
      procedure :: setpv   ! p,u,v,w,t,y1,y2,k,o
      procedure :: setdv   ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      procedure :: settv   ! vis,cond,emut
      procedure :: setdt   ! dt
      procedure :: setc
      procedure :: sett
      procedure(p_getx), deferred :: getx
      procedure, private :: eigen
  end type t_lhs
  
  type, extends(t_lhs) :: t_lhs_flowonly
    contains
      procedure :: getx => flowonly
  end type t_lhs_flowonly
  
  type, extends(t_lhs) :: t_lhs_flowturball
    contains
      procedure :: getx => flowturball
  end type t_lhs_flowturball

  type, extends(t_lhs) :: t_lhs_flowonly_ex
    contains
      procedure :: getx => flowonly_ex
  end type t_lhs_flowonly_ex
  
  type, extends(t_lhs) :: t_lhs_flowturball_ex
    contains
      procedure :: getx => flowturball_ex
  end type t_lhs_flowturball_ex

  abstract interface
    function p_getx(lhs,res) result(x)
      import t_lhs
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8), intent(in) :: res(lhs%npv)
      real(8) :: x(lhs%npv)
    end function p_getx
  end interface

  interface
    function p_getsndp2(lhs,uuu2) result(sndp2)
      import t_lhs
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8), intent(in) :: uuu2
      real(8) :: sndp2
    end function p_getsndp2
    
    function p_geteigenvis(lhs) result(eigenvis)
      import t_lhs
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: eigenvis
    end function p_geteigenvis

    function p_getenthalpy(lhs)
      import t_lhs
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: p_getenthalpy
    end function p_getenthalpy
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(lhs,config,grid)
      implicit none
      class(t_lhs), intent(out) :: lhs
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      
      lhs%ngrd = grid%getngrd()
      lhs%npv = config%getnpv()
      lhs%ndv = config%getndv()
      lhs%ntv = config%getntv()

      lhs%uref = config%geturef()
      lhs%str  = config%getstr()
      lhs%dt_phy = config%getdt_phy()
      lhs%nsteady = config%getnsteady()
      lhs%omega = config%getomega()
      
      select case(config%getprec())
      case(0)
        lhs%getsndp2 => no_prec
      case(1)
        lhs%getsndp2 => steady_prec
      case(2)
        lhs%getsndp2 => unsteady_prec
      end select

      select case(config%getiturb())
      case(-1,0)
        lhs%geteigenvis => turbulent
      case(-2)
        lhs%geteigenvis => laminar
      case(-3)
        lhs%geteigenvis => euler
      end select

      select case(config%getrotation())
      case(0)
        lhs%getenthalpy => enthalpy
      case(-1,1,-2,2,-3,3)
        lhs%getenthalpy => rothalpy
      case default
      end select

      allocate(lhs%c0(4))
      allocate(lhs%t0(4))
      lhs%c0 = 0.d0
      lhs%t0 = 0.d0
      lhs%c => lhs%c0
      lhs%t => lhs%t0
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(lhs)
      implicit none
      class(t_lhs), intent(inout) :: lhs

      if(associated(lhs%cx1))         nullify(lhs%cx1)          
      if(associated(lhs%cx2))         nullify(lhs%cx2)          
      if(associated(lhs%ex1))         nullify(lhs%ex1)          
      if(associated(lhs%ex2))         nullify(lhs%ex2)
      if(associated(lhs%tx1))         nullify(lhs%tx1)          
      if(associated(lhs%tx2))         nullify(lhs%tx2)  
      if(associated(lhs%grd))         nullify(lhs%grd)          
      if(associated(lhs%pv))          nullify(lhs%pv)           
      if(associated(lhs%dv))          nullify(lhs%dv)           
      if(associated(lhs%tv))          nullify(lhs%tv)           
      if(associated(lhs%dt))          nullify(lhs%dt)           
      if(associated(lhs%c))           nullify(lhs%c)            
      if(associated(lhs%t))           nullify(lhs%t)            
      if(associated(lhs%getsndp2))    nullify(lhs%getsndp2)     
      if(associated(lhs%geteigenvis)) nullify(lhs%geteigenvis)
      if(associated(lhs%getenthalpy)) nullify(lhs%getenthalpy)
      deallocate(lhs%c0,lhs%t0)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setnorm(lhs,cx1,cx2,ex1,ex2,tx1,tx2)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: cx1(3),cx2(3),ex1(3),ex2(3),tx1(3),tx2(3)
      
      lhs%cx1 => cx1
      lhs%cx2 => cx2
      lhs%ex1 => ex1
      lhs%ex2 => ex2
      lhs%tx1 => tx1
      lhs%tx2 => tx2
      
    end subroutine setnorm
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(lhs,grd)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: grd(lhs%ngrd)
      
      lhs%grd => grd
      
    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(lhs,pv)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: pv(lhs%npv)
      
      lhs%pv => pv
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(lhs,dv)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: dv(lhs%ndv)
      
      lhs%dv => dv
      
    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine settv(lhs,tv)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: tv(lhs%ntv)
      
      lhs%tv => tv
      
    end subroutine settv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdt(lhs,dt)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: dt
      
      lhs%dt => dt
      
    end subroutine setdt
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setc(lhs,c)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: c(4)
      
      lhs%c => c
      
    end subroutine setc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine sett(lhs,t)
      implicit none
      class(t_lhs), intent(inout) :: lhs
      real(8), intent(in), target :: t(4)
      
      lhs%t => t
      
    end subroutine sett
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function flowonly_ex(lhs,res) result(x)
      implicit none
      class(t_lhs_flowonly_ex), intent(in) :: lhs
      real(8), intent(in) :: res(lhs%npv)
      real(8) :: x(lhs%npv)
      real(8) :: a,b
      real(8) :: b1,i_b1,uv2,ss,rr,hst,cross,cross1,cross2,dhdp1
      real(8) :: ssy1,ssy2,ssh,aaa1,aaa2,uuu,coeff1,coeff2,coeff3,mm

      a = lhs%dt/lhs%grd(1)
      b = dble(lhs%nsteady)*1.5d0*lhs%grd(1)/lhs%dt_phy*a
      
      b1 = 1.d0+b
      i_b1 = 1.d0/(b1*lhs%dv(1))
      uv2 = lhs%pv(2)**2+lhs%pv(3)**2+lhs%pv(4)**2
      ss = lhs%dv(1)*lhs%dv(12)*(1.d0/lhs%getsndp2(uv2)+b/lhs%dv(6))
      rr = 1.d0/lhs%getsndp2(uv2)-1.d0/lhs%dv(6)+b1*lhs%dv(7)
      hst = lhs%getenthalpy()-uv2
      cross = lhs%dv(14)*lhs%dv(9)-lhs%dv(13)*lhs%dv(10)
      cross1 = lhs%dv(13)*lhs%dv(8)-lhs%dv(12)*lhs%dv(9)
      cross2 = lhs%dv(14)*lhs%dv(8)-lhs%dv(12)*lhs%dv(10)
      dhdp1 = 1.d0 - lhs%dv(11)*lhs%dv(1)
      ssy1 = lhs%dv(13)*lhs%dv(1)*rr+b1*lhs%dv(9)*dhdp1
      ssy2 = lhs%dv(14)*lhs%dv(1)*rr+b1*lhs%dv(10)*dhdp1
      ssh = lhs%dv(1)*rr*hst-b1*dhdp1*lhs%dv(1)
      aaa1 = lhs%dv(9)*hst+lhs%dv(13)*lhs%dv(1)-cross*lhs%pv(7)
      aaa2 = lhs%dv(8)*hst+lhs%dv(12)*lhs%dv(1)-cross2*lhs%pv(7)


      uuu = res(2)*lhs%pv(2)+res(3)*lhs%pv(3)+res(4)*lhs%pv(4)-res(5)

      coeff1 = res(1)*aaa1 + lhs%dv(9)*uuu + cross*res(7)
      coeff2 = res(1)*aaa2 + lhs%dv(8)*uuu + cross2*res(7)
      coeff3 = res(1)*(ssh-ssy2*lhs%pv(7)) + lhs%dv(1)*rr*uuu + res(7)*ssy2

      mm = 1.d0/ss

      x(1) = (coeff2+cross1*res(6)-cross1*lhs%pv(6)*res(1))*mm
      x(2) = (res(2)-res(1)*lhs%pv(2))*i_b1
      x(3) = (res(3)-res(1)*lhs%pv(3))*i_b1
      x(4) = (res(4)-res(1)*lhs%pv(4))*i_b1
      x(5) = -(coeff3-ssy1*lhs%pv(6)*res(1)+res(6)*ssy1)*mm*i_b1
      x(6) = (res(6)-res(1)*lhs%pv(6))*i_b1
      x(7) = (res(7)-res(1)*lhs%pv(7))*i_b1

      x = x*a
    end function flowonly_ex
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function flowturball_ex(lhs,res) result(x)
      implicit none
      class(t_lhs_flowturball_ex), intent(in) :: lhs
      real(8), intent(in) :: res(lhs%npv)
      real(8) :: x(lhs%npv)
      real(8) :: a,b
      real(8) :: b1,i_b1,uv2,ss,rr,hst,cross,cross1,cross2,dhdp1
      real(8) :: ssy1,ssy2,ssh,aaa1,aaa2,uuu,coeff1,coeff2,coeff3,mm

      a = lhs%dt/lhs%grd(1)
      b = dble(lhs%nsteady)*1.5d0*lhs%grd(1)/lhs%dt_phy*a
      
      b1 = 1.d0+b
      i_b1 = 1.d0/(b1*lhs%dv(1))
      uv2 = lhs%pv(2)**2+lhs%pv(3)**2+lhs%pv(4)**2
      ss = lhs%dv(1)*lhs%dv(12)*(1.d0/lhs%getsndp2(uv2)+b/lhs%dv(6))
      rr = 1.d0/lhs%getsndp2(uv2)-1.d0/lhs%dv(6)+b1*lhs%dv(7)
      hst = lhs%getenthalpy()-uv2
      cross = lhs%dv(14)*lhs%dv(9)-lhs%dv(13)*lhs%dv(10)
      cross1 = lhs%dv(13)*lhs%dv(8)-lhs%dv(12)*lhs%dv(9)
      cross2 = lhs%dv(14)*lhs%dv(8)-lhs%dv(12)*lhs%dv(10)
      dhdp1 = 1.d0 - lhs%dv(11)*lhs%dv(1)
      ssy1 = lhs%dv(13)*lhs%dv(1)*rr+b1*lhs%dv(9)*dhdp1
      ssy2 = lhs%dv(14)*lhs%dv(1)*rr+b1*lhs%dv(10)*dhdp1
      ssh = lhs%dv(1)*rr*hst-b1*dhdp1*lhs%dv(1)
      aaa1 = lhs%dv(9)*hst+lhs%dv(13)*lhs%dv(1)-cross*lhs%pv(7)
      aaa2 = lhs%dv(8)*hst+lhs%dv(12)*lhs%dv(1)-cross2*lhs%pv(7)


      uuu = res(2)*lhs%pv(2)+res(3)*lhs%pv(3)+res(4)*lhs%pv(4)-res(5)

      coeff1 = res(1)*aaa1 + lhs%dv(9)*uuu + cross*res(7)
      coeff2 = res(1)*aaa2 + lhs%dv(8)*uuu + cross2*res(7)
      coeff3 = res(1)*(ssh-ssy2*lhs%pv(7)) + lhs%dv(1)*rr*uuu + res(7)*ssy2

      mm = 1.d0/ss

      x(1) = (coeff2+cross1*res(6)-cross1*lhs%pv(6)*res(1))*mm
      x(2) = (res(2)-res(1)*lhs%pv(2))*i_b1
      x(3) = (res(3)-res(1)*lhs%pv(3))*i_b1
      x(4) = (res(4)-res(1)*lhs%pv(4))*i_b1
      x(5) = -(coeff3-ssy1*lhs%pv(6)*res(1)+res(6)*ssy1)*mm*i_b1
      x(6) = (res(6)-res(1)*lhs%pv(6))*i_b1
      x(7) = (res(7)-res(1)*lhs%pv(7))*i_b1
      x(8) = (res(8)-res(1)*lhs%pv(8))*i_b1
      x(9) = (res(9)-res(1)*lhs%pv(9))*i_b1

      x = x*a

    end function flowturball_ex
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function flowonly(lhs,res) result(x)
      implicit none
      class(t_lhs_flowonly), intent(in) :: lhs
      real(8), intent(in) :: res(lhs%npv)
      real(8) :: x(lhs%npv)
      real(8) :: a,b,c(4)
      real(8) :: b1,i_b1,uv2,ss,rr,hst,cross,cross1,cross2,dhdp1
      real(8) :: ssy1,ssy2,ssh,aaa1,aaa2,uuu,coeff1,coeff2,coeff3,mm

      a = 1.d0/(lhs%grd(1)/lhs%dt + lhs%eigen())
      b = (dble(lhs%nsteady)*1.5d0*lhs%grd(1)/lhs%dt_phy + 2.d0*lhs%geteigenvis()/lhs%grd(1))*a      
      c = lhs%c*lhs%grd(1)*a
      
      b1 = 1.d0+b
      i_b1 = 1.d0/(b1*lhs%dv(1))
      uv2 = lhs%pv(2)**2+lhs%pv(3)**2+lhs%pv(4)**2
      ss = lhs%dv(1)*lhs%dv(12)*(1.d0/lhs%getsndp2(uv2)+b/lhs%dv(6))
      rr = 1.d0/lhs%getsndp2(uv2)-1.d0/lhs%dv(6)+b1*lhs%dv(7)
      hst = lhs%getenthalpy()-uv2
      cross = lhs%dv(14)*lhs%dv(9)-lhs%dv(13)*lhs%dv(10)
      cross1 = lhs%dv(13)*lhs%dv(8)-lhs%dv(12)*lhs%dv(9)
      cross2 = lhs%dv(14)*lhs%dv(8)-lhs%dv(12)*lhs%dv(10)
      dhdp1 = 1.d0 - lhs%dv(11)*lhs%dv(1)
      ssy1 = lhs%dv(13)*lhs%dv(1)*rr+b1*lhs%dv(9)*dhdp1
      ssy2 = lhs%dv(14)*lhs%dv(1)*rr+b1*lhs%dv(10)*dhdp1
      ssh = lhs%dv(1)*rr*hst-b1*dhdp1*lhs%dv(1)
      aaa1 = lhs%dv(9)*hst+lhs%dv(13)*lhs%dv(1)-cross*lhs%pv(7)
      aaa2 = lhs%dv(8)*hst+lhs%dv(12)*lhs%dv(1)-cross2*lhs%pv(7)


      uuu = res(2)*lhs%pv(2)+res(3)*lhs%pv(3)+res(4)*lhs%pv(4)-res(5)

      coeff1 = res(1)*aaa1 + lhs%dv(9)*uuu + cross*res(7)
      coeff2 = res(1)*aaa2 + lhs%dv(8)*uuu + cross2*res(7)
      coeff3 = res(1)*(ssh-ssy2*lhs%pv(7)) + lhs%dv(1)*rr*uuu + res(7)*ssy2

      mm = b1*lhs%dv(1)*(ss + c(1)*cross1) - c(2)*ssy1 + c(3)*ss
      mm = 1.d0/mm

      x(1) = ( c(4)*cross1*(res(1)*lhs%pv(7)-res(7)) - c(2)*coeff1 + c(3)*coeff2 &
           + lhs%dv(1)*b1*(coeff2+cross1*res(6)-cross1*lhs%pv(6)*res(1)) )*mm
      x(2) = (res(2)-res(1)*lhs%pv(2))*i_b1
      x(3) = (res(3)-res(1)*lhs%pv(3))*i_b1
      x(4) = (res(4)-res(1)*lhs%pv(4))*i_b1
      x(5) = (c(4)*i_b1*ssy1*(res(7)-res(1)*lhs%pv(7)) + c(1)*coeff1 -c(3)*i_b1*coeff3 &
              -(coeff3-ssy1*lhs%pv(6)*res(1)+res(6)*ssy1) )*mm
      x(6) = ( ss*(res(6)-res(1)*lhs%pv(6)) -c(4)*ss*i_b1*(res(7)-res(1)*lhs%pv(7)) -c(1)*coeff2 &
             + c(2)*i_b1*coeff3 )*mm
      x(7) = (res(7)-res(1)*lhs%pv(7))*i_b1

      x = x*a

    end function flowonly
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function flowturball(lhs,res) result(x)
      implicit none
      class(t_lhs_flowturball), intent(in) :: lhs
      real(8), intent(in) :: res(lhs%npv)
      real(8) :: x(lhs%npv)
      real(8) :: a,b,c(4),t(4)
      real(8) :: b1,i_b1,uv2,ss,rr,hst,cross,cross1,cross2,dhdp1
      real(8) :: ssy1,ssy2,ssh,aaa1,aaa2,uuu,coeff1,coeff2,coeff3,mm,tt

      a = 1.d0/(lhs%grd(1)/lhs%dt + lhs%eigen())
      b = (dble(lhs%nsteady)*1.5d0*lhs%grd(1)/lhs%dt_phy + 2.d0*lhs%geteigenvis()/lhs%grd(1))*a      
      c = lhs%c*lhs%grd(1)*a
      t = lhs%t*lhs%grd(1)*a
      
      b1 = 1.d0+b
      i_b1 = 1.d0/(b1*lhs%dv(1))
      uv2 = lhs%pv(2)**2+lhs%pv(3)**2+lhs%pv(4)**2
      ss = lhs%dv(1)*lhs%dv(12)*(1.d0/lhs%getsndp2(uv2)+b/lhs%dv(6))
      rr = 1.d0/lhs%getsndp2(uv2)-1.d0/lhs%dv(6)+b1*lhs%dv(7)
      hst = lhs%getenthalpy()-uv2
      cross = lhs%dv(14)*lhs%dv(9)-lhs%dv(13)*lhs%dv(10)
      cross1 = lhs%dv(13)*lhs%dv(8)-lhs%dv(12)*lhs%dv(9)
      cross2 = lhs%dv(14)*lhs%dv(8)-lhs%dv(12)*lhs%dv(10)
      dhdp1 = 1.d0 - lhs%dv(11)*lhs%dv(1)
      ssy1 = lhs%dv(13)*lhs%dv(1)*rr+b1*lhs%dv(9)*dhdp1
      ssy2 = lhs%dv(14)*lhs%dv(1)*rr+b1*lhs%dv(10)*dhdp1
      ssh = lhs%dv(1)*rr*hst-b1*dhdp1*lhs%dv(1)
      aaa1 = lhs%dv(9)*hst+lhs%dv(13)*lhs%dv(1)-cross*lhs%pv(7)
      aaa2 = lhs%dv(8)*hst+lhs%dv(12)*lhs%dv(1)-cross2*lhs%pv(7)


      uuu = res(2)*lhs%pv(2)+res(3)*lhs%pv(3)+res(4)*lhs%pv(4)-res(5)

      coeff1 = res(1)*aaa1 + lhs%dv(9)*uuu + cross*res(7)
      coeff2 = res(1)*aaa2 + lhs%dv(8)*uuu + cross2*res(7)
      coeff3 = res(1)*(ssh-ssy2*lhs%pv(7)) + lhs%dv(1)*rr*uuu + res(7)*ssy2

      mm = b1*lhs%dv(1)*(ss + c(1)*cross1) - c(2)*ssy1 + c(3)*ss
      mm = 1.d0/mm

      x(1) = ( c(4)*cross1*(res(1)*lhs%pv(7)-res(7)) - c(2)*coeff1 + c(3)*coeff2 &
           + lhs%dv(1)*b1*(coeff2+cross1*res(6)-cross1*lhs%pv(6)*res(1)) )*mm
      x(2) = (res(2)-res(1)*lhs%pv(2))*i_b1
      x(3) = (res(3)-res(1)*lhs%pv(3))*i_b1
      x(4) = (res(4)-res(1)*lhs%pv(4))*i_b1
      x(5) = (c(4)*i_b1*ssy1*(res(7)-res(1)*lhs%pv(7)) + c(1)*coeff1 -c(3)*i_b1*coeff3 &
              -(coeff3-ssy1*lhs%pv(6)*res(1)+res(6)*ssy1) )*mm
      x(6) = ( ss*(res(6)-res(1)*lhs%pv(6)) -c(4)*ss*i_b1*(res(7)-res(1)*lhs%pv(7)) -c(1)*coeff2 &
             + c(2)*i_b1*coeff3 )*mm
      x(7) = (res(7)-res(1)*lhs%pv(7))*i_b1

      tt = 1.d0/( t(2)*t(3)-(t(1)+b1*lhs%dv(1))*(t(4)+b1*lhs%dv(1)) )
      x(8) = ((res(9) - res(1)*lhs%pv(9))*t(2) - (b1*lhs%dv(1)+t(4))*(res(8)-res(1)*lhs%pv(8)))*tt
      x(9) = ((res(8) - res(1)*lhs%pv(8))*t(3) - (b1*lhs%dv(1)+t(1))*(res(9)-res(1)*lhs%pv(9)))*tt

      x = x*a

    end function flowturball
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function eigen(lhs) result(lamda)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: lamda
      real(8) :: uv2,sndp2,ox,oy,oz,ds,u,up,d
      real(8), parameter :: cappa = 1.05d0
      
      uv2 = lhs%pv(2)**2+lhs%pv(3)**2+lhs%pv(4)**2
      sndp2 = lhs%getsndp2(uv2)
      
      ox = 0.5d0*(lhs%cx1(1)+lhs%cx2(1))
      oy = 0.5d0*(lhs%cx1(2)+lhs%cx2(2))
      oz = 0.5d0*(lhs%cx1(3)+lhs%cx2(3))
      ds = ox**2+oy**2+oz**2
      u  = ox*lhs%pv(2)+oy*lhs%pv(3)+oz*lhs%pv(4)
      up = 0.5d0*(1.d0+sndp2/lhs%dv(6))*u
      d  = 0.5d0*dsqrt(u**2*(1.d0-sndp2/lhs%dv(6))**2+4.d0*sndp2*ds)
      lamda = cappa*(dabs(up)+d)

      ox = 0.5d0*(lhs%ex1(1)+lhs%ex2(1))
      oy = 0.5d0*(lhs%ex1(2)+lhs%ex2(2))
      oz = 0.5d0*(lhs%ex1(3)+lhs%ex2(3))
      ds = ox**2+oy**2+oz**2
      u  = ox*lhs%pv(2)+oy*lhs%pv(3)+oz*lhs%pv(4)
      up = 0.5d0*(1.d0+sndp2/lhs%dv(6))*u
      d  = 0.5d0*dsqrt(u**2*(1.d0-sndp2/lhs%dv(6))**2+4.d0*sndp2*ds)
      lamda = lamda+cappa*(dabs(up)+d)
      
      ox = 0.5d0*(lhs%tx1(1)+lhs%tx2(1))
      oy = 0.5d0*(lhs%tx1(2)+lhs%tx2(2))
      oz = 0.5d0*(lhs%tx1(3)+lhs%tx2(3))
      ds = ox**2+oy**2+oz**2
      u  = ox*lhs%pv(2)+oy*lhs%pv(3)+oz*lhs%pv(4)
      up = 0.5d0*(1.d0+sndp2/lhs%dv(6))*u
      d  = 0.5d0*dsqrt(u**2*(1.d0-sndp2/lhs%dv(6))**2+4.d0*sndp2*ds)
      lamda = lamda+cappa*(dabs(up)+d)
      
    end function eigen
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function euler(lhs) result(eigenvis)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: eigenvis
      
      eigenvis = 0.d0
      
    end function euler
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function laminar(lhs) result(eigenvis)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: eigenvis

      eigenvis = dmax1(lhs%dv(7)*lhs%tv(2)/(lhs%dv(8)     &
                       +lhs%dv(12)*lhs%dv(7)*lhs%dv(1)    &
                       -lhs%dv(11)*lhs%dv(8)*lhs%dv(1)),  &
                       4.d0/3.d0*lhs%tv(1)/lhs%dv(1))
      eigenvis = eigenvis*( 0.25d0*(lhs%cx1(1)+lhs%cx2(1))**2 + 0.25d0*(lhs%cx1(2)+lhs%cx2(2))**2 + 0.25d0*(lhs%cx1(3)+lhs%cx2(3))**2 &
                          + 0.25d0*(lhs%ex1(1)+lhs%ex2(1))**2 + 0.25d0*(lhs%ex1(2)+lhs%ex2(2))**2 + 0.25d0*(lhs%ex1(3)+lhs%ex2(3))**2 &
                          + 0.25d0*(lhs%tx1(1)+lhs%tx2(1))**2 + 0.25d0*(lhs%tx1(2)+lhs%tx2(2))**2 + 0.25d0*(lhs%tx1(3)+lhs%tx2(3))**2 )
    end function laminar
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function turbulent(lhs) result(eigenvis)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: eigenvis
      real(8), parameter :: pr_t = 0.9d0
   
      eigenvis = dmax1(lhs%dv(7)*(lhs%tv(2)+lhs%dv(12)*lhs%tv(3)/pr_t) &
                      /(lhs%dv(8)+lhs%dv(12)*lhs%dv(7)*lhs%dv(1)       &
                       -lhs%dv(11)*lhs%dv(8)*lhs%dv(1)),               &
                       4.d0/3.d0*(lhs%tv(1)+lhs%tv(3))/lhs%dv(1))
      eigenvis = eigenvis*( 0.25d0*(lhs%cx1(1)+lhs%cx2(1))**2 + 0.25d0*(lhs%cx1(2)+lhs%cx2(2))**2 + 0.25d0*(lhs%cx1(3)+lhs%cx2(3))**2 &
                          + 0.25d0*(lhs%ex1(1)+lhs%ex2(1))**2 + 0.25d0*(lhs%ex1(2)+lhs%ex2(2))**2 + 0.25d0*(lhs%ex1(3)+lhs%ex2(3))**2 &
                          + 0.25d0*(lhs%tx1(1)+lhs%tx2(1))**2 + 0.25d0*(lhs%tx1(2)+lhs%tx2(2))**2 + 0.25d0*(lhs%tx1(3)+lhs%tx2(3))**2 )
    end function turbulent
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function no_prec(lhs,uuu2) result(sndp2)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8), intent(in) :: uuu2
      real(8) :: sndp2
      
      sndp2 = lhs%dv(6)
      
    end function no_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function steady_prec(lhs,uuu2) result(sndp2)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8), intent(in) :: uuu2
      real(8) :: sndp2  
      
      sndp2 = dmin1(lhs%dv(6),dmax1(uuu2,0.09d0*lhs%uref**2))
      
    end function steady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteady_prec(lhs,uuu2) result(sndp2)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8), intent(in) :: uuu2
      real(8) :: sndp2  
      
      sndp2 = dmin1(lhs%dv(6),dmax1(uuu2,0.09d0*lhs%uref**2,lhs%str**2))
      
    end function unsteady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function enthalpy(lhs)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: enthalpy
      enthalpy = lhs%dv(2) + 0.5d0*(lhs%pv(2)**2 + lhs%pv(3)**2 + lhs%pv(4)**2)
    end function enthalpy
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rothalpy(lhs)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: rothalpy
      rothalpy = lhs%dv(2) + 0.5d0*(lhs%pv(2)**2 + lhs%pv(3)**2 + lhs%pv(4)**2)
      rothalpy = rothalpy -0.5d0*(lhs%omega(1)**2*(lhs%grd(3)**2+lhs%grd(4)**2) &
                                + lhs%omega(2)**2*(lhs%grd(2)**2+lhs%grd(4)**2) &
                                + lhs%omega(3)**2*(lhs%grd(2)**2+lhs%grd(3)**2) &
                           -2.d0*(lhs%omega(1)*lhs%omega(2)*lhs%grd(2)*lhs%grd(3) &
                                + lhs%omega(2)*lhs%omega(3)*lhs%grd(3)*lhs%grd(4) &
                                + lhs%omega(3)*lhs%omega(1)*lhs%grd(4)*lhs%grd(2)))
    end function rothalpy
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module lhs_module
