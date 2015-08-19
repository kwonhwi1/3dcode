module cav_module
  use config_module
  use eos_module
  implicit none
  private
  public :: t_cav,t_merkle,t_kunz,t_singhal,t_cav_result
  
  type t_cav_result
     real(8) :: cavsource,icav(4)
  end type
  
  type, abstract :: t_cav
    private
    real(8) :: pref
    real(8) :: uref,c_c,c_v,dp_ref,t_ref
    real(8), pointer, public :: pv(:),dv(:),grd(:)
    contains
      procedure :: construct
      procedure :: destruct
      procedure(p_cavsource), deferred :: cavsource
  end type t_cav
  
  type, extends(t_cav) :: t_merkle
    contains
      procedure :: cavsource => merkle
  end type t_merkle
  
  type, extends(t_cav) :: t_kunz
    contains
      procedure :: cavsource => kunz
  end type t_kunz
  
  type, extends(t_cav) :: t_singhal
    contains
      procedure :: cavsource => singhal
  end type t_singhal
  
  abstract interface
    function p_cavsource(cav,eos) result(cav_result)
      import t_cav_result
      import t_cav
      import t_eos
      implicit none
      class(t_cav), intent(in) :: cav
      type(t_eos),intent(in) :: eos
      type(t_cav_result) :: cav_result
    end function p_cavsource
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(cav,config)
      implicit none
      class(t_cav), intent(out) :: cav
      type(t_config), intent(in) :: config
        
      cav%pref         = config%getpref()
      cav%uref         = config%geturef()
      cav%c_c          = config%getc_c()
      cav%c_v          = config%getc_v()
      cav%dp_ref = 2.d0/(config%getrhoref()*cav%uref**2)
      cav%t_ref = cav%uref/config%getl_chord()
         
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(cav)
      implicit none
      class(t_cav), intent(inout) :: cav
      
      if(associated(cav%grd))  nullify(cav%grd) 
      if(associated(cav%pv))   nullify(cav%pv) 
      if(associated(cav%dv))   nullify(cav%dv)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function merkle(cav,eos) result(cav_result)
      implicit none
      class(t_merkle), intent(in) :: cav
      type(t_eos), intent(in) :: eos
      type(t_cav_result) :: cav_result
      real(8) :: r_c,r_v,pww,lam,rho1
      real(8), parameter :: phi = 1.d0
      
      cav_result = t_cav_result(0.d0,(/0.d0,0.d0,0.d0,0.d0/))
      pww = eos%get_pww(cav%pv(5))
      rho1 = 1.d0/cav%dv(1)
      
      r_c = cav%c_c*cav%dv(1)*cav%pv(6)*(1.d0-cav%pv(7))*dmax1(cav%pv(1)+cav%pref-pww,0.d0)*cav%dp_ref*cav%t_ref
      r_v = cav%c_v*cav%dv(1)*(1.d0-cav%pv(6))*dmax1(pww-cav%pv(1)-cav%pref,0.d0)*cav%dp_ref*cav%t_ref
           
      cav_result%cavsource = (r_v - r_c)*cav%grd(1)
      
      if(r_c.ne.0.d0) then
        cav_result%icav(1) = - r_c*(1.d0/(cav%pv(1)+cav%pref-pww)+cav%dv(7)*rho1)
        cav_result%icav(2) = - r_c*cav%dv(8)*rho1
        cav_result%icav(3) = - r_c*(1.d0/cav%pv(6)+cav%dv(9)*rho1)
        cav_result%icav(4) = - r_c*(-1.d0/(1.d0-cav%pv(7))+cav%dv(10)*rho1)
      end if
      if(r_v.ne.0.d0) then
        cav_result%icav(1) = r_v*(-1.d0/(pww-cav%pv(1)-cav%pref)+cav%dv(7)*rho1)
        cav_result%icav(2) = r_v*cav%dv(8)*rho1
        cav_result%icav(3) = r_v*(cav%dv(9)*rho1-1.d0/(1.d0-cav%pv(6)))
        cav_result%icav(4) = r_v*cav%dv(10)*rho1
      end if
      lam = dsign(1.d0,cav_result%icav(3))
      cav_result%icav(:) = cav_result%icav(:)*(phi*(1.d0-0.5d0*(1.d0-lam))-0.5d0*(1.d0-lam))
    end function merkle
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function kunz(cav,eos) result(cav_result)
      implicit none
      class(t_kunz), intent(in) :: cav
      type(t_eos), intent(in) :: eos
      type(t_cav_result) :: cav_result
      real(8) :: r_c,r_v,al,av,ag,pww,av1,lam,rho1
      real(8), parameter :: phi = 1.d0

      cav_result = t_cav_result(0.d0,(/0.d0,0.d0,0.d0,0.d0/))
      pww = eos%get_pww(cav%pv(5))
      rho1 = 1.d0/cav%dv(1)
      
      av = cav%dv(1)*cav%pv(6)*(1.d0-cav%pv(7))/cav%dv(4)
      ag = cav%dv(1)*cav%pv(6)*cav%pv(7)/cav%dv(5)
      al = 1.d0 - av - ag

      r_c = cav%c_c*cav%dv(4)*av*(al+ag)**2*cav%t_ref
      r_v = cav%c_v*cav%dv(1)*(1.d0-cav%pv(6))*dmax1(pww-cav%pv(1)-cav%pref,0.d0)*cav%dp_ref*cav%t_ref

      cav_result%cavsource = (r_v - r_c)*cav%grd(1)
      
      av1 = 1.d0-4.d0*av+3.d0*av**2
      if(r_c.ne.0.d0 ) then
        cav_result%icav(1) = - cav%c_c*cav%t_ref*(cav%dv(7)*cav%pv(6)*(1.d0-cav%pv(7))*av1+2.d0*av**2*cav%dv(15)*(al+ag))
        cav_result%icav(2) = - cav%c_c*cav%t_ref*(cav%dv(8)*cav%pv(6)*(1.d0-cav%pv(7))*av1+2.d0*av**2*cav%dv(16)*(al+ag))
        cav_result%icav(3) = - cav%c_c*cav%t_ref*(cav%dv(9)*cav%pv(6)+cav%dv(1))*av1*(1.d0-cav%pv(7))
        cav_result%icav(4) = - cav%c_c*cav%t_ref*(cav%dv(10)*(1.d0-cav%pv(7))-cav%dv(1))*cav%pv(6)*av1
      end if
      if(r_v.ne.0.d0) then
        cav_result%icav(1) = r_v*(-1.d0/(pww-cav%pv(1)-cav%pref)+cav%dv(7)*rho1)
        cav_result%icav(2) = r_v*cav%dv(8)*rho1
        cav_result%icav(3) = r_v*(cav%dv(9)*rho1-1.d0/(1.d0-cav%pv(6)))
        cav_result%icav(4) = r_v*cav%dv(10)*rho1
      end if
      lam = dsign(1.d0,cav_result%icav(3))
      cav_result%icav(:) = cav_result%icav(:)*(phi*(1.d0-0.5d0*(1.d0-lam))-0.5d0*(1.d0-lam))
    end function kunz  
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function singhal(cav,eos) result(cav_result)
      implicit none
      class(t_singhal), intent(in) :: cav
      type(t_eos), intent(in) :: eos
      type(t_cav_result) :: cav_result
      real(8) :: r_c,r_v,pww,sigma,lam,rho1,rho2
      real(8), parameter :: phi = 1.d0
      
      cav_result = t_cav_result(0.d0,(/0.d0,0.d0,0.d0,0.d0/))
      pww = eos%get_pww(cav%pv(5))
      sigma = eos%get_sigma(cav%pv(5))
      rho1 = 1.d0/cav%dv(3)
      rho2 = 1.d0/cav%dv(4)
      
      r_c = cav%c_c*cav%uref/sigma*cav%dv(3)*cav%dv(4) &
          *dsqrt(2.d0/3.d0*dmax1(cav%pv(1)+cav%pref-pww,0.d0)*rho1)*cav%pv(6)*(1.d0-cav%pv(7))
      r_v = cav%c_v*cav%uref/sigma*cav%dv(3)*cav%dv(4) &
          *dsqrt(2.d0/3.d0*dmax1(pww-cav%pv(1)-cav%pref,0.d0)*rho1)*(1.d0-cav%pv(6))

      cav_result%cavsource = (r_v - r_c)*cav%grd(1)
      
      if(r_c.ne.0.d0 ) then
        cav_result%icav(1) = - r_c*(cav%dv(15)*rho2+0.5d0*cav%dv(17)*rho1+0.5d0/(cav%pv(1)+cav%pref-pww))
        cav_result%icav(2) = - r_c*(cav%dv(16)*rho2+0.5d0*cav%dv(18)*rho1)
        cav_result%icav(3) = - r_c/cav%pv(6)
        cav_result%icav(4) = + r_c/(1.d0-cav%pv(7))
      end if
      if(r_v.ne.0.d0 ) then
        cav_result%icav(1) = r_v*(cav%dv(15)*rho2+0.5d0*cav%dv(17)*rho1-0.5d0/(pww-cav%pv(1)-cav%pref))
        cav_result%icav(2) = r_v*(cav%dv(16)*rho2+0.5d0*cav%dv(18)*rho1)
        cav_result%icav(3) = - r_v/(1.d0-cav%pv(6))
        cav_result%icav(4) = 0.d0
      end if
      lam = dsign(1.d0,cav_result%icav(3))
      cav_result%icav(:) = cav_result%icav(:)*(phi*(1.d0-0.5d0*(1.d0-lam))-0.5d0*(1.d0-lam))
    end function singhal
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module cav_module 