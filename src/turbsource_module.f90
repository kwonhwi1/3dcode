module turbsource_module
  use config_module
  implicit none
  private
  public :: t_turbsource,t_kepsilon,t_kwsst, t_turb_result
  
  type t_turb_result
    real(8) :: source(2),itt(4),omega_cut
  end type
  
  type, abstract :: t_turbsource
    private
    logical :: tcomp
    real(8), pointer, public :: cx1(:),cx2(:),ex1(:),ex2(:),tx1(:),tx2(:)
    real(8), pointer, public :: pv(:,:)
    real(8), pointer, public :: dv(:),grd(:),tv(:)
    contains
      procedure :: construct
      procedure :: destruct
      procedure(p_calturbsource), deferred :: calturbsource
      procedure, private :: calbigf
  end type t_turbsource
  
  type, extends(t_turbsource) :: t_kepsilon
    contains
      procedure ::  calturbsource => kepsilon
  end type t_kepsilon
  
  type, extends(t_turbsource) :: t_kwsst
    contains
      procedure ::  calturbsource => kwsst
  end type t_kwsst
  
  abstract interface
    function p_calturbsource(ts) result(turb_result)
      import t_turb_result
      import t_turbsource
      implicit none
      class(t_turbsource), intent(in) :: ts
      type(t_turb_result) :: turb_result
    end function
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(turbsource,config)
      implicit none
      class(t_turbsource), intent(out) :: turbsource
      type(t_config), intent(in) :: config
     
      select case(config%gettcomp())
      case(0)
        turbsource%tcomp = .false.
      case(1)
        turbsource%tcomp = .true.
      end select

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(turbsource)
      implicit none
      class(t_turbsource), intent(inout) :: turbsource

      if(associated(turbsource%cx1))     nullify(turbsource%cx1)     
      if(associated(turbsource%cx2))     nullify(turbsource%cx2)     
      if(associated(turbsource%ex1))     nullify(turbsource%ex1)     
      if(associated(turbsource%ex2))     nullify(turbsource%ex2)
      if(associated(turbsource%tx1))     nullify(turbsource%tx1)     
      if(associated(turbsource%tx2))     nullify(turbsource%tx2)
      if(associated(turbsource%grd))     nullify(turbsource%grd)           
      if(associated(turbsource%pv))      nullify(turbsource%pv)      
      if(associated(turbsource%dv))      nullify(turbsource%dv)      
      if(associated(turbsource%tv))      nullify(turbsource%tv)      
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function kepsilon(ts) result(turb_result)
      implicit none
      class(t_kepsilon), intent(in) :: ts
      type(t_turb_result) :: turb_result
      real(8) :: k1,k2,k3,k4,k5,k6
      real(8) :: u1,u2,u3,u4,u5,u6
      real(8) :: v1,v2,v3,v4,v5,v6
      real(8) :: w1,w2,w3,w4,w5,w6
      real(8) :: dkdx,dkdy,dkdz,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vol
      real(8) :: kkyy,prod,re_t,fe1,fe2,mt2,ec,pd,dpdk,dpdo
      real(8) :: lam1,lam2,t,d
      real(8), parameter :: phi = 1.d0, c23=2.d0/3.d0
      
      k1 = 0.5d0*(ts%pv(1,8)+ts%pv(2,8))
      k2 = 0.5d0*(ts%pv(3,8)+ts%pv(2,8))
      k3 = 0.5d0*(ts%pv(4,8)+ts%pv(2,8))
      k4 = 0.5d0*(ts%pv(5,8)+ts%pv(2,8))
      k5 = 0.5d0*(ts%pv(6,8)+ts%pv(2,8))
      k6 = 0.5d0*(ts%pv(7,8)+ts%pv(2,8))
      
      u1 = 0.5d0*(ts%pv(1,2)+ts%pv(2,2))
      u2 = 0.5d0*(ts%pv(3,2)+ts%pv(2,2))
      u3 = 0.5d0*(ts%pv(4,2)+ts%pv(2,2))
      u4 = 0.5d0*(ts%pv(5,2)+ts%pv(2,2))
      u5 = 0.5d0*(ts%pv(6,2)+ts%pv(2,2))
      u6 = 0.5d0*(ts%pv(7,2)+ts%pv(2,2))
      
      v1 = 0.5d0*(ts%pv(1,3)+ts%pv(2,3))
      v2 = 0.5d0*(ts%pv(3,3)+ts%pv(2,3))
      v3 = 0.5d0*(ts%pv(4,3)+ts%pv(2,3))
      v4 = 0.5d0*(ts%pv(5,3)+ts%pv(2,3))
      v5 = 0.5d0*(ts%pv(6,3)+ts%pv(2,3))
      v6 = 0.5d0*(ts%pv(7,3)+ts%pv(2,3))

      w1 = 0.5d0*(ts%pv(1,4)+ts%pv(2,4))
      w2 = 0.5d0*(ts%pv(3,4)+ts%pv(2,4))
      w3 = 0.5d0*(ts%pv(4,4)+ts%pv(2,4))
      w4 = 0.5d0*(ts%pv(5,4)+ts%pv(2,4))
      w5 = 0.5d0*(ts%pv(6,4)+ts%pv(2,4))
      w6 = 0.5d0*(ts%pv(7,4)+ts%pv(2,4))
      
      vol = 1.d0/ts%grd(1)
      
      dkdx = (dsqrt(k2)*ts%cx2(1)+dsqrt(k4)*ts%ex2(1)+dsqrt(k6)*ts%tx2(1)-dsqrt(k1)*ts%cx1(1)-dsqrt(k3)*ts%ex1(1)-dsqrt(k5)*ts%tx1(1))*vol
      dkdy = (dsqrt(k2)*ts%cx2(2)+dsqrt(k4)*ts%ex2(2)+dsqrt(k6)*ts%tx2(2)-dsqrt(k1)*ts%cx1(2)-dsqrt(k3)*ts%ex1(2)-dsqrt(k5)*ts%tx1(2))*vol
      dkdz = (dsqrt(k2)*ts%cx2(3)+dsqrt(k4)*ts%ex2(3)+dsqrt(k6)*ts%tx2(3)-dsqrt(k1)*ts%cx1(3)-dsqrt(k3)*ts%ex1(3)-dsqrt(k5)*ts%tx1(3))*vol
      
      dudx = (u2*ts%cx2(1)+u4*ts%ex2(1)+u6*ts%tx2(1)-u1*ts%cx1(1)-u3*ts%ex1(1)-u5*ts%tx1(1))*vol
      dudy = (u2*ts%cx2(2)+u4*ts%ex2(2)+u6*ts%tx2(2)-u1*ts%cx1(2)-u3*ts%ex1(2)-u5*ts%tx1(2))*vol
      dudz = (u2*ts%cx2(3)+u4*ts%ex2(3)+u6*ts%tx2(3)-u1*ts%cx1(3)-u3*ts%ex1(3)-u5*ts%tx1(3))*vol
      
      dvdx = (v2*ts%cx2(1)+v4*ts%ex2(1)+v6*ts%tx2(1)-v1*ts%cx1(1)-v3*ts%ex1(1)-v5*ts%tx1(1))*vol
      dvdy = (v2*ts%cx2(2)+v4*ts%ex2(2)+v6*ts%tx2(2)-v1*ts%cx1(2)-v3*ts%ex1(2)-v5*ts%tx1(2))*vol
      dvdz = (v2*ts%cx2(3)+v4*ts%ex2(3)+v6*ts%tx2(3)-v1*ts%cx1(3)-v3*ts%ex1(3)-v5*ts%tx1(3))*vol
      
      dwdx = (w2*ts%cx2(1)+w4*ts%ex2(1)+w6*ts%tx2(1)-w1*ts%cx1(1)-w3*ts%ex1(1)-w5*ts%tx1(1))*vol
      dwdy = (w2*ts%cx2(2)+w4*ts%ex2(2)+w6*ts%tx2(2)-w1*ts%cx1(2)-w3*ts%ex1(2)-w5*ts%tx1(2))*vol
      dwdz = (w2*ts%cx2(3)+w4*ts%ex2(3)+w6*ts%tx2(3)-w1*ts%cx1(3)-w3*ts%ex1(3)-w5*ts%tx1(3))*vol
      
      kkyy = dkdx**2+dkdy**2+dkdz**2

      prod = ts%tv(3)*(2.d0*(dudx**2+dvdy**2+dwdz**2)+(dudy+dvdx)**2+(dudz+dwdx)**2+(dvdz+dwdy)**2         &
           - c23*(dudx+dvdy+dwdz)**2) -c23*ts%dv(1)*ts%pv(2,8)*(dudx+dvdy+dwdz)
      
      re_t = ts%dv(1)*ts%pv(2,8)**2/ts%pv(2,9)/ts%tv(1)
      
      fe1 = 1.d0 - dexp(-(re_t/40.d0)**2)
      fe2 = 1.d0 - 2.d0/9.d0*dexp(-(re_t/6.d0)**2)
      
      prod = dmin1(prod,30.d0*ts%dv(1)*ts%pv(2,9))
      
      if(ts%tcomp) then
        mt2 = 2.d0*ts%pv(2,8)/ts%dv(6)
        ec = mt2*ts%pv(2,9)
        pd = -0.4d0*prod*mt2+0.2d0*ts%dv(1)*ts%pv(2,9)*mt2
      else
        mt2 = 0.d0
        ec = 0.d0
        pd = 0.d0      
      end if
      
      turb_result%omega_cut = 0.d0
      turb_result%source(1) = (prod - ts%dv(1)*(ts%pv(2,9)+ec) + pd - 2.d0*ts%tv(1)*kkyy)*ts%grd(1)
      turb_result%source(2) = (1.5d0*fe1*prod*ts%pv(2,9) - 1.9d0*fe2*ts%dv(1)*ts%pv(2,9)**2 &
                            + 2.9556d0*ts%tv(1)*ts%pv(2,9)*kkyy)/ts%pv(2,8)*ts%grd(1)
     
      
      dpdk = 2.d0*(prod+c23*ts%dv(1)*ts%pv(2,8)*(dudx+dvdy+dwdz))/ts%pv(2,8)-c23*ts%dv(1)*(dudx+dvdy+dwdz)
      dpdo = - (prod+c23*ts%dv(1)*ts%pv(2,8)*(dudx+dvdy+dwdz))/ts%pv(2,9)
      
      turb_result%itt(1) = dpdk
      turb_result%itt(2) = dpdo - ts%dv(1)
      turb_result%itt(3) = (1.9d0*fe2*ts%dv(1)*ts%pv(2,9)**2 - 2.9556d0*ts%tv(1)*ts%pv(2,9)*kkyy )/ts%pv(2,8)**2
      turb_result%itt(4) = (1.5d0*fe1*dpdo*ts%pv(2,9) + 1.5d0*fe1*prod - 3.8d0*fe2*ts%dv(1)*ts%pv(2,9) + 2.9556d0*ts%tv(1)*kkyy)/ts%pv(2,8)
      
      if(ts%tcomp) then
        turb_result%itt(1) = turb_result%itt(1) + (- 1.6d0*ts%dv(1)*ts%pv(2,9) - 0.8d0*dpdk*ts%pv(2,8)- 0.8d0*prod)/ts%dv(6)
        turb_result%itt(2) = turb_result%itt(2) + (- ts%dv(1) - 0.4d0*dpdo + 0.2d0*ts%dv(1))*mt2
      end if
      
      t = turb_result%itt(1)+turb_result%itt(4)
      d = turb_result%itt(1)*turb_result%itt(4)-turb_result%itt(2)*turb_result%itt(3)
      lam1 = 0.5d0*t + dsqrt(dmax1(0.25d0*t**2-d,0.d0))
      lam2 = 0.5d0*t - dsqrt(dmax1(0.25d0*t**2-d,0.d0))
      lam1 = dsign(1.d0,lam1)
      lam2 = dsign(1.d0,lam2)
      
      turb_result%itt(1:2) = turb_result%itt(1:2)*(phi*(1.d0-0.5d0*(1.d0-lam1))-0.5d0*(1.d0-lam1))
      turb_result%itt(3:4) = turb_result%itt(3:4)*(phi*(1.d0-0.5d0*(1.d0-lam2))-0.5d0*(1.d0-lam2))
    end function kepsilon
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function kwsst(ts) result(turb_result)
      implicit none
      class(t_kwsst), intent(in) :: ts
      type(t_turb_result) :: turb_result
      real(8) :: k1,k2,k3,k4,k5,k6
      real(8) :: o1,o2,o3,o4,o5,o6
      real(8) :: u1,u2,u3,u4,u5,u6
      real(8) :: v1,v2,v3,v4,v5,v6
      real(8) :: w1,w2,w3,w4,w5,w6
      real(8) :: dkdx,dkdy,dkdz,dodx,dody,dodz
      real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vol
      real(8) :: prod,bigf,talpha,tbeta,mt2,pd,dpdk,dpdo
      real(8) :: lam1,lam2,t,d
      real(8), parameter :: phi = 1.d0, c23=2.d0/3.d0, c59 = 5.d0/9.d0

      k1 = 0.5d0*(ts%pv(1,8)+ts%pv(2,8))
      k2 = 0.5d0*(ts%pv(3,8)+ts%pv(2,8))
      k3 = 0.5d0*(ts%pv(4,8)+ts%pv(2,8))
      k4 = 0.5d0*(ts%pv(5,8)+ts%pv(2,8))
      k5 = 0.5d0*(ts%pv(6,8)+ts%pv(2,8))
      k6 = 0.5d0*(ts%pv(7,8)+ts%pv(2,8))

      o1 = 0.5d0*(ts%pv(1,9)+ts%pv(2,9))
      o2 = 0.5d0*(ts%pv(3,9)+ts%pv(2,9))
      o3 = 0.5d0*(ts%pv(4,9)+ts%pv(2,9))
      o4 = 0.5d0*(ts%pv(5,9)+ts%pv(2,9))
      o5 = 0.5d0*(ts%pv(6,9)+ts%pv(2,9))
      o6 = 0.5d0*(ts%pv(7,9)+ts%pv(2,9))
      
      u1 = 0.5d0*(ts%pv(1,2)+ts%pv(2,2))
      u2 = 0.5d0*(ts%pv(3,2)+ts%pv(2,2))
      u3 = 0.5d0*(ts%pv(4,2)+ts%pv(2,2))
      u4 = 0.5d0*(ts%pv(5,2)+ts%pv(2,2))
      u5 = 0.5d0*(ts%pv(6,2)+ts%pv(2,2))
      u6 = 0.5d0*(ts%pv(7,2)+ts%pv(2,2))
      
      v1 = 0.5d0*(ts%pv(1,3)+ts%pv(2,3))
      v2 = 0.5d0*(ts%pv(3,3)+ts%pv(2,3))
      v3 = 0.5d0*(ts%pv(4,3)+ts%pv(2,3))
      v4 = 0.5d0*(ts%pv(5,3)+ts%pv(2,3))
      v5 = 0.5d0*(ts%pv(6,3)+ts%pv(2,3))
      v6 = 0.5d0*(ts%pv(7,3)+ts%pv(2,3))

      w1 = 0.5d0*(ts%pv(1,4)+ts%pv(2,4))
      w2 = 0.5d0*(ts%pv(3,4)+ts%pv(2,4))
      w3 = 0.5d0*(ts%pv(4,4)+ts%pv(2,4))
      w4 = 0.5d0*(ts%pv(5,4)+ts%pv(2,4))
      w5 = 0.5d0*(ts%pv(6,4)+ts%pv(2,4))
      w6 = 0.5d0*(ts%pv(7,4)+ts%pv(2,4))
      
      vol = 1.d0/ts%grd(1)
      dkdx = (k2*ts%cx2(1)+k4*ts%ex2(1)+k6*ts%tx2(1)-k1*ts%cx1(1)-k3*ts%ex1(1)-k5*ts%tx1(1))*vol
      dkdy = (k2*ts%cx2(2)+k4*ts%ex2(2)+k6*ts%tx2(2)-k1*ts%cx1(2)-k3*ts%ex1(2)-k5*ts%tx1(2))*vol
      dkdz = (k2*ts%cx2(3)+k4*ts%ex2(3)+k6*ts%tx2(3)-k1*ts%cx1(3)-k3*ts%ex1(3)-k5*ts%tx1(3))*vol

      dodx = (o2*ts%cx2(1)+o4*ts%ex2(1)+o6*ts%tx2(1)-o1*ts%cx1(1)-o3*ts%ex1(1)-o5*ts%tx1(1))*vol
      dody = (o2*ts%cx2(2)+o4*ts%ex2(2)+o6*ts%tx2(2)-o1*ts%cx1(2)-o3*ts%ex1(2)-o5*ts%tx1(2))*vol
      dodz = (o2*ts%cx2(3)+o4*ts%ex2(3)+o6*ts%tx2(3)-o1*ts%cx1(3)-o3*ts%ex1(3)-o5*ts%tx1(3))*vol
      
      dudx = (u2*ts%cx2(1)+u4*ts%ex2(1)+u6*ts%tx2(1)-u1*ts%cx1(1)-u3*ts%ex1(1)-u5*ts%tx1(1))*vol
      dudy = (u2*ts%cx2(2)+u4*ts%ex2(2)+u6*ts%tx2(2)-u1*ts%cx1(2)-u3*ts%ex1(2)-u5*ts%tx1(2))*vol
      dudz = (u2*ts%cx2(3)+u4*ts%ex2(3)+u6*ts%tx2(3)-u1*ts%cx1(3)-u3*ts%ex1(3)-u5*ts%tx1(3))*vol
      
      dvdx = (v2*ts%cx2(1)+v4*ts%ex2(1)+v6*ts%tx2(1)-v1*ts%cx1(1)-v3*ts%ex1(1)-v5*ts%tx1(1))*vol
      dvdy = (v2*ts%cx2(2)+v4*ts%ex2(2)+v6*ts%tx2(2)-v1*ts%cx1(2)-v3*ts%ex1(2)-v5*ts%tx1(2))*vol
      dvdz = (v2*ts%cx2(3)+v4*ts%ex2(3)+v6*ts%tx2(3)-v1*ts%cx1(3)-v3*ts%ex1(3)-v5*ts%tx1(3))*vol
      
      dwdx = (w2*ts%cx2(1)+w4*ts%ex2(1)+w6*ts%tx2(1)-w1*ts%cx1(1)-w3*ts%ex1(1)-w5*ts%tx1(1))*vol
      dwdy = (w2*ts%cx2(2)+w4*ts%ex2(2)+w6*ts%tx2(2)-w1*ts%cx1(2)-w3*ts%ex1(2)-w5*ts%tx1(2))*vol
      dwdz = (w2*ts%cx2(3)+w4*ts%ex2(3)+w6*ts%tx2(3)-w1*ts%cx1(3)-w3*ts%ex1(3)-w5*ts%tx1(3))*vol
      
      prod = ts%tv(3)*(2.d0*(dudx**2+dvdy**2+dwdz**2)+(dudy+dvdx)**2+(dudz+dwdx)**2+(dvdz+dwdy)**2         &
           - c23*(dudx+dvdy+dwdz)**2) - c23*ts%dv(1)*ts%pv(2,8)*(dudx+dvdy+dwdz) 
      
      prod = dmin1(prod,0.9d0*ts%dv(1)*ts%pv(2,8)*ts%pv(2,9))
      
      turb_result%omega_cut = dsqrt(2.d0*(dudx**2+dvdy**2+dwdz**2)+(dudy+dvdx)**2+(dudz+dwdx)**2+(dvdz+dwdy)**2     &
                              - c23*(dudx+dvdy+dwdz)**2)
      
      bigf = ts%calbigf(dkdx*dodx+dkdy*dody+dkdz*dodz)
      
      talpha = bigf*c59 +(1.d0-bigf)*0.44d0
      tbeta  = bigf*0.075d0   +(1.d0-bigf)*0.0828d0
      
      if(ts%tcomp) then
        mt2 = 2.d0*ts%pv(2,8)/ts%dv(6)
        pd = -0.4d0*prod*mt2+0.018d0*ts%dv(1)*ts%pv(2,8)*ts%pv(2,9)*mt2
      else
        mt2 = 0.d0
        pd = 0.d0      
      end if
      
      turb_result%source(1) = (prod - 0.09d0*ts%dv(1)*ts%pv(2,8)*ts%pv(2,9)*(1.d0+c59*mt2*(1.d0-bigf)) &
                            + (1.d0-bigf)*pd)*ts%grd(1)
      turb_result%source(2) = (talpha*ts%dv(1)/ts%tv(3)*prod - tbeta*ts%dv(1)*ts%pv(2,9)**2 &
                            + 1.712d0*(1.d0-bigf)*ts%dv(1)/ts%pv(2,9)*(dkdx*dodx+dkdy*dody+dkdz*dodz) &
                            + (1.d0-bigf)*0.05d0*mt2*ts%dv(1)*ts%pv(2,9)**2 &
                            - (1.d0-bigf)*ts%dv(1)/ts%tv(3)*pd )*ts%grd(1)

      dpdk = prod/ts%pv(2,8)   
      dpdo = - (prod+c23*ts%dv(1)*ts%pv(2,8)*(dudx+dvdy+dwdz))/ts%pv(2,9)

      turb_result%itt(1) = dpdk - 0.09d0*ts%dv(1)*ts%pv(2,8)
      turb_result%itt(2) = dpdo - 0.09d0*ts%dv(1)*ts%pv(2,9)
      turb_result%itt(3) = 0.d0
      turb_result%itt(4) = talpha*ts%dv(1)/ts%tv(3)*dpdo + talpha*ts%dv(1)/ts%tv(3)*prod/ts%pv(2,9) - 2.d0*tbeta*ts%dv(1)*ts%pv(2,9) &
                         - 1.712d0*(1.d0-bigf)*ts%dv(1)/ts%pv(2,9)**2*(dkdx*dodx+dkdy*dody+dkdz*dodz)
                 
      if(ts%tcomp) then
        turb_result%itt(1) = turb_result%itt(1) + (- 0.128d0*ts%dv(1)*ts%pv(2,8)*ts%pv(2,9)/ts%dv(6) - 0.4d0*dpdk*mt2 - 0.8d0*prod/ts%dv(6) )*(1.d0-bigf)
        turb_result%itt(2) = turb_result%itt(2) + (- 0.032d0*ts%dv(1)*ts%pv(2,8)*mt2 - 0.4d0*dpdo*mt2 )*(1.d0-bigf)
        turb_result%itt(3) = turb_result%itt(3) + (0.1d0*ts%dv(1)*ts%pv(2,9)**2/ts%dv(6) + ts%dv(1)/ts%tv(3)*pd/ts%pv(2,8) &
                                                - ts%dv(1)/ts%tv(3)*(-0.4d0*dpdk*mt2-0.8d0*prod/ts%dv(6)+0.072d0*ts%dv(1)*ts%pv(2,8)*ts%pv(2,9)/ts%dv(6)))*(1.d0-bigf)
        turb_result%itt(4) = turb_result%itt(4) + (0.2d0*ts%dv(1)*ts%pv(2,9)*ts%pv(2,8)/ts%dv(6) - ts%dv(1)/ts%tv(3)*pd/ts%pv(2,9) &
                                                - ts%dv(1)/ts%tv(3)*(-0.4d0*dpdo*mt2+0.018d0*ts%dv(1)*ts%pv(2,8)*mt2))*(1.d0-bigf)
      end if
                 
      t = turb_result%itt(1)+turb_result%itt(4)
      d = turb_result%itt(1)*turb_result%itt(4)-turb_result%itt(2)*turb_result%itt(3)
      lam1 = 0.5d0*t + dsqrt(dmax1(0.25d0*t**2-d,0.d0))
      lam2 = 0.5d0*t - dsqrt(dmax1(0.25d0*t**2-d,0.d0))
      lam1 = dsign(1.d0,lam1)
      lam2 = dsign(1.d0,lam2)
      
      turb_result%itt(1:2) = turb_result%itt(1:2)*(phi*(1.d0-0.5d0*(1.d0-lam1))-0.5d0*(1.d0-lam1))
      turb_result%itt(3:4) = turb_result%itt(3:4)*(phi*(1.d0-0.5d0*(1.d0-lam2))-0.5d0*(1.d0-lam2))
      
    end function kwsst
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function calbigf(ts,pkw)
      implicit none
      class(t_turbsource), intent(in) :: ts
      real(8), intent(in) :: pkw
      real(8) :: calbigf
      real(8) :: term0,term1,term2,term3,cdkw,arg1,arg2
      
      term0 = 1.712d0*ts%dv(1)/ts%pv(2,9)*pkw
      cdkw = dmax1(term0,1.d-10)
      
      term1 = dsqrt(ts%pv(2,8))/(0.09d0*ts%pv(2,9)*ts%grd(5))
      term2 = 500.d0*ts%tv(1)/ts%dv(1)/(ts%pv(2,9)*ts%grd(5)**2)
      term3 = (3.424d0*ts%dv(1)*ts%pv(2,8))/(cdkw*ts%grd(5)**2)
      
      arg2 = dmax1(term1,term2)
      arg1 = dmin1(arg2,term3)
      
      calbigf = dtanh(arg1**4)
      
    end function calbigf
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module turbsource_module
