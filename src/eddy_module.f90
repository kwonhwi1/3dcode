module eddy_module
  use grid_module
  use variable_module
  implicit none
  private
  public :: t_eddy,t_eddy_ke,t_eddy_kwsst
  
  type, abstract :: t_eddy
    private
    integer :: ngrd,npv,ndv,ntv
    real(8), pointer :: pv(:,:)
    real(8), pointer :: cx1(:),cx2(:),ex1(:),ex2(:),tx1(:),tx2(:)
    real(8), pointer :: dv(:),tv(:),grd(:)
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setnorm ! (cx1,cx2,ex1,ex2,tx1,tx2)
      procedure :: setgrd  ! vol,xcen,ycen,zcen,ydns
      procedure :: setpv   ! p,u,v,w,t,y1,y2,k,o
      procedure :: setdv   ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      procedure :: settv   ! vis,cond,emut
      procedure(p_caleddy), deferred :: caleddy
  end type t_eddy
  
  abstract interface
    function p_caleddy(eddy) result(emut)
      import t_eddy
      class(t_eddy), intent(in) :: eddy
      real(8) :: emut
    end function p_caleddy
  end interface
  
  type, extends(t_eddy) :: t_eddy_ke
    contains
      procedure :: caleddy => kepsilon
  end type t_eddy_ke

  type, extends(t_eddy) :: t_eddy_kwsst
    contains
      procedure :: caleddy => kwsst
  end type t_eddy_kwsst

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(eddy,grid,variable)
      implicit none
      class(t_eddy), intent(out) :: eddy
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
    
      eddy%ngrd = grid%getngrd()
      eddy%npv  = variable%getnpv()
      eddy%ndv  = variable%getndv()
      eddy%ntv  = variable%getntv()      

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(eddy)
      implicit none
      class(t_eddy), intent(inout) :: eddy
      
      if(associated(eddy%cx1))     nullify(eddy%cx1)
      if(associated(eddy%cx2))     nullify(eddy%cx2)
      if(associated(eddy%ex1))     nullify(eddy%ex1)
      if(associated(eddy%ex2))     nullify(eddy%ex2)
      if(associated(eddy%tx1))     nullify(eddy%tx1)
      if(associated(eddy%tx2))     nullify(eddy%tx2)
      if(associated(eddy%grd))     nullify(eddy%grd)
      if(associated(eddy%pv))      nullify(eddy%pv) 
      if(associated(eddy%dv))      nullify(eddy%dv) 
      if(associated(eddy%tv))      nullify(eddy%tv)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setnorm(eddy,cx1,cx2,ex1,ex2,tx1,tx2)
      implicit none
      class(t_eddy), intent(inout) :: eddy
      real(8), intent(in), target :: cx1(3),cx2(3),ex1(3),ex2(3),tx1(3),tx2(3)
      
      eddy%cx1 => cx1
      eddy%cx2 => cx2
      eddy%ex1 => ex1
      eddy%ex2 => ex2
      eddy%tx1 => tx1
      eddy%tx2 => tx2
      
    end subroutine setnorm
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(eddy,grd)
      implicit none
      class(t_eddy), intent(inout) :: eddy
      real(8), intent(in), target :: grd(eddy%ngrd)
      
      eddy%grd => grd
    
    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(eddy,pv)
      implicit none
      class(t_eddy), intent(inout) :: eddy
      real(8), intent(in), target :: pv(7,eddy%npv)
      
      eddy%pv => pv
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(eddy,dv)
      implicit none
      class(t_eddy), intent(inout) :: eddy
      real(8), intent(in), target :: dv(eddy%ndv)
      
      eddy%dv => dv
    
    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine settv(eddy,tv)
      implicit none
      class(t_eddy), intent(inout) :: eddy
      real(8), intent(in), target :: tv(eddy%ntv)
      
      eddy%tv => tv
    
    end subroutine settv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function kepsilon(eddy) result(emut)
      class(t_eddy_ke), intent(in) :: eddy
      real(8) :: emut
      real(8) :: re_t,re_e,f_mu
      
      re_t = eddy%dv(1)*eddy%pv(2,8)**2/eddy%pv(2,9)/eddy%tv(1)
      re_e = eddy%dv(1)*(eddy%tv(1)/eddy%dv(1)*eddy%pv(2,9))**0.25d0*eddy%grd(5)/eddy%tv(1)
      f_mu = (1.d0+3.d0/re_t**(3.d0/4.d0))*(1.d0+80.d0*dexp(-re_e))*(1.d0-dexp(-re_e/43.d0-re_e**2/330.d0))**2
      emut = 0.09d0*eddy%dv(1)*eddy%pv(2,8)**2/eddy%pv(2,9)*f_mu    
    end function kepsilon
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function kwsst(eddy) result(emut)
      class(t_eddy_kwsst), intent(in) :: eddy
      real(8) :: emut
      real(8) :: u1,u2,u3,u4,u5,u6
      real(8) :: v1,v2,v3,v4,v5,v6
      real(8) :: w1,w2,w3,w4,w5,w6
      real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vol
      real(8) :: term1,term2,arg2,bigf2
      real(8) :: sijsij,sss,term4,term5,term6

      u1 = 0.5d0*(eddy%pv(1,2)+eddy%pv(2,2))
      u2 = 0.5d0*(eddy%pv(3,2)+eddy%pv(2,2))
      u3 = 0.5d0*(eddy%pv(4,2)+eddy%pv(2,2))
      u4 = 0.5d0*(eddy%pv(5,2)+eddy%pv(2,2))
      u5 = 0.5d0*(eddy%pv(6,2)+eddy%pv(2,2))
      u6 = 0.5d0*(eddy%pv(7,2)+eddy%pv(2,2))

      v1 = 0.5d0*(eddy%pv(1,3)+eddy%pv(2,3))
      v2 = 0.5d0*(eddy%pv(3,3)+eddy%pv(2,3))
      v3 = 0.5d0*(eddy%pv(4,3)+eddy%pv(2,3))
      v4 = 0.5d0*(eddy%pv(5,3)+eddy%pv(2,3))
      v5 = 0.5d0*(eddy%pv(6,3)+eddy%pv(2,3))
      v6 = 0.5d0*(eddy%pv(7,3)+eddy%pv(2,3))

      w1 = 0.5d0*(eddy%pv(1,4)+eddy%pv(2,4))
      w2 = 0.5d0*(eddy%pv(3,4)+eddy%pv(2,4))
      w3 = 0.5d0*(eddy%pv(4,4)+eddy%pv(2,4))
      w4 = 0.5d0*(eddy%pv(5,4)+eddy%pv(2,4))
      w5 = 0.5d0*(eddy%pv(6,4)+eddy%pv(2,4))
      w6 = 0.5d0*(eddy%pv(7,4)+eddy%pv(2,4))
      
      vol = 1.d0/eddy%grd(1)
      dudx = (u2*eddy%cx2(1)+u4*eddy%ex2(1)+u6*eddy%tx2(1)-u1*eddy%cx1(1)-u3*eddy%ex1(1)-u5*eddy%tx1(1))*vol
      dudy = (u2*eddy%cx2(2)+u4*eddy%ex2(2)+u6*eddy%tx2(2)-u1*eddy%cx1(2)-u3*eddy%ex1(2)-u5*eddy%tx1(2))*vol
      dudz = (u2*eddy%cx2(3)+u4*eddy%ex2(3)+u6*eddy%tx2(3)-u1*eddy%cx1(3)-u3*eddy%ex1(3)-u5*eddy%tx1(3))*vol
      dvdx = (v2*eddy%cx2(1)+v4*eddy%ex2(1)+v6*eddy%tx2(1)-v1*eddy%cx1(1)-v3*eddy%ex1(1)-v5*eddy%tx1(1))*vol
      dvdy = (v2*eddy%cx2(2)+v4*eddy%ex2(2)+v6*eddy%tx2(2)-v1*eddy%cx1(2)-v3*eddy%ex1(2)-v5*eddy%tx1(2))*vol
      dvdz = (v2*eddy%cx2(3)+v4*eddy%ex2(3)+v6*eddy%tx2(3)-v1*eddy%cx1(3)-v3*eddy%ex1(3)-v5*eddy%tx1(3))*vol
      dwdx = (w2*eddy%cx2(1)+w4*eddy%ex2(1)+w6*eddy%tx2(1)-w1*eddy%cx1(1)-w3*eddy%ex1(1)-w5*eddy%tx1(1))*vol
      dwdy = (w2*eddy%cx2(2)+w4*eddy%ex2(2)+w6*eddy%tx2(2)-w1*eddy%cx1(2)-w3*eddy%ex1(2)-w5*eddy%tx1(2))*vol
      dwdz = (w2*eddy%cx2(3)+w4*eddy%ex2(3)+w6*eddy%tx2(3)-w1*eddy%cx1(3)-w3*eddy%ex1(3)-w5*eddy%tx1(3))*vol
      
      term1 = dsqrt(eddy%pv(2,8))/(0.09d0*eddy%pv(2,9)*eddy%grd(5))
      term2 = 500.d0*eddy%tv(1)/eddy%dv(1)/(eddy%pv(2,9)*eddy%grd(5)**2)
      arg2 = dmax1(2.d0*term1,term2)
      bigf2 = dtanh(arg2**2)   
      sijsij = dudx**2+dvdy**2+dwdz**2+0.5d0*((dudy+dvdx)**2+(dudz+dwdx)**2+(dvdz+dwdy)**2)
      sss = dsqrt(2.d0*sijsij)
      term4 = 0.31d0*eddy%pv(2,9)
      term5 = sss*bigf2
      term6 = dmax1(term4,term5)
      emut = 0.31d0*eddy%dv(1)*eddy%pv(2,8)/term6
    end function kwsst
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module eddy_module
