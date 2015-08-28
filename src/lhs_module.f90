module lhs_module
  use config_module
  use grid_module
  use variable_module
  implicit none
  private
  public :: t_lhs,t_lhs_flowonly,t_lhs_flowturball
  
  type, abstract :: t_lhs
    private
    integer :: npv,ndv,ntv,ngrd,nsteady
    real(8) :: uref,str,dt_phy
    real(8), pointer :: cx1(:),cx2(:),ex1(:),ex2(:),tx1(:),tx2(:)
    real(8), pointer :: pv(:),tv(:),dv(:),grd(:),dt
    real(8), pointer :: c(:),t(:)
    real(8), pointer :: c0(:),t0(:)
    real(8), dimension(:,:), allocatable :: x
    procedure(p_getsndp2), pointer :: getsndp2
    procedure(p_geteigenvis), pointer :: geteigenvis
    procedure(p_eigen), pointer :: eigen
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
      procedure :: getx
      procedure(p_inverse), deferred :: inverse
  end type t_lhs
  
  type, extends(t_lhs) :: t_lhs_flowonly
    contains
      procedure :: inverse => flowonly
  end type t_lhs_flowonly
  
  type, extends(t_lhs) :: t_lhs_flowturball
    contains
      procedure :: inverse => flowturball
  end type t_lhs_flowturball
  
  abstract interface
    subroutine p_inverse(lhs)
      import t_lhs
      implicit none
      class(t_lhs), intent(inout) :: lhs
    end subroutine p_inverse  
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
    
    function p_eigen(lhs) result(lamda)
      import t_lhs
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: lamda
    end function p_eigen
    
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(lhs,config,grid,variable)
      implicit none
      class(t_lhs), intent(out) :: lhs
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      
      lhs%ngrd = grid%getngrd()
      lhs%npv = variable%getnpv()
      lhs%ndv = variable%getndv()
      lhs%ntv = variable%getntv()

      lhs%uref = config%geturef()
      lhs%str  = config%getstr()
      lhs%dt_phy = config%getdt_phy()
      lhs%nsteady = config%getnsteady()
      
      select case(config%getprec())
      case(0)
        lhs%getsndp2 => no_prec
      case(1)
        lhs%getsndp2 => steady_prec
      case(2)
        lhs%getsndp2 => unsteady_prec
      end select
      

      select case(config%gettimemethod())
      case(1,2)
        lhs%eigen => eigen_explicit
        lhs%geteigenvis => euler
      case(3)
        lhs%eigen => eigen_implicit
        select case(config%getiturb())
        case(-1,0)
          lhs%geteigenvis => turbulent
        case(-2)
          lhs%geteigenvis => laminar
        case(-3)
          lhs%geteigenvis => euler
        end select
      end select
      
      allocate(lhs%x(lhs%npv,lhs%npv))
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
      if(associated(lhs%eigen))       nullify(lhs%eigen)        

      deallocate(lhs%x,lhs%c0,lhs%t0)
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
    subroutine flowonly(lhs)
      implicit none
      class(t_lhs_flowonly), intent(inout) :: lhs
      real(8) :: uv2,ss,rr,hst,b1,dr,dhdp1,dr1,dr2,hst1
      real(8) :: cross1,cross2,cross3,a,b
      real(8) :: mm,ky1,ky2,x15,x55,x65,i_b1
      real(8) :: c(4)

      a = 1.d0/(lhs%grd(1)/lhs%dt + lhs%eigen())
      b = (dble(lhs%nsteady)*1.5d0*lhs%grd(1)/lhs%dt_phy + 2.d0*lhs%geteigenvis()/lhs%grd(1))*a      
      c = lhs%c*lhs%grd(1)*a
      
      b1 = 1.d0+b
      i_b1 = 1.d0/(b1*lhs%dv(1))
      uv2 = lhs%pv(2)**2+lhs%pv(3)**2+lhs%pv(4)**2
      ss = lhs%dv(1)*lhs%dv(12)*(1.d0/lhs%getsndp2(uv2)+b/lhs%dv(6))
      rr = 1.d0/lhs%getsndp2(uv2)-1.d0/lhs%dv(6)+b1*lhs%dv(7)
      hst1 = lhs%dv(2)-0.5d0*uv2-lhs%dv(14)*lhs%pv(7)
      hst = hst1 - lhs%dv(13)*lhs%pv(6)
      dr1 = lhs%dv(1)+lhs%dv(10)*lhs%pv(7)
      dr2 = lhs%dv(1)+lhs%dv(9)*lhs%pv(6)
      dr = dr1 + lhs%dv(9)*lhs%pv(6)
      dhdp1 = -1.d0 + lhs%dv(11)*lhs%dv(1)
      cross1 = lhs%dv(14)*lhs%dv(8)-lhs%dv(12)*lhs%dv(10)
      cross2 = lhs%dv(13)*lhs%dv(10)-lhs%dv(14)*lhs%dv(9)
      cross3 = lhs%dv(13)*lhs%dv(8)-lhs%dv(12)*lhs%dv(9)
      ky1 = lhs%dv(13)*lhs%dv(1)*rr - lhs%dv(9)*b1*dhdp1
      ky2 = lhs%dv(14)*lhs%dv(1)*rr - lhs%dv(10)*b1*dhdp1
      mm = 1.d0/(ss*b1*lhs%dv(1)+c(1)*b1*lhs%dv(1)*cross3-c(2)*ky1+c(3)*ss)
      x15 = ( c(2)*lhs%dv(9)-lhs%dv(8)*(c(3)+b1*lhs%dv(1)))*mm
      x55 = -(b1*c(1)*lhs%dv(9)-(c(3)+lhs%dv(1)*b1)*rr)*mm/b1
      x65 =  (b1*c(1)*lhs%dv(8)-c(2)*rr)*mm/b1
      
      lhs%x(1,1) = (b1*lhs%dv(1)*( lhs%dv(8)*hst + lhs%dv(12)*dr ) &
                 - c(2)*( lhs%dv(9)*hst1 + lhs%dv(13)*dr1 ) &
                 + c(3)*( lhs%dv(8)*hst1 + lhs%dv(12)*dr1 ) &
                 + c(4)*cross3*lhs%pv(7) )*mm
      lhs%x(1,2) = -x15*lhs%pv(2)
      lhs%x(1,3) = -x15*lhs%pv(3)
      lhs%x(1,4) = -x15*lhs%pv(4)
      lhs%x(1,5) = x15
      lhs%x(1,6) = b1*cross3*lhs%dv(1)*mm
      lhs%x(1,7) = ( b1*lhs%dv(1)*cross1 + c(2)*cross2 + c(3)*cross1 - c(4)*cross3 )*mm

      lhs%x(2,1) = -lhs%pv(2)*i_b1
      lhs%x(2,2) = i_b1
      lhs%x(2,3) = 0.d0
      lhs%x(2,4) = 0.d0
      lhs%x(2,5) = 0.d0
      lhs%x(2,6) = 0.d0
      lhs%x(2,7) = 0.d0

      lhs%x(3,1) = -lhs%pv(3)*i_b1
      lhs%x(3,2) = 0.d0
      lhs%x(3,3) = i_b1
      lhs%x(3,4) = 0.d0
      lhs%x(3,5) = 0.d0
      lhs%x(3,6) = 0.d0
      lhs%x(3,7) = 0.d0

      lhs%x(4,1) = -lhs%pv(4)*i_b1
      lhs%x(4,2) = 0.d0
      lhs%x(4,3) = 0.d0
      lhs%x(4,4) = i_b1
      lhs%x(4,5) = 0.d0
      lhs%x(4,6) = 0.d0
      lhs%x(4,7) = 0.d0
      
      lhs%x(5,1) = ( - b1*lhs%dv(1)*( hst*lhs%dv(1)*rr+b1*dhdp1*dr )   &
                 + c(1)*b1*lhs%dv(1)*(lhs%dv(9)*hst1 + lhs%dv(13)*dr2) &
                 - c(3)*(lhs%dv(1)*rr*hst1+b1*dhdp1*dr2)               & 
                 - c(4)*ky1*lhs%pv(7) )*mm*i_b1
      lhs%x(5,2) = -x55*lhs%pv(2)
      lhs%x(5,3) = -x55*lhs%pv(3)
      lhs%x(5,4) = -x55*lhs%pv(4)
      lhs%x(5,5) = x55
      lhs%x(5,6) = -ky1*mm
      lhs%x(5,7) = ( -lhs%dv(1)*b1*(ky2 + c(1)*cross2) - c(3)*ky2 + c(4)*ky1 )*mm*i_b1

      
      lhs%x(6,1) = ( - b1*lhs%dv(1)*ss*lhs%pv(6) - c(1)*b1*lhs%dv(1)*(lhs%dv(8)*hst1+lhs%dv(12)*dr2) &
                 + c(2)*(lhs%dv(1)*rr*hst1 + b1*dhdp1*dr2) + c(4)*ss*lhs%pv(7) )*mm*i_b1
      lhs%x(6,2) = - x65*lhs%pv(2)
      lhs%x(6,3) = - x65*lhs%pv(3)
      lhs%x(6,4) = - x65*lhs%pv(4)
      lhs%x(6,5) = x65
      lhs%x(6,6) = ss*mm
      lhs%x(6,7) = ( -c(1)*b1*lhs%dv(1)*cross1 + c(2)*ky2 - c(4)*ss)*mm*i_b1

      lhs%x(7,1) = -lhs%pv(7)*i_b1
      lhs%x(7,2) = 0.d0
      lhs%x(7,3) = 0.d0
      lhs%x(7,4) = 0.d0
      lhs%x(7,5) = 0.d0
      lhs%x(7,6) = 0.d0
      lhs%x(7,7) = i_b1
      
      lhs%x(:,:) = lhs%x(:,:)*a
    end subroutine flowonly
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine flowturball(lhs)
      implicit none
      class(t_lhs_flowturball), intent(inout) :: lhs
      real(8) :: uv2,ss,rr,hst,b1,dr,dhdp1,dr1,dr2,hst1
      real(8) :: cross1,cross2,cross3,a,b
      real(8) :: mm,ky1,ky2,x15,x55,x65,tt,i_b1
      real(8) :: c(4),t(4)

      a = 1.d0/(lhs%grd(1)/lhs%dt + lhs%eigen())
      b = (dble(lhs%nsteady)*1.5d0*lhs%grd(1)/lhs%dt_phy + 2.d0*lhs%geteigenvis()/lhs%grd(1))*a      
      c = lhs%c*lhs%grd(1)*a
      t = lhs%t*lhs%grd(1)*a
      
      b1 = 1.d0+b
      i_b1 = 1.d0/(b1*lhs%dv(1))
      uv2 = lhs%pv(2)**2+lhs%pv(3)**2+lhs%pv(4)**2
      ss = lhs%dv(1)*lhs%dv(12)*(1.d0/lhs%getsndp2(uv2)+b/lhs%dv(6))
      rr = 1.d0/lhs%getsndp2(uv2)-1.d0/lhs%dv(6)+b1*lhs%dv(7)
      hst1 = lhs%dv(2)-0.5d0*uv2-lhs%dv(14)*lhs%pv(7)
      hst = hst1 - lhs%dv(13)*lhs%pv(6)
      dr1 = lhs%dv(1)+lhs%dv(10)*lhs%pv(7)
      dr2 = lhs%dv(1)+lhs%dv(9)*lhs%pv(6)
      dr = dr1 + lhs%dv(9)*lhs%pv(6)
      dhdp1 = -1.d0 + lhs%dv(11)*lhs%dv(1)
      cross1 = lhs%dv(14)*lhs%dv(8)-lhs%dv(12)*lhs%dv(10)
      cross2 = lhs%dv(13)*lhs%dv(10)-lhs%dv(14)*lhs%dv(9)
      cross3 = lhs%dv(13)*lhs%dv(8)-lhs%dv(12)*lhs%dv(9)
      ky1 = lhs%dv(13)*lhs%dv(1)*rr - lhs%dv(9)*b1*dhdp1
      ky2 = lhs%dv(14)*lhs%dv(1)*rr - lhs%dv(10)*b1*dhdp1
      mm = 1.d0/(ss*b1*lhs%dv(1)+c(1)*b1*lhs%dv(1)*cross3-c(2)*ky1+c(3)*ss)
      x15 = ( c(2)*lhs%dv(9)-lhs%dv(8)*(c(3)+b1*lhs%dv(1)))*mm
      x55 = -(b1*c(1)*lhs%dv(9)-(c(3)+lhs%dv(1)*b1)*rr)*mm/b1
      x65 =  (b1*c(1)*lhs%dv(8)-c(2)*rr)*mm/b1
      tt = 1.d0/( t(2)*t(3)-(t(1)+b1*lhs%dv(1))*(t(4)+b1*lhs%dv(1)) )

      lhs%x(1,1) = (b1*lhs%dv(1)*( lhs%dv(8)*hst + lhs%dv(12)*dr ) &
                 - c(2)*( lhs%dv(9)*hst1 + lhs%dv(13)*dr1 ) &
                 + c(3)*( lhs%dv(8)*hst1 + lhs%dv(12)*dr1 ) &
                 + c(4)*cross3*lhs%pv(7) )*mm
      lhs%x(1,2) = -x15*lhs%pv(2)
      lhs%x(1,3) = -x15*lhs%pv(3)
      lhs%x(1,4) = -x15*lhs%pv(4)
      lhs%x(1,5) = x15
      lhs%x(1,6) = b1*cross3*lhs%dv(1)*mm
      lhs%x(1,7) = ( b1*lhs%dv(1)*cross1 + c(2)*cross2 + c(3)*cross1 - c(4)*cross3 )*mm
      lhs%x(1,8) = 0.d0
      lhs%x(1,9) = 0.d0
      
      lhs%x(2,1) = -lhs%pv(2)*i_b1
      lhs%x(2,2) = i_b1
      lhs%x(2,3) = 0.d0
      lhs%x(2,4) = 0.d0
      lhs%x(2,5) = 0.d0
      lhs%x(2,6) = 0.d0
      lhs%x(2,7) = 0.d0
      lhs%x(2,8) = 0.d0
      lhs%x(2,9) = 0.d0
      
      lhs%x(3,1) = -lhs%pv(3)*i_b1
      lhs%x(3,2) = 0.d0
      lhs%x(3,3) = i_b1
      lhs%x(3,4) = 0.d0
      lhs%x(3,5) = 0.d0
      lhs%x(3,6) = 0.d0
      lhs%x(3,7) = 0.d0
      lhs%x(3,8) = 0.d0
      lhs%x(3,9) = 0.d0
      
      lhs%x(4,1) = -lhs%pv(4)*i_b1
      lhs%x(4,2) = 0.d0
      lhs%x(4,3) = 0.d0
      lhs%x(4,4) = i_b1
      lhs%x(4,5) = 0.d0
      lhs%x(4,6) = 0.d0
      lhs%x(4,7) = 0.d0
      lhs%x(4,8) = 0.d0
      lhs%x(4,9) = 0.d0
      
      lhs%x(5,1) = ( - b1*lhs%dv(1)*( hst*lhs%dv(1)*rr+b1*dhdp1*dr )   &
                 + c(1)*b1*lhs%dv(1)*(lhs%dv(9)*hst1 + lhs%dv(13)*dr2) &
                 - c(3)*(lhs%dv(1)*rr*hst1+b1*dhdp1*dr2)               & 
                 - c(4)*ky1*lhs%pv(7) )*mm*i_b1
      lhs%x(5,2) = -x55*lhs%pv(2)
      lhs%x(5,3) = -x55*lhs%pv(3)
      lhs%x(5,4) = -x55*lhs%pv(4)
      lhs%x(5,5) = x55
      lhs%x(5,6) = -ky1*mm
      lhs%x(5,7) = ( -lhs%dv(1)*b1*(ky2 + c(1)*cross2) - c(3)*ky2 + c(4)*ky1 )*mm*i_b1
      lhs%x(5,8) = 0.d0
      lhs%x(5,9) = 0.d0
        
      lhs%x(6,1) = ( - b1*lhs%dv(1)*ss*lhs%pv(6) - c(1)*b1*lhs%dv(1)*(lhs%dv(8)*hst1+lhs%dv(12)*dr2) &
                 + c(2)*(lhs%dv(1)*rr*hst1 + b1*dhdp1*dr2) + c(4)*ss*lhs%pv(7) )*mm*i_b1
      lhs%x(6,2) = - x65*lhs%pv(2)
      lhs%x(6,3) = - x65*lhs%pv(3)
      lhs%x(6,4) = - x65*lhs%pv(4)
      lhs%x(6,5) = x65
      lhs%x(6,6) = ss*mm
      lhs%x(6,7) = ( -c(1)*b1*lhs%dv(1)*cross1 + c(2)*ky2 - c(4)*ss)*mm*i_b1
      lhs%x(6,8) = 0.d0
      lhs%x(6,9) = 0.d0
        
      lhs%x(7,1) = -lhs%pv(7)*i_b1
      lhs%x(7,2) = 0.d0
      lhs%x(7,3) = 0.d0
      lhs%x(7,4) = 0.d0
      lhs%x(7,5) = 0.d0
      lhs%x(7,6) = 0.d0
      lhs%x(7,7) = i_b1
      lhs%x(7,8) = 0.d0
      lhs%x(7,9) = 0.d0
       
      lhs%x(8,1) = ( - lhs%pv(9)*t(2) + lhs%pv(8)*(b1*lhs%dv(1)+t(4)) )*tt
      lhs%x(8,2) = 0.d0
      lhs%x(8,3) = 0.d0
      lhs%x(8,4) = 0.d0
      lhs%x(8,5) = 0.d0
      lhs%x(8,6) = 0.d0
      lhs%x(8,7) = 0.d0
      lhs%x(8,8) = -(b1*lhs%dv(1)+t(4))*tt
      lhs%x(8,9) = t(2)*tt

      lhs%x(9,1) = ( lhs%pv(9)*(b1*lhs%dv(1)+t(1)) - lhs%pv(8)*t(3) )*tt
      lhs%x(9,2) = 0.d0
      lhs%x(9,3) = 0.d0
      lhs%x(9,4) = 0.d0
      lhs%x(9,5) = 0.d0
      lhs%x(9,6) = 0.d0
      lhs%x(9,7) = 0.d0
      lhs%x(9,8) = t(3)*tt
      lhs%x(9,9) = -(b1*lhs%dv(1)+t(1))*tt
      
      lhs%x(:,:) = lhs%x(:,:)*a      
    end subroutine flowturball
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function eigen_explicit(lhs) result(lamda)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8) :: lamda

      lamda = 0.d0
      
    end function eigen_explicit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function eigen_implicit(lhs) result(lamda)
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
      up = (1.d0+sndp2/lhs%dv(6))*u
      d  = dsqrt(u**2*(1.d0-sndp2/lhs%dv(6))**2+4.d0*sndp2*ds)
      lamda = cappa*(dabs(up)+d)

      ox = 0.5d0*(lhs%ex1(1)+lhs%ex2(1))
      oy = 0.5d0*(lhs%ex1(2)+lhs%ex2(2))
      oz = 0.5d0*(lhs%ex1(3)+lhs%ex2(3))
      ds = ox**2+oy**2+oz**2
      u  = ox*lhs%pv(2)+oy*lhs%pv(3)+oz*lhs%pv(4)
      up = (1.d0+sndp2/lhs%dv(6))*u
      d  = dsqrt(u**2*(1.d0-sndp2/lhs%dv(6))**2+4.d0*sndp2*ds)
      lamda = lamda+cappa*(dabs(up)+d)
      
      ox = 0.5d0*(lhs%tx1(1)+lhs%tx2(1))
      oy = 0.5d0*(lhs%tx1(2)+lhs%tx2(2))
      oz = 0.5d0*(lhs%tx1(3)+lhs%tx2(3))
      ds = ox**2+oy**2+oz**2
      u  = ox*lhs%pv(2)+oy*lhs%pv(3)+oz*lhs%pv(4)
      up = (1.d0+sndp2/lhs%dv(6))*u
      d  = dsqrt(u**2*(1.d0-sndp2/lhs%dv(6))**2+4.d0*sndp2*ds)
      lamda = lamda+cappa*(dabs(up)+d)
      
    end function eigen_implicit
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
      
      sndp2 = dmin1(lhs%dv(6),dmax1(uuu2,lhs%uref**2))
      
    end function steady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteady_prec(lhs,uuu2) result(sndp2)
      implicit none
      class(t_lhs), intent(in) :: lhs
      real(8), intent(in) :: uuu2
      real(8) :: sndp2  
      
      sndp2 = dmin1(lhs%dv(6),dmax1(uuu2,lhs%uref**2,lhs%str**2*uuu2))
      
    end function unsteady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getx(lhs,i,j)
      implicit none
      class(t_lhs), intent(in) :: lhs
      integer, intent(in) :: i,j
      real(8) :: getx
      
      getx = lhs%x(i,j)
      
    end function getx
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module lhs_module
