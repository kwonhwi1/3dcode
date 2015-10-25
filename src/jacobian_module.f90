module jacobian_module
  use config_module
  use grid_module
  use variable_module
  implicit none
  private
  public :: t_jac,t_jac_flowonly,t_jac_flowturball
  
  type, abstract :: t_jac
    private
    integer :: ngrd,npv,ndv,ntv
    real(8) :: uref,str
    real(8) :: omega(3)
    real(8), pointer :: nx(:)
    real(8), pointer :: pv(:),tv(:),dv(:),grd(:)
    real(8), dimension(:,:), allocatable :: a
    procedure(p_getsndp2), pointer :: getsndp2
    procedure(p_geteigenvis), pointer :: geteigenvis
    procedure(p_getenthalpy), pointer :: getenthalpy
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setnorm ! (nx)
      procedure :: setgrd  ! vol,xcen,ycen,zcen,ydns
      procedure :: setpv   ! p,u,v,w,t,y1,y2,k,o
      procedure :: setdv   ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      procedure :: settv   ! vis,cond,emut
      procedure :: geta
      procedure(p_caljac), deferred :: caljac 
  end type t_jac
  
  type, extends(t_jac) :: t_jac_flowonly
    contains
      procedure :: caljac => flowonly
  end type t_jac_flowonly

  type, extends(t_jac) :: t_jac_flowturball
    contains
      procedure :: caljac => flowturball
  end type t_jac_flowturball
  
  abstract interface
    subroutine p_caljac(jac,sign)
      import t_jac
      implicit none
      class(t_jac), intent(inout) :: jac
      integer, intent(in) :: sign
    end subroutine p_caljac
  end interface
  
  interface
    function p_getsndp2(jac,u2) result(sndp2)
      import t_jac
      implicit none
      class(t_jac), intent(in) :: jac
      real(8), intent(in) :: u2
      real(8) :: sndp2
    end function p_getsndp2
    
    function p_geteigenvis(jac) result(eigenvis)
      import t_jac
      implicit none
      class(t_jac), intent(in) :: jac
      real(8) :: eigenvis
    end function p_geteigenvis

    function p_getenthalpy(jac)
      import t_jac
      implicit none
      class(t_jac), intent(in) :: jac
      real(8) :: p_getenthalpy
    end function p_getenthalpy
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(jac,config,grid,variable)
      implicit none
      class(t_jac), intent(out) :: jac
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
            
      jac%uref = config%geturef()
      jac%str  = config%getstr()
      jac%omega = config%getomega()
      
      select case(config%getprec())
      case(0)
        jac%getsndp2 => no_prec
      case(1)
        jac%getsndp2 => steady_prec
      case(2)
        jac%getsndp2 => unsteady_prec
      end select
      
      select case(config%getiturb())
      case(-1,0)
        jac%geteigenvis => turbulent
      case(-2)
        jac%geteigenvis => laminar
      case(-3)
        jac%geteigenvis => euler
      end select

      select case(config%getrotation())
      case(0)
        jac%getenthalpy => enthalpy
      case(-1,1,-2,2,-3,3)
        jac%getenthalpy => rothalpy
      case default
      end select

      jac%ngrd = grid%getngrd()
      jac%npv  = variable%getnpv()
      jac%ntv  = variable%getntv()
      jac%ndv  = variable%getndv()   
      allocate(jac%a(jac%npv,jac%npv))
    
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(jac)
      implicit none
      class(t_jac), intent(inout) :: jac

      if(associated(jac%nx))           nullify(jac%nx)                    
      if(associated(jac%grd))          nullify(jac%grd)         
      if(associated(jac%pv))           nullify(jac%pv)
      if(associated(jac%dv))           nullify(jac%dv) 
      if(associated(jac%tv))           nullify(jac%tv)                   
      if(associated(jac%getsndp2))     nullify(jac%getsndp2)    
      if(associated(jac%geteigenvis))  nullify(jac%geteigenvis)
      if(associated(jac%getenthalpy))  nullify(jac%getenthalpy)
      
      deallocate(jac%a)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setnorm(jac,nx)
      implicit none
      class(t_jac), intent(inout) :: jac
      real(8), intent(in), target :: nx(3)
      
      jac%nx => nx

    end subroutine setnorm
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(jac,grd)
      implicit none
      class(t_jac), intent(inout) :: jac
      real(8), intent(in), target :: grd(jac%ngrd)
      
      jac%grd => grd
      
    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setpv(jac,pv)
      implicit none
      class(t_jac), intent(inout) :: jac
      real(8), intent(in), target :: pv(jac%npv)
      
      jac%pv => pv
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(jac,dv)
      implicit none
      class(t_jac), intent(inout) :: jac
      real(8), intent(in), target :: dv(jac%ndv)
      
      jac%dv => dv
      
    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine settv(jac,tv)
      implicit none
      class(t_jac), intent(inout) :: jac
      real(8), intent(in), target :: tv(jac%ntv)
      
      jac%tv => tv
      
    end subroutine settv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine flowonly(jac,sign)
      implicit none
      class(t_jac_flowonly), intent(inout) :: jac
      integer, intent(in) :: sign
      real(8) :: ds,u,up,d,beta,sndp2,lamda,vis
      real(8), parameter :: cappa = 1.05d0
      real(8) :: a1,a2,a3,a4,a5,a6,a7
      real(8) :: g1,g5,g6,g7,ge1,ge5,ge6,ge7,uuu
      
      ds = jac%nx(1)**2+jac%nx(2)**2+jac%nx(3)**2
      u = jac%nx(1)*jac%pv(2) + jac%nx(2)*jac%pv(3)+jac%nx(3)*jac%pv(4)
      sndp2 = jac%getsndp2(jac%pv(2)**2+jac%pv(3)**2+jac%pv(4)**2)
      up = 0.5d0*(1.d0+sndp2/jac%dv(6))*u
      d  = 0.5d0*dsqrt(u**2*(1.d0-sndp2/jac%dv(6))**2+4.d0*sndp2*ds)
      lamda = dble(sign)*cappa*(dabs(up)+d)
      vis = dble(sign)*2.d0*jac%geteigenvis()*ds/jac%grd(1)
      beta = 1.d0/sndp2-1.d0/jac%dv(6)+jac%dv(7)
      
      a1 = u*jac%dv(7)
      a2 = jac%nx(1)*jac%dv(1)
      a3 = jac%nx(2)*jac%dv(1)
      a4 = jac%nx(3)*jac%dv(1)
      a5 = u*jac%dv(8)
      a6 = u*jac%dv(9)
      a7 = u*jac%dv(10)
      g1 = lamda*beta
      g5 = lamda*jac%dv(8)
      g6 = lamda*jac%dv(9)
      g7 = lamda*jac%dv(10)
      ge1 = vis*jac%dv(7)
      ge5 = vis*jac%dv(8)
      ge6 = vis*jac%dv(9)
      ge7 = vis*jac%dv(10)
      uuu = u+lamda+vis

      jac%a(1,1) = a1+g1+ge1
      jac%a(1,2) = a2
      jac%a(1,3) = a3
      jac%a(1,4) = a4
      jac%a(1,5) = a5+g5+ge5
      jac%a(1,6) = a6+g6+ge6
      jac%a(1,7) = a7+g7+ge7

      jac%a(2,1) = jac%pv(2)*jac%a(1,1) + jac%nx(1)
      jac%a(2,2) = jac%pv(2)*jac%a(1,2) + jac%dv(1)*uuu
      jac%a(2,3) = jac%pv(2)*jac%a(1,3)
      jac%a(2,4) = jac%pv(2)*jac%a(1,4)
      jac%a(2,5) = jac%pv(2)*jac%a(1,5)
      jac%a(2,6) = jac%pv(2)*jac%a(1,6)
      jac%a(2,7) = jac%pv(2)*jac%a(1,7)
        
      jac%a(3,1) = jac%pv(3)*jac%a(1,1) + jac%nx(2)
      jac%a(3,2) = jac%pv(3)*jac%a(1,2) 
      jac%a(3,3) = jac%pv(3)*jac%a(1,3) + jac%dv(1)*uuu
      jac%a(3,4) = jac%pv(3)*jac%a(1,4)
      jac%a(3,5) = jac%pv(3)*jac%a(1,5)
      jac%a(3,6) = jac%pv(3)*jac%a(1,6)
      jac%a(3,7) = jac%pv(3)*jac%a(1,7)

      jac%a(4,1) = jac%pv(4)*jac%a(1,1) + jac%nx(3)
      jac%a(4,2) = jac%pv(4)*jac%a(1,2) 
      jac%a(4,3) = jac%pv(4)*jac%a(1,3) 
      jac%a(4,4) = jac%pv(4)*jac%a(1,4) + jac%dv(1)*uuu
      jac%a(4,5) = jac%pv(4)*jac%a(1,5)
      jac%a(4,6) = jac%pv(4)*jac%a(1,6)
      jac%a(4,7) = jac%pv(4)*jac%a(1,7)
      
      jac%a(5,1) = jac%dv(1)*uuu*jac%dv(11)+jac%getenthalpy()*jac%a(1,1)-lamda-vis
      jac%a(5,2) = jac%dv(1)*uuu*jac%pv(2) +jac%getenthalpy()*jac%a(1,2)
      jac%a(5,3) = jac%dv(1)*uuu*jac%pv(3) +jac%getenthalpy()*jac%a(1,3)
      jac%a(5,4) = jac%dv(1)*uuu*jac%pv(4) +jac%getenthalpy()*jac%a(1,4)
      jac%a(5,5) = jac%dv(1)*uuu*jac%dv(12)+jac%getenthalpy()*jac%a(1,5)
      jac%a(5,6) = jac%dv(1)*uuu*jac%dv(13)+jac%getenthalpy()*jac%a(1,6)
      jac%a(5,7) = jac%dv(1)*uuu*jac%dv(14)+jac%getenthalpy()*jac%a(1,7)
        
      jac%a(6,1) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,1)
      jac%a(6,2) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,2)
      jac%a(6,3) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,3)
      jac%a(6,4) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,4)
      jac%a(6,5) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,5)
      jac%a(6,6) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,6) + jac%dv(1)*(1.d0-jac%pv(7))*uuu
      jac%a(6,7) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,7) - jac%dv(1)*jac%pv(6)*uuu
        
      jac%a(7,1) = jac%pv(6)*jac%pv(7)*jac%a(1,1)
      jac%a(7,2) = jac%pv(6)*jac%pv(7)*jac%a(1,2)
      jac%a(7,3) = jac%pv(6)*jac%pv(7)*jac%a(1,3)
      jac%a(7,4) = jac%pv(6)*jac%pv(7)*jac%a(1,4)
      jac%a(7,5) = jac%pv(6)*jac%pv(7)*jac%a(1,5)
      jac%a(7,6) = jac%pv(6)*jac%pv(7)*jac%a(1,6) + jac%dv(1)*jac%pv(7)*uuu
      jac%a(7,7) = jac%pv(6)*jac%pv(7)*jac%a(1,7) + jac%dv(1)*jac%pv(6)*uuu

      jac%a = jac%a*0.5d0
      
    end subroutine flowonly
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine flowturball(jac,sign)
      implicit none
      class(t_jac_flowturball), intent(inout) :: jac
      integer, intent(in) :: sign
      real(8) :: ds,u,up,d,beta,sndp2,lamda,vis
      real(8), parameter :: cappa = 1.05d0
      real(8) :: a1,a2,a3,a4,a5,a6,a7
      real(8) :: g1,g5,g6,g7,ge1,ge5,ge6,ge7,uuu
      
      ds = jac%nx(1)**2+jac%nx(2)**2+jac%nx(3)**2
      u = jac%nx(1)*jac%pv(2) + jac%nx(2)*jac%pv(3)+jac%nx(3)*jac%pv(4)
      sndp2 = jac%getsndp2(jac%pv(2)**2+jac%pv(3)**2+jac%pv(4)**2)
      up = 0.5d0*(1.d0+sndp2/jac%dv(6))*u
      d  = 0.5d0*dsqrt(u**2*(1.d0-sndp2/jac%dv(6))**2+4.d0*sndp2*ds)
      lamda = dble(sign)*cappa*(dabs(up)+d)
      vis = dble(sign)*2.d0*jac%geteigenvis()*ds/jac%grd(1)
      beta = 1.d0/sndp2-1.d0/jac%dv(6)+jac%dv(7)
      
      a1 = u*jac%dv(7)
      a2 = jac%nx(1)*jac%dv(1)
      a3 = jac%nx(2)*jac%dv(1)
      a4 = jac%nx(3)*jac%dv(1)
      a5 = u*jac%dv(8)
      a6 = u*jac%dv(9)
      a7 = u*jac%dv(10)
      g1 = lamda*beta
      g5 = lamda*jac%dv(8)
      g6 = lamda*jac%dv(9)
      g7 = lamda*jac%dv(10)
      ge1 = vis*jac%dv(7)
      ge5 = vis*jac%dv(8)
      ge6 = vis*jac%dv(9)
      ge7 = vis*jac%dv(10)
      uuu = u+lamda+vis

      jac%a(1,1) = a1+g1+ge1
      jac%a(1,2) = a2
      jac%a(1,3) = a3
      jac%a(1,4) = a4
      jac%a(1,5) = a5+g5+ge5
      jac%a(1,6) = a6+g6+ge6
      jac%a(1,7) = a7+g7+ge7
      jac%a(1,8) = 0.d0
      jac%a(1,9) = 0.d0
        
      jac%a(2,1) = jac%pv(2)*jac%a(1,1) + jac%nx(1)
      jac%a(2,2) = jac%pv(2)*jac%a(1,2) + jac%dv(1)*uuu
      jac%a(2,3) = jac%pv(2)*jac%a(1,3)
      jac%a(2,4) = jac%pv(2)*jac%a(1,4)
      jac%a(2,5) = jac%pv(2)*jac%a(1,5)
      jac%a(2,6) = jac%pv(2)*jac%a(1,6)
      jac%a(2,7) = jac%pv(2)*jac%a(1,7)
      jac%a(2,8) = 0.d0
      jac%a(2,9) = 0.d0
        
      jac%a(3,1) = jac%pv(3)*jac%a(1,1) + jac%nx(2)
      jac%a(3,2) = jac%pv(3)*jac%a(1,2) 
      jac%a(3,3) = jac%pv(3)*jac%a(1,3) + jac%dv(1)*uuu
      jac%a(3,4) = jac%pv(3)*jac%a(1,4)
      jac%a(3,5) = jac%pv(3)*jac%a(1,5)
      jac%a(3,6) = jac%pv(3)*jac%a(1,6)
      jac%a(3,7) = jac%pv(3)*jac%a(1,7)
      jac%a(3,8) = 0.d0
      jac%a(3,9) = 0.d0

      jac%a(4,1) = jac%pv(4)*jac%a(1,1) + jac%nx(3)
      jac%a(4,2) = jac%pv(4)*jac%a(1,2) 
      jac%a(4,3) = jac%pv(4)*jac%a(1,3) 
      jac%a(4,4) = jac%pv(4)*jac%a(1,4) + jac%dv(1)*uuu
      jac%a(4,5) = jac%pv(4)*jac%a(1,5)
      jac%a(4,6) = jac%pv(4)*jac%a(1,6)
      jac%a(4,7) = jac%pv(4)*jac%a(1,7)
      jac%a(4,8) = 0.d0
      jac%a(4,9) = 0.d0
      
      jac%a(5,1) = jac%dv(1)*uuu*jac%dv(11)+jac%getenthalpy()*jac%a(1,1)-lamda-vis
      jac%a(5,2) = jac%dv(1)*uuu*jac%pv(2) +jac%getenthalpy()*jac%a(1,2)
      jac%a(5,3) = jac%dv(1)*uuu*jac%pv(3) +jac%getenthalpy()*jac%a(1,3)
      jac%a(5,4) = jac%dv(1)*uuu*jac%pv(4) +jac%getenthalpy()*jac%a(1,4)
      jac%a(5,5) = jac%dv(1)*uuu*jac%dv(12)+jac%getenthalpy()*jac%a(1,5)
      jac%a(5,6) = jac%dv(1)*uuu*jac%dv(13)+jac%getenthalpy()*jac%a(1,6)
      jac%a(5,7) = jac%dv(1)*uuu*jac%dv(14)+jac%getenthalpy()*jac%a(1,7)
      jac%a(5,8) = 0.d0
      jac%a(5,9) = 0.d0
        
      jac%a(6,1) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,1)
      jac%a(6,2) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,2)
      jac%a(6,3) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,3)
      jac%a(6,4) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,4)
      jac%a(6,5) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,5)
      jac%a(6,6) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,6) + jac%dv(1)*(1.d0-jac%pv(7))*uuu
      jac%a(6,7) = jac%pv(6)*(1.d0-jac%pv(7))*jac%a(1,7) - jac%dv(1)*jac%pv(6)*uuu
      jac%a(6,8) = 0.d0
      jac%a(6,9) = 0.d0
        
      jac%a(7,1) = jac%pv(6)*jac%pv(7)*jac%a(1,1)
      jac%a(7,2) = jac%pv(6)*jac%pv(7)*jac%a(1,2)
      jac%a(7,3) = jac%pv(6)*jac%pv(7)*jac%a(1,3)
      jac%a(7,4) = jac%pv(6)*jac%pv(7)*jac%a(1,4)
      jac%a(7,5) = jac%pv(6)*jac%pv(7)*jac%a(1,5)
      jac%a(7,6) = jac%pv(6)*jac%pv(7)*jac%a(1,6) + jac%dv(1)*jac%pv(7)*uuu
      jac%a(7,7) = jac%pv(6)*jac%pv(7)*jac%a(1,7) + jac%dv(1)*jac%pv(6)*uuu
      jac%a(7,8) = 0.d0
      jac%a(7,9) = 0.d0

      jac%a(8,1) = jac%pv(8)*jac%a(1,1)
      jac%a(8,2) = jac%pv(8)*jac%a(1,2)
      jac%a(8,3) = jac%pv(8)*jac%a(1,3)
      jac%a(8,4) = jac%pv(8)*jac%a(1,4)
      jac%a(8,5) = jac%pv(8)*jac%a(1,5)
      jac%a(8,6) = jac%pv(8)*jac%a(1,6)
      jac%a(8,7) = jac%pv(8)*jac%a(1,7)
      jac%a(8,8) = uuu*jac%dv(1)
      jac%a(8,9) = 0.d0

      jac%a(9,1) = jac%pv(9)*jac%a(1,1)
      jac%a(9,2) = jac%pv(9)*jac%a(1,2)
      jac%a(9,3) = jac%pv(9)*jac%a(1,3)
      jac%a(9,4) = jac%pv(9)*jac%a(1,4)
      jac%a(9,5) = jac%pv(9)*jac%a(1,5)
      jac%a(9,6) = jac%pv(9)*jac%a(1,6)
      jac%a(9,7) = jac%pv(9)*jac%a(1,7)
      jac%a(9,8) = 0.d0
      jac%a(9,9) = uuu*jac%dv(1)
      
      jac%a = jac%a*0.5d0
    end subroutine flowturball
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function euler(jac) result(eigenvis)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8) :: eigenvis
      
      eigenvis = 0.d0
      
    end function euler
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function laminar(jac) result(eigenvis)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8) :: eigenvis

      eigenvis = dmax1(jac%dv(7)*jac%tv(2)/(jac%dv(8)     &
                       +jac%dv(12)*jac%dv(7)*jac%dv(1)    &
                       -jac%dv(11)*jac%dv(8)*jac%dv(1)),  &
                       4.d0/3.d0*jac%tv(1)/jac%dv(1)) 
    end function laminar
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function turbulent(jac) result(eigenvis)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8) :: eigenvis
      real(8), parameter :: pr_t = 0.9d0
   
      eigenvis = dmax1(jac%dv(7)*(jac%tv(2)+jac%dv(12)*jac%tv(3)/pr_t) &
                      /(jac%dv(8)+jac%dv(12)*jac%dv(7)*jac%dv(1)       &
                       -jac%dv(11)*jac%dv(8)*jac%dv(1)),               &
                       4.d0/3.d0*(jac%tv(1)+jac%tv(3))/jac%dv(1))
    end function turbulent
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function no_prec(jac,u2) result(sndp2)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8), intent(in) :: u2
      real(8) :: sndp2
      
      sndp2 = jac%dv(6)
      
    end function no_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function steady_prec(jac,u2) result(sndp2)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8), intent(in) :: u2
      real(8) :: sndp2  
      
      sndp2 = dmin1(jac%dv(6),dmax1(u2,jac%uref**2))
      
    end function steady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteady_prec(jac,u2) result(sndp2)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8), intent(in) :: u2
      real(8) :: sndp2  
      
      sndp2 = dmin1(jac%dv(6),dmax1(u2,jac%uref**2,jac%str**2*u2))
      
    end function unsteady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function geta(jac,i,j)
      implicit none
      class(t_jac), intent(in) :: jac
      integer, intent(in) :: i,j
      real(8) :: geta
      
      geta = jac%a(i,j)
      
    end function geta
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function enthalpy(jac)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8) :: enthalpy
      enthalpy = jac%dv(2) + 0.5d0*(jac%pv(2)**2 + jac%pv(3)**2 + jac%pv(4)**2)
    end function enthalpy
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rothalpy(jac)
      implicit none
      class(t_jac), intent(in) :: jac
      real(8) :: rothalpy
      rothalpy = jac%dv(2) + 0.5d0*(jac%pv(2)**2 + jac%pv(3)**2 + jac%pv(4)**2)
      rothalpy = rothalpy -0.5d0*(jac%omega(1)**2*(jac%grd(3)**2+jac%grd(4)**2) &
                                + jac%omega(2)**2*(jac%grd(2)**2+jac%grd(4)**2) &
                                + jac%omega(3)**2*(jac%grd(2)**2+jac%grd(3)**2) &
                           -2.d0*(jac%omega(1)*jac%omega(2)*jac%grd(2)*jac%grd(3) &
                                + jac%omega(2)*jac%omega(3)*jac%grd(3)*jac%grd(4) &
                                + jac%omega(3)*jac%omega(1)*jac%grd(4)*jac%grd(2)))
    end function rothalpy
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module jacobian_module
