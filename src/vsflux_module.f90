module vsflux_module
  use config_module
  use grid_module
  use variable_module
  implicit none
  private
  public :: t_vsflux,t_vsflux_laminar,t_vsflux_turbulent

  type, abstract :: t_vsflux
    private
    integer :: stencil,ngrd,npv,ndv,ntv
    real(8) :: tdk1,tdk2,tdo1,tdo2
    real(8) :: pr_t
    real(8), pointer :: pv(:,:)
    real(8), pointer :: nx(:)
    real(8), pointer :: ex1(:),ex2(:),ex3(:),ex4(:)
    real(8), pointer :: tx1(:),tx2(:),tx3(:),tx4(:)
    real(8), pointer :: grdl(:),grdr(:)
    real(8), pointer :: dvl(:),dvr(:)
    real(8), pointer :: tvl(:),tvr(:)
    procedure(p_calbigf), pointer :: calbigf
    contains
      procedure :: construct  
      procedure :: destruct
      procedure :: setnorm    ! (nx,ex1,ex2,ex3,ex4,tx1,tx2,tx3,tx4)
      procedure :: setgrd     ! (grdl,grdr) vol,xcen,ycen,zcen,ydns
      procedure :: setpv      ! (pv(14,npv)) p,u,v,w,t,y1,y2,k,o
      procedure :: setdv      ! (dvl,dvr) rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      procedure :: settv      ! (tvl,tvr) vis,cond,emut
      procedure(p_calflux), deferred :: calflux
  end type t_vsflux
 
  abstract interface
    subroutine p_calflux(vsflux,fx)
      import t_vsflux
      implicit none
      class(t_vsflux), intent(in) :: vsflux
      real(8), intent(out) :: fx(vsflux%npv)
    end subroutine p_calflux
  end interface
  
  type, extends(t_vsflux) :: t_vsflux_laminar
    contains
    procedure :: calflux => vsflux_laminar
  end type t_vsflux_laminar
  
  type, extends(t_vsflux) :: t_vsflux_turbulent
    contains
    procedure :: calflux => vsflux_turbulent
  end type t_vsflux_turbulent
  
  interface
    function p_calbigf(vsflux,pkw) result(bigf)
      import t_vsflux
      implicit none
      class(t_vsflux), intent(in) :: vsflux
      real(8), intent(in) :: pkw
      real(8) :: bigf
    end function p_calbigf
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    subroutine construct(vsflux,config,grid,variable)
      implicit none
      class(t_vsflux), intent(out) :: vsflux
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      
      select case(config%getiturb())
      case(-1)
        vsflux%tdk1 = 0.d0
        vsflux%tdk2 = 1.d0
        vsflux%tdo1 = 0.d0
        vsflux%tdo2 = 1.d0/1.3d0
        vsflux%calbigf => kepsilon
        vsflux%pr_t=0.9d0
      case(0)
        vsflux%tdk1 = 0.85d0
        vsflux%tdk2 = 1.d0
        vsflux%tdo1 = 0.5d0
        vsflux%tdo2 = 0.856d0
        vsflux%calbigf => kwsst
        vsflux%pr_t=0.9d0
      case(-2)
        vsflux%tdk1 = 0.d0
        vsflux%tdk2 = 0.d0
        vsflux%tdo1 = 0.d0
        vsflux%tdo2 = 0.d0
        vsflux%calbigf => null()
        vsflux%pr_t=0.d0
      end select
      
      vsflux%stencil = config%getstencil()
      vsflux%ngrd = grid%getngrd()
      vsflux%npv = variable%getnpv()
      vsflux%ndv = variable%getndv()
      vsflux%ntv = variable%getntv()
      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine destruct(vsflux)
      implicit none
      class(t_vsflux), intent(inout) :: vsflux

      if(associated(vsflux%nx))      nullify(vsflux%nx)
      if(associated(vsflux%ex1))     nullify(vsflux%ex1)    
      if(associated(vsflux%ex2))     nullify(vsflux%ex2)    
      if(associated(vsflux%ex3))     nullify(vsflux%ex3)    
      if(associated(vsflux%ex4))     nullify(vsflux%ex4)
      if(associated(vsflux%tx1))     nullify(vsflux%tx1)    
      if(associated(vsflux%tx2))     nullify(vsflux%tx2)    
      if(associated(vsflux%tx3))     nullify(vsflux%tx3)    
      if(associated(vsflux%tx4))     nullify(vsflux%tx4)
      if(associated(vsflux%grdl))    nullify(vsflux%grdl)   
      if(associated(vsflux%grdr))    nullify(vsflux%grdr) 
      if(associated(vsflux%pv))      nullify(vsflux%pv)     
      if(associated(vsflux%dvl))     nullify(vsflux%dvl)    
      if(associated(vsflux%dvr))     nullify(vsflux%dvr)     
      if(associated(vsflux%tvl))     nullify(vsflux%tvl)    
      if(associated(vsflux%tvr))     nullify(vsflux%tvr)    
      if(associated(vsflux%calbigf)) nullify(vsflux%calbigf)
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setnorm(vsflux,nx,ex1,ex2,ex3,ex4,tx1,tx2,tx3,tx4)
      implicit none
      class(t_vsflux), intent(inout) :: vsflux
      real(8), intent(in), target :: nx(3),ex1(3),ex2(3),ex3(3),ex4(3),tx1(3),tx2(3),tx3(3),tx4(3)
      
      vsflux%nx  => nx
      vsflux%ex1 => ex1
      vsflux%ex2 => ex2
      vsflux%ex3 => ex3
      vsflux%ex4 => ex4
      vsflux%tx1 => tx1
      vsflux%tx2 => tx2
      vsflux%tx3 => tx3
      vsflux%tx4 => tx4    

    end subroutine setnorm 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(vsflux,grdl,grdr)
      implicit none
      class(t_vsflux), intent(inout) :: vsflux
      real(8), intent(in), target :: grdl(vsflux%ngrd),grdr(vsflux%ngrd)
      
      vsflux%grdl => grdl
      vsflux%grdr => grdr
      
    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(vsflux,pv)
      implicit none
      class(t_vsflux), intent(inout) :: vsflux
      real(8), intent(in), target :: pv(vsflux%stencil,vsflux%npv)
      ! pv(19,:)~pv(38,:) are not used !! plz check
      
      vsflux%pv => pv
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
    subroutine setdv(vsflux,dvl,dvr)
      implicit none
      class(t_vsflux), intent(inout) :: vsflux
      real(8), intent(in), target :: dvl(vsflux%ndv),dvr(vsflux%ndv)
      
      vsflux%dvl => dvl
      vsflux%dvr => dvr
      
    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine settv(vsflux,tvl,tvr)
      implicit none
      class(t_vsflux), intent(inout) :: vsflux
      real(8), intent(in), target :: tvl(vsflux%ntv),tvr(vsflux%ntv)
      
      vsflux%tvl => tvl
      vsflux%tvr => tvr
      
    end subroutine settv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine vsflux_laminar(vsflux,fx)
      implicit none
      class(t_vsflux_laminar), intent(in) :: vsflux
      real(8), intent(out) :: fx(vsflux%npv)
      real(8) :: tx,ty,tz,ex,ey,ez
      real(8) :: vol,vis,con,uc,vc,wc
      real(8) :: dudn,dvdn,dwdn,dtdn
      real(8) :: dude,dvde,dwde,dtde
      real(8) :: dudt,dvdt,dwdt,dtdt
      real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dtdx,dtdy,dtdz
      real(8) :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
      real(8), parameter :: c23=2.d0/3.d0

      ex = 0.25d0*(vsflux%ex1(1)+vsflux%ex2(1)+vsflux%ex3(1)+vsflux%ex4(1)) 
      ey = 0.25d0*(vsflux%ex1(2)+vsflux%ex2(2)+vsflux%ex3(2)+vsflux%ex4(2))
      ez = 0.25d0*(vsflux%ex1(3)+vsflux%ex2(3)+vsflux%ex3(3)+vsflux%ex4(3))
      
      tx = 0.25d0*(vsflux%tx1(1)+vsflux%tx2(1)+vsflux%tx3(1)+vsflux%tx4(1)) 
      ty = 0.25d0*(vsflux%tx1(2)+vsflux%tx2(2)+vsflux%tx3(2)+vsflux%tx4(2))
      tz = 0.25d0*(vsflux%tx1(3)+vsflux%tx2(3)+vsflux%tx3(3)+vsflux%tx4(3))
      
      vol  = 2.d0/(vsflux%grdl(1) + vsflux%grdr(1))
      vis  = 0.5d0*(vsflux%tvl(1) + vsflux%tvr(1))
      con  = 0.5d0*(vsflux%tvl(2) + vsflux%tvr(2))
      uc   = 0.5d0*(vsflux%pv(9,2) + vsflux%pv(10,2))
      vc   = 0.5d0*(vsflux%pv(9,3) + vsflux%pv(10,3))
      wc   = 0.5d0*(vsflux%pv(9,4) + vsflux%pv(10,4))
      
      dudn = vsflux%pv(10,2) - vsflux%pv(9,2)
      dvdn = vsflux%pv(10,3) - vsflux%pv(9,3)
      dwdn = vsflux%pv(10,4) - vsflux%pv(9,4)
      dtdn = vsflux%pv(10,5) - vsflux%pv(9,5)
      
      dude = 0.0625d0*(- vsflux%pv(1,2)  - vsflux%pv(2,2)  + vsflux%pv(5,2)  + vsflux%pv(6,2)   &
               + 2.d0*(- vsflux%pv(7,2)  - vsflux%pv(8,2)  + vsflux%pv(11,2) + vsflux%pv(12,2)) &
                       - vsflux%pv(13,2) - vsflux%pv(14,2) + vsflux%pv(17,2) + vsflux%pv(18,2))
      dvde = 0.0625d0*(- vsflux%pv(1,3)  - vsflux%pv(2,3)  + vsflux%pv(5,3)  + vsflux%pv(6,3)   &
               + 2.d0*(- vsflux%pv(7,3)  - vsflux%pv(8,3)  + vsflux%pv(11,3) + vsflux%pv(12,3)) &
                       - vsflux%pv(13,3) - vsflux%pv(14,3) + vsflux%pv(17,3) + vsflux%pv(18,3))
      dwde = 0.0625d0*(- vsflux%pv(1,4)  - vsflux%pv(2,4)  + vsflux%pv(5,4)  + vsflux%pv(6,4)   &
               + 2.d0*(- vsflux%pv(7,4)  - vsflux%pv(8,4)  + vsflux%pv(11,4) + vsflux%pv(12,4)) &
                       - vsflux%pv(13,4) - vsflux%pv(14,4) + vsflux%pv(17,4) + vsflux%pv(18,4))
      dtde = 0.0625d0*(- vsflux%pv(1,5)  - vsflux%pv(2,5)  + vsflux%pv(5,5)  + vsflux%pv(6,5)   &
               + 2.d0*(- vsflux%pv(7,5)  - vsflux%pv(8,5)  + vsflux%pv(11,5) + vsflux%pv(12,5)) &
                       - vsflux%pv(13,5) - vsflux%pv(14,5) + vsflux%pv(17,5) + vsflux%pv(18,5))
               
      dudt = 0.0625d0*(- vsflux%pv(1,2)  - vsflux%pv(2,2)  + vsflux%pv(13,2) + vsflux%pv(14,2)   &
               + 2.d0*(- vsflux%pv(3,2)  - vsflux%pv(4,2)  + vsflux%pv(15,2) + vsflux%pv(16,2))  &
                       - vsflux%pv(5,2)  - vsflux%pv(6,2)  + vsflux%pv(17,2) + vsflux%pv(18,2))
      dvdt = 0.0625d0*(- vsflux%pv(1,3)  - vsflux%pv(2,3)  + vsflux%pv(13,3) + vsflux%pv(14,3)   &
               + 2.d0*(- vsflux%pv(3,3)  - vsflux%pv(4,3)  + vsflux%pv(15,3) + vsflux%pv(16,3))  &
                       - vsflux%pv(5,3)  - vsflux%pv(6,3)  + vsflux%pv(17,3) + vsflux%pv(18,3))
      dwdt = 0.0625d0*(- vsflux%pv(1,4)  - vsflux%pv(2,4)  + vsflux%pv(13,4) + vsflux%pv(14,4)   &
               + 2.d0*(- vsflux%pv(3,4)  - vsflux%pv(4,4)  + vsflux%pv(15,4) + vsflux%pv(16,4))  &
                       - vsflux%pv(5,4)  - vsflux%pv(6,4)  + vsflux%pv(17,4) + vsflux%pv(18,4))
      dtdt = 0.0625d0*(- vsflux%pv(1,5)  - vsflux%pv(2,5)  + vsflux%pv(13,5) + vsflux%pv(14,5)   &
               + 2.d0*(- vsflux%pv(3,5)  - vsflux%pv(4,5)  + vsflux%pv(15,5) + vsflux%pv(16,5))  &
                       - vsflux%pv(5,5)  - vsflux%pv(6,5)  + vsflux%pv(17,5) + vsflux%pv(18,5))
      
      dudx = (dudn*vsflux%nx(1)+dude*ex+dudt*tx)*vol
      dvdx = (dvdn*vsflux%nx(1)+dvde*ex+dvdt*tx)*vol
      dwdx = (dwdn*vsflux%nx(1)+dwde*ex+dwdt*tx)*vol
      dtdx = (dtdn*vsflux%nx(1)+dtde*ex+dtdt*tx)*vol
 
      dudy = (dudn*vsflux%nx(2)+dude*ey+dudt*ty)*vol
      dvdy = (dvdn*vsflux%nx(2)+dvde*ey+dvdt*ty)*vol
      dwdy = (dwdn*vsflux%nx(2)+dwde*ey+dwdt*ty)*vol
      dtdy = (dtdn*vsflux%nx(2)+dtde*ey+dtdt*ty)*vol

      dudz = (dudn*vsflux%nx(3)+dude*ez+dudt*tz)*vol
      dvdz = (dvdn*vsflux%nx(3)+dvde*ez+dvdt*tz)*vol
      dwdz = (dwdn*vsflux%nx(3)+dwde*ez+dwdt*tz)*vol
      dtdz = (dtdn*vsflux%nx(3)+dtde*ez+dtdt*tz)*vol
      
      tauxx = vis*(2.d0*dudx - c23*(dudx+dvdy+dwdz))
      tauyy = vis*(2.d0*dvdy - c23*(dudx+dvdy+dwdz))
      tauzz = vis*(2.d0*dwdz - c23*(dudx+dvdy+dwdz))
      tauxy = vis*(dudy+dvdx)
      tauxz = vis*(dudz+dwdx)
      tauyz = vis*(dvdz+dwdy)
      
      fx(1) = 0.d0
      fx(2) = (tauxx*vsflux%nx(1)+tauxy*vsflux%nx(2)+tauxz*vsflux%nx(3))
      fx(3) = (tauxy*vsflux%nx(1)+tauyy*vsflux%nx(2)+tauyz*vsflux%nx(3))
      fx(4) = (tauxz*vsflux%nx(1)+tauyz*vsflux%nx(2)+tauzz*vsflux%nx(3))
      fx(5) = (vsflux%nx(1)*(uc*tauxx+vc*tauxy+wc*tauxz+con*dtdx) &
            +  vsflux%nx(2)*(uc*tauxy+vc*tauyy+wc*tauyz+con*dtdy) &
            +  vsflux%nx(3)*(uc*tauxz+vc*tauyz+wc*tauzz+con*dtdz) )
      fx(6) = 0.d0
      fx(7) = 0.d0

    end subroutine vsflux_laminar
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    subroutine vsflux_turbulent(vsflux,fx)
      implicit none
      class(t_vsflux_turbulent), intent(in) :: vsflux
      real(8), intent(out) :: fx(vsflux%npv)
      real(8) :: tx,ty,tz,ex,ey,ez
      real(8) :: vol,vis,con,emut,tcon,bigf,uc,vc,wc
      real(8) :: dudn,dvdn,dwdn,dtdn,dkdn,dodn
      real(8) :: dude,dvde,dwde,dtde,dkde,dode
      real(8) :: dudt,dvdt,dwdt,dtdt,dkdt,dodt
      real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dtdx,dtdy,dtdz
      real(8) :: dkdx,dkdy,dkdz,dodx,dody,dodz
      real(8) :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
      real(8) :: tdk,tdo
      real(8), parameter :: c23=2.d0/3.d0
      
      ex = 0.25d0*(vsflux%ex1(1)+vsflux%ex2(1)+vsflux%ex3(1)+vsflux%ex4(1)) 
      ey = 0.25d0*(vsflux%ex1(2)+vsflux%ex2(2)+vsflux%ex3(2)+vsflux%ex4(2))
      ez = 0.25d0*(vsflux%ex1(3)+vsflux%ex2(3)+vsflux%ex3(3)+vsflux%ex4(3))
      
      tx = 0.25d0*(vsflux%tx1(1)+vsflux%tx2(1)+vsflux%tx3(1)+vsflux%tx4(1)) 
      ty = 0.25d0*(vsflux%tx1(2)+vsflux%tx2(2)+vsflux%tx3(2)+vsflux%tx4(2))
      tz = 0.25d0*(vsflux%tx1(3)+vsflux%tx2(3)+vsflux%tx3(3)+vsflux%tx4(3))
      
      vol  = 2.d0/(vsflux%grdl(1) + vsflux%grdr(1))
      vis  = 0.5d0*(vsflux%tvl(1) + vsflux%tvr(1))
      con  = 0.5d0*(vsflux%tvl(2) + vsflux%tvr(2))
      emut = 0.5d0*(vsflux%tvl(3) + vsflux%tvr(3))
      tcon = 0.5d0*(vsflux%tvl(3)*vsflux%dvl(12)+ vsflux%tvr(3)*vsflux%dvr(12))/vsflux%pr_t
      uc   = 0.5d0*(vsflux%pv(9,2) + vsflux%pv(10,2))
      vc   = 0.5d0*(vsflux%pv(9,3) + vsflux%pv(10,3))
      wc   = 0.5d0*(vsflux%pv(9,4) + vsflux%pv(10,4))
      
      dudn = vsflux%pv(10,2) - vsflux%pv(9,2)
      dvdn = vsflux%pv(10,3) - vsflux%pv(9,3)
      dwdn = vsflux%pv(10,4) - vsflux%pv(9,4)
      dtdn = vsflux%pv(10,5) - vsflux%pv(9,5)
      dkdn = vsflux%pv(10,8) - vsflux%pv(9,8)
      dodn = vsflux%pv(10,9) - vsflux%pv(9,9) 
      
      dude = 0.0625d0*(- vsflux%pv(1,2)  - vsflux%pv(2,2)  + vsflux%pv(5,2)  + vsflux%pv(6,2)   &
               + 2.d0*(- vsflux%pv(7,2)  - vsflux%pv(8,2)  + vsflux%pv(11,2) + vsflux%pv(12,2)) &
                       - vsflux%pv(13,2) - vsflux%pv(14,2) + vsflux%pv(17,2) + vsflux%pv(18,2))
      dvde = 0.0625d0*(- vsflux%pv(1,3)  - vsflux%pv(2,3)  + vsflux%pv(5,3)  + vsflux%pv(6,3)   &
               + 2.d0*(- vsflux%pv(7,3)  - vsflux%pv(8,3)  + vsflux%pv(11,3) + vsflux%pv(12,3)) &
                       - vsflux%pv(13,3) - vsflux%pv(14,3) + vsflux%pv(17,3) + vsflux%pv(18,3))
      dwde = 0.0625d0*(- vsflux%pv(1,4)  - vsflux%pv(2,4)  + vsflux%pv(5,4)  + vsflux%pv(6,4)   &
               + 2.d0*(- vsflux%pv(7,4)  - vsflux%pv(8,4)  + vsflux%pv(11,4) + vsflux%pv(12,4)) &
                       - vsflux%pv(13,4) - vsflux%pv(14,4) + vsflux%pv(17,4) + vsflux%pv(18,4))
      dtde = 0.0625d0*(- vsflux%pv(1,5)  - vsflux%pv(2,5)  + vsflux%pv(5,5)  + vsflux%pv(6,5)   &
               + 2.d0*(- vsflux%pv(7,5)  - vsflux%pv(8,5)  + vsflux%pv(11,5) + vsflux%pv(12,5)) &
                       - vsflux%pv(13,5) - vsflux%pv(14,5) + vsflux%pv(17,5) + vsflux%pv(18,5))
      dkde = 0.0625d0*(- vsflux%pv(1,8)  - vsflux%pv(2,8)  + vsflux%pv(5,8)  + vsflux%pv(6,8)   &
               + 2.d0*(- vsflux%pv(7,8)  - vsflux%pv(8,8)  + vsflux%pv(11,8) + vsflux%pv(12,8)) &
                       - vsflux%pv(13,8) - vsflux%pv(14,8) + vsflux%pv(17,8) + vsflux%pv(18,8))
      dode = 0.0625d0*(- vsflux%pv(1,9)  - vsflux%pv(2,9)  + vsflux%pv(5,9)  + vsflux%pv(6,9)   &
               + 2.d0*(- vsflux%pv(7,9)  - vsflux%pv(8,9)  + vsflux%pv(11,9) + vsflux%pv(12,9)) &
                       - vsflux%pv(13,9) - vsflux%pv(14,9) + vsflux%pv(17,9) + vsflux%pv(18,9))
                       
      dudt = 0.0625d0*(- vsflux%pv(1,2)  - vsflux%pv(2,2)  + vsflux%pv(13,2) + vsflux%pv(14,2)   &
               + 2.d0*(- vsflux%pv(3,2)  - vsflux%pv(4,2)  + vsflux%pv(15,2) + vsflux%pv(16,2))  &
                       - vsflux%pv(5,2)  - vsflux%pv(6,2)  + vsflux%pv(17,2) + vsflux%pv(18,2))
      dvdt = 0.0625d0*(- vsflux%pv(1,3)  - vsflux%pv(2,3)  + vsflux%pv(13,3) + vsflux%pv(14,3)   &
               + 2.d0*(- vsflux%pv(3,3)  - vsflux%pv(4,3)  + vsflux%pv(15,3) + vsflux%pv(16,3))  &
                       - vsflux%pv(5,3)  - vsflux%pv(6,3)  + vsflux%pv(17,3) + vsflux%pv(18,3))
      dwdt = 0.0625d0*(- vsflux%pv(1,4)  - vsflux%pv(2,4)  + vsflux%pv(13,4) + vsflux%pv(14,4)   &
               + 2.d0*(- vsflux%pv(3,4)  - vsflux%pv(4,4)  + vsflux%pv(15,4) + vsflux%pv(16,4))  &
                       - vsflux%pv(5,4)  - vsflux%pv(6,4)  + vsflux%pv(17,4) + vsflux%pv(18,4))
      dtdt = 0.0625d0*(- vsflux%pv(1,5)  - vsflux%pv(2,5)  + vsflux%pv(13,5) + vsflux%pv(14,5)   &
               + 2.d0*(- vsflux%pv(3,5)  - vsflux%pv(4,5)  + vsflux%pv(15,5) + vsflux%pv(16,5))  &
                       - vsflux%pv(5,5)  - vsflux%pv(6,5)  + vsflux%pv(17,5) + vsflux%pv(18,5))
      dkdt = 0.0625d0*(- vsflux%pv(1,8)  - vsflux%pv(2,8)  + vsflux%pv(13,8) + vsflux%pv(14,8)   &
               + 2.d0*(- vsflux%pv(3,8)  - vsflux%pv(4,8)  + vsflux%pv(15,8) + vsflux%pv(16,8))  &
                       - vsflux%pv(5,8)  - vsflux%pv(6,8)  + vsflux%pv(17,8) + vsflux%pv(18,8))
      dodt = 0.0625d0*(- vsflux%pv(1,9)  - vsflux%pv(2,9)  + vsflux%pv(13,9) + vsflux%pv(14,9)   &
               + 2.d0*(- vsflux%pv(3,9)  - vsflux%pv(4,9)  + vsflux%pv(15,9) + vsflux%pv(16,9))  &
                       - vsflux%pv(5,9)  - vsflux%pv(6,9)  + vsflux%pv(17,9) + vsflux%pv(18,9))

                       
      dudx = (dudn*vsflux%nx(1)+dude*ex+dudt*tx)*vol
      dvdx = (dvdn*vsflux%nx(1)+dvde*ex+dvdt*tx)*vol
      dwdx = (dwdn*vsflux%nx(1)+dwde*ex+dwdt*tx)*vol
      dtdx = (dtdn*vsflux%nx(1)+dtde*ex+dtdt*tx)*vol
      dkdx = (dkdn*vsflux%nx(1)+dkde*ex+dkdt*tx)*vol
      dodx = (dodn*vsflux%nx(1)+dode*ex+dodt*tx)*vol
      
      dudy = (dudn*vsflux%nx(2)+dude*ey+dudt*ty)*vol
      dvdy = (dvdn*vsflux%nx(2)+dvde*ey+dvdt*ty)*vol
      dwdy = (dwdn*vsflux%nx(2)+dwde*ey+dwdt*ty)*vol
      dtdy = (dtdn*vsflux%nx(2)+dtde*ey+dtdt*ty)*vol
      dkdy = (dkdn*vsflux%nx(2)+dkde*ey+dkdt*ty)*vol
      dody = (dodn*vsflux%nx(2)+dode*ey+dodt*ty)*vol
      
      dudz = (dudn*vsflux%nx(3)+dude*ez+dudt*tz)*vol
      dvdz = (dvdn*vsflux%nx(3)+dvde*ez+dvdt*tz)*vol
      dwdz = (dwdn*vsflux%nx(3)+dwde*ez+dwdt*tz)*vol
      dtdz = (dtdn*vsflux%nx(3)+dtde*ez+dtdt*tz)*vol
      dkdz = (dkdn*vsflux%nx(3)+dkde*ez+dkdt*tz)*vol
      dodz = (dodn*vsflux%nx(3)+dode*ez+dodt*tz)*vol
      
      bigf = vsflux%calbigf(dkdx*dodx+dkdy*dody+dkdz*dodz)
      
      tdk = bigf*vsflux%tdk1 + (1.d0-bigf)*vsflux%tdk2
      tdo = bigf*vsflux%tdo1 + (1.d0-bigf)*vsflux%tdo2

      tauxx = (vis+emut)*(2.d0*dudx - c23*(dudx+dvdy+dwdz))
      tauyy = (vis+emut)*(2.d0*dvdy - c23*(dudx+dvdy+dwdz))
      tauzz = (vis+emut)*(2.d0*dwdz - c23*(dudx+dvdy+dwdz))
      tauxy = (vis+emut)*(dudy+dvdx)
      tauxz = (vis+emut)*(dudz+dwdx)
      tauyz = (vis+emut)*(dvdz+dwdy)
            
      fx(1) = 0.d0
      fx(2) = (tauxx*vsflux%nx(1)+tauxy*vsflux%nx(2)+tauxz*vsflux%nx(3))
      fx(3) = (tauxy*vsflux%nx(1)+tauyy*vsflux%nx(2)+tauyz*vsflux%nx(3))
      fx(4) = (tauxz*vsflux%nx(1)+tauyz*vsflux%nx(2)+tauzz*vsflux%nx(3))
      fx(5) = (vsflux%nx(1)*(uc*tauxx+vc*tauxy+wc*tauxz+(con+tcon)*dtdx) &
            +  vsflux%nx(2)*(uc*tauxy+vc*tauyy+wc*tauyz+(con+tcon)*dtdy) &
            +  vsflux%nx(3)*(uc*tauxz+vc*tauyz+wc*tauzz+(con+tcon)*dtdz) )
      fx(6) = 0.d0
      fx(7) = 0.d0
      fx(8) = (vis+emut*tdk)*(dkdx*vsflux%nx(1)+dkdy*vsflux%nx(2)+dkdz*vsflux%nx(3))
      fx(9) = (vis+emut*tdo)*(dodx*vsflux%nx(1)+dody*vsflux%nx(2)+dodz*vsflux%nx(3))
      
    end subroutine vsflux_turbulent
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function kwsst(vsflux,pkw) result(bigf)
      implicit none
      class(t_vsflux), intent(in) :: vsflux
      real(8), intent(in) :: pkw
      real(8) :: bigf
      real(8) :: term0,term1,term2,term3,cdkw,arg1,arg2
      real(8) :: rho,ydns,vis,k,o
      
      rho = 0.5d0*(vsflux%dvl(1)+vsflux%dvr(1))
      ydns = 0.5d0*(vsflux%grdl(5)+vsflux%grdr(5))
      vis = 0.5d0*(vsflux%tvl(1)+vsflux%tvr(1))
      k = 0.5d0*(vsflux%pv(4,8)+vsflux%pv(3,8)) 
      o = 0.5d0*(vsflux%pv(4,9)+vsflux%pv(3,9))
      
      term0 = 1.712d0*rho/o*pkw
      cdkw = dmax1(term0,1.d-10)
      
      term1 = dsqrt(k)/(0.09d0*o*ydns)
      term2 = 500.d0*vis/rho/(o*ydns**2)
      term3 = (3.424d0*rho*k)/(cdkw*ydns**2)
      
      arg2 = dmax1(term1,term2)
      arg1 = dmin1(arg2,term3)
      
      bigf = dtanh(arg1**4)
      
    end function kwsst
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function kepsilon(vsflux,pkw) result(bigf)
      implicit none
      class(t_vsflux), intent(in) :: vsflux
      real(8), intent(in) :: pkw
      real(8) :: bigf

      bigf = 0.d0
      
    end function kepsilon
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module vsflux_module
