module flux_module
  use config_module
  use grid_module
  use variable_module
  use eos_module
  implicit none
  private
  public :: t_flux,t_roe,t_roem,t_ausmpwp,t_ausmpup

  type, abstract :: t_flux
    private
    integer :: npv,ndv,ngrd
    real(8) :: omega(3)
    real(8) :: uref,str,pref
    real(8), pointer :: pvl(:),pvr(:),dvl(:),dvr(:),sdst(:)
    real(8), pointer :: nx(:),grdl(:),grdr(:)
    procedure(p_getsndp2), pointer :: getsndp2
    procedure(p_getenthalpy_l), pointer :: getenthalpy_l
    procedure(p_getenthalpy_r), pointer :: getenthalpy_r
    procedure(p_getenthalpy_c), pointer :: getenthalpy_c
    contains
      procedure :: construct       
      procedure :: destruct
      procedure :: setnorm         ! (nx)
      procedure :: setgrd          ! (grdl,grdr) vol,x,y,z,ydns
      procedure :: setpv           ! (pvl,pvr) p,u,v,w,t,y1,y2,k,o
      procedure :: setdv           ! (dvl,dvr) rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      procedure :: setsdst         ! (sdst)
      procedure(p_calflux), deferred :: calflux
  end type t_flux

  abstract interface
    subroutine p_calflux(flux,eos,fx)
      import t_flux
      import t_eos
      implicit none
      class(t_flux), intent(in) :: flux
      type(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
    end subroutine p_calflux
  end interface
  
  type, extends(t_flux) :: t_roe
    contains
      procedure :: calflux => roe
  end type t_roe

  type, extends(t_flux) :: t_roem
    contains
      procedure :: calflux => roem
  end type t_roem
  
  type, extends(t_flux) :: t_ausmpwp
    contains
      procedure :: calflux => ausmpwp
  end type t_ausmpwp
  
  type, extends(t_flux) :: t_ausmpup
    contains
      procedure :: calflux => ausmpup
  end type t_ausmpup
  
  interface
    function p_getsndp2(flux,snd2,uuu2,cut) result(sndp2)
      import t_flux
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      integer, intent(in) :: cut
      real(8) :: sndp2
    end function p_getsndp2
    function p_getenthalpy_l(flux)
      import t_flux
      implicit none
      class(t_flux), intent(in) :: flux
      real(8) :: p_getenthalpy_l
    end function p_getenthalpy_l
    function p_getenthalpy_r(flux)
      import t_flux
      implicit none
      class(t_flux), intent(in) :: flux
      real(8) :: p_getenthalpy_r
    end function p_getenthalpy_r
    function p_getenthalpy_c(flux,h,uv2)
      import t_flux
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: h,uv2
      real(8) :: p_getenthalpy_c
    end function p_getenthalpy_c
  end interface
        
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    subroutine construct(flux,config,grid,variable)
      implicit none
      class(t_flux), intent(out) :: flux
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      
      flux%uref = config%geturef()
      flux%str  = config%getstr()
      flux%pref = config%getpref()
      flux%omega = config%getomega()
            
      select case(config%getprecd())
      case(0)
        flux%getsndp2 => no_prec
      case(1)
        flux%getsndp2 => steady_prec
      case(2)
        flux%getsndp2 => unsteady_prec
      end select

      select case(config%getrotation())
      case(0)
        flux%getenthalpy_l => enthalpy_l
        flux%getenthalpy_r => enthalpy_r
        flux%getenthalpy_c => enthalpy_c
      case(-1,1,-2,2,-3,3)
        flux%getenthalpy_l => rothalpy_l
        flux%getenthalpy_r => rothalpy_r
        flux%getenthalpy_c => rothalpy_c
      case default
      end select

      flux%ngrd = grid%getngrd()
      flux%npv = variable%getnpv()
      flux%ndv = variable%getndv()
      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine destruct(flux)
      implicit none
      class(t_flux), intent(inout) :: flux

      if(associated(flux%nx))            nullify(flux%nx)
      if(associated(flux%grdl))          nullify(flux%grdl)
      if(associated(flux%grdr))          nullify(flux%grdr)
      if(associated(flux%pvl))           nullify(flux%pvl)
      if(associated(flux%pvr))           nullify(flux%pvr)
      if(associated(flux%dvl))           nullify(flux%dvl)
      if(associated(flux%dvr))           nullify(flux%dvr)
      if(associated(flux%sdst))          nullify(flux%sdst)
      if(associated(flux%getsndp2))      nullify(flux%getsndp2)
      if(associated(flux%getenthalpy_l)) nullify(flux%getenthalpy_l)
      if(associated(flux%getenthalpy_r)) nullify(flux%getenthalpy_r)
      if(associated(flux%getenthalpy_c)) nullify(flux%getenthalpy_c)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setnorm(flux,nx)
      implicit none
      class(t_flux), intent(inout) :: flux
      real(8), intent(in), target :: nx(3)
      
      flux%nx => nx

    end subroutine setnorm
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(flux,grdl,grdr)
      implicit none
      class(t_flux), intent(inout) :: flux
      real(8), intent(in), target :: grdl(flux%ngrd),grdr(flux%ngrd)

      flux%grdl => grdl
      flux%grdr => grdr

    end subroutine setgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(flux,pvl,pvr)
      implicit none
      class(t_flux), intent(inout) :: flux
      real(8), intent(in), target :: pvl(flux%npv),pvr(flux%npv)
      
      flux%pvl => pvl
      flux%pvr => pvr
      
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(flux,dvl,dvr)
      implicit none
      class(t_flux), intent(inout) :: flux
      real(8), intent(in), target :: dvl(flux%ndv),dvr(flux%ndv)
      
      flux%dvl => dvl
      flux%dvr => dvr
      
    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setsdst(flux,sdst)
      implicit none
      class(t_flux), intent(inout) :: flux
      real(8), intent(in), target :: sdst(18)
      
      flux%sdst => sdst
      
    end subroutine setsdst
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine roe(flux,eos,fx)
      implicit none
      class(t_roe), intent(in) :: flux
      type(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull,uv2
      real(8) :: ravg(flux%npv),rdv(flux%ndv),ravg_d
      real(8) :: sndp2,sndp2_cut
      real(8) :: uuu,uup,ddd,ddd_cut,c_star,c_star_cut,m_star,du,dp
      real(8) :: df(flux%npv)
      
      dl = dsqrt(flux%nx(1)**2+flux%nx(2)**2+flux%nx(3)**2)
      
      if(dl.eq.0.d0) then
        fx = 0.d0
        return
      end if
      
      nx = flux%nx(1)/dl
      ny = flux%nx(2)/dl
      nz = flux%nx(3)/dl
      
      uull = nx*flux%pvl(2) + ny*flux%pvl(3) + nz*flux%pvl(4)
      uurr = nx*flux%pvr(2) + ny*flux%pvr(3) + nz*flux%pvr(4)

      ! roe average - 1/2 values
      ravg(1) = 0.5d0*(flux%pvr(1)+flux%pvl(1))+flux%pref
      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      do k=2,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do

      call eos%deteos(ravg(1),ravg(5),ravg(6),ravg(7),rdv)
      
      uv2 = ravg(2)**2+ravg(3)**2+ravg(4)**2
      uuu = nx*ravg(2) + ny*ravg(3) + nz*ravg(4)
      
      sndp2     = flux%getsndp2(rdv(6),uv2,0)
      sndp2_cut = flux%getsndp2(rdv(6),uv2,1)
      
      uup = 0.5d0*(1.d0+sndp2/rdv(6))*uuu
      ddd = 0.5d0*dsqrt((1.d0-sndp2/rdv(6))**2*uuu**2+4.d0*sndp2)
      ddd_cut = 0.5d0*dsqrt((1.d0-sndp2_cut/rdv(6))**2*uuu**2+4.d0*sndp2_cut)

      c_star = 0.5d0*(dabs(uup+ddd)+dabs(uup-ddd))
      c_star_cut = 0.5d0*(dabs(uup+ddd_cut)+dabs(uup-ddd_cut))
      m_star = 0.5d0*(dabs(uup+ddd_cut)-dabs(uup-ddd_cut))/ddd_cut

      du = m_star*(uurr-uull)+(c_star_cut-sndp2/rdv(6)*dabs(uuu)-0.5d0*(1.d0-sndp2/rdv(6))*uuu*m_star)*(flux%pvr(1)-flux%pvl(1))/rdv(1)/sndp2_cut
      dp = m_star*(flux%pvr(1)-flux%pvl(1))+(c_star-dabs(uuu)+0.5d0*(1.d0-sndp2/rdv(6))*uuu*m_star)*rdv(1)*(uurr-uull)
      
      df(1) = dabs(uuu)*(flux%dvr(1)             - flux%dvl(1)            ) + du*rdv(1)
      df(2) = dabs(uuu)*(flux%dvr(1)*flux%pvr(2) - flux%dvl(1)*flux%pvl(2)) + du*rdv(1)*ravg(2)  + dp*nx
      df(3) = dabs(uuu)*(flux%dvr(1)*flux%pvr(3) - flux%dvl(1)*flux%pvl(3)) + du*rdv(1)*ravg(3)  + dp*ny
      df(4) = dabs(uuu)*(flux%dvr(1)*flux%pvr(4) - flux%dvl(1)*flux%pvl(4)) + du*rdv(1)*ravg(4)  + dp*nz
      df(5) = dabs(uuu)*(flux%dvr(1)*flux%getenthalpy_r()-flux%pvr(1)     &
                      - (flux%dvl(1)*flux%getenthalpy_l()-flux%pvl(1))    ) + du*rdv(1)*flux%getenthalpy_c(rdv(2),uv2)  + dp*uuu
      do k=6,flux%npv
        df(k) = dabs(uuu)*(flux%dvr(1)*flux%pvr(k) - flux%dvl(1)*flux%pvl(k)) + du*rdv(1)*ravg(k)
      end do
      
      fx(1) = 0.5d0*(flux%dvl(1)*uull + flux%dvr(1)*uurr - df(1))*dl
      fx(2) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(2)+nx*(flux%pvl(1)+flux%pref) + flux%dvr(1)*uurr*flux%pvr(2)+nx*(flux%pvr(1)+flux%pref) - df(2))*dl
      fx(3) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(3)+ny*(flux%pvl(1)+flux%pref) + flux%dvr(1)*uurr*flux%pvr(3)+ny*(flux%pvr(1)+flux%pref) - df(3))*dl
      fx(4) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(4)+nz*(flux%pvl(1)+flux%pref) + flux%dvr(1)*uurr*flux%pvr(4)+nz*(flux%pvr(1)+flux%pref) - df(4))*dl
      fx(5) = 0.5d0*(flux%dvl(1)*uull*flux%getenthalpy_l() + flux%dvr(1)*uurr*flux%getenthalpy_r() - df(5))*dl
      do k=6,flux%npv
        fx(k) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(k) + flux%dvr(1)*uurr*flux%pvr(k) - df(k))*dl
      end do    
    end subroutine roe
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine roem(flux,eos,fx)
      implicit none
      class(t_roem), intent(in) :: flux
      type(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull,uv2
      real(8) :: ravg(flux%npv),ravg_d
      real(8) :: sndp2,sndp2_cut
      real(8) :: uuu,uup,ddd,ddd_cut,c_star,c_star_cut,m_star
      real(8) :: aaa,add,add1,b1,b2,b1b2,rrr,rrr1,ff,gg,sdst(18),pp_l,pp_r
      real(8) :: dqp(flux%npv),fl(flux%npv),fr(flux%npv),bdq(flux%npv),bdq1(flux%npv),dq(flux%npv)
      real(8) :: rdv(flux%ndv)
      
      dl = dsqrt(flux%nx(1)**2+flux%nx(2)**2+flux%nx(3)**2)
      
      if(dl.eq.0.d0) then
        fx = 0.d0
        return
      end if
      
      nx = flux%nx(1)/dl
      ny = flux%nx(2)/dl
      nz = flux%nx(3)/dl
      
      uull = nx*flux%pvl(2) + ny*flux%pvl(3) + nz*flux%pvl(4)
      uurr = nx*flux%pvr(2) + ny*flux%pvr(3) + nz*flux%pvr(4)
    
      ! roe average - 1/2 values
      ravg(1) = 0.5d0*(flux%pvr(1)+flux%pvl(1))+flux%pref
      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      do k=2,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do

      call eos%deteos(ravg(1),ravg(5),ravg(6),ravg(7),rdv)
      
      uv2 = ravg(2)**2+ravg(3)**2+ravg(4)**2
      uuu = nx*ravg(2) + ny*ravg(3) + nz*ravg(4)
      uv2_1 = 0.5d0*(flux%pvl(2)**2 + flux%pvl(3)**2 + flux%pvl(4)**2  &
                   + flux%pvr(2)**2 + flux%pvr(3)**2 + flux%pvr(4)**2 )
      
      sndp2     = flux%getsndp2(rdv(6),uv2,0)
      sndp2_1   = flux%getsndp2(rdv(6),uv2_1,0)
      sndp2_cut = flux%getsndp2(rdv(6),uv2,1)
      
      uup = 0.5d0*(1.d0+sndp2/rdv(6))*uuu
      ddd = 0.5d0*dsqrt((1.d0-sndp2_1/rdv(6))**2*uuu**2+4.d0*sndp2_1)
      ddd_cut = 0.5d0*dsqrt((1.d0-sndp2_cut/rdv(6))**2*uuu**2+4.d0*sndp2_cut)

      c_star = 0.5d0*(dabs(uup+ddd)+dabs(uup-ddd))
      c_star_cut = 0.5d0*(dabs(uup+ddd_cut)+dabs(uup-ddd_cut))
      m_star = 0.5d0*(dabs(uup+ddd_cut)-dabs(uup-ddd_cut))/ddd_cut

      b1 = dmax1(0.d0,uup + ddd_cut,0.5d0*(1.d0+sndp2/rdv(6))*uurr + ddd_cut)
      b2 = dmin1(0.d0,uup - ddd_cut,0.5d0*(1.d0+sndp2/rdv(6))*uull - ddd_cut)
      
      do k=1,flux%npv
        dqp(k) = flux%pvr(k)-flux%pvl(k)
      end do
      
      dq(1) = flux%dvr(1)             - flux%dvl(1) 
      dq(2) = flux%dvr(1)*flux%pvr(2) - flux%dvl(1)*flux%pvl(2)
      dq(3) = flux%dvr(1)*flux%pvr(3) - flux%dvl(1)*flux%pvl(3)
      dq(4) = flux%dvr(1)*flux%pvr(4) - flux%dvl(1)*flux%pvl(4)
      dq(5) = flux%dvr(1)*flux%getenthalpy_r() - flux%dvl(1)*flux%getenthalpy_l()
      do k=6,flux%npv
        dq(k) = flux%dvr(1)*flux%pvr(k) - flux%dvl(1)*flux%pvl(k)
      end do
      
      fl(1) = flux%dvl(1)*uull
      fl(2) = flux%dvl(1)*uull*flux%pvl(2)+nx*(flux%pvl(1)+flux%pref)
      fl(3) = flux%dvl(1)*uull*flux%pvl(3)+ny*(flux%pvl(1)+flux%pref)
      fl(4) = flux%dvl(1)*uull*flux%pvl(4)+nz*(flux%pvl(1)+flux%pref)
      fl(5) = flux%dvl(1)*uull*flux%getenthalpy_l()
      do k=6,flux%npv
        fl(k) = flux%dvl(1)*uull*flux%pvl(k)
      end do
      
      fr(1) = flux%dvr(1)*uurr
      fr(2) = flux%dvr(1)*uurr*flux%pvr(2)+nx*(flux%pvr(1)+flux%pref)
      fr(3) = flux%dvr(1)*uurr*flux%pvr(3)+ny*(flux%pvr(1)+flux%pref)
      fr(4) = flux%dvr(1)*uurr*flux%pvr(4)+nz*(flux%pvr(1)+flux%pref)
      fr(5) = flux%dvr(1)*uurr*flux%getenthalpy_r()
      do k=6,flux%npv
        fr(k) = flux%dvr(1)*uurr*flux%pvr(k)
      end do
      
      do k=1,18
        sdst(k) = flux%sdst(k)+0.5d0*flux%pref+rdv(1)*rdv(6)
      end do
      
      ff = 1.d0 - dmin1(sdst(9)/sdst(10),sdst(10)/sdst(9) &
                       ,sdst(11)/sdst(9),sdst(9)/sdst(11),sdst(9)/sdst(7),sdst(7)/sdst(9)     &
                       ,sdst(9)/sdst(3),sdst(3)/sdst(9),sdst(9)/sdst(15),sdst(15)/sdst(9)     &
                       ,sdst(10)/sdst(12),sdst(12)/sdst(10),sdst(10)/sdst(8),sdst(8)/sdst(10) &
                       ,sdst(4)/sdst(10),sdst(10)/sdst(4),sdst(10)/sdst(16),sdst(16)/sdst(10) )
      
      if(uuu .ne. 0.d0) then
        ff = (dabs(uup)/ddd_cut)**ff
      else
        ff = 1.d0
      end if
      
      pp_l = flux%pvl(1)+flux%pref+0.5d0*rdv(1)*rdv(6)
      pp_r = flux%pvr(1)+flux%pref+0.5d0*rdv(1)*rdv(6)
      gg = 1.d0 - dmin1(pp_l/pp_r,pp_r/pp_l)
      
      if(uuu .ne. 0.d0) then
        gg = (dabs(uup)/ddd_cut)**gg
      else
        gg = 1.d0
      end if
      
      add1 = c_star-dabs(uuu)+0.5d0*(1.d0-sndp2/rdv(6))*uuu*m_star
      add = c_star_cut-dabs(uuu)+0.5d0*(1.d0-sndp2/rdv(6))*uuu*m_star
      aaa = (c_star_cut-sndp2/rdv(6)*dabs(uuu)-0.5d0*(1.d0-sndp2/rdv(6))*uuu*m_star)
      
      if((m_star*uup-c_star_cut).eq.0.d0) then
        rrr = 0.d0
      else
        rrr = -1.d0/(m_star*uup-c_star_cut)
      end if

      if((m_star*uup-c_star).eq.0.d0) then
        rrr1 = 0.d0
      else
        rrr1 = -1.d0/(m_star*uup-c_star)
      end if

      bdq(1) = (add*dq(1)-ff*aaa*dqp(1)/sndp2_cut) 
      bdq(2) = bdq(1)*ravg(2)  + rdv(1)*add*dqp(2)
      bdq(3) = bdq(1)*ravg(3)  + rdv(1)*add*dqp(3)
      bdq(4) = bdq(1)*ravg(4)  + rdv(1)*add*dqp(4)
      bdq(5) = bdq(1)*flux%getenthalpy_c(rdv(2),uv2)  + rdv(1)*add*(flux%getenthalpy_r()-flux%getenthalpy_l())
      do k=6,flux%npv
        bdq(k) = bdq(1)*ravg(k) + rdv(1)*add*dqp(k)
      end do
      bdq1 = 0.d0
      bdq1(2) = -rdv(1)*add1*nx*(uurr-uull)
      bdq1(3) = -rdv(1)*add1*ny*(uurr-uull)
      bdq1(4) = -rdv(1)*add1*nz*(uurr-uull)
      b1b2 = b1*b2 + ddd_cut**2 - ddd*ddd_cut
      do k=1,flux%npv
        fx(k) = (b1*fl(k)-b2*fr(k)+b1*b2*dq(k)-gg*b1*b2*rrr*bdq(k)-gg*b1b2*rrr1*bdq1(k))/(b1-b2)*dl
      end do
    
    end subroutine roem
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
    subroutine ausmpwp(flux,eos,fx)
      implicit none
      class(t_ausmpwp), intent(in) :: flux
      type(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull
      real(8) :: ravg(flux%npv),rdv(flux%ndv),ravg_d
      real(8) :: amid,zml,zmr,am2mid
      real(8) :: am2rmid,am2rmid1,fmid,fmid1,alpha
      real(8) :: zmmr,pmr,zmpl,ppl,pmid,zmid
      real(8) :: ww1,ww2,ww,sdst(18),pp_l,pp_r
      real(8) :: pmt,pwl,pwr,zmpl1,zmmr1
      real(8), parameter :: beta = 0.125d0, ku = 0.25d0
      
      dl = dsqrt(flux%nx(1)**2+flux%nx(2)**2+flux%nx(3)**2)
      
      if(dl.eq.0.d0) then
        fx = 0.d0
        return
      end if
      
      nx = flux%nx(1)/dl
      ny = flux%nx(2)/dl
      nz = flux%nx(3)/dl
      
      uull = nx*flux%pvl(2) + ny*flux%pvl(3) + nz*flux%pvl(4)
      uurr = nx*flux%pvr(2) + ny*flux%pvr(3) + nz*flux%pvr(4)
    
      ! roe average - 1/2 values
      ravg(1) = 0.5d0*(flux%pvr(1)+flux%pvl(1))+flux%pref
      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      do k=2,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do

      call eos%deteos(ravg(1),ravg(5),ravg(6),ravg(7),rdv)
      
      amid = dsqrt(rdv(6))
      
      zmr = uurr/amid
      zml = uull/amid

      am2mid = (ravg(2)**2+ravg(3)**2+ravg(4)**2)/rdv(6)
      am2rmid  = flux%getsndp2(rdv(6),am2mid*rdv(6),1)/rdv(6)
      am2mid = 0.5d0*(flux%pvl(2)**2 + flux%pvl(3)**2 + flux%pvl(4)**2  &
                    + flux%pvr(2)**2 + flux%pvr(3)**2 + flux%pvr(4)**2 )/rdv(6)
      am2rmid1 = flux%getsndp2(rdv(6),am2mid*rdv(6),0)/rdv(6)

            
      fmid = dsqrt(am2rmid)*(2.d0-dsqrt(am2rmid))
      fmid1 = dsqrt(am2rmid1)*(2.d0-dsqrt(am2rmid1))
      
      alpha = 3.d0/16.d0*(-4.d0+5.d0*fmid1**2)
      alpha = dmax1(dmin1(alpha,3.d0/16.d0),-3.d0/4.d0)
      
      if(dabs(zmr).ge.1.d0) then
        zmmr = 0.5d0*(zmr-dabs(zmr))
        pmr = 0.5d0*(1.d0-zmr/dabs(zmr))
      else
        zmmr = -0.25d0*(zmr-1.d0)**2 - beta*(zmr**2-1.d0)**2
        pmr = 0.25d0*(zmr-1.d0)**2*(2.d0+zmr) - alpha*zmr*(zmr**2-1.d0)**2
      end if
      
      if(dabs(zml).ge.1.d0) then
        zmpl = 0.5d0*(zml+dabs(zml))
        ppl = 0.5d0*(1.d0+zml/dabs(zml))
      else
        zmpl = 0.25d0*(zml+1.d0)**2 + beta*(zml**2-1.d0)**2  
        ppl = 0.25d0*(zml+1.d0)**2*(2.d0-zml) + alpha*zml*(zml**2-1.d0)**2
      end if
      
      zmid = zmpl + zmmr
      pmid = ppl*(flux%pvl(1)+flux%pref) + pmr*(flux%pvr(1)+flux%pref) - ku*2.d0*ppl*pmr*rdv(1)*fmid1*amid*(uurr-uull)
      
      do k=1,18
        sdst(k) = flux%sdst(k)+flux%pref+0.5d0*rdv(1)*rdv(6)
      end do
      
      pp_l = flux%pvl(1)+flux%pref+0.5d0*rdv(1)*rdv(6)
      pp_r = flux%pvr(1)+flux%pref+0.5d0*rdv(1)*rdv(6)
      ww1 = 1.d0 - dmin1(pp_l/pp_r,pp_r/pp_l)**3
      ww2 = 1.d0-dmin1(1.d0,dmin1(sdst(7),sdst(8),sdst(11),sdst(12),sdst(3),sdst(4),sdst(15),sdst(16))/dmax1(sdst(7),sdst(8),sdst(11),sdst(12),sdst(3),sdst(4),sdst(15),sdst(16)))**2
      !ww2 = (1.d0-dmin1(1.d0,4.d0*(sdst(4)-sdst(3))/(sdst(6)+sdst(5)-sdst(1)-sdst(2)+1.d-12)))**2*(1.d0-dmin1(sdst(3)/sdst(4),sdst(4)/sdst(3)))**2
      ww = dmax1(ww1,ww2)
      pmt = pmid + 0.5d0*rdv(1)*rdv(6)
      
      if(pmt.ne.0.d0) then
        pwl = ((flux%pvl(1)+flux%pref+0.5d0*rdv(1)*rdv(6))/pmt-1.d0)*(1.d0-ww2)/fmid
        pwr = ((flux%pvr(1)+flux%pref+0.5d0*rdv(1)*rdv(6))/pmt-1.d0)*(1.d0-ww2)/fmid
      else
        pwl = 0.d0
        pwr = 0.d0
      end if
      
      if( zmid .ge. 0.d0) then
        zmpl1 = zmpl + zmmr*((1.d0-ww)*(1.d0+pwr)-pwl)
        zmmr1 = zmmr*ww*(1.d0+pwr)
      else
        zmpl1 = zmpl*ww*(1.d0+pwl)
        zmmr1 = zmmr + zmpl*((1.d0-ww)*(1.d0+pwl)-pwr)
      end if
      
      fx(1) = (zmpl1*amid*flux%dvl(1)              + zmmr1*amid*flux%dvr(1)                      )*dl
      fx(2) = (zmpl1*amid*flux%dvl(1)*flux%pvl(2)  + zmmr1*amid*flux%dvr(1)*flux%pvr(2) + pmid*nx)*dl
      fx(3) = (zmpl1*amid*flux%dvl(1)*flux%pvl(3)  + zmmr1*amid*flux%dvr(1)*flux%pvr(3) + pmid*ny)*dl
      fx(4) = (zmpl1*amid*flux%dvl(1)*flux%pvl(4)  + zmmr1*amid*flux%dvr(1)*flux%pvr(4) + pmid*nz)*dl
      fx(5) = (zmpl1*amid*flux%dvl(1)*flux%getenthalpy_l() + zmmr1*amid*flux%dvr(1)*flux%getenthalpy_r())*dl
      do k=6,flux%npv
        fx(k) = (zmpl1*amid*flux%dvl(1)*flux%pvl(k)+ zmmr1*amid*flux%dvr(1)*flux%pvr(k))*dl
      end do
        
    end subroutine ausmpwp
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine ausmpup(flux,eos,fx)
      implicit none
      class(t_ausmpup), intent(in) :: flux
      type(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull
      real(8) :: ravg(flux%npv),rdv(flux%ndv),ravg_d
      real(8) :: amid,zml,zmr,am2mid
      real(8) :: am2rmid,am2rmid1,fmid,fmid1,alpha
      real(8) :: zmmr,pmr,zmpl,ppl,pmid,zmid
      real(8) :: zmpl1,zmmr1  
      real(8), parameter :: beta = 0.125d0, kp = 0.25d0, ku = 0.25d0
      

      dl = dsqrt(flux%nx(1)**2+flux%nx(2)**2+flux%nx(3)**2)
      
      if(dl.eq.0.d0) then
        fx = 0.d0
        return
      end if
      
      nx = flux%nx(1)/dl
      ny = flux%nx(2)/dl
      nz = flux%nx(3)/dl
      
      uull = nx*flux%pvl(2) + ny*flux%pvl(3) + nz*flux%pvl(4)
      uurr = nx*flux%pvr(2) + ny*flux%pvr(3) + nz*flux%pvr(4)
    
      ! roe average - 1/2 values
      ravg(1) = 0.5d0*(flux%pvr(1)+flux%pvl(1))+flux%pref
      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      do k=2,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do

      call eos%deteos(ravg(1),ravg(5),ravg(6),ravg(7),rdv)
      
      amid = dsqrt(rdv(6))
      
      zmr = uurr/amid
      zml = uull/amid

      am2mid = (ravg(2)**2+ravg(3)**2+ravg(4)**2)/rdv(6)
      am2rmid  = flux%getsndp2(rdv(6),am2mid*rdv(6),1)/rdv(6)
      am2mid = 0.5d0*(flux%pvl(2)**2 + flux%pvl(3)**2 + flux%pvl(4)**2  &
                    + flux%pvr(2)**2 + flux%pvr(3)**2 + flux%pvr(4)**2 )/rdv(6)
      am2rmid1 = flux%getsndp2(rdv(6),am2mid*rdv(6),0)/rdv(6)

            
      fmid = dsqrt(am2rmid)*(2.d0-dsqrt(am2rmid))
      fmid1 = dsqrt(am2rmid1)*(2.d0-dsqrt(am2rmid1))
      
      alpha = 3.d0/16.d0*(-4.d0+5.d0*fmid1**2)
      alpha = dmax1(dmin1(alpha,3.d0/16.d0),-3.d0/4.d0)
      
      if(dabs(zmr).ge.1.d0) then
        zmmr = 0.5d0*(zmr-dabs(zmr))
        pmr = 0.5d0*(1.d0-zmr/dabs(zmr))
      else
        zmmr = -0.25d0*(zmr-1.d0)**2 - beta*(zmr**2-1.d0)**2
        pmr = 0.25d0*(zmr-1.d0)**2*(2.d0+zmr) - alpha*zmr*(zmr**2-1.d0)**2
      end if
      
      if(dabs(zml).ge.1.d0) then
        zmpl = 0.5d0*(zml+dabs(zml))
        ppl = 0.5d0*(1.d0+zml/dabs(zml))
      else
        zmpl = 0.25d0*(zml+1.d0)**2 + beta*(zml**2-1.d0)**2  
        ppl = 0.25d0*(zml+1.d0)**2*(2.d0-zml) + alpha*zml*(zml**2-1.d0)**2
      end if
      
      zmid = zmpl + zmmr - kp*dmax1(1.d0-am2mid,0.d0)*(flux%pvr(1)-flux%pvl(1))/rdv(1)/rdv(6)/fmid
      pmid = ppl*(flux%pvl(1)+flux%pref) + pmr*(flux%pvr(1)+flux%pref) - ku*2.d0*ppl*pmr*rdv(1)*fmid1*amid*(uurr-uull)
      
      if(zmid.gt.0.d0) then
        zmpl1 = zmid 
        zmmr1 = 0.d0
      else
        zmmr1 = zmid
        zmpl1 = 0.d0
      end if
      
      fx(1) = (zmpl1*amid*flux%dvl(1)              + zmmr1*amid*flux%dvr(1)                      )*dl
      fx(2) = (zmpl1*amid*flux%dvl(1)*flux%pvl(2)  + zmmr1*amid*flux%dvr(1)*flux%pvr(2) + pmid*nx)*dl
      fx(3) = (zmpl1*amid*flux%dvl(1)*flux%pvl(3)  + zmmr1*amid*flux%dvr(1)*flux%pvr(3) + pmid*ny)*dl
      fx(4) = (zmpl1*amid*flux%dvl(1)*flux%pvl(4)  + zmmr1*amid*flux%dvr(1)*flux%pvr(4) + pmid*nz)*dl
      fx(5) = (zmpl1*amid*flux%dvl(1)*flux%getenthalpy_l() + zmmr1*amid*flux%dvr(1)*flux%getenthalpy_r())*dl
      do k=6,flux%npv
        fx(k) = (zmpl1*amid*flux%dvl(1)*flux%pvl(k)+ zmmr1*amid*flux%dvr(1)*flux%pvr(k))*dl
      end do
     
    end subroutine ausmpup
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function no_prec(flux,snd2,uuu2,cut) result(sndp2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      integer, intent(in) :: cut
      real(8) :: sndp2
      
      sndp2 = snd2
      
    end function no_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function steady_prec(flux,snd2,uuu2,cut) result(sndp2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      integer, intent(in) :: cut
      real(8) :: sndp2  
      
      sndp2 = dmin1(snd2,dmax1(uuu2,dble(cut)*flux%uref**2))
      
    end function steady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteady_prec(flux,snd2,uuu2,cut) result(sndp2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      integer, intent(in) :: cut
      real(8) :: sndp2  
      
      sndp2 = dmin1(snd2,dmax1(uuu2,dble(cut)*flux%uref**2,dble(cut)*flux%str**2*uuu2))
      
    end function unsteady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function enthalpy_l(flux)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8) :: enthalpy_l
      enthalpy_l = flux%dvl(2) + 0.5d0*(flux%pvl(2)**2 + flux%pvl(3)**2 + flux%pvl(4)**2)
    end function enthalpy_l
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function enthalpy_r(flux)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8) :: enthalpy_r
      enthalpy_r = flux%dvr(2) + 0.5d0*(flux%pvr(2)**2 + flux%pvr(3)**2 + flux%pvr(4)**2)
    end function enthalpy_r
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function enthalpy_c(flux,h,uv2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: h,uv2
      real(8) :: enthalpy_c
      enthalpy_c = h+0.5d0*uv2
    end function enthalpy_c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rothalpy_l(flux)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8) :: rothalpy_l
      rothalpy_l = flux%dvl(2) + 0.5d0*(flux%pvl(2)**2 + flux%pvl(3)**2 + flux%pvl(4)**2)
      rothalpy_l = rothalpy_l -0.5d0*(flux%omega(1)**2*(flux%grdl(3)**2+flux%grdl(4)**2) &
                                     +flux%omega(2)**2*(flux%grdl(2)**2+flux%grdl(4)**2) &
                                     +flux%omega(3)**2*(flux%grdl(2)**2+flux%grdl(3)**2) &
                               -2.d0*(flux%omega(1)*flux%omega(2)*flux%grdl(2)*flux%grdl(3) &
                                    + flux%omega(2)*flux%omega(3)*flux%grdl(3)*flux%grdl(4) &
                                    + flux%omega(3)*flux%omega(1)*flux%grdl(4)*flux%grdl(2)))
    end function rothalpy_l
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rothalpy_r(flux)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8) :: rothalpy_r
      rothalpy_r = flux%dvr(2) + 0.5d0*(flux%pvr(2)**2 + flux%pvr(3)**2 + flux%pvr(4)**2)
      rothalpy_r = rothalpy_r -0.5d0*(flux%omega(1)**2*(flux%grdr(3)**2+flux%grdr(4)**2) &
                                     +flux%omega(2)**2*(flux%grdr(2)**2+flux%grdr(4)**2) &
                                     +flux%omega(3)**2*(flux%grdr(2)**2+flux%grdr(3)**2) &
                               -2.d0*(flux%omega(1)*flux%omega(2)*flux%grdr(2)*flux%grdr(3) &
                                    + flux%omega(2)*flux%omega(3)*flux%grdr(3)*flux%grdr(4) &
                                    + flux%omega(3)*flux%omega(1)*flux%grdr(4)*flux%grdr(2)))
    end function rothalpy_r
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function rothalpy_c(flux,h,uv2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: h,uv2
      real(8) :: rothalpy_c
      real(8) :: x(3)
      x(1) = 0.5d0*(flux%grdl(2)+flux%grdr(2))
      x(2) = 0.5d0*(flux%grdl(3)+flux%grdr(3))
      x(3) = 0.5d0*(flux%grdl(4)+flux%grdr(4))

      rothalpy_c = h+0.5d0*uv2
      rothalpy_c = rothalpy_c -0.5d0*(flux%omega(1)**2*(x(2)**2+x(3)**2) &
                                     +flux%omega(2)**2*(x(1)**2+x(3)**2) &
                                     +flux%omega(3)**2*(x(1)**2+x(2)**2) &
                               -2.d0*(flux%omega(1)*flux%omega(2)*x(1)*x(2) &
                                    + flux%omega(2)*flux%omega(3)*x(2)*x(3) &
                                    + flux%omega(3)*flux%omega(1)*x(3)*x(1)))
    end function rothalpy_c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module flux_module
