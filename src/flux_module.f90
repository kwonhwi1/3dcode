module flux_module
  use config_module
  use grid_module
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
    procedure(p_getsndp2), pointer :: getsndp2_du
    procedure(p_getenthalpy_l), pointer :: getenthalpy_l
    procedure(p_getenthalpy_r), pointer :: getenthalpy_r
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
      class(t_eos), intent(in) :: eos
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
    function p_getsndp2(flux,snd2,uuu2) result(sndp2)
      import t_flux
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
    end function p_getsndp2
    function p_getsndp2_du(flux,snd2,uuu2) result(sndp2)
      import t_flux
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
    end function p_getsndp2_du
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
  end interface
        
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    subroutine construct(flux,config,grid)
      implicit none
      class(t_flux), intent(out) :: flux
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid

      flux%npv = config%getnpv()
      flux%ndv = config%getndv()
      flux%uref = config%geturef()
      flux%str  = config%getstr()
      flux%pref = config%getpref()
      flux%omega = config%getomega()
            
      select case(config%getprec())
      case(0)
        flux%getsndp2 => no_prec
        flux%getsndp2_du => no_prec
      case(1)
        flux%getsndp2 => steady_prec
        flux%getsndp2_du => no_cut
      case(2)
        flux%getsndp2 => unsteady_prec
        flux%getsndp2_du => no_cut
      end select

      select case(config%getrotation())
      case(0)
        flux%getenthalpy_l => enthalpy_l
        flux%getenthalpy_r => enthalpy_r
      case(-1,1,-2,2,-3,3)
        flux%getenthalpy_l => rothalpy_l
        flux%getenthalpy_r => rothalpy_r
      case default
      end select

      flux%ngrd = grid%getngrd()

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
      if(associated(flux%getsndp2_du)) nullify(flux%getsndp2_du)
      if(associated(flux%getenthalpy_l)) nullify(flux%getenthalpy_l)
      if(associated(flux%getenthalpy_r)) nullify(flux%getenthalpy_r)
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
      class(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull
      real(8) :: ravg(flux%npv),rdv(flux%ndv),ravg_d,ravg_ht
      real(8) :: sndp2,u2
      real(8) :: uuu,uup,ddd,c_star,m_star,du,dp
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
      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      ravg(1) = 0.5d0*(flux%pvl(1)+flux%pvr(1))
      ravg(5) = 0.5d0*(flux%pvl(5)+flux%pvr(5))
      do k=2,4
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do
      do k=6,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do
      ravg_ht = (dsqrt(flux%dvl(1))*flux%getenthalpy_l()+dsqrt(flux%dvr(1))*flux%getenthalpy_r())*ravg_d

      call eos%deteos_simple(ravg(1)+flux%pref,ravg(5),ravg(6),ravg(7),rdv)
      u2 = ravg(2)**2+ravg(3)**2+ravg(4)**2
      rdv(1) = dsqrt(flux%dvl(1)*flux%dvr(1))

      uuu = nx*ravg(2) + ny*ravg(3) + nz*ravg(4)

      sndp2 = flux%getsndp2(rdv(6),u2)
      
      uup = 0.5d0*(1.d0+sndp2/rdv(6))*uuu
      ddd = 0.5d0*dsqrt((1.d0-sndp2/rdv(6))**2*uuu**2+4.d0*sndp2)

      c_star = 0.5d0*(dabs(uup+ddd)+dabs(uup-ddd))
      m_star = 0.5d0*(dabs(uup+ddd)-dabs(uup-ddd))/ddd

      du = m_star*(uurr-uull)+(c_star-sndp2/rdv(6)*dabs(uuu)-0.5d0*(1.d0-sndp2/rdv(6))*uuu*m_star)*(flux%pvr(1)-flux%pvl(1))/rdv(1)/sndp2
      dp = m_star*(flux%pvr(1)-flux%pvl(1))+(c_star-dabs(uuu)+0.5d0*(1.d0-sndp2/rdv(6))*uuu*m_star)*rdv(1)*(uurr-uull)
      
      df(1) = dabs(uuu)*(flux%dvr(1)             - flux%dvl(1)            ) + du*rdv(1)
      df(2) = dabs(uuu)*(flux%dvr(1)*flux%pvr(2) - flux%dvl(1)*flux%pvl(2)) + du*rdv(1)*ravg(2)  + dp*nx
      df(3) = dabs(uuu)*(flux%dvr(1)*flux%pvr(3) - flux%dvl(1)*flux%pvl(3)) + du*rdv(1)*ravg(3)  + dp*ny
      df(4) = dabs(uuu)*(flux%dvr(1)*flux%pvr(4) - flux%dvl(1)*flux%pvl(4)) + du*rdv(1)*ravg(4)  + dp*nz
      df(5) = dabs(uuu)*(flux%dvr(1)*flux%getenthalpy_r() - flux%pvr(1)   &
                       - flux%dvl(1)*flux%getenthalpy_l() + flux%pvl(1)   ) + du*rdv(1)*ravg_ht  + dp*uuu
      do k=6,flux%npv
        df(k) = dabs(uuu)*(flux%dvr(1)*flux%pvr(k) - flux%dvl(1)*flux%pvl(k)) + du*rdv(1)*ravg(k)
      end do
      
      fx(1) = 0.5d0*(flux%dvl(1)*uull + flux%dvr(1)*uurr                     - df(1))*dl
      fx(2) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(2)+nx*(flux%pvl(1)+flux%pref) &
                   + flux%dvr(1)*uurr*flux%pvr(2)+nx*(flux%pvr(1)+flux%pref) - df(2))*dl
      fx(3) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(3)+ny*(flux%pvl(1)+flux%pref) &
                   + flux%dvr(1)*uurr*flux%pvr(3)+ny*(flux%pvr(1)+flux%pref) - df(3))*dl
      fx(4) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(4)+nz*(flux%pvl(1)+flux%pref) &
                   + flux%dvr(1)*uurr*flux%pvr(4)+nz*(flux%pvr(1)+flux%pref) - df(4))*dl
      fx(5) = 0.5d0*(flux%dvl(1)*uull*flux%getenthalpy_l()                   &
                   + flux%dvr(1)*uurr*flux%getenthalpy_r()                   - df(5))*dl
      do k=6,flux%npv
        fx(k) = 0.5d0*(flux%dvl(1)*uull*flux%pvl(k) + flux%dvr(1)*uurr*flux%pvr(k) - df(k))*dl
      end do    
    end subroutine roe
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine roem(flux,eos,fx)
      implicit none
      class(t_roem), intent(in) :: flux
      class(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull
      real(8) :: ravg(flux%npv),ravg_d,ravg_ht
      real(8) :: am2rmid_du,am2rmid,fmid,fmid_du,amid
      real(8) :: uuu,c_star,c_star_du,m_star,u2
      real(8) :: aaa,add,b1,b2,b1_du,b2_du,ff,gg,sdst(18),pp_l,pp_r,rhom
      real(8) :: dqp(flux%npv),fl(flux%npv),fr(flux%npv),bdq(flux%npv),dq(flux%npv)
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
      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      ravg(1) = 0.5d0*(flux%pvl(1)+flux%pvr(1))
      ravg(5) = 0.5d0*(flux%pvl(5)+flux%pvr(5))
      do k=2,4
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do
      do k=6,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do
      ravg_ht = (dsqrt(flux%dvl(1))*flux%getenthalpy_l()+dsqrt(flux%dvr(1))*flux%getenthalpy_r())*ravg_d

      call eos%deteos_simple(ravg(1)+flux%pref,ravg(5),ravg(6),ravg(7),rdv)
      u2 = ravg(2)**2+ravg(3)**2+ravg(4)**2
      rdv(1) = dsqrt(flux%dvl(1)*flux%dvr(1))


      uuu = nx*ravg(2) + ny*ravg(3) + nz*ravg(4)

      am2rmid_du = flux%getsndp2_du(rdv(6),u2)/rdv(6)
      am2rmid  = flux%getsndp2(rdv(6),u2)/rdv(6)

      fmid = dsqrt(am2rmid)*(2.d0-dsqrt(am2rmid))
      fmid_du = dsqrt(am2rmid_du)*(2.d0-dsqrt(am2rmid_du))

      amid = dsqrt(rdv(6))

      b1 = dmax1(uuu+amid,uurr+amid)
      b2 = dmin1(uuu-amid,uull-amid)
      b1_du = dmax1(uuu+fmid_du*amid,uurr+fmid_du*amid)
      b2_du = dmin1(uuu-fmid_du*amid,uull-fmid_du*amid)

      c_star = 0.5d0*(dabs(b1)+dabs(b2))
      m_star = 0.5d0*(dabs(b1)-dabs(b2))/amid
      c_star_du = 0.5d0*(dabs(b1_du)+dabs(b2_du))
      
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
      
      rhom = dmin1(flux%dvl(1),flux%dvr(1))

      do k=1,18
        sdst(k) = flux%sdst(k)+0.5d0*flux%pref+rhom*rdv(6)
      end do
      
      ff = 1.d0 - dmin1(sdst(9)/sdst(10),sdst(10)/sdst(9) &
                       ,sdst(11)/sdst(9),sdst(9)/sdst(11),sdst(9)/sdst(7),sdst(7)/sdst(9)     &
                       ,sdst(9)/sdst(3),sdst(3)/sdst(9),sdst(9)/sdst(15),sdst(15)/sdst(9)     &
                       ,sdst(10)/sdst(12),sdst(12)/sdst(10),sdst(10)/sdst(8),sdst(8)/sdst(10) &
                       ,sdst(4)/sdst(10),sdst(10)/sdst(4),sdst(10)/sdst(16),sdst(16)/sdst(10) )
      
      if(uuu .ne. 0.d0) then
        ff = (dabs(uuu)/amid)**ff
      else
        ff = 1.d0
      end if
      
      pp_l = flux%pvl(1)+flux%pref+rhom*rdv(6)
      pp_r = flux%pvr(1)+flux%pref+rhom*rdv(6)
      gg = 1.d0 - dmin1(pp_l/pp_r,pp_r/pp_l)
      
      if(uuu .ne. 0.d0) then
        gg = (dabs(uuu)/amid)**gg
      else
        gg = 1.d0
      end if

      add = c_star-dabs(uuu)
      aaa = c_star_du-dabs(uuu)

      bdq(1) = add*(dq(1)-ff*dqp(1)/rdv(6)/fmid)
      bdq(2) = bdq(1)*ravg(2)  + rdv(1)*(add*dqp(2)-aaa*nx*(uurr-uull))
      bdq(3) = bdq(1)*ravg(3)  + rdv(1)*(add*dqp(3)-aaa*ny*(uurr-uull))
      bdq(4) = bdq(1)*ravg(4)  + rdv(1)*(add*dqp(4)-aaa*nz*(uurr-uull))
      bdq(5) = bdq(1)*ravg_ht  + rdv(1)*add*(flux%getenthalpy_r()-flux%getenthalpy_l())
      do k=6,flux%npv
        bdq(k) = bdq(1)*ravg(k) + rdv(1)*add*dqp(k)
      end do

      do k=1,flux%npv
        fx(k) = 0.5d0*(fl(k)+fr(k)-m_star*(fr(k)-fl(k))+(m_star*uuu-c_star)*dq(k)+gg*bdq(k))*dl
      end do
    
    end subroutine roem
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
    subroutine ausmpwp(flux,eos,fx)
      implicit none
      class(t_ausmpwp), intent(in) :: flux
      class(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull
      real(8) :: ravg(flux%npv),rdv(flux%ndv),ravg_d
      real(8) :: amid,zml,zmr,am2mid,rhom
      real(8) :: am2rmid,am2rmid1,fmid,fmid1,alpha
      real(8) :: zmmr,pmr,zmpl,ppl,pmid,zmid
      real(8) :: ww1,ww2,ww,sdst(18),pp_l,pp_r
      real(8) :: pmt,pwl,pwr,zmpl1,zmmr1
      real(8), parameter :: beta = 0.125d0,ku=1.d0
      
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

      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      ravg(1) = 0.5d0*(flux%pvl(1)+flux%pvr(1))
      ravg(5) = 0.5d0*(flux%pvl(5)+flux%pvr(5))
      do k=2,4
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do
      do k=6,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do

      call eos%deteos_simple(ravg(1)+flux%pref,ravg(5),ravg(6),ravg(7),rdv)

      rdv(1) = dsqrt(flux%dvl(1)*flux%dvr(1))

      amid = dsqrt(rdv(6))
      
      zmr = uurr/amid
      zml = uull/amid

      am2mid = ravg(2)**2+ravg(3)**2+ravg(4)**2
      am2rmid1 = flux%getsndp2_du(rdv(6),am2mid)/rdv(6)
      am2rmid  = flux%getsndp2(rdv(6),am2mid)/rdv(6)

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
      
      rhom = dmin1(flux%dvl(1),flux%dvr(1))
      do k=1,18
        sdst(k) = flux%sdst(k)+flux%pref+rhom*rdv(6)
      end do
      
      pp_l = flux%pvl(1)+flux%pref+rhom*rdv(6)
      pp_r = flux%pvr(1)+flux%pref+rhom*rdv(6)
      ww1 = 1.d0 - dmin1(pp_l/pp_r,pp_r/pp_l)**3
      ww2 = 1.d0-dmin1(1.d0,dmin1(sdst(7),sdst(8),sdst(11),sdst(12),sdst(3),sdst(4),sdst(15),sdst(16))/dmax1(sdst(7),sdst(8),sdst(11),sdst(12),sdst(3),sdst(4),sdst(15),sdst(16)))**2
      !ww2 = (1.d0-dmin1(1.d0,4.d0*(sdst(4)-sdst(3))/(sdst(6)+sdst(5)-sdst(1)-sdst(2)+1.d-12)))**2*(1.d0-dmin1(sdst(3)/sdst(4),sdst(4)/sdst(3)))**2
      ww = dmax1(ww1,ww2)
      pmt = pmid + rhom*rdv(6)
      
      if(pmt.ne.0.d0) then
        pwl = ((flux%pvl(1)+flux%pref+rhom*rdv(6))/pmt-1.d0)*(1.d0-ww2)/fmid
        pwr = ((flux%pvr(1)+flux%pref+rhom*rdv(6))/pmt-1.d0)*(1.d0-ww2)/fmid
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
      class(t_eos), intent(in) :: eos
      real(8), intent(out) :: fx(flux%npv)
      integer :: k
      real(8) :: nx,ny,nz,dl
      real(8) :: uurr,uull
      real(8) :: ravg(flux%npv),rdv(flux%ndv),ravg_d
      real(8) :: amid,zml,zmr,am2mid
      real(8) :: am2rmid,fmid,alpha
      real(8) :: zmmr,pmr,zmpl,ppl,pmid,zmid
      real(8) :: zmpl1,zmmr1  
      real(8), parameter :: beta = 0.125d0, kp = 0.25d0, ku = 0.5d0
      

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

      ravg_d = 1.d0/(dsqrt(flux%dvl(1))+dsqrt(flux%dvr(1)))
      ravg(1) = 0.5d0*(flux%pvl(1)+flux%pvr(1))
      ravg(5) = 0.5d0*(flux%pvl(5)+flux%pvr(5))
      do k=2,4
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do
      do k=6,flux%npv
        ravg(k) = (dsqrt(flux%dvl(1))*flux%pvl(k)+dsqrt(flux%dvr(1))*flux%pvr(k))*ravg_d
      end do

      call eos%deteos_simple(ravg(1)+flux%pref,ravg(5),ravg(6),ravg(7),rdv)

      rdv(1) = dsqrt(flux%dvl(1)*flux%dvr(1))
      
      amid = dsqrt(rdv(6))
      
      zmr = uurr/amid
      zml = uull/amid

      am2mid =  (nx*ravg(2) + ny*ravg(3) + nz*ravg(4))**2
      am2rmid = flux%getsndp2(rdv(6),am2mid)/rdv(6)
            
      fmid = dsqrt(am2rmid)*(2.d0-dsqrt(am2rmid))
      
      alpha = 3.d0/16.d0*(-4.d0+5.d0*fmid**2)
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
      
      zmid = zmpl + zmmr - kp*dmax1(1.d0-am2mid/rdv(6),0.d0)*(flux%pvr(1)-flux%pvl(1))/rdv(1)/rdv(6)/fmid
      pmid = ppl*(flux%pvl(1)+flux%pref) + pmr*(flux%pvr(1)+flux%pref) - ku*2.d0*ppl*pmr*rdv(1)*fmid*amid*(uurr-uull)
      
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
    function no_prec(flux,snd2,uuu2) result(sndp2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
      
      sndp2 = snd2
      
    end function no_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function steady_prec(flux,snd2,uuu2) result(sndp2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2  
      
      sndp2 = dmin1(snd2,dmax1(uuu2,flux%uref**2))
      
    end function steady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteady_prec(flux,snd2,uuu2) result(sndp2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2  
      
      sndp2 = dmin1(snd2,dmax1(uuu2,flux%uref**2,flux%str**2))
      
    end function unsteady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function no_cut(flux,snd2,uuu2) result(sndp2)
      implicit none
      class(t_flux), intent(in) :: flux
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2  
      
      sndp2 = dmin1(snd2,uuu2)
      
    end function no_cut
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
end module flux_module
