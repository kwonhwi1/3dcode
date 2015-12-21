module timestep_module
  use mpi
  use config_module
  use grid_module
  use variable_module
  implicit none
  private
  public :: t_timestep,t_localtime,t_mintime,t_fixedtime

  type, abstract :: t_timestep
    private
    integer :: rank,nprint
    integer :: npv,ndv,ntv,ngrd,imax,jmax,kmax
    real(8) :: cfl,uref,str
    real(8), dimension(:,:,:), allocatable :: dt
    procedure(p_getsndp2), pointer :: getsndp2
    procedure(p_geteigenvis), pointer :: geteigenvis
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: getdt
      procedure(p_caltimestep), deferred :: caltimestep
  end type t_timestep

  abstract interface
    subroutine p_caltimestep(timestep,grid,variable,nt_phy,nt)
      import t_timestep
      import t_grid
      import t_variable
      implicit none
      class(t_timestep), intent(inout) :: timestep
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      integer, intent(in) :: nt_phy,nt
      
    end subroutine p_caltimestep
  end interface
  
  type, extends(t_timestep) :: t_localtime
    contains
      procedure :: caltimestep => localtime
  end type t_localtime
  
  type, extends(t_timestep) :: t_mintime
    contains
      procedure :: caltimestep => mintime
  end type t_mintime
  
  type, extends(t_timestep) :: t_fixedtime
    contains
      procedure :: caltimestep => fixedtime
  end type t_fixedtime
  
  interface
    function p_getsndp2(timestep,snd2,uuu2) result(sndp2)
      import t_timestep
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
    end function p_getsndp2
    function p_geteigenvis(timestep,dv,tv) result(eigenvis)
      import t_timestep
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: dv(timestep%ndv),tv(timestep%ntv)
      real(8) :: eigenvis
    end function p_geteigenvis
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(timestep,config,grid)
      implicit none
      class(t_timestep), intent(out) :: timestep
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid

      timestep%npv = config%getnpv()
      timestep%ndv = config%getndv()
      timestep%ntv = config%getntv()
      timestep%rank = config%getrank()
      timestep%nprint = config%getnprint()
      timestep%cfl  = config%getcfl()
      timestep%uref = config%geturef()
      timestep%str  =  config%getstr()
      
      select case(config%getiturb())
      case(-1,0)
        timestep%geteigenvis => turbulent
      case(-2)
        timestep%geteigenvis => laminar
      case(-3)
        timestep%geteigenvis => euler
      end select

      select case(config%getprec())
      case(0)
        timestep%getsndp2 => no_prec
      case(1)
        timestep%getsndp2 => steady_prec
      case(2)
        timestep%getsndp2 => unsteady_prec
      end select

      timestep%ngrd = grid%getngrd()
      timestep%imax = grid%getimax()
      timestep%jmax = grid%getjmax()
      timestep%kmax = grid%getkmax()
      
      allocate(timestep%dt(2:timestep%imax,2:timestep%jmax,2:timestep%kmax))
      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(timestep)
      implicit none
      class(t_timestep), intent(inout) :: timestep
      
      if(associated(timestep%getsndp2))     nullify(timestep%getsndp2)
      if(associated(timestep%geteigenvis))  nullify(timestep%geteigenvis)
      
      deallocate(timestep%dt)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine localtime(timestep,grid,variable,nt_phy,nt)
      implicit none
      class(t_localtime), intent(inout) :: timestep
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      integer, intent(in) :: nt_phy,nt
      integer :: i,j,k    
      real(8) :: cx1(3),cx2(3),ex1(3),ex2(3),tx1(3),tx2(3)
      real(8) :: pv(timestep%npv)
      real(8) :: tv(timestep%ntv)
      real(8) :: dv(timestep%ndv)
      real(8) :: grd(timestep%ngrd)
      real(8) :: c1,c2,c3,e1,e2,e3,t1,t2,t3,s1,s2,s3
      real(8) :: uc,vc,wc,uv2,sndp2
      real(8) :: up,d,eigenx,eigeny,eigenz,eigenvis
      
      do k = 2,timestep%kmax
        do j = 2,timestep%jmax 
          do i = 2,timestep%imax            
            pv = variable%getpv(i,j,k)
            tv = variable%gettv(i,j,k)
            dv = variable%getdv(i,j,k)
            grd = grid%getgrd(i,j,k)
            cx1 = grid%getcx(i-1,j,k)
            cx2 = grid%getcx(i,j,k)
            ex1 = grid%getex(i,j-1,k)
            ex2 = grid%getex(i,j,k)
            tx1 = grid%gettx(i,j,k-1)
            tx2 = grid%gettx(i,j,k)
        
            c1 = 0.5d0*(cx1(1)+cx2(1))
            c2 = 0.5d0*(cx1(2)+cx2(2))
            c3 = 0.5d0*(cx1(3)+cx2(3))
            e1 = 0.5d0*(ex1(1)+ex2(1))
            e2 = 0.5d0*(ex1(2)+ex2(2))
            e3 = 0.5d0*(ex1(3)+ex2(3))
            t1 = 0.5d0*(tx1(1)+tx2(1))
            t2 = 0.5d0*(tx1(2)+tx2(2))
            t3 = 0.5d0*(tx1(3)+tx2(3))
            
            uc = dabs(c1*pv(2)+c2*pv(3)+c3*pv(4))
            vc = dabs(e1*pv(2)+e2*pv(3)+e3*pv(4))
            wc = dabs(t1*pv(2)+t2*pv(3)+t3*pv(4))
            
            s1 = c1**2+c2**2+c3**2
            s2 = e1**2+e2**2+e3**2
            s3 = t1**2+t2**2+t3**2
            
            uv2 = pv(2)**2+pv(3)**2+pv(4)**2
            sndp2 = timestep%getsndp2(dv(6),uv2)
            
            up = uc*(1.d0+sndp2/dv(6))
            d = dsqrt(uc**2*(1.d0-sndp2/dv(6))**2+4.d0*sndp2*s1)
            eigenx = 0.5d0*(up+d)

            up = vc*(1.d0+sndp2/dv(6))
            d = dsqrt(vc**2*(1.d0-sndp2/dv(6))**2+4.d0*sndp2*s2)
            eigeny = 0.5d0*(up+d)

            up = wc*(1.d0+sndp2/dv(6))
            d = dsqrt(wc**2*(1.d0-sndp2/dv(6))**2+4.d0*sndp2*s3)
            eigenz = 0.5d0*(up+d)
            
            eigenvis = timestep%geteigenvis(dv,tv)*(s1+s2+s3)/grd(1)
            
            timestep%dt(i,j,k) = timestep%cfl*grd(1)/(eigenx + eigeny + eigenz + 4.d0*eigenvis)
           end do
        end do
      end do
    end subroutine localtime
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine mintime(timestep,grid,variable,nt_phy,nt)
      implicit none
      class(t_mintime), intent(inout) :: timestep
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      integer, intent(in) :: nt_phy,nt
      integer :: i,j,k    
      real(8) :: cx1(3),cx2(3),ex1(3),ex2(3),tx1(3),tx2(3)
      real(8) :: pv(timestep%npv)
      real(8) :: tv(timestep%ntv)
      real(8) :: dv(timestep%ndv)
      real(8) :: grd(timestep%ngrd)
      real(8) :: c1,c2,c3,e1,e2,e3,t1,t2,t3,s1,s2,s3
      real(8) :: uc,vc,wc,uv2,sndp2
      real(8) :: up,d,eigenx,eigeny,eigenz,eigenvis
      real(8) :: dtmin,mpi_dtmin
      real(8), save :: time = 0.d0
      integer :: ierr
      
      dtmin = 1.d10
      
      do k = 2,timestep%kmax
        do j = 2,timestep%jmax 
          do i = 2,timestep%imax            
            pv = variable%getpv(i,j,k)
            tv = variable%gettv(i,j,k)
            dv = variable%getdv(i,j,k)
            grd = grid%getgrd(i,j,k)
            cx1 = grid%getcx(i-1,j,k)
            cx2 = grid%getcx(i,j,k)
            ex1 = grid%getex(i,j-1,k)
            ex2 = grid%getex(i,j,k)
            tx1 = grid%gettx(i,j,k-1)
            tx2 = grid%gettx(i,j,k)
        
            c1 = 0.5d0*(cx1(1)+cx2(1))
            c2 = 0.5d0*(cx1(2)+cx2(2))
            c3 = 0.5d0*(cx1(3)+cx2(3))
            e1 = 0.5d0*(ex1(1)+ex2(1))
            e2 = 0.5d0*(ex1(2)+ex2(2))
            e3 = 0.5d0*(ex1(3)+ex2(3))
            t1 = 0.5d0*(tx1(1)+tx2(1))
            t2 = 0.5d0*(tx1(2)+tx2(2))
            t3 = 0.5d0*(tx1(3)+tx2(3))
            
            uc = dabs(c1*pv(2)+c2*pv(3)+c3*pv(4))
            vc = dabs(e1*pv(2)+e2*pv(3)+e3*pv(4))
            wc = dabs(t1*pv(2)+t2*pv(3)+t3*pv(4))
            
            s1 = c1**2+c2**2+c3**2
            s2 = e1**2+e2**2+e3**2
            s3 = t1**2+t2**2+t3**2
            
            uv2 = pv(2)**2+pv(3)**2+pv(4)**2
            sndp2 = timestep%getsndp2(dv(6),uv2)
            
            up = uc*(1.d0+sndp2/dv(6))
            d = dsqrt(uc**2*(1.d0-sndp2/dv(6))**2+4.d0*sndp2*s1)
            eigenx = 0.5d0*(up+d)

            up = vc*(1.d0+sndp2/dv(6))
            d = dsqrt(vc**2*(1.d0-sndp2/dv(6))**2+4.d0*sndp2*s2)
            eigeny = 0.5d0*(up+d)

            up = wc*(1.d0+sndp2/dv(6))
            d = dsqrt(wc**2*(1.d0-sndp2/dv(6))**2+4.d0*sndp2*s3)
            eigenz = 0.5d0*(up+d)
            
            eigenvis = timestep%geteigenvis(dv,tv)*(s1+s2+s3)/grd(1)
            
            timestep%dt(i,j,k) = timestep%cfl*grd(1)/(eigenx + eigeny + eigenz + 4.d0*eigenvis)

            if(dtmin.ge.timestep%dt(i,j,k)) then
              dtmin = timestep%dt(i,j,k)
            end if
           end do
        end do
      end do

      call mpi_reduce(dtmin,mpi_dtmin,1,mpi_real8,mpi_min,0,mpi_comm_world,ierr)
      call mpi_bcast(mpi_dtmin,1,mpi_real8,0,mpi_comm_world,ierr)
      timestep%dt = mpi_dtmin
      time = time + mpi_dtmin
      if((timestep%rank.eq.0).and.(mod(nt,timestep%nprint).eq.0)) write(*,*) 'Solution time=',time, mpi_dtmin
    end subroutine mintime
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine fixedtime(timestep,grid,variable,nt_phy,nt)
      implicit none
      class(t_fixedtime), intent(inout) :: timestep
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      integer, intent(in) :: nt_phy,nt

      timestep%dt = timestep%cfl

    end subroutine fixedtime
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function no_prec(timestep,snd2,uuu2) result(sndp2)
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
      
      sndp2 = snd2
      
    end function no_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function steady_prec(timestep,snd2,uuu2) result(sndp2)
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
      
      sndp2 = dmin1(snd2,dmax1(uuu2,timestep%uref**2))
      
    end function steady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteady_prec(timestep,snd2,uuu2) result(sndp2)
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
      
      sndp2 = dmin1(snd2,dmax1(uuu2,timestep%uref**2,timestep%str**2))
      
    end function unsteady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function euler(timestep,dv,tv) result(eigenvis)
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: dv(timestep%ndv),tv(timestep%ntv)
      real(8) :: eigenvis
      
      eigenvis = 0.d0
      
    end function euler
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function laminar(timestep,dv,tv) result(eigenvis)
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: dv(timestep%ndv),tv(timestep%ntv)
      real(8) :: eigenvis

      eigenvis = dmax1(dv(7)*tv(2)/(dv(8)+dv(12)*dv(7)*dv(1)-dv(11)*dv(8)*dv(1)),  &
                       4.d0/3.d0*tv(1)/dv(1)) 
    end function laminar
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function turbulent(timestep,dv,tv) result(eigenvis)
      implicit none
      class(t_timestep), intent(in) :: timestep
      real(8), intent(in) :: dv(timestep%ndv),tv(timestep%ntv)
      real(8) :: eigenvis
      real(8), parameter :: pr_t = 0.9d0
   
      eigenvis = dmax1(dv(7)*(tv(2)+dv(12)*tv(3)/pr_t)/(dv(8)+dv(12)*dv(7)*dv(1)    &
                       -dv(11)*dv(8)*dv(1)),4.d0/3.d0*(tv(1)+tv(3))/dv(1))
    end function turbulent
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getdt(timestep,i,j,k)
      implicit none
      class(t_timestep), intent(in) :: timestep
      integer, intent(in) :: i,j,k
      real(8) :: getdt
      
      getdt = timestep%dt(i,j,k)
      
    end function getdt
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module timestep_module
