module rhs_module
  use config_module
  use variable_module
  use grid_module
  use eos_module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! subordinate to rhs module  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  use flux_module
  use vsflux_module
  use muscl_module
  use cav_module
  use turbsource_module
  use unsteady_module
  use gravity_module
  use rotation_module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
  implicit none
  private
  public :: t_rhs
  
  type t_rhs
    private
    integer :: stencil,npv,ndv,ntv,ngrd,imax,jmax,kmax
    logical :: l_flux,l_muscl,l_cav,l_turbsource,l_vsflux,l_unsteady,l_gravity,l_rotation
    real(8) :: pref
    real(8) :: omega(3)
    real(8), dimension(:,:,:,:), allocatable :: res,icav,itt
    real(8), dimension(:,:,:), allocatable ::omega_cut
    real(8), dimension(:,:), allocatable :: ea,fa,ga,eva,fva,gva
    class(t_flux),       allocatable :: flux
    class(t_muscl),      allocatable :: muscl
    class(t_cav),        allocatable :: cav
    class(t_turbsource), allocatable :: turbsource
    class(t_vsflux),     allocatable :: vsflux
    class(t_unsteady),   allocatable :: unsteady
    class(t_gravity),    allocatable :: gravity
    class(t_rotation),   allocatable :: rotation
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: calrhs
      procedure :: getres
      procedure :: geticav
      procedure :: getitt
      procedure :: getomega_cut
  end type t_rhs
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(rhs,config,grid,variable)
      implicit none
      class(t_rhs), intent(out) :: rhs
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      
      rhs%stencil = config%getstencil()
      rhs%pref = config%getpref()
      rhs%omega = config%getomega()

      select case(config%getiturb())
      case(-3)

      case(-2)
        allocate(t_vsflux_laminar::rhs%vsflux)

      case(-1)
        allocate(t_vsflux_turbulent::rhs%vsflux)
        allocate(t_kepsilon::rhs%turbsource)

      case(0)
        allocate(t_vsflux_turbulent::rhs%vsflux)
        allocate(t_kwsst::rhs%turbsource)

      end select
      
      selectcase(config%getnscheme())
      case(1)
        allocate(t_roe::rhs%flux)
      case(2)
        allocate(t_roem::rhs%flux)
      case(3)
        allocate(t_ausmpwp::rhs%flux)
      case(4)
        allocate(t_ausmpup::rhs%flux)
      end select
            
      select case(config%getnmuscl())
      case(0)
      case(1)
        allocate(t_tvd::rhs%muscl)
      case(2)
        allocate(t_mlp::rhs%muscl)
      end select
                
      select case(config%getncav())
      case(0)
      case(1)
        allocate(t_merkle::rhs%cav)
      case(2)
        allocate(t_kunz::rhs%cav)
      case(3)
        allocate(t_singhal::rhs%cav)
      end select

      select case(config%getnsteady())
      case(0)
      case(1)
        allocate(t_unsteady::rhs%unsteady)
      end select

      select case(config%getgravity())
      case(0)
      case(-1,1,-2,2,-3,3)
        allocate(t_gravity::rhs%gravity)
      case default
      end select

      select case(config%getrotation())
      case(0)
      case(-1,1,-2,2,-3,3)
        allocate(t_rotation::rhs%rotation)
      case default
      end select

      rhs%l_flux       = .false.
      rhs%l_muscl      = .false.
      rhs%l_cav        = .false.
      rhs%l_turbsource = .false.
      rhs%l_vsflux     = .false.
      rhs%l_unsteady   = .false.
      rhs%l_gravity    = .false.
      rhs%l_rotation   = .false.

      rhs%ngrd = grid%getngrd()
      rhs%imax = grid%getimax()
      rhs%jmax = grid%getjmax()
      rhs%kmax = grid%getkmax()
      rhs%npv = variable%getnpv()
      rhs%ndv = variable%getndv()
      rhs%ntv = variable%getntv()
      
      
      if(allocated(rhs%flux)) then
        call rhs%flux%construct(config,grid,variable)
        rhs%l_flux = .true.
        allocate(rhs%ea(rhs%npv,rhs%imax),rhs%fa(rhs%npv,rhs%jmax),rhs%ga(rhs%npv,rhs%kmax))
      end if
      
      if(allocated(rhs%vsflux)) then
        call rhs%vsflux%construct(config,grid,variable)
        rhs%l_vsflux = .true.
        allocate(rhs%eva(rhs%npv,rhs%imax),rhs%fva(rhs%npv,rhs%jmax),rhs%gva(rhs%npv,rhs%kmax))
      end if
      
      if(allocated(rhs%muscl)) then
        call rhs%muscl%construct(config,variable)
        rhs%l_muscl = .true.
      end if
   
      if(allocated(rhs%cav)) then
        call rhs%cav%construct(config,grid,variable)
        if(config%gettimemethod().eq.3) then
          allocate(rhs%icav(4,2:rhs%imax,2:rhs%jmax,2:rhs%kmax))
        end if
        rhs%l_cav = .true.
      end if
      
      if(allocated(rhs%turbsource)) then
        call rhs%turbsource%construct(config,grid,variable)
        if(config%gettimemethod().eq.3) then
          allocate(rhs%itt(4,2:rhs%imax,2:rhs%jmax,2:rhs%kmax))
        end if
        allocate(rhs%omega_cut(2:rhs%imax,2:rhs%jmax,2:rhs%kmax))
        rhs%l_turbsource = .true.
      end if
      
      if(allocated(rhs%unsteady)) then
        call rhs%unsteady%construct(config,grid,variable)
        rhs%l_unsteady   = .true.
      end if

      if(allocated(rhs%gravity)) then
        call rhs%gravity%construct(config,grid,variable)
        rhs%l_gravity = .true.
      end if

      if(allocated(rhs%rotation)) then
        call rhs%rotation%construct(config,grid,variable)
        rhs%l_rotation = .true.
      end if

      allocate(rhs%res(rhs%npv,2:rhs%imax,2:rhs%jmax,2:rhs%kmax))

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(rhs)
      implicit none
      class(t_rhs), intent(inout) :: rhs
      
      if(rhs%l_flux) then
        call rhs%flux%destruct()  
        deallocate(rhs%flux)
        deallocate(rhs%ea,rhs%fa,rhs%ga)
      end if
      if(rhs%l_vsflux) then
        call rhs%vsflux%destruct()
        deallocate(rhs%vsflux)
        deallocate(rhs%eva,rhs%fva,rhs%gva)
      end if
      if(rhs%l_muscl) then
        call rhs%muscl%destruct()
        deallocate(rhs%muscl)
      end if
      if(rhs%l_cav) then
        call rhs%cav%destruct()
        deallocate(rhs%cav)
        if(allocated(rhs%icav)) deallocate(rhs%icav)
      end if
      if(rhs%l_turbsource) then
        call rhs%turbsource%destruct()
        deallocate(rhs%turbsource,rhs%omega_cut)
        if(allocated(rhs%itt)) deallocate(rhs%itt)
      end if
      if(rhs%l_unsteady) then
        call rhs%unsteady%destruct()
        deallocate(rhs%unsteady)
      end if
      if(rhs%l_gravity) then
        call rhs%gravity%destruct()
        deallocate(rhs%gravity)
      end if
      if(rhs%l_rotation) then
        call rhs%rotation%destruct()
        deallocate(rhs%rotation)
      end if
      deallocate(rhs%res)
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calrhs(rhs,grid,variable,eos)
      implicit none
      class(t_rhs), intent(inout) :: rhs
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      type(t_eos), intent(in) :: eos
      integer :: i,j,k
      integer :: ii,jj,kk,ll
      real(8) :: nx(3)
      real(8) :: ex1(3),ex2(3),ex3(3),ex4(3)
      real(8) :: tx1(3),tx2(3),tx3(3),tx4(3)
      real(8) :: grdl(rhs%ngrd),grdr(rhs%ngrd)
      real(8) :: x(rhs%stencil,rhs%npv)
      real(8) :: pvl(rhs%npv),pvr(rhs%npv)
      real(8) :: dvl(rhs%ndv),dvr(rhs%ndv)
      real(8) :: tvl(rhs%ntv),tvr(rhs%ntv)
      type(t_cav_result) :: cav_result
      type(t_turb_result) :: turb_result
      
      do k=2,rhs%kmax
        do j=2,rhs%jmax
          do i=1,rhs%imax
          
            nx = grid%getcx(i,j,k)
            
            ll = 0
            do kk = -1,1
              do jj = -1,1
                do ii = 0,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+ii,j+jj,k+kk)
                end do
              end do
            end do
            
            if(rhs%l_muscl) then
              do kk=-1,1
                do jj= -1,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i-1,j+jj,k+kk)
                end do
              end do
              do kk=-1,1
                do jj= -1,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+2,j+jj,k+kk)
                end do
              end do
              x(ll+1,:) = variable%getpv(i-2,j,k)
              x(ll+2,:) = variable%getpv(i+3,j,k)
              
              call rhs%muscl%setpv(x)
              call rhs%muscl%interpolation(pvl,pvr)

              call eos%deteos(pvl(1)+rhs%pref,pvl(5),pvl(6),pvl(7),dvl)
              call eos%deteos(pvr(1)+rhs%pref,pvr(5),pvr(6),pvr(7),dvr)
            else
              pvl = variable%getpv(i,j,k) 
              pvr = variable%getpv(i+1,j,k)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i+1,j,k)
            end if
            
            grdl = grid%getgrd(i,j,k)
            grdr = grid%getgrd(i+1,j,k)

            call rhs%flux%setnorm(nx)
            call rhs%flux%setgrd(grdl,grdr)
            call rhs%flux%setpv(pvl,pvr)
            call rhs%flux%setdv(dvl,dvr)
            call rhs%flux%setsdst(x(1:18,1))
            call rhs%flux%calflux(eos,rhs%ea(:,i))
            
            if(rhs%l_vsflux) then
              ex1 = grid%getex(i,j-1,k)
              ex2 = grid%getex(i+1,j-1,k)
              ex3 = grid%getex(i,j,k)
              ex4 = grid%getex(i+1,j,k)
              tx1 = grid%gettx(i,j,k-1)
              tx2 = grid%gettx(i+1,j,k-1)
              tx3 = grid%gettx(i,j,k)
              tx4 = grid%gettx(i+1,j,k)              
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i+1,j,k)
              tvl = variable%gettv(i,j,k)
              tvr = variable%gettv(i+1,j,k)
              
              call rhs%vsflux%setnorm(nx,ex1,ex2,ex3,ex4,tx1,tx2,tx3,tx4)
              call rhs%vsflux%setgrd(grdl,grdr)
              call rhs%vsflux%setpv(x)
              call rhs%vsflux%setdv(dvl,dvr)
              call rhs%vsflux%settv(tvl,tvr)
              call rhs%vsflux%calflux(rhs%eva(:,i))
            end if
          end do
          do i=2,rhs%imax
            rhs%res(:,i,j,k) = -( rhs%ea(:,i) - rhs%ea(:,i-1) )
            if(rhs%l_vsflux) rhs%res(:,i,j,k) = rhs%res(:,i,j,k) + (rhs%eva(:,i) - rhs%eva(:,i-1) )
          end do
        end do
      end do
      
      do i=2,rhs%imax
        do k=2,rhs%kmax
          do j=1,rhs%jmax
            
            nx = grid%getex(i,j,k)
            
            ll = 0
            do ii = -1,1
              do kk = -1,1
                do jj = 0,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+ii,j+jj,k+kk)
                end do
              end do
            end do
            
            if(rhs%l_muscl) then
              do ii = -1,1
                do kk = -1,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+ii,j-1,k+kk)
                end do
              end do
              do ii = -1,1
                do kk = -1,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+ii,j+2,k+kk)
                end do
              end do
              x(ll+1,:) = variable%getpv(i,j-2,k)
              x(ll+2,:) = variable%getpv(i,j+3,k)
              
              call rhs%muscl%setpv(x)
              call rhs%muscl%interpolation(pvl,pvr)

              call eos%deteos(pvl(1)+rhs%pref,pvl(5),pvl(6),pvl(7),dvl)
              call eos%deteos(pvr(1)+rhs%pref,pvr(5),pvr(6),pvr(7),dvr)
            else
              pvl = variable%getpv(i,j,k) 
              pvr = variable%getpv(i,j+1,k)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j+1,k)
            end if
            
            grdl = grid%getgrd(i,j,k)
            grdr = grid%getgrd(i,j+1,k)

            call rhs%flux%setnorm(nx)
            call rhs%flux%setgrd(grdl,grdr)
            call rhs%flux%setpv(pvl,pvr)
            call rhs%flux%setdv(dvl,dvr)
            call rhs%flux%setsdst(x(1:18,1))
            call rhs%flux%calflux(eos,rhs%fa(:,j))
            
            if(rhs%l_vsflux) then
              ex1 = grid%gettx(i,j,k-1)
              ex2 = grid%gettx(i,j+1,k-1)
              ex3 = grid%gettx(i,j,k)
              ex4 = grid%gettx(i,j+1,k)
              tx1 = grid%getcx(i-1,j,k)
              tx2 = grid%getcx(i,j,k)
              tx3 = grid%getcx(i-1,j+1,k)
              tx4 = grid%getcx(i,j+1,k)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j+1,k)
              tvl = variable%gettv(i,j,k)
              tvr = variable%gettv(i,j+1,k)

              call rhs%vsflux%setnorm(nx,ex1,ex2,ex3,ex4,tx1,tx2,tx3,tx4)
              call rhs%vsflux%setgrd(grdl,grdr)
              call rhs%vsflux%setpv(x)
              call rhs%vsflux%setdv(dvl,dvr)
              call rhs%vsflux%settv(tvl,tvr)
              call rhs%vsflux%calflux(rhs%fva(:,j))
            end if
          end do
          do j=2,rhs%jmax
            rhs%res(:,i,j,k) = rhs%res(:,i,j,k) -( rhs%fa(:,j) - rhs%fa(:,j-1) )
            if(rhs%l_vsflux) rhs%res(:,i,j,k) = rhs%res(:,i,j,k) + (rhs%fva(:,j) - rhs%fva(:,j-1) )
          end do
        end do
      end do

      do j=2,rhs%jmax
        do i=2,rhs%imax
          do k=1,rhs%kmax
          
            nx = grid%gettx(i,j,k)
            
            ll = 0
            do jj = -1,1
              do ii = -1,1
                do kk = 0,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+ii,j+jj,k+kk)
                end do
              end do
            end do
            
            if(rhs%l_muscl) then
              do jj = -1,1
                do ii = -1,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+ii,j+jj,k-1)
                end do
              end do
              do jj = -1,1
                do ii = -1,1
                  ll = ll+1
                  x(ll,:) = variable%getpv(i+ii,j+jj,k+2)
                end do
              end do
              x(ll+1,:) = variable%getpv(i,j,k-2)
              x(ll+2,:) = variable%getpv(i,j,k+3)
              
              call rhs%muscl%setpv(x)
              call rhs%muscl%interpolation(pvl,pvr)

              call eos%deteos(pvl(1)+rhs%pref,pvl(5),pvl(6),pvl(7),dvl)
              call eos%deteos(pvr(1)+rhs%pref,pvr(5),pvr(6),pvr(7),dvr)
            else
              pvl = variable%getpv(i,j,k) 
              pvr = variable%getpv(i,j,k+1)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j,k+1)
            end if
            
            grdl = grid%getgrd(i,j,k)
            grdr = grid%getgrd(i,j,k+1)

            call rhs%flux%setnorm(nx)
            call rhs%flux%setgrd(grdl,grdr)
            call rhs%flux%setpv(pvl,pvr)
            call rhs%flux%setdv(dvl,dvr)
            call rhs%flux%setsdst(x(1:18,1))
            call rhs%flux%calflux(eos,rhs%ga(:,k))
            
            if(rhs%l_vsflux) then
              ex1 = grid%getcx(i-1,j,k)
              ex2 = grid%getcx(i-1,j,k+1)
              ex3 = grid%getcx(i,j,k)
              ex4 = grid%getcx(i,j,k+1)
              tx1 = grid%getex(i,j-1,k)
              tx2 = grid%getex(i,j,k)
              tx3 = grid%getex(i,j-1,k+1)
              tx4 = grid%getex(i,j,k+1)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j,k+1)
              tvl = variable%gettv(i,j,k)
              tvr = variable%gettv(i,j,k+1)

              call rhs%vsflux%setnorm(nx,ex1,ex2,ex3,ex4,tx1,tx2,tx3,tx4)
              call rhs%vsflux%setgrd(grdl,grdr)
              call rhs%vsflux%setpv(x)
              call rhs%vsflux%setdv(dvl,dvr)
              call rhs%vsflux%settv(tvl,tvr)
              call rhs%vsflux%calflux(rhs%gva(:,k))
            end if
          end do
          do k=2,rhs%kmax
            rhs%res(:,i,j,k) = rhs%res(:,i,j,k) - ( rhs%ga(:,k) - rhs%ga(:,k-1) ) 
            if(rhs%l_vsflux) rhs%res(:,i,j,k) = rhs%res(:,i,j,k) + (rhs%gva(:,k) - rhs%gva(:,k-1) )
          end do
        end do
      end do
      
      do k=2,rhs%kmax
        do j=2,rhs%jmax
          do i=2,rhs%imax
            ex1 = grid%getcx(i-1,j,k)
            ex2 = grid%getcx(i,j,k)
            ex3 = grid%getex(i,j-1,k)
            ex4 = grid%getex(i,j,k)
            tx1 = grid%gettx(i,j,k-1)
            tx2 = grid%gettx(i,j,k)
            grdl = grid%getgrd(i,j,k)
            x(1,:) = variable%getpv(i-1,j,k)
            x(2,:) = variable%getpv(i,j,k)
            x(3,:) = variable%getpv(i+1,j,k)
            x(4,:) = variable%getpv(i,j-1,k)
            x(5,:) = variable%getpv(i,j+1,k)
            x(6,:) = variable%getpv(i,j,k-1)
            x(7,:) = variable%getpv(i,j,k+1)
            dvl = variable%getdv(i,j,k)
            tvl = variable%gettv(i,j,k)
        
            if(rhs%l_cav) then
              call rhs%cav%setgrd(grdl)
              call rhs%cav%setpv(x)
              call rhs%cav%setdv(dvl)
              cav_result = rhs%cav%cavsource(eos)
              rhs%res(6,i,j,k) = rhs%res(6,i,j,k) + cav_result%cavsource
              if(allocated(rhs%icav)) rhs%icav(:,i,j,k) = cav_result%icav(:)
            end if
            
            if(rhs%l_turbsource) then
              call rhs%turbsource%setnorm(ex1,ex2,ex3,ex4,tx1,tx2)
              call rhs%turbsource%setgrd(grdl)
              call rhs%turbsource%setpv(x)
              call rhs%turbsource%setdv(dvl)
              call rhs%turbsource%settv(tvl)            
              turb_result = rhs%turbsource%calturbsource()
              rhs%omega_cut(i,j,k) = turb_result%omega_cut
              rhs%res(8:9,i,j,k) = rhs%res(8:9,i,j,k) + turb_result%source(:)
              if(allocated(rhs%itt)) rhs%itt(:,i,j,k) = turb_result%itt(:)
            end if
            
            if(rhs%l_unsteady) then
              pvl = variable%getqq(1,i,j,k)
              pvr = variable%getqq(2,i,j,k)
              call rhs%unsteady%setgrd(grdl)
              call rhs%unsteady%setpv(x)
              call rhs%unsteady%setdv(dvl)
              call rhs%unsteady%setqq(pvl,pvr)
              rhs%res(:,i,j,k) = rhs%res(:,i,j,k) - rhs%unsteady%unsteadysource()
            end if

            if(rhs%l_gravity) then
              call rhs%gravity%setdv(dvl)
              call rhs%gravity%setgrd(grdl)
              rhs%res(:,i,j,k) = rhs%res(:,i,j,k) + rhs%gravity%getgravitysource()
            end if

            if(rhs%l_rotation) then
              call rhs%rotation%setpv(x)
              call rhs%rotation%setdv(dvl)
              call rhs%rotation%setgrd(grdl)
              rhs%res(:,i,j,k) = rhs%res(:,i,j,k) - rhs%rotation%getrotationsource()
            end if
          end do
        end do
      end do
    end subroutine calrhs
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getres(rhs,n,i,j,k)
      implicit none
      class(t_rhs), intent(in) :: rhs
      integer, intent(in) :: n,i,j,k
      real(8) :: getres
      
      getres = rhs%res(n,i,j,k)
      
    end function getres
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function geticav(rhs,i,j,k)
      implicit none
      class(t_rhs), intent(in) :: rhs
      integer, intent(in) :: i,j,k
      real(8) :: geticav(4)
      
      geticav = rhs%icav(:,i,j,k)
      
    end function geticav
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getitt(rhs,i,j,k)
      implicit none
      class(t_rhs), intent(in) :: rhs
      integer, intent(in) :: i,j,k
      real(8) :: getitt(4)
      
      getitt= rhs%itt(:,i,j,k)
      
    end function getitt
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getomega_cut(rhs,i,j,k)
      implicit none
      class(t_rhs), intent(in) :: rhs
      integer, intent(in) :: i,j,k
      real(8) :: getomega_cut
      
      getomega_cut= rhs%omega_cut(i,j,k)
      
    end function getomega_cut
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module rhs_module
