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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
  implicit none
  private
  public :: t_rhs
  
  type t_rhs
    private
    integer :: npv,ndv,ntv,ngrd,imax,jmax,kmax
    logical :: l_flux,l_muscl,l_cav,l_turbsource,l_vsflux,l_unsteady
    real(8) :: pref
    real(8), dimension(:,:,:,:), allocatable :: res,icav,itt
    real(8), dimension(:,:,:), allocatable ::omega_cut
    real(8), dimension(:,:), allocatable :: ea,fa,ga,eva,fva,gva
    class(t_flux),       allocatable :: flux
    class(t_muscl),      allocatable :: muscl
    class(t_cav),        allocatable :: cav
    class(t_turbsource), allocatable :: turbsource
    class(t_vsflux),     allocatable :: vsflux
    class(t_unsteady),   allocatable :: unsteady
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
      
      rhs%pref = config%getpref()
      
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

      rhs%l_flux       = .false.
      rhs%l_muscl      = .false.
      rhs%l_cav        = .false.
      rhs%l_turbsource = .false.
      rhs%l_vsflux     = .false.
      rhs%l_unsteady   = .false.

      rhs%ngrd = grid%getngrd()
      rhs%imax = grid%getimax()
      rhs%jmax = grid%getjmax()
      rhs%kmax = grid%getkmax()
      rhs%npv = variable%getnpv()
      rhs%ndv = variable%getndv()
      rhs%ntv = variable%getntv()
      
      
      if(allocated(rhs%flux)) then
        call rhs%flux%construct(config,variable)
        rhs%l_flux = .true.
        allocate(rhs%ea(rhs%npv,rhs%imax),rhs%fa(rhs%npv,rhs%jmax),rhs%ga(rhs%npv,rhs%kmax))
      end if
      
      if(allocated(rhs%vsflux)) then
        call rhs%vsflux%construct(config,variable)
        rhs%l_vsflux = .true.
        allocate(rhs%eva(rhs%npv,rhs%imax),rhs%fva(rhs%npv,rhs%jmax),rhs%gva(rhs%npv,rhs%kmax))
      end if
      
      if(allocated(rhs%muscl)) then
        call rhs%muscl%construct(config,variable)
        rhs%l_muscl = .true.
      end if
   
      if(allocated(rhs%cav)) then
        call rhs%cav%construct(config)
        if(config%gettimemethod().eq.3) then
          allocate(rhs%icav(4,2:rhs%imax,2:rhs%jmax,2:rhs%kmax))
        end if
        rhs%l_cav = .true.
      end if
      
      if(allocated(rhs%turbsource)) then
        call rhs%turbsource%construct(config)
        if(config%gettimemethod().eq.3) then
          allocate(rhs%itt(4,2:rhs%imax,2:rhs%jmax,2:rhs%kmax))
        end if
        allocate(rhs%omega_cut(2:rhs%imax,2:rhs%jmax,2:rhs%kmax))
        rhs%l_turbsource = .true.
      end if
      
      if(allocated(rhs%unsteady)) then
        call rhs%unsteady%construct(config,variable)
        rhs%l_unsteady   = .true.
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
      real(8), target :: nx(3)
      real(8), target :: ex1(3),ex2(3),ex3(3),ex4(3)
      real(8), target :: tx1(3),tx2(3),tx3(3),tx4(3)
      real(8), target :: grdl(rhs%ngrd),grdr(rhs%ngrd)
      real(8), target :: x(38,rhs%npv)
      real(8), target :: pvl(rhs%npv),pvr(rhs%npv)
      real(8), target :: dvl(rhs%ndv),dvr(rhs%ndv)
      real(8), target :: tvl(rhs%ntv),tvr(rhs%ntv)
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
              
              rhs%muscl%x => x
              call rhs%muscl%interpolation(pvl,pvr)
              if((pvl(6).lt.0.d0).or.(pvl(6).gt.1.d0)) pvl(6) = x(9,6)
              if((pvr(6).lt.0.d0).or.(pvr(6).gt.1.d0)) pvr(6) = x(10,6)
              if((pvl(7).lt.0.d0).or.(pvl(7).gt.1.d0)) pvl(7) = x(9,7)
              if((pvr(7).lt.0.d0).or.(pvr(7).gt.1.d0)) pvr(7) = x(10,7)
              
              if(rhs%l_turbsource) then
                if((pvl(8).lt.0.d0).and.(x(9,8).gt.0.d0))  pvl(8) = x(9,8)
                if((pvr(8).lt.0.d0).and.(x(10,8).gt.0.d0)) pvr(8) = x(10,8)
                if((pvl(9).lt.0.d0).and.(x(9,9).gt.0.d0))  pvl(9) = x(9,9)
                if((pvr(9).lt.0.d0).and.(x(10,9).gt.0.d0)) pvr(9) = x(10,9)
              end if
              
              call eos%deteos(pvl(1)+rhs%pref,pvl(5),pvl(6),pvl(7),dvl)
              call eos%deteos(pvr(1)+rhs%pref,pvr(5),pvr(6),pvr(7),dvr)
            else
              pvl = variable%getpv(i,j,k) 
              pvr = variable%getpv(i+1,j,k)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i+1,j,k)
            end if
            
            rhs%flux%nx => nx  
            rhs%flux%pvl => pvl 
            rhs%flux%pvr => pvr
            rhs%flux%dvl => dvl 
            rhs%flux%dvr => dvr
            rhs%flux%sdst => x(1:18,1)
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
              grdl = grid%getgrd(i,j,k)
              grdr = grid%getgrd(i+1,j,k)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i+1,j,k)
              tvl = variable%gettv(i,j,k)
              tvr = variable%gettv(i+1,j,k)
              
              rhs%vsflux%nx => nx 
              rhs%vsflux%ex1 => ex1 
              rhs%vsflux%ex2 => ex2 
              rhs%vsflux%ex3 => ex3 
              rhs%vsflux%ex4 => ex4 
              rhs%vsflux%tx1 => tx1 
              rhs%vsflux%tx2 => tx2 
              rhs%vsflux%tx3 => tx3 
              rhs%vsflux%tx4 => tx4 
              rhs%vsflux%grdl => grdl
              rhs%vsflux%grdr => grdr
              rhs%vsflux%dvl => dvl
              rhs%vsflux%dvr => dvr
              rhs%vsflux%tvl => tvl
              rhs%vsflux%tvr => tvr
              rhs%vsflux%pv => x(1:18,:)
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
              
              rhs%muscl%x => x
              call rhs%muscl%interpolation(pvl,pvr)
              if((pvl(6).lt.0.d0).or.(pvl(6).gt.1.d0)) pvl(6) = x(9,6)
              if((pvr(6).lt.0.d0).or.(pvr(6).gt.1.d0)) pvr(6) = x(10,6)
              if((pvl(7).lt.0.d0).or.(pvl(7).gt.1.d0)) pvl(7) = x(9,7)
              if((pvr(7).lt.0.d0).or.(pvr(7).gt.1.d0)) pvr(7) = x(10,7)
              
              if(rhs%l_turbsource) then
                if((pvl(8).lt.0.d0).and.(x(9,8).gt.0.d0))  pvl(8) = x(9,8)
                if((pvr(8).lt.0.d0).and.(x(10,8).gt.0.d0)) pvr(8) = x(10,8)
                if((pvl(9).lt.0.d0).and.(x(9,9).gt.0.d0))  pvl(9) = x(9,9)
                if((pvr(9).lt.0.d0).and.(x(10,9).gt.0.d0)) pvr(9) = x(10,9)
              end if
              
              call eos%deteos(pvl(1)+rhs%pref,pvl(5),pvl(6),pvl(7),dvl)
              call eos%deteos(pvr(1)+rhs%pref,pvr(5),pvr(6),pvr(7),dvr)
            else
              pvl = variable%getpv(i,j,k) 
              pvr = variable%getpv(i,j+1,k)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j+1,k)
            end if


            rhs%flux%nx => nx  
            rhs%flux%pvl => pvl 
            rhs%flux%pvr => pvr
            rhs%flux%dvl => dvl 
            rhs%flux%dvr => dvr
            rhs%flux%sdst => x(1:18,1)
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
              grdl = grid%getgrd(i,j,k)
              grdr = grid%getgrd(i,j+1,k)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j+1,k)
              tvl = variable%gettv(i,j,k)
              tvr = variable%gettv(i,j+1,k)

              
              
              rhs%vsflux%nx => nx 
              rhs%vsflux%ex1 => ex1 
              rhs%vsflux%ex2 => ex2 
              rhs%vsflux%ex3 => ex3 
              rhs%vsflux%ex4 => ex4 
              rhs%vsflux%tx1 => tx1 
              rhs%vsflux%tx2 => tx2 
              rhs%vsflux%tx3 => tx3 
              rhs%vsflux%tx4 => tx4 
              rhs%vsflux%grdl => grdl
              rhs%vsflux%grdr => grdr
              rhs%vsflux%dvl => dvl
              rhs%vsflux%dvr => dvr
              rhs%vsflux%tvl => tvl
              rhs%vsflux%tvr => tvr
              rhs%vsflux%pv => x(1:18,:)
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
              
              rhs%muscl%x => x
              call rhs%muscl%interpolation(pvl,pvr)
              if((pvl(6).lt.0.d0).or.(pvl(6).gt.1.d0)) pvl(6) = x(9,6)
              if((pvr(6).lt.0.d0).or.(pvr(6).gt.1.d0)) pvr(6) = x(10,6)
              if((pvl(7).lt.0.d0).or.(pvl(7).gt.1.d0)) pvl(7) = x(9,7)
              if((pvr(7).lt.0.d0).or.(pvr(7).gt.1.d0)) pvr(7) = x(10,7)
              
              if(rhs%l_turbsource) then
                if((pvl(8).lt.0.d0).and.(x(9,8).gt.0.d0))  pvl(8) = x(9,8)
                if((pvr(8).lt.0.d0).and.(x(10,8).gt.0.d0)) pvr(8) = x(10,8)
                if((pvl(9).lt.0.d0).and.(x(9,9).gt.0.d0))  pvl(9) = x(9,9)
                if((pvr(9).lt.0.d0).and.(x(10,9).gt.0.d0)) pvr(9) = x(10,9)
              end if
              
              call eos%deteos(pvl(1)+rhs%pref,pvl(5),pvl(6),pvl(7),dvl)
              call eos%deteos(pvr(1)+rhs%pref,pvr(5),pvr(6),pvr(7),dvr)
            else
              pvl = variable%getpv(i,j,k) 
              pvr = variable%getpv(i,j,k+1)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j,k+1)
            end if

            
            
            rhs%flux%nx => nx  
            rhs%flux%pvl => pvl 
            rhs%flux%pvr => pvr
            rhs%flux%dvl => dvl 
            rhs%flux%dvr => dvr
            rhs%flux%sdst => x(1:18,1)
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
              grdl = grid%getgrd(i,j,k)
              grdr = grid%getgrd(i,j,k+1)
              dvl = variable%getdv(i,j,k)
              dvr = variable%getdv(i,j,k+1)
              tvl = variable%gettv(i,j,k)
              tvr = variable%gettv(i,j,k+1)

              rhs%vsflux%nx => nx 
              rhs%vsflux%ex1 => ex1 
              rhs%vsflux%ex2 => ex2 
              rhs%vsflux%ex3 => ex3 
              rhs%vsflux%ex4 => ex4 
              rhs%vsflux%tx1 => tx1 
              rhs%vsflux%tx2 => tx2 
              rhs%vsflux%tx3 => tx3 
              rhs%vsflux%tx4 => tx4 
              rhs%vsflux%grdl => grdl
              rhs%vsflux%grdr => grdr
              rhs%vsflux%dvl => dvl
              rhs%vsflux%dvr => dvr
              rhs%vsflux%tvl => tvl
              rhs%vsflux%tvr => tvr
              rhs%vsflux%pv => x(1:18,:)
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
              rhs%cav%grd => grdl
              rhs%cav%pv => x(2,:)
              rhs%cav%dv => dvl
              cav_result = rhs%cav%cavsource(eos)
              rhs%res(6,i,j,k) = rhs%res(6,i,j,k) + cav_result%cavsource
              if(allocated(rhs%icav)) rhs%icav(:,i,j,k) = cav_result%icav(:)
            end if
            
            if(rhs%l_turbsource) then
              rhs%turbsource%cx1 => ex1
              rhs%turbsource%cx2 => ex2
              rhs%turbsource%ex1 => ex3
              rhs%turbsource%ex2 => ex4
              rhs%turbsource%tx1 => tx1
              rhs%turbsource%tx2 => tx2
              rhs%turbsource%grd => grdl
              rhs%turbsource%pv => x(1:7,:)
              rhs%turbsource%dv => dvl
              rhs%turbsource%tv => tvl
              turb_result = rhs%turbsource%calturbsource()
              rhs%omega_cut(i,j,k) = turb_result%omega_cut
              rhs%res(8:9,i,j,k) = rhs%res(8:9,i,j,k) + turb_result%source(:)
              if(allocated(rhs%itt)) rhs%itt(:,i,j,k) = turb_result%itt(:)
            end if
            
            if(rhs%l_unsteady) then
              pvl = variable%getqq(1,i,j,k)
              pvr = variable%getqq(2,i,j,k)
              rhs%unsteady%grd => grdl
              rhs%unsteady%pv => x(2,:)
              rhs%unsteady%dv => dvl
              rhs%unsteady%qq1 => pvl
              rhs%unsteady%qq2 => pvr
              rhs%res(:,i,j,k) = rhs%res(:,i,j,k) - rhs%unsteady%unsteadysource()
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
