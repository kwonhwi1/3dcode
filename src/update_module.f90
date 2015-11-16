module update_module
  use config_module 
  use grid_module
  use variable_module
  use eos_module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! subordinate to update module  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
  use rhs_module
  use timestep_module
  use lhs_module
  use jacobian_module
  use eddy_module
  use bc_module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
  implicit none
  private
  public :: t_update,t_eulerex,t_rk3rd,t_lusgs
  
  type, abstract :: t_update
    private
    integer :: npv,ndv,ntv,ngrd,imax,jmax,kmax
    logical :: l_timestep,l_lhs,l_eddy,l_jac,l_tv,l_turb,l_cav
    real(8) :: pref,kref,oref
    class(t_lhs), allocatable :: lhs
    class(t_timestep), allocatable :: timestep
    class(t_eddy), allocatable :: eddy
    class(t_jac), allocatable :: jac
    class(t_rhs), allocatable :: rhs
    class(t_bc), allocatable :: bc
    contains
      procedure :: construct => construct_eulerex
      procedure :: destruct  => destruct_eulerex
      procedure(p_timeinteg), deferred :: timeinteg
  end type t_update
  
  abstract interface
    subroutine p_timeinteg(update,grid,variable,eos,nt_phy,nt)
      import t_update
      import t_grid
      import t_variable
      import t_eos
      implicit none
      class(t_update), intent(inout) :: update
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(in) :: nt_phy,nt
    end subroutine p_timeinteg
  end interface
  
  type, extends(t_update) :: t_eulerex
    contains
      procedure :: timeinteg => eulerex
  end type t_eulerex
  
  type, extends(t_update) :: t_rk3rd
    private
    real(8) :: a1(3),a2(3),a3(3)
    real(8), dimension(:,:,:,:), allocatable :: rk
    contains
      procedure :: construct => construct_rk3rd
      procedure :: destruct  => destruct_rk3rd
      procedure :: timeinteg => rk3rd
  end type t_rk3rd
  
  type, extends(t_update) :: t_lusgs
    private
    real(8), dimension(:,:,:,:), allocatable  :: dqs,dcv
    contains
      procedure :: construct => construct_lusgs
      procedure :: destruct  => destruct_lusgs
      procedure :: timeinteg => lusgs
  end type t_lusgs
  

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct_eulerex(update,config,grid,eos)
      implicit none
      class(t_update), intent(out) :: update
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      class(t_eos), intent(in) :: eos

      update%npv = config%getnpv()
      update%ndv = config%getndv()
      update%ntv = config%getntv()
      update%pref = config%getpref()
      update%kref = config%getkref()*1.d-12
      update%oref = config%getoref()*1.d-12
      
      update%l_timestep = .false.
      update%l_lhs = .false.
      update%l_eddy = .false.
      update%l_jac = .false.
      
      select case(config%getlocal())
      case(-1)
        allocate(t_fixedtime::update%timestep)
      case(0)
        allocate(t_mintime::update%timestep)
      case(1)
        allocate(t_localtime::update%timestep)
      end select
            
      select case(config%getiturb())
      case(0)
        update%l_tv = .true.
        update%l_turb = .true.
        allocate(t_eddy_kwsst::update%eddy)
        allocate(t_lhs_flowturball_ex::update%lhs)
      case(-1)
        update%l_tv = .true.
        update%l_turb = .true.
        allocate(t_eddy_ke::update%eddy)
        allocate(t_lhs_flowturball_ex::update%lhs)
      case(-2)
        update%l_tv = .true.
        update%l_turb = .false.
        allocate(t_lhs_flowonly_ex::update%lhs)
      case(-3)
        update%l_tv = .false.
        update%l_turb = .false.
        allocate(t_lhs_flowonly_ex::update%lhs)
      end select
      
      select case(config%getncav())
      case(0)
        update%l_cav = .false.
      case(1,2,3)
        update%l_cav = .true.
      end select
      
      if(allocated(update%timestep)) then
        call update%timestep%construct(config,grid)
        update%l_timestep = .true.
      end if

      if(allocated(update%lhs)) then
        call update%lhs%construct(config,grid)
        update%l_lhs = .true.
      end if
      
      if(allocated(update%eddy)) then
        call update%eddy%construct(config,grid)
        update%l_eddy = .true.
      end if
      
      allocate(update%rhs)
      call update%rhs%construct(config,grid)

      allocate(update%bc)
      call update%bc%construct(config,grid,eos)

      update%ngrd = grid%getngrd()
      update%imax = grid%getimax()
      update%jmax = grid%getjmax()
      update%kmax = grid%getkmax()
      
    end subroutine construct_eulerex
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct_rk3rd(update,config,grid,eos)
      implicit none
      class(t_rk3rd), intent(out) :: update
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      class(t_eos), intent(in) :: eos

      update%npv = config%getnpv()
      update%ndv = config%getndv()
      update%ntv = config%getntv()
      update%pref = config%getpref()
      update%kref = config%getkref()*1.d-12
      update%oref = config%getoref()*1.d-12
      
      update%l_timestep = .false.
      update%l_lhs = .false.
      update%l_eddy = .false.
      update%l_jac = .false.
      
      select case(config%getlocal())
      case(-1)
        allocate(t_fixedtime::update%timestep)
      case(0)
        allocate(t_mintime::update%timestep)
      case(1)
        allocate(t_localtime::update%timestep)
      end select
            
      select case(config%getiturb())
      case(0)
        update%l_tv = .true.
        update%l_turb = .true.
        allocate(t_eddy_kwsst::update%eddy)
        allocate(t_lhs_flowturball_ex::update%lhs)
      case(-1)
        update%l_tv = .true.
        update%l_turb = .true.
        allocate(t_eddy_ke::update%eddy)
        allocate(t_lhs_flowturball_ex::update%lhs)
      case(-2)
        update%l_tv = .true.
        update%l_turb = .false.
        allocate(t_lhs_flowonly_ex::update%lhs)
      case(-3)
        update%l_tv = .false.
        update%l_turb = .false.
        allocate(t_lhs_flowonly_ex::update%lhs)
      end select
      
      select case(config%getncav())
      case(0)
        update%l_cav = .false.
      case(1,2,3)
        update%l_cav = .true.
      end select
      
      if(allocated(update%timestep)) then
        call update%timestep%construct(config,grid)
        update%l_timestep = .true.
      end if

      if(allocated(update%lhs)) then
        call update%lhs%construct(config,grid)
        update%l_lhs = .true.
      end if
      
      if(allocated(update%eddy)) then
        call update%eddy%construct(config,grid)
        update%l_eddy = .true.
      end if
      
      allocate(update%rhs)
      call update%rhs%construct(config,grid)

      allocate(update%bc)
      call update%bc%construct(config,grid,eos)

      update%ngrd = grid%getngrd()
      update%imax = grid%getimax()
      update%jmax = grid%getjmax()
      update%kmax = grid%getkmax()
      
      update%a1=(/0.d0,0.75d0,1.d0/3.d0/)
      update%a2=(/1.d0,0.25d0,2.d0/3.d0/)
      update%a3=(/1.d0,0.25d0,2.d0/3.d0/)
      
      allocate(update%rk(update%npv,update%imax,update%jmax,update%kmax))
    end subroutine construct_rk3rd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct_lusgs(update,config,grid,eos)
      implicit none
      class(t_lusgs), intent(out) :: update
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      class(t_eos), intent(in) :: eos

      update%npv = config%getnpv()
      update%ndv = config%getndv()
      update%ntv = config%getntv()
      update%pref = config%getpref()
      update%kref = config%getkref()*1.d-12
      update%oref = config%getoref()*1.d-12
      
      update%l_timestep = .false.
      update%l_lhs = .false.
      update%l_eddy = .false.
      update%l_jac = .false.
      
      select case(config%getlocal())
      case(-1)
        allocate(t_fixedtime::update%timestep)
      case(0)
        allocate(t_mintime::update%timestep)
      case(1)
        allocate(t_localtime::update%timestep)
      end select
            
      select case(config%getiturb())
      case(0)
        update%l_tv = .true.
        update%l_turb = .true.
        allocate(t_eddy_kwsst::update%eddy)
        allocate(t_lhs_flowturball::update%lhs) 
        allocate(t_jac_flowturball::update%jac)
      case(-1)
        update%l_tv = .true.
        update%l_turb = .true.
        allocate(t_eddy_ke::update%eddy)
        allocate(t_lhs_flowturball::update%lhs)
        allocate(t_jac_flowturball::update%jac)
      case(-2)
        update%l_tv = .true.
        update%l_turb = .false.
        allocate(t_lhs_flowonly::update%lhs)
        allocate(t_jac_flowonly::update%jac)
      case(-3)
        update%l_tv = .false.
        update%l_turb = .false.
        allocate(t_lhs_flowonly::update%lhs)
        allocate(t_jac_flowonly::update%jac)
      end select
      
      select case(config%getncav())
      case(0)
        update%l_cav = .false.
      case(1,2,3)
        update%l_cav = .true.
      end select
      
      if(allocated(update%timestep)) then
        call update%timestep%construct(config,grid)
        update%l_timestep = .true.
      end if

      if(allocated(update%lhs)) then
        call update%lhs%construct(config,grid)
        update%l_lhs = .true.
      end if
      
      if(allocated(update%eddy)) then
        call update%eddy%construct(config,grid)
        update%l_eddy = .true.
      end if
      
      if(allocated(update%jac)) then
        call update%jac%construct(config,grid)
        update%l_jac = .true.
      end if
      
      allocate(update%rhs)
      call update%rhs%construct(config,grid)

      allocate(update%bc)
      call update%bc%construct(config,grid,eos)

      update%ngrd = grid%getngrd()
      update%imax = grid%getimax()
      update%jmax = grid%getjmax()
      update%kmax = grid%getkmax()
      
      allocate(update%dqs(update%npv,update%imax+1,update%jmax+1,update%kmax+1))
      allocate(update%dcv(update%npv,update%imax+1,update%jmax+1,update%kmax+1))
      
    end subroutine construct_lusgs
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct_eulerex(update)
      implicit none
      class(t_update), intent(inout) :: update
          
      if(update%l_timestep) then
        call update%timestep%destruct()
        deallocate(update%timestep)
      end if
      
      if(update%l_lhs) then
        call update%lhs%destruct()
        deallocate(update%lhs)
      end if
      
      if(update%l_eddy) then
        call update%eddy%destruct()
        deallocate(update%eddy)
      end if

      call update%rhs%destruct()
      deallocate(update%rhs)

      call update%bc%destruct()
      deallocate(update%bc)
      
    end subroutine destruct_eulerex
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct_rk3rd(update)
      implicit none
      class(t_rk3rd), intent(inout) :: update
          
      if(update%l_timestep) then
        call update%timestep%destruct()
        deallocate(update%timestep)
      end if
      
      if(update%l_lhs) then
        call update%lhs%destruct()
        deallocate(update%lhs)
      end if
      
      if(update%l_eddy) then
        call update%eddy%destruct()
        deallocate(update%eddy)
      end if

      call update%rhs%destruct()
      deallocate(update%rhs)

      call update%bc%destruct()
      deallocate(update%bc)
      
      deallocate(update%rk)
    end subroutine destruct_rk3rd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct_lusgs(update)
      implicit none
      class(t_lusgs), intent(inout) :: update
          
      if(update%l_timestep) then
        call update%timestep%destruct()
        deallocate(update%timestep)
      end if
      
      if(update%l_lhs) then
        call update%lhs%destruct()
        deallocate(update%lhs)
      end if
      
      if(update%l_eddy) then
        call update%eddy%destruct()
        deallocate(update%eddy)
      end if

      if(update%l_jac) then
        call update%jac%destruct()
        deallocate(update%jac)
      end if
      
      call update%rhs%destruct()
      deallocate(update%rhs)

      call update%bc%destruct()
      deallocate(update%bc)

      deallocate(update%dqs,update%dcv)
      
    end subroutine destruct_lusgs
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine eulerex(update,grid,variable,eos,nt_phy,nt)
      implicit none
      class(t_eulerex), intent(inout) :: update
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(in) :: nt_phy,nt
      integer :: i,j,k,n,l
      real(8) :: nx1(3),nx2(3),nx3(3),nx4(3),nx5(3),nx6(3)
      real(8) :: pv(update%npv),dv(update%ndv),tv(update%ntv)
      real(8) :: grd(update%ngrd),dt
      real(8) :: x(7,update%npv)
      real(8) :: res(update%npv)

      
      call update%bc%setbc(grid,variable,eos)
      call update%timestep%caltimestep(grid,variable,nt_phy,nt)
      call update%rhs%calrhs(grid,variable,eos)
      
      do k=2,update%kmax
        do j=2,update%jmax
          do i=2,update%imax
          
            nx1 = grid%getcx(i-1,j,k)
            nx2 = grid%getcx(i,j,k)
            nx3 = grid%getex(i,j-1,k)
            nx4 = grid%getex(i,j,k)
            nx5 = grid%gettx(i,j,k-1)
            nx6 = grid%gettx(i,j,k)
            grd = grid%getgrd(i,j,k)
            pv = variable%getpv(i,j,k)
            dv = variable%getdv(i,j,k)
            tv = variable%gettv(i,j,k)
            dt = update%timestep%getdt(i,j,k)

            call update%lhs%setnorm(nx1,nx2,nx3,nx4,nx5,nx6)
            call update%lhs%setgrd(grd)
            call update%lhs%setpv(pv)
            call update%lhs%setdv(dv)
            call update%lhs%settv(tv)
            call update%lhs%setdt(dt)

            do l=1,update%npv
              res(l) = update%rhs%getres(l,i,j,k)
            end do
            pv = pv + update%lhs%getx(res)

            pv(1) = dmax1(-update%pref+1.d1,pv(1))
            pv(6) = dmin1(1.d0,dmax1(0.d0,pv(6)))
            pv(7) = dmin1(1.d0,dmax1(0.d0,pv(7)))
            if(update%l_turb) then
              pv(8) = dmax1(pv(8),update%kref)
              pv(9) = dmax1(pv(9),update%oref,update%rhs%getomega_cut(i,j,k))
            end if
            do n=1,update%npv
              call variable%setpv(n,i,j,k,pv(n))
            end do
            
            call eos%deteos(pv(1)+update%pref,pv(5),pv(6),pv(7),dv,tv)
            
            do n=1,update%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            do n=1,update%ntv
              call variable%settv(n,i,j,k,tv(n))
            end do
          end do
        end do
      end do
      
      if(update%l_eddy) then
        do k=2,update%kmax
          do j=2,update%jmax
            do i=2,update%imax
              nx1 = grid%getcx(i-1,j,k)
              nx2 = grid%getcx(i,j,k)
              nx3 = grid%getex(i,j-1,k)
              nx4 = grid%getex(i,j,k)
              nx5 = grid%gettx(i,j,k-1)
              nx6 = grid%gettx(i,j,k)
              grd = grid%getgrd(i,j,k)
              dv = variable%getdv(i,j,k)
              tv = variable%gettv(i,j,k)
              x(1,:) = variable%getpv(i-1,j,k)
              x(2,:) = variable%getpv(i,j,k)
              x(3,:) = variable%getpv(i+1,j,k)
              x(4,:) = variable%getpv(i,j-1,k)
              x(5,:) = variable%getpv(i,j+1,k)
              x(6,:) = variable%getpv(i,j,k-1)
              x(7,:) = variable%getpv(i,j,k+1)

              call update%eddy%setnorm(nx1,nx2,nx3,nx4,nx5,nx6)
              call update%eddy%setgrd(grd)
              call update%eddy%setpv(x)
              call update%eddy%setdv(dv)
              call update%eddy%settv(tv)
              tv(3) = update%eddy%caleddy()
              call variable%settv(3,i,j,k,tv(3))
            end do
          end do
        end do
      end if
    end subroutine eulerex
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine rk3rd(update,grid,variable,eos,nt_phy,nt)
      implicit none
      class(t_rk3rd), intent(inout) :: update
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(in) :: nt_phy,nt
      integer :: i,j,k,n,l,m
      real(8) :: nx1(3),nx2(3),nx3(3),nx4(3),nx5(3),nx6(3)
      real(8) :: pv(update%npv),dv(update%ndv),tv(update%ntv)
      real(8) :: grd(update%ngrd),dt
      real(8) :: x(7,update%npv)
      real(8) :: res(update%npv)

      
      
      do m=1,3
        call update%bc%setbc(grid,variable,eos)
        if(m.eq.1) call update%timestep%caltimestep(grid,variable,nt_phy,nt)
        call update%rhs%calrhs(grid,variable,eos)
        do k=2,update%kmax
          do j=2,update%jmax
            do i=2,update%imax
              if(m.eq.1) update%rk(:,i,j,k) = variable%getpv(i,j,k)
              nx1 = grid%getcx(i-1,j,k)
              nx2 = grid%getcx(i,j,k)
              nx3 = grid%getex(i,j-1,k)
              nx4 = grid%getex(i,j,k)
              nx5 = grid%gettx(i,j,k-1)
              nx6 = grid%gettx(i,j,k)
              grd = grid%getgrd(i,j,k)
              pv = variable%getpv(i,j,k)
              dv = variable%getdv(i,j,k)
              tv = variable%gettv(i,j,k)
              dt = update%timestep%getdt(i,j,k)
            
              call update%lhs%setnorm(nx1,nx2,nx3,nx4,nx5,nx6)
              call update%lhs%setgrd(grd)
              call update%lhs%setpv(pv)
              call update%lhs%setdv(dv)
              call update%lhs%settv(tv)
              call update%lhs%setdt(dt)

              do l=1,update%npv
                res(l) = update%rhs%getres(l,i,j,k)
              end do
            
              pv = update%a1(m)*update%rk(:,i,j,k) + update%a2(m)*pv + update%a3(m)*update%lhs%getx(res)
          
              pv(1) = dmax1(-update%pref+1.d1,pv(1))
              pv(6) = dmin1(1.d0,dmax1(0.d0,pv(6)))
              pv(7) = dmin1(1.d0,dmax1(0.d0,pv(7)))
              if(update%l_turb) then
                pv(8) = dmax1(pv(8),update%kref)
                pv(9) = dmax1(pv(9),update%oref,update%rhs%getomega_cut(i,j,k))
              end if
              do n=1,update%npv
                call variable%setpv(n,i,j,k,pv(n))
              end do
                     
              call eos%deteos(pv(1)+update%pref,pv(5),pv(6),pv(7),dv,tv)
            
              do n=1,update%ndv
                call variable%setdv(n,i,j,k,dv(n))
              end do
            
              do n=1,update%ntv
                call variable%settv(n,i,j,k,tv(n))
              end do
              
            end do
          end do
        end do
        if(update%l_eddy) then
          do k=2,update%kmax
            do j=2,update%jmax
              do i=2,update%imax
                nx1 = grid%getcx(i-1,j,k)
                nx2 = grid%getcx(i,j,k)
                nx3 = grid%getex(i,j-1,k)
                nx4 = grid%getex(i,j,k)
                nx5 = grid%gettx(i,j,k-1)
                nx6 = grid%gettx(i,j,k)
                grd = grid%getgrd(i,j,k)
                dv = variable%getdv(i,j,k)
                tv = variable%gettv(i,j,k)
                x(1,:) = variable%getpv(i-1,j,k)
                x(2,:) = variable%getpv(i,j,k)
                x(3,:) = variable%getpv(i+1,j,k)
                x(4,:) = variable%getpv(i,j-1,k)
                x(5,:) = variable%getpv(i,j+1,k)
                x(6,:) = variable%getpv(i,j,k-1)
                x(7,:) = variable%getpv(i,j,k+1)

                call update%eddy%setnorm(nx1,nx2,nx3,nx4,nx5,nx6)
                call update%eddy%setgrd(grd)
                call update%eddy%setpv(x)
                call update%eddy%setdv(dv)
                call update%eddy%settv(tv)
                tv(3) = update%eddy%caleddy()
                call variable%settv(3,i,j,k,tv(3))
              end do
            end do
          end do
        end if
      end do
      
    end subroutine rk3rd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine lusgs(update,grid,variable,eos,nt_phy,nt)
      implicit none
      class(t_lusgs), intent(inout) :: update
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(in) :: nt_phy,nt
      integer :: i,j,k,n,l
      real(8) :: nx1(3),nx2(3),nx3(3),nx4(3),nx5(3),nx6(3)
      real(8) :: pv(update%npv),dv(update%ndv),tv(update%ntv)
      real(8) :: grd(update%ngrd),dt,c(4),t(4)
      real(8) :: x(7,update%npv)
            
      call update%bc%setbc(grid,variable,eos)
      call update%timestep%caltimestep(grid,variable,nt_phy,nt)
      call update%rhs%calrhs(grid,variable,eos)
      
      update%dcv = 0.d0
      update%dqs = 0.d0
      do k=2,update%kmax
        do j=2,update%jmax
          do i=2,update%imax
            
            nx1 = grid%getcx(i-1,j,k)
            grd = grid%getgrd(i-1,j,k)
            pv = variable%getpv(i-1,j,k)
            dv = variable%getdv(i-1,j,k)
            tv = variable%gettv(i-1,j,k)
            call update%jac%setnorm(nx1)
            call update%jac%setgrd(grd)
            call update%jac%setpv(pv)
            call update%jac%setdv(dv)
            call update%jac%settv(tv)
            call update%jac%caljac(1)
            
            do n=1,update%npv
              update%dcv(n,i,j,k) = update%rhs%getres(n,i,j,k)
              do l=1,update%npv
                update%dcv(n,i,j,k) = update%dcv(n,i,j,k) + update%jac%geta(n,l)*update%dqs(l,i-1,j,k)
              end do
            end do
            
            nx1 = grid%getex(i,j-1,k)
            grd = grid%getgrd(i,j-1,k)
            pv = variable%getpv(i,j-1,k)
            dv = variable%getdv(i,j-1,k)
            tv = variable%gettv(i,j-1,k)          
            call update%jac%setnorm(nx1)
            call update%jac%setgrd(grd)
            call update%jac%setpv(pv)
            call update%jac%setdv(dv)
            call update%jac%settv(tv)
            call update%jac%caljac(1)
        
            do n=1,update%npv
              do l=1,update%npv
                update%dcv(n,i,j,k) = update%dcv(n,i,j,k) + update%jac%geta(n,l)*update%dqs(l,i,j-1,k)
              end do
            end do

            nx1 = grid%gettx(i,j,k-1)
            grd = grid%getgrd(i,j,k-1)
            pv = variable%getpv(i,j,k-1)
            dv = variable%getdv(i,j,k-1)
            tv = variable%gettv(i,j,k-1)          
            call update%jac%setnorm(nx1)
            call update%jac%setgrd(grd)
            call update%jac%setpv(pv)
            call update%jac%setdv(dv)
            call update%jac%settv(tv)
            call update%jac%caljac(1)
        
            do n=1,update%npv
              do l=1,update%npv
                update%dcv(n,i,j,k) = update%dcv(n,i,j,k) + update%jac%geta(n,l)*update%dqs(l,i,j,k-1)
              end do
            end do

            nx1 = grid%getcx(i-1,j,k)
            nx2 = grid%getcx(i,j,k)
            nx3 = grid%getex(i,j-1,k)
            nx4 = grid%getex(i,j,k)
            nx5 = grid%gettx(i,j,k-1)
            nx6 = grid%gettx(i,j,k)
            grd = grid%getgrd(i,j,k)
            pv = variable%getpv(i,j,k)
            dv = variable%getdv(i,j,k)
            tv = variable%gettv(i,j,k)
            dt = update%timestep%getdt(i,j,k)
            
            call update%lhs%setnorm(nx1,nx2,nx3,nx4,nx5,nx6)
            call update%lhs%setgrd(grd)
            call update%lhs%setpv(pv)
            call update%lhs%setdv(dv)
            call update%lhs%settv(tv)
            call update%lhs%setdt(dt)
            
            if(update%l_cav) then
              c = update%rhs%geticav(i,j,k)
              call update%lhs%setc(c)
            end if
            if(update%l_turb) then
              t = update%rhs%getitt(i,j,k)
              call update%lhs%sett(t)
            end if

            update%dqs(:,i,j,k) = update%lhs%getx(update%dcv(:,i,j,k))
          end do
        end do
      end do
      
      update%dqs = 0.d0
      do k=update%kmax,2,-1
        do j=update%jmax,2,-1
          do i=update%imax,2,-1
            
            nx1 = grid%getcx(i,j,k)
            grd = grid%getgrd(i+1,j,k)
            pv = variable%getpv(i+1,j,k)
            dv = variable%getdv(i+1,j,k)
            tv = variable%gettv(i+1,j,k)
            call update%jac%setnorm(nx1)
            call update%jac%setgrd(grd)
            call update%jac%setpv(pv)
            call update%jac%setdv(dv)
            call update%jac%settv(tv)
            call update%jac%caljac(-1)
        
            do n=1,update%npv
              do l=1,update%npv
                update%dcv(n,i,j,k) = update%dcv(n,i,j,k) - update%jac%geta(n,l)*update%dqs(l,i+1,j,k)
              end do
            end do
            
            nx1 = grid%getex(i,j,k)
            grd = grid%getgrd(i,j+1,k)
            pv = variable%getpv(i,j+1,k)
            dv = variable%getdv(i,j+1,k)
            tv = variable%gettv(i,j+1,k)          
            call update%jac%setnorm(nx1)
            call update%jac%setgrd(grd)
            call update%jac%setpv(pv)
            call update%jac%setdv(dv)
            call update%jac%settv(tv)
            call update%jac%caljac(-1)
            
            do n=1,update%npv
              do l=1,update%npv
                update%dcv(n,i,j,k) = update%dcv(n,i,j,k) - update%jac%geta(n,l)*update%dqs(l,i,j+1,k)
              end do
            end do

            nx1 = grid%gettx(i,j,k)
            grd = grid%getgrd(i,j,k+1)
            pv = variable%getpv(i,j,k+1)
            dv = variable%getdv(i,j,k+1)
            tv = variable%gettv(i,j,k+1)          
            call update%jac%setnorm(nx1)
            call update%jac%setgrd(grd)
            call update%jac%setpv(pv)
            call update%jac%setdv(dv)
            call update%jac%settv(tv)
            call update%jac%caljac(-1)
            
            do n=1,update%npv
              do l=1,update%npv
                update%dcv(n,i,j,k) = update%dcv(n,i,j,k) - update%jac%geta(n,l)*update%dqs(l,i,j,k+1)
              end do
            end do

            nx1 = grid%getcx(i-1,j,k)
            nx2 = grid%getcx(i,j,k)
            nx3 = grid%getex(i,j-1,k)
            nx4 = grid%getex(i,j,k)
            nx5 = grid%gettx(i,j,k-1)
            nx6 = grid%gettx(i,j,k)
            grd = grid%getgrd(i,j,k)
            pv = variable%getpv(i,j,k)
            dv = variable%getdv(i,j,k)
            tv = variable%gettv(i,j,k)
            dt = update%timestep%getdt(i,j,k)
            
            call update%lhs%setnorm(nx1,nx2,nx3,nx4,nx5,nx6)
            call update%lhs%setgrd(grd)
            call update%lhs%setpv(pv)
            call update%lhs%setdv(dv)
            call update%lhs%settv(tv)
            call update%lhs%setdt(dt)
            
            if(update%l_cav) then
              c = update%rhs%geticav(i,j,k)
              call update%lhs%setc(c)
            end if
            if(update%l_turb) then
              t = update%rhs%getitt(i,j,k)
              call update%lhs%sett(t)
            end if

            update%dqs(:,i,j,k) = update%lhs%getx(update%dcv(:,i,j,k))
          end do
        end do
      end do
      
      do k=2,update%kmax
        do j=2,update%jmax
          do i=2,update%imax
            pv = variable%getpv(i,j,k) + update%dqs(:,i,j,k)
            
            pv(1) = dmax1(-update%pref+1.d1,pv(1))
            pv(6) = dmin1(1.d0,dmax1(0.d0,pv(6)))
            pv(7) = dmin1(1.d0,dmax1(0.d0,pv(7)))
            if(update%l_turb) then
              pv(8) = dmax1(pv(8),update%kref)
              pv(9) = dmax1(pv(9),update%oref,update%rhs%getomega_cut(i,j,k))
            end if
            
            do n=1,update%npv
              call variable%setpv(n,i,j,k,pv(n))
            end do
            
            call eos%deteos(pv(1)+update%pref,pv(5),pv(6),pv(7),dv,tv)
            
            do n=1,update%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            do n=1,update%ntv
              call variable%settv(n,i,j,k,tv(n))
            end do
          end do
        end do
      end do
      if(update%l_eddy) then
        do k=2,update%kmax
          do j=2,update%jmax
            do i=2,update%imax
              nx1 = grid%getcx(i-1,j,k)
              nx2 = grid%getcx(i,j,k)
              nx3 = grid%getex(i,j-1,k)
              nx4 = grid%getex(i,j,k)
              nx5 = grid%gettx(i,j,k-1)
              nx6 = grid%gettx(i,j,k)
              grd = grid%getgrd(i,j,k)
              dv = variable%getdv(i,j,k)
              tv = variable%gettv(i,j,k)
              x(1,:) = variable%getpv(i-1,j,k)
              x(2,:) = variable%getpv(i,j,k)
              x(3,:) = variable%getpv(i+1,j,k)
              x(4,:) = variable%getpv(i,j-1,k)
              x(5,:) = variable%getpv(i,j+1,k)
              x(6,:) = variable%getpv(i,j,k-1)
              x(7,:) = variable%getpv(i,j,k+1)

              call update%eddy%setnorm(nx1,nx2,nx3,nx4,nx5,nx6)
              call update%eddy%setgrd(grd)
              call update%eddy%setpv(x)
              call update%eddy%setdv(dv)
              call update%eddy%settv(tv)
              tv(3) = update%eddy%caleddy()
              call variable%settv(3,i,j,k,tv(3))
            end do
          end do
        end do
      end if
    end subroutine lusgs
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module update_module
