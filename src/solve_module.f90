module solve_module
  use config_module
  use eos_module
  use grid_module
  use variable_module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! subordinate to update module  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  use initial_module
  use update_module
  use residual_module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  implicit none
  private
  public :: t_solve
  
  type t_solve
    private
    logical :: l_update,l_nsteady,l_ini
    integer :: npv,ndv,imax,jmax,kmax
    integer :: npmax,ntmax,nexport
    real(8) :: pref
    class(t_update), allocatable :: update
    class(t_ini), allocatable :: ini
    class(t_resi), allocatable :: resi
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: solve_equation
      procedure,private :: unsteadyupdate
  end type t_solve

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(solve,config,grid,eos)
      implicit none
      class(t_solve), intent(out) :: solve
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      class(t_eos), intent(in) :: eos

      solve%npv = config%getnpv()
      solve%ndv = config%getndv()
      solve%npmax = config%getnpmax()
      solve%ntmax = config%getntmax()
      solve%nexport = config%getnexport()
      
      solve%pref = config%getpref()

      solve%l_update = .false.
      solve%l_ini = .false.
      
      select case(config%gettimemethod())
      case(1)
        allocate(t_eulerex::solve%update)
      case(2)
        allocate(t_rk3rd::solve%update)
      case(3)
        allocate(t_lusgs::solve%update)
      end select
      
      select case(config%getiread())
      case(0)
        select case(config%getrotation())
        case(0)
          allocate(t_ini_initial::solve%ini)
        case default
          allocate(t_ini_initial_rot::solve%ini)
        end select
      case(1)
        allocate(t_ini_restart::solve%ini)
      end select
      
      select case(config%getnsteady())
      case(0)
        solve%l_nsteady = .false.
      case(1)
        solve%l_nsteady = .true.
      end select
      
      if(allocated(solve%update)) then
        call solve%update%construct(config,grid,eos)
        solve%l_update = .true.
      end if
      
      if(allocated(solve%ini)) then
        call solve%ini%construct(config,grid)
        solve%l_ini = .true.
      end if
      
      allocate(solve%resi)
      call solve%resi%construct(config,grid)
     
      solve%imax = grid%getimax()
      solve%jmax = grid%getjmax()
      solve%kmax = grid%getkmax()
      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(solve)
      implicit none
      class(t_solve), intent(inout) :: solve
      
      if(solve%l_update) then
        call solve%update%destruct()
        deallocate(solve%update)
      end if
      
      if(solve%l_ini) then
        call solve%ini%destruct()
        deallocate(solve%ini)
      end if
      
      call solve%resi%destruct()
      deallocate(solve%resi)
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine solve_equation(solve,grid,variable,eos)
      implicit none
      class(t_solve), intent(inout) :: solve
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer :: nt_phy,nt,nps,nts
      
      call solve%ini%initialize(grid,variable,eos,nps,nts)

      do nt_phy=nps,solve%npmax
        call solve%resi%setl_converge(.false.)
        do nt=nts,solve%ntmax
          call solve%resi%setqres(variable)
          call solve%update%timeinteg(grid,variable,eos,nt_phy,nt)
          call solve%resi%residual(variable,nt_phy,nt)        
          if(solve%resi%getconverge()) then
            if(.not.solve%l_nsteady) call variable%export_variable(nt_phy,nt)
            exit
          end if
          if((mod(nt,solve%nexport).eq.0).and.(.not.solve%l_nsteady)) then
            call variable%export_variable(nt_phy,nt)
          end if
        end do
        if(solve%l_nsteady) then
          call solve%unsteadyupdate(variable)
          if(mod(nt_phy,solve%nexport).eq.0) then
            call variable%export_variable(nt_phy,nt-1)
          end if
        end if
      end do
      
    end subroutine solve_equation
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine unsteadyupdate(solve,variable)
      implicit none
      class(t_solve), intent(inout) :: solve
      type(t_variable), intent(inout) :: variable
      integer :: i,j,k,n
      real(8) :: qq(solve%npv),dv(solve%ndv),pv(solve%npv)
      
      do k=2,solve%kmax
        do j=2,solve%jmax
          do i=2,solve%imax
            qq = variable%getqq(1,i,j,k)
            call variable%setqq(2,i,j,k,qq)
            dv = variable%getdv(i,j,k)
            pv = variable%getpv(i,j,k)
            qq(1) = dv(1)
            qq(2) = dv(1)*pv(2)
            qq(3) = dv(1)*pv(3)
            qq(4) = dv(1)*pv(4)
            qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-solve%pref
            do n=6,solve%npv
              qq(n) = dv(1)*pv(n)
            end do
            
            call variable%setqq(1,i,j,k,qq)
            
          end do
        end do
      end do
    end subroutine unsteadyupdate
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module solve_module
