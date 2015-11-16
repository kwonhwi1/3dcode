program main
  use config_module
  use eos_module
  use grid_module
  use variable_module
  use solve_module
  implicit none
  type(t_config) :: config
  class(t_eos), allocatable :: eos
  type(t_grid) :: grid
  type(t_variable) :: variable
  type(t_solve) :: solve
  real(8), dimension(:), allocatable :: dv,tv

  call config%construct()

  select case(config%getiturb())
  case(0,-1,-2)
    allocate(t_eos_prop::eos)
  case(-3)
    allocate(t_eosonly::eos)
  case default
  end select

  allocate(dv(config%getndv()),tv(config%getntv()))

  call eos%construct(config)
  call grid%construct(config)
  call variable%construct(config,grid)

  call eos%deteos(config%getpref(),config%gettref(),config%gety1ref(),config%gety2ref(),dv,tv)
  call config%setref(dv,tv)

  call solve%construct(config,grid,eos)

  call solve%solve_equation(grid,variable,eos)

  call solve%destruct()
  call variable%destruct()
  call grid%destruct()
  call eos%destruct()

  deallocate(eos,dv,tv)

  call config%destruct()
end program main
