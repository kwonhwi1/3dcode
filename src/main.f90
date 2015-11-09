program main
  use config_module
  use eos_module
  use grid_module
  use variable_module
  use solve_module
  implicit none
  type(t_config) :: config
  type(t_eos) :: eos
  type(t_grid) :: grid
  type(t_variable) :: variable
  type(t_solve) :: solve

  call config%construct(eos) !eos is also constructed
  call grid%construct(config)
  call variable%construct(config,grid)
  call solve%construct(config,grid,variable,eos)
  
  call solve%solve_equation(grid,variable,eos)
  
  call solve%destruct()
  call variable%destruct()
  call grid%destruct()
  call config%destruct(eos) !eos is destructed
end program main   
