program main
  use config_module
  use eos_module
  use prop_module
  use grid_module
  use variable_module
  use solve_module
  implicit none
  type(t_config) :: config
  type(t_eos) :: eos
  type(t_prop) :: prop
  type(t_grid) :: grid
  type(t_variable) :: variable
  type(t_solve) :: solve

  call config%construct(eos,prop) !eos & prop are also constructed
  call grid%construct(config)
  call variable%construct(config,grid)
  call solve%construct(config,grid,variable,eos,prop)
  
  call solve%solve_equation(grid,variable,eos,prop)
  
  call solve%destruct()
  call variable%destruct()
  call grid%destruct()
  call config%destruct(eos,prop) !eos & prop are destructed
end program main   
