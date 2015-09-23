program postprocess
  use config_module
  use eos_module
  use prop_module
  use postgrid_module
  use postvariable_module
  
  implicit none
  type(t_config) :: config
  type(t_eos) :: eos
  type(t_prop) :: prop
  type(t_grid) :: grid
  type(t_variable) :: variable
  
  integer :: io,istart,iend,nsolution,mode
  real(8) :: area

  open(newunit=io,file='./input_post.inp',status='old',action='read')
  read(io,*); read(io,*) istart
  read(io,*); read(io,*) iend
  read(io,*); read(io,*) nsolution
  read(io,*); read(io,*) mode
  read(io,*); read(io,*) area
  close(io)
      
  call config%construct(eos,prop) !eos & prop are also constructed
  call grid%construct(config)
  call variable%construct(config,grid,eos,istart,iend,nsolution)
  
  select case(mode)
  case(1)
    call variable%cgnswriting(config,grid)
  case(2)
    call variable%surface_writing(config,grid)
  case(3)
    call variable%clcd_writing(config,grid,area)
  case default
  end select
  
  call variable%destruct(grid)
  call grid%destruct()
  call config%destruct(eos,prop) !eos & prop are destructed
end program postprocess
