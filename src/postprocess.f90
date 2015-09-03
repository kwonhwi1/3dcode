program postprocess
  use config_module
  use eos_module
  use prop_module
  use postgrid_module
  use postvariable_module
  use datawriting_module
  
  implicit none
  type(t_config) :: config
  type(t_eos) :: eos
  type(t_prop) :: prop
  type(t_grid) :: grid
  type(t_variable) :: variable
  type(t_datawriting) :: datawriting
  
  integer :: io,istart,iend,nsolution
  real(8) :: area,zposition

  open(newunit=io,file='./input_post.inp',status='old',action='read')
  read(io,*); read(io,*) istart
  read(io,*); read(io,*) iend
  read(io,*); read(io,*) nsolution
  read(io,*); read(io,*) mode
  read(io,*); read(io,*) area
  close(io)
      
  call config%construct(eos,prop) !eos & prop are also constructed
  call grid%construct(config)
  call variable%construct(config,eos,istart,iend,nsolution)
  
  select case(mode)
  case(1)
    call datawriting%cgnswriting(config,variable)
  case(2)
    call datawriting%surface_writing(config,grid,variable)
  case(3)
    call datawriting%clcd_writing(config,grid,variable,area)
  case default
  end select
  
  call variable%destruct()
  call grid%destruct()
  call config%destruct(eos,prop) !eos & prop are destructed
end program postprocess
