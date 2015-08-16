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
  type(t_grod) :: grid
  type(t_variable) :: variable
  type(t_datawriting) :: datawriting
  
  integer :: io,istart,iend,nsolution
  
  open(newunit=io,file='./input_post.inp',status='old',action='read')
  read(io,*); read(io,*) istart
  read(io,*); read(io,*) iend
  read(io,*); read(io,*) nsolution
  close(io)
      
  call config%construct(eos,prop) !eos & prop are also constructed
  call grid%construct(config)
  call variable%construct(config,eos,istart,iend,nsolution)
      
  call datawriting%cgnswriting(config,variable)
  call datawriting%clcd_writing(config,grid,variable)
  
  call variable%destruct()
  call grid%destruct()
  call config%destruct(eos,prop) !eos & prop are destructed
end program postprocess
