program postprocess
  use config_module
  use eos_module
  use postgrid_module
  use postvariable_module

  implicit none
  type(t_config) :: config
  class(t_eos), allocatable :: eos
  type(t_grid) :: grid
  type(t_variable) :: variable
  real(8), dimension(:), allocatable :: dv,tv

  integer :: io,istart,iend,nsolution,mode,nsolname
  real(8) :: area
  integer :: l
  character(len=4) :: extravar_tmp
  character(len=:), allocatable :: extravar

  open(newunit=io,file='./input_post.inp',status='old',action='read')
    read(io,*); read(io,*) istart
    read(io,*); read(io,*) iend
    read(io,*); read(io,*) nsolution
    read(io,*); read(io,*) mode
    read(io,*); read(io,*) area
    read(io,*); read(io,*) nsolname
    if(mode.eq.1) read(io,*); read(io,*) extravar_tmp
  close(io)

  call config%construct()

  extravar = trim(extravar_tmp)
  if(extravar.ne.'0') then
    do l=1,len(extravar)
      if(extravar(l:l).eq.'3') then
        if(config%getncav().eq.0) write(*,*) '!!! wrong extravar; check ncav !!!'
      end if
      if(extravar(l:l).eq.'4') then
        if(config%getcsf().eq.0) write(*,*) '!!! wrong extravar; check csf !!!'
      end if
   end do
  end if

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
  call variable%construct(config,grid,eos,istart,iend,nsolution)

  call eos%deteos(config%getpref(),config%gettref(),config%gety1ref(),config%gety2ref(),dv,tv)
  call config%setref(dv,tv)

  select case(mode)
  case(1)
    call variable%cgnswriting(config,grid,nsolname,extravar)
  case(2)
    call variable%surface_writing(config,grid)
  case(3)
    call variable%clcd_writing(config,grid,area)
  case default
  end select

  call variable%destruct(grid)
  call grid%destruct()
  call eos%destruct()

  deallocate(eos,dv,tv)

  call config%destruct()
end program postprocess
