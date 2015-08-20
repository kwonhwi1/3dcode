module config_module
  use eos_module
  use prop_module
  implicit none
  private
  public :: t_config
  
  type t_config
    private
    logical :: status
    character(30) :: name
    integer :: rank,size,stencil
    integer :: iread,rstnum
    integer :: nsteady,npmax,ntmax,bond
    real(8) :: dt_phy
    integer :: nexport,nprint
    integer :: iturb,tcomp
    integer :: nscheme,precd,nmuscl,nlim
    integer :: timemethod,local,prec
    real(8) :: cfl
    integer :: fluid,fluid_eostype,ngas,gas_eostype,mixingrule
    integer :: ncav
    real(8) :: c_v,c_c
    real(8) :: pref,uref,aoa,tref,y1ref,y2ref
    real(8) :: l_chord,l_character,scale
    real(8) :: str,pi
    real(8) :: dvref(18),tvref(2)
    real(8) :: kref,oref,emutref
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: getrank
      procedure :: getsize
      procedure :: getstencil
      procedure :: getname
      procedure :: getiread
      procedure :: getnsteady
      procedure :: getnpmax
      procedure :: getntmax
      procedure :: getbond
      procedure :: getdt_phy
      procedure :: getnexport
      procedure :: getnprint
      procedure :: getrstnum
      procedure :: getiturb
      procedure :: getnscheme
      procedure :: getprecd
      procedure :: getnmuscl
      procedure :: getnlim
      procedure :: gettimemethod
      procedure :: getlocal
      procedure :: getprec
      procedure :: getncav
      procedure :: gettcomp
      procedure :: getc_v
      procedure :: getc_c
      procedure :: getcfl
      procedure :: getpref
      procedure :: geturef
      procedure :: getaoa
      procedure :: gettref
      procedure :: gety1ref
      procedure :: gety2ref
      procedure :: getrhoref
      procedure :: getvisref
      procedure :: getkref
      procedure :: getoref
      procedure :: getemutref
      procedure :: getl_chord
      procedure :: getl_character
      procedure :: getscale
      procedure :: getstr
      procedure :: getpi
  end type t_config
    
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(config,eos,prop)
      implicit none
      include 'mpif.h'
      class(t_config), intent(out) :: config
      type(t_eos), intent(out) :: eos
      type(t_prop), intent(out) :: prop
      integer :: n,io,ierr
      
      config%stencil = 38
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world,config%size,ierr)
      call mpi_comm_rank(mpi_comm_world,config%rank,ierr)
    
      do n = 0,config%size-1
        if(n.eq.config%rank) then
          open(newunit=io,file='./input.inp',status='old',action='read')
          read(io,*); read(io,*) config%name
          read(io,*); read(io,*) config%iread,config%rstnum
          read(io,*); read(io,*) config%nsteady,config%npmax,config%ntmax,config%bond,config%dt_phy
          read(io,*); read(io,*) config%nexport,config%nprint
          read(io,*); read(io,*) config%iturb,config%tcomp
          read(io,*); read(io,*) config%nscheme,config%precd,config%nmuscl,config%nlim
          read(io,*); read(io,*) config%timemethod,config%local,config%prec,config%cfl
          read(io,*); read(io,*) config%fluid,config%fluid_eostype,config%ngas,config%gas_eostype,config%mixingrule
          read(io,*); read(io,*) config%ncav,config%c_v,config%c_c
          read(io,*); read(io,*) config%pref,config%uref,config%aoa,config%tref,config%y1ref,config%y2ref
          read(io,*); read(io,*) config%l_chord,config%l_character,config%scale
          close(io)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
      end do
      
      call eos%construct(config%fluid,config%fluid_eostype,config%ngas,config%gas_eostype,config%mixingrule)
      call prop%construct(config%fluid,config%fluid_eostype,config%ngas)
      
      
      config%pi  = datan(1.d0)*4.d0
      config%aoa = config%aoa*config%pi/180.d0
      config%str = config%l_character/config%pi/config%dt_phy/config%uref
      
      call eos%deteos(config%pref,config%tref,config%y1ref,config%y2ref,config%dvref)
      call prop%detprop(config%dvref(3),config%dvref(4),config%dvref(5),config%tref,config%y1ref,config%y2ref,config%tvref)
      
      
      if(config%iturb.eq.0) then
        config%oref = 10.d0*config%uref/config%l_character
        config%emutref = 10.d0**(-5)*config%tvref(1)
        config%kref = config%emutref/config%dvref(1)*config%oref
      else if(config%iturb.eq.-1) then
        config%kref = 1.5d0*(0.0001d0*config%uref)**2
        config%emutref =  10.d0**(-5)*config%tvref(1)
        config%oref = 0.09d0*config%dvref(1)*config%kref**2/config%emutref
      else
        config%kref = 0.d0
        config%oref = 0.d0
        config%emutref =  0.d0    
      end if
      
      config%status = .true.
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(config,eos,prop)
      implicit none
      include 'mpif.h'
      class(t_config), intent(inout) :: config
      type(t_eos), intent(inout) :: eos
      type(t_prop), intent(inout) :: prop
      integer :: ierr
      
      config%status = .false.      
      call eos%destruct()
      call prop%destruct()
      
      call mpi_finalize(ierr)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    pure function getrank(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getrank
      
      getrank = config%rank
      
    end function getrank
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getsize(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getsize
      
      getsize = config%size
      
    end function getsize
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getstencil(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getstencil

      getstencil = config%stencil

    end function getstencil
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getname(config)
      implicit none
      class(t_config), intent(in) :: config
      character(30) :: getname
      
      getname = config%name
      
    end function getname
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
    pure function getiread(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getiread
      
      getiread = config%iread
      
    end function getiread
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getnsteady(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getnsteady
      
      getnsteady = config%nsteady
      
    end function getnsteady
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getnpmax(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getnpmax
      
      getnpmax = config%npmax
      
    end function getnpmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getntmax(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getntmax
      
      getntmax = config%ntmax
      
    end function getntmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getbond(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getbond
      
      getbond = config%bond
      
    end function getbond
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getdt_phy(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getdt_phy
      
      getdt_phy = config%dt_phy
      
    end function getdt_phy
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getnexport(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getnexport
      
      getnexport = config%nexport
      
    end function getnexport
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getnprint(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getnprint
      
      getnprint = config%nprint
      
    end function getnprint
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getrstnum(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getrstnum
      
      getrstnum = config%rstnum
      
    end function getrstnum
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getiturb(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getiturb
      
      getiturb = config%iturb
      
    end function getiturb
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getnscheme(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getnscheme
      
      getnscheme = config%nscheme
      
    end function getnscheme
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getprecd(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getprecd
      
      getprecd = config%precd
      
    end function getprecd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getnmuscl(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getnmuscl
      
      getnmuscl = config%nmuscl
      
    end function getnmuscl
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getnlim(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getnlim
      
      getnlim = config%nlim
    end function getnlim
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function gettimemethod(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: gettimemethod
      
      gettimemethod = config%timemethod
      
    end function gettimemethod
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getlocal(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getlocal
      
      getlocal = config%local
      
    end function getlocal
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getprec(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getprec
      
      getprec = config%prec
      
    end function getprec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getncav(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: getncav
      
      getncav = config%ncav
      
    end function getncav
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function gettcomp(config)
      implicit none
      class(t_config), intent(in) :: config
      integer :: gettcomp
      
      gettcomp = config%tcomp
      
    end function gettcomp
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getc_v(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getc_v
      
      getc_v = config%c_v
      
    end function getc_v
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getc_c(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getc_c
      
      getc_c = config%c_c
      
    end function getc_c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getcfl(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getcfl
      
      getcfl = config%cfl
      
    end function getcfl
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getpref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getpref
      
      getpref = config%pref
      
    end function getpref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function geturef(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: geturef
      
      geturef = config%uref
      
    end function geturef
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getaoa(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getaoa
      
      getaoa = config%aoa
      
    end function getaoa
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function gettref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: gettref
      
      gettref = config%tref
      
    end function gettref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function gety1ref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: gety1ref
      
      gety1ref = config%y1ref
    
    end function gety1ref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function gety2ref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: gety2ref
      
      gety2ref = config%y2ref
      
    end function gety2ref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getrhoref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getrhoref
      
      getrhoref = config%dvref(1)
      
    end function getrhoref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getvisref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getvisref
      
      getvisref = config%tvref(1)
      
    end function getvisref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getkref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getkref
      
      getkref = config%kref
      
    end function getkref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getoref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getoref
      
      getoref = config%oref
      
    end function getoref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getemutref(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getemutref
      
      getemutref = config%emutref
      
    end function getemutref
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getl_chord(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getl_chord
      
      getl_chord = config%l_chord
      
    end function getl_chord
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getl_character(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getl_character
      
      getl_character = config%l_character
      
    end function getl_character
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getscale(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getscale
      
      getscale = config%scale
      
    end function getscale
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getstr(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getstr
      
      getstr = config%str
      
    end function getstr
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getpi(config)
      implicit none
      class(t_config), intent(in) :: config
      real(8) :: getpi
      
      getpi = config%pi
      
    end function getpi
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module config_module
