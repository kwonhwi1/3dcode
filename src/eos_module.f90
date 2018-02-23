module eos_module
  use mpi
  use config_module
  implicit none
  private
  public :: t_eos,t_eosonly,t_eos_prop

  type t_iapws_coeff_eos
    integer, dimension(:), allocatable :: i,j,ir,jr
    real(8), dimension(:), allocatable :: n,nr
  end type t_iapws_coeff_eos

  type t_iapws_coeff_prop
    real(8), dimension(:), allocatable :: h,l
    real(8), dimension(:,:), allocatable :: hh,ll
  end type t_iapws_coeff_prop

  type t_srk_property
    real(8) :: tc,pc,tb,mw,w,k,cv0 ! fluid properties
    real(8) :: s,a,b ! coefficients for srk
  end type t_srk_property

  type t_db_data
    real(8) :: t
    integer :: p_nstart,p_nend,p_nsub
    integer, allocatable :: p_iboundary(:)
    real(8), allocatable :: p_boundary(:),p_delta(:)
    real(8), allocatable :: p(:)
    real(8), allocatable :: rho(:,:),h(:,:),vis(:,:),cond(:,:)
  end type t_db_data

  type t_db_const
    integer :: t_ndata,t_nsub
    integer, allocatable :: t_iboundary(:)
    real(8), allocatable :: t_boundary(:),t_delta(:)
    type(t_db_data), allocatable :: db(:)
  end type t_db_const

  type, abstract :: t_eos
    private
    logical :: prop
    integer :: ndv,ntv,rank,size
    real(8) :: gamma_f,gamma_g
    procedure(p_pww), public, pointer :: get_pww
    procedure(p_tww), public, pointer :: get_tww
    procedure(p_sigma), public, pointer :: get_sigma ! sigma of fluid with its vapor
    procedure(p_eos_l), pointer :: eos_l
    procedure(p_eos_v), pointer :: eos_v,eos_vs
    procedure(p_eos_g), pointer :: eos_g,eos_gs
    procedure(p_eos_l_prop), pointer :: eos_l_prop
    procedure(p_eos_v_prop), pointer :: eos_v_prop,eos_vs_prop
    procedure(p_eos_g_prop), pointer :: eos_g_prop,eos_gs_prop
    type(t_iapws_coeff_eos)  :: iapws_liquid_coeff,iapws_vapor_coeff,iapws_m_vapor_coeff
    type(t_iapws_coeff_prop) :: iapws_prop_coeff
    type(t_srk_property) :: srk_vapor_property,srk_gas_property
    type(t_db_const) :: db_const(3) ! 1=liquid,2=vapor,3=gas
    contains
      procedure(p_construct), deferred :: construct
      procedure :: destruct
      procedure :: deteos_simple
      procedure :: getgamma_f
      procedure :: getgamma_g
      procedure, private :: set_iapws97_eos
      procedure, private :: set_iapws97_prop
      procedure, private :: set_srk_property
      procedure, private :: set_database
      procedure(p_deteos), deferred :: deteos
  end type t_eos


  abstract interface
    subroutine p_deteos(eos,p,t,y1,y2,dv,tv)
      import t_eos
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(eos%ndv) ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      real(8), intent(inout) :: tv(eos%ntv)
    end subroutine p_deteos

    subroutine p_construct(eos,config)
      import t_config
      import t_eos
      implicit none
      class(t_eos), intent(inout) :: eos
      type(t_config), intent(in) :: config
    end subroutine p_construct
  end interface

  type, extends(t_eos) :: t_eosonly
    contains
      procedure :: construct => construct_eos
      procedure :: deteos => deteosonly
  end type t_eosonly

  type, extends(t_eos) :: t_eos_prop
    contains
      procedure :: construct => construct_eos_prop
      procedure :: deteos => deteos_prop
  end type t_eos_prop

  type t_eos2
    real(8) :: rho,h,drdp,drdt,dhdp,dhdt
  end type t_eos2

  type t_prop2
    real(8) :: vis,cond
  end type t_prop2

  interface
    function p_pww(eos,t) result(pww)
      import t_eos
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
    end function p_pww

    function p_tww(eos,p) result(tww)
      import t_eos
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p
      real(8) :: tww
    end function p_tww

   function p_sigma(eos,t) result(sigma)
      import t_eos
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: sigma
    end function p_sigma

    subroutine p_eos_l(eos,p,t,phase,eos2)
      import t_eos
      import t_eos2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
    end subroutine p_eos_l

    subroutine p_eos_v(eos,p,t,phase,eos2)
      import t_eos
      import t_eos2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
    end subroutine p_eos_v

    subroutine p_eos_g(eos,p,t,phase,eos2)
      import t_eos
      import t_eos2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
    end subroutine p_eos_g

    subroutine p_eos_l_prop(eos,p,t,phase,eos2,prop2)
      import t_eos
      import t_eos2
      import t_prop2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
    end subroutine p_eos_l_prop

    subroutine p_eos_v_prop(eos,p,t,phase,eos2,prop2)
      import t_eos
      import t_eos2
      import t_prop2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
    end subroutine p_eos_v_prop

    subroutine p_eos_g_prop(eos,p,t,phase,eos2,prop2)
      import t_eos
      import t_eos2
      import t_prop2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
    end subroutine p_eos_g_prop
  end interface

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct_eos(eos,config)
      implicit none
      class(t_eosonly), intent(inout) :: eos
      type(t_config), intent(in) :: config

      eos%ndv = config%getndv()
      eos%ntv = config%getntv()
      !eos%size  = config%getsize()
      !eos%rank  = config%getrank()

      select case(config%getfluid())
      case(1) ! stiffened for water
        eos%eos_l  => stiffened
        eos%eos_v  => ideal
        eos%get_pww => h2o_pww
        eos%get_tww => h2o_tww
        eos%get_sigma => h2o_sigma
        eos%gamma_f = 4.4d0
        eos%eos_vs => null()
      case(2) ! iapws97 for water
        call eos%set_iapws97_eos()
        eos%eos_l  => iapws97_l
        eos%eos_v  => iapws97_v
        eos%get_pww => h2o_pww
        eos%get_tww => h2o_tww
        eos%get_sigma => h2o_sigma
        eos%gamma_f = 1.33d0
        eos%eos_vs => null()
      case(3)! database water
        call eos%set_database(1,1)
        call eos%set_database(1,2)
        eos%eos_l  => database
        eos%eos_v  => database
        eos%get_pww => h2o_pww
        eos%get_tww => h2o_tww
        eos%get_sigma => h2o_sigma
        eos%gamma_f = 1.33d0
        eos%eos_vs => null()
      case(4)! database nitrogen
        call eos%set_database(2,1)
        call eos%set_database(2,2)
        eos%eos_l  => database
        eos%eos_v  => database
        eos%get_pww => n2_pww
        eos%get_tww => n2_tww
        eos%get_sigma => n2_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(2,eos%srk_vapor_property)
        eos%eos_vs  => srk_g
      case(5)! database oxygen
        call eos%set_database(3,1)
        call eos%set_database(3,2)
        eos%eos_l  => database
        eos%eos_v  => database
        eos%get_pww => o2_pww
        eos%get_tww => o2_tww
        eos%get_sigma => o2_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(3,eos%srk_vapor_property)
        eos%eos_vs  => srk_g
      case(6)! database hydrogen
        call eos%set_database(4,1)
        call eos%set_database(4,2)
        eos%eos_l  => database
        eos%eos_v  => database
        eos%get_pww => h2_pww
        eos%get_tww => h2_tww
        eos%get_sigma => h2_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(4,eos%srk_vapor_property)
        eos%eos_vs  => srk_g
      case(7)! database DME
        call eos%set_database(7,1)
        call eos%set_database(7,2)
        eos%eos_l  => database
        eos%eos_v  => database
        eos%get_pww => DME_pww
        eos%get_tww => DME_tww
        eos%get_sigma => DME_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(4,eos%srk_vapor_property)
        eos%eos_vs  => srk_g
      case default
      end select

      select case(config%getngas())
      case(1) ! ideal gas
        eos%eos_g  => ideal
        eos%gamma_g = 1.4d0
        eos%eos_gs => null()
      case(2,3,4) ! nitrogen,oxygen,hydrogen
        call eos%set_database(config%getngas(),3)
        eos%eos_g  => database
        eos%gamma_g = 1.4d0
        call eos%set_srk_property(config%getngas(),eos%srk_gas_property)
        eos%eos_gs => srk_g
      case(5) ! helium
        call eos%set_database(config%getngas(),3)
        eos%eos_g  => database
        eos%gamma_g = 1.66d0
        call eos%set_srk_property(config%getngas(),eos%srk_gas_property)
        eos%eos_gs => srk_g
     case(6,7,8) ! nitrogen,oxygen,hydrogen
        call eos%set_srk_property(config%getngas()-4,eos%srk_gas_property)
        eos%eos_g  => srk_g
        eos%gamma_g = 1.4d0
        eos%eos_gs => null()
      case(9) ! helium
        call eos%set_srk_property(config%getngas()-4,eos%srk_gas_property)
        eos%eos_g  => srk_g
        eos%gamma_g = 1.66d0
        eos%eos_gs => null()
      case default
      end select

      eos%eos_l_prop => null()
      eos%eos_v_prop => null()
      eos%eos_g_prop => null()
      eos%eos_vs_prop => null()
      eos%eos_gs_prop => null()

    end subroutine construct_eos
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct_eos_prop(eos,config)
      implicit none
      class(t_eos_prop), intent(inout) :: eos
      type(t_config), intent(in) :: config

      eos%ndv = config%getndv()
      eos%ntv = config%getntv()
      !eos%size  = config%getsize()
      !eos%rank  = config%getrank()

      select case(config%getfluid())
      case(1) ! stiffened for water
        eos%eos_l  => stiffened
        eos%eos_v  => ideal
        eos%eos_l_prop  => stiffened_prop
        eos%eos_v_prop  => ideal_prop
        eos%get_pww => h2o_pww
        eos%get_tww => h2o_tww
        eos%get_sigma => h2o_sigma
        eos%gamma_f = 4.4d0
        eos%eos_vs => null()
        eos%eos_vs_prop => null()
      case(2) ! iapws97 for water
        call eos%set_iapws97_eos()
        call eos%set_iapws97_prop()
        eos%eos_l  => iapws97_l
        eos%eos_v  => iapws97_v
        eos%eos_l_prop  => iapws97_l_prop
        eos%eos_v_prop  => iapws97_v_prop
        eos%get_pww => h2o_pww
        eos%get_tww => h2o_tww
        eos%get_sigma => h2o_sigma
        eos%gamma_f = 1.33d0
        eos%eos_vs => null()
        eos%eos_vs_prop => null()
      case(3)! database water
        call eos%set_database(1,1)
        call eos%set_database(1,2)
        eos%eos_l  => database
        eos%eos_l_prop  => database_prop
        eos%eos_v  => database
        eos%eos_v_prop  => database_prop
        eos%get_pww => h2o_pww
        eos%get_tww => h2o_tww
        eos%get_sigma => h2o_sigma
        eos%gamma_f = 1.33d0
        eos%eos_vs => null()
        eos%eos_vs_prop => null()
      case(4)! database nitrogen
        call eos%set_database(2,1)
        call eos%set_database(2,2)
        eos%eos_l  => database
        eos%eos_l_prop  => database_prop
        eos%eos_v  => database
        eos%eos_v_prop  => database_prop
        eos%get_pww => n2_pww
        eos%get_tww => n2_tww
        eos%get_sigma => n2_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(2,eos%srk_vapor_property)
        eos%eos_vs => srk_g
        eos%eos_vs_prop  => srk_g_prop
     case(5)! database oxygen
        call eos%set_database(3,1)
        call eos%set_database(3,2)
        eos%eos_l  => database
        eos%eos_l_prop  => database_prop
        eos%eos_v  => database
        eos%eos_v_prop  => database_prop
        eos%get_pww => o2_pww
        eos%get_tww => o2_tww
        eos%get_sigma => o2_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(3,eos%srk_vapor_property)
        eos%eos_vs => srk_g
        eos%eos_vs_prop  => srk_g_prop
     case(6)! database hydrogen
        call eos%set_database(4,1)
        call eos%set_database(4,2)
        eos%eos_l  => database
        eos%eos_l_prop  => database_prop
        eos%eos_v  => database
        eos%eos_v_prop  => database_prop
        eos%get_pww => h2_pww
        eos%get_tww => h2_tww
        eos%get_sigma => h2_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(4,eos%srk_vapor_property)
        eos%eos_vs => srk_g
        eos%eos_vs_prop  => srk_g_prop
     case(7)! database DME
        call eos%set_database(7,1)
        call eos%set_database(7,2)
        eos%eos_l  => database
        eos%eos_l_prop  => database_prop
        eos%eos_v  => database
        eos%eos_v_prop  => database_prop
        eos%get_pww => DME_pww
        eos%get_tww => DME_tww
        eos%get_sigma => DME_sigma
        eos%gamma_f = 1.4d0
        call eos%set_srk_property(7,eos%srk_vapor_property)
        eos%eos_vs => srk_g
        eos%eos_vs_prop  => srk_g_prop
      case default
      end select

      select case(config%getngas())
      case(1) ! ideal gas
        eos%eos_g  => ideal
        eos%eos_g_prop  => ideal_prop
        eos%gamma_g = 1.4d0
        eos%eos_gs => null()
        eos%eos_gs_prop => null()
      case(2,3,4) ! nitrogen,oxygen,hydrogen
        call eos%set_database(config%getngas(),3)
        eos%eos_g  => database
        eos%eos_g_prop  => database_prop
        eos%gamma_g = 1.4d0
        call eos%set_srk_property(config%getngas(),eos%srk_gas_property)
        eos%eos_gs => srk_g
        eos%eos_gs_prop => srk_g_prop
     case(5) ! helium
        call eos%set_database(config%getngas(),3)
        eos%eos_g  => database
        eos%eos_g_prop  => database_prop
        eos%gamma_g = 1.66d0
        call eos%set_srk_property(config%getngas(),eos%srk_gas_property)
        eos%eos_gs => srk_g
        eos%eos_gs_prop => srk_g_prop
      case(6,7,8) ! nitrogen,oxygen,hydrogen
        call eos%set_srk_property(config%getngas()-4,eos%srk_gas_property)
        call eos%set_database(config%getngas()-4,3)
        eos%eos_g  => srk_g
        eos%eos_g_prop  => srk_g_prop
        eos%gamma_g = 1.4d0
        eos%eos_gs => null()
        eos%eos_gs_prop => null()
      case(9) ! helium
        call eos%set_srk_property(config%getngas()-4,eos%srk_gas_property)
        call eos%set_database(config%getngas()-4,3)
        eos%eos_g  => srk_g
        eos%eos_g_prop  => srk_g_prop
        eos%gamma_g = 1.66d0
        eos%eos_gs => null()
        eos%eos_gs_prop => null()
      case default
      end select

    end subroutine construct_eos_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(eos)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer :: m,n

      if(associated(eos%get_pww))    nullify(eos%get_pww)
      if(associated(eos%get_tww))    nullify(eos%get_tww)
      if(associated(eos%get_sigma))  nullify(eos%get_sigma)
      if(associated(eos%eos_l))      nullify(eos%eos_l)
      if(associated(eos%eos_v))      nullify(eos%eos_v)
      if(associated(eos%eos_g))      nullify(eos%eos_g)
      if(associated(eos%eos_l_prop)) nullify(eos%eos_l_prop)
      if(associated(eos%eos_v_prop)) nullify(eos%eos_v_prop)
      if(associated(eos%eos_g_prop)) nullify(eos%eos_g_prop)
      if(associated(eos%eos_vs))     nullify(eos%eos_vs)
      if(associated(eos%eos_gs))     nullify(eos%eos_gs)
      if(associated(eos%eos_vs_prop)) nullify(eos%eos_vs_prop)
      if(associated(eos%eos_gs_prop)) nullify(eos%eos_gs_prop)

      do n=1,3
        if(allocated(eos%db_const(n)%db)) then
          do m=1,eos%db_const(n)%t_ndata
            if(allocated(eos%db_const(n)%db(m)%p_boundary))  deallocate(eos%db_const(n)%db(m)%p_boundary)
            if(allocated(eos%db_const(n)%db(m)%p_iboundary)) deallocate(eos%db_const(n)%db(m)%p_iboundary)
            if(allocated(eos%db_const(n)%db(m)%p_delta))     deallocate(eos%db_const(n)%db(m)%p_delta)
            if(allocated(eos%db_const(n)%db(m)%p))           deallocate(eos%db_const(n)%db(m)%p)
            if(allocated(eos%db_const(n)%db(m)%rho))         deallocate(eos%db_const(n)%db(m)%rho)
            if(allocated(eos%db_const(n)%db(m)%h))           deallocate(eos%db_const(n)%db(m)%h)
            if(allocated(eos%db_const(n)%db(m)%vis))         deallocate(eos%db_const(n)%db(m)%vis)
            if(allocated(eos%db_const(n)%db(m)%cond))        deallocate(eos%db_const(n)%db(m)%cond)
          end do
          if(allocated(eos%db_const(n)%t_boundary))  deallocate(eos%db_const(n)%t_boundary)
          if(allocated(eos%db_const(n)%t_iboundary)) deallocate(eos%db_const(n)%t_iboundary)
          if(allocated(eos%db_const(n)%t_delta))     deallocate(eos%db_const(n)%t_delta)
          deallocate(eos%db_const(n)%db)
        end if
      end do

      if(allocated(eos%iapws_liquid_coeff%i))   deallocate(eos%iapws_liquid_coeff%i)
      if(allocated(eos%iapws_liquid_coeff%j))   deallocate(eos%iapws_liquid_coeff%j)
      if(allocated(eos%iapws_liquid_coeff%n))   deallocate(eos%iapws_liquid_coeff%n)
      if(allocated(eos%iapws_vapor_coeff%j))    deallocate(eos%iapws_vapor_coeff%j)
      if(allocated(eos%iapws_vapor_coeff%n))    deallocate(eos%iapws_vapor_coeff%n)
      if(allocated(eos%iapws_vapor_coeff%ir))   deallocate(eos%iapws_vapor_coeff%ir)
      if(allocated(eos%iapws_vapor_coeff%jr))   deallocate(eos%iapws_vapor_coeff%jr)
      if(allocated(eos%iapws_vapor_coeff%nr))   deallocate(eos%iapws_vapor_coeff%nr)
      if(allocated(eos%iapws_m_vapor_coeff%j))  deallocate(eos%iapws_m_vapor_coeff%j)
      if(allocated(eos%iapws_m_vapor_coeff%n))  deallocate(eos%iapws_m_vapor_coeff%n)
      if(allocated(eos%iapws_m_vapor_coeff%ir)) deallocate(eos%iapws_m_vapor_coeff%ir)
      if(allocated(eos%iapws_m_vapor_coeff%jr)) deallocate(eos%iapws_m_vapor_coeff%jr)
      if(allocated(eos%iapws_m_vapor_coeff%nr)) deallocate(eos%iapws_m_vapor_coeff%nr)
      if(allocated(eos%iapws_prop_coeff%h))     deallocate(eos%iapws_prop_coeff%h)
      if(allocated(eos%iapws_prop_coeff%l))     deallocate(eos%iapws_prop_coeff%l)
      if(allocated(eos%iapws_prop_coeff%hh))    deallocate(eos%iapws_prop_coeff%hh)
      if(allocated(eos%iapws_prop_coeff%ll))    deallocate(eos%iapws_prop_coeff%ll)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_iapws97_eos(eos)
      implicit none
      class(t_eos), intent(inout) :: eos

      allocate(eos%iapws_liquid_coeff%i(34),eos%iapws_liquid_coeff%j(34),eos%iapws_liquid_coeff%n(34))

      allocate(eos%iapws_vapor_coeff%j(9),eos%iapws_vapor_coeff%n(9))
      allocate(eos%iapws_vapor_coeff%ir(43),eos%iapws_vapor_coeff%jr(43),eos%iapws_vapor_coeff%nr(43))

      allocate(eos%iapws_m_vapor_coeff%j(9),eos%iapws_m_vapor_coeff%n(9))
      allocate(eos%iapws_m_vapor_coeff%ir(13),eos%iapws_m_vapor_coeff%jr(13),eos%iapws_m_vapor_coeff%nr(13))

      eos%iapws_liquid_coeff%i(1)  =  0;     eos%iapws_liquid_coeff%j(1)  =   -2;    eos%iapws_liquid_coeff%n(1)  =   0.14632971213167d0
      eos%iapws_liquid_coeff%i(2)  =  0;     eos%iapws_liquid_coeff%j(2)  =   -1;    eos%iapws_liquid_coeff%n(2)  = - 0.84548187169114d0
      eos%iapws_liquid_coeff%i(3)  =  0;     eos%iapws_liquid_coeff%j(3)  =    0;    eos%iapws_liquid_coeff%n(3)  = - 0.37563603672040d1
      eos%iapws_liquid_coeff%i(4)  =  0;     eos%iapws_liquid_coeff%j(4)  =    1;    eos%iapws_liquid_coeff%n(4)  =   0.33855169168385d1
      eos%iapws_liquid_coeff%i(5)  =  0;     eos%iapws_liquid_coeff%j(5)  =    2;    eos%iapws_liquid_coeff%n(5)  = - 0.95791963387872d0
      eos%iapws_liquid_coeff%i(6)  =  0;     eos%iapws_liquid_coeff%j(6)  =    3;    eos%iapws_liquid_coeff%n(6)  =   0.15772038513228d0
      eos%iapws_liquid_coeff%i(7)  =  0;     eos%iapws_liquid_coeff%j(7)  =    4;    eos%iapws_liquid_coeff%n(7)  = - 0.16616417199501d-1
      eos%iapws_liquid_coeff%i(8)  =  0;     eos%iapws_liquid_coeff%j(8)  =    5;    eos%iapws_liquid_coeff%n(8)  =   0.81214629983568d-3
      eos%iapws_liquid_coeff%i(9)  =  1;     eos%iapws_liquid_coeff%j(9)  =   -9;    eos%iapws_liquid_coeff%n(9)  =   0.28319080123804d-3
      eos%iapws_liquid_coeff%i(10) =  1;     eos%iapws_liquid_coeff%j(10)  =  -7;    eos%iapws_liquid_coeff%n(10) = - 0.60706301565874d-3
      eos%iapws_liquid_coeff%i(11) =  1;     eos%iapws_liquid_coeff%j(11)  =  -1;    eos%iapws_liquid_coeff%n(11) = - 0.18990068218419d-1
      eos%iapws_liquid_coeff%i(12) =  1;     eos%iapws_liquid_coeff%j(12)  =   0;    eos%iapws_liquid_coeff%n(12) = - 0.32529748770505d-1
      eos%iapws_liquid_coeff%i(13) =  1;     eos%iapws_liquid_coeff%j(13)  =   1;    eos%iapws_liquid_coeff%n(13) = - 0.21841717175414d-1
      eos%iapws_liquid_coeff%i(14) =  1;     eos%iapws_liquid_coeff%j(14)  =   3;    eos%iapws_liquid_coeff%n(14) = - 0.52838357969930d-4
      eos%iapws_liquid_coeff%i(15) =  2;     eos%iapws_liquid_coeff%j(15)  =  -3;    eos%iapws_liquid_coeff%n(15) = - 0.47184321073267d-3
      eos%iapws_liquid_coeff%i(16) =  2;     eos%iapws_liquid_coeff%j(16)  =   0;    eos%iapws_liquid_coeff%n(16) = - 0.30001780793026d-3
      eos%iapws_liquid_coeff%i(17) =  2;     eos%iapws_liquid_coeff%j(17)  =   1;    eos%iapws_liquid_coeff%n(17) =   0.47661393906987d-4
      eos%iapws_liquid_coeff%i(18) =  2;     eos%iapws_liquid_coeff%j(18)  =   3;    eos%iapws_liquid_coeff%n(18) = - 0.44141845330846d-5
      eos%iapws_liquid_coeff%i(19) =  2;     eos%iapws_liquid_coeff%j(19)  =  17;    eos%iapws_liquid_coeff%n(19) = - 0.72694996297594d-15
      eos%iapws_liquid_coeff%i(20) =  3;     eos%iapws_liquid_coeff%j(20)  =  -4;    eos%iapws_liquid_coeff%n(20) = - 0.31679644845054d-4
      eos%iapws_liquid_coeff%i(21) =  3;     eos%iapws_liquid_coeff%j(21)  =   0;    eos%iapws_liquid_coeff%n(21) = - 0.28270797985312d-5
      eos%iapws_liquid_coeff%i(22) =  3;     eos%iapws_liquid_coeff%j(22)  =   6;    eos%iapws_liquid_coeff%n(22) = - 0.85205128120103d-9
      eos%iapws_liquid_coeff%i(23) =  4;     eos%iapws_liquid_coeff%j(23)  =  -5;    eos%iapws_liquid_coeff%n(23) = - 0.22425281908000d-5
      eos%iapws_liquid_coeff%i(24) =  4;     eos%iapws_liquid_coeff%j(24)  =  -2;    eos%iapws_liquid_coeff%n(24) = - 0.65171222895601d-6
      eos%iapws_liquid_coeff%i(25) =  4;     eos%iapws_liquid_coeff%j(25)  =  10;    eos%iapws_liquid_coeff%n(25) = - 0.14341729937924d-12
      eos%iapws_liquid_coeff%i(26) =  5;     eos%iapws_liquid_coeff%j(26)  =  -8;    eos%iapws_liquid_coeff%n(26) = - 0.40516996860117d-6
      eos%iapws_liquid_coeff%i(27) =  8;     eos%iapws_liquid_coeff%j(27)  = -11;    eos%iapws_liquid_coeff%n(27) = - 0.12734301741641d-8
      eos%iapws_liquid_coeff%i(28) =  8;     eos%iapws_liquid_coeff%j(28)  =  -6;    eos%iapws_liquid_coeff%n(28) = - 0.17424871230634d-9
      eos%iapws_liquid_coeff%i(29) = 21;     eos%iapws_liquid_coeff%j(29)  = -29;    eos%iapws_liquid_coeff%n(29) = - 0.68762131295531d-18
      eos%iapws_liquid_coeff%i(30) = 23;     eos%iapws_liquid_coeff%j(30)  = -31;    eos%iapws_liquid_coeff%n(30) =   0.14478307828521d-19
      eos%iapws_liquid_coeff%i(31) = 29;     eos%iapws_liquid_coeff%j(31)  = -38;    eos%iapws_liquid_coeff%n(31) =   0.26335781662795d-22
      eos%iapws_liquid_coeff%i(32) = 30;     eos%iapws_liquid_coeff%j(32)  = -39;    eos%iapws_liquid_coeff%n(32) = - 0.11947622640071d-22
      eos%iapws_liquid_coeff%i(33) = 31;     eos%iapws_liquid_coeff%j(33)  = -40;    eos%iapws_liquid_coeff%n(33) =   0.18228094581404d-23
      eos%iapws_liquid_coeff%i(34) = 32;     eos%iapws_liquid_coeff%j(34)  = -41;    eos%iapws_liquid_coeff%n(34) = - 0.93537087292458d-25

      eos%iapws_vapor_coeff%j(1)  =    0;    eos%iapws_vapor_coeff%n(1)  = - 0.96927686500217d1
      eos%iapws_vapor_coeff%j(2)  =    1;    eos%iapws_vapor_coeff%n(2)  =   0.10086655968018d2
      eos%iapws_vapor_coeff%j(3)  =   -5;    eos%iapws_vapor_coeff%n(3)  = - 0.56087911283020d-2
      eos%iapws_vapor_coeff%j(4)  =   -4;    eos%iapws_vapor_coeff%n(4)  =   0.71452738081455d-1
      eos%iapws_vapor_coeff%j(5)  =   -3;    eos%iapws_vapor_coeff%n(5)  = - 0.40710498223928d0
      eos%iapws_vapor_coeff%j(6)  =   -2;    eos%iapws_vapor_coeff%n(6)  =   0.14240819171444d1
      eos%iapws_vapor_coeff%j(7)  =   -1;    eos%iapws_vapor_coeff%n(7)  = - 0.43839511319450d1
      eos%iapws_vapor_coeff%j(8)  =    2;    eos%iapws_vapor_coeff%n(8)  = - 0.28408632460772d0
      eos%iapws_vapor_coeff%j(9)  =    3;    eos%iapws_vapor_coeff%n(9)  =   0.21268463753307d-1

      eos%iapws_m_vapor_coeff%j(1)  =    0;    eos%iapws_m_vapor_coeff%n(1)  = - 0.96937268393049d1
      eos%iapws_m_vapor_coeff%j(2)  =    1;    eos%iapws_m_vapor_coeff%n(2)  =   0.10087275970006d2
      eos%iapws_m_vapor_coeff%j(3)  =   -5;    eos%iapws_m_vapor_coeff%n(3)  = - 0.56087911283020d-2
      eos%iapws_m_vapor_coeff%j(4)  =   -4;    eos%iapws_m_vapor_coeff%n(4)  =   0.71452738081455d-1
      eos%iapws_m_vapor_coeff%j(5)  =   -3;    eos%iapws_m_vapor_coeff%n(5)  = - 0.40710498223928d0
      eos%iapws_m_vapor_coeff%j(6)  =   -2;    eos%iapws_m_vapor_coeff%n(6)  =   0.14240819171444d1
      eos%iapws_m_vapor_coeff%j(7)  =   -1;    eos%iapws_m_vapor_coeff%n(7)  = - 0.43839511319450d1
      eos%iapws_m_vapor_coeff%j(8)  =    2;    eos%iapws_m_vapor_coeff%n(8)  = - 0.28408632460772d0
      eos%iapws_m_vapor_coeff%j(9)  =    3;    eos%iapws_m_vapor_coeff%n(9)  =   0.21268463753307d-1

      eos%iapws_vapor_coeff%ir(1)  =  1;     eos%iapws_vapor_coeff%jr(1)  =    0;    eos%iapws_vapor_coeff%nr(1)  = - 0.17731742473213d-2
      eos%iapws_vapor_coeff%ir(2)  =  1;     eos%iapws_vapor_coeff%jr(2)  =    1;    eos%iapws_vapor_coeff%nr(2)  = - 0.17834862292358d-1
      eos%iapws_vapor_coeff%ir(3)  =  1;     eos%iapws_vapor_coeff%jr(3)  =    2;    eos%iapws_vapor_coeff%nr(3)  = - 0.45996013696365d-1
      eos%iapws_vapor_coeff%ir(4)  =  1;     eos%iapws_vapor_coeff%jr(4)  =    3;    eos%iapws_vapor_coeff%nr(4)  = - 0.57581259083432d-1
      eos%iapws_vapor_coeff%ir(5)  =  1;     eos%iapws_vapor_coeff%jr(5)  =    6;    eos%iapws_vapor_coeff%nr(5)  = - 0.50325278727930d-1
      eos%iapws_vapor_coeff%ir(6)  =  2;     eos%iapws_vapor_coeff%jr(6)  =    1;    eos%iapws_vapor_coeff%nr(6)  = - 0.33032641670203d-4
      eos%iapws_vapor_coeff%ir(7)  =  2;     eos%iapws_vapor_coeff%jr(7)  =    2;    eos%iapws_vapor_coeff%nr(7)  = - 0.18948987516315d-3
      eos%iapws_vapor_coeff%ir(8)  =  2;     eos%iapws_vapor_coeff%jr(8)  =    4;    eos%iapws_vapor_coeff%nr(8)  = - 0.39392777243355d-2
      eos%iapws_vapor_coeff%ir(9)  =  2;     eos%iapws_vapor_coeff%jr(9)  =    7;    eos%iapws_vapor_coeff%nr(9)  = - 0.43797295650573d-1
      eos%iapws_vapor_coeff%ir(10) =  2;     eos%iapws_vapor_coeff%jr(10)  =  36;    eos%iapws_vapor_coeff%nr(10) = - 0.26674547914087d-4
      eos%iapws_vapor_coeff%ir(11) =  3;     eos%iapws_vapor_coeff%jr(11)  =   0;    eos%iapws_vapor_coeff%nr(11) =   0.20481737692309d-7
      eos%iapws_vapor_coeff%ir(12) =  3;     eos%iapws_vapor_coeff%jr(12)  =   1;    eos%iapws_vapor_coeff%nr(12) =   0.43870667284435d-6
      eos%iapws_vapor_coeff%ir(13) =  3;     eos%iapws_vapor_coeff%jr(13)  =   3;    eos%iapws_vapor_coeff%nr(13) = - 0.32277677238570d-4
      eos%iapws_vapor_coeff%ir(14) =  3;     eos%iapws_vapor_coeff%jr(14)  =   6;    eos%iapws_vapor_coeff%nr(14) = - 0.15033924542148d-2
      eos%iapws_vapor_coeff%ir(15) =  3;     eos%iapws_vapor_coeff%jr(15)  =  35;    eos%iapws_vapor_coeff%nr(15) = - 0.40668253562649d-1
      eos%iapws_vapor_coeff%ir(16) =  4;     eos%iapws_vapor_coeff%jr(16)  =   1;    eos%iapws_vapor_coeff%nr(16) = - 0.78847309559367d-9
      eos%iapws_vapor_coeff%ir(17) =  4;     eos%iapws_vapor_coeff%jr(17)  =   2;    eos%iapws_vapor_coeff%nr(17) =   0.12790717852285d-7
      eos%iapws_vapor_coeff%ir(18) =  4;     eos%iapws_vapor_coeff%jr(18)  =   3;    eos%iapws_vapor_coeff%nr(18) =   0.48225372718507d-6
      eos%iapws_vapor_coeff%ir(19) =  5;     eos%iapws_vapor_coeff%jr(19)  =   7;    eos%iapws_vapor_coeff%nr(19) =   0.22922076337661d-5
      eos%iapws_vapor_coeff%ir(20) =  6;     eos%iapws_vapor_coeff%jr(20)  =   3;    eos%iapws_vapor_coeff%nr(20) = - 0.16714766451061d-10
      eos%iapws_vapor_coeff%ir(21) =  6;     eos%iapws_vapor_coeff%jr(21)  =  16;    eos%iapws_vapor_coeff%nr(21) = - 0.21171472321355d-2
      eos%iapws_vapor_coeff%ir(22) =  6;     eos%iapws_vapor_coeff%jr(22)  =  35;    eos%iapws_vapor_coeff%nr(22) = - 0.23895741934104d2
      eos%iapws_vapor_coeff%ir(23) =  7;     eos%iapws_vapor_coeff%jr(23)  =   0;    eos%iapws_vapor_coeff%nr(23) = - 0.59059564324270d-17
      eos%iapws_vapor_coeff%ir(24) =  7;     eos%iapws_vapor_coeff%jr(24)  =  11;    eos%iapws_vapor_coeff%nr(24) = - 0.12621808899101d-5
      eos%iapws_vapor_coeff%ir(25) =  7;     eos%iapws_vapor_coeff%jr(25)  =  25;    eos%iapws_vapor_coeff%nr(25) = - 0.38946842435739d-1
      eos%iapws_vapor_coeff%ir(26) =  8;     eos%iapws_vapor_coeff%jr(26)  =   8;    eos%iapws_vapor_coeff%nr(26) =   0.11256211360459d-10
      eos%iapws_vapor_coeff%ir(27) =  8;     eos%iapws_vapor_coeff%jr(27)  =  36;    eos%iapws_vapor_coeff%nr(27) = - 0.82311340897998d1
      eos%iapws_vapor_coeff%ir(28) =  9;     eos%iapws_vapor_coeff%jr(28)  =  13;    eos%iapws_vapor_coeff%nr(28) =   0.19809712802088d-7
      eos%iapws_vapor_coeff%ir(29) = 10;     eos%iapws_vapor_coeff%jr(29)  =   4;    eos%iapws_vapor_coeff%nr(29) =   0.10406965210174d-18
      eos%iapws_vapor_coeff%ir(30) = 10;     eos%iapws_vapor_coeff%jr(30)  =  10;    eos%iapws_vapor_coeff%nr(30) = - 0.10234747095929d-12
      eos%iapws_vapor_coeff%ir(31) = 10;     eos%iapws_vapor_coeff%jr(31)  =  14;    eos%iapws_vapor_coeff%nr(31) = - 0.10018179379511d-8
      eos%iapws_vapor_coeff%ir(32) = 16;     eos%iapws_vapor_coeff%jr(32)  =  29;    eos%iapws_vapor_coeff%nr(32) = - 0.80882908646985d-10
      eos%iapws_vapor_coeff%ir(33) = 16;     eos%iapws_vapor_coeff%jr(33)  =  50;    eos%iapws_vapor_coeff%nr(33) =   0.10693031879409d0
      eos%iapws_vapor_coeff%ir(34) = 18;     eos%iapws_vapor_coeff%jr(34)  =  57;    eos%iapws_vapor_coeff%nr(34) = - 0.33662250574171d0
      eos%iapws_vapor_coeff%ir(35) = 20;     eos%iapws_vapor_coeff%jr(35)  =  20;    eos%iapws_vapor_coeff%nr(35) =   0.89185845355421d-24
      eos%iapws_vapor_coeff%ir(36) = 20;     eos%iapws_vapor_coeff%jr(36)  =  35;    eos%iapws_vapor_coeff%nr(36) =   0.30629316876232d-12
      eos%iapws_vapor_coeff%ir(37) = 20;     eos%iapws_vapor_coeff%jr(37)  =  48;    eos%iapws_vapor_coeff%nr(37) = - 0.42002467698208d-5
      eos%iapws_vapor_coeff%ir(38) = 21;     eos%iapws_vapor_coeff%jr(38)  =  21;    eos%iapws_vapor_coeff%nr(38) = - 0.59056029685639d-25
      eos%iapws_vapor_coeff%ir(39) = 22;     eos%iapws_vapor_coeff%jr(39)  =  53;    eos%iapws_vapor_coeff%nr(39) =   0.37826947613457d-5
      eos%iapws_vapor_coeff%ir(40) = 23;     eos%iapws_vapor_coeff%jr(40)  =  39;    eos%iapws_vapor_coeff%nr(40) = - 0.12768608934681d-14
      eos%iapws_vapor_coeff%ir(41) = 24;     eos%iapws_vapor_coeff%jr(41)  =  26;    eos%iapws_vapor_coeff%nr(41) =   0.73087610595061d-28
      eos%iapws_vapor_coeff%ir(42) = 24;     eos%iapws_vapor_coeff%jr(42)  =  40;    eos%iapws_vapor_coeff%nr(42) =   0.55414715350778d-16
      eos%iapws_vapor_coeff%ir(43) = 24;     eos%iapws_vapor_coeff%jr(43)  =  58;    eos%iapws_vapor_coeff%nr(43) = - 0.94369707241210d-6

      eos%iapws_m_vapor_coeff%ir(1)  =  1;     eos%iapws_m_vapor_coeff%jr(1)  =    0;    eos%iapws_m_vapor_coeff%nr(1)  = - 0.73362260186506d-2
      eos%iapws_m_vapor_coeff%ir(2)  =  1;     eos%iapws_m_vapor_coeff%jr(2)  =    2;    eos%iapws_m_vapor_coeff%nr(2)  = - 0.88223831943146d-1
      eos%iapws_m_vapor_coeff%ir(3)  =  1;     eos%iapws_m_vapor_coeff%jr(3)  =    5;    eos%iapws_m_vapor_coeff%nr(3)  = - 0.72334555213245d-1
      eos%iapws_m_vapor_coeff%ir(4)  =  1;     eos%iapws_m_vapor_coeff%jr(4)  =   11;    eos%iapws_m_vapor_coeff%nr(4)  = - 0.40813178534455d-2
      eos%iapws_m_vapor_coeff%ir(5)  =  2;     eos%iapws_m_vapor_coeff%jr(5)  =    1;    eos%iapws_m_vapor_coeff%nr(5)  =   0.20097803380207d-2
      eos%iapws_m_vapor_coeff%ir(6)  =  2;     eos%iapws_m_vapor_coeff%jr(6)  =    7;    eos%iapws_m_vapor_coeff%nr(6)  = - 0.53045921898642d-1
      eos%iapws_m_vapor_coeff%ir(7)  =  2;     eos%iapws_m_vapor_coeff%jr(7)  =   16;    eos%iapws_m_vapor_coeff%nr(7)  = - 0.76190409086970d-2
      eos%iapws_m_vapor_coeff%ir(8)  =  3;     eos%iapws_m_vapor_coeff%jr(8)  =    4;    eos%iapws_m_vapor_coeff%nr(8)  = - 0.63498037657313d-2
      eos%iapws_m_vapor_coeff%ir(9)  =  3;     eos%iapws_m_vapor_coeff%jr(9)  =   16;    eos%iapws_m_vapor_coeff%nr(9)  = - 0.86043093028588d-1
      eos%iapws_m_vapor_coeff%ir(10) =  4;     eos%iapws_m_vapor_coeff%jr(10)  =   7;    eos%iapws_m_vapor_coeff%nr(10) =   0.75321581522770d-2
      eos%iapws_m_vapor_coeff%ir(11) =  4;     eos%iapws_m_vapor_coeff%jr(11)  =  10;    eos%iapws_m_vapor_coeff%nr(11) = - 0.79238375446139d-2
      eos%iapws_m_vapor_coeff%ir(12) =  5;     eos%iapws_m_vapor_coeff%jr(12)  =   9;    eos%iapws_m_vapor_coeff%nr(12) = - 0.22888160778447d-3
      eos%iapws_m_vapor_coeff%ir(13) =  5;     eos%iapws_m_vapor_coeff%jr(13)  =  10;    eos%iapws_m_vapor_coeff%nr(13) = - 0.26456501482810d-2

    end subroutine set_iapws97_eos
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_iapws97_prop(eos)
      implicit none
      class(t_eos), intent(inout) :: eos

      allocate(eos%iapws_prop_coeff%h(0:3),eos%iapws_prop_coeff%l(0:3))
      allocate(eos%iapws_prop_coeff%hh(0:5,0:6),eos%iapws_prop_coeff%ll(0:4,0:5))

      eos%iapws_prop_coeff%h(0) =   1.d0          ;    eos%iapws_prop_coeff%l(0) =   1.d0
      eos%iapws_prop_coeff%h(1) =   0.978197d0    ;    eos%iapws_prop_coeff%l(1) =   6.978267d0
      eos%iapws_prop_coeff%h(2) =   0.579829d0    ;    eos%iapws_prop_coeff%l(2) =   2.599096d0
      eos%iapws_prop_coeff%h(3) = - 0.202354d0    ;    eos%iapws_prop_coeff%l(3) = - 0.998254d0


      eos%iapws_prop_coeff%hh(0,0) =   0.5132047d0;    eos%iapws_prop_coeff%hh(1,0) =   0.3205656d0
      eos%iapws_prop_coeff%hh(0,1) =   0.2151778d0;    eos%iapws_prop_coeff%hh(1,1) =   0.7317883d0
      eos%iapws_prop_coeff%hh(0,2) = - 0.2818107d0;    eos%iapws_prop_coeff%hh(1,2) = - 1.070786d0
      eos%iapws_prop_coeff%hh(0,3) =   0.1778064d0;    eos%iapws_prop_coeff%hh(1,3) =   0.460504d0
      eos%iapws_prop_coeff%hh(0,4) = - 0.0417661d0;    eos%iapws_prop_coeff%hh(1,4) =   0.d0
      eos%iapws_prop_coeff%hh(0,5) =   0.d0       ;    eos%iapws_prop_coeff%hh(1,5) = - 0.01578386d0
      eos%iapws_prop_coeff%hh(0,6) =   0.d0       ;    eos%iapws_prop_coeff%hh(1,6) =   0.d0

      eos%iapws_prop_coeff%hh(2,0) =   0.d0       ;    eos%iapws_prop_coeff%hh(3,0) =   0.d0
      eos%iapws_prop_coeff%hh(2,1) =   1.241044d0 ;    eos%iapws_prop_coeff%hh(3,1) =   1.476783d0
      eos%iapws_prop_coeff%hh(2,2) = - 1.263184d0 ;    eos%iapws_prop_coeff%hh(3,2) =   0.d0
      eos%iapws_prop_coeff%hh(2,3) =   0.2340379d0;    eos%iapws_prop_coeff%hh(3,3) = - 0.4924179d0
      eos%iapws_prop_coeff%hh(2,4) =   0.d0       ;    eos%iapws_prop_coeff%hh(3,4) =   0.1600435d0
      eos%iapws_prop_coeff%hh(2,5) =   0.d0       ;    eos%iapws_prop_coeff%hh(3,5) =   0.d0
      eos%iapws_prop_coeff%hh(2,6) =   0.d0       ;    eos%iapws_prop_coeff%hh(3,6) = - 0.003629481d0

      eos%iapws_prop_coeff%hh(4,0) = - 0.7782567d0;    eos%iapws_prop_coeff%hh(5,0) =   0.1885447d0
      eos%iapws_prop_coeff%hh(4,1) =   0.d0       ;    eos%iapws_prop_coeff%hh(5,1) =   0.d0
      eos%iapws_prop_coeff%hh(4,2) =   0.d0       ;    eos%iapws_prop_coeff%hh(5,2) =   0.d0
      eos%iapws_prop_coeff%hh(4,3) =   0.d0       ;    eos%iapws_prop_coeff%hh(5,3) =   0.d0
      eos%iapws_prop_coeff%hh(4,4) =   0.d0       ;    eos%iapws_prop_coeff%hh(5,4) =   0.d0
      eos%iapws_prop_coeff%hh(4,5) =   0.d0       ;    eos%iapws_prop_coeff%hh(5,5) =   0.d0
      eos%iapws_prop_coeff%hh(4,6) =   0.d0       ;    eos%iapws_prop_coeff%hh(5,6) =   0.d0

      eos%iapws_prop_coeff%ll(0,0) =   1.3293046d0  ;    eos%iapws_prop_coeff%ll(1,0) =   1.7018363d0
      eos%iapws_prop_coeff%ll(0,1) = - 0.40452437d0 ;    eos%iapws_prop_coeff%ll(1,1) = - 2.2156845d0
      eos%iapws_prop_coeff%ll(0,2) =   0.24409490d0 ;    eos%iapws_prop_coeff%ll(1,2) =   1.6511057d0
      eos%iapws_prop_coeff%ll(0,3) =   0.018660751d0;    eos%iapws_prop_coeff%ll(1,3) = - 0.76736002d0
      eos%iapws_prop_coeff%ll(0,4) = - 0.12961068d0 ;    eos%iapws_prop_coeff%ll(1,4) =   0.37283344d0
      eos%iapws_prop_coeff%ll(0,5) =   0.044809953d0;    eos%iapws_prop_coeff%ll(1,5) = - 0.11203160d0

      eos%iapws_prop_coeff%ll(2,0) =   5.2246158d0  ;    eos%iapws_prop_coeff%ll(3,0) =   8.7127675d0
      eos%iapws_prop_coeff%ll(2,1) = - 10.124111d0  ;    eos%iapws_prop_coeff%ll(3,1) = - 9.5000611d0
      eos%iapws_prop_coeff%ll(2,2) =   4.9874687d0  ;    eos%iapws_prop_coeff%ll(3,2) =   4.3786606d0
      eos%iapws_prop_coeff%ll(2,3) = - 0.27297694d0 ;    eos%iapws_prop_coeff%ll(3,3) = - 0.91783782d0
      eos%iapws_prop_coeff%ll(2,4) = - 0.43083393d0 ;    eos%iapws_prop_coeff%ll(3,4) =   0.d0
      eos%iapws_prop_coeff%ll(2,5) =   0.13333849d0 ;    eos%iapws_prop_coeff%ll(3,5) =   0.d0

      eos%iapws_prop_coeff%ll(4,0) = - 1.8525999d0
      eos%iapws_prop_coeff%ll(4,1) =   0.93404690d0
      eos%iapws_prop_coeff%ll(4,2) =   0.d0
      eos%iapws_prop_coeff%ll(4,3) =   0.d0
      eos%iapws_prop_coeff%ll(4,4) =   0.d0
      eos%iapws_prop_coeff%ll(4,5) =   0.d0
    end subroutine set_iapws97_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_srk_property(eos,fluid,srk_property)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer, intent(in) :: fluid
      type(t_srk_property),intent(out) :: srk_property
      real(8),parameter :: ru=8.314d0

      select case(fluid)
      case(2) ! nitrogen
        srk_property%tc  = 126.19d0
        srk_property%pc  = 3.3958d6
        srk_property%tb  = 77.355d0
        srk_property%mw  = 28.013d-3
        srk_property%w   = 0.0372d0
        srk_property%cv0 = 742.2d0

      case(3) ! oxygen
        srk_property%tc  = 154.58d0
        srk_property%pc  = 5.043d6
        srk_property%tb  = 90.188d0
        srk_property%mw  = 31.999d-3
        srk_property%w   = 0.0222d0
        srk_property%cv0 = 650.d0

      case(4) ! hydrogen
        srk_property%tc  = 33.145d0
        srk_property%pc  = 1.2964d6
        srk_property%tb  = 20.369d0      ! boiling temp
        srk_property%mw  = 2.0159d-3  ! Molar mass
        srk_property%w   = -0.219d0     ! Acentric factor
        srk_property%cv0 = 6186.6d0     ! should be filled   ! 

      case(5) ! helium
        srk_property%tc  = 5.1953d0
        srk_property%pc  = 0.22746d6
        srk_property%tb  = 4.23d0
        srk_property%mw  = 4.0026d-3
        srk_property%w   = -0.382d0
        srk_property%cv0 = 0.d0         ! should be filled

      case(7) ! dimethylether
        srk_property%tc  = 400.38d0
        srk_property%pc  = 5.3368d6
        srk_property%tb  = 4.23d0
        srk_property%mw  = 4.0026d-3
        srk_property%w   = -0.382d0
        srk_property%cv0 = 0.d0         ! should be filled

      case default
      end select

      srk_property%s     = 0.48508d0+1.5517d0*srk_property%w-0.15613d0*srk_property%w**2
      srk_property%a     = 0.42747d0*(ru*srk_property%tc)**2/srk_property%pc
      srk_property%b     = 0.08664d0*ru*srk_property%tc/srk_property%pc

    end subroutine set_srk_property
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_database(eos,fluid,phase)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer, intent(in) :: fluid,phase
      integer :: n,l,k
      integer :: m,io,ier
      integer :: intsize,realsize
      integer(kind=mpi_offset_kind) :: disp
      character(3) :: c_phase,c_fluid

      select case(fluid)
      case(1) ! water
        c_fluid = 'h2o'
      case(2) ! nitrogen
        c_fluid = 'n2'
      case(3) ! oxygen
        c_fluid = 'o2'
      case(4) ! hydrogen
        c_fluid = 'h2'
      case(5) ! helium
        c_fluid = 'he'
      case(7) ! dimethylether
        c_fluid = 'DME'
      case default
      end select

      select case(phase)
      case(1) ! liquid
        c_phase = 'liq'
      case(2,3) ! vapor, gas
        c_phase = 'gas'
      case default
      end select

!! database set start
 !      do k=0,eos%size-1
 !        if(k.eq.eos%rank) then
 !          open(unit=10,file='./../fld/'//trim(c_fluid)//'_'//trim(c_phase)//'_old.fld',form='binary')
 !
 !          read(10) eos%db_const(phase)%t_ndata
 !          allocate(eos%db_const(phase)%db(eos%db_const(phase)%t_ndata))
 !
 !          read(10) eos%db_const(phase)%t_nsub
 !          allocate(eos%db_const(phase)%t_boundary(eos%db_const(phase)%t_nsub+1))
 !          allocate(eos%db_const(phase)%t_iboundary(eos%db_const(phase)%t_nsub+1))
 !          allocate(eos%db_const(phase)%t_delta(eos%db_const(phase)%t_nsub))
 !
 !          do n=1,eos%db_const(phase)%t_nsub+1
 !            read(10) eos%db_const(phase)%t_boundary(n)
 !            read(10) eos%db_const(phase)%t_iboundary(n)
 !          end do
 !
 !          do n=1,eos%db_const(phase)%t_nsub
 !            read(10) eos%db_const(phase)%t_delta(n)
 !          end do
 !
 !          do m=1,eos%db_const(phase)%t_ndata
 !            read(10) eos%db_const(phase)%db(m)%t
 !            read(10) eos%db_const(phase)%db(m)%p_nstart,eos%db_const(phase)%db(m)%p_nend
 !            allocate(eos%db_const(phase)%db(m)%p(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend))
 !
 !            read(10) eos%db_const(phase)%db(m)%p_nsub
 !            allocate(eos%db_const(phase)%db(m)%p_boundary(eos%db_const(phase)%db(m)%p_nsub+1))
 !            allocate(eos%db_const(phase)%db(m)%p_iboundary(eos%db_const(phase)%db(m)%p_nsub+1))
 !            allocate(eos%db_const(phase)%db(m)%p_delta(eos%db_const(phase)%db(m)%p_nsub))
 !
 !            do n=1,eos%db_const(phase)%db(m)%p_nsub+1
 !              read(10) eos%db_const(phase)%db(m)%p_boundary(n)
 !              read(10) eos%db_const(phase)%db(m)%p_iboundary(n)
 !            end do
 !
 !            do n=1,eos%db_const(phase)%db(m)%p_nsub
 !              read(10) eos%db_const(phase)%db(m)%p_delta(n)
 !            end do
 !
 !            do n=eos%db_const(phase)%db(m)%p_nstart,eos%db_const(phase)%db(m)%p_nend
 !              read(10) eos%db_const(phase)%db(m)%p(n)
 !            end do
 !
 !          end do
 !
 !          do m=1,eos%db_const(phase)%t_ndata-1
 !            allocate(eos%db_const(phase)%db(m)%rho(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,9))
 !            allocate(eos%db_const(phase)%db(m)%h(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,9))
 !            allocate(eos%db_const(phase)%db(m)%vis(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,4))
 !            allocate(eos%db_const(phase)%db(m)%cond(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,4))
 !
 !            do n=eos%db_const(phase)%db(m)%p_nstart,eos%db_const(phase)%db(m)%p_nend-1
 !              do l=1,9
 !                read(10) eos%db_const(phase)%db(m)%rho(n,l)
 !              end do
 !              do l=1,9
 !                read(10) eos%db_const(phase)%db(m)%h(n,l)
 !              end do
 !              do l=1,4
 !                read(10) eos%db_const(phase)%db(m)%vis(n,l)
 !              end do
 !              do l=1,4
 !                read(10) eos%db_const(phase)%db(m)%cond(n,l)
 !              end do
 !            end do
 !          end do
 !          close(10)
 !
 !        end if
 !      end do
 !
 !      call mpi_type_size(mpi_integer,intsize,ier)
 !      call mpi_type_size(mpi_real8,realsize,ier)
 !
 !      call mpi_file_open(mpi_comm_world,'./../fld/'//trim(c_fluid)//'_'//trim(c_phase)//'.fld',mpi_mode_wronly+mpi_mode_create,mpi_info_null,io,ier)
 !
 !      disp = 0
 !      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
 !      call mpi_file_write_all(io,eos%db_const(phase)%t_ndata,1,mpi_integer,mpi_status_ignore,ier)
 !      disp = disp + intsize
 !
 !      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
 !      call mpi_file_write_all(io,eos%db_const(phase)%t_nsub,1,mpi_integer,mpi_status_ignore,ier)
 !      disp = disp + intsize
 !
 !      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
 !      call mpi_file_write_all(io,eos%db_const(phase)%t_iboundary,eos%db_const(phase)%t_nsub+1,mpi_integer,mpi_status_ignore,ier)
 !      disp = disp + intsize*(eos%db_const(phase)%t_nsub+1)
 !
 !      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !      call mpi_file_write_all(io,eos%db_const(phase)%t_boundary,eos%db_const(phase)%t_nsub+1,mpi_real8,mpi_status_ignore,ier)
 !      disp = disp + realsize*(eos%db_const(phase)%t_nsub+1)
 !
 !      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !      call mpi_file_write_all(io,eos%db_const(phase)%t_delta,eos%db_const(phase)%t_nsub,mpi_real8,mpi_status_ignore,ier)
 !      disp = disp + realsize*eos%db_const(phase)%t_nsub
 !
 !      do m=1,eos%db_const(phase)%t_ndata
 !
 !        call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%t,1,mpi_real8,mpi_status_ignore,ier)
 !        disp = disp + realsize
 !
 !        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%p_nstart,1,mpi_integer,mpi_status_ignore,ier)
 !        disp = disp + intsize
 !
 !        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%p_nend,1,mpi_integer,mpi_status_ignore,ier)
 !        disp = disp + intsize
 !
 !        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%p_nsub,1,mpi_integer,mpi_status_ignore,ier)
 !        disp = disp + intsize
 !
 !        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%p_iboundary,eos%db_const(phase)%db(m)%p_nsub+1,mpi_integer,mpi_status_ignore,ier)
 !        disp = disp + intsize*(eos%db_const(phase)%db(m)%p_nsub+1)
 !
 !        call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%p_boundary,eos%db_const(phase)%db(m)%p_nsub+1,mpi_real8,mpi_status_ignore,ier)
 !        disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nsub+1)
 !
 !        call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%p_delta,eos%db_const(phase)%db(m)%p_nsub,mpi_real8,mpi_status_ignore,ier)
 !        disp = disp + realsize*eos%db_const(phase)%db(m)%p_nsub
 !
 !        call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !        call mpi_file_write_all(io,eos%db_const(phase)%db(m)%p,eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart+1,mpi_real8,mpi_status_ignore,ier)
 !        disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart+1)
 !
 !        if (m.lt.eos%db_const(phase)%t_ndata) then
 !          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !          call mpi_file_write_all(io,eos%db_const(phase)%db(m)%rho,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9,mpi_real8,mpi_status_ignore,ier)
 !          disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9
 !
 !          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !          call mpi_file_write_all(io,eos%db_const(phase)%db(m)%h,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9,mpi_real8,mpi_status_ignore,ier)
 !          disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9
 !
 !          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !          call mpi_file_write_all(io,eos%db_const(phase)%db(m)%vis,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4,mpi_real8,mpi_status_ignore,ier)
 !          disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4
 !
 !          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
 !          call mpi_file_write_all(io,eos%db_const(phase)%db(m)%cond,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4,mpi_real8,mpi_status_ignore,ier)
 !          disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4
 !        end if
 !
 !      end do
 !
 !      call mpi_file_close(io,ier)
!! database set end

        call mpi_type_size(mpi_integer,intsize,ier)
        call mpi_type_size(mpi_real8,realsize,ier)
  
        call mpi_file_open(mpi_comm_world,'./../fld/'//trim(c_fluid)//'_'//trim(c_phase)//'.fld',mpi_mode_rdonly,mpi_info_null,io,ier)
  
        disp = 0
        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,eos%db_const(phase)%t_ndata,1,mpi_integer,mpi_status_ignore,ier)
        disp = disp + intsize
        allocate(eos%db_const(phase)%db(eos%db_const(phase)%t_ndata))
  
        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,eos%db_const(phase)%t_nsub,1,mpi_integer,mpi_status_ignore,ier)
        disp = disp + intsize
        allocate(eos%db_const(phase)%t_boundary(eos%db_const(phase)%t_nsub+1))
        allocate(eos%db_const(phase)%t_iboundary(eos%db_const(phase)%t_nsub+1))
        allocate(eos%db_const(phase)%t_delta(eos%db_const(phase)%t_nsub))
  
        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,eos%db_const(phase)%t_iboundary,eos%db_const(phase)%t_nsub+1,mpi_integer,mpi_status_ignore,ier)
        disp = disp + intsize*(eos%db_const(phase)%t_nsub+1)
  
        call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,eos%db_const(phase)%t_boundary,eos%db_const(phase)%t_nsub+1,mpi_real8,mpi_status_ignore,ier)
        disp = disp + realsize*(eos%db_const(phase)%t_nsub+1)
  
        call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,eos%db_const(phase)%t_delta,eos%db_const(phase)%t_nsub,mpi_real8,mpi_status_ignore,ier)
        disp = disp + realsize*eos%db_const(phase)%t_nsub
  
        do m=1,eos%db_const(phase)%t_ndata
  
          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%t,1,mpi_real8,mpi_status_ignore,ier)
          disp = disp + realsize
  
          call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%p_nstart,1,mpi_integer,mpi_status_ignore,ier)
          disp = disp + intsize
  
          call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%p_nend,1,mpi_integer,mpi_status_ignore,ier)
          disp = disp + intsize
          allocate(eos%db_const(phase)%db(m)%p(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend))
  
          call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%p_nsub,1,mpi_integer,mpi_status_ignore,ier)
          disp = disp + intsize
          allocate(eos%db_const(phase)%db(m)%p_boundary(eos%db_const(phase)%db(m)%p_nsub+1))
          allocate(eos%db_const(phase)%db(m)%p_iboundary(eos%db_const(phase)%db(m)%p_nsub+1))
          allocate(eos%db_const(phase)%db(m)%p_delta(eos%db_const(phase)%db(m)%p_nsub))
  
          call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%p_iboundary,eos%db_const(phase)%db(m)%p_nsub+1,mpi_integer,mpi_status_ignore,ier)
          disp = disp + intsize*(eos%db_const(phase)%db(m)%p_nsub+1)
  
          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%p_boundary,eos%db_const(phase)%db(m)%p_nsub+1,mpi_real8,mpi_status_ignore,ier)
          disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nsub+1)
  
          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%p_delta,eos%db_const(phase)%db(m)%p_nsub,mpi_real8,mpi_status_ignore,ier)
          disp = disp + realsize*eos%db_const(phase)%db(m)%p_nsub
  
          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
          call mpi_file_read_all(io,eos%db_const(phase)%db(m)%p,eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart+1,mpi_real8,mpi_status_ignore,ier)
          disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart+1)
  
          if(m.lt.eos%db_const(phase)%t_ndata) then
            allocate(eos%db_const(phase)%db(m)%rho(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,9))
            allocate(eos%db_const(phase)%db(m)%h(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,9))
            allocate(eos%db_const(phase)%db(m)%vis(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,4))
            allocate(eos%db_const(phase)%db(m)%cond(eos%db_const(phase)%db(m)%p_nstart:eos%db_const(phase)%db(m)%p_nend-1,4))
  
            call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
            call mpi_file_read_all(io,eos%db_const(phase)%db(m)%rho,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9,mpi_real8,mpi_status_ignore,ier)
            disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9
  
            call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
            call mpi_file_read_all(io,eos%db_const(phase)%db(m)%h,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9,mpi_real8,mpi_status_ignore,ier)
            disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*9
  
            call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
            call mpi_file_read_all(io,eos%db_const(phase)%db(m)%vis,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4,mpi_real8,mpi_status_ignore,ier)
            disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4
  
            call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
            call mpi_file_read_all(io,eos%db_const(phase)%db(m)%cond,(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4,mpi_real8,mpi_status_ignore,ier)
            disp = disp + realsize*(eos%db_const(phase)%db(m)%p_nend-eos%db_const(phase)%db(m)%p_nstart)*4
          end if
  
        end do
  
        call mpi_file_close(io,ier)

    end subroutine set_database
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine deteos_simple(eos,p,t,y1,y2,dv)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(eos%ndv) ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      type(t_eos2)  :: l_eos,v_eos,g_eos
      real(8) :: irhol,irhov,irhog,irho,yl
      real(8) :: ic2l,ic2v,ic2g

      call eos%eos_l(p,t,1,l_eos)
      call eos%eos_v(p,t,2,v_eos)
      call eos%eos_g(p,t,3,g_eos)

      dv(3)  = l_eos%rho
      dv(4)  = v_eos%rho
      dv(5)  = g_eos%rho

      irhol = 1.d0/dv(3)
      irhov = 1.d0/dv(4)
      irhog = 1.d0/dv(5)
      yl = 1.d0-y1-y2
      irho = yl*irhol + y1*irhov + y2*irhog
      dv(1)  = 1.d0/irho
      dv(2)  = yl*l_eos%h + y1*v_eos%h + y2*g_eos%h

      dv(7)  = dv(1)**2*(yl*irhol**2*l_eos%drdp   &
                       + y1*irhov**2*v_eos%drdp   &
                       + y2*irhog**2*g_eos%drdp)
      dv(8)  = dv(1)**2*(yl*irhol**2*l_eos%drdt   &
                       + y1*irhov**2*v_eos%drdt   &
                       + y2*irhog**2*g_eos%drdt)
      dv(9)  = dv(1)**2*(irhol - irhov)
      dv(10) = dv(1)**2*(irhol - irhog)

      dv(11) = yl*l_eos%dhdp + y1*v_eos%dhdp + y2*g_eos%dhdp
      dv(12) = yl*l_eos%dhdt + y1*v_eos%dhdt + y2*g_eos%dhdt
      dv(13) = v_eos%h - l_eos%h
      dv(14) = g_eos%h - l_eos%h

      dv(6)  = dv(1)*dv(12)/(dv(1)*dv(7)*dv(12)+dv(8)*(1.d0-dv(1)*dv(11)))
      if(dv(6).le.0.d0) then
        write(*,*) 'alternative calculation for snd2'

        ic2l = (dv(3)*l_eos%drdp*l_eos%dhdt + l_eos%drdt*(1.d0-dv(3)*l_eos%dhdp))/(dv(3)*l_eos%dhdt)
        ic2v = (dv(4)*v_eos%drdp*v_eos%dhdt + v_eos%drdt*(1.d0-dv(4)*v_eos%dhdp))/(dv(4)*v_eos%dhdt)
        ic2g = (dv(5)*g_eos%drdp*g_eos%dhdt + g_eos%drdt*(1.d0-dv(5)*g_eos%dhdp))/(dv(5)*g_eos%dhdt)
        dv(6) = 1.d0/(dv(1)**2*(yl*ic2l*irhol**2+y1*ic2v*irhov**2+y2*ic2g*irhog**2))
      end if

      dv(15) = v_eos%drdp
      dv(16) = v_eos%drdt
      dv(17) = l_eos%drdp
      dv(18) = l_eos%drdt

    end subroutine deteos_simple
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine deteosonly(eos,p,t,y1,y2,dv,tv)
      implicit none
      class(t_eosonly), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(eos%ndv) ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      real(8), intent(inout) :: tv(eos%ntv)
      type(t_eos2)  :: l_eos,v_eos,g_eos
      real(8) :: irhol,irhov,irhog,irho,yl
      real(8) :: ic2l,ic2v,ic2g

      call eos%eos_l(p,t,1,l_eos)
      call eos%eos_v(p,t,2,v_eos)
      call eos%eos_g(p,t,3,g_eos)

      dv(3)  = l_eos%rho
      dv(4)  = v_eos%rho
      dv(5)  = g_eos%rho

      irhol = 1.d0/dv(3)
      irhov = 1.d0/dv(4)
      irhog = 1.d0/dv(5)
      yl = 1.d0-y1-y2
      irho = yl*irhol + y1*irhov + y2*irhog
      dv(1)  = 1.d0/irho
      dv(2)  = yl*l_eos%h + y1*v_eos%h + y2*g_eos%h

      dv(7)  = dv(1)**2*(yl*irhol**2*l_eos%drdp   &
                       + y1*irhov**2*v_eos%drdp   &
                       + y2*irhog**2*g_eos%drdp)
      dv(8)  = dv(1)**2*(yl*irhol**2*l_eos%drdt   &
                       + y1*irhov**2*v_eos%drdt   &
                       + y2*irhog**2*g_eos%drdt)
      dv(9)  = dv(1)**2*(irhol - irhov)
      dv(10) = dv(1)**2*(irhol - irhog)

      dv(11) = yl*l_eos%dhdp + y1*v_eos%dhdp + y2*g_eos%dhdp
      dv(12) = yl*l_eos%dhdt + y1*v_eos%dhdt + y2*g_eos%dhdt
      dv(13) = v_eos%h - l_eos%h
      dv(14) = g_eos%h - l_eos%h

      dv(6)  = dv(1)*dv(12)/(dv(1)*dv(7)*dv(12)+dv(8)*(1.d0-dv(1)*dv(11)))
      if(dv(6).le.0.d0) then
        write(*,*) 'alternative calculation for snd2'
        ic2l = (dv(3)*l_eos%drdp*l_eos%dhdt + l_eos%drdt*(1.d0-dv(3)*l_eos%dhdp))/(dv(3)*l_eos%dhdt)
        ic2v = (dv(4)*v_eos%drdp*v_eos%dhdt + v_eos%drdt*(1.d0-dv(4)*v_eos%dhdp))/(dv(4)*v_eos%dhdt)
        ic2g = (dv(5)*g_eos%drdp*g_eos%dhdt + g_eos%drdt*(1.d0-dv(5)*g_eos%dhdp))/(dv(5)*g_eos%dhdt)
        dv(6) = 1.d0/(dv(1)**2*(yl*ic2l*irhol**2+y1*ic2v*irhov**2+y2*ic2g*irhog**2))
      end if

      dv(15) = v_eos%drdp
      dv(16) = v_eos%drdt
      dv(17) = l_eos%drdp
      dv(18) = l_eos%drdt

    end subroutine deteosonly
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine deteos_prop(eos,p,t,y1,y2,dv,tv)
      implicit none
      class(t_eos_prop), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(eos%ndv) ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      real(8), intent(inout) :: tv(eos%ntv)
      type(t_eos2)  :: l_eos,v_eos,g_eos
      type(t_prop2) :: l_prop,v_prop,g_prop
      real(8) :: irhol,irhov,irhog,irho,yl
      real(8) :: ic2l,ic2v,ic2g

      call eos%eos_l_prop(p,t,1,l_eos,l_prop)
      call eos%eos_v_prop(p,t,2,v_eos,v_prop)
      call eos%eos_g_prop(p,t,3,g_eos,g_prop)

      dv(3)  = l_eos%rho
      dv(4)  = v_eos%rho
      dv(5)  = g_eos%rho

      irhol = 1.d0/dv(3)
      irhov = 1.d0/dv(4)
      irhog = 1.d0/dv(5)
      yl = 1.d0-y1-y2
      irho = yl*irhol + y1*irhov + y2*irhog
      dv(1)  = 1.d0/irho
      dv(2)  = yl*l_eos%h + y1*v_eos%h + y2*g_eos%h

      dv(7)  = dv(1)**2*(yl*irhol**2*l_eos%drdp   &
                       + y1*irhov**2*v_eos%drdp   &
                       + y2*irhog**2*g_eos%drdp)
      dv(8)  = dv(1)**2*(yl*irhol**2*l_eos%drdt   &
                       + y1*irhov**2*v_eos%drdt   &
                       + y2*irhog**2*g_eos%drdt)
      dv(9)  = dv(1)**2*(irhol - irhov)
      dv(10) = dv(1)**2*(irhol - irhog)

      dv(11) = yl*l_eos%dhdp + y1*v_eos%dhdp + y2*g_eos%dhdp
      dv(12) = yl*l_eos%dhdt + y1*v_eos%dhdt + y2*g_eos%dhdt
      dv(13) = v_eos%h - l_eos%h
      dv(14) = g_eos%h - l_eos%h

      dv(6)  = dv(1)*dv(12)/(dv(1)*dv(7)*dv(12)+dv(8)*(1.d0-dv(1)*dv(11)))
      if(dv(6).le.0.d0) then
        write(*,*) 'alternative calculation for snd2'
        ic2l = (dv(3)*l_eos%drdp*l_eos%dhdt + l_eos%drdt*(1.d0-dv(3)*l_eos%dhdp))/(dv(3)*l_eos%dhdt)
        ic2v = (dv(4)*v_eos%drdp*v_eos%dhdt + v_eos%drdt*(1.d0-dv(4)*v_eos%dhdp))/(dv(4)*v_eos%dhdt)
        ic2g = (dv(5)*g_eos%drdp*g_eos%dhdt + g_eos%drdt*(1.d0-dv(5)*g_eos%dhdp))/(dv(5)*g_eos%dhdt)
        dv(6) = 1.d0/(dv(1)**2*(yl*ic2l*irhol**2+y1*ic2v*irhov**2+y2*ic2g*irhog**2))
      end if

      dv(15) = v_eos%drdp
      dv(16) = v_eos%drdt
      dv(17) = l_eos%drdp
      dv(18) = l_eos%drdt

      tv(1) = dv(1)*(l_prop%vis*yl*irhol  + v_prop%vis*y1*irhov + g_prop%vis*y2*irhog)
      tv(2) = dv(1)*(l_prop%cond*yl*irhol + v_prop%cond*y1*irhov + g_prop%cond*y2*irhog)

    end subroutine deteos_prop
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine ideal(eos,p,t,phase,eos2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      real(8), parameter :: r1 = 0.00348383500557d0 !1.4d0/0.4d0/1004.64
      real(8) :: tau

      tau = 1.d0/t

      ! ideal gas law
      eos2%rho  = p*r1*tau
      eos2%h    = 1004.64d0*t
      eos2%drdp = r1*tau
      eos2%drdt = -p*r1*tau**2
      eos2%dhdp = 0.d0
      eos2%dhdt = 1004.64d0

    end subroutine ideal
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine ideal_prop(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      real(8), parameter :: r1 = 0.00348383500557d0 !1.4d0/0.4d0/1004.64
      real(8) :: tau

      tau = 1.d0/t

      ! ideal gas law
      eos2%rho  = p*r1*tau
      eos2%h    = 1004.64d0*t
      eos2%drdp = r1*tau
      eos2%drdt = -p*r1*tau**2
      eos2%dhdp = 0.d0
      eos2%dhdt = 1004.64d0

      !sutherland
      prop2%vis  = 1.716d-5*(t/273.11d0)**1.5d0*383.67d0/(t+110.56d0)
      !prop2%cond = 0.0241d0*(t/273.d0)**1.5d0*467.d0/(t+194.d0)
      prop2%cond = 1395.333333d0*prop2%vis !pr = 0.72

    end subroutine ideal_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine stiffened(eos,p,t,phase,eos2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      real(8), parameter :: r1 = 3.0915376183918d-4 !0.00037160906726d0 !4.4d0/3.4d0/4186.d0
      real(8) :: tau

      tau = 1.d0/t

      eos2%rho  = (p+6.d8)*r1*tau
      eos2%h    = 4186.d0*t!-917822.36d0  ! 917822.36??
      eos2%drdp = r1*tau
      eos2%drdt = -(p+6.d8)*r1*tau**2
      eos2%dhdp = 0.d0
      eos2%dhdt = 4186.d0

    end subroutine stiffened
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine stiffened_prop(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out):: prop2
      real(8), parameter :: r1 = 0.00037160906726d0 !2.8d0/1.8d0/4186.d0
      real(8) :: tau,temp

      tau = 1.d0/t

      eos2%rho  = (p+8.5d8)*r1*tau
      eos2%h    = 4186.d0*t-917822.36d0
      eos2%drdp = r1*tau
      eos2%drdt = -(p+8.5d8)*r1*tau**2
      eos2%dhdp = 0.d0
      eos2%dhdt = 4186.d0

      temp = dmax1(273.d0,t)

      prop2%vis  = 1.788d-3*dexp(-1.704d0 - 5.306d0*273.d0/temp + 7.003d0*(273.d0/temp)**2)
      !prop2%cond = -0.000009438d0*t**2 + 0.007294d0*t - 0.7284d0
      prop2%cond = 4598.d0*prop2%vis !pr = 7

    end subroutine stiffened_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97_l(eos,p,t,phase,eos2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      real(8) :: pi,tau,gt,gp,gpp,gtt,gpt
      real(8) :: g,g1,g2,g3,g4,g5,pi1,tau1,gp1
      integer :: k

      pi = p*0.00000006049607d0!/16.53d6
      tau = 1386d0/t

      pi1 = 1.d0/(7.1d0-pi)
      tau1 = 1.d0/(tau-1.222d0)
      gt = 0.d0
      gp = 0.d0
      gpp = 0.d0
      gtt = 0.d0
      gpt = 0.d0
      do k = 1,34
        g = eos%iapws_liquid_coeff%n(k)*(7.1d0-pi)**eos%iapws_liquid_coeff%i(k)*(tau-1.222d0)**eos%iapws_liquid_coeff%j(k)
        g1 = dble(eos%iapws_liquid_coeff%j(k))*g*tau1
        g2 = - dble(eos%iapws_liquid_coeff%i(k))*g*pi1
        g3 = - dble(eos%iapws_liquid_coeff%i(k)-1)*g2*pi1
        g4 = dble(eos%iapws_liquid_coeff%j(k)-1)*g1*tau1
        g5 = dble(eos%iapws_liquid_coeff%j(k))*g2*tau1
        gt = gt + g1
        gp = gp + g2
        gpp = gpp + g3
        gtt = gtt + g4
        gpt = gpt + g5
      end do

      gp1 = 1.d0/gp

      eos2%rho  = 25.841246054191803727*tau*gp1
      eos2%h    = 639675.036d0*gt
      eos2%drdp = -gpp*1.5632937721834122d-6*tau*gp1**2
      eos2%drdt = -eos2%rho*7.215007215007215d-4*tau + 0.01864447767257706*gpt*tau**3*gp1**2
      eos2%dhdp = 0.03869782431941923775*gpt
      eos2%dhdt = -461.526d0*tau**2*gtt

    end subroutine iapws97_l
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97_v(eos,p,t,phase,eos2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      real(8) :: pww
      real(8) :: pi,tau,g0t,g0p,g0pp,g0tt,g0pt,grt,grp,grpp,grtt,grpt
      real(8) :: g0,g01,g02
      real(8) :: gr,gr1,gr2,gr3,gr4,gr5,pi1,tau1,tau2,gp1
      integer :: k

      pww = eos%get_pww(t)

      pi = p*1d-6
      tau = 540d0/t

      pi1 = 1.d6/p
      tau1 = 1.d0/tau
      tau2 = 1.d0/(tau-0.5d0)

      g0t = 0.d0
      g0p = pi1
      g0pp = -pi1**2
      g0tt = 0.d0
      g0pt = 0.d0
      grt = 0.d0
      grp = 0.d0
      grpp = 0.d0
      grtt = 0.d0
      grpt = 0.d0

      if(p.gt.pww) then
        do k = 1,13
          if(k.le.9) then
            g0 = eos%iapws_m_vapor_coeff%n(k)*tau**eos%iapws_m_vapor_coeff%j(k)
            g01 = dble(eos%iapws_m_vapor_coeff%j(k))*g0*tau1
            g02 = dble(eos%iapws_m_vapor_coeff%j(k)-1)*g01*tau1
            g0t = g0t + g01
            g0tt = g0tt + g02
          end if
          gr = eos%iapws_m_vapor_coeff%nr(k)*pi**eos%iapws_m_vapor_coeff%ir(k)*(tau-0.5d0)**eos%iapws_m_vapor_coeff%jr(k)
          gr1 = dble(eos%iapws_m_vapor_coeff%jr(k))*gr*tau2
          gr2 = dble(eos%iapws_m_vapor_coeff%ir(k))*gr*pi1
          gr3 = dble(eos%iapws_m_vapor_coeff%ir(k)-1)*gr2*pi1
          gr4 = dble(eos%iapws_m_vapor_coeff%jr(k)-1)*gr1*tau2
          gr5 = dble(eos%iapws_m_vapor_coeff%jr(k))*gr2*tau2
          grt = grt + gr1
          grp = grp + gr2
          grpp = grpp + gr3
          grtt = grtt + gr4
          grpt = grpt + gr5
        end do
      else
        do k = 1,43
          if(k.le.9) then
            g0 = eos%iapws_vapor_coeff%n(k)*tau**eos%iapws_vapor_coeff%j(k)
            g01 = dble(eos%iapws_vapor_coeff%j(k))*g0*tau1
            g02 = dble(eos%iapws_vapor_coeff%j(k)-1)*g01*tau1
            g0t = g0t + g01
            g0tt = g0tt + g02
          end if
          gr = eos%iapws_vapor_coeff%nr(k)*pi**eos%iapws_vapor_coeff%ir(k)*(tau-0.5d0)**eos%iapws_vapor_coeff%jr(k)
          gr1 = dble(eos%iapws_vapor_coeff%jr(k))*gr*tau2
          gr2 = dble(eos%iapws_vapor_coeff%ir(k))*gr*pi1
          gr3 = dble(eos%iapws_vapor_coeff%ir(k)-1)*gr2*pi1
          gr4 = dble(eos%iapws_vapor_coeff%jr(k)-1)*gr1*tau2
          gr5 = dble(eos%iapws_vapor_coeff%jr(k))*gr2*tau2
          grt = grt + gr1
          grp = grp + gr2
          grpp = grpp + gr3
          grtt = grtt + gr4
          grpt = grpt + gr5
        end do
      end if

      gp1 = 1.d0/(g0p+grp)

      eos2%rho  = 4.01245401527075799d0*tau*gp1
      eos2%h    = 249224.04d0*(g0t+grt)
      eos2%drdp = -(g0pp+grpp)*4.01245401527075799d-6*tau*gp1**2
      eos2%drdt = -eos2%rho*0.00185185185185185185*tau + 0.0074304703986495518*(g0pt+grpt)*tau**3*gp1**2
      eos2%dhdp = 0.24922404d0*(g0pt+grpt)
      eos2%dhdt = -461.526d0*tau**2*(g0tt+grtt)

    end subroutine iapws97_v
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97_l_prop(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      real(8) :: pi,tau,gt,gp,gpp,gtt,gpt
      real(8) :: g,g1,g2,g3,g4,g5,pi1,tau1,gp1
      real(8) :: tt,rr,mu0,mu1,mu2,lam0,lam1,lam2,tt1
      integer :: i,j,k

      pi = p*0.00000006049607d0!/16.53d6
      tau = 1386d0/t

      pi1 = 1.d0/(7.1d0-pi)
      tau1 = 1.d0/(tau-1.222d0)
      gt = 0.d0
      gp = 0.d0
      gpp = 0.d0
      gtt = 0.d0
      gpt = 0.d0
      do k = 1,34
        g = eos%iapws_liquid_coeff%n(k)*(7.1d0-pi)**eos%iapws_liquid_coeff%i(k)*(tau-1.222d0)**eos%iapws_liquid_coeff%j(k)
        g1 = dble(eos%iapws_liquid_coeff%j(k))*g*tau1
        g2 = - dble(eos%iapws_liquid_coeff%i(k))*g*pi1
        g3 = - dble(eos%iapws_liquid_coeff%i(k)-1)*g2*pi1
        g4 = dble(eos%iapws_liquid_coeff%j(k)-1)*g1*tau1
        g5 = dble(eos%iapws_liquid_coeff%j(k))*g2*tau1
        gt = gt + g1
        gp = gp + g2
        gpp = gpp + g3
        gtt = gtt + g4
        gpt = gpt + g5
      end do

      gp1 = 1.d0/gp

      eos2%rho  = 25.841246054191803727*tau*gp1
      eos2%h    = 639675.036d0*gt
      eos2%drdp = -gpp*1.5632937721834122d-6*tau*gp1**2
      eos2%drdt = -eos2%rho*7.215007215007215d-4*tau + 0.01864447767257706*gpt*tau**3*gp1**2
      eos2%dhdp = 0.03869782431941923775*gpt
      eos2%dhdt = -461.526d0*tau**2*gtt

      tt = t*0.00154495032985d0!/647.27d0
      rr = dmax1(1.d-3,dmin1(1332.4d0,eos2%rho))*0.00314699949333d0!/317.763d0
      tt1 = 647.27d0/t

      mu0 = dsqrt(tt)/(eos%iapws_prop_coeff%h(0) + eos%iapws_prop_coeff%h(1)*tt1  &
                  + eos%iapws_prop_coeff%h(2)*tt1**2 + eos%iapws_prop_coeff%h(3)*tt1**3)

      mu1 = 0.d0
      do i = 0,5
        do j = 0,6
          mu1 = mu1 + eos%iapws_prop_coeff%hh(i,j)*(tt1-1.d0)**i*(rr-1.d0)**j
        end do
      end do
      mu1 = dexp(rr*mu1)
      mu2 = 1.d0

      prop2%vis = 55.071d-6*mu0*mu1*mu2

      lam0 = dsqrt(tt)/(eos%iapws_prop_coeff%l(0) + eos%iapws_prop_coeff%l(1)*tt1  &
                   + eos%iapws_prop_coeff%l(2)*tt1**2 + eos%iapws_prop_coeff%l(3)*tt1**3)
      lam1 = 0.d0
      do i = 0,4
        do j = 0,5
          lam1 = lam1 + eos%iapws_prop_coeff%ll(i,j)*(tt1-1.d0)**i*(rr-1.d0)**j
        end do
      end do
      lam1 = dexp(rr*lam1)
      lam2 = 0.d0

      prop2%cond = 0.4945d0*(lam0*lam1 + lam2)

    end subroutine iapws97_l_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97_v_prop(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      real(8) :: pww
      real(8) :: pi,tau,g0t,g0p,g0pp,g0tt,g0pt,grt,grp,grpp,grtt,grpt
      real(8) :: g0,g01,g02
      real(8) :: gr,gr1,gr2,gr3,gr4,gr5,pi1,tau1,tau2,gp1
      real(8) :: tt,rr,mu0,mu1,mu2,lam0,lam1,lam2,tt1
      integer :: i,j,k

      pww = eos%get_pww(t)

      pi = p*1d-6
      tau = 540d0/t

      pi1 = 1.d6/p
      tau1 = 1.d0/tau
      tau2 = 1.d0/(tau-0.5d0)

      g0t = 0.d0
      g0p = pi1
      g0pp = -pi1**2
      g0tt = 0.d0
      g0pt = 0.d0
      grt = 0.d0
      grp = 0.d0
      grpp = 0.d0
      grtt = 0.d0
      grpt = 0.d0

      if(p.gt.pww) then
        do k = 1,13
          if(k.le.9) then
            g0 = eos%iapws_m_vapor_coeff%n(k)*tau**eos%iapws_m_vapor_coeff%j(k)
            g01 = dble(eos%iapws_m_vapor_coeff%j(k))*g0*tau1
            g02 = dble(eos%iapws_m_vapor_coeff%j(k)-1)*g01*tau1
            g0t = g0t + g01
            g0tt = g0tt + g02
          end if
          gr = eos%iapws_m_vapor_coeff%nr(k)*pi**eos%iapws_m_vapor_coeff%ir(k)*(tau-0.5d0)**eos%iapws_m_vapor_coeff%jr(k)
          gr1 = dble(eos%iapws_m_vapor_coeff%jr(k))*gr*tau2
          gr2 = dble(eos%iapws_m_vapor_coeff%ir(k))*gr*pi1
          gr3 = dble(eos%iapws_m_vapor_coeff%ir(k)-1)*gr2*pi1
          gr4 = dble(eos%iapws_m_vapor_coeff%jr(k)-1)*gr1*tau2
          gr5 = dble(eos%iapws_m_vapor_coeff%jr(k))*gr2*tau2
          grt = grt + gr1
          grp = grp + gr2
          grpp = grpp + gr3
          grtt = grtt + gr4
          grpt = grpt + gr5
        end do
      else
        do k = 1,43
          if(k.le.9) then
            g0 = eos%iapws_vapor_coeff%n(k)*tau**eos%iapws_vapor_coeff%j(k)
            g01 = dble(eos%iapws_vapor_coeff%j(k))*g0*tau1
            g02 = dble(eos%iapws_vapor_coeff%j(k)-1)*g01*tau1
            g0t = g0t + g01
            g0tt = g0tt + g02
          end if
          gr = eos%iapws_vapor_coeff%nr(k)*pi**eos%iapws_vapor_coeff%ir(k)*(tau-0.5d0)**eos%iapws_vapor_coeff%jr(k)
          gr1 = dble(eos%iapws_vapor_coeff%jr(k))*gr*tau2
          gr2 = dble(eos%iapws_vapor_coeff%ir(k))*gr*pi1
          gr3 = dble(eos%iapws_vapor_coeff%ir(k)-1)*gr2*pi1
          gr4 = dble(eos%iapws_vapor_coeff%jr(k)-1)*gr1*tau2
          gr5 = dble(eos%iapws_vapor_coeff%jr(k))*gr2*tau2
          grt = grt + gr1
          grp = grp + gr2
          grpp = grpp + gr3
          grtt = grtt + gr4
          grpt = grpt + gr5
        end do
      end if

      gp1 = 1.d0/(g0p+grp)

      eos2%rho  = 4.01245401527075799d0*tau*gp1
      eos2%h    = 249224.04d0*(g0t+grt)
      eos2%drdp = -(g0pp+grpp)*4.01245401527075799d-6*tau*gp1**2
      eos2%drdt = -eos2%rho*0.00185185185185185185*tau + 0.0074304703986495518*(g0pt+grpt)*tau**3*gp1**2
      eos2%dhdp = 0.24922404d0*(g0pt+grpt)
      eos2%dhdt = -461.526d0*tau**2*(g0tt+grtt)

      tt = t*0.00154495032985d0!/647.27d0
      rr = dmax1(1.d-3,dmin1(1332.4d0,eos2%rho))*0.00314699949333d0!/317.763d0
      tt1 = 647.27d0/t

      mu0 = dsqrt(tt)/(eos%iapws_prop_coeff%h(0) + eos%iapws_prop_coeff%h(1)*tt1  &
                  + eos%iapws_prop_coeff%h(2)*tt1**2 + eos%iapws_prop_coeff%h(3)*tt1**3)

      mu1 = 0.d0
      do i = 0,5
        do j = 0,6
          mu1 = mu1 + eos%iapws_prop_coeff%hh(i,j)*(tt1-1.d0)**i*(rr-1.d0)**j
        end do
      end do
      mu1 = dexp(rr*mu1)
      mu2 = 1.d0

      prop2%vis = 55.071d-6*mu0*mu1*mu2

      lam0 = dsqrt(tt)/(eos%iapws_prop_coeff%l(0) + eos%iapws_prop_coeff%l(1)*tt1  &
                   + eos%iapws_prop_coeff%l(2)*tt1**2 + eos%iapws_prop_coeff%l(3)*tt1**3)
      lam1 = 0.d0
      do i = 0,4
        do j = 0,5
          lam1 = lam1 + eos%iapws_prop_coeff%ll(i,j)*(tt1-1.d0)**i*(rr-1.d0)**j
        end do
      end do
      lam1 = dexp(rr*lam1)
      lam2 = 0.d0

      prop2%cond = 0.4945d0*(lam0*lam1 + lam2)

    end subroutine iapws97_v_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine database(eos,p,t,phase,eos2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      integer :: t_index,p_index,l
      integer :: n1,n2,m1,m2,m3,m4,t_subn,p_subn1,p_subn2
      real(8) :: t1,t2,p1,p2
      real(8) :: x(2),y(4)
      real(8) :: xi,eta,dp,dt,output(3),coeff(9)

      t_subn = 0
      p_subn1 = 0
      p_subn2 = 0


      ! defining the temperature index
      do l=1,eos%db_const(phase)%t_nsub
        t1 = eos%db_const(phase)%t_boundary(l)
        t2 = eos%db_const(phase)%t_boundary(l+1)
        if ((t-t1)*(t-t2).le.0.d0) then
          t_subn = l
          exit
        end if
      end do

      if(t_subn.eq.0) then
        if(t.le.eos%db_const(phase)%t_boundary(1)) then
          n1 = 1
          n2 = 2
        else if(t.ge.eos%db_const(phase)%t_boundary(eos%db_const(phase)%t_nsub+1)) then
          n1 = eos%db_const(phase)%t_ndata-1
          n2 = eos%db_const(phase)%t_ndata
        end if
      else
        t_index = int((t-eos%db_const(phase)%t_boundary(t_subn))/eos%db_const(phase)%t_delta(t_subn)) &
                + eos%db_const(phase)%t_iboundary(t_subn)
        n1 = min0(t_index,eos%db_const(phase)%t_iboundary(t_subn+1)-1)
        n2 = min0(t_index+1,eos%db_const(phase)%t_iboundary(t_subn+1))
      end if

      x(1) = eos%db_const(phase)%db(n1)%t
      x(2) = eos%db_const(phase)%db(n2)%t


      ! defining the pressure index 1
      do l=1,eos%db_const(phase)%db(n1)%p_nsub
        p1 = eos%db_const(phase)%db(n1)%p_boundary(l)
        p2 = eos%db_const(phase)%db(n1)%p_boundary(l+1)
        if ((p-p1)*(p-p2).le.0.d0) then
          p_subn1 = l
          exit
        end if
      end do

      if(p_subn1.eq.0) then
        if(p.le.eos%db_const(phase)%db(n1)%p_boundary(1)) then
          m1 = eos%db_const(phase)%db(n1)%p_nstart
          m2 = eos%db_const(phase)%db(n1)%p_nstart+1
          p_subn1 = 1
        else if(p.ge.eos%db_const(phase)%db(n1)%p_boundary(eos%db_const(phase)%db(n1)%p_nsub+1)) then
          m1 = eos%db_const(phase)%db(n1)%p_iboundary(eos%db_const(phase)%db(n1)%p_nsub+1)-1
          m2 = eos%db_const(phase)%db(n1)%p_iboundary(eos%db_const(phase)%db(n1)%p_nsub+1)
          p_subn1 = eos%db_const(phase)%db(n1)%p_nsub
        end if
      else
        p_index = int((p-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) &
                + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
        m1 = min0(p_index,eos%db_const(phase)%db(n1)%p_iboundary(p_subn1+1)-1)
        m2 = min0(p_index+1,eos%db_const(phase)%db(n1)%p_iboundary(p_subn1+1))
      end if

      y(1) = eos%db_const(phase)%db(n1)%p(m1)
      y(2) = eos%db_const(phase)%db(n1)%p(m2)


      ! defining the pressure index 2
      do l=1,eos%db_const(phase)%db(n2)%p_nsub
        p1 = eos%db_const(phase)%db(n2)%p_boundary(l)
        p2 = eos%db_const(phase)%db(n2)%p_boundary(l+1)
        if ((p-p1)*(p-p2).le.0.d0) then
          p_subn2 = l
          exit
        end if
      end do

      if(p_subn2.eq.0) then
        if(p.le.eos%db_const(phase)%db(n2)%p_boundary(1)) then
          m4 = eos%db_const(phase)%db(n2)%p_nstart
          m3 = eos%db_const(phase)%db(n2)%p_nstart+1
          p_subn2 = 1
        else if(p.ge.eos%db_const(phase)%db(n2)%p_boundary(eos%db_const(phase)%db(n2)%p_nsub+1)) then
          m4 = eos%db_const(phase)%db(n2)%p_iboundary(eos%db_const(phase)%db(n2)%p_nsub+1)-1
          m3 = eos%db_const(phase)%db(n2)%p_iboundary(eos%db_const(phase)%db(n2)%p_nsub+1)
          p_subn2 = eos%db_const(phase)%db(n2)%p_nsub
        end if
      else
        p_index = int((p-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) &
                + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
        m4 = min0(p_index,eos%db_const(phase)%db(n2)%p_iboundary(p_subn2+1)-1)
        m3 = min0(p_index+1,eos%db_const(phase)%db(n2)%p_iboundary(p_subn2+1))
      end if

      y(4) = eos%db_const(phase)%db(n2)%p(m4)
      y(3) = eos%db_const(phase)%db(n2)%p(m3)


      ! construct regular grid
      ! m1~m4: clockwise from the lower left point (m1)
      if (y(1) .lt. y(4)) then
        m4 = int((y(1)-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
      else
        m1 = int((y(4)-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
      end if

      if (y(2) .gt. y(3)) then
        m3 = int((y(2)-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
      else
        m2 = int((y(3)-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
      end if

      y(1) = eos%db_const(phase)%db(n1)%p(m1)
      y(2) = eos%db_const(phase)%db(n2)%p(m3)


      ! biquadratic interpolation
      dt = x(2)-x(1)
      dp = y(2)-y(1)

      xi  = dmin1(dmax1((t-x(1))/(x(2)-x(1)),0.d0),1.d0)
      eta = dmin1(dmax1((p-y(1))/(y(2)-y(1)),0.d0),1.d0)

      coeff = eos%db_const(phase)%db(n1)%rho(m1,:)
      call biquadratic_interpolation(coeff,xi,eta,dt,dp,output)
      eos2%rho  = output(1)
      eos2%drdt = output(2)
      eos2%drdp = output(3)

      coeff = eos%db_const(phase)%db(n1)%h(m1,:)
      call biquadratic_interpolation(coeff,xi,eta,dt,dp,output)
      eos2%h    = output(1)
      eos2%dhdt = output(2)
      eos2%dhdp = output(3)

    end subroutine database
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine database_prop(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      integer :: t_index,p_index,l
      integer :: n1,n2,m1,m2,m3,m4,t_subn,p_subn1,p_subn2
      real(8) :: t1,t2,p1,p2
      real(8) :: x(2),y(4)
      real(8) :: xi,eta,dp,dt,output(3),coeff(9)

      t_subn = 0
      p_subn1 = 0
      p_subn2 = 0


      ! defining the temperature index
      do l=1,eos%db_const(phase)%t_nsub
        t1 = eos%db_const(phase)%t_boundary(l)
        t2 = eos%db_const(phase)%t_boundary(l+1)
        if ((t-t1)*(t-t2).le.0.d0) then
          t_subn = l
          exit
        end if
      end do

      if(t_subn.eq.0) then
        if(t.le.eos%db_const(phase)%t_boundary(1)) then
          n1 = 1
          n2 = 2
        else if(t.ge.eos%db_const(phase)%t_boundary(eos%db_const(phase)%t_nsub+1)) then
          n1 = eos%db_const(phase)%t_ndata-1
          n2 = eos%db_const(phase)%t_ndata
        end if
      else
        t_index = int((t-eos%db_const(phase)%t_boundary(t_subn))/eos%db_const(phase)%t_delta(t_subn)) &
                + eos%db_const(phase)%t_iboundary(t_subn)
        n1 = min0(t_index,eos%db_const(phase)%t_iboundary(t_subn+1)-1)
        n2 = min0(t_index+1,eos%db_const(phase)%t_iboundary(t_subn+1))
      end if

      x(1) = eos%db_const(phase)%db(n1)%t
      x(2) = eos%db_const(phase)%db(n2)%t


      ! defining the pressure index 1
      do l=1,eos%db_const(phase)%db(n1)%p_nsub
        p1 = eos%db_const(phase)%db(n1)%p_boundary(l)
        p2 = eos%db_const(phase)%db(n1)%p_boundary(l+1)
        if ((p-p1)*(p-p2).le.0.d0) then
          p_subn1 = l
          exit
        end if
      end do

      if(p_subn1.eq.0) then
        if(p.le.eos%db_const(phase)%db(n1)%p_boundary(1)) then
          m1 = eos%db_const(phase)%db(n1)%p_nstart
          m2 = eos%db_const(phase)%db(n1)%p_nstart+1
          p_subn1 = 1
        else if(p.ge.eos%db_const(phase)%db(n1)%p_boundary(eos%db_const(phase)%db(n1)%p_nsub+1)) then
          m1 = eos%db_const(phase)%db(n1)%p_iboundary(eos%db_const(phase)%db(n1)%p_nsub+1)-1
          m2 = eos%db_const(phase)%db(n1)%p_iboundary(eos%db_const(phase)%db(n1)%p_nsub+1)
          p_subn1 = eos%db_const(phase)%db(n1)%p_nsub
        end if
      else
        p_index = int((p-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) &
                + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
        m1 = min0(p_index,eos%db_const(phase)%db(n1)%p_iboundary(p_subn1+1)-1)
        m2 = min0(p_index+1,eos%db_const(phase)%db(n1)%p_iboundary(p_subn1+1))
      end if

      y(1) = eos%db_const(phase)%db(n1)%p(m1)
      y(2) = eos%db_const(phase)%db(n1)%p(m2)


      ! defining the pressure index 2
      do l=1,eos%db_const(phase)%db(n2)%p_nsub
        p1 = eos%db_const(phase)%db(n2)%p_boundary(l)
        p2 = eos%db_const(phase)%db(n2)%p_boundary(l+1)
        if ((p-p1)*(p-p2).le.0.d0) then
          p_subn2 = l
          exit
        end if
      end do

      if(p_subn2.eq.0) then
        if(p.le.eos%db_const(phase)%db(n2)%p_boundary(1)) then
          m4 = eos%db_const(phase)%db(n2)%p_nstart
          m3 = eos%db_const(phase)%db(n2)%p_nstart+1
          p_subn2 = 1
        else if(p.ge.eos%db_const(phase)%db(n2)%p_boundary(eos%db_const(phase)%db(n2)%p_nsub+1)) then
          m4 = eos%db_const(phase)%db(n2)%p_iboundary(eos%db_const(phase)%db(n2)%p_nsub+1)-1
          m3 = eos%db_const(phase)%db(n2)%p_iboundary(eos%db_const(phase)%db(n2)%p_nsub+1)
          p_subn2 = eos%db_const(phase)%db(n2)%p_nsub
        end if
      else
        p_index = int((p-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) &
                + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
        m4 = min0(p_index,eos%db_const(phase)%db(n2)%p_iboundary(p_subn2+1)-1)
        m3 = min0(p_index+1,eos%db_const(phase)%db(n2)%p_iboundary(p_subn2+1))
      end if

      y(4) = eos%db_const(phase)%db(n2)%p(m4)
      y(3) = eos%db_const(phase)%db(n2)%p(m3)

      ! construct regular grid
      ! m1~m4: clockwise from the lower left point (m1)
      if (y(1) .lt. y(4)) then
        m4 = int((y(1)-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
      else
        m1 = int((y(4)-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
      end if

      if (y(2) .gt. y(3)) then
        m3 = int((y(2)-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
      else
        m2 = int((y(3)-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
      end if

      y(1) = eos%db_const(phase)%db(n1)%p(m1)
      y(2) = eos%db_const(phase)%db(n2)%p(m3)


      ! biquadratic interpolation
      dt = x(2)-x(1)
      dp = y(2)-y(1)

      xi  = dmin1(dmax1((t-x(1))/(x(2)-x(1)),0.d0),1.d0)
      eta = dmin1(dmax1((p-y(1))/(y(2)-y(1)),0.d0),1.d0)

      coeff = eos%db_const(phase)%db(n1)%rho(m1,:)
      call biquadratic_interpolation(coeff,xi,eta,dt,dp,output)
      eos2%rho  = output(1)
      eos2%drdt = output(2)
      eos2%drdp = output(3)

      coeff = eos%db_const(phase)%db(n1)%h(m1,:)
      call biquadratic_interpolation(coeff,xi,eta,dt,dp,output)
      eos2%h    = output(1)
      eos2%dhdt = output(2)
      eos2%dhdp = output(3)


      ! bilinear interpolation
      coeff(1:4) = eos%db_const(phase)%db(n1)%vis(m1,:)
      call bilinear_interpolation(coeff(1:4),xi,eta,dt,dp,output)
      prop2%vis = output(1)

      coeff(1:4) = eos%db_const(phase)%db(n1)%cond(m1,:)
      call bilinear_interpolation(coeff(1:4),xi,eta,dt,dp,output)
      prop2%cond = output(1)

    end subroutine database_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine srk_g(eos,p,t,phase,eos2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      real(8),parameter :: ru=8.314d0
      real(8) :: a,alpha,b,mw,s,tc,pc,tb,cv0
      real(8) :: rho,e,h,pv
      real(8) :: dpdt,dpdr
      real(8) :: daadt,dsaadt,daidt
      real(8) :: ddaaddt,ddsaaddt,ddaiddt
      real(8) :: cv,cp,bp,bt,ap,ar

      if(phase.eq.2) then
        tc   = eos%srk_vapor_property%tc
        pc   = eos%srk_vapor_property%pc
        tb   = eos%srk_vapor_property%tb
        mw   = eos%srk_vapor_property%mw
        cv0  = eos%srk_vapor_property%cv0
        s    = eos%srk_vapor_property%s
        a    = eos%srk_vapor_property%a
        b    = eos%srk_vapor_property%b
      else if(phase.eq.3) then
        tc   = eos%srk_gas_property%tc
        pc   = eos%srk_gas_property%pc
        tb   = eos%srk_gas_property%tb
        mw   = eos%srk_gas_property%mw
        cv0  = eos%srk_gas_property%cv0
        s    = eos%srk_gas_property%s
        a    = eos%srk_gas_property%a
        b    = eos%srk_gas_property%b
      end if

      alpha = (1.d0 + s*(1.d0-dsqrt(t/tc)))**2

      pv = computepv(tc,pc,tb,t)
      !if ((p.le.pv).and.(t.le.tc)) then
      if ((p.le.pv)) then
        alpha = pv/pc
      end if

      rho = cubicsolver(p,t,a*alpha,b,mw)

      ! find daadt for dpdt
      daidt  = -s/dsqrt(t*tc)*(1.d0+s*(1.d0-dsqrt(t/tc)))
      dsaadt = daidt
      daadt  = a*dsaadt

      ! dpdt, dpdr
      dpdt = rho*ru/(mw-b*rho) - daadt*rho**2/(mw*(mw+b*rho))
      dpdr = mw*ru*t/(mw-b*rho)**2 - a*alpha*rho*(2.d0*mw+b*rho)/(mw*(mw+b*rho)**2)

      ! find ddaaddt for cv
      ddaiddt  = 0.5d0*s* (s/(t*tc) + (1.d0+s*(1.d0-dsqrt(t/tc)))/dsqrt(t**3*tc))
      ddsaaddt = ddaiddt
      ddaaddt  = a*ddsaaddt

      ! cv for dhdt & e
      cv = cv0 + t*ddaaddt*dlog(1.d0+b*rho/mw)/(b*mw)
      cp = cv + t*dpdt**2/(rho**2*dpdr)

      ! dhdp, dhdt
      bp = 1.d0/rho - p/(rho**2*dpdr)
      bt = cv + dpdt*p/(rho**2*dpdr)

      ap = cv/dpdt+1.d0/rho
      ar = -cv/dpdt*dpdr - p/rho**2

      e = cv0*t + (t*daadt-a*alpha)*dlog(1.d0+b*rho/mw)/(b*mw)
      h = e + dpdr*(1.d0-rho*bp)

      eos2%rho  = rho
      eos2%h    = h
      eos2%drdp = 1.d0/dpdr
      eos2%drdt = -dpdt/dpdr
      eos2%dhdp = bp
      eos2%dhdt = cp

    end subroutine srk_g
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine srk_g_prop(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      real(8),parameter :: ru=8.314d0
      real(8) :: a,alpha,b,mw,s,tc,pc,tb,cv0
      real(8) :: rho,e,h,pv
      real(8) :: dpdt,dpdr
      real(8) :: daadt,dsaadt,daidt
      real(8) :: ddaaddt,ddsaaddt,ddaiddt
      real(8) :: cv,cp,bp,bt,ap,ar
      integer :: t_index,p_index,l
      integer :: n1,n2,m1,m2,m3,m4,t_subn,p_subn1,p_subn2
      real(8) :: t1,t2,p1,p2
      real(8) :: x(2),y(4)
      real(8) :: xi,eta,dp,dt,output(3),coeff(4)

      ! srk for eos2

      if(phase.eq.2) then
        tc   = eos%srk_vapor_property%tc
        pc   = eos%srk_vapor_property%pc
        tb   = eos%srk_vapor_property%tb
        mw   = eos%srk_vapor_property%mw
        cv0  = eos%srk_vapor_property%cv0
        s    = eos%srk_vapor_property%s
        a    = eos%srk_vapor_property%a
        b    = eos%srk_vapor_property%b
     else if(phase.eq.3) then
        tc   = eos%srk_gas_property%tc
        pc   = eos%srk_gas_property%pc
        tb   = eos%srk_gas_property%tb
        mw   = eos%srk_gas_property%mw
        cv0  = eos%srk_gas_property%cv0
        s    = eos%srk_gas_property%s
        a    = eos%srk_gas_property%a
        b    = eos%srk_gas_property%b
      end if

      alpha = (1.d0 + s*(1.d0-dsqrt(t/tc)))**2

      pv = computepv(tc,pc,tb,t)
      !if ((p.le.pv).and.(t.le.tc)) then
      if ((p.le.pv)) then
        alpha = pv/pc
      end if

      rho = cubicsolver(p,t,a*alpha,b,mw)

      ! find daadt for dpdt
      daidt  = -s/dsqrt(t*tc)*(1.d0+s*(1.d0-dsqrt(t/tc)))
      dsaadt = daidt
      daadt  = a*dsaadt

      ! dpdt, dpdr
      dpdt = rho*ru/(mw-b*rho) - daadt*rho**2/(mw*(mw+b*rho))
      dpdr = mw*ru*t/(mw-b*rho)**2 - a*alpha*rho*(2.d0*mw+b*rho)/(mw*(mw+b*rho)**2)

      ! find ddaaddt for cv
      ddaiddt  = 0.5d0*s* (s/(t*tc) + (1.d0+s*(1.d0-dsqrt(t/tc)))/dsqrt(t**3*tc))
      ddsaaddt = ddaiddt
      ddaaddt  = a*ddsaaddt

      ! cv for dhdt & e
      cv = cv0 + t*ddaaddt*dlog(1.d0+b*rho/mw)/(b*mw)
      cp = cv + t*dpdt**2/(rho**2*dpdr)

      ! dhdp, dhdt
      bp = 1.d0/rho - p/(rho**2*dpdr)
      bt = cv + dpdt*p/(rho**2*dpdr)

      ap = cv/dpdt+1.d0/rho
      ar = -cv/dpdt*dpdr - p/rho**2

      e = cv0*t + (t*daadt-a*alpha)*dlog(1.d0+b*rho/mw)/(b*mw)
      h = e + dpdr*(1.d0-rho*bp)

      eos2%rho  = rho
      eos2%h    = h
      eos2%drdp = 1.d0/dpdr
      eos2%drdt = -dpdt/dpdr
      eos2%dhdp = bp
      eos2%dhdt = cp


      ! database for prop2

      t_subn = 0
      p_subn1 = 0
      p_subn2 = 0


      ! defining the temperature index
      do l=1,eos%db_const(phase)%t_nsub
        t1 = eos%db_const(phase)%t_boundary(l)
        t2 = eos%db_const(phase)%t_boundary(l+1)
        if ((t-t1)*(t-t2).le.0.d0) then
          t_subn = l
          exit
        end if
      end do

      if(t_subn.eq.0) then
        if(t.le.eos%db_const(phase)%t_boundary(1)) then
          n1 = 1
          n2 = 2
        else if(t.ge.eos%db_const(phase)%t_boundary(eos%db_const(phase)%t_nsub+1)) then
          n1 = eos%db_const(phase)%t_ndata-1
          n2 = eos%db_const(phase)%t_ndata
        end if
      else
        t_index = int((t-eos%db_const(phase)%t_boundary(t_subn))/eos%db_const(phase)%t_delta(t_subn)) &
                + eos%db_const(phase)%t_iboundary(t_subn)
        n1 = min0(t_index,eos%db_const(phase)%t_iboundary(t_subn+1)-1)
        n2 = min0(t_index+1,eos%db_const(phase)%t_iboundary(t_subn+1))
      end if

      x(1) = eos%db_const(phase)%db(n1)%t
      x(2) = eos%db_const(phase)%db(n2)%t


      ! defining the pressure index 1
      do l=1,eos%db_const(phase)%db(n1)%p_nsub
        p1 = eos%db_const(phase)%db(n1)%p_boundary(l)
        p2 = eos%db_const(phase)%db(n1)%p_boundary(l+1)
        if ((p-p1)*(p-p2).le.0.d0) then
          p_subn1 = l
          exit
        end if
      end do

      if(p_subn1.eq.0) then
        if(p.le.eos%db_const(phase)%db(n1)%p_boundary(1)) then
          m1 = eos%db_const(phase)%db(n1)%p_nstart
          m2 = eos%db_const(phase)%db(n1)%p_nstart+1
          p_subn1 = 1
        else if(p.ge.eos%db_const(phase)%db(n1)%p_boundary(eos%db_const(phase)%db(n1)%p_nsub+1)) then
          m1 = eos%db_const(phase)%db(n1)%p_iboundary(eos%db_const(phase)%db(n1)%p_nsub+1)-1
          m2 = eos%db_const(phase)%db(n1)%p_iboundary(eos%db_const(phase)%db(n1)%p_nsub+1)
          p_subn1 = eos%db_const(phase)%db(n1)%p_nsub
        end if
      else
        p_index = int((p-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) &
                + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
        m1 = min0(p_index,eos%db_const(phase)%db(n1)%p_iboundary(p_subn1+1)-1)
        m2 = min0(p_index+1,eos%db_const(phase)%db(n1)%p_iboundary(p_subn1+1))
      end if

      y(1) = eos%db_const(phase)%db(n1)%p(m1)
      y(2) = eos%db_const(phase)%db(n1)%p(m2)


      ! defining the pressure index 2
      do l=1,eos%db_const(phase)%db(n2)%p_nsub
        p1 = eos%db_const(phase)%db(n2)%p_boundary(l)
        p2 = eos%db_const(phase)%db(n2)%p_boundary(l+1)
        if ((p-p1)*(p-p2).le.0.d0) then
          p_subn2 = l
          exit
        end if
      end do

      if(p_subn2.eq.0) then
        if(p.le.eos%db_const(phase)%db(n2)%p_boundary(1)) then
          m4 = eos%db_const(phase)%db(n2)%p_nstart
          m3 = eos%db_const(phase)%db(n2)%p_nstart+1
          p_subn2 = 1
        else if(p.ge.eos%db_const(phase)%db(n2)%p_boundary(eos%db_const(phase)%db(n2)%p_nsub+1)) then
          m4 = eos%db_const(phase)%db(n2)%p_iboundary(eos%db_const(phase)%db(n2)%p_nsub+1)-1
          m3 = eos%db_const(phase)%db(n2)%p_iboundary(eos%db_const(phase)%db(n2)%p_nsub+1)
          p_subn2 = eos%db_const(phase)%db(n2)%p_nsub
        end if
      else
        p_index = int((p-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) &
                + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
        m4 = min0(p_index,eos%db_const(phase)%db(n2)%p_iboundary(p_subn2+1)-1)
        m3 = min0(p_index+1,eos%db_const(phase)%db(n2)%p_iboundary(p_subn2+1))
      end if

      y(4) = eos%db_const(phase)%db(n2)%p(m4)
      y(3) = eos%db_const(phase)%db(n2)%p(m3)


      ! construct regular grid
      if (y(1) .lt. y(4)) then
        m4 = int((y(1)-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
      else
        m1 = int((y(4)-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
      end if

      if (y(2) .gt. y(3)) then
        m3 = int((y(2)-eos%db_const(phase)%db(n2)%p_boundary(p_subn2))/eos%db_const(phase)%db(n2)%p_delta(p_subn2)) + eos%db_const(phase)%db(n2)%p_iboundary(p_subn2)
      else
        m2 = int((y(3)-eos%db_const(phase)%db(n1)%p_boundary(p_subn1))/eos%db_const(phase)%db(n1)%p_delta(p_subn1)) + eos%db_const(phase)%db(n1)%p_iboundary(p_subn1)
      end if

      y(1) = eos%db_const(phase)%db(n1)%p(m1)
      y(2) = eos%db_const(phase)%db(n2)%p(m3)


      dt = x(2)-x(1)
      dp = y(2)-y(1)

      xi  = dmin1(dmax1((t-x(1))/(x(2)-x(1)),0.d0),1.d0)
      eta = dmin1(dmax1((p-y(1))/(y(2)-y(1)),0.d0),1.d0)

      ! bilinear interpolation
      coeff(1:4) = eos%db_const(phase)%db(n1)%vis(m1,:)
      call bilinear_interpolation(coeff(1:4),xi,eta,dt,dp,output)
      prop2%vis = output(1)

      coeff(1:4) = eos%db_const(phase)%db(n1)%cond(m1,:)
      call bilinear_interpolation(coeff(1:4),xi,eta,dt,dp,output)
      prop2%cond = output(1)

    end subroutine srk_g_prop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2o_pww(eos,t) result(pww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
      real(8) :: tt,theta,a,b,c

      tt = dmin1(dmax1(t,273.15d0),647.096d0)
      theta = tt - 0.23855557567849d0/(tt - 0.65017534844798d3)
      a = theta**2 + 0.11670521452767d4*theta - 0.72421316703206d6
      b = - 0.17073846940092d2*theta**2 + 0.12020824702470d5*theta - 0.32325550322333d7
      c = 0.14915108613530d2*theta**2 - 0.48232657361591d4*theta + 0.40511340542057d6
      pww = 1.d6*(2.d0*c/(-b+dsqrt(b**2-4.d0*a*c)))**4
    end function h2o_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function n2_pww(eos,t) result(pww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
      real(8) :: theta

      theta = 1.d0-t/126.192d0
      pww = dexp( (126.192d0/t) * (-6.12445284d0*theta+1.2632722d0*theta**1.5d0               &
                   -0.765910082d0*theta**2.5d0-1.77570564d0*theta**5.d0) ) * 3395800.d0

    end function n2_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function o2_pww(eos,t) result(pww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
      real(8) :: theta

      theta = 1.d0-t/154.581d0
      pww = dexp ( (154.581d0/t) * (-6.043938d0*theta+1.175627d0*theta**1.5d0                 &
                   -0.994086d0*theta**3.d0-3.456781d0*theta**7.d0+3.361499d0*theta**9.d0) )   &
            * 5043000.d0

    end function o2_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2_pww(eos,t) result(pww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
      real(8) :: theta

      theta = 1.d0-t/33.145d0
      pww = dexp( (33.145d0/t) * (-4.89789d0*theta+0.988558d0*theta**1.5d0                    &
                   +0.349689d0*theta**2.d0 + 0.499356d0*theta**2.85d0) ) * 1296400.d0

    end function h2_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function DME_pww(eos,t) result(pww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww

      if(t.le.239d0) then
        pww = 2.51936448014985d-26 * t**1.27973987394184d1
      else if((t.gt.240d0).and.(t.le.279d0)) then
        pww = 2.75551246362286d-20 * t**1.02562411384623d1
      else if((t.gt.280d0).and.(t.le.319d0)) then
        pww = 3.40630759109831d-16 * t**8.58195169665804d0
      else if((t.gt.320d0).and.(t.le.359d0)) then
        pww = 1.83248328175938d-13 * t**7.49068550945385d0
      else
        pww = 5.89759662712948d-12 * t**6.89997409467499d0
      end if
    end function DME_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2o_tww(eos,p) result(tww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p
      real(8) :: tww
      real(8) :: beta,e,f,g,d

      beta = (p/1.d6)**0.25d0
      e = beta**2 - 0.17073846940092d2*beta + 0.14915108613530d2
      f = 0.11670521452767d4*beta**2 + 0.12020824702470d5*beta - 0.48232657361591d4
      g = - 0.72421316703206d6*beta**2 - 0.32325550322333d7*beta + 0.40511340542057d6
      d = 2.d0*g/(-f-dsqrt(f**2-4.d0*e*g))
      tww = 0.5d0*(0.65017534844798d3 + d - dsqrt((0.65017534844798d3 + d)**2                 &
                   - 4.d0*(-0.23855557567849d0+0.65017534844798d3*d)))
    end function h2o_tww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function n2_tww(eos,p) result(tww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p
      real(8) :: tww

      if(p.le.100000.d0) then
        tww = 3.186d0*p**0.217d0 + 38.5d0
      else if((p.gt.100000.d0).and.(p.le.300000.d0)) then
        tww = 1.899d0*p**0.2499d0 + 43.51d0
      else if((p.gt.300000.d0).and.(p.le.1000000.d0)) then
        tww = 1.245d0*p**0.2302d0 + 47.67d0
      else
        tww = 1.591*p**0.2622d0 + 44.17d0
      end if
    end function n2_tww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function o2_tww(eos,p) result(tww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p
      real(8) :: tww

      if(p.lt.1.d4) then
        tww = 9.153d0*p**0.1546d0 + 34.65d0
      else if((p.ge.1.d4).and.(p.lt.1.d5)) then
        tww = 4.41d0*p**0.2042d0 + 43.79d0
      else if((p.ge.1.d5).and.(p.lt.1.d6)) then
        tww = 1.899d0*p**0.2569d0 + 53.53d0
      else
        tww = 1.51d0*p**0.2709d0 + 55.87d0
      end if
    end function o2_tww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2_tww(eos,p) result(tww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p
      real(8) :: tww

      if(p.le.1.d5) then
        tww = 0.682d0*p**0.2562d0 + 7.289d0
      else if((p.gt.1.d5).and.(p.le.3.d5)) then
        tww = 0.4102d0*p**0.2903d0 + 8.724d0
      else if((p.gt.3.d5).and.(p.le.7.d5)) then
        tww = 0.366d0*p**0.2975d0 + 9.091d0
      else if((p.gt.7.d5).and.(p.le.1.d6)) then
        tww = 0.5396d0*p**0.2746d0 + 7.4216d0
      else
        tww = 3.377d0*p**0.1748d0 - 6.394d0
      end if
    end function h2_tww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function DME_tww(eos,p) result(tww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p
      real(8) :: tww

      if(p.le.1.d5) then
        tww = 3.32986432802911d1 + 1.85534379486282d1*LOG(p)
      else if((p.gt.1.d5).and.(p.le.3.d5)) then
        tww = -5.57196453839085d1 + 2.63140232108003d1*LOG(p)
      else if((p.gt.3.d5).and.(p.le.7.d5)) then
        tww = -1.38577037376019d2 + 3.28692957589036d1*LOG(p)
      else if((p.gt.7.d5).and.(p.le.1.d6)) then
        tww = -2.07517520436400d2 + 3.80114376290005d1*LOG(p)
      else if((p.gt.1.d6).and.(p.le.3.d6)) then
        tww = -3.25691065803963d2 + 4.64639077454879d1*LOG(p)
      else
        tww = -4.63045141170944d2 + 5.57207292900318d1*LOG(p)
      end if
    end function DME_tww
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2o_sigma(eos,t) result(sigma) ! surface tension
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: sigma

      sigma = 0.0881130429998226d0 + 0.0000616475914044782d0*t - 4.4860167748009d-7*t**2      &
            + 2.53682468718936d-10*t**3 - 2.62005668908913d-13*t**4 + 2.97878955533757d-16*t**5

    end function h2o_sigma
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function n2_sigma(eos,t) result(sigma)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: sigma

      sigma = 0.0277139593881926d0 - 0.0001871381651395d0*t - 2.47277506113986d-6*t**2        &
            + 3.73212791877442d-8*t**3 - 2.44521143856208d-10*t**4 + 6.90956437271856d-13*t**5

    end function n2_sigma
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function o2_sigma(eos,t) result(sigma)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: sigma

      sigma = 0.0383957397934739d0 - 0.000292004596094666d0*t - 1.63784793251168d-7*t**2      &
            + 5.57366451318535d-9*t**3 - 3.42769625844841d-11*t**4 + 1.05411729447465d-13*t**5

    end function o2_sigma
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2_sigma(eos,t) result(sigma)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: sigma

      sigma = 0.00520006147168522d0 - 0.000129664670257188d0*t -  4.14551445475162d-6*t**2    &
            + 2.18364493308481d-7*t**3 -  5.43761039654375d-9*t**4 + 5.62476539114984d-11*t**5

    end function h2_sigma
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function DME_sigma(eos,t) result(sigma)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: sigma

      if(t.le.200.d0) then
        sigma = 7.674692739168850d12 * t**(-4.460029937536576d0)
      else if((t.gt.200.d0).and.(t.le.300.d0)) then
        sigma = 3.360744776376887d9 * t**(-3.004965901433068d0)
      else
        sigma = 7.719551553645169d13 * t**(-4.726125971944006d0)
      end if

    end function DME_sigma
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function cubicsolver(pin,tin,aalpha,bin,mw)
      implicit none
      real(8), intent(in)  :: pin,tin,aalpha,bin,mw
      real(8) :: cubicsolver
      real(8) :: p,q,r
      real(8) :: a,b,d,m,n
      real(8) :: phi,k,pi
      real(8),parameter :: ru=8.314d0
      pi = datan(1.d0)*4.d0

      p = (bin**2*pin*mw - aalpha*mw + ru*tin*bin*mw)/(aalpha*bin)
      q = mw**2*ru*tin/(aalpha*bin)
      r = -mw**3*pin/(aalpha*bin)

      a = (3.d0*q - p**2)/3.d0
      b = (2.d0*p**3 - 9.d0*p*q + 27.d0*r)/27.d0
      d = a**3/27.d0 + b**2/4.d0

      if (d.gt.0.d0) then
        m = (-0.5d0*b+dsqrt(d))
        n = (-0.5d0*b-dsqrt(d))
        m = dsign(1.d0,m)*dabs(m)**(1.d0/3.d0)
        n = dsign(1.d0,n)*dabs(n)**(1.d0/3.d0)
        cubicsolver = m+n-p/3.d0
      else if(d.eq.0.d0) then
        m = (-0.5d0*b+dsqrt(d))
        n = (-0.5d0*b-dsqrt(d))
        m = dsign(1.d0,m)*dabs(m)**(1.d0/3.d0)
        n = dsign(1.d0,n)*dabs(n)**(1.d0/3.d0)
        if((m+n - p/3.d0).ge.0.d0) then
          cubicsolver = m+n - p/3.d0
        else
          cubicsolver = -0.5d0*(m+n) - p/3.d0
        end if
      else if (d.lt.0.d0) then
        k = 0.d0
        if (b.gt.0.d0) then
          phi = dacos(-dsqrt(-27.d0*b**2/(4.d0*a**3)))
        else
          phi = dacos(dsqrt(-27.d0*b**2/(4.d0*a**3)))
        end if
        cubicsolver = 2.d0*dsqrt(-a/3.d0)*dcos((phi+120.d0*k)*pi/180.d0) - p/3.d0
      end if

    end function cubicsolver
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function computepv(tc,pc,tb,t)
      implicit none
      real(8), intent(in) :: tc,pc,tb,t
      real(8) :: computepv

      ! using reidel's formula

      real(8) :: a,b,c,d
      real(8) :: tr,tbr,pbr
      real(8) :: k1,k2,psi,q,alphac,tt

      k1 = 0.0838d0
      k2 = 3.758d0

      tt=dmin1(tc-2.d0,t)

      tr  = tt/tc
      tbr = tb/tc
      pbr = 101325.d0/pc

      psi = -35.d0 + 36.d0/tbr + 42.d0*dlog(tbr) - tbr**6
      alphac = (k1*k2*psi-dlog(pbr))/(k1*psi-dlog(tbr))
      q = k1*(k2-alphac)

      a = -35.d0*q
      b = -36.d0*q
      c = 42.d0*q + alphac
      d = -q
      computepv = pc*dexp(a - b/tr + c*dlog(tr) + d*tr**6)

    end function computepv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getgamma_f(eos)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8) :: getgamma_f

      getgamma_f = eos%gamma_f

    end function getgamma_f
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getgamma_g(eos)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8) :: getgamma_g

      getgamma_g = eos%gamma_g

    end function getgamma_g
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine biquadratic_interpolation(coeff,xi,eta,dt,dp,output)
      real(8), intent(in) :: coeff(9)
      real(8), intent(in) :: xi,eta,dp,dt
      real(8), intent(inout) :: output(3)

      output(1) = coeff(1)        + coeff(2)*xi        + coeff(3)*xi**2 +        &
                  coeff(4)*eta    + coeff(5)*xi*eta    + coeff(6)*xi**2*eta +    &
                  coeff(7)*eta**2 + coeff(8)*xi*eta**2 + coeff(9)*xi**2*eta**2

      output(2) = coeff(2) + 2.d0*coeff(3)*xi + coeff(5)*eta + 2.d0*coeff(6)*xi*eta + coeff(8)*eta**2 + 2.d0*coeff(9)*xi*eta**2

      output(3) = coeff(4) + 2.d0*coeff(7)*eta + coeff(5)*xi  + 2.d0*coeff(8)*xi*eta + coeff(6)*xi**2  + 2.d0*coeff(9)*xi**2*eta

      output(2) = output(2)/dt
      output(3) = output(3)/dp

    end subroutine biquadratic_interpolation
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bilinear_interpolation(coeff,xi,eta,dt,dp,output)
      real(8), intent(in) :: coeff(4)
      real(8), intent(in) :: xi,eta,dt,dp
      real(8), intent(inout) :: output(3)

      output(1) = coeff(1) + coeff(2)*xi + coeff(3)*eta + coeff(4)*xi*eta
      output(2) = (coeff(2) + coeff(4)*eta)/dt
      output(3) = (coeff(3) + coeff(4)*xi )/dp

    end subroutine bilinear_interpolation
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module eos_module
