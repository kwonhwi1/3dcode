module eos_module
  use mpi
  implicit none
  private
  public :: t_eos
    
  type t_iapws_coeff_eos
    integer, dimension(:), allocatable :: i,j,ir,jr
    real(8), dimension(:), allocatable :: n,nr
  end type t_iapws_coeff_eos
  
  type t_iapws_coeff_prop
    real(8), dimension(:), allocatable :: h,l
    real(8), dimension(:,:), allocatable :: hh,ll
  end type t_iapws_coeff_prop

  type t_srk_property
    real(8) :: tc, pc, tb, mw, w, k, cv0
  end type t_srk_property
  
  type t_db_data
    integer :: p_ndata
    real(8) :: t
    real(8),allocatable :: dbv(:,:) !p,rho,h,drdp,drdt,dhdp,dhdt,vis,cond
  end type t_db_data

  type t_db_const
    integer :: tndata_db
    real(8) :: delta_t_db,delta_p_db
    type(t_db_data),allocatable :: db(:)
  end type t_db_const

  type t_eos
    private
    integer :: ntv,rank,size
    procedure(p_pww), public, pointer :: get_pww
    procedure(p_sigma), public, pointer :: get_sigma
    procedure(p_eos_l), pointer :: eos_l
    procedure(p_eos_v), pointer :: eos_v
    procedure(p_eos_g), pointer :: eos_g
    type(t_iapws_coeff_eos) :: iapws_liquid_coeff,iapws_vapor_coeff,iapws_m_vapor_coeff
    type(t_iapws_coeff_prop) :: iapws_prop_coeff
    type(t_srk_property) :: srk_gas_property
    type(t_db_const) :: db_const(3) ! 1=liquid,2=vapor,3=gas
    contains
      procedure :: construct !(fluid,fluid_eostype,ngas,gas_eostype,size,rank)
      procedure :: deteos
      procedure, private :: set_iapws97_eos
      procedure, private :: set_iapws97_prop
      procedure, private :: set_srk_property
      procedure, private :: set_database
      procedure :: destruct
  end type t_eos

  type t_eos2
    real(8) :: rho,h,drdp,drdt,dhdp,dhdt
  end type t_eos2

  type t_prop2
    real(8) :: vis,cond
  end type t_prop2

  type t_srk_coeff
    real(8) :: a,alpha,b,s
  end type t_srk_coeff
  
  interface
    function p_pww(eos,t) result(pww)
      import t_eos
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
    end function p_pww
    
    function p_sigma(eos,t) result(sigma)
      import t_eos
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: sigma
    end function p_sigma
    
    subroutine p_eos_l(eos,p,t,phase,eos2,prop2)
      import t_eos
      import t_eos2
      import t_prop2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
    end subroutine p_eos_l
    
    subroutine p_eos_v(eos,p,t,phase,eos2,prop2)
      import t_eos
      import t_eos2
      import t_prop2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
    end subroutine p_eos_v
    
    subroutine p_eos_g(eos,p,t,phase,eos2,prop2)
      import t_eos
      import t_eos2
      import t_prop2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
    end subroutine p_eos_g
  end interface

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine construct(eos,fluid,fluid_eostype,ngas,gas_eostype,ntv,size,rank)
      implicit none
      class(t_eos), intent(out) :: eos
      integer, intent(in) :: fluid,fluid_eostype,ngas,gas_eostype,ntv,size,rank

      eos%ntv = ntv
      eos%size  = size
      eos%rank  = rank

      select case(ngas)
      case(1) ! ideal gas
        select case(gas_eostype)
        case(1) ! ideal gas law
          eos%eos_g  => ideal
        case default
        end select
      case(2,3,4,5) ! nitrogen,oxygen,hydrogen,helium
        select case(gas_eostype)
        case(1) ! database
          call eos%set_database(ngas,3)
          eos%eos_g  => database
        case(2) ! srk
          call eos%set_srk_property(ngas)
          if(ntv.ge.2) call eos%set_database(ngas,3)
          eos%eos_g  => srk_g
        case default
        end select
      case default
      end select

      select case(fluid_eostype)
      case(1) ! database
        call eos%set_database(fluid,1)
        eos%eos_l  => database
        call eos%set_database(fluid,2)
        eos%eos_v  => database
      case(2) ! iapws97 for water
        call eos%set_iapws97_eos()
        if(ntv.ge.2) call eos%set_iapws97_prop()
        eos%eos_l  => iapws97_l
        eos%eos_v  => iapws97_v
      case(3) ! stiffened for water
        eos%eos_l  => stiffened
        eos%eos_v  => ideal
      case default
      end select

      select case(fluid)
      case(1) ! water
        eos%get_pww => h2o_pww
        eos%get_sigma => h2o_sigma
      case(2) ! nitrogen
        eos%get_pww => n2_pww
        eos%get_sigma => n2_sigma
      case(3) ! oxygen
        eos%get_pww => o2_pww
        eos%get_sigma => o2_sigma
      case(4) ! hydrogen
        eos%get_pww => h2_pww
        eos%get_sigma => h2_sigma
      case default
      end select
    
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine destruct(eos)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer :: m,n

      if(associated(eos%get_pww))   nullify(eos%get_pww)
      if(associated(eos%get_sigma)) nullify(eos%get_sigma)
      if(associated(eos%eos_l))     nullify(eos%eos_l)
      if(associated(eos%eos_v))     nullify(eos%eos_v)
      if(associated(eos%eos_g))     nullify(eos%eos_g)

      do n=1,3
        do m=1,eos%db_const(n)%tndata_db
          if(allocated(eos%db_const(n)%db(m)%dbv)) deallocate(eos%db_const(n)%db(m)%dbv)
        end do
        if(allocated(eos%db_const(n)%db)) deallocate(eos%db_const(n)%db)
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
    subroutine deteos(eos,p,t,y1,y2,dv,tv)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(18) ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      real(8), intent(inout) :: tv(eos%ntv)
      type(t_eos2)  :: l_eos,v_eos,g_eos
      type(t_prop2) :: l_prop,v_prop,g_prop

      call eos%eos_l(p,t,1,l_eos,l_prop)
      call eos%eos_v(p,t,2,v_eos,v_prop)
      call eos%eos_g(p,t,3,g_eos,g_prop)
      
      dv(3)  = l_eos%rho
      dv(4)  = v_eos%rho
      dv(5)  = g_eos%rho

      dv(1)  = 1.d0/((1.d0-y1-y2)/l_eos%rho + y1/v_eos%rho + y2/g_eos%rho)
      dv(2)  = (1.d0-y1-y2)*l_eos%h + y1*v_eos%h + y2*g_eos%h

      dv(7)  = dv(1)**2*((1.d0-y1-y2)/l_eos%rho**2*l_eos%drdp &
                                 + y1/v_eos%rho**2*v_eos%drdp   &
                                 + y2/g_eos%rho**2*g_eos%drdp)
      dv(8)  = dv(1)**2*((1.d0-y1-y2)/l_eos%rho**2*l_eos%drdt &
                                 + y1/v_eos%rho**2*v_eos%drdt   &
                                 + y2/g_eos%rho**2*g_eos%drdt)
      dv(9)  = dv(1)**2*(1.d0/l_eos%rho - 1.d0/v_eos%rho)
      dv(10) = dv(1)**2*(1.d0/l_eos%rho - 1.d0/g_eos%rho)
    
      dv(11) = (1.d0-y1-y2)*l_eos%dhdp + y1*v_eos%dhdp + y2*g_eos%dhdp
      dv(12) = (1.d0-y1-y2)*l_eos%dhdt + y1*v_eos%dhdt + y2*g_eos%dhdt
      dv(13) = v_eos%h - l_eos%h
      dv(14) = g_eos%h - l_eos%h
      
      dv(6)  = dv(1)*dv(12)/(dv(1)*dv(7)*dv(12)+dv(8)*(1.d0-dv(1)*dv(11)))
      
      dv(15) = v_eos%drdp
      dv(16) = v_eos%drdt
      dv(17) = l_eos%drdp
      dv(18) = l_eos%drdt

      if(eos%ntv.ge.2) then
        tv(1) = dv(1)*(l_prop%vis*(1.d0-y1-y2)/dv(3)  + v_prop%vis*y1/dv(4)  &
                       + g_prop%vis*y2/dv(5))
        tv(2) = dv(1)*(l_prop%cond*(1.d0-y1-y2)/dv(3) + v_prop%cond*y1/dv(4) &
                       + g_prop%cond*y2/dv(5))
      end if

    end subroutine deteos
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine ideal(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      real(8), parameter :: r1 = 1.4d0/0.4d0/1004.64
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
      if(eos%ntv.ge.2) then
#ifdef shocktube
      prop2%vis = 0.006d0
#else
      prop2%vis  = 1.716d-5*(t/273.11d0)**1.5d0*383.67d0/(t+110.56d0)
      !prop2%cond = 0.0241d0*(t/273.d0)**1.5d0*467.d0/(t+194.d0)
#endif
      prop2%cond = 1004.64d0/0.72d0*prop2%vis
      end if

    end subroutine ideal
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_database(eos,fluid,phase)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer, intent(in) :: fluid,phase
      integer :: n,m,l,k,ierr

      do k=0,eos%size-1
        if(k.eq.eos%rank) then

          select case(phase)
          case(1) ! liquid
            select case(fluid)
            case(1) ! water
              open(unit = 10, file = "./../database/DATABASE_WATER_LIQ.DAT", FORM = "BINARY")
            case(2) ! nitrogen
              open(unit = 10, file = "./../database/DATABASE_NITROGEN_LIQ.DAT", FORM = "BINARY")
            case(3) ! oxygen
              open(unit = 10, file = "./../database/DATABASE_OXYGEN_LIQ.DAT", FORM = "BINARY")
            case(4) ! hydrogen
              open(unit = 10, file = "./../database/DATABASE_HYDROGEN_LIQ.DAT", FORM = "BINARY")
            end select

          case(2,3) ! vapor, gas
            select case(fluid)
            case(1) ! water
              open(unit = 10, file = "./../database/DATABASE_WATER_GAS.DAT", FORM = "BINARY")
            case(2) ! nitrogen
              open(unit = 10, file = "./../database/DATABASE_NITROGEN_GAS.DAT", FORM = "BINARY")
            case(3) ! oxygen
              open(unit = 10, file = "./../database/DATABASE_OXYGEN_GAS.DAT", FORM = "BINARY")
            case(4) ! hydrogen
              open(unit = 10, file = "./../database/DATABASE_HYDROGEN_GAS.DAT", FORM = "BINARY")
            case(5) ! helium
              open(unit = 10, file = "./../database/DATABASE_HELIUM_GAS.DAT", FORM = "BINARY")
            end select
          end select

          read(10) eos%db_const(phase)%tndata_db
          read(10) eos%db_const(phase)%delta_t_db
          read(10) eos%db_const(phase)%delta_p_db

          allocate(eos%db_const(phase)%db(eos%db_const(phase)%tndata_db))

          do m=1,eos%db_const(phase)%tndata_db
            read(10) eos%db_const(phase)%db(m)%t
            read(10) eos%db_const(phase)%db(m)%p_ndata

            allocate(eos%db_const(phase)%db(m)%dbv(eos%db_const(phase)%db(m)%p_ndata,9))

            do n=1,eos%db_const(phase)%db(m)%p_ndata
              do l=1,9
                read(10) eos%db_const(phase)%db(m)%dbv(n,l)
              end do
            end do
          end do

          close(10)

        end if
      end do

    end subroutine set_database
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine database(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      integer :: t_index,p_index
      integer :: n1,n2,m1,m2
      real(8) :: x(2),y(2)
      real(8) :: a0,a1,a2,a3,a4

      ! defining the temperature index
      t_index = min0(max0(int((t-eos%db_const(phase)%db(1)%t)/eos%db_const(phase)%delta_t_db)+1,1),eos%db_const(phase)%tndata_db)
      !t_index = min0(max0(nint((t-eos%db_const(phase)%db(1)%t)/eos%db_const(phase)%delta_t_db+1),1),eos%db_const(phase)%tndata_db)
      if(t.ge.eos%db_const(phase)%db(t_index)%t) then
        n1 = min0(max0(t_index,1),eos%db_const(phase)%tndata_db-1)
        n2 = min0(max0(t_index+1,2),eos%db_const(phase)%tndata_db)
      else
        n1 = max0(t_index-1,1)
        n2 = max0(t_index,2)
      end if
      
      x(1) = eos%db_const(phase)%db(n1)%t
      x(2) = eos%db_const(phase)%db(n2)%t
      
      ! defining the pressure index
      
      p_index = min0(max0(int(p/eos%db_const(phase)%delta_p_db)+1,1),eos%db_const(phase)%db(n1)%p_ndata)
      !p_index = nint(p/eos%db_const(phase)%delta_p_db)
      if(p.ge.eos%db_const(phase)%db(n1)%dbv(p_index,1)) then
        m1 = min0(max0(p_index,1),eos%db_const(phase)%db(n1)%p_ndata-1)
        m2 = min0(max0(p_index+1,2),eos%db_const(phase)%db(n1)%p_ndata)
      else
        m1 = max0(p_index-1,1)
        m2 = max0(p_index,2)
      end if
      
      y(1) = eos%db_const(phase)%db(n1)%dbv(m1,1)
      y(2) = eos%db_const(phase)%db(n2)%dbv(m2,1)
      
      a0 = 1.d0/((x(2)-x(1))*(y(2)-y(1)))
      a1 = dabs((x(2)-t)*(y(2)-p))
      a2 = dabs((t-x(1))*(y(2)-p))
      a3 = dabs((x(2)-t)*(p-y(1)))
      a4 = dabs((t-x(1))*(p-y(1)))


      eos2%rho  = (eos%db_const(phase)%db(n1)%dbv(m1,2)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,2)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,2)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,2)*a4 )*a0
      eos2%h    = (eos%db_const(phase)%db(n1)%dbv(m1,3)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,3)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,3)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,3)*a4 )*a0
      eos2%drdp = (eos%db_const(phase)%db(n1)%dbv(m1,4)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,4)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,4)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,4)*a4 )*a0
      eos2%drdt = (eos%db_const(phase)%db(n1)%dbv(m1,5)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,5)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,5)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,5)*a4 )*a0
      eos2%dhdp = (eos%db_const(phase)%db(n1)%dbv(m1,6)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,6)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,6)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,6)*a4 )*a0
      eos2%dhdt = (eos%db_const(phase)%db(n1)%dbv(m1,7)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,7)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,7)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,7)*a4 )*a0

      if(eos%ntv.ge.2) then
      prop2%vis  = (eos%db_const(phase)%db(n1)%dbv(m1,8)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,8)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,8)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,8)*a4 )*a0
      prop2%cond = (eos%db_const(phase)%db(n1)%dbv(m1,9)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,9)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,9)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,9)*a4 )*a0
      end if

    end subroutine database
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_srk_property(eos,fluid)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer, intent(in) :: fluid
      type(t_srk_property) :: srk_property

      select case(fluid)
      case(2) ! nitrogen
        srk_property%tc  = 126.19d0
        srk_property%pc  = 3.3958d6
        srk_property%tb  = 77.355d0
        srk_property%mw  = 28.013d-3
        srk_property%w   = 0.0372d0
        srk_property%cv0 = 742.2d0
        
      case(3) ! oxyten
        srk_property%tc  = 154.58d0
        srk_property%pc  = 5.043d6
        srk_property%tb  = 90.188d0
        srk_property%mw  = 31.999d-3
        srk_property%w   = 0.0222d0
        srk_property%cv0 = 650.d0

      case(4) ! hydrogen
        srk_property%tc  = 33.145d0
        srk_property%pc  = 1.2964d6
        srk_property%tb  = 20.369d0
        srk_property%mw  = 2.0159d-3
        srk_property%w   = -0.219d0
        srk_property%cv0 = 0.d0         ! should be filled
        
      case(5) ! helium
        srk_property%tc  = 5.1953d0
        srk_property%pc  = 0.22746d6
        srk_property%tb  = 4.23d0
        srk_property%mw  = 4.0026d-3
        srk_property%w   = -0.382d0
        srk_property%cv0 = 0.d0         ! should be filled
      case default
      end select

      eos%srk_gas_property = srk_property

    end subroutine set_srk_property
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_srk_coeff(p,t,srk_property,srk_coeff)
      implicit none
      real(8), intent(in) :: p,t
      type(t_srk_property), intent(in) :: srk_property
      type(t_srk_coeff), intent(out) :: srk_coeff
      real(8),parameter :: ru=8.314d0

      srk_coeff%s     = 0.48508d0 + 1.5517d0*srk_property%w - 0.15613d0*srk_property%w**2
      srk_coeff%a     = 0.42747d0*(ru*srk_property%tc)**2/srk_property%pc
      srk_coeff%alpha = (1.d0 + srk_coeff%s*(1.d0-dsqrt(t/srk_property%tc)))**2
      srk_coeff%b     = 0.08664d0*ru*srk_property%tc/srk_property%pc

    end subroutine set_srk_coeff
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine srk_g(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out) :: prop2
      type(t_srk_coeff) :: srk_coeff
      real(8),parameter :: ru=8.314d0
      real(8) :: a,alpha,b,mw,s,tc,pc,cv0
      real(8) :: rho,e,h,pv
      real(8) :: dpdt,dpdr
      real(8) :: daadt,dsaadt,daidt
      real(8) :: ddaaddt,ddsaaddt,ddaiddt
      real(8) :: cv,cp,bp,bt,ap,ar
      integer :: t_index,p_index
      integer :: n1,n2,m1,m2
      real(8) :: x(2),y(2)
      real(8) :: a0,a1,a2,a3,a4

      call set_srk_coeff(p,t,eos%srk_gas_property,srk_coeff)

      a     = srk_coeff%a
      alpha = srk_coeff%alpha
      b     = srk_coeff%b
      mw    = eos%srk_gas_property%mw
      s     = srk_coeff%s
      tc    = eos%srk_gas_property%tc
      pc    = eos%srk_gas_property%pc
      cv0   = eos%srk_gas_property%cv0

      pv = computepv(tc,pc,eos%srk_gas_property%tb,t)
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

      if(eos%ntv.ge.2) then
      ! defining the temperature index
      t_index = min0(max0(int((t-eos%db_const(phase)%db(1)%t)/eos%db_const(phase)%delta_t_db)+1,1),eos%db_const(phase)%tndata_db)
      !t_index = min0(max0(nint((t-eos%db_const(phase)%db(1)%t)/eos%db_const(phase)%delta_t_db+1),1),eos%db_const(phase)%tndata_db)
      if(t.ge.eos%db_const(phase)%db(t_index)%t) then
        n1 = min0(max0(t_index,1),eos%db_const(phase)%tndata_db-1)
        n2 = min0(max0(t_index+1,2),eos%db_const(phase)%tndata_db)
      else
        n1 = max0(t_index-1,1)
        n2 = max0(t_index,2)
      end if

      x(1) = eos%db_const(phase)%db(n1)%t
      x(2) = eos%db_const(phase)%db(n2)%t

      ! defining the pressure index
      p_index = min0(max0(int(p/eos%db_const(phase)%delta_p_db)+1,1),eos%db_const(phase)%db(n1)%p_ndata)
      !p_index = nint(p/eos%db_const(phase)%delta_p_db)
      if(p.ge.eos%db_const(phase)%db(n1)%dbv(p_index,1)) then
        m1 = min0(max0(p_index,1),eos%db_const(phase)%db(n1)%p_ndata-1)
        m2 = min0(max0(p_index+1,2),eos%db_const(phase)%db(n1)%p_ndata)
      else
        m1 = max0(p_index-1,1)
        m2 = max0(p_index,2)
      end if

      y(1) = eos%db_const(phase)%db(n1)%dbv(m1,1)
      y(2) = eos%db_const(phase)%db(n2)%dbv(m2,1)

      a0 = 1.d0/((x(2)-x(1))*(y(2)-y(1)))
      a1 = dabs((x(2)-t)*(y(2)-p))
      a2 = dabs((t-x(1))*(y(2)-p))
      a3 = dabs((x(2)-t)*(p-y(1)))
      a4 = dabs((t-x(1))*(p-y(1)))

      prop2%vis  = (eos%db_const(phase)%db(n1)%dbv(m1,8)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,8)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,8)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,8)*a4 )*a0
      prop2%cond = (eos%db_const(phase)%db(n1)%dbv(m1,9)*a1 + eos%db_const(phase)%db(n2)%dbv(m1,9)*a2 + eos%db_const(phase)%db(n1)%dbv(m2,9)*a3 + eos%db_const(phase)%db(n2)%dbv(m2,9)*a4 )*a0
      end if

    end subroutine srk_g
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
    subroutine iapws97_l(eos,p,t,phase,eos2,prop2)
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

      pi = p/16.53d6
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

      if(eos%ntv.ge.2) then
      tt = t/647.27d0
      rr = dmax1(1.d-3,dmin1(1332.4d0,eos2%rho))/317.763d0
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
      end if
    end subroutine iapws97_l
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97_v(eos,p,t,phase,eos2,prop2)
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
      
      pi = p/1d6
      tau = 540d0/t
      
      pi1 = 1.d0/pi
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
      
      if(eos%ntv.ge.2) then
      tt = t/647.27d0
      rr = dmax1(1.d-3,dmin1(1332.4d0,eos2%rho))/317.763d0
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
      end if
    end subroutine iapws97_v
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine stiffened(eos,p,t,phase,eos2,prop2)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      integer, intent(in) :: phase
      type(t_eos2), intent(out) :: eos2
      type(t_prop2), intent(out):: prop2
      real(8), parameter :: r1 = 1.d0/2691.d0
      real(8) :: tau,temp

      tau = 1.d0/t

      eos2%rho  = (p+8.5d8)*r1*tau
      eos2%h    = 4186.d0*t-917822.36d0
      eos2%drdp = r1*tau
      eos2%drdt = -(p+8.5d8)*r1*tau**2
      eos2%dhdp = 0.d0
      eos2%dhdt = 4186.d0

      if(eos%ntv.ge.2) then
      temp = dmax1(273.d0,t)

      prop2%vis  = 1.788d-3*dexp(-1.704d0 - 5.306d0*273.d0/temp + 7.003d0*(273.d0/temp)**2)
      prop2%cond = -0.000009438d0*t**2 + 0.007294d0*t - 0.7284d0
      end if
    end subroutine stiffened
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

      pww = 0.03326d0*t**4 + 0.1123d0*t**3 - 680.d0*t**2 + 53550.d0*t - 1215000.d0

    end function n2_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function o2_pww(eos,t) result(pww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
      pww = 0.d0
    end function o2_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2_pww(eos,t) result(pww)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: t
      real(8) :: pww
      pww = 0.d0
    end function h2_pww
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function h2o_sigma(eos,t) result(sigma)
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
        m = (-0.5d0*b+dsqrt(d))**(1.d0/3.d0)
        n = (-0.5d0*b-dsqrt(d))**(1.d0/3.d0)
        cubicsolver = m+n-p/3.d0  
      else if(d.eq.0.d0) then
        m = (-0.5d0*b+dsqrt(d))**(1.d0/3.d0)
        n = (-0.5d0*b-dsqrt(d))**(1.d0/3.d0)
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
end module eos_module
