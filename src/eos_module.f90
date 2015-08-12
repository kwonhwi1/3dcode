module eos_module
  implicit none
  private
  public :: t_eos
    
  type t_iapws_coeff
    integer, dimension(:), allocatable :: i,j,ir,jr
    real(8), dimension(:), allocatable :: n,nr
  end type t_iapws_coeff
  
  type t_srk_property
    real(8) :: tc, pc, tb, mw, w, k, cv0
  end type t_srk_property
  
  type t_eos
    private
    procedure(p_pww), public, pointer :: get_pww
    procedure(p_sigma), public, pointer :: get_sigma
    procedure(p_deteos), public, pointer :: deteos
    procedure(p_eos_l), pointer :: eos_l
    procedure(p_eos_v), pointer :: eos_v
    procedure(p_eos_g), pointer :: eos_g
    type(t_iapws_coeff) :: iapws_liquid_coeff,iapws_vapor_coeff,iapws_m_vapor_coeff
    type(t_srk_property) :: srk_vapor_property, srk_gas_property
    contains
      procedure :: construct !(fluid,eos_type,n_gas)
      procedure, private :: set_iapws97
      procedure, private :: set_srk_property
      procedure, private :: set_srk_binary
      procedure, private :: srk_m
      procedure :: destruct
  end type t_eos

  type t_eos2
    real(8) :: rho,h,drdp,drdt,dhdp,dhdt,cv
  end type t_eos2

  type t_srk_coeff
    real(8) :: a, alpha, b, s
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
    
    subroutine p_deteos(eos,p,t,y1,y2,dv)
      import t_eos
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(18)
    end subroutine p_deteos
    
    subroutine p_eos_l(eos,p,t,liquid)
      import t_eos
      import t_eos2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid
    end subroutine p_eos_l
    
    subroutine p_eos_v(eos,p,t,vapor)
      import t_eos
      import t_eos2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: vapor
    end subroutine p_eos_v
    
    subroutine p_eos_g(eos,p,t,gas)
      import t_eos
      import t_eos2
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: gas
    end subroutine p_eos_g
  end interface

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine construct(eos,fluid,fluid_eostype,ngas,gas_eostype,mixingrule)
      implicit none
      class(t_eos), intent(out) :: eos
      integer, intent(in) :: fluid,fluid_eostype,ngas,gas_eostype,mixingrule
      
      select case(mixingrule) ! mixing rule for ullage
        case(1) ! amagat's law
          eos%deteos => deteos_amagat
        case(2) ! dalton's law
          eos%deteos => deteos_dalton
      end select
      
      select case(gas_eostype)
      case(1) ! fitting or ideal gas law  
        select case(ngas)
        case(1) ! ideal gas
          eos%eos_g => ideal
        case(2) ! nitrogen
          eos%eos_g => n2_v_fit
        case(3) ! oxygen
          eos%eos_g => o2_v_fit
        case(4) ! hydrogen
          eos%eos_g => h2_v_fit
        case(5) ! helium
          eos%eos_g => he_v_fit
        case default  
        end select
      case(2) ! srk
        call eos%set_srk_property(ngas,2)
        call eos%set_srk_binary(fluid,ngas)
        eos%eos_g => srk_g
      case default  
      end select
      
      select case(fluid)
      case(1) ! water
        eos%get_pww => h2o_pww
        eos%get_sigma => h2o_sigma
        select case(fluid_eostype)
        case(1) ! fitting
          eos%eos_l  => h2o_l_fit
          eos%eos_v  => h2o_v_fit
        case(2) ! iapws97
          call eos%set_iapws97()
          eos%eos_l  => iapws97_l
          eos%eos_v  => iapws97_v
        case(3) ! stiffened
          eos%eos_l  => stiffened
          eos%eos_v  => ideal
        end select

      case(2) ! nitrogen
        eos%get_pww => n2_pww
        eos%get_sigma => n2_sigma
        eos%eos_l  => n2_l_fit
        select case(fluid_eostype)
        case(1) ! fitting
          eos%eos_v  => n2_v_fit
        case(2) ! srk
          call eos%set_srk_property(fluid,1)
          eos%eos_v => srk_v
        case default
        end select
      
      case(3) ! oxygen
        eos%get_pww => o2_pww
        eos%get_sigma => o2_sigma
        !eos%eos_l  => srk_l
        eos%eos_l => o2_l_fit
        select case(fluid_eostype)
        case(1) ! fitting
          eos%eos_v  => o2_v_fit
        case(2) ! srk
          call eos%set_srk_property(fluid,1)
          eos%eos_v => srk_v
        case default
        end select
        
      case(4) ! hydrogen
        eos%get_pww => h2_pww
        eos%get_sigma => h2_sigma
        eos%eos_l  => h2_l_fit
        select case(fluid_eostype)
        case(1) ! fitting
          eos%eos_v  => h2_v_fit
        case(2) ! srk
          call eos%set_srk_property(fluid,1)
          eos%eos_v => srk_v     
        case default
        end select
      end select
    
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine destruct(eos)
      implicit none
      class(t_eos), intent(inout) :: eos
    
      if(associated(eos%get_pww))   nullify(eos%get_pww)
      if(associated(eos%get_sigma)) nullify(eos%get_sigma)
      if(associated(eos%deteos))    nullify(eos%deteos)
      if(associated(eos%eos_l))     nullify(eos%eos_l)
      if(associated(eos%eos_v))     nullify(eos%eos_v)
      if(associated(eos%eos_g))     nullify(eos%eos_g)
      
      if(allocated(eos%iapws_liquid_coeff%i)) deallocate(eos%iapws_liquid_coeff%i)
      if(allocated(eos%iapws_liquid_coeff%j)) deallocate(eos%iapws_liquid_coeff%j)
      if(allocated(eos%iapws_liquid_coeff%n)) deallocate(eos%iapws_liquid_coeff%n)
      if(allocated(eos%iapws_vapor_coeff%j)) deallocate(eos%iapws_vapor_coeff%j)
      if(allocated(eos%iapws_vapor_coeff%n)) deallocate(eos%iapws_vapor_coeff%n)
      if(allocated(eos%iapws_vapor_coeff%ir)) deallocate(eos%iapws_vapor_coeff%ir)
      if(allocated(eos%iapws_vapor_coeff%jr)) deallocate(eos%iapws_vapor_coeff%jr)
      if(allocated(eos%iapws_vapor_coeff%nr)) deallocate(eos%iapws_vapor_coeff%nr)
      if(allocated(eos%iapws_m_vapor_coeff%j)) deallocate(eos%iapws_m_vapor_coeff%j)
      if(allocated(eos%iapws_m_vapor_coeff%n)) deallocate(eos%iapws_m_vapor_coeff%n)
      if(allocated(eos%iapws_m_vapor_coeff%ir)) deallocate(eos%iapws_m_vapor_coeff%ir)
      if(allocated(eos%iapws_m_vapor_coeff%jr)) deallocate(eos%iapws_m_vapor_coeff%jr)
      if(allocated(eos%iapws_m_vapor_coeff%nr)) deallocate(eos%iapws_m_vapor_coeff%nr)
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine deteos_amagat(eos,p,t,y1,y2,dv)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(18) ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      type(t_eos2) :: liquid,vapor,gas

      call eos%eos_l(p,t,liquid)
      call eos%eos_v(p,t,vapor)
      call eos%eos_g(p,t,gas)    
      
      dv(3)  = liquid%rho
      dv(4)  = vapor%rho
      dv(5)  = gas%rho

      dv(1)  = 1.d0/((1.d0-y1)/liquid%rho + y1*((1.d0-y2)/vapor%rho + y2/gas%rho))
      dv(2)  = (1.d0-y1)*liquid%h + y1*((1.d0-y2)*vapor%h + y2*gas%h)
      
      dv(7)  = dv(1)**2*((1.d0-y1)/liquid%rho**2*liquid%drdp + &
                     y1*((1.d0-y2)/vapor%rho**2*vapor%drdp  + y2/gas%rho**2*gas%drdp))
      dv(8)  = dv(1)**2*((1.d0-y1)/liquid%rho**2*liquid%drdt + &
                     y1*((1.d0-y2)/vapor%rho**2*vapor%drdt  + y2/gas%rho**2*gas%drdt))
      dv(9)  = dv(1)**2*(1.d0/liquid%rho - ((1.d0-y2)/vapor%rho + y2/gas%rho))
    
      dv(11) = (1.d0-y1)*liquid%dhdp + y1*((1.d0-y2)*vapor%dhdp + y2*gas%dhdp)
      dv(12) = (1.d0-y1)*liquid%dhdt + y1*((1.d0-y2)*vapor%dhdt + y2*gas%dhdt)
      dv(13) = (1.d0-y2)*vapor%h + y2*gas%h - liquid%h
      dv(14) = y1*(gas%h - vapor%h)
      
      dv(6)  = dv(1)*dv(12)/(dv(1)*dv(7)*dv(12)+dv(8)*(1.d0-dv(1)*dv(11)))
      
      dv(15) = vapor%drdp
      dv(16) = vapor%drdt
      dv(17) = liquid%drdp
      dv(18) = liquid%drdt
    
    end subroutine deteos_amagat
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine deteos_dalton(eos,p,t,y1,y2,dv)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t,y1,y2
      real(8), intent(out) :: dv(18) ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
      type(t_eos2) :: liquid,ullage,vapor,gas
      real(8) :: cv
      
      call eos%eos_l(p,t,liquid)
      call eos%eos_v(p,t,vapor)
      call eos%eos_g(p,t,gas)
      
      dv(3)  = liquid%rho
      dv(4)  = vapor%rho
      dv(5)  = gas%rho
      
      dv(13) = (1.d0-y2)*vapor%h +y2*gas%h - liquid%h !warning
      dv(14) = y1*(gas%h - vapor%h)                   !warning
      
      dv(15) = vapor%drdp
      dv(16) = vapor%drdt
      dv(17) = liquid%drdp
      dv(18) = liquid%drdt
      
      if (y1.eq.0.d0) then
        dv(1) = dv(3)
        dv(2) = liquid%h
        
        dv(7) = liquid%drdp
        dv(8) = liquid%drdt
        dv(11) = liquid%dhdp
        dv(12) = liquid%dhdt
        
        dv(9)  = dv(1)**2*(1.d0/dv(3) - 1.d0/dv(4))   !warning 
        dv(10) = dv(1)**2*(1.d0/dv(3) - 1.d0/dv(5))   !warning
        
        dv(6) = dv(12)/(liquid%cv*dv(7))
      else
        call eos%srk_m(p,t,y2,ullage)
        dv(1)  = 1.d0/((1.d0-y1)/liquid%rho + y1/ullage%rho)
        dv(2)  = (1.d0-y1)*liquid%h + y1*ullage%h
        
        dv(7)  = dv(1)**2*((1.d0-y1)/liquid%rho**2*liquid%drdp + y1/ullage%rho**2*ullage%drdp )  
        dv(8)  = dv(1)**2*((1.d0-y1)/liquid%rho**2*liquid%drdt + y1/ullage%rho**2*ullage%drdt )  
        dv(11) = (1.d0-y1)*liquid%dhdp + y1*ullage%dhdp
        dv(12) = (1.d0-y1)*liquid%dhdt + y1*ullage%dhdt
                
        dv(9)  = dv(1)**2*(1.d0/dv(3) - 1.d0/ullage%rho) 
        dv(10) = dv(1)**2*y1*(1.d0/dv(4) - 1.d0/dv(5))      
        
        cv = (1.d0-y1)*liquid%cv + y1*ullage%cv
        dv(6) = dv(12)/(cv*dv(7))
      end if

    end subroutine deteos_dalton
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_iapws97(eos)
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
     
    end subroutine set_iapws97
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_srk_property(eos,fluid,gas)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer, intent(in) :: fluid,gas
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

      select case(gas)
      case(1)
        eos%srk_vapor_property = srk_property
      case(2)
        eos%srk_gas_property = srk_property
      case default
      end select
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
    subroutine set_srk_binary(eos,fluid,ngas)
      implicit none
      class(t_eos), intent(inout) :: eos
      integer, intent(in) :: fluid, ngas
      
      ! n2-h2,n2-he,o2-he,h2-he
      ! binary interaction parameter k_ij for calculating the second cross-virial coefficients of mixtures
      ! long meng, yuan-yuan duan, xiao-dong wang 
      ! above are for virial, actually!!
      
      ! n2-o2
      ! vapor-liquid equilibria of mixtures containing n2, o2, co2, and ethane.
      ! j.stoll, j.vrabec, h.hasse
      
      ! o2-h2
      ! do not use this combination
           
      select case(ngas)
      case(2) ! nitrogen
        select case(fluid)
        case(2) ! nitrogen
          eos%srk_gas_property%k = 0.d0
        case(3) ! oxygen
          eos%srk_gas_property%k = -0.00978d0 
        case(4) ! hydrogen 
          eos%srk_gas_property%k = 0.037d0
        end select
        
      case(3) ! oxygen
        select case(fluid)
        case(2) ! nitrogen
          eos%srk_gas_property%k = -0.00978d0 
        case(3) ! oxygen
          eos%srk_gas_property%k = 0.d0
        case(4) ! hydrogen 
          eos%srk_gas_property%k = 0.d0
        end select
        
      case(4) ! hydrogen
        select case(fluid)
        case(2) ! nitrogen
          eos%srk_gas_property%k = 0.037d0
        case(3) ! oxygen
          eos%srk_gas_property%k = 0.d0
        case(4) ! hydrogen 
          eos%srk_gas_property%k = 0.d0
        end select
        
      case(5) ! helium
        select case(fluid)
        case(2) ! nitrogen
          eos%srk_gas_property%k = 0.281d0
        case(3) ! oxygen
          eos%srk_gas_property%k = 0.248d0
        case(4) ! hydrogen 
          eos%srk_gas_property%k = 0.328d0   
        end select
      end select
      
      eos%srk_vapor_property%k = 0.d0
    end subroutine set_srk_binary
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    subroutine srk_v(eos,p,t,vapor)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: vapor
      type(t_srk_coeff) :: srk_coeff
      real(8) :: a,alpha,b,mw,s,tc,pc,cv0
      real(8) :: rho,e,h,pv
      real(8) :: dpdt,dpdr
      real(8) :: daadt,dsaadt,daidt
      real(8) :: ddaaddt,ddsaaddt,ddaiddt
      real(8) :: cv,cp,bp,bt
      real(8),parameter :: ru=8.314d0
      
      call set_srk_coeff(p,t,eos%srk_vapor_property,srk_coeff)
            
      a     = srk_coeff%a
      alpha = srk_coeff%alpha
      b     = srk_coeff%b
      mw    = eos%srk_vapor_property%mw 
      s     = srk_coeff%s
      tc    = eos%srk_vapor_property%tc
      pc    = eos%srk_vapor_property%pc
      cv0   = eos%srk_vapor_property%cv0
      
      pv = computepv(tc,pc,eos%srk_vapor_property%tb,t)
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

      e = cv0*t + (t*daadt-a*alpha)*dlog(1.d0+b*rho/mw)/(b*mw)
      h = e + dpdr*(1.d0-rho*bp)       
      
      vapor%rho  = rho
      vapor%h    = h
      vapor%drdp = 1.d0/dpdr
      vapor%drdt = -dpdt/dpdr
      vapor%dhdp = bp
      vapor%dhdt = cp
      vapor%cv   = cv
      
    end subroutine srk_v
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine srk_g(eos,p,t,gas)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: gas
      type(t_srk_coeff) :: srk_coeff
      real(8),parameter :: ru=8.314d0
      real(8) :: a,alpha,b,mw,s,tc,pc,cv0
      real(8) :: rho,e,h,pv
      real(8) :: dpdt,dpdr
      real(8) :: daadt,dsaadt,daidt
      real(8) :: ddaaddt,ddsaaddt,ddaiddt
      real(8) :: cv,cp,bp,bt,ap,ar
      
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
      
      gas%rho  = rho
      gas%h    = h
      gas%drdp = 1.d0/dpdr
      gas%drdt = -dpdt/dpdr
      gas%dhdp = bp
      gas%dhdt = cp
      gas%cv   = cv
        
    end subroutine srk_g
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
    subroutine srk_l(eos,p,t,liquid)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid
      type(t_srk_coeff) :: srk_coeff
      real(8),parameter :: ru=8.314d0
      real(8) :: a,alpha,b,mw,s,tc,pc,cv0
      real(8) :: rho,e,h,pv
      real(8) :: dpdt,dpdr
      real(8) :: daadt,dsaadt,daidt
      real(8) :: ddaaddt,ddsaaddt,ddaiddt
      real(8) :: cv,cp,bp,bt,ap,ar
      
      call set_srk_coeff(p,t,eos%srk_vapor_property,srk_coeff)
            
      a     = srk_coeff%a
      alpha = srk_coeff%alpha
      b     = srk_coeff%b
      mw    = eos%srk_vapor_property%mw 
      s     = srk_coeff%s
      tc    = eos%srk_vapor_property%tc
      pc    = eos%srk_vapor_property%pc
      cv0   = eos%srk_vapor_property%cv0
      
      pv = computepv(tc,pc,eos%srk_vapor_property%tb,t)
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
      
      liquid%rho  = rho
      liquid%h    = h
      liquid%drdp = 1.d0/dpdr
      liquid%drdt = -dpdt/dpdr
      liquid%dhdp = bp
      liquid%dhdt = cp
      liquid%cv   = cv
        
    end subroutine srk_l
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
    subroutine srk_m(eos,p,t,y2,ullage)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t,y2
      type(t_eos2), intent(out) :: ullage
      type(t_srk_coeff) :: srk_gas_coeff,srk_vapor_coeff
      real(8),parameter :: ru=8.314d0
      real(8) :: ai,alphai,bi,mwi,si,tci,pci,cv0i,xi
      real(8) :: aj,alphaj,bj,mwj,sj,tcj,pcj,cv0j,xj,kij
      real(8) :: alphaaij,aalpha,b,mw,cv0
      real(8) :: rho,e,h,pv
      real(8) :: dpdt,dpdr
      real(8) :: daadt,dsaadt,daidt,dajdt
      real(8) :: ddaaddt,ddsaaddt,ddaiddt,ddajddt
      real(8) :: daaijdt,e_p_den_i,e_p_den_j,h_p_mass_i,h_p_mass_j,h_p_den_i,h_p_den_j,dpdri,dpdrj
      real(8) :: cv,cp,bp,bt
      
      ! get coefficients of each constituent gas
      call set_srk_coeff(p,t,eos%srk_gas_property,srk_gas_coeff)
      call set_srk_coeff(p,t,eos%srk_vapor_property,srk_vapor_coeff)

      ai     = srk_vapor_coeff%a
      alphai = srk_vapor_coeff%alpha
      bi     = srk_vapor_coeff%b
      mwi    = eos%srk_vapor_property%mw 
      si     = srk_vapor_coeff%s
      tci    = eos%srk_vapor_property%tc
      pci    = eos%srk_vapor_property%pc
      cv0i   = eos%srk_vapor_property%cv0
      
      aj     = srk_gas_coeff%a
      alphaj = srk_gas_coeff%alpha
      bj     = srk_gas_coeff%b
      mwj    = eos%srk_gas_property%mw 
      sj     = srk_gas_coeff%s
      tcj    = eos%srk_gas_property%tc
      pcj    = eos%srk_gas_property%pc
      cv0j   = eos%srk_gas_property%cv0

      kij    = eos%srk_gas_property%k
  
      pv = computepv(tci,pci,eos%srk_vapor_property%tb,t)
      !if ((p.le.pv).and.(t.le.tci)) then
      if ((p.le.pv)) then
        alphai = pv/pci
      end if
      
      pv = computepv(tcj,pcj,eos%srk_gas_property%tb,t)
      !if ((p.le.pv).and.(t.le.tcj)) then
      if ((p.le.pv)) then
        alphaj = pv/pcj
      end if
      
      ! calculate mole fraction
      mw = 1.d0/((1.d0-y2)/mwi + y2/mwj)
      xi = (1.d0-y2)*mw/mwi
      xj = y2*mw/mwj
      
      ! find mixture property
      alphaaij = dsqrt(alphai*alphaj*ai*aj)*(1.d0-kij)
      aalpha   = xi**2*alphai*ai + 2.d0*xi*xj*alphaaij + xj**2*alphaj*aj
      b        = xi*bi + xj*bj
      cv0      = xi*cv0i + xj*cv0j

      rho = cubicsolver(p,t,aalpha,b,mw)
      
      ! find daadt for dpdt
      daidt  = -si/dsqrt(t*tci)*(1.d0+si*(1.d0-dsqrt(t/tci)))
      dajdt  = -sj/dsqrt(t*tcj)*(1.d0+sj*(1.d0-dsqrt(t/tcj))) 
      dsaadt = 0.5d0*(dsqrt(alphai/alphaj)*dajdt + dsqrt(alphaj/alphai)*daidt)
      daadt  = xi**2*ai*daidt + 2.d0*xi*xj*dsqrt(ai*aj)*dsaadt + xj**2*aj*dajdt
      
      ! dpdt, dpdr
      dpdt = rho*ru/(mw-b*rho) - daadt*rho**2/(mw*(mw+b*rho))
      dpdr = mw*ru*t/(mw-b*rho)**2 - aalpha*rho*(2.d0*mw+b*rho)/(mw*(mw+b*rho)**2)
      dpdri = mw*ru*t/(mwi*(mw-b*rho)**2)*(mw+rho*(bi-b)) &
              - 2.d0*rho*xj*alphaaij/(mwi*(mw+b*rho)) &
              + aalpha*rho**2*bi/(mwi*(mw+b*rho)**2)
      dpdrj = mw*ru*t/(mwj*(mw-b*rho)**2)*(mw+rho*(bj-b)) &
              - 2.d0*rho*xi*alphaaij/(mwj*(mw+b*rho)) &
              + aalpha*rho**2*bj/(mwj*(mw+b*rho)**2)
      
      ! find ddaaddt for cv
      ddaiddt  = 0.5d0*si * (si/(t*tci) + (1.d0+si*(1.d0-dsqrt(t/tci)))/dsqrt(t**3*tci))
      ddajddt  = 0.5d0*sj * (sj/(t*tcj) + (1.d0+sj*(1.d0-dsqrt(t/tcj)))/dsqrt(t**3*tcj))
      ddsaaddt = 0.5d0 * (daidt*dajdt/dsqrt(alphai*alphaj)                &
                         - 0.5d0*dsqrt(alphai/alphaj**3)*dajdt**2         &
                         - 0.5d0*dsqrt(alphaj/alphai**3)*daidt**2         &
                         + dsqrt(alphai/alphaj)*ddajddt                   &
                         + dsqrt(alphaj/alphai)*ddaiddt )
      ddaaddt  = xi**2*ai*ddaiddt + 2.d0*xi*xj*dsqrt(ai*aj)*ddsaaddt + xj**2*aj*ddajddt
      
      ! cv for dhdt & e
      cv = cv0 + t*ddaaddt*dlog(1.d0+b*rho/mw)/(b*mw)
      cp = cv + t*dpdt**2/(rho**2*dpdr)
      
      ! dhdp, dhdt
      bp = 1.d0/rho - p/(rho**2*dpdr)
      bt = cv + dpdt*p/(rho**2*dpdr)

      ! daaijdt for partial density energy(e_p_den)
      daaijdt = dsqrt(ai*aj)*dsaadt
      
      ! e_p_den, h_p_den, h_p_mass
      e_p_den_i = cv0i*t + ( 2.d0*(xj*(t*daaijdt-alphaaij))*dlog(1.d0+b*rho/mw) &
                           + bi*(t*daadt-aalpha)*(rho/(mw+b*rho)-dlog(1.d0+b*rho/mw)/b) ) / (b*mwi)
      e_p_den_j = cv0j*t + ( 2.d0*(xi*(t*daaijdt-alphaaij))*dlog(1.d0+b*rho/mw) &
                           + bj*(t*daadt-aalpha)*(rho/(mw+b*rho)-dlog(1.d0+b*rho/mw)/b) ) / (b*mwj)  
      h_p_den_i = e_p_den_i + dpdri
      h_p_den_j = e_p_den_j + dpdrj
      h_p_mass_i = e_p_den_i + dpdri*(1.d0-rho*bp)
      h_p_mass_j = e_p_den_j + dpdrj*(1.d0-rho*bp)
      
      e = cv0*t + (t*daadt-aalpha)*dlog(1.d0+b*rho/mw)/(b*mw)
      h = e + dpdr*(1.d0-rho*bp)       
            
      ullage%rho  = rho
      ullage%h    = h
      ullage%drdp = 1.d0/dpdr
      ullage%drdt = -dpdt/dpdr
      ullage%dhdp = bp
      ullage%dhdt = cp
      ullage%cv   = cv
      
    end subroutine srk_m   
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
    subroutine ideal(eos,p,t,gas)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: gas
      real(8), parameter :: r1 = 1.4d0/0.4d0/1004.64
      real(8) :: tau
      
      tau = 1.d0/t
      
      gas%rho  = p*r1*tau
      gas%h    = 1004.64d0*t
      gas%drdp = r1*tau
      gas%drdt = -p*r1*tau**2
      gas%dhdp = 0.d0
      gas%dhdt = 1004.64d0
      
    end subroutine ideal
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2o_l_fit(eos,p,t,liquid)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid
      
      liquid%rho = - 6311.024687022091d0+105.96719226625858d0*t-0.6180387750363456d0*t**2+0.0018192120809256687d0*t**3- &
                 2.7046566485973865d-6*t**4 + 1.61722220645795d-9*t**5 + 0.000029080581316291686d0*p -                  &
                 3.3309022145970827d-7*t*p +1.4640381037355473d-9*t**2*p -                                              &
                 2.8865979129253496d-12*t**3*p + 2.156808271941056d-15*t**4*p -                                         &
                 2.6280680541055985d-14*p**2 + 2.297216572824232d-16*t*p**2 -                                           &
                 6.747504264295466d-19*t**2*p**2 + 6.525442466633386d-22*t**3*p**2
                 
      liquid%h = 1.3557612872909496d9 - 1.363266289856961d7*t + 40681.30627818539d0*t**2 - 41.77799452978375d0*t**3 +                          &
               0.016078938509377864d0*t**4 - 1.6731139882275142d-6*t**5 - 3.2351874957958134d6*liquid%rho + 31394.61628884811d0*t*liquid%rho - &
               83.9426871170273d0*t**2*liquid%rho + 0.06306446061998854d0*t**3*liquid%rho - 0.000013263064759017577d0*t**4*liquid%rho +        &
               2553.4878230058544d0*liquid%rho**2 - 23.88678218373101d0*t*liquid%rho**2 + 0.055714123710032606d0*t**2*liquid%rho**2 -          &
               0.000023226414329920397d0*t**3*liquid%rho**2 - 0.6662015780399999d0*liquid%rho**3 + 0.006005285813282444d0*t*liquid%rho**3 -    &
               0.000011766537107697452d0*t**2*liquid%rho**3
      liquid%drdp = 0.028344984598699954d0 - 0.0002652384870719236d0*t + 8.069978640361292d-7*t**2 -                                        &
                  8.975533323360278d-10*t**3 + 4.2088568723050893d-13*t**4 - 7.187943097341636d-17*t**5 -                                   &
                  0.00006936637304995329d0*liquid%rho + 6.119809907043226d-7*t*liquid%rho - 1.6429994317814744d-9*t**2*liquid%rho +         &
                  1.2949079244023172d-12*t**3*liquid%rho - 3.117135505170058d-16*t**4*liquid%rho + 5.665107930759331d-8*liquid%rho**2 -     &
                  4.69641704021055d-10*t*liquid%rho**2+1.0885462450668403d-12*t**2*liquid%rho**2-4.63961788687922d-16*t**3*liquid%rho**2 -  &
                  1.543349645193048d-11*liquid%rho**3 + 1.1974984976879823d-13*t*liquid%rho**3 - 2.3212446257875464d-16*t**2*liquid%rho**3
      liquid%drdt = 121879.9584706434d0-1139.2925141420403d0*t+3.4530261160544415d0*t**2 - 0.003775156632058203d0*t**3+                    &
                  1.712735110674499d-6*t**4 - 2.7555948221273537d-10*t**5 - 299.6006492602664d0*liquid%rho +                               &
                  2.646621331302561d0*t*liquid%rho -0.007107401898351704d0*t**2*liquid%rho + 5.532737555143375d-6*t**3*liquid%rho -        &
                  1.300382894492162d-9*t**4*liquid%rho + 0.24544365927014916d0*liquid%rho**2 - 0.002041380849491402d0*t*liquid%rho**2 +    &
                  4.7506749483188065d-6*t**2*liquid%rho**2 - 2.0054026031199562d-9*t**3*liquid%rho**2 -                                    &
                  0.00006702275389171266d0*liquid%rho**3 + 5.22658483578724d-7*t*liquid%rho**3 - 1.021429064669065d-9*t**2*liquid%rho**3
      liquid%dhdp = 28.256990502588454d0 - 0.258386234711587d0*t + 0.0007698968883972372d0*t**2 -                                         &
                  8.299018578366212d-7*t**3 + 3.7225797441168453d-10*t**4 - 5.932327720843982d-14*t**5 -                                  &
                  0.0699434126351822d0*liquid%rho +  0.0006037175382305189d0*t*liquid%rho - 1.5925877833800108d-6*t**2*liquid%rho +       &
                  1.2213368187673786d-9*t**3*liquid%rho - 2.836880571012265d-13*t**4*liquid%rho + 0.00005771616450097861d0*liquid%rho**2- &
                  4.6843157322323967d-7*t*liquid%rho**2 + 1.0699805810683837d-9*t**2*liquid%rho**2 -                                      &
                  4.445120270945741d-13*t**3*liquid%rho**2 - 1.588026958029046d-8*liquid%rho**3 +                                         &
                  1.2067728903792023d-10*t*liquid%rho**3 - 2.313101693533085d-13*t**2*liquid%rho**3                      
      liquid%dhdt = 2.418822477923761d7+35108.447996488154d0*t - 873.5689347951068d0*t**2 + 1.7230272675151386d0*t**3-                &
                  0.0010893940589922265d0*t**4 + 2.177699794807627d-7*t**5 - 81178.10072386597d0*liquid%rho +                         &
                  90.42654563938547d0*t*liquid%rho + 1.4007102727849796d0*t**2*liquid%rho - 0.0022661796428748962d0*t**3*liquid%rho + &
                  7.777074422878771d-7*t**4*liquid%rho + 85.17031674149585d0*liquid%rho**2 - 0.20952117580112675d0*t*liquid%rho**2 -  &
                  0.000655595931910761d0*t**2*liquid%rho**2 + 7.213406726853327d-7*t**3*liquid%rho**2 -                               &
                  0.028586662251862786d0*liquid%rho**3 + 0.00009122145304724242d0*t*liquid%rho**3 +                                   &
                  7.761232305034003d-8*t**2*liquid%rho**3  
    end subroutine h2o_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2o_v_fit(eos,p,t,vapor)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: vapor
      
      if(t.le.347.d0) then
        vapor%rho = 0.0046086722990453105d0 - 0.00004151388632156331d0*t + 1.2413429087175857d-7*t**2 - 1.2306907168991188d-10*t**3 &
                  + 0.000020178322982098646d0*p - 6.247722798143032d-8*t*p                                                          &
                  + 6.431744010624924d-11*t**2*p+ 1.18294435757289d-11*p**2                                                         &
                  - 3.0114803183388726d-14*t*p**2
        vapor%h = 1.6140385439321147d6 + 8052.693799326647d0*t - 40.01547777953337d0*t**2 + 0.1282203670472617d0*t**3 -                       &
                0.00020429579076161367d0*t**4 + 1.2999110502072536d-7*t**5 - 4.6229047884475544d7*vapor%rho + 556220.6222673394d0*t*vapor%rho &  
                -2518.9396226736376d0*t**2*vapor%rho + 5.080373140251466d0*t**3*vapor%rho - 0.0038468216423550375d0*t**4*vapor%rho -          &
                1.4102180069315466d8*vapor%rho**2 + 1.260420439726543d6*t*vapor%rho**2 - 3753.392034834953d0*t**2*vapor%rho**2 +              &
                3.7239906172249513d0*t**3*vapor%rho**2 - 5.175977413572329d7*vapor%rho**3 + 300383.7348311928d0*t*vapor%rho**3 -              &
                435.84557503111483d0*t**2*vapor%rho**3
        vapor%drdp = 0.000052009410802462945d0 - 5.004533737095439d-7*t + 2.496157019891614d-9*t**2 -                                         &
                        6.850414374885032d-12*t**3 + 9.850424673683412d-15*t**4 - 5.81580779044606d-18*t**5 +                                 &
                        0.0007468297964490092d0*vapor%rho - 8.945009619536416d-6*t*vapor%rho + 4.035865049281818d-8*t**2*vapor%rho -          &
                        8.11581415175705d-11*t**3*vapor%rho + 6.131099812799595d-14*t**4*vapor%rho +  0.0018361131184468887d0*vapor%rho**2 -  &
                        0.000016425813360828023d0*t*vapor%rho**2 + 4.896642737873156d-8*t**2*vapor%rho**2 -                                   &
                        4.864021269344103d-11*t**3*vapor%rho**2 + 0.0006498859974798217d0*vapor%rho**3 -                                      &
                        3.7770305631556957d-6*t*vapor%rho**3 + 5.488322941633847d-9*t**2*vapor%rho**3
        vapor%dhdt = 69740.2862537548d0 - 1106.8606451381106d0*t + 7.173138161430398d0*t**2 -                                                 &
                    0.023112000968498456d0*t**3 + 0.00003704343835272604d0*t**4 - 2.3632785520084313d-8*t**5 +                                &
                    7.06467480366046d6*vapor%rho - 86531.28933536631d0*t*vapor%rho + 397.27353974813946d0*t**2*vapor%rho -                    &
                    0.8100829096462536d0*t**3*vapor%rho + 0.0006189421808640733d0*t**4*vapor%rho + 2.4331535736487016d7*vapor%rho**2 -        &
                    217697.33317293136d0*t*vapor%rho**2 + 648.8848821623319d0*t**2*vapor%rho**2 -                                             &
                    0.6443472958123013d0*t**3*vapor%rho**2 + 8.957917762373075d6*vapor%rho**3 -                                               &
                    52014.21123588152d0*t*vapor%rho**3 + 75.50955185252435d0*t**2*vapor%rho**3 
      else if((t.gt.347.d0).and.(t.le.497.d0)) then
        vapor%rho = 0.0304345716654541d0 -0.000217065175507298d0*t + 5.17227679271671d-07*t**2 - 4.07615934825363d-10*t**3 +        &
                    0.0000206853043135603d0*p - 7.37792204298515d-08*t*p                                                            &
                    + 1.16426003692796d-10*t**2*p - 6.85486255365125d-14*t**3*p                                                     &
                    + 3.02040748441321d-11*p**2 - 1.78331822460195d-13*t*p**2                                                       &
                    + 3.57437719439424d-16*t**2*p**2 - 2.41822055967549d-19*t**3*p**2                                               &
                    +5.17216589243678d-18*p**3 - 2.10471393810734d-20*t*p**3                                                        &
                    +2.15338023622516d-23*t**2*p**3 + 1.97005590122501d-25*p**4 -                                                   &
                    4.05774798589835d-28*t*p**4 + 2.12852285936374d-33*p**5
        vapor%h = 1672315.67433477d0 + 4955.3297331043d0*t - 11.1328027215681d0*t**2 + 0.0175990883625111d0*t**3 - &
                  0.0000101449179466104d0*t**4 - 574582.768022701d0*vapor%rho + 3426.70214256579d0*t*vapor%rho - 7.02840680468889d0*t**2*vapor%rho + &
                  0.00487707943643386d0*t**3*vapor%rho - 21816.9278267704d0*vapor%rho**2 + 89.2301619577232d0*t*vapor%rho**2 -&
                  0.0915244176282538d0*t**2*vapor%rho**2 - 124.509070972886d0*vapor%rho**3 + 0.248625694075255d0*t*vapor%rho**3
        vapor%drdp = 0.0000214068495783794d0 - 7.89132492686781d-08*t + 1.28600002435269d-10*t**2 - 7.81509720014333d-14*t**3 +        &
                     9.70273812432014d-06*vapor%rho - 5.78351768072687d-08*t*vapor%rho + 1.17602375083603d-10*t**2*vapor%rho -                        &
                     8.08411915570631d-14*t**3*vapor%rho + 5.40779755406623d-07*vapor%rho**2 - 2.21646856027985d-09*t*vapor%rho**2 +                  &
                     2.2828969468728d-12*t**2*vapor%rho**2 + 4.99113366855691d-09*vapor%rho**3 - 1.02123051912891d-11*t*vapor%rho**3 +                &
                     9.44914641203704d-12*vapor%rho**4
        vapor%dhdt = -3033.80716978799d0 + 45.3409779082821d0*t - 0.157465860914466d0*t**2 + 0.000243173696121016d0*t**3 -            &
                     1.39392770211301d-07*t**4 + 98937.3024214889d0*vapor%rho - 891.054167343861d0*t*vapor%rho + 3.02076705048526d0*t**2*vapor%rho -&
                     0.00455702585759261d0*t**3*vapor%rho + 2.57741825292208d-06*t**4*vapor%rho + 15936.2394936362d0*vapor%rho**2 -                 &
                     104.969845431369d0*t*vapor%rho**2 + 0.230176196829332d0*t**2*vapor%rho**2 - 0.000167962640251061d0*t**3*vapor%rho**2 +         &
                     479.361786337656d0*vapor%rho**3 - 2.03444799195589d0*t*vapor%rho**3 + 0.00215388853083512d0*t**2*vapor%rho**3 +                &
                     2.99492096741046d0*vapor%rho** 4 - 0.00609386126945215d0*t*vapor%rho**4 + 0.00214502154063786d0*vapor%rho**5
      else
        vapor%rho = 13.197225283259d0 - 0.0781811172267488d0*t + 0.000154184911560128d0*t**2 - 1.01228091394607d-07*t**3 -          &
                    0.0000425337828013887d0*p + 3.04665239562157d-07*t*p                                                            &
                    - 6.38086984502479d-10*t**2*p + 4.32474213995152d-13*t**3*p                                                     &
                    + 4.73157952125924d-11*p**2 - 2.78968106119848d-13*t*p**2                                                       &
                    + 5.52729713074471d-16*t**2*p**2 - 3.66680108954259d-19*t**3*p**2                                               &
                    + 3.76386771510989d-18*p**3 - 1.45063542714493d-20*t*p**3                                                       &
                    + 1.40177910612357d-23*t**2*p**3 + 7.82536215892866d-26 *p**4                                                   &
                    - 1.49391037531458d-28*t*p**4   + 4.57130226214232d-34*p**5
        vapor%h = 7139235.2818126d0 - 39461.1047270727d0*t + 123.016364051649d0*t**2 - 0.161216479038991d0*t**3 +                       &
                  0.0000787398895903491d0*t**4 - 2133191.8212778d0*vapor%rho + 16504.5657953198d0*t*vapor%rho - 48.4107131761524d0*t**2*vapor%rho +    &
                  0.0634070302544731d0*t**3*vapor%rho - 0.0000312157773910099d0*t**4*vapor%rho - 54671.1142939832d0*vapor%rho**2 +                     &
                  312.54941433458d0*t*vapor%rho**2 - 0.597723026381932d0*t**2*vapor%rho**2 + 0.000381926351635443d0*t**3*vapor%rho**2 -                &
                  301.061606537303d0*vapor%rho**3 + 1.15970794808449d0*t*vapor%rho**3 - 0.00111190863641621d0*t**2*vapor%rho**3 -                      &
                  0.651906676977929d0*vapor%rho**4 + 0.00118075876509709d0*t*vapor%rho**4
        vapor%drdp = 0.000027646722992399d0 - 1.11926434493766d-07*t + 1.85279270247462d-10*t**2 - 1.09363229043887d-13*t**3 -          &
                     0.0000070118253298229d0*vapor%rho + 4.41078803695323d-08*t*vapor%rho - 8.99397396544076d-11*t**2*vapor%rho +                      &
                     6.02028263785701d-14*t**3*vapor%rho + 1.60429170974966d-06*vapor%rho**2 - 9.31573281384837d-09*t*vapor%rho**2 +                   &
                     1.80830036391113d-11*t**2*vapor%rho**2 - 1.17202189573276d-14*t**3*vapor%rho**2 + 1.57622617451259d-08*vapor%rho**3 -             &
                     5.99883628772789d-11*t*vapor%rho**3 + 5.70918830851557d-14*t**2*vapor%rho**3 + 4.1341541526741d-11*vapor%rho**4 -                 &
                     7.71586993101989d-14*t*vapor%rho**4 + 2.58357488221122d-14*vapor%rho**5
        vapor%dhdt = 54122.0046918165d0 - 302.437279261426d0*t + 0.582869886344912d0*t**2 - 0.000373348944239481d0*t**3 -             &
                     38696.5616301691d0*vapor%rho + 227.607080757579d0*t*vapor%rho - 0.444233634160479d0*t**2*vapor%rho +                            &
                     0.000288112945516371d0*t**3*vapor%rho + 5420.14763683506d0*vapor%rho**2 - 31.902217596953d0*t*vapor%rho**2 +                    &
                     0.0625549751764669d0*t**2*vapor%rho**2 - 0.0000408485848158127d0*t**3*vapor%rho**2 + 49.0212719421652d0*vapor%rho**3 -          &
                     0.187189705044773d0*t*vapor%rho**3 + 0.000178485823529128d0*t**2*vapor%rho**3 + 0.110318082427605d0*vapor%rho**4 -              &
                     0.00020498071909394d0*t*vapor%rho**4 + 0.0000507096802350709d0*vapor%rho**5
      end if
      if (t.lt.296.d0) then                
        vapor%drdt = 0.0015480993389557069d0 - 0.00002204315012405441d0*t + 1.1761132674188987d-7*t**2 - 2.786961066751213d-10*t**3     &
                + 2.474814735018812d-13*t**4 - 0.8001981115028577d0*vapor%rho + 0.011427865347741394d0*t*vapor%rho                           &
                - 0.00006157522268531719d0*t**2*vapor%rho + 1.475339596758024d-7*t**3*vapor%rho - 1.3248641132424813d-10*t**4*vapor%rho           &
                - 18.823422434851317d0*vapor%rho**2 + 0.19747534922119145d0*t*vapor%rho**2 - 0.0006904434493780153d0*t**2*vapor%rho**2            &
                + 8.044532992903289d-7*t**3*vapor%rho**2 - 94.28088325558488d0*vapor%rho**3 + 0.6459714529858793d0*t*vapor%rho**3                 &
                - 0.001106113187641709d0*t**2*vapor%rho**3 - 99.90491610726272d0*vapor%rho**4 + 0.33224232928388087d0*t*vapor%rho**4              &
                + 12.558583095426068d0*vapor%rho**5
        vapor%dhdp = -42948.92139971102d0 + 774.1474629248053d0*t - 5.574937750585217d0*t**2 +                                        &
                    0.020049561101358117d0*t**3 - 0.00003601017265315555d0*t**4 + 2.5840711650456145d-8*t**5 -                     &
                    3.8795311376638315d6*vapor%rho + 53841.27611225954d0*t*vapor%rho - 280.16500438649547d0*t**2*vapor%rho +                      &
                    0.647830494592685d0*t**3*vapor%rho - 0.0005616580339274063d0*t**4*vapor%rho - 2.529432983615679d7*vapor%rho**2+               &
                    257585.75930428825d0*t*vapor%rho**2 - 874.587299851378d0*t**2*vapor%rho**2 +                                             &
                    0.9900592360302818d0*t**3*vapor%rho**2     
        
      else if((t.ge.296.d0).and.(t.lt.314.d0)) then
        vapor%drdt = 0.00045711848406281485 - 6.574332974047608d-6*t + 3.517329163462288d-8*t**2 - 8.306194636795476d-11*t**3 +         &
                7.31206542355063d-14*t**4 - 0.27583036049813975d0*vapor%rho + 0.003717655093135968d0*t*vapor%rho -                           &
                0.000019128186158450722d0*t**2*vapor%rho + 4.384516304283248d-8*t**3*vapor%rho - 3.766467591730544d-11*t**4*vapor%rho -           &
                4.775668651062799d0*vapor%rho**2 + 0.04712225492912923d0*t*vapor%rho**2 - 0.00015503387171534767d0*t**2*vapor%rho**2 +            &
                1.7003113271227604d-7*t**3*vapor%rho**2 - 10.931939154387278d0*vapor%rho**3 + 0.07047793818002629d0*t*vapor%rho**3 -              &
                0.00011358413099484183d0*t**2*vapor%rho**3 - 4.537303275769487d0*vapor%rho**4 + 0.014330905251429124d0*t*vapor%rho**4 -           &
                0.02782522311180477d0*vapor%rho**5
      else if((t.ge.314.d0).and.(t.lt.332.d0)) then
        vapor%drdt = 0.002479474833613407d0 - 0.00003187828249623505d0*t + 1.5348441095810086d-7*t**2 - 3.2801531309942245d-10*t**3     &
                + 2.62559269178044d-13*t**4 - 0.2941904846101775d0*vapor%rho + 0.0037224237724600763d0*t*vapor%rho                           &
                - 0.000017955519612122393d0*t**2*vapor%rho + 3.860049618527247d-8*t**3*vapor%rho - 3.112086584570705d-11*t**4*vapor%rho           &
                - 1.7531392341348808d0*vapor%rho**2 + 0.016337040762004425d0*t*vapor%rho**2 - 0.00005080207895538633d0*t**2*vapor%rho**2          &
                + 5.2690873875885054d-8*t**3*vapor%rho**2 - 1.6918721657663482d0*vapor%rho**3 + 0.010352644150310147d0*t*vapor%rho**3             &
                - 0.00001583783934071373d0*t**2*vapor%rho**3 - 0.34622902695411484d0*vapor%rho**4 + 0.001042069442089262d0*t*vapor%rho**4         &
                - 0.007278966282372564d0*vapor%rho**5
      else if((t.ge.332.d0).and.(t.le.347.d0)) then
        vapor%drdt = 0.0027349989184879747d0 - 0.00003418171573287271d0*t + 1.5952183417150958d-7*t**2 - 3.296335553972185d-10*t**3     &
                + 2.545725331783d-13*t**4 - 0.21908046643019907d0*vapor%rho + 0.0026310597719816555d0*t*vapor%rho                            &
                - 0.000012103906318237557d0*t**2*vapor%rho + 2.4840583043625168d-8*t**3*vapor%rho - 1.9122677235546308d-11*t**4*vapor%rho         &
                - 0.6965991472905452d0*vapor%rho**2 + 0.0061551518242377095d0*t*vapor%rho**2 - 0.00001817241916521105d0*t**2*vapor%rho**2         &
                + 1.791163543979636d-8*t**3*vapor%rho**2 - 0.32467128630581d0*vapor%rho**3 + 0.0018957809746163372d0*t*vapor%rho**3               &
                - 2.768681600604109d-6*t**2*vapor%rho**3 - 0.03814270926478833d0*vapor%rho**4 + 0.00011037670500046392d0*t*vapor%rho**4           &
                - 0.0010845296436386492d0*vapor%rho**5
      else if((t.gt.347.d0).and.(t.le.397.d0)) then
        vapor%drdt = -0.0000491213869994715d0 + 4.3699109475449d-7*t - 1.27801403663623d-9*t**2 + 1.23702103316548d-12*t**3 -           &
                 0.00745673136638664d0*vapor%rho + 0.0000155921168719041d0*t*vapor%rho + 1.11675572608361d-10*t**2*vapor%rho -                    &
                 2.03129818205504d-11*t**3*vapor%rho - 0.0727116915503993d0*vapor%rho**2 + 0.000559946132162153d0*t*vapor%rho**2 -                &
                 1.45506600290402d-6*t**2*vapor%rho**2 + 1.27101008348442d-9*t**3*vapor%rho**2 - 0.0178975505425514d0*vapor%rho**3 +              &
                 0.0000943938898436783d0*t*vapor%rho**3 - 1.24666727989567d-7*t**2*vapor%rho**3 - 0.00143272552687376d0*vapor%rho**4 +            &
                 3.73822031270469d-6*t*vapor%rho**4-0.0000272850579250176d0*vapor%rho**5
      else if((t.gt.397.d0).and.(t.le.447.d0)) then
        vapor%drdt = -0.0116146111379786d0 + 0.0000834812107294513d0*t - 1.99614022037687d-07*t**2 + 1.58811953545974d-10*t**3 +        &
                 0.0552686133839977d0*vapor%rho - 0.000435201052381534d0*t*vapor%rho + 1.07610493924041d-06*t**2*vapor%rho -                      &
                 8.73145090591731d-10*t**3*vapor%rho - 0.0709149440418435d0*vapor%rho**2 + 0.000510724229940245d0*t*vapor%rho**2-                 &
                 1.23091425806775d-06*t**2*vapor%rho**2 + 9.90242686781605d-10*t**3*vapor%rho**2 - 0.00643663811831822*vapor%rho**3+              &
                 0.0000299908146192532d0*t*vapor%rho**3 - 3.49445390054649d-08*t**2*vapor%rho**3 - 0.000126284826162668d0*vapor%rho**4+           &
                 2.88102955025302d-07*t*vapor%rho**4 - 5.55731072635301d-07*vapor%rho**5
      else if((t.gt.447.d0).and.(t.le.497.d0)) then
        vapor%drdt = -0.0864611400507011d0 + 0.000564503319948696d0*t - 1.22681430912334d-06*t**2 + 8.87456451441192d-10*t**3 +         &
                 0.198845609823498d0*vapor%rho - 0.00133530111502184d0*t*vapor%rho + 2.94049042029591d-06*t**2*vapor%rho -                        &
                 2.14699720902986d-09*t**3*vapor%rho - 0.0745875895526571d0*vapor%rho**2 + 0.000496395999451785d0*t*vapor%rho**2 -                &
                 1.10225985771488d-06*t**2*vapor%rho**2 + 8.15329494815788d-10*t**3*vapor%rho**2 - 0.00272923523424554d0*vapor%rho**3 +           &
                 0.0000116460569402577d0*t*vapor%rho**3 - 1.24218143895751d-08*t**2*vapor%rho**3 - 0.0000220275512608702d0*vapor%rho**4 +         &
                 4.57228531090116d-08*t*vapor%rho**4 - 3.62049533964305d-08*vapor%rho**5
        
      else !if((t.gt.497.d0).and.(t.le.552.d0)) then
        vapor%drdt = -0.00971360853462954d0 + 0.0000366376246703219d0*t - 3.4479458553232d-08*t**2 - 0.0000432642035811531d0*vapor%rho -     &
                 0.000010274561150715d0*t*vapor%rho + 1.2766404543913d-08*t**2*vapor%rho - 0.000742479359962173d0*vapor%rho**2 +                  &
                 0.0000019451020297395d0*t*vapor%rho**2 - 1.22586476060632d-09*t**2*vapor%rho**2 - 0.00029193963090101d0*vapor%rho**3 +           &
                 1.09811167895505d-06*t*vapor%rho**3 - 1.03935342359949d-09*t**2*vapor%rho**3 - 1.61008564793834d-06*vapor%rho**4 +               &
                 3.03605870434151d-09*t*vapor%rho**4 - 1.81412095488845e-09*vapor%rho**5
      end if
    
      if (t.lt.296.d0) then
        vapor%dhdp = -42948.92139971102d0 + 774.1474629248053d0*t - 5.574937750585217d0*t**2 +                                        &
                    0.020049561101358117d0*t**3 - 0.00003601017265315555d0*t**4 + 2.5840711650456145d-8*t**5 -                     &
                    3.8795311376638315d6*vapor%rho + 53841.27611225954d0*t*vapor%rho - 280.16500438649547d0*t**2*vapor%rho +                      &
                    0.647830494592685d0*t**3*vapor%rho - 0.0005616580339274063d0*t**4*vapor%rho - 2.529432983615679d7*vapor%rho**2+               &
                    257585.75930428825d0*t*vapor%rho**2 - 874.587299851378d0*t**2*vapor%rho**2 +                                             &
                    0.9900592360302818d0*t**3*vapor%rho**2            
      else if((t.ge.296.d0).and.(t.lt.314.d0)) then        
        vapor%dhdp = -15586.728701871567d0 + 264.0264835470356d0*t - 1.7875045343117237d0*t**2 +                                      &
                    0.006045307503336912d0*t**3 - 0.000010212748709652075d0*t**4 + 6.894542736674656d-9*t**5 -                     &
                    527744.9289123041d0*vapor%rho + 6900.184295261724d0*t*vapor%rho - 33.827285188196655d0*t**2*vapor%rho +                       &
                    0.07369354996628959d0*t**3*vapor%rho - 0.000060195142611365495d0*t**4*vapor%rho -                                        &
                    1.2692441139548968d6*vapor%rho**2 + 12179.77837626357d0*t*vapor%rho**2 -                                                 &
                    38.9691482898955d0*t**2*vapor%rho**2 + 0.04157034917642339d0*t**3*vapor%rho**2
      else if((t.ge.314.d0).and.(t.lt.332.d0)) then
        vapor%dhdp = -5985.323520493672d0 + 95.44189489466376d0*t - 0.6086603702954034d0*t**2 +                                       &
                    0.0019399724785255999d0*t**3 - 3.0898699024029392d-6*t**4 + 1.967274233030398d-9*t**5 -                        &
                    78428.828592011d0*vapor%rho+970.2347151449338d0*t*vapor%rho - 4.500508099945776d0*t**2*vapor%rho +                            &
                    0.00927714313504291d0*t**3*vapor%rho - 7.170466057187251d-6*t**4*vapor%rho - 81033.69763174784d0*vapor%rho**2 +               &
                    736.6191168105244d0*t*vapor%rho**2 - 2.2326112086981085d0*t**2*vapor%rho**2 +                                            &
                    0.0022561595659607156d0*t**3*vapor%rho**2 - 1845.1443282315129d0*vapor%rho**3 +                                          &
                    11.233381579730002d0*t*vapor%rho**3 - 0.017101053012734537d0*t**2*vapor%rho**3
      else if((t.ge.332.d0).and.(t.le.347.d0)) then
        vapor%dhdp = -2685.651624312814d0 + 40.408089706028505d0*t - 0.2433658728309133d0*t**2 +                                      &
                    0.0007330697580489127d0*t**3 - 1.1041140047375149d-6*t**4 + 6.650974354118622d-10*t**5 -                       &
                    13687.143431442617d0*vapor%rho + 161.11583016560508d0*t*vapor%rho - 0.7111991512467304d0*t**2*vapor%rho +                     &
                    0.0013952479260520127d0*t**3*vapor%rho-1.0264259793980518d-6*t**4*vapor%rho-6823.122132832979d0*vapor%rho**2+                 & 
                    59.14067793196761d0*t*vapor%rho**2 - 0.17092804592023148d0*t**2*vapor%rho**2 +                                           &
                    0.00016472384887903586d0*t**3*vapor%rho**2 - 253.26035577631035d0*vapor%rho**3 +                                         &
                    1.4874538670237314d0*t*vapor%rho**3 - 0.002182372555770298d0*t**2*vapor%rho**3 -                                         &
                    36.98478579604281d0*vapor%rho**4 + 0.10488084893954559d0*t*vapor%rho**4
      else if((t.gt.347.d0).and.(t.le.397.d0)) then
        vapor%dhdp = -13.7891645369465d0 + 0.101162496581722d0*t - 0.000251834570684114d0*t**2 + 2.11403618349564d-07*t**3 -          &
                   3.48932771355333d0*vapor%rho + 0.0177668288455497d0*t*vapor%rho - 0.0000226806539879019d0*t**2*vapor%rho -                     &
                   0.213028973710571d0*vapor%rho**2 + 0.000527329802151844d0*t*vapor%rho**2
      else if((t.gt.397.d0).and.(t.le.447.d0)) then
        vapor%dhdp = -2.30175009407832d0 + 0.0120585839248053d0*t - 0.0000212164847528242d0*t**2 + 1.22563572857418d-08*t**3 -        &
                   6.12857017700577d0*vapor%rho + 0.0433971794477965d0*t*vapor%rho - 0.000102787728553932d0*t**2*vapor%rho +                      &
                   8.13319224771949d-08*t**3*vapor%rho - 0.556896790690874d0*vapor%rho**2 + 0.00261694430608853d0*t*vapor%rho**2 -                &
                   3.06968998459052d-06*t**2*vapor%rho**2 - 0.0124048270971314d0*vapor%rho**3 + 0.0000281589576227285d0*t*vapor%rho**3 -          &
                   0.0000388702957448261*vapor%rho**4
      else if((t.gt.447.d0).and.(t.le.497.d0)) then
        vapor%dhdp = -3.53105327478957d0 + 0.0205069948762115d0*t - 0.0000406757484955379d0*t**2 + 2.72722604602877d-08*t**3 +        &
                   2.37467611198915d0*vapor%rho - 0.0156394857378528d0*t*vapor%rho + 0.0000341112896898702d0*t**2*vapor%rho -                     &
                   2.46756737993989d-08*t**3*vapor%rho - 1.09371605112861d0*vapor%rho**2 + 0.00709746784239992d0*t*vapor%rho**2 -                 &
                   0.0000153269087339451d0*t**2*vapor%rho**2 + 1.10142491313411d-08*t**3*vapor%rho**2 - 0.0204302542017631d0*vapor%rho**3 +       &
                   0.0000855009278132256d0*t*vapor%rho**3 - 8.93387442246377d-08*t**2*vapor%rho**3 - 0.0000765698781615174d0*vapor%rho**4 +       &
                   1.54780324969024d-07*t*vapor%rho**4-2.89819007617669d-08*vapor%rho**5
        
      else if((t.gt.497.d0).and.(t.le.517.d0)) then
    
        vapor%dhdp = -4.97816126755044d0 + 0.028377508774646d0*t - 0.0000547211830793747d0*t**2 + 3.54570677399166d-08*t**3 +         &
                   2.86357345093025d0*vapor%rho - 0.0174229907386162d0*t*vapor%rho + 0.0000352485407932868d0*t**2*vapor%rho -                     &
                   2.37258046341918d-08*t**3*vapor%rho - 0.423789001925988d0*vapor%rho**2 + 0.00259334071331075d0*t*vapor%rho**2 -                &
                   5.28716718390658d-06*t**2*vapor%rho**2 + 3.59102899347635d-09*t**3*vapor%rho**2 - 0.00393075761935776d0*vapor%rho**3 +         &
                   0.0000157503582818176d0*t*vapor%rho**3 - 1.5771273186733d-08*t**2*vapor%rho**3 - 0.0000109960071167684d0*vapor%rho**4 +        &
                   2.16567800282055d-08*t*vapor%rho**4 - 7.12630336161619d-09*vapor%rho**5
        
      else if((t.gt.517.d0).and.(t.le.537.d0)) then
        vapor%dhdp = -4.22428894926394d0 + 0.0231923752936463d0*t - 0.0000430905939952741d0*t**2 + 2.69041685271792d-08*t**3 +        &
                   1.72798881908645d0*vapor%rho - 0.0101426340594348d0*t*vapor%rho + 0.0000197878211815207d0*t**2*vapor%rho -                     &
                   1.2840856893509d-08*t**3*vapor%rho - 0.184734743817539d0*vapor%rho**2 + 0.00108817934315952d0*t*vapor%rho**2 -                 &
                   2.13593136453837d-06*t**2*vapor%rho**2 + 1.39692920624806d-09*t**3*vapor%rho**2 - 0.00124849540237478d0*vapor%rho**3 +         &
                   4.82126022692123d-06*t*vapor%rho**3 - 4.65321504153324d-09 *t**2*vapor%rho**3 - 2.61806209683763d-06*vapor%rho**4 +            &
                   4.97932038812081d-09*t*vapor%rho**4 - 1.37762262774754d-09*vapor%rho**5
      else !if((t.gt.537.d0).and.(t.le.552.d0)) then
        vapor%dhdp = -0.276039525042864d0 + 0.000744183648679953d0*t - 5.31879113542593d-07*t**2 - 0.0878128596698185d0*vapor%rho +        &
                   0.000310472439579764d0*t*vapor%rho - 2.76518971688037d-07*t**2*vapor%rho + 0.00793839664309268d0*vapor%rho**2 -                &
                   0.0000298325227021426d0*t*vapor%rho**2 + 2.79885337100082d-08*t**2*vapor%rho**2 - 0.000424804774961827d0*vapor%rho**3 +        &
                   1.59897629495813d-06*t*vapor%rho**3 - 1.50456834308976d-09*t**2*vapor%rho**3 - 8.57195694781748d-07*vapor%rho**4 +             &
                   1.59041253984922d-09*t*vapor%rho**4 - 4.2109909723164d-10 *vapor%rho**5
        
      end if
    end subroutine h2o_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97_l(eos,p,t,liquid)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid
      real(8) :: pi,tau,gt,gp,gpp,gtt,gpt
      real(8) :: g,g1,g2,g3,g4,g5,pi1,tau1,gp1
      integer :: k

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
      
      liquid%rho  = 25.841246054191803727*tau*gp1
      liquid%h    = 639675.036d0*gt
      liquid%drdp = -gpp*1.5632937721834122d-6*tau*gp1**2
      liquid%drdt = -liquid%rho*7.215007215007215d-4*tau + 0.01864447767257706*gpt*tau**3*gp1**2
      liquid%dhdp = 0.03869782431941923775*gpt
      liquid%dhdt = -461.526d0*tau**2*gtt
      
    end subroutine iapws97_l
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97_v(eos,p,t,vapor)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: vapor
      real(8) :: pww
      real(8) :: pi,tau,g0t,g0p,g0pp,g0tt,g0pt,grt,grp,grpp,grtt,grpt
      real(8) :: g0,g01,g02
      real(8) :: gr,gr1,gr2,gr3,gr4,gr5,pi1,tau1,tau2,gp1
      integer :: k
      
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
      
      vapor%rho  = 4.01245401527075799d0*tau*gp1
      vapor%h    = 249224.04d0*(g0t+grt)
      vapor%drdp = -(g0pp+grpp)*4.01245401527075799d-6*tau*gp1**2
      vapor%drdt = -vapor%rho*0.00185185185185185185*tau + 0.0074304703986495518*(g0pt+grpt)*tau**3*gp1**2
      vapor%dhdp = 0.24922404d0*(g0pt+grpt)
      vapor%dhdt = -461.526d0*tau**2*(g0tt+grtt)
      
    end subroutine iapws97_v
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine stiffened(eos,p,t,liquid)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid
      real(8), parameter :: r1 = 1.d0/2691.d0
      real(8) :: tau
      
      tau = 1.d0/t

      liquid%rho  = (p+8.5d8)*r1*tau
      liquid%h    = 4186.d0*t
      liquid%drdp = r1*tau
      liquid%drdt = -(p+8.5d8)*r1*tau**2
      liquid%dhdp = 0.d0
      liquid%dhdt = 4186.d0
      
    end subroutine stiffened
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine n2_l_fit(eos,p,t,liquid)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid
      real(8) :: rhol

      rhol = 3278.7479320710095d0 -140.53866400032032d0*t + 3.4391996294042255d0*t**2 - 0.04313918674235664d0*t**3+                &

             0.0002701512499721705d0*t**4 - 6.808567943865866d-7*t**5 + 0.0001263733642858279d0*p -                                &
             6.922008184428884d-6*t*p +1.4177982228082052d-7*t**2*p - 1.27888626669828d-9*t**3*p +                                 &
             4.337249290759174d-12*t**4*p + 5.963833700360859d-12*p**2 -                                                           &
             2.288595626069022d-13*t*p**2 +2.9357774279657827d-15*t**2*p**2 -                                                      &
             1.2665399523355606d-17*t**3*p**2  
      liquid%rho = rhol             
      liquid%h = 4.168304457241471d7 - 1.1042002477173102d6*t + 11052.049281367019d0*t**2 - 52.23483148623568d0*t**3 +                     &
            0.11846481736761054d0*t**4 - 0.00010643417382698117d0*t**5 - 110824.88074946059d0*rhol +                               &
            2598.250999713307d0*t*rhol - 21.36040545710846d0*t**2*rhol + 0.07224212430642923d0*t**3*rhol -                         &
            0.00008569087822637289d0*t**4*rhol + 95.28876405538512d0*rhol**2 - 1.9625537317149264d0*t*rhol**2 +                    &
            0.012463055910504148d0*t**2*rhol**2 - 0.000023537822079719658d0*t**3*rhol**2 - 0.02654233996233127d0*rhol**3 +         &
            0.00047324415207110907d0*t*rhol**3 - 2.0239668871820323d-6*t**2*rhol**3
      liquid%dhdt = 108819.94371947719d0 - 727.3413290413118d0*t + 1.269400197297674d0*t**2 - 189.0073325742833d0*rhol+                  &
                0.6388683793616025d0*t*rhol + 0.08373418130165292d0*rhol**2

      if(t.lt.75.d0) then
        liquid%drdp = 0.00020259945009638386d0 - 8.557700251478422d-7*t + 9.866811451399918d-10*t**2 -                                   &
                    3.8236563833158784d-7*liquid%rho + 8.16705219276515d-10*t*liquid%rho + 1.809590831406454d-10*liquid%rho**2
        liquid%drdt = -247.86002574696704 + 1.461427245574839*t - 0.002193220588925308*t**2 + 0.4220428895914477*liquid%rho -                  &
                    0.0012658720516194725*t*liquid%rho - 0.00018145431872668984*liquid%rho**2
        liquid%dhdp = -0.028313020074434314d0 + 0.00011804310572403998d0*t - 1.3397056490623622d-7*t**2 +                              &
                    0.000056694982040525745d0*rhol - 1.1453228268267254d-7*t*rhol - 2.751528916063015d-8*rhol**2

      else if((t.ge.75.d0).and.(t.lt.85.d0)) then
        liquid%drdp = 0.0003107714792534878d0-1.4580031144969927d-6*t+1.878326893643818d-9*t**2 -                                        &
                    5.913408276569343d-7*liquid%rho + 1.3886778703445793d-9*t*liquid%rho + 2.823329207920792d-10*liquid%rho**2
        liquid%drdt = -297.385479939352d0 + 1.7627825395290624d0*t - 0.00264790794680055d0*t**2 +                                        &
                    0.5153374297778337d0*liquid%rho -0.0015504587543153292d0*t*liquid%rho - 0.0002253532129203019d0*liquid%rho**2
        liquid%dhdp = -0.04672291067295588d0 + 0.0002145536800198554d0*t - 2.6721340030384323d-7*t**2 +                                &
                    0.00009282878073359741d0*rhol - 2.0802911493195356d-7*t*rhol - 4.530039947064013d-8*rhol**2

      else if((t.ge.85.d0).and.(t.lt.95.d0)) then
        liquid%drdp = 0.0005381464209109027d0-2.844760393759367d-6*t + 4.078531373621749d-9*t**2 -                                       &
                    1.027004966405762d-6*liquid%rho + 2.6986099993842745d-9*t*liquid%rho + 4.920323909333415d-10*liquid%rho**2
        liquid%drdt = -407.9014901431543d0 + 2.5256977364739486d0*t - 0.0039453749684370005d0*t**2 +                                     &
                    0.7170907625038794d0*liquid%rho - 0.002251246844406133d0*t*liquid%rho - 0.0003171779048991788d0*liquid%rho**2
        liquid%dhdp = -0.08594702942195079d0 + 0.0004496906307035824d0*t - 6.386575730157394d-7*t**2 +                                 &
                    0.00016853194316053486d0*rhol - 4.3084759764382327d-7*t*rhol - 8.204710498526094d-8*rhol**2

      else
        liquid%drdp = 0.0011931950341395573d0 - 7.632317529105058d-6*t + 1.3195922863204093d-8*t**2 -                                    &
                    2.2091813560313297d-6*liquid%rho + 6.920022440956056d-9*t*liquid%rho + 1.0319826829467008d-9*liquid%rho**2
        liquid%drdt = -664.661596560617d0 + 4.5231316387311296d0*t - 0.00775764679239673d0*t**2 +                                        &
                    1.164397373989908d0*liquid%rho - 0.004010605364921429d0*t*liquid%rho - 0.0005106823684308862d0*liquid%rho**2
        liquid%dhdp = -0.8447822231307364d0 + 0.005724722263525644d0*t - 9.854988248711547d-6*t**2 +                                   &
                    0.002619758754575089d0*rhol - 0.00001264121304115033d0*t*rhol + 1.2918556348269813d-8*t**2*rhol -              &
                    2.6615532549020144d-6*rhol**2 + 6.751573958788486d-9*t*rhol**2 + 8.893807771614017d-10*rhol**3

      end if
      liquid%cv = 0.935d0

    end subroutine n2_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine n2_v_fit(eos,p,t,vapor)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: vapor
      real(8) :: rhov

      if(t.lt.75.d0) then
        rhov  = 0.5166d0 - 0.02249d0*t + 0.0001501d0*p + 0.000325d0*t**2 - 0.000002214d0*t*p +                                     &
                0.0000000001156d0*p**2 - 0.000001559d0*t**3 + 0.00000001081d0*t**2*p -                                             &
                0.000000000001246d0*t*p**2
        vapor%drdp = 0.0002002d0 - 0.000004451d0*t + 0.00001142d0*rhov + 0.0000000439d0*t**2 - 0.0000002484d0*t*rhov                    &
                + 0.0000001744d0*rhov**2 - 0.0000000001621d0*t**3 + 0.000000001439d0*t**2*rhov                                     &
                - 0.000000002071d0*t*rhov**2 + 0.0000000005365d0*rhov**3
        vapor%drdt = 0.001042d0 - 0.00002954d0*t - 0.04749d0*rhov + 0.000000209d0*t**2 + 0.0007431d0*t*rhov-                            &
                0.002029d0*rhov**2 - 0.000003837d0*t**2*rhov + 0.00002203d0*t*rhov**2 - 0.00001267d0*rhov**3         
        vapor%dhdp = - 0.2506d0 + 0.005177d0*t - 0.006404d0*rhov - 0.00002963d0*t**2 + 0.00007627d0*t*rhov - 0.00003413d0*rhov**2

      else if((t.ge.75.d0).and.(t.lt.85.d0)) then
        rhov  = - 0.3101d0+0.00782d0*t+0.0001439d0*p - 0.00004921d0*t**2 - 0.000002023d0*t*p +                                     &
                0.00000000008892d0*p**2  + 0.000000009375d0*t**2*p - 0.0000000000009134d0*t*p**2 +                                 &
                1.234d-17*p**3
        vapor%drdp = 0.0001269d0 - 0.000001592d0*t + 0.000007968d0*rhov + 0.00000000665d0*t**2 - 0.000000154d0*t*rhov +                 &
                0.00000009901d0*rhov**2 + 0.0000000007947d0*t**2*rhov - 0.000000001059d0*t*rhov**2 +                               &
                0.000000000257d0*rhov**3
        vapor%drdt = - 0.01857d0 + 0.0008489d0*t - 0.0497d0*rhov - 0.00001259d0*t**2 + 0.0007831d0*t*rhov - 0.001607d0*rhov**2          &
                + 0.00000006098d0*t**3 - 0.000003973d0*t**2*rhov + 0.00001609d0*t*rhov**2 - 0.000007701d0*rhov**3
        vapor%dhdp = - 0.3547d0 + 0.009983d0*t - 0.01229d0*rhov - 0.0001031d0*t**2 + 0.0002587d0*t*rhov -                             &
                    0.00008049d0*rhov**2 + 0.0000003721d0*t**3 - 0.000001397d0*t**2*rhov + 0.0000008312d0*t*rhov**2

      else if((t.ge.85.d0).and.(t.lt.95.d0)) then
        rhov  = - 2.669d0 + 0.06128d0*t + 0.0001518d0*p - 0.0003514d0*t**2 - 0.00000217d0*t*p +                                    &
                0.00000000007329d0*p**2 + 0.000000009994d0*t**2*p - 0.0000000000007113d0*t*p**2 +                                  &
                7.852d-18*p**3
        vapor%drdt = - 0.003125d0 + 0.0000724d0*t - 0.02908d0*rhov - 0.0000004186d0*t**2 + 0.0002712d0*t*rhov -                         &
                0.00351d0*rhov**2 - 0.0000007957d0*t**2*rhov + 0.00006725d0*t*rhov**2 - 0.0000391d0*rhov**3 -                      &
                0.0000003415d0*t**2*rhov**2 + 0.0000003997d0*t*rhov**3 - 0.00000009137d0*rhov**4
        vapor%drdp = 0.000134d0 - 0.000001751d0*t + 0.000002705d0*rhov + 0.000000007533d0*t**2 -                                        &
                0.00000002436d0*t*rhov + 0.000000009563d0*rhov**2
        vapor%dhdp = - 0.276d0 + 0.006999d0*t - 0.007838d0*rhov - 0.00006547d0*t**2 + 0.0001495d0*t*rhov -                            &
                    0.00004313d0*rhov**2 + 0.0000002143d0*t**3 - 0.0000007297d0*t**2*rhov + 0.0000003972d0*t*rhov**2

      else 
        rhov  = 1.034d0 - 0.07884d0*t + 0.0001568d0*p + 0.001336d0*t**2 - 0.000002235d0*t*p +                                      &
                0.0000000000651d0*p**2 - 0.000006516d0*t**3 + 0.00000001007d0*t**2*p -                                             &
                0.0000000000006035d0*t*p**2 + 5.537d-18*p**3
        vapor%drdt = - 0.001017d0 + 0.00001877d0*t - 0.03033d0*rhov - 0.00000008603d0*t**2 + 0.0003087d0*t*rhov -                       &
                0.001221d0*rhov**2 - 0.000001054d0*t**2*rhov + 0.00001559d0*t*rhov**2 - 0.00009703d0*rhov**3 -                     &
                0.00000005037d0*t**2*rhov**2 + 0.000001885d0*t*rhov**3 - 0.0000005862d0*rhov**4 -                                  &
                0.000000009349d0*t**2*rhov**3 + 0.00000000576d0*t*rhov**4 - 0.000000000703d0*rhov**5 
        vapor%drdp = 0.0001016d0 - 0.000001022d0*t + 0.000003396d0*rhov + 0.000000003427d0*t**2 - 0.00000004841d0*t*rhov +              &
                0.0000001622d0*rhov**2 + 0.0000000001831d0*t**2*rhov - 0.000000003039d0*t*rhov**2 +                                &
                0.0000000009155d0*rhov**3 +  0.00000000001459d0*t**2*rhov**2 - 0.000000000008902d0*t*rhov**3 +                     &
                0.000000000001154d0*rhov**4
        vapor%dhdp = - 0.08365d0 + 0.001028d0*t - 0.006387d0*rhov - 0.000003573d0*t**2 + 0.0001179d0*t*rhov -                         &
                    0.00004007d0*rhov**2 - 0.0000005597d0*t**2*rhov + 0.0000003848d0*t*rhov**2 -                                   &
                    0.00000005636d0*rhov**3

      end if

      vapor%rho = rhov
      vapor%h = -503.2d0 + 1041.d0*t - 996.9d0*rhov + 5.045d0*t*rhov - 0.1174d0*rhov**2
      vapor%dhdt = 1230.d0 - 7.078d0*t + 125.2d0*rhov + 0.08643d0*t**2 - 2.239d0*t*rhov + 1.453d0*rhov**2 - 0.000348d0*t**3             &
                + 0.0108d0*t**2*rhov - 0.01384d0*t*rhov**2 + 0.002898d0*rhov**3


    end subroutine n2_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine o2_l_fit(eos,p,t,liquid)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid

      if(t.lt.70.d0) then
        liquid%rho  = 1501.d0 - 2.926d0*t - 8.138d-8*p - 0.01211d0*t**2  + 2.167d-8*t*p - 4.511d-15*p**2
        liquid%drdp = 0.00008776915721929953d0 - 2.2342524348669326d-7*t + 1.8141443035116113d-10*t**2 - 1.1754769451884666d-7*liquid%rho +       &
                      1.5048472539207272d-10*t*liquid%rho + 3.973954268541276d-11*liquid%rho**2                                                   
        liquid%drdt = -169.9573359197483d0 + 0.47922054939184056d0*t + 0.0005039707930294704d0*t**2 + 0.2164405353757457d0*liquid%rho -           &
                      0.00038108000586276937d0*t*liquid%rho - 0.00006886510035061562d0*liquid%rho**2                                              
        liquid%dhdp = -0.00730870470152617d0 + 3.8347472507992255d-6*t + 5.200986111633859d-8*t**2 +                                              &
                      0.000011968196827258478d0*liquid%rho - 8.536610900563244d-9*t*liquid%rho - 4.366271010728139d-9*liquid%rho**2               
        liquid%dhdt = 47052.378307313425d0 - 279.3877042029916d0*t + 0.37778990552899994d0*t**2 - 55.69792372925869d0*liquid%rho +                &
                      0.1775162796540966d0*t*liquid%rho + 0.016899439713916296d0*liquid%rho**2                                                    
      else if((t.ge.70.d0).and.(t.lt.85.d0)) then                                                                                                 
        liquid%rho  = 1518.d0 - 3.465d0*t - 9.49d-7*p - 0.007947d0*t**2 + 0.0000000342d0*t*p - 0.000000000000009018d0*p**2                        
        liquid%drdp = 0.00016901241652512783d0 - 5.961402721697739d-7*t + 6.239467083109434d-10*t**2 - 2.2779828363949149d-7*liquid%rho +         &
                      4.0160879711071303d-10*t*liquid%rho + 7.719525323215468d-11*liquid%rho**2                                                   
        liquid%drdt = -247.03936232218197d0 + 1.1179364828310698d0*t - 0.0013471142999006913d0*t**2 + 0.3050897707058698d0*liquid%rho -           &
                      0.0006893860625919242d0*t*liquid%rho - 0.00009600058090309282d0*liquid%rho**2                                               
        liquid%dhdp = -0.018644466497095738d0 + 0.00006609888609863794d0*t - 6.44935741312229d-8*t**2 +                                           &
                      0.000026780418705627872d0*liquid%rho - 4.571303726965722d-8*t*liquid%rho - 9.303966708048961d-9*liquid%rho**2               
        liquid%dhdt = 35585.60979044708d0 - 151.58475414505295d0*t + 0.1749340495961691d0*t**2 - 44.43367072628075d0*liquid%rho +                 &
                      0.09726527126267216d0*t*liquid%rho + 0.01463184914003078d0*liquid%rho**2                                                    
      else if((t.ge.85.d0).and.(t.lt.100.d0)) then                                                                                                
        liquid%rho  = 1462.d0 - 2.119d0*t - 0.000002796d0*p - 0.01602d0*t**2 + 0.00000005613d0*t*p                                                &
                      - 0.00000000000001847d0*p**2                                                                                                
        liquid%drdp = 0.0015512976090951017d0 - 7.3255697316179285d-6*t + 8.988677750931408d-9*t**2 - 3.309605354370924d-6*liquid%rho +           &
                      1.0894923700822548d-8*t*liquid%rho - 7.492435840777116d-12*t**2*liquid%rho + 2.3452488471876264d-9*liquid%rho**2 -          &
                      4.002938350210652d-12*t*liquid%rho**2 - 5.520668106679057d-13*liquid%rho**3                                                 
        liquid%drdt = -335.9704458503552d0 + 1.6190358088963328d0*t - 0.0020589065325389362d0*t**2 + 0.42083182086310705d0*liquid%rho -           &
                      0.0010150989461190684d0*t*liquid%rho - 0.00013365894511512824d0*liquid%rho**2                                               
        liquid%dhdp = -0.03716579503864091d0 + 0.00015778688708237198d0*t - 1.909804586556369d-7*t**2 +                                           &
                      0.000051827247269057006d0*liquid%rho - 1.0584792429617636d-7*t*liquid%rho - 1.7838567910501354d-8*liquid%rho**2             
        liquid%dhdt = 39335.27694086211d0 - 166.66275558328675d0*t + 0.19599865974318578d0*t**2 - 49.728784845470486d0*liquid%rho +               &
                      0.10722163369619506d0*t*liquid%rho + 0.016520946417112048d0*liquid%rho**2                                                   
      else                                                                                                                                        
        liquid%rho  = 1317.d0 + 0.8469d0*t - 0.000007636d0*p  - 0.03114d0*t**2 + 0.0000001051d0*t*p - 0.00000000000004438d0*p**2                  
        liquid%drdp = 0.0006670590554266884d0 - 3.2916420186119954d-6*t +  4.47420442330802d-9*t**2 - 8.888294773549778d-7*liquid%rho +           &
                      2.1561506434117898d-9*t*liquid%rho + 2.980042793744716d-10*liquid%rho**2                                                    
        liquid%drdt = -496.0107160500555d0 + 2.5628145689358757d0*t - 0.003390485799877293d0*t**2 + 0.6271693191024084d0*liquid%rho -             &
                      0.0016345119848304563d0*t*liquid%rho - 0.00019965225365294045d0*liquid%rho**2                                               
        liquid%dhdp = -0.08046060269748674d0 + 0.0003991073945728477d0*t - 5.483787602922069d-7*t**2 +                                            &
                      0.00010899415071437275d0*liquid%rho - 2.612833026845765d-7*t*liquid%rho - 3.688790604999396d-8*liquid%rho**2                
        liquid%dhdt = 60126.07415718087d0 - 289.79951116169536d0*t + 0.3696713613029599d0*t**2 - 76.4887566760527d0*liquid%rho +                  &
                      0.1880417883602597d0*t*liquid%rho + 0.025059718783963285d0*liquid%rho**2     
      end if
      
      if (t.lt.65.d0) then     
        liquid%h = 188344.28513978166d0 + 33189.920824275454d0*t - 158.71274463422736d0*t**2 - 0.1971006940102005d0*t**3-                         &
                   1825.3244113713697d0*liquid%rho - 32.32874078304518d0*t*liquid%rho + 0.14488139179892803d0*t**2*liquid%rho +                   &
                   1.061222115988811d0*liquid%rho**2 + 0.0067357264348732685d0*t*liquid%rho**2                                                          
      else if((t.ge.65.d0).and.(t.lt.75.d0)) then                                                                                                 
        liquid%h = 2.693820010849507d6-17842.324834953455d0*t+14.860712553743063d0*t**2 + 0.021985362402576526d0*t**3-                            &
                   5159.614157640176d0*liquid%rho + 30.40327343109178d0*t*liquid%rho - 0.02413596888059153d0*t**2*liquid%rho +                    &
                   2.0997897334718143d0*liquid%rho**2 - 0.009425307668614939d0*t*liquid%rho**2                                                          
      else if((t.ge.75.d0).and.(t.lt.85.d0)) then                                                                                                 
        liquid%h = 3.2677256466973964d6 - 28881.16763356495d0*t + 61.41339331937471d0*t**2 - 0.001902770015740886d0*t**3-                         &
                   5897.999164193667d0*liquid%rho + 42.809773910643116d0*t*liquid%rho - 0.057739210010128696d0*t**2*liquid%rho +                  &
                   2.32683825318106d0*liquid%rho**2 - 0.012465065325212965d0*t*liquid%rho**2                                                            
      else if((t.ge.85.d0).and.(t.lt.95.d0)) then                                                                                                 
        liquid%h = 1.062884607310819d6 + 595.0398216397284d0*t - 7.057940616995017d0*t**2 - 2610.003238100057d0*liquid%rho +                      &
                   3.1249835040989016d0*t*liquid%rho + 1.1244431312261893d0*liquid%rho**2                                                               
      else                                                                                                                                        
        liquid%h = 3.9863701370553784d6 - 41832.792234712186d0*t + 125.32397983730809d0*t**2 -                                                    &
                   0.10315557591569982d0*t**3 - 6761.866796744369d0*liquid%rho + 55.347333056645034d0*t*liquid%rho -                              &
                   0.08957261017260913d0*t**2*liquid%rho + 2.5904940126370963d0*liquid%rho**2 - 0.015432565994687501d0*t*liquid%rho**2
      end if   
      liquid%cv = 0.935d0
    end subroutine o2_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine o2_v_fit(eos,p,t,vapor)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: vapor
     
      if(t.lt.70.d0) then
        vapor%rho  = 0.09839d0 - 0.002916d0*t + 0.0001142d0*p + 0.00002159d0*t**2 - 0.0000008471d0*t*p + 0.0000000000452d0*p**2
        vapor%drdp = 0.0007707348924580752d0 - 0.000041958503188899d0*t + 1.0448016877818711d-6*t**2 - 1.3484224840948894d-8*t**3 +               &
                     8.811905122709702d-11*t**4 - 2.313579674114593d-13*t**5 + 0.00142912831370235d0*vapor%rho -                                  &
                     0.0000726312592119253d0*t*vapor%rho + 1.383068888346053d-6*t**2*vapor%rho - 1.1686878601551056d-8*t**3*vapor%rho +           &
                     3.696798836829787d-11*t**4*vapor%rho + 0.0003750168074693673d0*vapor%rho**2 - 0.000013553811902321636d0*t*vapor%rho**2 +     &
                     1.6337841688517373d-7*t**2*vapor%rho**2 - 6.567433738348047d-10*t**3*vapor%rho**2
        vapor%drdt = 0.014906237909174724d0 - 0.0006838002599396501d0*t + 0.000010434222091002844d0*t**2 -                                        &
                     5.2969519762551634d-8*t**3 - 0.24062214550091177d0*vapor%rho + 0.010104368078388902d0*t*vapor%rho -                          &
                     0.0001536639947518677d0*t**2*vapor%rho + 7.930110953343673d-7*t**3*vapor%rho - 2.222024436521593d0*vapor%rho**2 +            &
                     0.09127920840035736d0*t*vapor%rho**2 - 0.001242149652136588d0*t**2*vapor%rho**2 + 5.594359313170392d-6*t**3*vapor%rho**2 -   &
                     0.7255066544058085d0*vapor%rho**3 + 0.018444541277570482d0*t*vapor%rho**3 - 0.00011601069430877057d0*t**2*vapor%rho**3 +     &
                     0.08053602556802765d0*vapor%rho**4 - 0.0011710426076103977d0*t*vapor%rho**4 + 0.0011449298831236335d0*vapor%rho**5
        vapor%dhdp = -5099.172498105065d0 + 380.528083883449d0*t - 11.35334176169662d0*t**2 + 0.1692958799114193d0*t**3 -                         &
                     0.0012617517364290234d0*t**4 + 3.760191059655844d-6*t**5 - 9827.159407957586d0*vapor%rho +                                   &
                     560.2763183150105d0*t*vapor%rho - 11.998275105009657d0*t**2*vapor%rho + 0.11436575036862953d0*t**3*vapor%rho -               &
                     0.00040934334971020687d0*t**4*vapor%rho - 182.06655056471268d0*vapor%rho**2 + 7.842132808398672d0*t*vapor%rho**2 -           &
                     0.11254015090812083d0*t**2*vapor%rho**2 + 0.0005381105463875461d0*t**3*vapor%rho**2 + 6.642745261937848d0*vapor%rho**3 -     &
                     0.20393126775621448d0*t*vapor%rho**3 + 0.0015547916446733549d0*t**2*vapor%rho**3
        vapor%dhdt = 563309.187628491d0 - 45731.317018270536d0*t + 1475.414069624006d0*t**2 - 23.637030947078387d0*t**3 +                         &
                     0.1882222367295729d0*t**4 - 0.0005964480354189977d0*t**5 + 5.116730037432892d6*vapor%rho -                                   &
                     303394.85512502067d0*t*vapor%rho + 6743.570827164565d0*t**2*vapor%rho - 66.59537385825031d0*t**3*vapor%rho +                 &
                     0.2465460467461267d0*t**4*vapor%rho + 3.9001650223352215d6*vapor%rho**2 - 165926.51516490726d0*t*vapor%rho**2 +              &
                     2357.039442948709d0*t**2*vapor%rho**2 - 11.178055933876037d0*t**3*vapor%rho**2 + 87173.77547421321d0*vapor%rho**3 -          &
                     2517.441352677458d0*t*vapor%rho**3 + 18.173234911213374d0*t**2*vapor%rho**3
      else if((t.ge.70.d0).and.(t.lt.85.d0)) then
        vapor%rho  = - 0.0004756d0 + 0.00001172d0*t + 0.000152d0*p - 0.00000007195d0*t**2 - 0.000001994d0*t*p                                     &
                     + 0.0000000001147d0*p**2 + 0.000000008694d0*t**2*p - 0.000000000001166d0*t*p**2 + 2.484d-17*p**3
        vapor%drdp = 0.0007707348924580752d0 - 0.000041958503188899d0*t + 1.0448016877818711d-6*t**2 - 1.3484224840948894d-8*t**3 +               &
                     8.811905122709702d-11*t**4 - 2.313579674114593d-13*t**5 + 0.00142912831370235d0*vapor%rho -                                  &
                     0.0000726312592119253d0*t*vapor%rho + 1.383068888346053d-6*t**2*vapor%rho - 1.1686878601551056d-8*t**3*vapor%rho +           &
                     3.696798836829787d-11*t**4*vapor%rho + 0.0003750168074693673d0*vapor%rho**2 - 0.000013553811902321636d0*t*vapor%rho**2 +     &
                     1.6337841688517373d-7*t**2*vapor%rho**2 - 6.567433738348047d-10*t**3*vapor%rho**2
        if ((t.ge.70.d0).and.(t.lt.76.d0)) then              
          vapor%drdt = 0.014906237909174724d0 - 0.0006838002599396501d0*t + 0.000010434222091002844d0*t**2 -                                      &
                       5.2969519762551634d-8*t**3 - 0.24062214550091177d0*vapor%rho + 0.010104368078388902d0*t*vapor%rho -                        &
                       0.0001536639947518677d0*t**2*vapor%rho + 7.930110953343673d-7*t**3*vapor%rho - 2.222024436521593d0*vapor%rho**2 +          &
                       0.09127920840035736d0*t*vapor%rho**2 - 0.001242149652136588d0*t**2*vapor%rho**2 + 5.594359313170392d-6*t**3*vapor%rho**2 - &
                       0.7255066544058085d0*vapor%rho**3 + 0.018444541277570482d0*t*vapor%rho**3 - 0.00011601069430877057d0*t**2*vapor%rho**3 +   &
                       0.08053602556802765d0*vapor%rho**4 - 0.0011710426076103977d0*t*vapor%rho**4 + 0.0011449298831236335d0*vapor%rho**5
          vapor%dhdp = -1820.3042938684684d0 + 124.79671006046507d0*t - 3.4205261714246307d0*t**2 + 0.046854912176703865d0*t**3-                  &
                       0.0003207867772214351d0*t**4 + 8.781864124784989d-7*t**5 - 1492.3555940657948d0*vapor%rho +                                &
                       78.16533675291102d0*t*vapor%rho - 1.53795291043518d0*t**2*vapor%rho + 0.01347029975329776d0*t**3*vapor%rho -               &
                       0.00004430685521187763d0*t**4*vapor%rho - 14.816104935984022d0*vapor%rho**2 + 0.6114246391170685d0*t*vapor%rho**2 -        &
                       0.008393988617545237d0*t**2*vapor%rho**2 + 0.00003834723318112214d0*t**3*vapor%rho**2 -                                    &
                       1.5348756113957471d0*vapor%rho**3 + 0.040272778836831336d0*t*vapor%rho**3 - 0.0002643750611815428d0*t**2*vapor%rho**3
          vapor%dhdt = 183273.5505750412d0 - 16025.112210677551d0*t + 535.4218239747133d0*t**2 - 8.647297438749607d0*t**3 +                       &
                       0.06815998250509395d0*t**4 - 0.00021102564140562663d0*t**5 + 2.1961869005671577d6*vapor%rho -                              &
                       120184.92414268502d0*t*vapor%rho + 2465.0309735757846d0*t**2*vapor%rho - 22.459508787390995d0*t**3*vapor%rho +             &
                       0.07670397495727141d0*t**4*vapor%rho + 660931.8075247495d0*vapor%rho**2 - 25868.53525710509d0*t*vapor%rho**2 +             &
                       338.0260411082621d0*t**2*vapor%rho**2 - 1.47444111627998d0*t**3*vapor%rho**2
        else if ((t.ge.76.d0).and.(t.lt.80.d0)) then       
          vapor%drdt = -0.04854331619304875d0 + 0.0018625637021031302d0*t - 0.000023817754142181705d0*t**2 +                                      &
                       1.0150815139525493d-7*t**3 + 0.40671541285969615d0*vapor%rho - 0.01658087479962605d0*t*vapor%rho +                         &
                       0.00021614163163704268d0*t**2*vapor%rho - 9.297911681877266d-7*t**3*vapor%rho - 1.027532307027257d0*vapor%rho**2 +         &
                       0.03949440839012544d0*t*vapor%rho**2 - 0.0005059457667093034d0*t**2*vapor%rho**2 +                                         &
                       2.1598348204722054d-6*t**3*vapor%rho**2 - 0.08354601271375685d0*vapor%rho**3 + 0.0020508124290226633d0*t*vapor%rho**3 -    &
                       0.000012607879342385479d0*t**2*vapor%rho**3 - 0.00021418885346537008d0*vapor%rho**4 +                                      &
                       2.8334581494197435d-6*t*vapor%rho**4 - 4.452886490115064d-6*vapor%rho**5
          vapor%dhdp = -163.57238613423138d0 + 8.423138131758563d0*t - 0.16254082860573565d0*t**2 +0.0013930576743292702d0*t**3-                  &
                       4.474580654124217d-6*t**4 - 53.178274889731014d0*vapor%rho + 1.957035440351718d0*t*vapor%rho -                             &
                       0.024051878786737528d0*t**2*vapor%rho + 0.00009869921174019753d0*t**3*vapor%rho - 0.13162960174507596d0*vapor%rho**2 +     &
                       0.0036178674197744325d0*t*vapor%rho**2 - 0.00002452880070910895d0*t**2*vapor%rho**2 -                                      &
                       0.012757780185688423d0*vapor%rho**3 + 0.0001556785951708372d0*t*vapor%rho**3                           
          vapor%dhdt = -206027.57432064012d0 + 10608.053734733474d0*t - 203.89896331947034d0*t**2 + 1.7416608969007021d0*t**3 -                   &
                       0.005578218830657945d0*t**4 + 987949.0438127739d0*vapor%rho - 50718.46526891996d0*t*vapor%rho +                            &
                       976.0046623476202d0*t**2*vapor%rho - 8.344697638500415d0*t**3*vapor%rho + 0.026747362042505118d0*t**4*vapor%rho +          &
                       156934.16976084435d0*vapor%rho**2 - 5773.509161557668d0*t*vapor%rho**2 + 70.93328911382892d0*t**2*vapor%rho**2 -           &
                       0.2909909179777248d0*t**3*vapor%rho**2 + 267.06540513392736d0*vapor%rho**3 - 7.375877673824645d0*t*vapor%rho**3 +          &
                       0.05020991727009204d0*t**2*vapor%rho**3 + 21.143045855507914d0*vapor%rho**4 - 0.2582839625268897d0*t*vapor%rho**4 +        &
                       0.01336769169446307d0*vapor%rho**5
        else
          vapor%drdt = 0.47974996138307924d0 - 0.024128849557415104d0*t + 0.0004543447918707363d0*t**2 -                                          &
                       3.7965661987300736d-6*t**3 + 1.1879741844521969d-8*t**4 - 1.7281404709285162d0*vapor%rho +                                 &
                       0.08720801502889175d0*t*vapor%rho - 0.0016648686751448687d0*t**2*vapor%rho + 0.000014123944835723762d0*t**3*vapor%rho -    &
                       4.485575916280877d-8*t**4*vapor%rho - 0.5029230161955298d0*vapor%rho**2 + 0.01827637234889779d0*t*vapor%rho**2 -           &
                       0.00022145318352130195d0*t**2*vapor%rho**2 + 8.94325375023116d-7*t**3*vapor%rho**2 - 0.02327833210996124d0*vapor%rho**3 +  &
                       0.0005424584096473847d0*t*vapor%rho**3 - 3.1665945399323005d-6*t**2*vapor%rho**3 - 0.00009058365937300747d0*vapor%rho**4 + &
                       1.1221225199509703d-6*t*vapor%rho**4 - 1.0746303747724054d-6*vapor%rho**5
          vapor%dhdp = -79.80330163853974d0 + 3.88757274267661d0*t - 0.07096908305532168d0*t**2 + 0.0005753414460866079d0*t**3 -                  &
                        1.7478310939415295d-6*t**4 - 15.295717171466867d0*vapor%rho + 0.5322982088854236d0*t*vapor%rho -                          &
                        0.006186682492240336d0*t**2*vapor%rho + 0.00002401032541190659d0*t**3*vapor%rho - 0.041804840416864913d0*vapor%rho**2 +   &
                        0.0011029569461511418d0*t*vapor%rho**2 - 7.149382002954615d-6*t**2*vapor%rho**2 - 0.0030109382295691813d0*vapor%rho**3+   &
                        0.000034570238053543616d0*t*vapor%rho**3 
          vapor%dhdt = 19897.845131710368d0 - 3177.986531591686d0*t + 126.09104424300202d0*t**2 - 2.121485612233824d0*t**3 +                      &
                       0.016441862971140482d0*t**4 - 0.00004851982372421387d0*t**5 + 481774.4839590495d0*vapor%rho -                              &
                       23379.789559965284d0*t*vapor%rho + 425.25286439212846d0*t**2*vapor%rho - 3.436269653408993d0*t**3*vapor%rho +              &
                       0.010408912589833242d0*t**4*vapor%rho + 44793.85318530364d0*vapor%rho**2 - 1558.2148713242539d0*t*vapor%rho**2 +           &
                       18.102863688070727d0*t**2*vapor%rho**2 - 0.07022730563376371d0*t**3*vapor%rho**2 + 94.01585127406705d0*vapor%rho**3 -      &
                       2.4470092032328252d0*t*vapor%rho**3 + 0.015694693462109874d0*t**2*vapor%rho**3 + 4.762111745531126d0*vapor%rho**4 -        &
                       0.05455320489374141d0*t*vapor%rho**4
        end if
      else if((t.ge.85.d0).and.(t.lt.100.d0)) then
        vapor%rho  = - 0.02582d0 + 0.0005927d0*t + 0.00013d0*p - 0.000003392d0*t**2 - 0.000001467d0*t*p                                           &
                     + 0.00000000008644d0*p**2+ 0.000000005532d0*t**2*p - 0.000000000000824d0*t*p**2 + 1.541d-17*p**3
        vapor%drdp = 0.000166895070887691d0 - 2.7150419042070515d-6*t + 1.9640040289688417d-8*t**2 - 5.331483083429449d-11*t**3 +                 &
                     5.644462745348359d-6*vapor%rho - 9.625145699824075d-8*t*vapor%rho + 4.469609834034141d-10*t**2*vapor%rho +                   &
                     3.3442983102396117d-7*vapor%rho**2 - 3.5901782914117336d-9*t*vapor%rho**2
        vapor%drdt = 0.0003597d0 - 0.00001197d0*t - 0.04671d0*vapor%rho + 0.0000001349d0*t**2 + 0.000822d0*t*vapor%rho -                          &
                     0.02823d0*vapor%rho**2 - 0.0000000005156d0*t**3 - 0.000006473d0*t**2*vapor%rho + 0.00006079d0*t*vapor%rho**2 -               &
                     0.0007027d0*vapor%rho**3 + 0.00000001929d0*t**3*vapor%rho - 0.0000003541d0*t**2*vapor%rho**2 +                               &
                     0.000007856d0*t*vapor%rho**3 - 0.000001556d0*vapor%rho**4
        vapor%dhdt = 2203.6112035650203d0 - 90.14473230568774d0*t + 2.390003485202703d0*t**2 - 0.030626759354033382d0*t**3 +                      &
                     0.0001914496462308313d0*t**4 - 4.697594073248713d-7*t**5 + 11342.80115819394d0*vapor%rho -                                   &
                     512.18589347155d0*t*vapor%rho + 8.597236776654182d0*t**2*vapor%rho - 0.06366942705073778d0*t**3*vapor%rho +                  &
                     0.00017578120806203588d0*t**4*vapor%rho + 4932.327616002995d0*vapor%rho**2 - 152.6602458978181d0*t*vapor%rho**2 +            &
                     1.576024450302989d0*t**2*vapor%rho**2 - 0.005426446413328457d0*t**3*vapor%rho**2 - 21.885639287056527d0*vapor%rho**3 +       &
                     0.4062498339942897d0*t*vapor%rho**3 - 0.0018331972154759458d0*t**2*vapor%rho**3 + 2.09504673362577d0*vapor%rho**4 -          &
                     0.023144795503526796d0*t*vapor%rho**4
        if (t.lt.91.d0) then    
          vapor%dhdp = 0.13000955507009873d0 - 0.0033513769515030258d0*t + 0.000019003813364007332d0*t**2 -                                       &
                       0.471368548100311d0*vapor%rho + 0.010212305523918915d0*t*vapor%rho - 0.00005543220186946809d0*t**2*vapor%rho +             &
                       0.028987751443963147d0*vapor%rho**2 - 0.0006213565521429516d0*t*vapor%rho**2 + 3.3383066141924628d-6*t**2*vapor%rho**2 -   &
                       0.0102120169716408d0*vapor%rho**3 + 0.00022256948265323033d0*t*vapor%rho**3 - 1.2144262025908764d-6*t**2*vapor%rho**3 -    &
                       0.000030039333916901116d0*vapor%rho**4 + 3.365652964801301d-7*t*vapor%rho**4 - 3.230151429203042d-8*vapor%rho**5
        else
          vapor%dhdp = 0.3247381399523399d0 - 0.010521663341191409d0*t + 0.00010619158851164091d0*t**2 -                                          &
                        3.506350858865523d-7*t**3 - 0.6602971091356331d0*vapor%rho + 0.019777545834480666d0*t*vapor%rho -                         &
                        0.00019791065907962655d0*t**2*vapor%rho + 6.61326614616879d-7*t**3*vapor%rho + 0.039409670704404626d0*vapor%rho**2 -      &
                        0.0011535794188021811d0*t*vapor%rho**2 + 0.000011284954798273413d0*t**2*vapor%rho**2 -                                    &
                        3.6886295058841984d-8*t**3*vapor%rho**2 - 0.0012008331972785885d0*vapor%rho**3 +                                          &
                        0.00002353474242400194d0*t*vapor%rho**3 - 1.1536369937566051d-7*t**2*vapor%rho**3 +                                       &
                        0.00001810207935028912d0*vapor%rho**4 - 1.8680258753120824d-7*t*vapor%rho**4 + 2.1166235190874365d-7*vapor%rho**5           
        end if 
      else
        vapor%rho  = 0.009632253324945106d0 - 0.00017573530285693376d0*t + 8.168869369962596d-7*t**2 + 0.00010777377953337692d0*p -               &
                     1.004547904883942d-6*t*p + 3.1157965295030853d-9*t**2*p + 2.9499227174247904d-11*p**2 -                                      &
                     2.1138410719274732d-13*t*p**2 + 1.5327208040241761d-18*p**3
        vapor%drdp = 0.0001803274568602044d0 - 3.3795391580549686d-6*t + 3.165897582126288d-8*t**2 - 1.4822850694723483d-10*t**3 +                &
                     2.774873306903886d-13*t**4 + 5.928950205356179d-6*vapor%rho - 1.1824844073292441d-7*t*vapor%rho +                            &
                     8.48477866401482d-10*t**2*vapor%rho - 2.1313160698610714d-12*t**3*vapor%rho + 9.333978024089662d-8*vapor%rho**2 -            &
                     1.4579445504397311d-9*t*vapor%rho**2 + 5.8282608584087095d-12*t**2*vapor%rho**2 + 3.222235852124598d-11*vapor%rho**3 -       &
                     1.3709718330677488d-13*t*vapor%rho**3
        vapor%drdt = -0.00002515d0 + 0.000000475d0*t - 0.02816d0*vapor%rho - 0.000000002236d0*t**2 + 0.0002639d0*t*vapor%rho -                    &
                     0.0004721d0*vapor%rho**2 - 0.000000823d0*t**2*vapor%rho + 0.000003115d0*t*vapor%rho**2 - 0.000001188d0*vapor%rho**3
        vapor%dhdp = 0.1759340420608009d0 - 0.007073771329378833d0*t + 0.00009218190230974735d0*t**2 -                                            &
                     5.099758487698818d-7*t**3 + 1.0345126472325508d-9*t**4 - 0.047240246573340444d0*vapor%rho +                                  &
                     0.0012229392313518933d0*t*vapor%rho - 0.000010635833919253268d0*t**2*vapor%rho + 3.0987308214774d-8*t**3*vapor%rho +         &
                     0.0005780905361482492d0*vapor%rho**2 - 0.000010005520490858424d0*t*vapor%rho**2 + 4.3323122926162755d-8*t**2*vapor%rho**2 -  &
                     1.0971249591802358d-6*vapor%rho**3 + 8.90315214224214d-9*t*vapor%rho**3
        vapor%dhdt = 859.5361003250456d0 + 1.4008644968401225d0*t - 0.012977767339376965d0*t**2 + 0.00003996054935106315d0*t**3 -                 &
                     164.20620494516544d0*vapor%rho + 4.12333889917936d0*t*vapor%rho - 0.032794402039010524d0*t**2*vapor%rho +                    &
                     0.00008506491677537112d0*t**3*vapor%rho + 13.509572745523892d0*vapor%rho**2 - 0.23281284728100154d0*t*vapor%rho**2 +         &
                     0.0010086105642078795d0*t**2*vapor%rho**2 - 0.030859911460888973d0*vapor%rho**3 + 0.0002471441196742016d0*t*vapor%rho**3 +   &
                     0.00012443689479352496d0*vapor%rho**4
      end if
      
      if (t.lt.65.d0) then     
        vapor%h = -668305.9182671085d0 + 52776.26639215887d0*t - 1594.8682285165742d0*t**2 + 24.31457194925873d0*t**3 -                           &
                  0.18401275450645885d0*t**4 + 0.0005535811365063484d0*t**5 - 5.465336679119291d6*vapor%rho +                                     &
                  307082.8631464915d0*t*vapor%rho - 6461.838613084738d0*t**2*vapor%rho + 60.363287751797884d0*t**3*vapor%rho -                    &
                  0.21124800496224108d0*t**4*vapor%rho - 3.1235172523588953d6*vapor%rho**2 + 124781.93304209826d0*t*vapor%rho**2 -                &
                  1664.2754203931222d0*t**2*vapor%rho**2 + 7.409539038225778d0*t**3*vapor%rho**2 + 5005.064772759815d0*vapor%rho**3 -             &
                  130.36528663262695d0*t*vapor%rho**3 + 0.8505776739383123d0*t**2*vapor%rho**3      
      else if((t.ge.65.d0).and.(t.lt.75.d0)) then        
        vapor%h = -668305.9182671085d0 + 52776.26639215887d0*t - 1594.8682285165742d0*t**2 + 24.31457194925873d0*t**3 -                           &
                  0.18401275450645885d0*t**4 + 0.0005535811365063484d0*t**5 - 5.465336679119291d6*vapor%rho +                                     &
                  307082.8631464915d0*t*vapor%rho - 6461.838613084738d0*t**2*vapor%rho + 60.363287751797884d0*t**3*vapor%rho -                    &
                  0.21124800496224108d0*t**4*vapor%rho - 3.1235172523588953d6*vapor%rho**2 + 124781.93304209826d0*t*vapor%rho**2 -                &
                  1664.2754203931222d0*t**2*vapor%rho**2 + 7.409539038225778d0*t**3*vapor%rho**2 + 5005.064772759815d0*vapor%rho**3 -             &
                  130.36528663262695d0*t*vapor%rho**3 + 0.8505776739383123d0*t**2*vapor%rho**3                     
      else if((t.ge.75.d0).and.(t.lt.85.d0)) then        
        vapor%h = 602607.0025649129d0 - 31434.64179498057d0*t + 674.9372956355897d0*t**2 - 6.766693571017301d0*t**3 +                             &
                  0.03188135094788187d0*t**4 - 0.000053797851534581746d0*t**5 - 1.9047616532657095d6*vapor%rho +                                  &
                  95510.17988726991d0*t*vapor%rho - 1793.3691527397891d0*t**2*vapor%rho + 14.944774059958775d0*t**3*vapor%rho -                   &
                  0.046644295226138235d0*t**4*vapor%rho - 282530.8553560337d0*vapor%rho**2 + 10070.164811021295d0*t*vapor%rho**2 -                &
                  119.81887920061938d0*t**2*vapor%rho**2 + 0.47583657123133827d0*t**3*vapor%rho**2 - 383.0601068638107d0*vapor%rho**3 +           &
                  9.733063155726166d0*t*vapor%rho**3 - 0.06142616242241004d0*t**2*vapor%rho**3       
      else if((t.ge.85.d0).and.(t.lt.95.d0)) then            
        vapor%h = -17335.847493337365d0 + 5590.0913077739715d0*t - 186.84942731521483d0*t**2 + 2.9933095250332964d0*t**3 -                        &
                  0.02169878512751981d0*t**4 + 0.00005940881417254163d0*t**5 - 458516.68404563505d0*vapor%rho +                                   &
                  20481.007671163396d0*t*vapor%rho - 342.78114918942157d0*t**2*vapor%rho + 2.5457651437802755d0*t**3*vapor%rho -                  &
                  0.0070800660231476265d0*t**4*vapor%rho - 21998.44897558312d0*vapor%rho**2 + 699.217785643861d0*t*vapor%rho**2 -                 &
                  7.422756054880894d0*t**2*vapor%rho**2 + 0.02631154307075599d0*t**3*vapor%rho**2 - 36.32333433750597d0*vapor%rho**3 +            &
                  0.9038822763369369d0*t*vapor%rho**3 - 0.005446312090879567d0*t**2*vapor%rho**3 - 1.2978404707087177d0*vapor%rho**4 +            &
                  0.013298294549614897d0*t*vapor%rho**4       
      else  
        vapor%h = -213927.8590721971d0 + 12325.48231787402d0*t - 241.06129329094887d0*t**2 + 2.5192572314707062d0*t**3 -                          &
                  0.013048433283000527d0*t**4 + 0.000026826088799777146d0*t**5 - 51745.19300574396d0*vapor%rho +                                  &
                  2003.0066953336134d0*t*vapor%rho - 29.25954922775434d0*t**2*vapor%rho + 0.18947160576149089d0*t**3*vapor%rho -                  &
                  0.0004585329660949181d0*t**4*vapor%rho - 1176.6344464267477d0*vapor%rho**2 + 33.030639959726315d0*t*vapor%rho**2 -              &
                  0.30893946169331915d0*t**2*vapor%rho**2 + 0.0009623949949494379d0*t**3*vapor%rho**2 - 7.452446482887528d0*vapor%rho**3 +        &
                  0.13844669394649298d0*t*vapor%rho**3 - 0.0006392895175965527d0*t**2*vapor%rho**3 - 0.024998007187044494d0*vapor%rho**4 +        &
                  0.0002130764320100439d0*t*vapor%rho**4        
      end if       

    end subroutine o2_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2_l_fit(eos,p,t,liquid)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: liquid
      liquid = t_eos2(0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
    end subroutine h2_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2_v_fit(eos,p,t,vapor)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: vapor
      vapor = t_eos2(0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
    end subroutine h2_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine he_v_fit(eos,p,t,gas)
      implicit none
      class(t_eos), intent(in) :: eos
      real(8), intent(in) :: p,t
      type(t_eos2), intent(out) :: gas
      gas = t_eos2(0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
    end subroutine he_v_fit
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
      pww = 0.d0
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
