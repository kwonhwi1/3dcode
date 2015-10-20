module bc_module
  use mpi
  use config_module
  use grid_module
  use variable_module
  use eos_module
  use prop_module
  implicit none
  private
  public :: t_bc
  
  type, extends(t_bcinfo) :: t_bcinfo2
    integer :: origin(3),dir(3)
    integer :: neighbor1(3),neighbor2(3),neighbor3(3)
    integer :: neighbor4(3),neighbor5(3),neighbor6(3)
    character(4) :: face
    real(8) :: massflowrate,total_pressure,total_temperature,omega(3)
    real(8), dimension(:), allocatable :: pv,tv,dv
    procedure(p_bctype), pointer :: bctype
  end type t_bcinfo2

  type t_prec
    real(8) :: str,uref
    procedure(p_getsndp2), pointer :: getsndp2
  end type t_prec

  type t_bc
    private
    integer :: rank,size,iturb,nbc,ncon
    integer :: npv,ntv,ndv,ngrd
    class(t_bcinfo2), dimension(:), allocatable :: bcinfo,corner,edge
    class(t_connectinfo), dimension(:), allocatable :: connectinfo
    class(t_mpitemp), dimension(:), allocatable :: mpitemp
    type(t_prec) :: prec
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setbc
  end type t_bc

  interface
    function p_getsndp2(prec,snd2,uuu2) result(sndp2)
      import t_prec
      implicit none
      class(t_prec), intent(in) :: prec
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
    end function p_getsndp2

    subroutine p_bctype(bcinfo,grid,variable,eos,prop,prec)
      import t_bcinfo2
      import t_grid
      import t_variable
      import t_eos
      import t_prop
      import t_prec
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
    end subroutine p_bctype
  end interface
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(bc,config,grid,variable,eos,prop)
      implicit none
      class(t_bc), intent(out) :: bc
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: l,n,m
      integer, parameter :: dim = 3
      logical :: isurf_edge(12),jsurf_edge(12),ksurf_edge(12)
      logical :: isurf_corner(8),jsurf_corner(8),ksurf_corner(8)

      bc%rank = config%getrank() 
      bc%size = config%getsize()
      bc%iturb = config%getiturb()
      
      bc%nbc = grid%getnbc()
      bc%ncon = grid%getncon()
      bc%npv = variable%getnpv()
      bc%ntv = variable%getntv()
      bc%ndv = variable%getndv()
      bc%ngrd = grid%getngrd()

      bc%prec%uref = config%geturef()
      bc%prec%str  = config%getstr()
      select case(config%getprec())
      case(0)
        bc%prec%getsndp2 => no_prec
      case(1)
        bc%prec%getsndp2 => steady_prec
      case(2)
        bc%prec%getsndp2 => unsteady_prec
      case default
        bc%prec%getsndp2 => null()
      end select

      allocate(bc%bcinfo(bc%nbc),bc%connectinfo(bc%ncon))

      do n=1,bc%nbc
        bc%bcinfo(n)%bcname = grid%getbcname(n)
        do m=1,dim
          bc%bcinfo(n)%istart(m) = grid%getbcistart(n,m)
          bc%bcinfo(n)%iend(m)   = grid%getbciend(n,m)
        end do
        
        if(bc%bcinfo(n)%istart(1).eq.bc%bcinfo(n)%iend(1)) then
          if(bc%bcinfo(n)%istart(1).eq.1) then
            bc%bcinfo(n)%face = 'imin'
            bc%bcinfo(n)%istart(1) = bc%bcinfo(n)%istart(1)-2
            bc%bcinfo(n)%origin(1) = 2*bc%bcinfo(n)%iend(1)+1
            bc%bcinfo(n)%origin(2) = 0
            bc%bcinfo(n)%origin(3) = 0
            bc%bcinfo(n)%dir(1) = -1
            bc%bcinfo(n)%dir(2) =  1
            bc%bcinfo(n)%dir(3) =  1
          else
            bc%bcinfo(n)%face = 'imax'
            bc%bcinfo(n)%iend(1) = bc%bcinfo(n)%iend(1)+2
            bc%bcinfo(n)%origin(1) = 2*bc%bcinfo(n)%istart(1)-1
            bc%bcinfo(n)%origin(2) = 0
            bc%bcinfo(n)%origin(3) = 0
            bc%bcinfo(n)%dir(1) = -1
            bc%bcinfo(n)%dir(2) =  1
            bc%bcinfo(n)%dir(3) =  1
          end if
        else if(bc%bcinfo(n)%istart(2).eq.bc%bcinfo(n)%iend(2)) then
          if(bc%bcinfo(n)%istart(2).eq.1) then
            bc%bcinfo(n)%face = 'jmin'
            bc%bcinfo(n)%istart(2) = bc%bcinfo(n)%istart(2)-2
            bc%bcinfo(n)%origin(1) = 0
            bc%bcinfo(n)%origin(2) = 2*bc%bcinfo(n)%iend(2)+1
            bc%bcinfo(n)%origin(3) = 0
            bc%bcinfo(n)%dir(1) =  1
            bc%bcinfo(n)%dir(2) = -1
            bc%bcinfo(n)%dir(3) =  1
         else
            bc%bcinfo(n)%face = 'jmax'
            bc%bcinfo(n)%iend(2) = bc%bcinfo(n)%iend(2)+2
            bc%bcinfo(n)%origin(1) = 0
            bc%bcinfo(n)%origin(2) = 2*bc%bcinfo(n)%istart(2)-1
            bc%bcinfo(n)%origin(3) = 0
            bc%bcinfo(n)%dir(1) =  1
            bc%bcinfo(n)%dir(2) = -1
            bc%bcinfo(n)%dir(3) =  1         
          end if
        else if(bc%bcinfo(n)%istart(3).eq.bc%bcinfo(n)%iend(3)) then
          if(bc%bcinfo(n)%istart(3).eq.1) then
            bc%bcinfo(n)%face = 'kmin'
            bc%bcinfo(n)%istart(3) = bc%bcinfo(n)%istart(3)-2
            bc%bcinfo(n)%origin(1) = 0
            bc%bcinfo(n)%origin(2) = 0
            bc%bcinfo(n)%origin(3) = 2*bc%bcinfo(n)%iend(3)+1
            bc%bcinfo(n)%dir(1) =  1
            bc%bcinfo(n)%dir(2) =  1
            bc%bcinfo(n)%dir(3) = -1         
          else
            bc%bcinfo(n)%face = 'kmax'
            bc%bcinfo(n)%iend(3) = bc%bcinfo(n)%iend(3)+2
            bc%bcinfo(n)%origin(1) = 0
            bc%bcinfo(n)%origin(2) = 0
            bc%bcinfo(n)%origin(3) =  2*bc%bcinfo(n)%istart(3)-1
            bc%bcinfo(n)%dir(1) =  1
            bc%bcinfo(n)%dir(2) =  1
            bc%bcinfo(n)%dir(3) = -1         
          end if 
        end if
        
        if(trim(bc%bcinfo(n)%bcname).eq.'BCWall') then 
          select case(config%getiturb())
          case(-1)
            bc%bcinfo(n)%bctype => bcwallviscouske
          case(0)
            bc%bcinfo(n)%bctype => bcwallviscouskw
          case(-2)
            bc%bcinfo(n)%bctype => bcwallviscous
          case default
            bc%bcinfo(n)%bctype => bcwallinviscid
          end select
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCInflow') then
          bc%bcinfo(n)%bctype => bcinflow
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCInflowSubsonic') then
          bc%bcinfo(n)%bctype => bcinflowsubsonic
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCInflowSupersonic') then
          bc%bcinfo(n)%bctype => bcinflowsupersonic
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCOutflow') then
          bc%bcinfo(n)%bctype => bcoutflow
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCOutflowSubsonic') then
          bc%bcinfo(n)%bctype => bcoutflowsubsonic
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCOutflowSupersonic') then
          bc%bcinfo(n)%bctype => bcoutflowsupersonic
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCExtrapolate') then
          bc%bcinfo(n)%bctype => bcoutflowsupersonic
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCFarfield') then
          bc%bcinfo(n)%bctype => bcfarfield
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCDegeneratePoint') then
          bc%bcinfo(n)%bctype => bcdegeneratepoint
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCDegenerateLine') then
          bc%bcinfo(n)%bctype => bcdegeneratepoint
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCSymmetryPlane') then
          bc%bcinfo(n)%bctype => bcsymmetryplane
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCShiftedPeriodic') then
          bc%bcinfo(n)%bctype => bcshiftedperiodic
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCTotalPressureIn') then
          bc%bcinfo(n)%total_temperature = 0.d0
          bc%bcinfo(n)%total_pressure = 0.d0
          bc%bcinfo(n)%bctype => bctotalpressurein
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCCounterRotatingWall') then
          select case(config%getiturb())
          case(-2,-1,0)
            bc%bcinfo(n)%bctype => bccounterrotatingwallviscous
          case default
            bc%bcinfo(n)%bctype => bccounterrotatingwallinviscid
          end select
        else
          bc%bcinfo(n)%bctype => null()
          write(*,*) 'error, check bc name',bc%bcinfo(n)%bcname
        end if
        bc%bcinfo(n)%neighbor1 = 0 ! meaningless
        bc%bcinfo(n)%neighbor2 = 0 ! meaningless
        bc%bcinfo(n)%neighbor3 = 0 ! meaningless
        bc%bcinfo(n)%neighbor4 = 0 ! meaningless
        bc%bcinfo(n)%neighbor5 = 0 ! meaningless
        bc%bcinfo(n)%neighbor6 = 0 ! meaningless

        allocate(bc%bcinfo(n)%pv(bc%npv),bc%bcinfo(n)%dv(bc%ndv),bc%bcinfo(n)%tv(bc%ntv))

        bc%bcinfo(n)%pv(1) = config%getpref()
        bc%bcinfo(n)%pv(2) = config%geturef()*dcos(config%getaos())*dcos(config%getaoa())
        bc%bcinfo(n)%pv(3) = config%geturef()*dcos(config%getaos())*dsin(config%getaoa())
        bc%bcinfo(n)%pv(4) = config%geturef()*dsin(config%getaos())
        bc%bcinfo(n)%pv(5) = config%gettref()
        bc%bcinfo(n)%pv(6) = config%gety1ref()
        bc%bcinfo(n)%pv(7) = config%gety2ref()

        call eos%deteos(bc%bcinfo(n)%pv(1),bc%bcinfo(n)%pv(5),bc%bcinfo(n)%pv(6),bc%bcinfo(n)%pv(7),bc%bcinfo(n)%dv)

        if(config%getiturb().ge.-2) then
          call prop%detprop(bc%bcinfo(n)%dv(3),bc%bcinfo(n)%dv(4),bc%bcinfo(n)%dv(5),bc%bcinfo(n)%pv(5) &
                           ,bc%bcinfo(n)%pv(6),bc%bcinfo(n)%pv(7),bc%bcinfo(n)%tv(1:2))
          if(config%getiturb().ge.-1) then
            bc%bcinfo(n)%tv(3) = config%getemutref()
            bc%bcinfo(n)%pv(8) = config%getkref()
            bc%bcinfo(n)%pv(9) = config%getoref()
          end if
        end if
      end do

      bc%bcinfo(n)%omega = config%getomega()

      do n=1,bc%ncon
        bc%connectinfo(n)%donor = grid%getconnectdonor(n)
        
        do m=1,dim
          bc%connectinfo(n)%istart(m) = grid%getconnectistart(n,m)
          bc%connectinfo(n)%iend(m)   = grid%getconnectiend(n,m)
          bc%connectinfo(n)%istart_donor(m) = grid%getconnectistart_donor(n,m)
          bc%connectinfo(n)%iend_donor(m)   = grid%getconnectiend_donor(n,m)
          do l=1,dim
            bc%connectinfo(n)%transmat(m,l) = grid%getconnecttransmat(n,m,l)
          end do
        end do

        if(bc%connectinfo(n)%istart(2).eq.bc%connectinfo(n)%iend(2)) then
          if(bc%connectinfo(n)%istart(2).eq.2) then
            bc%connectinfo(n)%iend(2) = bc%connectinfo(n)%iend(2) + 2
            if(bc%connectinfo(n)%transmat(1,2).eq.1) then
              bc%connectinfo(n)%iend_donor(1) = bc%connectinfo(n)%iend_donor(1) + 2 
            else if(bc%connectinfo(n)%transmat(1,2).eq.-1) then
              bc%connectinfo(n)%iend_donor(1) = bc%connectinfo(n)%iend_donor(1) - 2
            else if(bc%connectinfo(n)%transmat(2,2).eq.1) then
              bc%connectinfo(n)%iend_donor(2) = bc%connectinfo(n)%iend_donor(2) + 2 
            else if(bc%connectinfo(n)%transmat(2,2).eq.-1) then
              bc%connectinfo(n)%iend_donor(2) = bc%connectinfo(n)%iend_donor(2) - 2
            else if(bc%connectinfo(n)%transmat(3,2).eq.1) then
              bc%connectinfo(n)%iend_donor(3) = bc%connectinfo(n)%iend_donor(3) + 2 
            else if(bc%connectinfo(n)%transmat(3,2).eq.-1) then
              bc%connectinfo(n)%iend_donor(3) = bc%connectinfo(n)%iend_donor(3) - 2
            end if
          else
            bc%connectinfo(n)%istart(2) = bc%connectinfo(n)%istart(2) - 2
            if(bc%connectinfo(n)%transmat(1,2).eq.1) then
              bc%connectinfo(n)%istart_donor(1) = bc%connectinfo(n)%istart_donor(1) - 2 
            else if(bc%connectinfo(n)%transmat(1,2).eq.-1) then
              bc%connectinfo(n)%istart_donor(1) = bc%connectinfo(n)%istart_donor(1) + 2
            else if(bc%connectinfo(n)%transmat(2,2).eq.1) then
              bc%connectinfo(n)%istart_donor(2) = bc%connectinfo(n)%istart_donor(2) - 2
            else if(bc%connectinfo(n)%transmat(2,2).eq.-1) then
              bc%connectinfo(n)%istart_donor(2) = bc%connectinfo(n)%istart_donor(2) + 2
            else if(bc%connectinfo(n)%transmat(3,2).eq.1) then
              bc%connectinfo(n)%istart_donor(3) = bc%connectinfo(n)%istart_donor(3) - 2
            else if(bc%connectinfo(n)%transmat(3,2).eq.-1) then
              bc%connectinfo(n)%istart_donor(3) = bc%connectinfo(n)%istart_donor(3) + 2
            end if
          end if
        else if(bc%connectinfo(n)%istart(1).eq.bc%connectinfo(n)%iend(1)) then
          if(bc%connectinfo(n)%istart(1).eq.2) then
            bc%connectinfo(n)%iend(1) = bc%connectinfo(n)%iend(1) + 2
            if(bc%connectinfo(n)%transmat(1,1).eq.1) then
              bc%connectinfo(n)%iend_donor(1) = bc%connectinfo(n)%iend_donor(1) + 2 
            else if(bc%connectinfo(n)%transmat(1,1).eq.-1) then
              bc%connectinfo(n)%iend_donor(1) = bc%connectinfo(n)%iend_donor(1) - 2
            else if(bc%connectinfo(n)%transmat(2,1).eq.1) then
              bc%connectinfo(n)%iend_donor(2) = bc%connectinfo(n)%iend_donor(2) + 2 
            else if(bc%connectinfo(n)%transmat(2,1).eq.-1) then
              bc%connectinfo(n)%iend_donor(2) = bc%connectinfo(n)%iend_donor(2) - 2
            else if(bc%connectinfo(n)%transmat(3,1).eq.1) then
              bc%connectinfo(n)%iend_donor(3) = bc%connectinfo(n)%iend_donor(3) + 2 
            else if(bc%connectinfo(n)%transmat(3,1).eq.-1) then
              bc%connectinfo(n)%iend_donor(3) = bc%connectinfo(n)%iend_donor(3) - 2
            end if
          else
            bc%connectinfo(n)%istart(1) = bc%connectinfo(n)%istart(1) - 2
            if(bc%connectinfo(n)%transmat(1,1).eq.1) then
              bc%connectinfo(n)%istart_donor(1) = bc%connectinfo(n)%istart_donor(1) - 2 
            else if(bc%connectinfo(n)%transmat(1,1).eq.-1) then
              bc%connectinfo(n)%istart_donor(1) = bc%connectinfo(n)%istart_donor(1) + 2
            else if(bc%connectinfo(n)%transmat(2,1).eq.1) then
              bc%connectinfo(n)%istart_donor(2) = bc%connectinfo(n)%istart_donor(2) - 2
            else if(bc%connectinfo(n)%transmat(2,1).eq.-1) then
              bc%connectinfo(n)%istart_donor(2) = bc%connectinfo(n)%istart_donor(2) + 2
            else if(bc%connectinfo(n)%transmat(3,1).eq.1) then
              bc%connectinfo(n)%istart_donor(3) = bc%connectinfo(n)%istart_donor(3) - 2
            else if(bc%connectinfo(n)%transmat(3,1).eq.-1) then
              bc%connectinfo(n)%istart_donor(3) = bc%connectinfo(n)%istart_donor(3) + 2
            end if
          end if
        else if(bc%connectinfo(n)%istart(3).eq.bc%connectinfo(n)%iend(3)) then
          if(bc%connectinfo(n)%istart(3).eq.2) then
            bc%connectinfo(n)%iend(3) = bc%connectinfo(n)%iend(3) + 2
            if(bc%connectinfo(n)%transmat(1,3).eq.1) then
              bc%connectinfo(n)%iend_donor(1) = bc%connectinfo(n)%iend_donor(1) + 2 
            else if(bc%connectinfo(n)%transmat(1,3).eq.-1) then
              bc%connectinfo(n)%iend_donor(1) = bc%connectinfo(n)%iend_donor(1) - 2
            else if(bc%connectinfo(n)%transmat(2,3).eq.1) then
              bc%connectinfo(n)%iend_donor(2) = bc%connectinfo(n)%iend_donor(2) + 2 
            else if(bc%connectinfo(n)%transmat(2,3).eq.-1) then
              bc%connectinfo(n)%iend_donor(2) = bc%connectinfo(n)%iend_donor(2) - 2
            else if(bc%connectinfo(n)%transmat(3,3).eq.1) then
              bc%connectinfo(n)%iend_donor(3) = bc%connectinfo(n)%iend_donor(3) + 2 
            else if(bc%connectinfo(n)%transmat(3,3).eq.-1) then
              bc%connectinfo(n)%iend_donor(3) = bc%connectinfo(n)%iend_donor(3) - 2
            end if
          else
            bc%connectinfo(n)%istart(3) = bc%connectinfo(n)%istart(3) - 2
            if(bc%connectinfo(n)%transmat(1,3).eq.1) then
              bc%connectinfo(n)%istart_donor(1) = bc%connectinfo(n)%istart_donor(1) - 2 
            else if(bc%connectinfo(n)%transmat(1,3).eq.-1) then
              bc%connectinfo(n)%istart_donor(1) = bc%connectinfo(n)%istart_donor(1) + 2
            else if(bc%connectinfo(n)%transmat(2,3).eq.1) then
              bc%connectinfo(n)%istart_donor(2) = bc%connectinfo(n)%istart_donor(2) - 2
            else if(bc%connectinfo(n)%transmat(2,3).eq.-1) then
              bc%connectinfo(n)%istart_donor(2) = bc%connectinfo(n)%istart_donor(2) + 2
            else if(bc%connectinfo(n)%transmat(3,3).eq.1) then
              bc%connectinfo(n)%istart_donor(3) = bc%connectinfo(n)%istart_donor(3) - 2
            else if(bc%connectinfo(n)%transmat(3,3).eq.-1) then
              bc%connectinfo(n)%istart_donor(3) = bc%connectinfo(n)%istart_donor(3) + 2
            end if
          end if
        end if
      end do
      
      allocate(bc%mpitemp(bc%ncon))
      do n=1,bc%ncon
        bc%mpitemp(n)%num = (abs(bc%connectinfo(n)%istart(1)-bc%connectinfo(n)%iend(1))+1) &
                           *(abs(bc%connectinfo(n)%istart(2)-bc%connectinfo(n)%iend(2))+1) &
                           *(abs(bc%connectinfo(n)%istart(3)-bc%connectinfo(n)%iend(3))+1) & 
                           *(bc%npv+bc%ntv+bc%ndv)
        allocate(bc%mpitemp(n)%sendbuf(bc%mpitemp(n)%num),bc%mpitemp(n)%recvbuf(bc%mpitemp(n)%num))
      end do
      
      allocate(bc%edge(12),bc%corner(8))

      bc%edge(1)%istart(1)  = -1                  ; bc%edge(1)%iend(1)  = 1
      bc%edge(1)%istart(2)  = -1                  ; bc%edge(1)%iend(2)  = 1
      bc%edge(1)%istart(3)  =   2                 ; bc%edge(1)%iend(3)  = grid%getkmax()
      bc%edge(1)%neighbor1(1) = 1                 ; bc%edge(1)%neighbor2(1) = 2
      bc%edge(1)%neighbor1(2) = 2                 ; bc%edge(1)%neighbor2(2) = 1
      bc%edge(1)%neighbor1(3) = 0                 ; bc%edge(1)%neighbor2(3) = 0
      bc%edge(1)%neighbor3(1) = 0                 ; bc%edge(1)%origin(1) = 2
      bc%edge(1)%neighbor3(2) = 0                 ; bc%edge(1)%origin(2) = 2
      bc%edge(1)%neighbor3(3) = 0                 ; bc%edge(1)%origin(3) = 0

      bc%edge(2)%istart(1)  = grid%getimax()+1    ; bc%edge(2)%iend(1)  = grid%getimax()+3
      bc%edge(2)%istart(2)  = -1                  ; bc%edge(2)%iend(2)  = 1
      bc%edge(2)%istart(3)  =  2                  ; bc%edge(2)%iend(3)  = grid%getkmax()
      bc%edge(2)%neighbor1(1) = grid%getimax()+1  ; bc%edge(2)%neighbor2(1) = grid%getimax()
      bc%edge(2)%neighbor1(2) = 2                 ; bc%edge(2)%neighbor2(2) = 1
      bc%edge(2)%neighbor1(3) = 0                 ; bc%edge(2)%neighbor2(3) = 0
      bc%edge(2)%neighbor3(1) = 0                 ; bc%edge(2)%origin(1) = grid%getimax() 
      bc%edge(2)%neighbor3(2) = 0                 ; bc%edge(2)%origin(2) = 2
      bc%edge(2)%neighbor3(3) = 0                 ; bc%edge(2)%origin(3) = 0

      bc%edge(3)%istart(1)  = -1                  ; bc%edge(3)%iend(1)  = 1
      bc%edge(3)%istart(2)  = grid%getjmax()+1    ; bc%edge(3)%iend(2)  = grid%getjmax()+3
      bc%edge(3)%istart(3)  = 2                   ; bc%edge(3)%iend(3)  = grid%getkmax()
      bc%edge(3)%neighbor1(1) = 1                 ; bc%edge(3)%neighbor2(1) = 2
      bc%edge(3)%neighbor1(2) = grid%getjmax()    ; bc%edge(3)%neighbor2(2) = grid%getjmax()+1
      bc%edge(3)%neighbor1(3) = 0                 ; bc%edge(3)%neighbor2(3) = 0
      bc%edge(3)%neighbor3(1) = 0                 ; bc%edge(3)%origin(1) = 2
      bc%edge(3)%neighbor3(2) = 0                 ; bc%edge(3)%origin(2) = grid%getjmax()
      bc%edge(3)%neighbor3(3) = 0                 ; bc%edge(3)%origin(3) = 0

      bc%edge(4)%istart(1)  = grid%getimax()+1    ; bc%edge(4)%iend(1)  = grid%getimax()+3
      bc%edge(4)%istart(2)  = grid%getjmax()+1    ; bc%edge(4)%iend(2)  = grid%getjmax()+3
      bc%edge(4)%istart(3)  = 2                   ; bc%edge(4)%iend(3)  = grid%getkmax()
      bc%edge(4)%neighbor1(1) = grid%getimax()+1  ; bc%edge(4)%neighbor2(1) = grid%getimax()
      bc%edge(4)%neighbor1(2) = grid%getjmax()    ; bc%edge(4)%neighbor2(2) = grid%getjmax()+1
      bc%edge(4)%neighbor1(3) = 0                 ; bc%edge(4)%neighbor2(3) = 0
      bc%edge(4)%neighbor3(1) = 0                 ; bc%edge(4)%origin(1) = grid%getimax()
      bc%edge(4)%neighbor3(2) = 0                 ; bc%edge(4)%origin(2) = grid%getjmax()
      bc%edge(4)%neighbor3(3) = 0                 ; bc%edge(4)%origin(3) = 0
      
      bc%edge(5)%istart(1)  = -1                  ; bc%edge(5)%iend(1)  = 1
      bc%edge(5)%istart(2)  =  2                  ; bc%edge(5)%iend(2)  = grid%getjmax()
      bc%edge(5)%istart(3)  = -1                  ; bc%edge(5)%iend(3)  = 1
      bc%edge(5)%neighbor1(1) = 1                 ; bc%edge(5)%neighbor2(1) = 0
      bc%edge(5)%neighbor1(2) = 0                 ; bc%edge(5)%neighbor2(2) = 0
      bc%edge(5)%neighbor1(3) = 2                 ; bc%edge(5)%neighbor2(3) = 0
      bc%edge(5)%neighbor3(1) = 2                 ; bc%edge(5)%origin(1) = 2
      bc%edge(5)%neighbor3(2) = 0                 ; bc%edge(5)%origin(2) = 0
      bc%edge(5)%neighbor3(3) = 1                 ; bc%edge(5)%origin(3) = 2

      bc%edge(6)%istart(1)  = grid%getimax()+1    ; bc%edge(6)%iend(1)  = grid%getimax()+3
      bc%edge(6)%istart(2)  =  2                  ; bc%edge(6)%iend(2)  = grid%getjmax()
      bc%edge(6)%istart(3)  = -1                  ; bc%edge(6)%iend(3)  = 1
      bc%edge(6)%neighbor1(1) = grid%getimax()+1  ; bc%edge(6)%neighbor2(1) = 0
      bc%edge(6)%neighbor1(2) = 0                 ; bc%edge(6)%neighbor2(2) = 0
      bc%edge(6)%neighbor1(3) = 2                 ; bc%edge(6)%neighbor2(3) = 0
      bc%edge(6)%neighbor3(1) = grid%getimax()    ; bc%edge(6)%origin(1) = grid%getimax()
      bc%edge(6)%neighbor3(2) = 0                 ; bc%edge(6)%origin(2) = 0
      bc%edge(6)%neighbor3(3) = 1                 ; bc%edge(6)%origin(3) = 2

      bc%edge(7)%istart(1)  = -1                  ; bc%edge(7)%iend(1)  = 1
      bc%edge(7)%istart(2)  =  2                  ; bc%edge(7)%iend(2)  = grid%getjmax()
      bc%edge(7)%istart(3)  = grid%getkmax()+1    ; bc%edge(7)%iend(3)  = grid%getkmax()+3
      bc%edge(7)%neighbor1(1) = 1                 ; bc%edge(7)%neighbor2(1) = 0
      bc%edge(7)%neighbor1(2) = 0                 ; bc%edge(7)%neighbor2(2) = 0
      bc%edge(7)%neighbor1(3) = grid%getkmax()    ; bc%edge(7)%neighbor2(3) = 0
      bc%edge(7)%neighbor3(1) = 2                 ; bc%edge(7)%origin(1) = 2
      bc%edge(7)%neighbor3(2) = 0                 ; bc%edge(7)%origin(2) = 0
      bc%edge(7)%neighbor3(3) = grid%getkmax()+1  ; bc%edge(7)%origin(3) = grid%getkmax()

      bc%edge(8)%istart(1)  = grid%getimax()+1    ; bc%edge(8)%iend(1)  = grid%getimax()+3
      bc%edge(8)%istart(2)  = 2                   ; bc%edge(8)%iend(2)  = grid%getjmax()
      bc%edge(8)%istart(3)  = grid%getkmax()+1    ; bc%edge(8)%iend(3)  = grid%getkmax()+3  
      bc%edge(8)%neighbor1(1) = grid%getimax()+1  ; bc%edge(8)%neighbor2(1) = 0
      bc%edge(8)%neighbor1(2) = 0                 ; bc%edge(8)%neighbor2(2) = 0
      bc%edge(8)%neighbor1(3) = grid%getkmax()    ; bc%edge(8)%neighbor2(3) = 0
      bc%edge(8)%neighbor3(1) = grid%getimax()    ; bc%edge(8)%origin(1) = grid%getimax()
      bc%edge(8)%neighbor3(2) = 0                 ; bc%edge(8)%origin(2) = 0
      bc%edge(8)%neighbor3(3) = grid%getkmax()+1  ; bc%edge(8)%origin(3) = grid%getkmax()
      
      bc%edge(9)%istart(1)  =  2                  ; bc%edge(9)%iend(1)  = grid%getimax()
      bc%edge(9)%istart(2)  = -1                  ; bc%edge(9)%iend(2)  = 1
      bc%edge(9)%istart(3)  = -1                  ; bc%edge(9)%iend(3)  = 1
      bc%edge(9)%neighbor1(1) = 0                 ; bc%edge(9)%neighbor2(1) = 0
      bc%edge(9)%neighbor1(2) = 0                 ; bc%edge(9)%neighbor2(2) = 1
      bc%edge(9)%neighbor1(3) = 0                 ; bc%edge(9)%neighbor2(3) = 2
      bc%edge(9)%neighbor3(1) = 0                 ; bc%edge(9)%origin(1) = 0
      bc%edge(9)%neighbor3(2) = 2                 ; bc%edge(9)%origin(2) = 2
      bc%edge(9)%neighbor3(3) = 1                 ; bc%edge(9)%origin(3) = 2

      bc%edge(10)%istart(1) =  2                  ; bc%edge(10)%iend(1) = grid%getimax()
      bc%edge(10)%istart(2) = grid%getjmax()+1    ; bc%edge(10)%iend(2) = grid%getjmax()+3
      bc%edge(10)%istart(3) = -1                  ; bc%edge(10)%iend(3) = 1
      bc%edge(10)%neighbor1(1) = 0                ; bc%edge(10)%neighbor2(1) = 0
      bc%edge(10)%neighbor1(2) = 0                ; bc%edge(10)%neighbor2(2) = grid%getjmax()+1
      bc%edge(10)%neighbor1(3) = 0                ; bc%edge(10)%neighbor2(3) = 2
      bc%edge(10)%neighbor3(1) = 0                ; bc%edge(10)%origin(1) = 0
      bc%edge(10)%neighbor3(2) = grid%getjmax()   ; bc%edge(10)%origin(2) = grid%getjmax()
      bc%edge(10)%neighbor3(3) = 1                ; bc%edge(10)%origin(3) = 2

      bc%edge(11)%istart(1) =  2                  ; bc%edge(11)%iend(1) = grid%getimax()
      bc%edge(11)%istart(2) = -1                  ; bc%edge(11)%iend(2) = 1
      bc%edge(11)%istart(3) = grid%getkmax()+1    ; bc%edge(11)%iend(3) = grid%getkmax()+3
      bc%edge(11)%neighbor1(1) = 0                ; bc%edge(11)%neighbor2(1) = 0
      bc%edge(11)%neighbor1(2) = 0                ; bc%edge(11)%neighbor2(2) = 1
      bc%edge(11)%neighbor1(3) = 0                ; bc%edge(11)%neighbor2(3) = grid%getkmax()
      bc%edge(11)%neighbor3(1) = 0                ; bc%edge(11)%origin(1) = 0
      bc%edge(11)%neighbor3(2) = 2                ; bc%edge(11)%origin(2) = 2
      bc%edge(11)%neighbor3(3) = grid%getkmax()+1 ; bc%edge(11)%origin(3) = grid%getkmax()

      bc%edge(12)%istart(1) = 2                   ; bc%edge(12)%iend(1) = grid%getimax()
      bc%edge(12)%istart(2) = grid%getjmax()+1    ; bc%edge(12)%iend(2) = grid%getjmax()+3
      bc%edge(12)%istart(3) = grid%getkmax()+1    ; bc%edge(12)%iend(3) = grid%getkmax()+3      
      bc%edge(12)%neighbor1(1) = 0                ; bc%edge(12)%neighbor2(1) = 0
      bc%edge(12)%neighbor1(2) = 0                ; bc%edge(12)%neighbor2(2) = grid%getjmax()+1
      bc%edge(12)%neighbor1(3) = 0                ; bc%edge(12)%neighbor2(3) = grid%getkmax()
      bc%edge(12)%neighbor3(1) = 0                ; bc%edge(12)%origin(1) = 0
      bc%edge(12)%neighbor3(2) = grid%getjmax()   ; bc%edge(12)%origin(2) = grid%getjmax()
      bc%edge(12)%neighbor3(3) = grid%getkmax()+1 ; bc%edge(12)%origin(3) = grid%getkmax()
      
      do m=1,12
        bc%edge(m)%dir = 0
        allocate(bc%edge(m)%pv(bc%npv),bc%edge(m)%dv(bc%ndv),bc%edge(m)%tv(bc%ntv))

        bc%edge(m)%pv(1) = config%getpref()
        bc%edge(m)%pv(2) = config%geturef()*dcos(config%getaos())*dcos(config%getaoa())
        bc%edge(m)%pv(3) = config%geturef()*dcos(config%getaos())*dsin(config%getaoa())
        bc%edge(m)%pv(4) = config%geturef()*dsin(config%getaos())
        bc%edge(m)%pv(5) = config%gettref()
        bc%edge(m)%pv(6) = config%gety1ref()
        bc%edge(m)%pv(7) = config%gety2ref()

        call eos%deteos(bc%edge(m)%pv(1),bc%edge(m)%pv(5),bc%edge(m)%pv(6),bc%edge(m)%pv(7),bc%edge(m)%dv)

        if(config%getiturb().ge.-2) then
          call prop%detprop(bc%edge(m)%dv(3),bc%edge(m)%dv(4),bc%edge(m)%dv(5),bc%edge(m)%pv(5) &
                           ,bc%edge(m)%pv(6),bc%edge(m)%pv(7),bc%edge(m)%tv(1:2))
          if(config%getiturb().ge.-1) then
            bc%edge(m)%tv(3) = config%getemutref()
            bc%edge(m)%pv(8) = config%getkref()
            bc%edge(m)%pv(9) = config%getoref()
          end if
        end if
      end do

      isurf_edge = .false.; jsurf_edge = .false.; ksurf_edge = .false.

      do m=1,4
        do n=1,bc%nbc
          if((trim(bc%edge(n)%bcname).eq.'BCWall').and. &
             (trim(bc%bcinfo(n)%bcname).eq.'BCSymmetryPlane')) then
            if((bc%bcinfo(n)%istart(1).le.bc%edge(m)%neighbor1(1)).and. &
               (bc%bcinfo(n)%iend(1).ge.bc%edge(m)%neighbor1(1)).and.   &
               (bc%bcinfo(n)%istart(2).le.bc%edge(m)%neighbor1(2)).and. &
               (bc%bcinfo(n)%iend(2).ge.bc%edge(m)%neighbor1(2)) ) then
              isurf_edge(m) = .true.
            end if
            if((bc%bcinfo(n)%istart(1).le.bc%edge(m)%neighbor2(1)).and. &
               (bc%bcinfo(n)%iend(1).ge.bc%edge(m)%neighbor2(1)).and.   &
               (bc%bcinfo(n)%istart(2).le.bc%edge(m)%neighbor2(2)).and. &
               (bc%bcinfo(n)%iend(2).ge.bc%edge(m)%neighbor2(2)) ) then
              jsurf_edge(m) = .true.
            end if
          end if
        end do
      end do
      do m=5,8
        do n=1,bc%nbc 
          if((trim(bc%bcinfo(n)%bcname).eq.'BCWall').and. &
             (trim(bc%bcinfo(n)%bcname).eq.'BCSymmetryPlane')) then
            if((bc%bcinfo(n)%istart(1).le.bc%edge(m)%neighbor1(1)).and. &
               (bc%bcinfo(n)%iend(1).ge.bc%edge(m)%neighbor1(1)).and.   &
               (bc%bcinfo(n)%istart(3).le.bc%edge(m)%neighbor1(3)).and. &
               (bc%bcinfo(n)%iend(3).ge.bc%edge(m)%neighbor1(3)) ) then
              isurf_edge(m) = .true.
            end if
            if((bc%bcinfo(n)%istart(1).le.bc%edge(m)%neighbor3(1)).and. &
               (bc%bcinfo(n)%iend(1).ge.bc%edge(m)%neighbor3(1)).and.   &
               (bc%bcinfo(n)%istart(3).le.bc%edge(m)%neighbor3(3)).and. &
               (bc%bcinfo(n)%iend(3).ge.bc%edge(m)%neighbor3(3)) ) then
              ksurf_edge(m) = .true.
            end if
          end if
        end do
      end do
      do m=9,12
        do n=1,bc%nbc 
          if((trim(bc%bcinfo(n)%bcname).eq.'BCWall').and. &
             (trim(bc%bcinfo(n)%bcname).eq.'BCSymmetryPlane')) then
            if((bc%bcinfo(n)%istart(2).le.bc%edge(m)%neighbor2(2)).and. &
               (bc%bcinfo(n)%iend(2).ge.bc%edge(m)%neighbor2(2)).and.   &
               (bc%bcinfo(n)%istart(3).le.bc%edge(m)%neighbor2(3)).and. &
               (bc%bcinfo(n)%iend(3).ge.bc%edge(m)%neighbor2(3)) ) then
              jsurf_edge(m) = .true.
            end if
            if((bc%bcinfo(n)%istart(2).le.bc%edge(m)%neighbor3(2)).and. &
               (bc%bcinfo(n)%iend(2).ge.bc%edge(m)%neighbor3(2)).and.   &
               (bc%bcinfo(n)%istart(3).le.bc%edge(m)%neighbor3(3)).and. &
               (bc%bcinfo(n)%iend(3).ge.bc%edge(m)%neighbor3(3)) ) then
              ksurf_edge(m) = .true.
            end if
          end if
        end do
      end do

      do m=1,4
        if(isurf_edge(m).and.jsurf_edge(m)) then
          if(bc%iturb.eq.0) then
            bc%edge(m)%bctype => edgewallwallkw
          else
            bc%edge(m)%bctype => edgewallwall
          end if
        else if(isurf_edge(m).and.(.not.jsurf_edge(m))) then
          select case(m)
          case(1,3)
            bc%edge(m)%face = 'imin'
          case(2,4)
            bc%edge(m)%face = 'imax'
          end select
          bc%edge(m)%origin(1) = bc%edge(m)%neighbor2(1)
          bc%edge(m)%origin(2) = bc%edge(m)%neighbor2(2)
          select case(bc%iturb)
          case(-1)
            bc%edge(m)%bctype => bcwallviscouske
          case(0)
            bc%edge(m)%bctype => bcwallviscouskw
          case(-2)
            bc%edge(m)%bctype => bcwallviscous
          case default
            bc%edge(m)%bctype => bcwallinviscid
          end select

        else if((.not.isurf_edge(m)).and.jsurf_edge(m)) then
          select case(m)
          case(1,2)
            bc%edge(m)%face = 'jmin'
          case(3,4)
            bc%edge(m)%face = 'jmax'
          end select
          bc%edge(m)%origin(1) = bc%edge(m)%neighbor1(1)
          bc%edge(m)%origin(2) = bc%edge(m)%neighbor1(2)
          select case(bc%iturb)
          case(-1)
            bc%edge(m)%bctype => bcwallviscouske
          case(0)
            bc%edge(m)%bctype => bcwallviscouskw
          case(-2)
            bc%edge(m)%bctype => bcwallviscous
          case default
            bc%edge(m)%bctype => bcwallinviscid
          end select
        else
          bc%edge(m)%bctype => edgenowall        
        end if
        bc%edge(m)%dir(3) = 1
      end do
      
      do m=5,8
        if(isurf_edge(m).and.ksurf_edge(m)) then
          if(bc%iturb.eq.0) then
            bc%edge(m)%bctype => edgewallwallkw
          else
            bc%edge(m)%bctype => edgewallwall
          end if
        else if(isurf_edge(m).and.(.not.ksurf_edge(m))) then
          select case(m)
          case(5,7)
            bc%edge(m)%face = 'imin'
          case(6,8)
            bc%edge(m)%face = 'imax'
          end select
          bc%edge(m)%origin(1) = bc%edge(m)%neighbor3(1)
          bc%edge(m)%origin(3) = bc%edge(m)%neighbor3(3)
          select case(bc%iturb)
          case(-1)
            bc%edge(m)%bctype => bcwallviscouske
          case(0)
            bc%edge(m)%bctype => bcwallviscouskw
          case(-2)
            bc%edge(m)%bctype => bcwallviscous
          case default
            bc%edge(m)%bctype => bcwallinviscid
          end select

        else if((.not.isurf_edge(m)).and.ksurf_edge(m)) then
          select case(m)
          case(5,6)
            bc%edge(m)%face = 'kmin'
          case(7,8)
            bc%edge(m)%face = 'kmax'
          end select
          bc%edge(m)%origin(1) = bc%edge(m)%neighbor1(1)
          bc%edge(m)%origin(3) = bc%edge(m)%neighbor1(3)
          select case(bc%iturb)
          case(-1)
            bc%edge(m)%bctype => bcwallviscouske
          case(0)
            bc%edge(m)%bctype => bcwallviscouskw
          case(-2)
            bc%edge(m)%bctype => bcwallviscous
          case default
            bc%edge(m)%bctype => bcwallinviscid
          end select
        else
          bc%edge(m)%bctype => edgenowall        
        end if
        bc%edge(m)%dir(2) = 1
      end do

      do m=9,12
        if(jsurf_edge(m).and.ksurf_edge(m)) then
          if(bc%iturb.eq.0) then
            bc%edge(m)%bctype => edgewallwallkw
          else
            bc%edge(m)%bctype => edgewallwall
          end if
        else if(jsurf_edge(m).and.(.not.ksurf_edge(m))) then
          select case(m)
          case(9,11)
            bc%edge(m)%face = 'jmin'
          case(10,12)
            bc%edge(m)%face = 'jmax'
          end select
          bc%edge(m)%origin(2) = bc%edge(m)%neighbor3(2)
          bc%edge(m)%origin(3) = bc%edge(m)%neighbor3(3)
          select case(bc%iturb)
          case(-1)
            bc%edge(m)%bctype => bcwallviscouske
          case(0)
            bc%edge(m)%bctype => bcwallviscouskw
          case(-2)
            bc%edge(m)%bctype => bcwallviscous
          case default
            bc%edge(m)%bctype => bcwallinviscid
          end select

        else if((.not.jsurf_edge(m)).and.ksurf_edge(m)) then
          select case(m)
          case(9,10)
            bc%edge(m)%face = 'kmin'
          case(11,12)
            bc%edge(m)%face = 'kmax'
          end select
          bc%edge(m)%origin(2) = bc%edge(m)%neighbor2(2)
          bc%edge(m)%origin(3) = bc%edge(m)%neighbor2(3)
          select case(bc%iturb)
          case(-1)
            bc%edge(m)%bctype => bcwallviscouske
          case(0)
            bc%edge(m)%bctype => bcwallviscouskw
          case(-2)
            bc%edge(m)%bctype => bcwallviscous
          case default
            bc%edge(m)%bctype => bcwallinviscid
          end select
        else
          bc%edge(m)%bctype => edgenowall        
        end if
        bc%edge(m)%dir(1) = 1
      end do
      
      bc%corner(1)%istart(1) = -1                   ; bc%corner(1)%iend(1) = 1
      bc%corner(1)%istart(2) = -1                   ; bc%corner(1)%iend(2) = 1
      bc%corner(1)%istart(3) = -1                   ; bc%corner(1)%iend(3) = 1
      bc%corner(1)%neighbor1(1) = 1                 ; bc%corner(1)%neighbor2(1) = 2 
      bc%corner(1)%neighbor1(2) = 2                 ; bc%corner(1)%neighbor2(2) = 1
      bc%corner(1)%neighbor1(3) = 2                 ; bc%corner(1)%neighbor2(3) = 2
      bc%corner(1)%neighbor3(1) = 2                 ; bc%corner(1)%origin(1) = 2
      bc%corner(1)%neighbor3(2) = 2                 ; bc%corner(1)%origin(2) = 2
      bc%corner(1)%neighbor3(3) = 1                 ; bc%corner(1)%origin(3) = 2
      bc%corner(1)%neighbor4(1) = 2                 ; bc%corner(1)%neighbor5(1) = 1 
      bc%corner(1)%neighbor4(2) = 1                 ; bc%corner(1)%neighbor5(2) = 2
      bc%corner(1)%neighbor4(3) = 1                 ; bc%corner(1)%neighbor5(3) = 1
      bc%corner(1)%neighbor6(1) = 1                 
      bc%corner(1)%neighbor6(2) = 1                 
      bc%corner(1)%neighbor6(3) = 2                

      bc%corner(2)%istart(1) = grid%getimax()+1     ; bc%corner(2)%iend(1) = grid%getimax()+3
      bc%corner(2)%istart(2) = -1                   ; bc%corner(2)%iend(2) = 1
      bc%corner(2)%istart(3) = -1                   ; bc%corner(2)%iend(3) = 1
      bc%corner(2)%neighbor1(1) = grid%getimax()+1  ; bc%corner(2)%neighbor2(1) = grid%getimax() 
      bc%corner(2)%neighbor1(2) = 2                 ; bc%corner(2)%neighbor2(2) = 1
      bc%corner(2)%neighbor1(3) = 2                 ; bc%corner(2)%neighbor2(3) = 2
      bc%corner(2)%neighbor3(1) = grid%getimax()    ; bc%corner(2)%origin(1) = grid%getimax()
      bc%corner(2)%neighbor3(2) = 2                 ; bc%corner(2)%origin(2) = 2
      bc%corner(2)%neighbor3(3) = 1                 ; bc%corner(2)%origin(3) = 2
      bc%corner(2)%neighbor4(1) = grid%getimax()    ; bc%corner(2)%neighbor5(1) = grid%getimax()+1 
      bc%corner(2)%neighbor4(2) = 1                 ; bc%corner(2)%neighbor5(2) = 2
      bc%corner(2)%neighbor4(3) = 1                 ; bc%corner(2)%neighbor5(3) = 1
      bc%corner(2)%neighbor6(1) = grid%getimax()+1    
      bc%corner(2)%neighbor6(2) = 1              
      bc%corner(2)%neighbor6(3) = 2            
      
      bc%corner(3)%istart(1) = -1                   ; bc%corner(3)%iend(1) = 1
      bc%corner(3)%istart(2) = grid%getjmax()+1     ; bc%corner(3)%iend(2) = grid%getjmax()+3
      bc%corner(3)%istart(3) = -1                   ; bc%corner(3)%iend(3) = 1
      bc%corner(3)%neighbor1(1) = 1                 ; bc%corner(3)%neighbor2(1) = 2 
      bc%corner(3)%neighbor1(2) = grid%getjmax()    ; bc%corner(3)%neighbor2(2) = grid%getjmax()+1
      bc%corner(3)%neighbor1(3) = 2                 ; bc%corner(3)%neighbor2(3) = 2
      bc%corner(3)%neighbor3(1) = 2                 ; bc%corner(3)%origin(1) = 2
      bc%corner(3)%neighbor3(2) = grid%getjmax()    ; bc%corner(3)%origin(2) = grid%getjmax()
      bc%corner(3)%neighbor3(3) = 1                 ; bc%corner(3)%origin(3) = 2
      bc%corner(3)%neighbor4(1) = 2                 ; bc%corner(3)%neighbor5(1) = 1 
      bc%corner(3)%neighbor4(2) = grid%getjmax()+1  ; bc%corner(3)%neighbor5(2) = grid%getjmax()
      bc%corner(3)%neighbor4(3) = 1                 ; bc%corner(3)%neighbor5(3) = 1
      bc%corner(3)%neighbor6(1) = 1                
      bc%corner(3)%neighbor6(2) = grid%getjmax()+1  
      bc%corner(3)%neighbor6(3) = 2               
      
      bc%corner(4)%istart(1) = grid%getimax()+1     ; bc%corner(4)%iend(1) = grid%getimax()+3
      bc%corner(4)%istart(2) = grid%getjmax()+1     ; bc%corner(4)%iend(2) = grid%getjmax()+3
      bc%corner(4)%istart(3) = -1                   ; bc%corner(4)%iend(3) = 1
      bc%corner(4)%neighbor1(1) = grid%getimax()+1  ; bc%corner(4)%neighbor2(1) = grid%getimax() 
      bc%corner(4)%neighbor1(2) = grid%getjmax()    ; bc%corner(4)%neighbor2(2) = grid%getjmax()+1
      bc%corner(4)%neighbor1(3) = 2                 ; bc%corner(4)%neighbor2(3) = 2
      bc%corner(4)%neighbor3(1) = grid%getimax()    ; bc%corner(4)%origin(1) = grid%getimax()
      bc%corner(4)%neighbor3(2) = grid%getjmax()    ; bc%corner(4)%origin(2) = grid%getjmax()
      bc%corner(4)%neighbor3(3) = 1                 ; bc%corner(4)%origin(3) = 2
      bc%corner(4)%neighbor4(1) = grid%getimax()    ; bc%corner(4)%neighbor5(1) = grid%getimax()+1 
      bc%corner(4)%neighbor4(2) = grid%getjmax()+1  ; bc%corner(4)%neighbor5(2) = grid%getjmax()
      bc%corner(4)%neighbor4(3) = 1                 ; bc%corner(4)%neighbor5(3) = 1
      bc%corner(4)%neighbor6(1) = grid%getimax()+1   
      bc%corner(4)%neighbor6(2) = grid%getjmax()+1   
      bc%corner(4)%neighbor6(3) = 2             
      
      bc%corner(5)%istart(1) = -1                   ; bc%corner(5)%iend(1) = 1
      bc%corner(5)%istart(2) = -1                   ; bc%corner(5)%iend(2) = 1
      bc%corner(5)%istart(3) = grid%getkmax()+1     ; bc%corner(5)%iend(3) = grid%getkmax()+3
      bc%corner(5)%neighbor1(1) = 1                 ; bc%corner(5)%neighbor2(1) = 2 
      bc%corner(5)%neighbor1(2) = 2                 ; bc%corner(5)%neighbor2(2) = 1
      bc%corner(5)%neighbor1(3) = grid%getkmax()    ; bc%corner(5)%neighbor2(3) = grid%getkmax()
      bc%corner(5)%neighbor3(1) = 2                 ; bc%corner(5)%origin(1) = 2
      bc%corner(5)%neighbor3(2) = 2                 ; bc%corner(5)%origin(2) = 2
      bc%corner(5)%neighbor3(3) = grid%getkmax()+1  ; bc%corner(5)%origin(3) = grid%getkmax()
      bc%corner(5)%neighbor4(1) = 2                 ; bc%corner(5)%neighbor5(1) = 1 
      bc%corner(5)%neighbor4(2) = 1                 ; bc%corner(5)%neighbor5(2) = 2
      bc%corner(5)%neighbor4(3) = grid%getkmax()+1  ; bc%corner(5)%neighbor5(3) = grid%getkmax()+1
      bc%corner(5)%neighbor6(1) = 1                
      bc%corner(5)%neighbor6(2) = 1                
      bc%corner(5)%neighbor6(3) = grid%getkmax() 
      
      bc%corner(6)%istart(1) = grid%getimax()+1     ; bc%corner(6)%iend(1) = grid%getimax()+3
      bc%corner(6)%istart(2) = -1                   ; bc%corner(6)%iend(2) = 1
      bc%corner(6)%istart(3) = grid%getkmax()+1     ; bc%corner(6)%iend(3) = grid%getkmax()+3
      bc%corner(6)%neighbor1(1) = grid%getimax()+1  ; bc%corner(6)%neighbor2(1) = grid%getimax() 
      bc%corner(6)%neighbor1(2) = 2                 ; bc%corner(6)%neighbor2(2) = 1
      bc%corner(6)%neighbor1(3) = grid%getkmax()    ; bc%corner(6)%neighbor2(3) = grid%getkmax()
      bc%corner(6)%neighbor3(1) = grid%getimax()    ; bc%corner(6)%origin(1) = grid%getimax()
      bc%corner(6)%neighbor3(2) = 1                 ; bc%corner(6)%origin(2) = 2
      bc%corner(6)%neighbor3(3) = grid%getkmax()+1  ; bc%corner(6)%origin(3) = grid%getkmax()
      bc%corner(6)%neighbor4(1) = grid%getimax()    ; bc%corner(6)%neighbor5(1) = grid%getimax()+1 
      bc%corner(6)%neighbor4(2) = 1                 ; bc%corner(6)%neighbor5(2) = 2
      bc%corner(6)%neighbor4(3) = grid%getkmax()+1  ; bc%corner(6)%neighbor5(3) = grid%getkmax()+1
      bc%corner(6)%neighbor6(1) = grid%getimax()+1   
      bc%corner(6)%neighbor6(2) = 1                
      bc%corner(6)%neighbor6(3) = grid%getkmax() 
      
      bc%corner(7)%istart(1) = -1                   ; bc%corner(7)%iend(1) = 1
      bc%corner(7)%istart(2) = grid%getjmax()+1     ; bc%corner(7)%iend(2) = grid%getjmax()+3
      bc%corner(7)%istart(3) = grid%getkmax()+1     ; bc%corner(7)%iend(3) = grid%getkmax()+3
      bc%corner(7)%neighbor1(1) = 1                 ; bc%corner(7)%neighbor2(1) = 2 
      bc%corner(7)%neighbor1(2) = grid%getjmax()    ; bc%corner(7)%neighbor2(2) = grid%getjmax()+1
      bc%corner(7)%neighbor1(3) = grid%getkmax()    ; bc%corner(7)%neighbor2(3) = grid%getkmax()
      bc%corner(7)%neighbor3(1) = 2                 ; bc%corner(7)%origin(1) = 2
      bc%corner(7)%neighbor3(2) = grid%getjmax()    ; bc%corner(7)%origin(2) = grid%getjmax()
      bc%corner(7)%neighbor3(3) = grid%getkmax()+1  ; bc%corner(7)%origin(3) = grid%getkmax()
      bc%corner(7)%neighbor4(1) = 2                 ; bc%corner(7)%neighbor5(1) = 1 
      bc%corner(7)%neighbor4(2) = grid%getjmax()+1  ; bc%corner(7)%neighbor5(2) = grid%getjmax()
      bc%corner(7)%neighbor4(3) = grid%getkmax()+1  ; bc%corner(7)%neighbor5(3) = grid%getkmax()+1
      bc%corner(7)%neighbor6(1) = 1                
      bc%corner(7)%neighbor6(2) = grid%getjmax()+1   
      bc%corner(7)%neighbor6(3) = grid%getkmax() 
      
      bc%corner(8)%istart(1) = grid%getimax()+1     ; bc%corner(8)%iend(1) = grid%getimax()+3
      bc%corner(8)%istart(2) = grid%getjmax()+1     ; bc%corner(8)%iend(2) = grid%getjmax()+3
      bc%corner(8)%istart(3) = grid%getkmax()+1     ; bc%corner(8)%iend(3) = grid%getkmax()+3
      bc%corner(8)%neighbor1(1) = grid%getimax()+1  ; bc%corner(8)%neighbor2(1) = grid%getimax() 
      bc%corner(8)%neighbor1(2) = grid%getjmax()    ; bc%corner(8)%neighbor2(2) = grid%getjmax()+1
      bc%corner(8)%neighbor1(3) = grid%getkmax()    ; bc%corner(8)%neighbor2(3) = grid%getkmax()
      bc%corner(8)%neighbor3(1) = grid%getimax()    ; bc%corner(8)%origin(1) = grid%getimax()
      bc%corner(8)%neighbor3(2) = grid%getjmax()    ; bc%corner(8)%origin(2) = grid%getjmax()
      bc%corner(8)%neighbor3(3) = grid%getkmax()+1  ; bc%corner(8)%origin(3) = grid%getkmax()
      bc%corner(8)%neighbor4(1) = grid%getimax()    ; bc%corner(8)%neighbor5(1) = grid%getimax()+1 
      bc%corner(8)%neighbor4(2) = grid%getjmax()+1  ; bc%corner(8)%neighbor5(2) = grid%getjmax()
      bc%corner(8)%neighbor4(3) = grid%getkmax()+1  ; bc%corner(8)%neighbor5(3) = grid%getkmax()+1
      bc%corner(8)%neighbor6(1) = grid%getimax()+1   
      bc%corner(8)%neighbor6(2) = grid%getjmax()+1   
      bc%corner(8)%neighbor6(3) = grid%getkmax()
      
      do n=1,8
        bc%corner(n)%dir = 0
        allocate(bc%corner(n)%pv(bc%npv),bc%corner(n)%dv(bc%ndv),bc%corner(n)%tv(bc%ntv))

        bc%corner(n)%pv(1) = config%getpref()
        bc%corner(n)%pv(2) = config%geturef()*dcos(config%getaos())*dcos(config%getaoa())
        bc%corner(n)%pv(3) = config%geturef()*dcos(config%getaos())*dsin(config%getaoa())
        bc%corner(n)%pv(4) = config%geturef()*dsin(config%getaos())
        bc%corner(n)%pv(5) = config%gettref()
        bc%corner(n)%pv(6) = config%gety1ref()
        bc%corner(n)%pv(7) = config%gety2ref()

        call eos%deteos(bc%corner(n)%pv(1),bc%corner(n)%pv(5),bc%corner(n)%pv(6),bc%corner(n)%pv(7),bc%corner(n)%dv)

        if(config%getiturb().ge.-2) then
          call prop%detprop(bc%corner(n)%dv(3),bc%corner(n)%dv(4),bc%corner(n)%dv(5),bc%corner(n)%pv(5) &
                           ,bc%corner(n)%pv(6),bc%corner(n)%pv(7),bc%corner(n)%tv(1:2))
          if(config%getiturb().ge.-1) then
            bc%corner(n)%tv(3) = config%getemutref()
            bc%corner(n)%pv(8) = config%getkref()
            bc%corner(n)%pv(9) = config%getoref()
          end if
        end if
      end do

      isurf_corner = .false.; jsurf_corner = .false.; ksurf_corner = .false.

      do m=1,8
        do n=1,bc%nbc
          if((trim(bc%bcinfo(n)%bcname).eq.'BCWall').and. &
             (trim(bc%bcinfo(n)%bcname).eq.'BCSymmetryPlane')) then
            if((bc%bcinfo(n)%istart(1).le.bc%corner(m)%neighbor1(1)).and. &
               (bc%bcinfo(n)%iend(1).ge.bc%corner(m)%neighbor1(1)).and.   &
               (bc%bcinfo(n)%istart(2).le.bc%corner(m)%neighbor1(2)).and. &
               (bc%bcinfo(n)%iend(2).ge.bc%corner(m)%neighbor1(2)).and.   &
               (bc%bcinfo(n)%istart(3).le.bc%corner(m)%neighbor1(3)).and. &
               (bc%bcinfo(n)%iend(3).ge.bc%corner(m)%neighbor1(3)) ) then
              isurf_corner(m) = .true.
            end if
            if((bc%bcinfo(n)%istart(1).le.bc%corner(m)%neighbor2(1)).and. &
               (bc%bcinfo(n)%iend(1).ge.bc%corner(m)%neighbor2(1)).and.   &
               (bc%bcinfo(n)%istart(2).le.bc%corner(m)%neighbor2(2)).and. &
               (bc%bcinfo(n)%iend(2).ge.bc%corner(m)%neighbor2(2)).and.   &
               (bc%bcinfo(n)%istart(3).le.bc%corner(m)%neighbor2(3)).and. &
               (bc%bcinfo(n)%iend(3).ge.bc%corner(m)%neighbor2(3)) ) then
              jsurf_corner(m) = .true.
            end if
            if((bc%bcinfo(n)%istart(1).le.bc%corner(m)%neighbor3(1)).and. &
               (bc%bcinfo(n)%iend(1).ge.bc%corner(m)%neighbor3(1)).and.   &
               (bc%bcinfo(n)%istart(2).le.bc%corner(m)%neighbor3(2)).and. &
               (bc%bcinfo(n)%iend(2).ge.bc%corner(m)%neighbor3(2)).and.   &
               (bc%bcinfo(n)%istart(3).le.bc%corner(m)%neighbor3(3)).and. &
               (bc%bcinfo(n)%iend(3).ge.bc%corner(m)%neighbor3(3)) ) then
              ksurf_corner(m) = .true.
            end if
          end if
        end do
      end do
      

      do m =1,8
        if(isurf_corner(m).and.jsurf_corner(m).and.ksurf_corner(m)) then
          if(bc%iturb.eq.0) then
            bc%corner(m)%bctype => cornerwallwallkw
          else
            bc%corner(m)%bctype => cornerwallwall
          end if
        else if((.not.isurf_corner(m)).and.jsurf_corner(m).and.ksurf_corner(m)) then
          bc%corner(m)%origin(1) = bc%corner(m)%neighbor1(1)
          bc%corner(m)%origin(2) = bc%corner(m)%neighbor1(2)
          bc%corner(m)%origin(3) = bc%corner(m)%neighbor1(3)
          bc%corner(m)%neighbor1(1) = 0
          bc%corner(m)%neighbor1(2) = 0
          bc%corner(m)%neighbor1(3) = 0
          bc%corner(m)%neighbor2(1) = bc%corner(m)%neighbor5(1)
          bc%corner(m)%neighbor2(2) = bc%corner(m)%neighbor5(2)
          bc%corner(m)%neighbor2(3) = bc%corner(m)%neighbor5(3)
          bc%corner(m)%neighbor3(1) = bc%corner(m)%neighbor6(1)
          bc%corner(m)%neighbor3(2) = bc%corner(m)%neighbor6(2)
          bc%corner(m)%neighbor3(3) = bc%corner(m)%neighbor6(3)
          if(bc%iturb.eq.0) then
            bc%corner(m)%bctype => edgewallwallkw
          else
            bc%corner(m)%bctype => edgewallwall
          end if
        else if(isurf_corner(m).and.(.not.jsurf_corner(m)).and.ksurf_corner(m)) then
          bc%corner(m)%origin(1) = bc%corner(m)%neighbor2(1)
          bc%corner(m)%origin(2) = bc%corner(m)%neighbor2(2)
          bc%corner(m)%origin(3) = bc%corner(m)%neighbor2(3)
          bc%corner(m)%neighbor1(1) = bc%corner(m)%neighbor4(1)
          bc%corner(m)%neighbor1(2) = bc%corner(m)%neighbor4(2)
          bc%corner(m)%neighbor1(3) = bc%corner(m)%neighbor4(3)
          bc%corner(m)%neighbor2(1) = 0
          bc%corner(m)%neighbor2(2) = 0
          bc%corner(m)%neighbor2(3) = 0
          bc%corner(m)%neighbor3(1) = bc%corner(m)%neighbor6(1)
          bc%corner(m)%neighbor3(2) = bc%corner(m)%neighbor6(2)
          bc%corner(m)%neighbor3(3) = bc%corner(m)%neighbor6(3)
          if(bc%iturb.eq.0) then
            bc%corner(m)%bctype => edgewallwallkw
          else
            bc%corner(m)%bctype => edgewallwall
          end if
        else if(isurf_corner(m).and.jsurf_corner(m).and.(.not.ksurf_corner(m))) then
          bc%corner(m)%origin(1) = bc%corner(m)%neighbor3(1)
          bc%corner(m)%origin(2) = bc%corner(m)%neighbor3(2)
          bc%corner(m)%origin(3) = bc%corner(m)%neighbor3(3)
          bc%corner(m)%neighbor1(1) = bc%corner(m)%neighbor4(1)
          bc%corner(m)%neighbor1(2) = bc%corner(m)%neighbor4(2)
          bc%corner(m)%neighbor1(3) = bc%corner(m)%neighbor4(3)
          bc%corner(m)%neighbor2(1) = bc%corner(m)%neighbor5(1)
          bc%corner(m)%neighbor2(2) = bc%corner(m)%neighbor5(2)
          bc%corner(m)%neighbor2(3) = bc%corner(m)%neighbor5(3)
          bc%corner(m)%neighbor3(1) = 0
          bc%corner(m)%neighbor3(2) = 0
          bc%corner(m)%neighbor3(3) = 0
          if(bc%iturb.eq.0) then
            bc%corner(m)%bctype => edgewallwallkw
          else
            bc%corner(m)%bctype => edgewallwall
          end if
        else if((.not.isurf_corner(m)).and.(.not.jsurf_corner(m)).and.ksurf_corner(m)) then
          select case(m)
          case(1,2,3,4)
            bc%corner(m)%face = 'Kmin'
          case(5,6,7,8)
            bc%corner(m)%face = 'Kmax'
          end select
          bc%corner(m)%origin(1) = bc%corner(m)%neighbor6(1)
          bc%corner(m)%origin(2) = bc%corner(m)%neighbor6(2)
          bc%corner(m)%origin(3) = bc%corner(m)%neighbor6(3)
          select case(bc%iturb)
          case(-1)
            bc%corner(m)%bctype => bcwallviscouske
          case(0)
            bc%corner(m)%bctype => bcwallviscouskw
          case(-2)
            bc%corner(m)%bctype => bcwallviscous
          case default
            bc%corner(m)%bctype => bcwallinviscid
          end select
        else if(isurf_corner(m).and.(.not.jsurf_corner(m)).and.(.not.ksurf_corner(m))) then
          select case(m)
          case(1,3,5,7)
            bc%corner(m)%face = 'imin'
          case(2,4,6,8)
            bc%corner(m)%face = 'imax'
          end select
          bc%corner(m)%origin(1) = bc%corner(m)%neighbor4(1)
          bc%corner(m)%origin(2) = bc%corner(m)%neighbor4(2)
          bc%corner(m)%origin(3) = bc%corner(m)%neighbor4(3)
          select case(bc%iturb)
          case(-1)
            bc%corner(m)%bctype => bcwallviscouske
          case(0)
            bc%corner(m)%bctype => bcwallviscouskw
          case(-2)
            bc%corner(m)%bctype => bcwallviscous
          case default
            bc%corner(m)%bctype => bcwallinviscid
          end select
        else if((.not.isurf_corner(m)).and.jsurf_corner(m).and.(.not.ksurf_corner(m))) then
          select case(m)
          case(1,2,5,6)
            bc%corner(m)%face = 'jmin'
          case(3,4,7,8)
            bc%corner(m)%face = 'jmax'
          end select
          bc%corner(m)%origin(1) = bc%corner(m)%neighbor5(1)
          bc%corner(m)%origin(2) = bc%corner(m)%neighbor5(2)
          bc%corner(m)%origin(3) = bc%corner(m)%neighbor5(3)
          select case(bc%iturb)
          case(-1)
            bc%corner(m)%bctype => bcwallviscouske
          case(0)
            bc%corner(m)%bctype => bcwallviscouskw
          case(-2)
            bc%corner(m)%bctype => bcwallviscous
          case default
            bc%corner(m)%bctype => bcwallinviscid
          end select
        else
          bc%corner(m)%bctype => cornernowall
        end if 
      end do

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(bc)
      implicit none
      class(t_bc), intent(inout) :: bc
      integer :: n

      if(associated(bc%prec%getsndp2)) nullify(bc%prec%getsndp2)

      do n=1,bc%nbc
        if(associated(bc%bcinfo(n)%bctype)) nullify(bc%bcinfo(n)%bctype)
        deallocate(bc%bcinfo(n)%pv,bc%bcinfo(n)%dv,bc%bcinfo(n)%tv)
      end do
      
      do n=1,12
        if(associated(bc%edge(n)%bctype)) nullify(bc%edge(n)%bctype)
        deallocate(bc%edge(n)%pv,bc%edge(n)%dv,bc%edge(n)%tv)
      end do
      
      do n=1,8
        if(associated(bc%corner(n)%bctype)) nullify(bc%corner(n)%bctype)
        deallocate(bc%corner(n)%pv,bc%corner(n)%dv,bc%corner(n)%tv)
      end do

      do n=1,bc%ncon
        if(allocated(bc%mpitemp(n)%sendbuf)) deallocate(bc%mpitemp(n)%sendbuf)
        if(allocated(bc%mpitemp(n)%recvbuf)) deallocate(bc%mpitemp(n)%recvbuf)
      end do

      deallocate(bc%bcinfo,bc%edge,bc%corner,bc%connectinfo)
      deallocate(bc%mpitemp)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setbc(bc,grid,variable,eos,prop)
      implicit none
      class(t_bc), intent(inout) :: bc
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,n,m,l
      integer :: ii,jj,kk
      integer :: ier
      integer :: request_s(bc%ncon),request_r(bc%ncon),request_sa(bc%ncon),request_ra(bc%ncon)
      integer :: status(mpi_status_size)
      real(8) :: pv(bc%npv),dv(bc%ndv),tv(bc%ntv)
      
      do n=1,bc%nbc
        call bc%bcinfo(n)%bctype(grid,variable,eos,prop,bc%prec)
      end do
 
      do n=1,bc%ncon
        if(bc%rank.ne.bc%connectinfo(n)%donor) then 
          call mpi_irecv(bc%mpitemp(n)%recvadress,21,mpi_integer,bc%connectinfo(n)%donor,bc%connectinfo(n)%donor+bc%size,mpi_comm_world,request_ra(n),ier)
          call mpi_irecv(bc%mpitemp(n)%recvbuf,bc%mpitemp(n)%num,mpi_real8,bc%connectinfo(n)%donor,bc%connectinfo(n)%donor,mpi_comm_world,request_r(n),ier)
        end if
      end do
     
      do n=1,bc%ncon
        l = 0
        do k=bc%connectinfo(n)%istart(3),bc%connectinfo(n)%iend(3)
          do j=bc%connectinfo(n)%istart(2),bc%connectinfo(n)%iend(2)
            do i=bc%connectinfo(n)%istart(1),bc%connectinfo(n)%iend(1)
              pv = variable%getpv(i,j,k)
              tv = variable%gettv(i,j,k)
              dv = variable%getdv(i,j,k)
              do m=1,bc%npv
                l = l + 1
                bc%mpitemp(n)%sendbuf(l) = pv(m)
              end do
              do m=1,bc%ntv
                l = l + 1
                bc%mpitemp(n)%sendbuf(l) = tv(m)
              end do
              do m=1,bc%ndv
                l = l + 1
                bc%mpitemp(n)%sendbuf(l) = dv(m)
              end do            
            end do
          end do
        end do
          
        if(bc%rank.ne.bc%connectinfo(n)%donor) then 
          bc%mpitemp(n)%sendadress = (/bc%connectinfo(n)%istart(1),bc%connectinfo(n)%iend(1), &
                                       bc%connectinfo(n)%istart(2),bc%connectinfo(n)%iend(2), &
                                       bc%connectinfo(n)%istart(3),bc%connectinfo(n)%iend(3), &
                                       bc%connectinfo(n)%istart_donor(1), &
                                       bc%connectinfo(n)%iend_donor(1),   &
                                       bc%connectinfo(n)%istart_donor(2), &
                                       bc%connectinfo(n)%iend_donor(2),   &
                                       bc%connectinfo(n)%istart_donor(3), &
                                       bc%connectinfo(n)%iend_donor(3),   &  
                                       bc%connectinfo(n)%transmat(1,1),   &
                                       bc%connectinfo(n)%transmat(1,2),   &
                                       bc%connectinfo(n)%transmat(1,3),   &
                                       bc%connectinfo(n)%transmat(2,1),   &
                                       bc%connectinfo(n)%transmat(2,2),   &
                                       bc%connectinfo(n)%transmat(2,3),   &
                                       bc%connectinfo(n)%transmat(3,1),   &
                                       bc%connectinfo(n)%transmat(3,2),   &
                                       bc%connectinfo(n)%transmat(3,3)/)
                     

          call mpi_isend(bc%mpitemp(n)%sendadress,21,mpi_integer,bc%connectinfo(n)%donor,bc%rank+bc%size,mpi_comm_world,request_sa(n),ier)
          call mpi_isend(bc%mpitemp(n)%sendbuf,bc%mpitemp(n)%num,mpi_real8,bc%connectinfo(n)%donor,bc%rank,mpi_comm_world,request_s(n),ier)
        else  
          l = 0
          do k=bc%connectinfo(n)%istart(3),bc%connectinfo(n)%iend(3)
            do j=bc%connectinfo(n)%istart(2),bc%connectinfo(n)%iend(2)
              do i=bc%connectinfo(n)%istart(1),bc%connectinfo(n)%iend(1)
                ii = bc%connectinfo(n)%transmat(1,1)*(i-bc%connectinfo(n)%istart(1)) &
                   + bc%connectinfo(n)%transmat(1,2)*(j-bc%connectinfo(n)%istart(2)) &
                   + bc%connectinfo(n)%transmat(1,3)*(k-bc%connectinfo(n)%istart(3)) &
                   + bc%connectinfo(n)%istart_donor(1)
                jj = bc%connectinfo(n)%transmat(2,1)*(i-bc%connectinfo(n)%istart(1)) &
                   + bc%connectinfo(n)%transmat(2,2)*(j-bc%connectinfo(n)%istart(2)) &
                   + bc%connectinfo(n)%transmat(2,3)*(k-bc%connectinfo(n)%istart(3)) &
                   + bc%connectinfo(n)%istart_donor(2)
                kk = bc%connectinfo(n)%transmat(3,1)*(i-bc%connectinfo(n)%istart(1)) &
                   + bc%connectinfo(n)%transmat(3,2)*(j-bc%connectinfo(n)%istart(2)) &
                   + bc%connectinfo(n)%transmat(3,3)*(k-bc%connectinfo(n)%istart(3)) &
                   + bc%connectinfo(n)%istart_donor(3)
                do m=1,bc%npv
                  l = l + 1
                  call variable%setpv(m,ii,jj,kk,bc%mpitemp(n)%sendbuf(l))
                end do
                do m=1,bc%ntv
                  l = l + 1
                  call variable%settv(m,ii,jj,kk,bc%mpitemp(n)%sendbuf(l))
                end do
                do m=1,bc%ndv
                  l = l + 1
                  call variable%setdv(m,ii,jj,kk,bc%mpitemp(n)%sendbuf(l))
                end do
              end do
            end do
          end do
        end if
      end do
      
      do n=1,bc%ncon
        if(bc%rank.ne.bc%connectinfo(n)%donor) then
          call mpi_wait(request_s(n),status,ier)
          call mpi_wait(request_r(n),status,ier)
          call mpi_wait(request_sa(n),status,ier)
          call mpi_wait(request_ra(n),status,ier)
        end if
      end do

      do n=1,bc%ncon
        if(bc%rank.ne.bc%connectinfo(n)%donor) then 
          l = 0
          do k=bc%mpitemp(n)%recvadress(5),bc%mpitemp(n)%recvadress(6)
            do j=bc%mpitemp(n)%recvadress(3),bc%mpitemp(n)%recvadress(4)
              do i=bc%mpitemp(n)%recvadress(1),bc%mpitemp(n)%recvadress(2)
                ii = bc%mpitemp(n)%recvadress(13)*(i-bc%mpitemp(n)%recvadress(1)) &
                   + bc%mpitemp(n)%recvadress(14)*(j-bc%mpitemp(n)%recvadress(3)) &
                   + bc%mpitemp(n)%recvadress(15)*(k-bc%mpitemp(n)%recvadress(5)) &
                   + bc%mpitemp(n)%recvadress(7)
                jj = bc%mpitemp(n)%recvadress(16)*(i-bc%mpitemp(n)%recvadress(1)) &
                   + bc%mpitemp(n)%recvadress(17)*(j-bc%mpitemp(n)%recvadress(3)) &
                   + bc%mpitemp(n)%recvadress(18)*(k-bc%mpitemp(n)%recvadress(5)) &
                   + bc%mpitemp(n)%recvadress(9)
                kk = bc%mpitemp(n)%recvadress(19)*(i-bc%mpitemp(n)%recvadress(1)) &
                   + bc%mpitemp(n)%recvadress(20)*(j-bc%mpitemp(n)%recvadress(3)) &
                   + bc%mpitemp(n)%recvadress(21)*(k-bc%mpitemp(n)%recvadress(5)) &
                   + bc%mpitemp(n)%recvadress(11)
                do m=1,bc%npv
                  l = l + 1
                  call variable%setpv(m,ii,jj,kk,bc%mpitemp(n)%recvbuf(l))
                end do
                do m=1,bc%ntv
                  l = l + 1
                  call variable%settv(m,ii,jj,kk,bc%mpitemp(n)%recvbuf(l))
                end do
                do m=1,bc%ndv
                  l = l + 1
                  call variable%setdv(m,ii,jj,kk,bc%mpitemp(n)%recvbuf(l))
                end do
              end do
            end do
          end do
        end if
      end do
      
      do n=1,12
        call bc%edge(n)%bctype(grid,variable,eos,prop,bc%prec)
      end do

      do n=1,8
        call bc%corner(n)%bctype(grid,variable,eos,prop,bc%prec)
      end do
      
    end subroutine setbc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine edgewallwall(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: ppv(variable%getnpv())

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            ppv = -variable%getpv(ii,jj,kk)
            if((bcinfo%neighbor1(1).ne.0).or.(bcinfo%neighbor1(2).ne.0).or.(bcinfo%neighbor1(3).ne.0)) then
              ii = bcinfo%neighbor1(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor1(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor1(3)+bcinfo%dir(3)*k
              ppv = ppv - variable%getpv(ii,jj,kk)
            end if
            if((bcinfo%neighbor2(1).ne.0).or.(bcinfo%neighbor2(2).ne.0).or.(bcinfo%neighbor2(3).ne.0)) then
              ii = bcinfo%neighbor2(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor2(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor2(3)+bcinfo%dir(3)*k
              ppv = ppv - variable%getpv(ii,jj,kk)
            end if
            if((bcinfo%neighbor3(1).ne.0).or.(bcinfo%neighbor3(2).ne.0).or.(bcinfo%neighbor3(3).ne.0)) then
              ii = bcinfo%neighbor3(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor3(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor3(3)+bcinfo%dir(3)*k
              ppv = ppv - variable%getpv(ii,jj,kk)
            end if
            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4,8,9)
                call variable%setpv(m,i,j,k,ppv(m))
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    end subroutine edgewallwall
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine edgewallwallkw(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: ppv(variable%getnpv()),opv(variable%getnpv()),grd(grid%getngrd())
      real(8) :: var

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            grd = grid%getgrd(ii,jj,kk)
            var = 800.d0*tv(1)/dv(1)/grd(5)**2
            ppv = -variable%getpv(ii,jj,kk)
            if((bcinfo%neighbor1(1).ne.0).or.(bcinfo%neighbor1(2).ne.0).or.(bcinfo%neighbor1(3).ne.0)) then
              ii = bcinfo%neighbor1(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor1(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor1(3)+bcinfo%dir(3)*k
              ppv = ppv - variable%getpv(ii,jj,kk)
            end if
            if((bcinfo%neighbor2(1).ne.0).or.(bcinfo%neighbor2(2).ne.0).or.(bcinfo%neighbor2(3).ne.0)) then
              ii = bcinfo%neighbor2(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor2(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor2(3)+bcinfo%dir(3)*k
              ppv = ppv - variable%getpv(ii,jj,kk)
            end if
            if((bcinfo%neighbor3(1).ne.0).or.(bcinfo%neighbor3(2).ne.0).or.(bcinfo%neighbor3(3).ne.0)) then
              ii = bcinfo%neighbor3(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor3(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor3(3)+bcinfo%dir(3)*k
              ppv = ppv - variable%getpv(ii,jj,kk)
            end if
            opv = 4.d0*var + ppv
            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4,8)
                call variable%setpv(m,i,j,k,ppv(m))
              case(9)
                call variable%setpv(m,i,j,k,opv(m))
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    end subroutine edgewallwallkw
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine cornerwallwall(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: ppv(variable%getnpv())

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            ppv = -variable%getpv(ii,jj,kk)
            
            ii = bcinfo%neighbor1(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor1(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor1(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)
            
            ii = bcinfo%neighbor2(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor2(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor2(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor3(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor3(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor3(3)+bcinfo%dir(3)*k          
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor4(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor4(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor4(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor5(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor5(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor5(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor6(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor6(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor6(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4,8,9)
                call variable%setpv(m,i,j,k,ppv(m))
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    end subroutine cornerwallwall
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine cornerwallwallkw(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: ppv(variable%getnpv()),opv(variable%getnpv()),grd(grid%getngrd())
      real(8) :: var

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            grd = grid%getgrd(ii,jj,kk)
            var = 800.d0*tv(1)/dv(1)/grd(5)**2
            ppv = -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor1(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor1(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor1(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)
            
            ii = bcinfo%neighbor2(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor2(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor2(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor3(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor3(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor3(3)+bcinfo%dir(3)*k          
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor4(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor4(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor4(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor5(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor5(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor5(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            ii = bcinfo%neighbor6(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor6(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor6(3)+bcinfo%dir(3)*k
            ppv = ppv -variable%getpv(ii,jj,kk)

            opv = 8.d0*var + ppv
            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4,8)
                call variable%setpv(m,i,j,k,ppv(m))
              case(9)
                call variable%setpv(m,i,j,k,opv(m))
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    end subroutine cornerwallwallkw
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine edgenowall(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
 
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            if((bcinfo%neighbor1(1).ne.0).or.(bcinfo%neighbor1(2).ne.0).or.(bcinfo%neighbor1(3).ne.0)) then
              ii = bcinfo%neighbor1(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor1(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor1(3)+bcinfo%dir(3)*k
              pv = pv + variable%getpv(ii,jj,kk)
              tv = tv + variable%gettv(ii,jj,kk)
            end if
            if((bcinfo%neighbor2(1).ne.0).or.(bcinfo%neighbor2(2).ne.0).or.(bcinfo%neighbor2(3).ne.0)) then
              ii = bcinfo%neighbor2(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor2(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor2(3)+bcinfo%dir(3)*k
              pv = pv + variable%getpv(ii,jj,kk)
              tv = tv + variable%gettv(ii,jj,kk)
            end if
            if((bcinfo%neighbor3(1).ne.0).or.(bcinfo%neighbor3(2).ne.0).or.(bcinfo%neighbor3(3).ne.0)) then
              ii = bcinfo%neighbor3(1)+bcinfo%dir(1)*i
              jj = bcinfo%neighbor3(2)+bcinfo%dir(2)*j
              kk = bcinfo%neighbor3(3)+bcinfo%dir(3)*k
              pv = pv + variable%getpv(ii,jj,kk)
              tv = tv + variable%gettv(ii,jj,kk)
            end if
            pv = 1.d0/3.d0*pv
            tv = 1.d0/3.d0*tv
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+bcinfo%pv(1),pv(5),pv(6),pv(7),dv)
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do

    end subroutine edgenowall
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine cornernowall(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
 
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            ii = bcinfo%neighbor1(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor1(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor1(3)+bcinfo%dir(3)*k
            pv = pv + variable%getpv(ii,jj,kk)
            tv = tv + variable%gettv(ii,jj,kk)
            ii = bcinfo%neighbor2(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor2(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor2(3)+bcinfo%dir(3)*k
            pv = pv + variable%getpv(ii,jj,kk)
            tv = tv + variable%gettv(ii,jj,kk)
            ii = bcinfo%neighbor3(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor3(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor3(3)+bcinfo%dir(3)*k
            pv = pv + variable%getpv(ii,jj,kk)
            tv = tv + variable%gettv(ii,jj,kk)
            ii = bcinfo%neighbor4(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor4(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor4(3)+bcinfo%dir(3)*k
            pv = pv + variable%getpv(ii,jj,kk)
            tv = tv + variable%gettv(ii,jj,kk)
            ii = bcinfo%neighbor5(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor5(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor5(3)+bcinfo%dir(3)*k
            pv = pv + variable%getpv(ii,jj,kk)
            tv = tv + variable%gettv(ii,jj,kk)
            ii = bcinfo%neighbor6(1)+bcinfo%dir(1)*i
            jj = bcinfo%neighbor6(2)+bcinfo%dir(2)*j
            kk = bcinfo%neighbor6(3)+bcinfo%dir(3)*k
            pv = pv + variable%getpv(ii,jj,kk)
            tv = tv + variable%gettv(ii,jj,kk)
            pv = 1.d0/7.d0*pv
            tv = 1.d0/7.d0*tv
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+bcinfo%pv(1),pv(5),pv(6),pv(7),dv)
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do

    end subroutine cornernowall
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcwallinviscid(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),pv_s(variable%getnpv())
      real(8) :: dv(variable%getndv()),tv(variable%getntv())
      real(8) :: nx(3),var
  
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getcx(ii-1,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getcx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getex(ii,jj-1,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getex(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              nx = grid%gettx(ii,jj,kk-1)
              pv_s = variable%getpv(ii,jj,kk)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              nx = - grid%gettx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            end select
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            do m=1,variable%getnpv()
              select case(m)
              case(2)
                var = 2.d0*(pv_s(2)-nx(1)*(pv_s(2)*nx(1)+pv_s(3)*nx(2)+pv_s(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)) - pv(2)
                call variable%setpv(m,i,j,k,var)
              case(3)
                var = 2.d0*(pv_s(3)-nx(2)*(pv_s(2)*nx(1)+pv_s(3)*nx(2)+pv_s(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)) - pv(3)
                call variable%setpv(m,i,j,k,var)
              case(4)
                var = 2.d0*(pv_s(4)-nx(3)*(pv_s(2)*nx(1)+pv_s(3)*nx(2)+pv_s(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)) - pv(4)
                call variable%setpv(m,i,j,k,var)
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
      
    end subroutine bcwallinviscid
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcwallviscous(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4)
                call variable%setpv(m,i,j,k,-pv(m))
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
   
    end subroutine bcwallviscous
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcwallviscouske(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4,8,9)
                call variable%setpv(m,i,j,k,-pv(m))
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              select case(m)
              case(3)
                call variable%settv(m,i,j,k,-tv(m))
              case default
                call variable%settv(m,i,j,k,tv(m))
              end select
            end do
          end do
        end do
      end do

    end subroutine bcwallviscouske
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcwallviscouskw(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv()),grd(grid%getngrd())
      real(8) :: dv_b(variable%getndv()),tv_b(variable%getntv())
      real(8) :: var

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
            end select
            dv_b = variable%getdv(ii,jj,kk)
            tv_b = variable%gettv(ii,jj,kk)
            grd = grid%getgrd(ii,jj,kk)
            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4,8)
                call variable%setpv(m,i,j,k,-pv(m))
              case(9)
                var = 1600.d0*tv_b(1)/dv_b(1)/grd(5)**2-pv(m)
                call variable%setpv(m,i,j,k,var)
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              select case(m)
              case(3)
                call variable%settv(m,i,j,k,-tv(m))
              case default
                call variable%settv(m,i,j,k,tv(m))
              end select
            end do
          end do
        end do
      end do     
    end subroutine bcwallviscouskw
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcinflow(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k
    
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            
          end do
        end do
      end do
     
    end subroutine bcinflow
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcinflowsubsonic(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: pv_s(variable%getnpv()),nx(3),var,dl,a

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getcx(ii-1,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getcx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getex(ii,jj-1,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getex(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              nx = - grid%gettx(ii,jj,kk-1)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              nx = grid%gettx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            end select
            dl = 1.d0/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            a = prec%getsndp2(dv(6),pv_s(2)**2+pv_s(3)**2+pv_s(4)**2)
            a = dsqrt(a)
            var = pv_s(1) - dv(1)*a*(nx(1)*(bcinfo%pv(2)-pv_s(2))  + nx(2)*(bcinfo%pv(3)-pv_s(3)) &
                                        + nx(3)*(bcinfo%pv(4)-pv_s(4)))*dl
            pv(1) = var - pv(1)
            pv(2) = 2.d0*( bcinfo%pv(2)+nx(1)*0.5d0*var/dv(1)/a*dl) - pv(2)
            pv(3) = 2.d0*( bcinfo%pv(3)+nx(2)*0.5d0*var/dv(1)/a*dl) - pv(3)
            pv(4) = 2.d0*( bcinfo%pv(4)+nx(3)*0.5d0*var/dv(1)/a*dl) - pv(4)
            pv(5:variable%getnpv()) = 2.d0*bcinfo%pv(5:variable%getnpv()) - pv(5:variable%getnpv())
            pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
            pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+bcinfo%pv(1),pv(5),pv(6),pv(7),dv)
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
     
    end subroutine bcinflowsubsonic
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcinflowsupersonic(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      integer :: i,j,k,ii,jj,kk,m
      
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            pv(1) = -pv(1)
            pv(2:variable%getnpv()) = 2.d0*bcinfo%pv(2:variable%getnpv()) - pv(2:variable%getnpv())
            pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
            pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+bcinfo%pv(1),pv(5),pv(6),pv(7),dv)
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
     
    end subroutine bcinflowsupersonic
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcoutflow(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k
      
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            
          end do
        end do
      end do
      
    end subroutine bcoutflow
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcoutflowsubsonic(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: pv_s(variable%getnpv()),nx(3),dl,a
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getcx(ii-1,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getcx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getex(ii,jj-1,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getex(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              nx = - grid%gettx(ii,jj,kk-1)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              nx = grid%gettx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
            end select
            dl = 1.d0/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            a = prec%getsndp2(dv(6),pv_s(2)**2+pv_s(3)**2+pv_s(4)**2)
            a = dsqrt(a)
            pv(1) = -pv(1)
            pv(2) = 2.d0*( pv_s(2)+nx(1)*pv_s(1)/dv(1)/a*dl) - pv(2)
            pv(3) = 2.d0*( pv_s(3)+nx(2)*pv_s(1)/dv(1)/a*dl) - pv(3)
            pv(4) = 2.d0*( pv_s(4)+nx(3)*pv_s(1)/dv(1)/a*dl) - pv(4)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+bcinfo%pv(1),pv(5),pv(6),pv(7),dv)
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
     
    end subroutine bcoutflowsubsonic
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcoutflowsupersonic(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    end subroutine bcoutflowsupersonic
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcfarfield(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: nx(3),vel,mach,dl,var,a
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: pv_b(variable%getnpv()),dv_b(variable%getndv())

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getcx(ii-1,jj,kk)
              pv_b = variable%getpv(ii,jj,kk)
              dv_b = variable%getdv(ii,jj,kk)
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getcx(ii,jj,kk)
              pv_b = variable%getpv(ii,jj,kk)
              dv_b = variable%getdv(ii,jj,kk)
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getex(ii,jj-1,kk)
              pv_b = variable%getpv(ii,jj,kk)
              dv_b = variable%getdv(ii,jj,kk)
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getex(ii,jj,kk)
              pv_b = variable%getpv(ii,jj,kk)
              dv_b = variable%getdv(ii,jj,kk)
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              nx = - grid%gettx(ii,jj,kk-1)
              pv_b = variable%getpv(ii,jj,kk)
              dv_b = variable%getdv(ii,jj,kk)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              nx = grid%gettx(ii,jj,kk)
              pv_b = variable%getpv(ii,jj,kk)
              dv_b = variable%getdv(ii,jj,kk)
            end select
            dl = 1.d0/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            vel = (pv_b(2)*nx(1)+pv_b(3)*nx(2)+pv_b(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
            mach = vel/dsqrt(dv_b(6))
            a = prec%getsndp2(dv_b(6),pv_b(2)**2+pv_b(3)**2+pv_b(4)**2)
            a = dsqrt(a)
            if(mach.ge.0.d0) then !out
              if(mach.ge.1.d0) then !super                
              else !sub
                pv(1) = -pv(1)
                pv(2) = 2.d0*( pv_b(2)+nx(1)*pv_b(1)/dv_b(1)/a*dl) - pv(2)
                pv(3) = 2.d0*( pv_b(3)+nx(2)*pv_b(1)/dv_b(1)/a*dl) - pv(3)
                pv(4) = 2.d0*( pv_b(4)+nx(3)*pv_b(1)/dv_b(1)/a*dl) - pv(4)
              end if
            else ! in
              if(mach.le.-1.d0) then !super
                pv(1) = -pv(1)
                pv(2:variable%getnpv()) = 2.d0*bcinfo%pv(2:variable%getnpv()) - pv(2:variable%getnpv())
                pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
                pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
              else !sub
                var = pv_b(1) - dv_b(1)*a*(nx(1)*(bcinfo%pv(2)-pv_b(2)) + nx(2)*(bcinfo%pv(3)-pv_b(3)) &
                                              + nx(3)*(bcinfo%pv(4)-pv_b(4)))*dl
                pv(1) = var - pv(1)
                pv(2) = 2.d0*( bcinfo%pv(2)+nx(1)*0.5d0*var/dv_b(1)/a*dl) - pv(2)
                pv(3) = 2.d0*( bcinfo%pv(3)+nx(2)*0.5d0*var/dv_b(1)/a*dl) - pv(3)
                pv(4) = 2.d0*( bcinfo%pv(4)+nx(3)*0.5d0*var/dv_b(1)/a*dl) - pv(4)
                pv(5:variable%getnpv()) = 2.d0*bcinfo%pv(5:variable%getnpv()) - pv(5:variable%getnpv())
                pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
                pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
              end if
            end if
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+bcinfo%pv(1),pv(5),pv(6),pv(7),dv)
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    
    end subroutine bcfarfield
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcdegeneratepoint(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    
    end subroutine bcdegeneratepoint
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcsymmetryplane(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),pv_s(variable%getnpv())
      real(8) :: dv(variable%getndv()),tv(variable%getntv())
      real(8) :: nx(3),var
  
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getcx(ii-1,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getcx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getex(ii,jj-1,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getex(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              nx = grid%gettx(ii,jj,kk-1)
              pv_s = variable%getpv(ii,jj,kk)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              nx = - grid%gettx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)
            end select
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            do m=1,variable%getnpv()
              select case(m)
              case(2)
                var = 2.d0*(pv_s(2)-nx(1)*(pv_s(2)*nx(1)+pv_s(3)*nx(2)+pv_s(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)) - pv(2)
                call variable%setpv(m,i,j,k,var)
              case(3)
                var = 2.d0*(pv_s(3)-nx(2)*(pv_s(2)*nx(1)+pv_s(3)*nx(2)+pv_s(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)) - pv(3)
                call variable%setpv(m,i,j,k,var)
              case(4)
                var = 2.d0*(pv_s(4)-nx(3)*(pv_s(2)*nx(1)+pv_s(3)*nx(2)+pv_s(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)) - pv(4)
                call variable%setpv(m,i,j,k,var)
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do
    end subroutine bcsymmetryplane
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcshiftedperiodic(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
 
      ! h=1/16
      select case(bcinfo%face)
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)-16
              pv = variable%getpv(i+16,j+16,k)
              dv = variable%getdv(i+16,j+16,k)
              tv = variable%gettv(i+16,j+16,k)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
            do i=bcinfo%iend(1)-15,bcinfo%iend(1)
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              pv = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
              tv = variable%gettv(ii,jj,kk)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do

      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)        
            do i=bcinfo%istart(1)+16,bcinfo%iend(1)
              pv = variable%getpv(i-16,j-16,k)
              dv = variable%getdv(i-16,j-16,k)
              tv = variable%gettv(i-16,j-16,k)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
            do i=bcinfo%istart(1),bcinfo%istart(1)+15
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              pv = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
              tv = variable%gettv(ii,jj,kk)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do

      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)-16
              pv = variable%getpv(i+16,j,k+16)
              dv = variable%getdv(i+16,j,k+16)
              tv = variable%gettv(i+16,j,k+16)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
            do i=bcinfo%iend(1)-15,bcinfo%iend(1)
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              pv = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
              tv = variable%gettv(ii,jj,kk)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do

      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)        
            do i=bcinfo%istart(1)+16,bcinfo%iend(1)
              pv = variable%getpv(i-16,j,k-16)
              dv = variable%getdv(i-16,j,k-16)
              tv = variable%gettv(i-16,j,k-16)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
            do i=bcinfo%istart(1),bcinfo%istart(1)+15
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              pv = variable%getpv(ii,jj,kk)
              dv = variable%getdv(ii,jj,kk)
              tv = variable%gettv(ii,jj,kk)
              do m=1,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do

      case default
      end select

    end subroutine bcshiftedperiodic
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bctotalpressurein(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: grd(grid%getngrd()),nx(3),gridvel(3)
      real(8) :: ua,va,wa,area,uvwa2,x(2),h_total,pa,ta
      real(8),parameter :: pi = 4.d0*datan(1.d0)
      logical :: check

      ua = 0.d0
      va = 0.d0
      wa = 0.d0
      pa = 0.d0
      ta = 0.d0
      area = 0.d0

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            select case(bcinfo%face)
            case('imin')
              if(i.eq.bcinfo%istart(1)) then
                ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
                jj = bcinfo%origin(2)+bcinfo%dir(2)*j
                kk = bcinfo%origin(3)+bcinfo%dir(3)*k
                nx = - grid%getcx(ii-1,jj,kk)
                grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii-1,jj,kk))

                gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
                gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
                gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)

                ua = ua + (pv(2)+gridvel(1))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                va = va + (pv(3)+gridvel(2))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                wa = wa + (pv(4)+gridvel(3))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                pa = pa + pv(1)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                ta = ta + pv(5)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                area = area + dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              end if
            case('imax')
              if(i.eq.bcinfo%iend(1)) then
                ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
                jj = bcinfo%origin(2)+bcinfo%dir(2)*j
                kk = bcinfo%origin(3)+bcinfo%dir(3)*k
                nx = grid%getcx(ii,jj,kk)
                grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii+1,jj,kk))

                gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
                gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
                gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)

                ua = ua + (pv(2)+gridvel(1))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                va = va + (pv(3)+gridvel(2))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                wa = wa + (pv(4)+gridvel(3))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                pa = pa + pv(1)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                ta = ta + pv(5)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                area = area + dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              end if
            case('jmin')
              if(j.eq.bcinfo%istart(2)) then
                ii = bcinfo%origin(1)+bcinfo%dir(1)*i
                jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
                kk = bcinfo%origin(3)+bcinfo%dir(3)*k
                nx = - grid%getex(ii,jj-1,kk)
                grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj-1,kk))

                gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
                gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
                gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)

                ua = ua + (pv(2)+gridvel(1))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                va = va + (pv(3)+gridvel(2))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                wa = wa + (pv(4)+gridvel(3))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                pa = pa + pv(1)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                ta = ta + pv(5)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                area = area + dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              end if
            case('jmax')
              if(j.eq.bcinfo%iend(2)) then
                ii = bcinfo%origin(1)+bcinfo%dir(1)*i
                jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
                kk = bcinfo%origin(3)+bcinfo%dir(3)*k
                nx = grid%getex(ii,jj,kk)
                grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj+1,kk))

                gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
                gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
                gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)

                ua = ua + (pv(2)+gridvel(1))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                va = va + (pv(3)+gridvel(2))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                wa = wa + (pv(4)+gridvel(3))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                pa = pa + pv(1)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                ta = ta + pv(5)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                area = area + dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              end if
            case('kmin')
              if(k.eq.bcinfo%istart(3)) then
                ii = bcinfo%origin(1)+bcinfo%dir(1)*i
                jj = bcinfo%origin(2)+bcinfo%dir(2)*j
                kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
                nx = - grid%gettx(ii,jj,kk-1)
                grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj,kk-1))

                gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
                gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
                gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)

                ua = ua + (pv(2)+gridvel(1))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                va = va + (pv(3)+gridvel(2))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                wa = wa + (pv(4)+gridvel(3))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                pa = pa + pv(1)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                ta = ta + pv(5)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                area = area + dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              end if
            case('kmax')
              if(k.eq.bcinfo%iend(3)) then
                ii = bcinfo%origin(1)+bcinfo%dir(1)*i
                jj = bcinfo%origin(2)+bcinfo%dir(2)*j
                kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
                nx = grid%gettx(ii,jj,kk)
                grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj,kk+1))

                gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
                gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
                gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)

                ua = ua + (pv(2)+gridvel(1))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                va = va + (pv(3)+gridvel(2))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                wa = wa + (pv(4)+gridvel(3))*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                pa = pa + pv(1)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                ta = ta + pv(5)*dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
                area = area + dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              end if
            end select
          end do
        end do
      end do

      ua = ua/area
      va = va/area
      pa = pa/area
      ta = ta/area
      uvwa2 = ua**2+va**2+wa**2

      ! reference pv(6), pv(7) should be given as boundary condition
      call eos%deteos(bcinfo%total_pressure,bcinfo%total_temperature,bcinfo%pv(6),bcinfo%pv(7),dv)
      h_total = dv(2)

      ! initial guess
      ! x(1)=static pressure, x(2)=static temperature
      x(1) = pa+bcinfo%pv(1)
      x(2) = ta

      call newt(eos,2,bcinfo%total_pressure,h_total,bcinfo%pv(6),bcinfo%pv(7),uvwa2,x,check)

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)

            pv(1) = 2.d0*(x(1)-bcinfo%pv(1))-pv(1)
            pv(5) = 2.d0*x(2)-pv(5)
            pv(6) = 2.d0*bcinfo%pv(6)-pv(6)
            pv(7) = 2.d0*bcinfo%pv(7)-pv(7)
            pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
            pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)

            do m=1,7
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+bcinfo%pv(1),pv(4),pv(6),pv(7),dv)
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do

    end subroutine bctotalpressurein
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bccounterrotatingwallviscous(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: grd(grid%getngrd()),gridvel(3)

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)

            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii-1,jj,kk))

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii+1,jj,kk))

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('jmin')
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj-1,kk))

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('jmax')
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj+1,kk))

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('kmin')
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj,kk-1))

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('kmax')
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj,kk+1))

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            end select

            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4)
                call variable%setpv(m,i,j,k,-pv(m)+2.d0*gridvel(m-1))
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do


    end subroutine bccounterrotatingwallviscous
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bccounterrotatingwallinviscid(bcinfo,grid,variable,eos,prop,prec)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      type(t_prec), intent(in) :: prec
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),pv_s(variable%getnpv())
      real(8) :: dv(variable%getndv()),tv(variable%getntv())
      real(8) :: grd(grid%getngrd()),gridvel(3),nx(3),var

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            dv = variable%getdv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)

            select case(bcinfo%face)
            case('imin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%iend(1)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii-1,jj,kk))
              nx = grid%getcx(ii-1,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii+1,jj,kk))
              nx = -grid%getcx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('jmin')
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj-1,kk))
              nx = grid%getex(ii,jj-1,kk)
              pv_s = variable%getpv(ii,jj,kk)

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('jmax')
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj+1,kk))
              nx = -grid%getex(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('kmin')
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj,kk-1))
              nx = grid%gettx(ii,jj,kk-1)
              pv_s = variable%getpv(ii,jj,kk)

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            case('kmax')
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              grd = 0.5d0*(grid%getgrd(ii,jj,kk)+grid%getgrd(ii,jj,kk+1))
              nx = -grid%gettx(ii,jj,kk)
              pv_s = variable%getpv(ii,jj,kk)

              gridvel(1) = bcinfo%omega(2)*grd(4)-bcinfo%omega(3)*grd(3)
              gridvel(2) = bcinfo%omega(3)*grd(2)-bcinfo%omega(1)*grd(4)
              gridvel(3) = bcinfo%omega(1)*grd(3)-bcinfo%omega(2)*grd(2)
            end select

            do m=1,variable%getnpv()
              select case(m)
              case(2,3,4)
                var = 2.d0*(pv_s(m) - nx(m-1)*(pv_s(2)*nx(1)+pv_s(3)*nx(2)+pv_s(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2) + 2.d0*gridvel(m-1)) - pv(m)
                call variable%setpv(m,i,j,k,var)
              case default
                call variable%setpv(m,i,j,k,pv(m))
              end select
            end do
            do m=1,variable%getndv()
              call variable%setdv(m,i,j,k,dv(m))
            end do
            do m=1,variable%getntv()
              call variable%settv(m,i,j,k,tv(m))
            end do
          end do
        end do
      end do

    end subroutine bccounterrotatingwallinviscid
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function no_prec(prec,snd2,uuu2) result(sndp2)
      implicit none
      class(t_prec), intent(in) :: prec
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
      
      sndp2 = snd2
      
    end function no_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function steady_prec(prec,snd2,uuu2) result(sndp2)
      implicit none
      class(t_prec), intent(in) :: prec
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2
      
      sndp2 = dmin1(snd2,dmax1(uuu2,prec%uref**2))
      
    end function steady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function unsteady_prec(prec,snd2,uuu2) result(sndp2)
      implicit none
      class(t_prec), intent(in) :: prec
      real(8), intent(in) :: snd2,uuu2
      real(8) :: sndp2

      sndp2 = dmin1(snd2,dmax1(uuu2,prec%uref**2,prec%str**2*uuu2))
      
    end function unsteady_prec
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine newt(eos,n,p_total,h_total,y1,y2,uvwa2,x,check)
      implicit none
      type(t_eos), intent(in) :: eos
      integer, intent(in) :: n
      real(8), intent(in) :: p_total,h_total,y1,y2,uvwa2
      logical, intent(out) :: check
      real(8), intent(inout) :: x(n)
      integer, parameter :: maxits=200
      real(8), parameter :: tolf=1.d-4,tolmin=1.d-6,tolx=1.d-7,stpmx=100.d0
      integer :: i,its,j,indx(n)
      real(8) :: fvec(n)
      real(8) :: d,den,f,fold,stpmax,sum,temp,test
      real(8) :: fjac(n,n),g(n),p(n),xold(n)

      f=fmin(eos,n,p_total,h_total,y1,y2,uvwa2,x,fvec)
      test=0.d0

      do i=1,n
        if(dabs(fvec(i)).gt.test) test=dabs(fvec(i))
      end do
      if(test.lt.0.01d0*tolf) then
        check=.false.
        return
      end if

      sum=0.d0
      do i=1,n
        sum=sum+x(i)**2
      end do
      stpmax=stpmx*dmax1(dsqrt(sum),dble(n))
      do its=1,maxits
        call fdjac(eos,n,p_total,h_total,y1,y2,uvwa2,x,fvec,fjac)
        do i=1,n
          sum=0.d0
          do j=1,n
            sum=sum+fjac(j,i)*fvec(j)
          end do
          g(i)=sum
        end do
        do i=1,n
          xold(i)=x(i)
        end do
        fold=f
        do i=1,n
          p(i)=-fvec(i)
        end do

        call ludcmp(fjac,n,indx,d)
        call lubksb(fjac,n,indx,p)
        call lnsrch(eos,n,p_total,h_total,y1,y2,uvwa2,xold,fold,g,p,x,f,stpmax,check)
        test=0.d0
        do i=1,n
          if(dabs(fvec(i)).gt.test) test=dabs(fvec(i))
        end do
        if (test.lt.tolf) then
          check=.false.
          return
        end if
        if(check) then
          test=0.d0
          den=dmax1(f,0.5d0*n)
          do i=1,n
            temp=dabs(g(i))*dmax1(dabs(x(i)),1.d0)/den
            if(temp.gt.test) test=temp
          end do
          if(test.lt.tolmin) then
            check=.true.
          else
            check=.false.
          end if
          return
        end if
        test=0.d0
        do i=1,n
          temp=(dabs(x(i)-xold(i)))/dmax1(dabs(x(i)),1.d0)
          if(temp.gt.test) test=temp
        end do
        if(test.lt.tolx) return
      end do

      write(*,*) "maxits exceeded in newt"

    end subroutine newt
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function fmin(eos,n,p_total,h_total,y1,y2,uvwa2,x,fvec)
      implicit none
      type(t_eos), intent(in) :: eos
      integer, intent(in) :: n
      real(8), intent(in) :: h_total,p_total,y1,y2,uvwa2,x(n)
      real(8), intent(out) :: fvec(n)
      integer :: i
      real(8) :: fmin,sum

      call funcv(eos,n,p_total,h_total,y1,y2,uvwa2,x,fvec)

      sum=0.d0
      do i=1,n
        sum=sum+fvec(i)**2
      enddo
      fmin=0.5d0*sum
      return
    end function fmin
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine funcv(eos,n,p_total,h_total,y1,y2,uvwa2,x,fvec)
      implicit none
      type(t_eos), intent(in) :: eos
      integer, intent(in) :: n
      real(8), intent(in) :: p_total,h_total,y1,y2,uvwa2,x(n)
      real(8), intent(out) :: fvec(n)
      real(8) :: dv(18)

      call eos%deteos(x(1),x(2),y1,y2,dv)

      fvec(1) = p_total - (x(1)+0.5d0*dv(1)*uvwa2)
      fvec(2) = h_total - (dv(2)+0.5d0*uvwa2)
    end subroutine funcv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine fdjac(eos,n,p_total,h_total,y1,y2,uvwa2,x,fvec,df)
      implicit none
      type(t_eos), intent(in) :: eos
      integer, intent(in) :: n
      real(8), intent(in) :: p_total,h_total,y1,y2,uvwa2,x(n),fvec(n)
      real(8), intent(out) :: df(n,n)
      real(8),parameter :: eps=1.0d-4
      integer :: i,j
      real(8) :: h,temp(n),f(n)

      do j=1,n
        temp = x
        h=eps*dabs(temp(j))
        if(h.eq.0.0d0) h=eps
        temp(j)=temp(j)+h
        h=temp(j)-x(j)
        call funcv(eos,n,p_total,h_total,y1,y2,uvwa2,temp,f)
        do i=1,n
          df(i,j)=(f(i)-fvec(i))/h
        enddo
      enddo

      return
    end subroutine fdjac
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine lnsrch(eos,n,p_total,h_total,y1,y2,uvwa2,xold,fold,g,p,x,f,stpmax,check)
      implicit none
      type(t_eos), intent(in) :: eos
      integer, intent(in) :: n
      real(8), intent(in) :: p_total,h_total,y1,y2,uvwa2
      logical :: check
      real(8) :: f,fold,stpmax,g(n),p(n),x(n),xold(n),fvec(n)
      real(8), parameter :: alf=1.d-4,tolx=1.d-7
      integer :: i
      real(8) :: a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2
      real(8) :: slope,sum,temp,test,tmplam

      check=.false.
      sum=0.0d0
      do i=1,n
        sum=sum+p(i)*p(i)
      enddo
      sum=dsqrt(sum)
      if(sum.gt.stpmax) then
        do i=1,n
          p(i)=p(i)*stpmax/sum
         enddo
      endif
      slope=0.0d0
      do i=1,n
        slope=slope+g(i)*p(i)
      enddo
      if(slope.ge.0.d0) write(*,*) "roundoff problem in lnsrch"
      test=0.0d0
      do i=1,n
        temp=dabs(p(i))/dmax1(dabs(xold(i)),1.0d0)
        if(temp.gt.test) test=temp
      enddo
      alamin=tolx/test
      alam=1.0d0
1     continue
      do i=1,n
        x(i)=xold(i)+alam*p(i)
      enddo
      f=fmin(eos,n,p_total,h_total,y1,y2,uvwa2,x,fvec)
      if(alam.lt.alamin) then
        do i=1,n
          x(i)=xold(i)
        enddo
        check=.true.
        return
      else if(f.le.fold+alf*alam*slope) then
        return
      else
        if(alam.eq.1.0d0) then
          tmplam=-slope/(2.0d0*(f-fold-slope))
        else
          rhs1=f-fold-alam*slope
          rhs2=f2-fold2-alam2*slope
          a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
          b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
          if(a.eq.0.0d0) then
            tmplam=-slope/(2.0d0*b)
          else
            disc=b*b-3.0d0*a*slope
            if(disc.lt.0.d0) then
              tmplam=0.5d0*alam
            else if(b.le.0.d0) then
              tmplam=(-b+dsqrt(disc))/(3.0d0*a)
            else
              tmplam=-slope/(b+dsqrt(disc))
            end if
          endif
          if(tmplam.gt.0.5d0*alam) tmplam=0.5d0*alam
        endif
      endif
      alam2=alam
      f2=f
      fold2=fold
      alam=dmax1(tmplam,0.1d0*alam)
      goto 1

    end subroutine lnsrch
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine ludcmp(a,n,indx,d)
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: indx(n)
      real(8), intent(out) :: d
      real(8), intent(inout) :: a(n,n)
      real(8), parameter :: tiny=1.d-20
      integer :: i,imax,j,k
      real(8) :: aamax,dum,sum,vv(n)

      d=1.d0
      do i=1,n
        aamax=0.d0
        do j=1,n
          if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
        end do
        if(aamax.eq.0.d0) write(*,*) "singular matrix in ludcmp"
        vv(i)=1.d0/aamax
      end do

      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          end do
          a(i,j)=sum
        end do
        aamax=0.d0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          end do
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if(dum.ge.aamax) then
            imax=i
            aamax=dum
          end if
        end do
        if(j.ne.imax) then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          end do
          d=-d
          vv(imax)=vv(j)
        end if
        indx(j)=imax
        if(a(j,j).eq.0.d0) a(j,j)=tiny
        if(j.ne.n) then
          dum=1.d0/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          end do
        end if
      end do
      return

    end subroutine ludcmp
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine lubksb(a,n,indx,b)
      implicit none
      integer, intent(in) :: n,indx(n)
      real(8), intent(in) :: a(n,n)
      real(8), intent(inout) :: b(n)
      integer :: i,ii,j,ll
      real(8) :: sum

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          end do
        else if(sum.ne.0.d0) then
          ii=i
        end if
        b(i)=sum
      end do
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        end do
        b(i)=sum/a(i,i)
      end do
      return

    end subroutine lubksb
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module bc_module
