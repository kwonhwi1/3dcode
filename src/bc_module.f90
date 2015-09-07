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
    procedure(p_bctype), pointer :: bctype
  end type t_bcinfo2
  
  type t_ref
    real(8), dimension(:), allocatable :: pv,tv,dv
  end type t_ref
  
  type t_bc
    private
    integer :: rank,size,iturb,nbc,ncon
    integer :: npv,ntv,ndv,ngrd
    type(t_ref) :: ref
    type(t_bcinfo2), dimension(:), allocatable :: bcinfo,corner,edge
    type(t_connectinfo), dimension(:), allocatable :: connectinfo
    type(t_mpitemp), dimension(:), allocatable :: mpitemp
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setbc
  end type t_bc

  interface
    subroutine p_bctype(bcinfo,ref,grid,variable,eos,prop)
      import t_bcinfo2
      import t_ref
      import t_grid
      import t_variable
      import t_eos
      import t_prop  
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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

      allocate(bc%ref%pv(bc%npv),bc%ref%dv(bc%ndv),bc%ref%tv(bc%ntv))
      
      bc%ref%pv(1) = config%getpref()
      bc%ref%pv(2) = config%geturef()*dcos(config%getaos())*dcos(config%getaoa())
      bc%ref%pv(3) = config%geturef()*dcos(config%getaos())*dsin(config%getaoa())
      bc%ref%pv(4) = config%geturef()*dsin(config%getaos())
      bc%ref%pv(5) = config%gettref()
      bc%ref%pv(6) = config%gety1ref()
      bc%ref%pv(7) = config%gety2ref()
      
      call eos%deteos(bc%ref%pv(1),bc%ref%pv(5),bc%ref%pv(6),bc%ref%pv(7),bc%ref%dv)
      
      if(config%getiturb().ge.-2) then
        call prop%detprop(bc%ref%dv(3),bc%ref%dv(4),bc%ref%dv(5),bc%ref%pv(5),bc%ref%pv(6),bc%ref%pv(7),bc%ref%tv(1:2))
        if(config%getiturb().ge.-1) then
          bc%ref%tv(3) = config%getemutref()
          bc%ref%pv(8) = config%getkref()
          bc%ref%pv(9) = config%getoref()
        end if
      end if
 
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
        else if(trim(bc%bcinfo(n)%bcname).eq.'UserDefined') then
          bc%bcinfo(n)%bctype => bcshiftedperiodic
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
      end do

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
      end do

      isurf_edge = .false.; jsurf_edge = .false.; ksurf_edge = .false.

      do m=1,4
        do n=1,bc%nbc 
          if((trim(bc%bcinfo(n)%bcname).eq.'BCWall').and. &
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
      
      do m=1,8
        bc%corner(m)%dir = 0
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
      
      do n=1,bc%nbc
        if(associated(bc%bcinfo(n)%bctype)) nullify(bc%bcinfo(n)%bctype)
      end do
      
      do n=1,12
        if(associated(bc%edge(n)%bctype)) nullify(bc%edge(n)%bctype)
      end do
      
      do n=1,8
        if(associated(bc%corner(n)%bctype)) nullify(bc%corner(n)%bctype)
      end do

      do n=1,bc%ncon
        if(allocated(bc%mpitemp(n)%sendbuf)) deallocate(bc%mpitemp(n)%sendbuf)
        if(allocated(bc%mpitemp(n)%recvbuf)) deallocate(bc%mpitemp(n)%recvbuf)
      end do
      
      deallocate(bc%ref%pv,bc%ref%dv,bc%ref%tv)
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
      integer :: status(mpi_status_size,bc%ncon)
      real(8) :: pv(bc%npv),dv(bc%ndv),tv(bc%ntv)
      
      do n=1,bc%nbc
        call bc%bcinfo(n)%bctype(bc%ref,grid,variable,eos,prop)
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
      
      if(bc%size.gt.1) then
        call mpi_waitall(bc%ncon,request_s ,status,ier)
        call mpi_waitall(bc%ncon,request_sa,status,ier)
        call mpi_waitall(bc%ncon,request_r ,status,ier)
        call mpi_waitall(bc%ncon,request_ra,status,ier)
      end if

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
        call bc%edge(n)%bctype(bc%ref,grid,variable,eos,prop)
      end do

      do n=1,8
        call bc%corner(n)%bctype(bc%ref,grid,variable,eos,prop)
      end do
      
    end subroutine setbc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine edgewallwall(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine edgewallwallkw(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine cornerwallwall(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine cornerwallwallkw(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine edgenowall(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
            call eos%deteos(pv(1)+ref%pv(1),pv(5),pv(6),pv(7),dv) 
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
    subroutine cornernowall(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
            call eos%deteos(pv(1)+ref%pv(1),pv(5),pv(6),pv(7),dv) 
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
    subroutine bcwallinviscid(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
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
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getcx(ii,jj,kk)
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getex(ii,jj-1,kk)
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getex(ii,jj,kk)
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              nx = grid%gettx(ii,jj,kk-1)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              nx = - grid%gettx(ii,jj,kk)
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
                var = pv(2)-2.d0*nx(1)*(pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)
                call variable%setpv(m,i,j,k,var)
              case(3)
                var = pv(3)-2.d0*nx(2)*(pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)
                call variable%setpv(m,i,j,k,var)
              case(4)
                var = pv(4)-2.d0*nx(3)*(pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)
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
    subroutine bcwallviscous(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine bcwallviscouske(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine bcwallviscouskw(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine bcinflow(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k
    
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            
          end do
        end do
      end do
     
    end subroutine bcinflow
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcinflowsubsonic(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())

      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            pv(2:variable%getnpv()) = 2.d0*ref%pv(2:variable%getnpv()) - pv(2:variable%getnpv())
            pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
            pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+ref%pv(1),pv(5),pv(6),pv(7),dv)
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
    subroutine bcinflowsupersonic(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
            pv(2:variable%getnpv()) = 2.d0*ref%pv(2:variable%getnpv()) - pv(2:variable%getnpv())
            pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
            pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+ref%pv(1),pv(5),pv(6),pv(7),dv)
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
    subroutine bcoutflow(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k
      
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            
          end do
        end do
      end do
      
    end subroutine bcoutflow
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bcoutflowsubsonic(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
    
      do k=bcinfo%istart(3),bcinfo%iend(3)
        do j=bcinfo%istart(2),bcinfo%iend(2)
          do i=bcinfo%istart(1),bcinfo%iend(1)
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            pv(1) = -pv(1)
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+ref%pv(1),pv(5),pv(6),pv(7),dv)
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
    subroutine bcoutflowsupersonic(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine bcfarfield(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,ii,jj,kk,m
      real(8) :: nx(3),vel,mach
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
            ii = bcinfo%origin(1)+bcinfo%dir(1)*i
            jj = bcinfo%origin(2)+bcinfo%dir(2)*j
            kk = bcinfo%origin(3)+bcinfo%dir(3)*k
            pv = variable%getpv(ii,jj,kk)
            tv = variable%gettv(ii,jj,kk)
            vel = (pv_b(2)*nx(1)+pv_b(3)*nx(2)+pv_b(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
            mach = vel/dsqrt(dv_b(6))
            if(mach.ge.0.d0) then !out
              if(mach.ge.1.d0) then !super                
              else !sub
                pv(1) = -pv(1)
              end if
            else ! in
              if(mach.le.-1.d0) then !super
                pv(1) = -pv(1)
                pv(2:variable%getnpv()) = 2.d0*ref%pv(2:variable%getnpv()) - pv(2:variable%getnpv())
                pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
                pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
              else !sub
                pv(2:variable%getnpv()) = 2.d0*ref%pv(2:variable%getnpv()) - pv(2:variable%getnpv())
                pv(6) = dmin1(dmax1(pv(6),0.d0),1.d0)
                pv(7) = dmin1(dmax1(pv(7),0.d0),1.d0)
              end if
            end if
            do m=1,variable%getnpv()
              call variable%setpv(m,i,j,k,pv(m))
            end do
            call eos%deteos(pv(1)+ref%pv(1),pv(5),pv(6),pv(7),dv)
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
    subroutine bcdegeneratepoint(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
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
    subroutine bcsymmetryplane(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
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
            case('imax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*bcinfo%istart(1)
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getcx(ii,jj,kk)
            case('jmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%iend(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = grid%getex(ii,jj-1,kk)
            case('jmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*bcinfo%istart(2)
              kk = bcinfo%origin(3)+bcinfo%dir(3)*k
              nx = - grid%getex(ii,jj,kk)
            case('kmin')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%iend(3)
              nx = grid%gettx(ii,jj,kk-1)
            case('kmax')
              ii = bcinfo%origin(1)+bcinfo%dir(1)*i
              jj = bcinfo%origin(2)+bcinfo%dir(2)*j
              kk = bcinfo%origin(3)+bcinfo%dir(3)*bcinfo%istart(3)
              nx = - grid%gettx(ii,jj,kk)
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
                var = pv(2)-2.d0*nx(1)*(pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)
                call variable%setpv(m,i,j,k,var)
              case(3)
                var = pv(3)-2.d0*nx(2)*(pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)
                call variable%setpv(m,i,j,k,var)
              case(4)
                var = pv(4)-2.d0*nx(3)*(pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/(nx(1)**2+nx(2)**2+nx(3)**2)
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
    subroutine bcshiftedperiodic(bcinfo,ref,grid,variable,eos,prop)
      implicit none
      class(t_bcinfo2), intent(in) :: bcinfo
      type(t_ref), intent(in) :: ref
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,m,ii,jj,kk
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      real(8) :: nx(3),var
 
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
end module bc_module
