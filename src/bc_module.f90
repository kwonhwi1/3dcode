module bc_module
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
    character(4) :: face
    procedure(p_bctype), pointer :: bctype
  end type t_bcinfo2
  
  type t_ref
    real(8), dimension(:), allocatable :: pv,tv,dv
  end type t_ref
  
  type t_mpitemp
    integer :: num
    integer :: sendadress(21),recvadress(21)
    real(8), dimension(:), allocatable :: sendbuf,recvbuf
  end type t_moitemp

  type t_bc
    private
    integer :: rank,size,iturb,nbc,ncon
    integer :: npv,ntv,ndv,ngrd
    type(t_ref) :: ref
    type(t_bcinfo2), dimension(:), allocatable :: bcinfo
    type(t_bcinfo2) :: corner(8),edge(12)
    type(t_connectinfo), dimension(:), allocatable :: connectinfo
    type(t_mpitemp), dimension(:), allocatable :: mpitemp
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setbc
      procedure, private :: bccorner
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
      bc%ref%pv(2) = config%geturef()*dcos(config%getaoa())
      bc%ref%pv(3) = config%geturef()*dsin(config%getaoa())
      bc%ref%pv(4) = 0.d0
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
        
        if(bc%bcinfo(n)%istart(2).eq.bc%bcinfo(n)%iend(2)) then
          if(bc%bcinfo(n)%istart(2).eq.1) then
            bc%bcinfo(n)%face = 'jmin'
            bc%bcinfo(n)%istart(2) = bc%bcinfo(n)%istart(2)-2
          else
            bc%bcinfo(n)%face = 'jmax'
            bc%bcinfo(n)%iend(2) = bc%bcinfo(n)%iend(2)+2
          end if
        else if(bc%bcinfo(n)%istart(1).eq.bc%bcinfo(n)%iend(1)) then
          if(bc%bcinfo(n)%istart(1).eq.1) then
            bc%bcinfo(n)%face = 'imin'
            bc%bcinfo(n)%istart(1) = bc%bcinfo(n)%istart(1)-2
          else
            bc%bcinfo(n)%face = 'imax'
            bc%bcinfo(n)%iend(1) = bc%bcinfo(n)%iend(1)+2
          end if
        else if(bc%bcinfo(n)%istart(3).eq.bc%bcinfo(n)%iend(3)) then
          if(bc%bcinfo(n)%istart(3).eq.1) then
            bc%bcinfo(n)%face = 'kmin'
            bc%bcinfo(n)%istart(3) = bc%bcinfo(n)%istart(3)-2
          else
            bc%bcinfo(n)%face = 'kmax'
            bc%bcinfo(n)%iend(3) = bc%bcinfo(n)%iend(3)+2
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
        else if(trim(bc%bcinfo(n)%bcname).eq.'BCO112utflow') then
          bc%bcinfo(n)%bctype => bcoutflow
        else if(trim(bc%bcinfo(n)%bcname).eq.'bcoutflowsubsonic') then
          bc%bcinfo(n)%bctype => bcoutflowsubsonic
        else if(trim(bc%bcinfo(n)%bcname).eq.'bcoutflowsupersonic') then
          bc%bcinfo(n)%bctype => bcoutflowsupersonic
        else if(trim(bc%bcinfo(n)%bcname).eq.'bcfarfield') then
          bc%bcinfo(n)%bctype => bcfarfield
        else if(trim(bc%bcinfo(n)%bcname).eq.'bcdegeneratepoint') then
          bc%bcinfo(n)%bctype => bcdegeneratepoint
        else if(trim(bc%bcinfo(n)%bcname).eq.'bcsymmetryplane') then
          bc%bcinfo(n)%bctype => bcdegeneratepoint
        else
          bc%bcinfo(n)%bctype => null()
          write(*,*) 'error, check bc name'
          stop
        end if
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

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(bc)
      implicit none
      class(t_bc), intent(inout) :: bc
      integer :: n
      
      do n=1,bc%nbc
        if(associated(bc%bcinfo(n)%bctype)) nullify(bc%bcinfo(n)%bctype)
      end do
      
      deallocate(bc%ref%pv,bc%ref%dv,bc%ref%tv)
      deallocate(bc%bcinfo,bc%connectinfo)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setbc(bc,grid,variable,eos,prop)
      implicit none
      include 'mpif.h'
      class(t_bc), intent(inout) :: bc
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer :: i,j,k,n,m,l
      integer :: ii,jj,kk,num
      integer :: ier,request1,request2,request3,request4
      integer :: status(mpi_status_size)
      integer :: sendadress(21),recvadress(21)
      real(8) :: pv(bc%npv),dv(bc%ndv),tv(bc%ntv)
      real(8), dimension(:), allocatable :: sendbuf,recvbuf
      
      do n=1,bc%nbc
        call bc%bcinfo(n)%bctype(bc%ref,grid,variable,eos,prop)
      end do
      
      do n=1,bc%ncon
        num = (abs(bc%connectinfo(n)%istart(2)-bc%connectinfo(n)%iend(2))+1) &
             *(abs(bc%connectinfo(n)%istart(1)-bc%connectinfo(n)%iend(1))+1) &
             *(abs(bc%connectinfo(n)%istart(3)-bc%connectinfo(n)%iend(3))+1)*(bc%npv+bc%ntv+bc%ndv)
        allocate(sendbuf(num))
        l = 0
        do k=bc%connectinfo(n)%istart(3),bc%connectinfo(n)%iend(3)
          do j=bc%connectinfo(n)%istart(2),bc%connectinfo(n)%iend(2)
            do i=bc%connectinfo(n)%istart(1),bc%connectinfo(n)%iend(1)
              pv = variable%getpv(i,j,k)
              tv = variable%gettv(i,j,k)
              dv = variable%getdv(i,j,k)
              do m=1,bc%npv
                l = l + 1
                sendbuf(l) = pv(m)
              end do
              do m=1,bc%ntv
                l = l + 1
                sendbuf(l) = tv(m)
              end do
              do m=1,bc%ndv
                l = l + 1
                sendbuf(l) = dv(m)
              end do            
            end do
          end do
        end do
          
        if(bc%rank.ne.bc%connectinfo(n)%donor) then 
          sendadress = (/bc%connectinfo(n)%istart(1),bc%connectinfo(n)%iend(1), &
                         bc%connectinfo(n)%istart(2),bc%connectinfo(n)%iend(2), &
                         bc%connectinfo(n)%istart(3),bc%connectinfo(n)%iend(3), &
                         bc%connectinfo(n)%istart_donor(1),bc%connectinfo(n)%iend_donor(1), &
                         bc%connectinfo(n)%istart_donor(2),bc%connectinfo(n)%iend_donor(2), &
                         bc%connectinfo(n)%istart_donor(3),bc%connectinfo(n)%iend_donor(3), &  
                         bc%connectinfo(n)%transmat(1,1),bc%connectinfo(n)%transmat(1,2),bc%connectinfo(n)%transmat(1,3), &
                         bc%connectinfo(n)%transmat(2,1),bc%connectinfo(n)%transmat(2,2),bc%connectinfo(n)%transmat(2,3), &
                         bc%connectinfo(n)%transmat(3,1),bc%connectinfo(n)%transmat(3,2),bc%connectinfo(n)%transmat(3,3)/)
                     

          call mpi_isend(sendadress,21,mpi_integer,bc%connectinfo(n)%donor,bc%rank+10000,mpi_comm_world,request1,ier)
          call mpi_isend(sendbuf,num,mpi_real8,bc%connectinfo(n)%donor,bc%rank,mpi_comm_world,request2,ier)
          call mpi_wait(request1,status,ier)
          call mpi_wait(request2,status,ier)
        else  
          l = 0
          do k=bc%connectinfo(n)%istart(3),bc%connectinfo(n)%iend(3)
            do j=bc%connectinfo(n)%istart(2),bc%connectinfo(n)%iend(2)
              do i=bc%connectinfo(n)%istart(1),bc%connectinfo(n)%iend(1)
                ii = bc%connectinfo(n)%transmat(1,1)*(i-bc%connectinfo(n)%istart(1)) &
                   + bc%connectinfo(n)%transmat(1,2)*(j-bc%connectinfo(n)%istart(2)) &
                   + bc%connectinfo(n)%transmat(1,3)*(k-bc%connectinfo(n)%istart(3)) + bc%connectinfo(n)%istart_donor(1)
                jj = bc%connectinfo(n)%transmat(2,1)*(i-bc%connectinfo(n)%istart(1)) &
                   + bc%connectinfo(n)%transmat(2,2)*(j-bc%connectinfo(n)%istart(2)) &
                   + bc%connectinfo(n)%transmat(2,3)*(k-bc%connectinfo(n)%istart(3)) + bc%connectinfo(n)%istart_donor(2)
                kk = bc%connectinfo(n)%transmat(3,1)*(i-bc%connectinfo(n)%istart(1)) &
                   + bc%connectinfo(n)%transmat(3,2)*(j-bc%connectinfo(n)%istart(2)) &
                   + bc%connectinfo(n)%transmat(3,3)*(k-bc%connectinfo(n)%istart(3)) + bc%connectinfo(n)%istart_donor(3)
                do m=1,bc%npv
                  l = l + 1
                  call variable%setpv(m,ii,jj,kk,sendbuf(l))
                end do
                do m=1,bc%ntv
                  l = l + 1
                  call variable%settv(m,ii,jj,kk,sendbuf(l))
                end do
                do m=1,bc%ndv
                  l = l + 1
                  call variable%setdv(m,ii,jj,kk,sendbuf(l))
                end do
              end do
            end do
          end do
        end if
        deallocate(sendbuf)
      end do

      do n=1,bc%ncon
        num = (abs(bc%connectinfo(n)%istart(2)-bc%connectinfo(n)%iend(2))+1) &
             *(abs(bc%connectinfo(n)%istart(1)-bc%connectinfo(n)%iend(1))+1) &
             *(abs(bc%connectinfo(n)%istart(3)-bc%connectinfo(n)%iend(3))+1)*(bc%npv+bc%ntv+bc%ndv)
        allocate(recvbuf(num))
        if(bc%rank.ne.bc%connectinfo(n)%donor) then 
          call mpi_irecv(recvadress,21,mpi_integer,bc%connectinfo(n)%donor,bc%connectinfo(n)%donor+10000,mpi_comm_world,request3,ier)
          call mpi_irecv(recvbuf,num,mpi_real8,bc%connectinfo(n)%donor,bc%connectinfo(n)%donor,mpi_comm_world,request4,ier)
          call mpi_wait(request3,status,ier)
          call mpi_wait(request4,status,ier)
      
          l = 0
          do k=recvadress(5),recvadress(6)
            do j=recvadress(3),recvadress(4)
              do i=recvadress(1),recvadress(2)
                ii = recvadress(13)*(i-recvadress(1)) + recvadress(14)*(j-recvadress(3)) + recvadress(15)*(k-recvadress(5)) + recvadress(7)
                jj = recvadress(16)*(i-recvadress(1)) + recvadress(17)*(j-recvadress(3)) + recvadress(18)*(k-recvadress(5)) + recvadress(9)
                kk = recvadress(19)*(i-recvadress(1)) + recvadress(20)*(j-recvadress(3)) + recvadress(21)*(k-recvadress(5)) + recvadress(11)
                do m=1,bc%npv
                  l = l + 1
                  call variable%setpv(m,ii,jj,kk,recvbuf(l))
                end do
                do m=1,bc%ntv
                  l = l + 1
                  call variable%settv(m,ii,jj,kk,recvbuf(l))
                end do
                do m=1,bc%ndv
                  l = l + 1
                  call variable%setdv(m,ii,jj,kk,recvbuf(l))
                end do
              end do
            end do
          end do
        end if
        deallocate(recvbuf)
      end do
      
      do i=2,grid%getimax()
        call bc%bccorner(grid,variable,eos,prop,i,2,2,0,-1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,1,2,0,-1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,2,1,0,-1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,1,1,0,-1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()  ,2,0,1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()+1,2,0,1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()  ,1,0,1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()+1,1,0,1,-1)
        call bc%bccorner(grid,variable,eos,prop,i,2,grid%getkmax()  ,0,-1,1)
        call bc%bccorner(grid,variable,eos,prop,i,1,grid%getkmax()  ,0,-1,1)
        call bc%bccorner(grid,variable,eos,prop,i,2,grid%getkmax()+1,0,-1,1)
        call bc%bccorner(grid,variable,eos,prop,i,1,grid%getkmax()+1,0,-1,1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()  ,grid%getkmax()  ,0,1,1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()+1,grid%getkmax()  ,0,1,1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()  ,grid%getkmax()+1,0,1,1)
        call bc%bccorner(grid,variable,eos,prop,i,grid%getjmax()+1,grid%getkmax()+1,0,1,1)
      end do
      do j=2,grid%getjmax()
        call bc%bccorner(grid,variable,eos,prop,2,j,2,-1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,1,j,2,-1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,2,j,1,-1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,1,j,1,-1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,j,2,1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,j,2,1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,j,1,1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,j,1,1,0,-1)
        call bc%bccorner(grid,variable,eos,prop,2,j,grid%getkmax()  ,-1,0,1)
        call bc%bccorner(grid,variable,eos,prop,1,j,grid%getkmax()  ,-1,0,1)
        call bc%bccorner(grid,variable,eos,prop,2,j,grid%getkmax()+1,-1,0,1)
        call bc%bccorner(grid,variable,eos,prop,1,j,grid%getkmax()+1,-1,0,1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,j,grid%getkmax()  ,1,0,1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,j,grid%getkmax()  ,1,0,1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,j,grid%getkmax()+1,1,0,1)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,j,grid%getkmax()+1,1,0,1)
      end do
      do k=2,grid%getkmax()
        call bc%bccorner(grid,variable,eos,prop,2,2,k,-1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,1,2,k,-1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,2,1,k,-1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,1,1,k,-1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,2,k,1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,2,k,1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,1,k,1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,1,k,1,-1,0)
        call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()  ,k,-1,1,0)
        call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()  ,k,-1,1,0)
        call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()+1,k,-1,1,0)
        call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()+1,k,-1,1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()  ,k,1,1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()  ,k,1,1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()+1,k,1,1,0)
        call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()+1,k,1,1,0)
      end do
      
      call bc%bccorner(grid,variable,eos,prop,2,2,2,-1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,2,2,-1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,2,1,2,-1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,2,2,1,-1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,1,2,-1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,2,1,1,-1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,2,1,-1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,1,1,-1,-1,-1)
      
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,2,2,1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,2,2,1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,1,2,1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,2,1,1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,1,2,1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,1,1,1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,2,1,1,-1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,1,1,1,-1,-1)
      
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()  ,2,-1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()  ,2,-1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()+1,2,-1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()  ,1,-1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()+1,2,-1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()+1,1,-1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()  ,1,-1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()+1,1,-1,1,-1)
      
      call bc%bccorner(grid,variable,eos,prop,2,2,grid%getkmax()  ,-1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,1,2,grid%getkmax()  ,-1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,2,1,grid%getkmax()  ,-1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,2,2,grid%getkmax()+1,-1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,1,1,grid%getkmax()  ,-1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,2,1,grid%getkmax()+1,-1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,1,2,grid%getkmax()+1,-1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,1,1,grid%getkmax()+1,-1,-1,1)
      
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()  ,2,1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()  ,2,1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()+1,2,1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()  ,1,1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()+1,2,1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()+1,1,1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()  ,1,1,1,-1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()+1,1,1,1,-1)
      
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()  ,grid%getkmax()  ,-1,1,1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()  ,grid%getkmax()  ,-1,1,1)
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()+1,grid%getkmax()  ,-1,1,1)
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()  ,grid%getkmax()+1,-1,1,1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()+1,grid%getkmax()  ,-1,1,1)
      call bc%bccorner(grid,variable,eos,prop,2,grid%getjmax()+1,grid%getkmax()+1,-1,1,1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()  ,grid%getkmax()+1,-1,1,1)
      call bc%bccorner(grid,variable,eos,prop,1,grid%getjmax()+1,grid%getkmax()+1,-1,1,1)
      
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,2,grid%getkmax()  ,1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,2,grid%getkmax()  ,1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,1,grid%getkmax()  ,1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,2,grid%getkmax()+1,1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,1,grid%getkmax()  ,1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,1,grid%getkmax()+1,1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,2,grid%getkmax()+1,1,-1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,1,grid%getkmax()+1,1,-1,1)
      
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()  ,grid%getkmax()  ,1,1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()  ,grid%getkmax()  ,1,1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()+1,grid%getkmax()  ,1,1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()  ,grid%getkmax()+1,1,1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()+1,grid%getkmax()  ,1,1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()  ,grid%getjmax()+1,grid%getkmax()+1,1,1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()  ,grid%getkmax()+1,1,1,1)
      call bc%bccorner(grid,variable,eos,prop,grid%getimax()+1,grid%getjmax()+1,grid%getkmax()+1,1,1,1)
      
    end subroutine setbc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine bccorner(bc,grid,variable,eos,prop,i,j,k,ii,jj,kk)
      implicit none
      class(t_bc), intent(in) :: bc
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer, intent(in) :: i,j,k,ii,jj,kk
      integer :: n
      real(8) :: pv(bc%npv),pvx(bc%npv),pvy(bc%npv),pvz(bc%npv)
      real(8) :: dv(bc%ndv),tv(bc%ntv)
      
      pv = variable%getpv(i,j,k)
      pvx = variable%getpv(i+ii,j,k)
      pvy = variable%getpv(i,j+jj,k)
      pvz = variable%getpv(i,j,k+kk)
    
      if((-pv(2).eq.pvx(2)).and.(-pv(3).eq.pvx(3)).and.(-pv(4).eq.pvx(4)) then
        if((-pv(2).eq.pvy(2)).and.(-pv(3).eq.pvy(3)).and.(-pv(4).eq.pvy(4))) then
          if((-pv(2).eq.pvz(2)).and.(-pv(3).eq.pvz(3)).and.(-pv(4).eq.pvz(4))) then
!----------- wall-wall-wall(x,y,z)
            pv = variable%getpv(i,j,k)
            pv(2) = -pv(2)
            pv(3) = -pv(3)
            pv(4) = -pv(4)
            if(bc%iturb.eq.0) then
              pv(8) = -pv(8)
            else if(bc%iturb.eq.-1) then
              pv(8) = -pv(8)
              pv(9) = -pv(9)
            end if
          else
!----------- wall-wall(x,y)
            if(kk.eq.0) then
              pv = variable%getpv(i,j,k)
            else
            
            end if
          end if
        else 
          if((-pv(2).eq.pvz(2)).and.(-pv(3).eq.pvz(3)).and.(-pv(4).eq.pvz(4))) then
!----------- wall-wall(x,z)
            if(jj.eq.0) then
              pv = variable%getpv(i,j,k)
            else
            
            end if
          else 
!----------- wall(x)
            if((jj.eq.0).or.(kk.eq.0)) then
              pv = variable%getpv(i,j+jj,k+kk)
            else
            
            end if
          end if
        end if
      else
        if((-pv(2).eq.pvy(2)).and.(-pv(3).eq.pvy(3)).and.(-pv(4).eq.pvy(4))) then
          if((-pv(2).eq.pvz(2)).and.(-pv(3).eq.pvz(3)).and.(-pv(4).eq.pvz(4))) then
!----------- wall-wall(y,z)
            if(ii.eq.0) then
              pv = variable%getpv(i,j,k) 
            else
            
            end if
          else
!----------- wall(y)
            pv = variable%getpv(i+ii,j,k+kk)
          end if
        else 
          if((-pv(2).eq.pvz(2)).and.(-pv(3).eq.pvz(3)).and.(-pv(4).eq.pvz(4))) then
!----------- wall(z)
            pv = variable%getpv(i,j,k)
          else
!----------- no wall
            pv = 1.d0/3.d0*(pvx+pvy+pvz)
          end if
        end if
      end if
      
      do n=1,bc%npv
        call variable%setpv(n,i+ii,j+jj,k+kk,pv(n))
      end do
      
      call eos%deteos(pv(1)+bc%ref%pv(1),pv(5),pv(6),pv(7),dv)
      
      do n=1,bc%ndv
        call variable%setdv(n,i+ii,j+jj,k+kk,dv(n))
      end do
      
      if(bc%ntv.ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
      
      do n=1,bc%ntv
        call variable%settv(n,i+ii,j+jj,k+kk,tv(n))
      end do
      

    end subroutine bccorner
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
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%getcx(bcinfo%iend(1),j,k)
              ii = bcinfo%iend(1)-i + bcinfo%iend(1)+1
              pv = variable%getpv(ii,j,k)
              dv = variable%getdv(ii,j,k)
              tv = variable%gettv(ii,j,k)
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
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%getcx(bcinfo%istart(1)-1,j,k)
              ii = bcinfo%istart(1)-i + bcinfo%istart(1)-1
              pv = variable%getpv(ii,j,k)
              dv = variable%getdv(ii,j,k)
              tv = variable%gettv(ii,j,k)
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
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%getex(i,bcinfo%iend(2),k)
              jj = bcinfo%iend(2)-j + bcinfo%iend(2)+1
              pv = variable%getpv(i,jj,k)
              dv = variable%getdv(i,jj,k)
              tv = variable%gettv(i,jj,k)
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
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%getex(i,bcinfo%istart(2)-1,k)
              jj = bcinfo%istart(2)-j + bcinfo%istart(2)-1
              pv = variable%getpv(i,jj,k)
              dv = variable%getdv(i,jj,k)
              tv = variable%gettv(i,jj,k)
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
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%gettx(i,j,bcinfo%iend(3))
              kk = bcinfo%iend(3)-k + bcinfo%iend(3)+1
              pv = variable%getpv(i,j,kk)
              dv = variable%getdv(i,j,kk)
              tv = variable%gettv(i,j,kk)
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
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%gettx(i,j,bcinfo%istart(3)-1)
              kk = bcinfo%istart(3)-k + bcinfo%istart(3)-1
              pv = variable%getpv(i,j,kk)
              dv = variable%getdv(i,j,kk)
              tv = variable%gettv(i,j,kk)
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
      end select
    
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
      integer :: i,j,k,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              dv = variable%getdv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
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
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              dv = variable%getdv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
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
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              dv = variable%getdv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
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
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              dv = variable%getdv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
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
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              dv = variable%getdv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
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
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              dv = variable%getdv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
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
      end select
      
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
      integer :: i,j,k,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())

      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              dv = variable%getdv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
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
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              dv = variable%getdv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
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
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              dv = variable%getdv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
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
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              dv = variable%getdv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
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
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              dv = variable%getdv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
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
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              dv = variable%getdv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
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
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      end select
      
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
      integer :: i,j,k,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv()),grd(grid%getngrd())
      real(8) :: x(3),x1(3),x2(3),x3(3),var,xcc,ycc,zcc

      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              dv = variable%getdv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
              do m=1,variable%getnpv()
                select case(m)
                case(2,3,4,8)
                  call variable%setpv(m,i,j,k,-pv(m))
                case(9)
                  x = grid%getx(bcinfo%iend(1)+1,j,k)
                  x1 = grid%getx(bcinfo%iend(1)+1,j+1,k)
                  x2 = grid%getx(bcinfo%iend(1)+1,j,k+1)
                  x3 = grid%getx(bcinfo%iend(1)+1,j+1,k+1)
                  grd = grid%getgrd(bcinfo%iend(1)+1,j,k)
                  xcc = 0.25d0*(x(1)+x1(1)+x2(1)+x3(1))
                  ycc = 0.25d0*(x(2)+x1(2)+x2(2)+x3(2))
                  zcc = 0.25d0*(x(3)+x1(3)+x2(3)+x3(3))
                  var = 1600.d0*tv(1)/dv(1)/((grd(2)-xcc)**2+(grd(3)-ycc)**2+(grd(4)-zcc)**2)-pv(m)
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
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              dv = variable%getdv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
              do m=1,variable%getnpv()
                select case(m)
                case(2,3,4,8)
                  call variable%setpv(m,i,j,k,-pv(m))
                case(9)
                  x = grid%getx(bcinfo%istart(1),j,k)
                  x1 = grid%getx(bcinfo%istart(1),j+1,k)
                  x2 = grid%getx(bcinfo%istart(1),j,k+1)
                  x3 = grid%getx(bcinfo%istart(1),j+1,k+1)
                  grd = grid%getgrd(bcinfo%istart(1)-1,j,k)
                  xcc = 0.25d0*(x(1)+x1(1)+x2(1)+x3(1))
                  ycc = 0.25d0*(x(2)+x1(2)+x2(2)+x3(2))
                  zcc = 0.25d0*(x(3)+x1(3)+x2(3)+x3(3))
                  var = 1600.d0*tv(1)/dv(1)/((grd(2)-xcc)**2+(grd(3)-ycc)**2+(grd(4)-zcc)**2)-pv(m)
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
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              dv = variable%getdv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
              do m=1,variable%getnpv()
                select case(m)
                case(2,3,4,8)
                  call variable%setpv(m,i,j,k,-pv(m))
                case(9)
                  x = grid%getx(i,bcinfo%iend(2)+1,k)
                  x1 = grid%getx(i+1,bcinfo%iend(2)+1,k)
                  x2 = grid%getx(i,bcinfo%iend(2)+1,k+1)
                  x3 = grid%getx(i+1,bcinfo%iend(2)+1,k+1)
                  grd = grid%getgrd(i,bcinfo%iend(2)+1,k)
                  xcc = 0.25d0*(x(1)+x1(1)+x2(1)+x3(1))
                  ycc = 0.25d0*(x(2)+x1(2)+x2(2)+x3(2))
                  zcc = 0.25d0*(x(3)+x1(3)+x2(3)+x3(3))
                  var = 1600.d0*tv(1)/dv(1)/((grd(2)-xcc)**2+(grd(3)-ycc)**2+(grd(4)-zcc)**2)-pv(m)
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
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              dv = variable%getdv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
              do m=1,variable%getnpv()
                select case(m)
                case(2,3,4,8)
                  call variable%setpv(m,i,j,k,-pv(m))
                case(9)
                  x = grid%getx(i,bcinfo%istart(2),k)
                  x1 = grid%getx(i+1,bcinfo%istart(2),k)
                  x2 = grid%getx(i,bcinfo%istart(2),k+1)
                  x3 = grid%getx(i+1,bcinfo%istart(2),k+1)
                  grd = grid%getgrd(i,bcinfo%istart(2)-1,k)
                  xcc = 0.25d0*(x(1)+x1(1)+x2(1)+x3(1))
                  ycc = 0.25d0*(x(2)+x1(2)+x2(2)+x3(2))
                  zcc = 0.25d0*(x(3)+x1(3)+x2(3)+x3(3))
                  var = 1600.d0*tv(1)/dv(1)/((grd(2)-xcc)**2+(grd(3)-ycc)**2+(grd(4)-zcc)**2)-pv(m)
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
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              dv = variable%getdv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
              do m=1,variable%getnpv()
                select case(m)
                case(2,3,4,8)
                  call variable%setpv(m,i,j,k,-pv(m))
                case(9)
                  x = grid%getx(i,j,bcinfo%iend(3)+1)
                  x1 = grid%getx(i+1,j,bcinfo%iend(3)+1)
                  x2 = grid%getx(i,j+1,bcinfo%iend(3)+1)
                  x3 = grid%getx(i+1,j+1,bcinfo%iend(3)+1)
                  grd = grid%getgrd(i,j,bcinfo%iend(3)+1)
                  xcc = 0.25d0*(x(1)+x1(1)+x2(1)+x3(1))
                  ycc = 0.25d0*(x(2)+x1(2)+x2(2)+x3(2))
                  zcc = 0.25d0*(x(3)+x1(3)+x2(3)+x3(3))
                  var = 1600.d0*tv(1)/dv(1)/((grd(2)-xcc)**2+(grd(3)-ycc)**2+(grd(4)-zcc)**2)-pv(m)
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
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              dv = variable%getdv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
              do m=1,variable%getnpv()
                select case(m)
                case(2,3,4,8)
                  call variable%setpv(m,i,j,k,-pv(m))
                case(9)
                  x = grid%getx(i,j,bcinfo%istart(3))
                  x1 = grid%getx(i+1,j,bcinfo%istart(3))
                  x2 = grid%getx(i,j+1,bcinfo%istart(3))
                  x3 = grid%getx(i+1,j+1,bcinfo%istart(3))
                  grd = grid%getgrd(i,j,bcinfo%istart(3)-1)
                  xcc = 0.25d0*(x(1)+x1(1)+x2(1)+x3(1))
                  ycc = 0.25d0*(x(2)+x1(2)+x2(2)+x3(2))
                  zcc = 0.25d0*(x(3)+x1(3)+x2(3)+x3(3))
                  var = 1600.d0*tv(1)/dv(1)/((grd(2)-xcc)**2+(grd(3)-ycc)**2+(grd(4)-zcc)**2)-pv(m)
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
      end select
      
      
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
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      end select
      
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
      integer :: i,j,k,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
              call variable%setpv(1,i,j,k,pv(1))
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
              call variable%setpv(1,i,j,k,pv(1))
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
              call variable%setpv(1,i,j,k,pv(1))
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
              call variable%setpv(1,i,j,k,pv(1))
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
              call variable%setpv(1,i,j,k,pv(1))
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
              call variable%setpv(1,i,j,k,pv(1))
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,dv(m))
              end do
              if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
              do m=1,variable%getntv()
                call variable%settv(m,i,j,k,tv(m))
              end do
            end do
          end do
        end do
      end select
      
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
      real(8) :: tv(variable%getntv())
      integer :: i,j,k,m
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,ref%dv(m))
              end do
              do m=1,variable%getntv()
                select case(m)
                case(3)
                  call variable%settv(m,i,j,k,tv(m))
                case default
                  call variable%settv(m,i,j,k,ref%tv(m))
                end select
              end do
            end do
          end do
        end do
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,ref%dv(m))
              end do
              do m=1,variable%getntv()
                select case(m)
                case(3)
                  call variable%settv(m,i,j,k,tv(m))
                case default
                  call variable%settv(m,i,j,k,ref%tv(m))
                end select
              end do
            end do
          end do
        end do
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,ref%dv(m))
              end do
              do m=1,variable%getntv()
                select case(m)
                case(3)
                  call variable%settv(m,i,j,k,tv(m))
                case default
                  call variable%settv(m,i,j,k,ref%tv(m))
                end select
              end do
            end do
          end do
        end do
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k) 
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,ref%dv(m))
              end do
              do m=1,variable%getntv()
                select case(m)
                case(3)
                  call variable%settv(m,i,j,k,tv(m))
                case default
                  call variable%settv(m,i,j,k,ref%tv(m))
                end select
              end do
            end do
          end do
        end do
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,ref%dv(m))
              end do
              do m=1,variable%getntv()
                select case(m)
                case(3)
                  call variable%settv(m,i,j,k,tv(m))
                case default
                  call variable%settv(m,i,j,k,ref%tv(m))
                end select
              end do
            end do
          end do
        end do
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1) 
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,ref%pv(m))
              end do
              do m=1,variable%getndv()
                call variable%setdv(m,i,j,k,ref%dv(m))
              end do
              do m=1,variable%getntv()
                select case(m)
                case(3)
                  call variable%settv(m,i,j,k,tv(m))
                case default
                  call variable%settv(m,i,j,k,ref%tv(m))
                end select
              end do
            end do
          end do
        end do
      end select

      
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
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
            
            end do
          end do
        end do
      end select
      
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
      integer :: i,j,k,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
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
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
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
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
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
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
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
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
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
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
              call variable%setpv(1,i,j,k,0.d0)
              do m=2,variable%getnpv()
                call variable%setpv(m,i,j,k,pv(m))
              end do
              call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
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
      end select
      
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
      integer :: i,j,k,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              dv = variable%getdv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
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
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              dv = variable%getdv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
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
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              dv = variable%getdv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
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
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              dv = variable%getdv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
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
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              dv = variable%getdv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
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
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              dv = variable%getdv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
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
      end select
    
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
      integer :: i,j,k,m
      real(8) :: nx(3),vel
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%getcx(bcinfo%iend(1),j,k)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              dv = variable%getdv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
              vel = (pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              if(vel/dsqrt(dv(6)).ge.0.d0) then !out
                if(vel/dsqrt(dv(6)).ge.1.d0) then !super                
                  do m=1,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              else ! in
                if(vel/dsqrt(dv(6)).le.-1.d0) then !super
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,ref%dv(m))
                  end do
                  do m=1,variable%getntv()
                    select case(m)
                    case(3)
                      call variable%settv(m,i,j,k,tv(m))
                    case default
                      call variable%settv(m,i,j,k,ref%tv(m))
                    end select
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,pv(1))
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              end if
            end do
          end do
        end do
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%getcx(bcinfo%istart(1)-1,j,k)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              dv = variable%getdv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
              vel = (pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              if(vel/dsqrt(dv(6)).ge.0.d0) then !out
                if(vel/dsqrt(dv(6)).ge.1.d0) then !super                
                  do m=1,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              else ! in
                if(vel/dsqrt(dv(6)).le.-1.d0) then !super
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,ref%dv(m))
                  end do
                  do m=1,variable%getntv()
                    select case(m)
                    case(3)
                      call variable%settv(m,i,j,k,tv(m))
                    case default
                      call variable%settv(m,i,j,k,ref%tv(m))
                    end select
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,pv(1))
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              end if
            end do
          end do
        end do
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%getex(i,bcinfo%iend(2),k)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              dv = variable%getdv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
              vel = (pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              if(vel/dsqrt(dv(6)).ge.0.d0) then !out
                if(vel/dsqrt(dv(6)).ge.1.d0) then !super                
                  do m=1,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              else ! in
                if(vel/dsqrt(dv(6)).le.-1.d0) then !super
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,ref%dv(m))
                  end do
                  do m=1,variable%getntv()
                    select case(m)
                    case(3)
                      call variable%settv(m,i,j,k,tv(m))
                    case default
                      call variable%settv(m,i,j,k,ref%tv(m))
                    end select
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,pv(1))
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              end if
            end do
          end do
        end do
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%getex(i,bcinfo%istart(2)-1,k)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              dv = variable%getdv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
              vel = (pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              if(vel/dsqrt(dv(6)).ge.0.d0) then !out
                if(vel/dsqrt(dv(6)).ge.1.d0) then !super                
                  do m=1,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              else ! in
                if(vel/dsqrt(dv(6)).le.-1.d0) then !super
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,ref%dv(m))
                  end do
                  do m=1,variable%getntv()
                    select case(m)
                    case(3)
                      call variable%settv(m,i,j,k,tv(m))
                    case default
                      call variable%settv(m,i,j,k,ref%tv(m))
                    end select
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,pv(1))
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              end if
            end do
          end do
        end do
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%getex(i,j,bcinfo%iend(3))
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              dv = variable%getdv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
              vel = (pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              if(vel/dsqrt(dv(6)).ge.0.d0) then !out
                if(vel/dsqrt(dv(6)).ge.1.d0) then !super                
                  do m=1,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              else ! in
                if(vel/dsqrt(dv(6)).le.-1.d0) then !super
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,ref%dv(m))
                  end do
                  do m=1,variable%getntv()
                    select case(m)
                    case(3)
                      call variable%settv(m,i,j,k,tv(m))
                    case default
                      call variable%settv(m,i,j,k,ref%tv(m))
                    end select
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,pv(1))
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              end if
            end do
          end do
        end do
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%getex(i,j,bcinfo%istart(3)-1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              dv = variable%getdv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
              vel = (pv(2)*nx(1)+pv(3)*nx(2)+pv(4)*nx(3))/dsqrt(nx(1)**2+nx(2)**2+nx(3)**2)
              if(vel/dsqrt(dv(6)).ge.0.d0) then !out
                if(vel/dsqrt(dv(6)).ge.1.d0) then !super                
                  do m=1,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,pv(m))
                  end do
                  call eos%deteos(ref%pv(1),pv(5),pv(6),pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              else ! in
                if(vel/dsqrt(dv(6)).le.-1.d0) then !super
                  call variable%setpv(1,i,j,k,0.d0)
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,ref%dv(m))
                  end do
                  do m=1,variable%getntv()
                    select case(m)
                    case(3)
                      call variable%settv(m,i,j,k,tv(m))
                    case default
                      call variable%settv(m,i,j,k,ref%tv(m))
                    end select
                  end do
                else !sub
                  call variable%setpv(1,i,j,k,pv(1))
                  do m=2,variable%getnpv()
                    call variable%setpv(m,i,j,k,ref%pv(m))
                  end do
                  call eos%deteos(pv(1)+ref%pv(1),ref%pv(5),ref%pv(6),ref%pv(7),dv)
                  do m=1,variable%getndv()
                    call variable%setdv(m,i,j,k,dv(m))
                  end do
                  if(variable%getntv().ne.0) call prop%detprop(dv(3),dv(4),dv(5),ref%pv(5),ref%pv(6),ref%pv(7),tv(1:2))
                  do m=1,variable%getntv()
                    call variable%settv(m,i,j,k,tv(m))
                  end do
                end if
              end if
            end do
          end do
        end do
      end select
    
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
      integer :: i,j,k,m
      real(8) :: pv(variable%getnpv()),dv(variable%getndv()),tv(variable%getntv())
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%iend(1)+1,j,k)
              dv = variable%getdv(bcinfo%iend(1)+1,j,k)
              tv = variable%gettv(bcinfo%iend(1)+1,j,k)
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
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(bcinfo%istart(1)-1,j,k)
              dv = variable%getdv(bcinfo%istart(1)-1,j,k)
              tv = variable%gettv(bcinfo%istart(1)-1,j,k)
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
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%iend(2)+1,k)
              dv = variable%getdv(i,bcinfo%iend(2)+1,k)
              tv = variable%gettv(i,bcinfo%iend(2)+1,k)
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
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,bcinfo%istart(2)-1,k)
              dv = variable%getdv(i,bcinfo%istart(2)-1,k)
              tv = variable%gettv(i,bcinfo%istart(2)-1,k)
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
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%iend(3)+1)
              dv = variable%getdv(i,j,bcinfo%iend(3)+1)
              tv = variable%gettv(i,j,bcinfo%iend(3)+1)
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
            do i=bcinfo%istart(1),bcinfo%iend(1)
              pv = variable%getpv(i,j,bcinfo%istart(3)-1)
              dv = variable%getdv(i,j,bcinfo%istart(3)-1)
              tv = variable%gettv(i,j,bcinfo%istart(3)-1)
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
      end select
    
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
      real(8) :: var,nx(3)
      
      select case(bcinfo%face)
      case('imin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%getcx(bcinfo%iend(1),j,k)
              ii = bcinfo%iend(1)-i + bcinfo%iend(1)+1
              pv = variable%getpv(ii,j,k)
              dv = variable%getdv(ii,j,k)
              tv = variable%gettv(ii,j,k)
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
      case('imax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%getcx(bcinfo%istart(1)-1,j,k)
              ii = bcinfo%istart(1)-i + bcinfo%istart(1)-1
              pv = variable%getpv(ii,j,k)
              dv = variable%getdv(ii,j,k)
              tv = variable%gettv(ii,j,k)
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
      case('jmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%getex(i,bcinfo%iend(2),k)
              jj = bcinfo%iend(2)-j + bcinfo%iend(2)+1
              pv = variable%getpv(i,jj,k)
              dv = variable%getdv(i,jj,k)
              tv = variable%gettv(i,jj,k)
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
      case('jmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%getex(i,bcinfo%istart(2)-1,k)
              jj = bcinfo%istart(2)-j + bcinfo%istart(2)-1
              pv = variable%getpv(i,jj,k)
              dv = variable%getdv(i,jj,k)
              tv = variable%gettv(i,jj,k)
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
      case('kmin')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = grid%gettx(i,j,bcinfo%iend(3))
              kk = bcinfo%iend(3)-k + bcinfo%iend(3)+1
              pv = variable%getpv(i,j,kk)
              dv = variable%getdv(i,j,kk)
              tv = variable%gettv(i,j,kk)
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
      case('kmax')
        do k=bcinfo%istart(3),bcinfo%iend(3)
          do j=bcinfo%istart(2),bcinfo%iend(2)
            do i=bcinfo%istart(1),bcinfo%iend(1)
              nx = - grid%gettx(i,j,bcinfo%istart(3)-1)
              kk = bcinfo%istart(3)-k + bcinfo%istart(3)-1
              pv = variable%getpv(i,j,kk)
              dv = variable%getdv(i,j,kk)
              tv = variable%gettv(i,j,kk)
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
      end select
    end subroutine bcsymmetryplane
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module bc_module
