module grid_module
  use mpi
  use config_module
  implicit none
#include <cgnslib_f.h>
#include <cgnstypes_f.h>
  private
  public :: t_grid,t_bcinfo,t_connectinfo,t_mpitemp

  type t_bcinfo
    integer :: istart(3),iend(3)
    character(32) :: bcname,famname
  end type t_bcinfo
  
  type t_connectinfo
    integer :: donor,transmat(3,3)
    integer :: istart(3),iend(3),istart_donor(3),iend_donor(3)
  end type t_connectinfo
  
  type t_mpitemp
    integer :: num
    integer :: sendadress(21),recvadress(21)
    real(8), dimension(:), allocatable :: sendbuf,recvbuf
  end type t_mpitemp

  type t_zoneinfo
    integer :: imax,jmax,kmax
  end type t_zoneinfo

  type t_grid
    private
    character(32) :: zonename
    integer :: rank,size
    integer :: imax,jmax,kmax,nbc,ncon,ngrd
    type(t_bcinfo), dimension(:), allocatable :: bcinfo
    type(t_connectinfo), dimension(:), allocatable ::connectinfo
    type(t_zoneinfo), dimension(:), allocatable :: zoneinfo
    real(8), dimension(:,:,:,:), allocatable :: x
    real(8), dimension(:,:,:,:), allocatable :: cx,ex,tx
    real(8), dimension(:,:,:,:), allocatable :: grd ! vol,xcen,ycen,zcen,ydns
    contains
      procedure :: construct
      procedure :: destruct
      procedure, private :: calydns
      procedure, private :: calnormal
      procedure, private :: calvolume
      procedure :: getzonename
      procedure :: getimax
      procedure :: getjmax
      procedure :: getkmax
      procedure :: getimax_zone
      procedure :: getjmax_zone
      procedure :: getkmax_zone
      procedure :: getnbc
      procedure :: getncon
      procedure :: getx
      procedure :: getcx
      procedure :: getex
      procedure :: gettx
      procedure :: getngrd
      procedure :: getgrd
      procedure :: getbcname
      procedure :: getfamname
      procedure :: getbcistart
      procedure :: getbciend
      procedure :: getconnectdonor
      procedure :: getconnecttransmat
      procedure :: getconnectistart
      procedure :: getconnectiend
      procedure :: getconnectistart_donor
      procedure :: getconnectiend_donor
  end type t_grid

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(grid,config)
      implicit none
      class(t_grid), intent(out) :: grid
      type(t_config), intent(in) :: config
      character(32) :: nodename,connectname,donorname
      character(32), dimension(:), allocatable :: famname,fambcname,zonename
      integer :: i,j,n,m
      integer :: ifile,ier
      integer :: nzone,nfam
      cgsize_t :: nmax(3),nstart(3)
      cgsize_t :: isize(9),ipnts(6),ipntsdonor(6)
      integer :: transvec(3) 
      integer :: nfambc,ngeo,ibctype
      integer :: normallist
      integer, parameter :: dim = 3
      real(8), dimension(:,:,:), allocatable :: xx,yy,zz
      
      grid%rank = config%getrank()
      grid%size = config%getsize()
      do n = 0,grid%size-1
        if(n.eq.grid%rank) then

          call cg_open_f('./'//trim(config%getname())//'.cgns',cg_mode_read,ifile,ier)
          if(ier.ne.cg_ok) call cg_error_exit_f
          
          call cg_nzones_f(ifile,1,nzone,ier)
          if(nzone.ne.config%getsize()) write(*,*) 'invalid size of domain! check cpu number'
          
          allocate(zonename(nzone),grid%zoneinfo(0:nzone-1))

          do j= 0,nzone-1
            call cg_zone_read_f(ifile,1,j+1,zonename(j+1),isize,ier)
            grid%zoneinfo(j)%imax = isize(1)
            grid%zoneinfo(j)%jmax = isize(2)
            grid%zoneinfo(j)%kmax = isize(3)
            if(j.eq.n) then
              grid%zonename = zonename(j+1)
              grid%imax = isize(1)
              grid%jmax = isize(2)
              grid%kmax = isize(3)
            end if
          end do
          
          allocate(grid%x(dim,2:grid%imax+1,2:grid%jmax+1,2:grid%kmax+1))
          allocate(xx(2:grid%imax+1,2:grid%jmax+1,2:grid%kmax+1),yy(2:grid%imax+1,2:grid%jmax+1,2:grid%kmax+1),zz(2:grid%imax+1,2:grid%jmax+1,2:grid%kmax+1))
          
          nmax(1) = grid%imax; nmax(2) = grid%jmax; nmax(3) = grid%kmax
          nstart  = 1
          
          call cg_coord_read_f(ifile,1,n+1,'CoordinateX',realdouble,nstart,nmax,xx,ier)
          call cg_coord_read_f(ifile,1,n+1,'CoordinateY',realdouble,nstart,nmax,yy,ier)
          call cg_coord_read_f(ifile,1,n+1,'CoordinateZ',realdouble,nstart,nmax,zz,ier)
          
          grid%x(1,:,:,:) = xx*config%getscale()
          grid%x(2,:,:,:) = yy*config%getscale()
          grid%x(3,:,:,:) = zz*config%getscale()
          
          call cg_nfamilies_f(ifile,1,nfam,ier)
          
          allocate(famname(nfam),fambcname(nfam))
          
          do i=1,nfam
            call cg_family_read_f(ifile,1,i,famname(i),nfambc,ngeo,ier)
            call cg_fambc_read_f(ifile,1,i,1,nodename,ibctype,ier)
            if(trim(nodename)=='FamBC') then
              fambcname(i) = bctypename(ibctype)
            else
              fambcname(i) = '....'
            end if 
          end do     
          
          call cg_nbocos_f(ifile,1,n+1,grid%nbc,ier)
          allocate(grid%bcinfo(grid%nbc))
          
          do i=1,grid%nbc
            call cg_goto_f(ifile,1,ier,'Zone_t',n+1,'ZoneBC_t',1,'BC_t',i,'end')
            call cg_famname_read_f(nodename,ier)
            call cg_boco_read_f(ifile,1,n+1,i,ipnts,normallist,ier)
            
            do j=1,dim
              if(ipnts(j).eq.ipnts(dim+j)) then
                if(ipnts(j).eq.1) then
                  select case(j)
                  case(1)
                    grid%bcinfo(i)%istart(1) = ipnts(1)
                    grid%bcinfo(i)%istart(2) = ipnts(2)+1
                    grid%bcinfo(i)%istart(3) = ipnts(3)+1
                    grid%bcinfo(i)%iend(1)   = ipnts(4)
                    grid%bcinfo(i)%iend(2)   = ipnts(5)
                    grid%bcinfo(i)%iend(3)   = ipnts(6)
                  case(2)
                    grid%bcinfo(i)%istart(1) = ipnts(1)+1
                    grid%bcinfo(i)%istart(2) = ipnts(2)
                    grid%bcinfo(i)%istart(3) = ipnts(3)+1
                    grid%bcinfo(i)%iend(1)   = ipnts(4)
                    grid%bcinfo(i)%iend(2)   = ipnts(5)
                    grid%bcinfo(i)%iend(3)   = ipnts(6)
                  case(3)
                    grid%bcinfo(i)%istart(1) = ipnts(1)+1
                    grid%bcinfo(i)%istart(2) = ipnts(2)+1
                    grid%bcinfo(i)%istart(3) = ipnts(3)
                    grid%bcinfo(i)%iend(1)   = ipnts(4)
                    grid%bcinfo(i)%iend(2)   = ipnts(5)
                    grid%bcinfo(i)%iend(3)   = ipnts(6)
                  end select
                else
                  select case(j)
                  case(1)
                    grid%bcinfo(i)%istart(1) = ipnts(1)+1
                    grid%bcinfo(i)%istart(2) = ipnts(2)+1
                    grid%bcinfo(i)%istart(3) = ipnts(3)+1
                    grid%bcinfo(i)%iend(1)   = ipnts(4)+1
                    grid%bcinfo(i)%iend(2)   = ipnts(5)
                    grid%bcinfo(i)%iend(3)   = ipnts(6)
                  case(2)
                    grid%bcinfo(i)%istart(1) = ipnts(1)+1
                    grid%bcinfo(i)%istart(2) = ipnts(2)+1
                    grid%bcinfo(i)%istart(3) = ipnts(3)+1
                    grid%bcinfo(i)%iend(1)   = ipnts(4)
                    grid%bcinfo(i)%iend(2)   = ipnts(5)+1
                    grid%bcinfo(i)%iend(3)   = ipnts(6)
                  case(3)
                    grid%bcinfo(i)%istart(1) = ipnts(1)+1
                    grid%bcinfo(i)%istart(2) = ipnts(2)+1
                    grid%bcinfo(i)%istart(3) = ipnts(3)+1
                    grid%bcinfo(i)%iend(1)   = ipnts(4)
                    grid%bcinfo(i)%iend(2)   = ipnts(5)
                    grid%bcinfo(i)%iend(3)   = ipnts(6)+1
                  end select
                end if
              end if
            end do
         
            do j=1,nfam
              if(trim(nodename).eq.trim(famname(j))) then
                grid%bcinfo(i)%bcname = fambcname(j)
                if(trim(grid%bcinfo(i)%bcname).eq.'UserDefined') then
                  grid%bcinfo(i)%bcname = famname(j)
                end if
                grid%bcinfo(i)%famname = famname(j)
              end if
            end do
          end do

          call cg_n1to1_f(ifile,1,n+1,grid%ncon,ier)
          allocate(grid%connectinfo(grid%ncon))
          
          do j=1,grid%ncon
            call cg_1to1_read_f(ifile,1,n+1,j,connectname,donorname,ipnts,ipntsdonor,transvec,ier)
            
            grid%connectinfo(j)%transmat = 0
            grid%connectinfo(j)%transmat(abs(transvec(1)),1)=transvec(1)/abs(transvec(1))
            grid%connectinfo(j)%transmat(abs(transvec(2)),2)=transvec(2)/abs(transvec(2))
            grid%connectinfo(j)%transmat(abs(transvec(3)),3)=transvec(3)/abs(transvec(3))
            
            do m=1,dim
              if(ipnts(m).eq.ipnts(dim+m)) then
                if(ipnts(m).eq.1) then
                  select case(m)
                  case(1)
                    grid%connectinfo(j)%istart(1) = ipnts(1)+1
                    grid%connectinfo(j)%istart(2) = ipnts(2)
                    grid%connectinfo(j)%istart(3) = ipnts(3)
                    grid%connectinfo(j)%iend(1)   = ipnts(4)+1
                    grid%connectinfo(j)%iend(2)   = ipnts(5)
                    grid%connectinfo(j)%iend(3)   = ipnts(6)
                  case(2)
                    grid%connectinfo(j)%istart(1) = ipnts(1)
                    grid%connectinfo(j)%istart(2) = ipnts(2)+1
                    grid%connectinfo(j)%istart(3) = ipnts(3)
                    grid%connectinfo(j)%iend(1)   = ipnts(4)
                    grid%connectinfo(j)%iend(2)   = ipnts(5)+1
                    grid%connectinfo(j)%iend(3)   = ipnts(6)
                  case(3)
                    grid%connectinfo(j)%istart(1) = ipnts(1)
                    grid%connectinfo(j)%istart(2) = ipnts(2)
                    grid%connectinfo(j)%istart(3) = ipnts(3)+1
                    grid%connectinfo(j)%iend(1)   = ipnts(4)
                    grid%connectinfo(j)%iend(2)   = ipnts(5)
                    grid%connectinfo(j)%iend(3)   = ipnts(6)+1
                  end select
                else
                  grid%connectinfo(j)%istart(1) = ipnts(1)
                  grid%connectinfo(j)%istart(2) = ipnts(2)
                  grid%connectinfo(j)%istart(3) = ipnts(3)
                  grid%connectinfo(j)%iend(1)   = ipnts(4)
                  grid%connectinfo(j)%iend(2)   = ipnts(5)
                  grid%connectinfo(j)%iend(3)   = ipnts(6)
                end if
              end if
              
              if(ipntsdonor(m).eq.ipntsdonor(dim+m)) then
                if(ipntsdonor(m).eq.1) then
                  grid%connectinfo(j)%istart_donor(1) = ipntsdonor(1)
                  grid%connectinfo(j)%istart_donor(2) = ipntsdonor(2)
                  grid%connectinfo(j)%istart_donor(3) = ipntsdonor(3)
                  grid%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)
                  grid%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)
                  grid%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)
                else
                  select case(m)
                  case(1)
                    grid%connectinfo(j)%istart_donor(1) = ipntsdonor(1)+1
                    grid%connectinfo(j)%istart_donor(2) = ipntsdonor(2)
                    grid%connectinfo(j)%istart_donor(3) = ipntsdonor(3)
                    grid%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)+1
                    grid%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)
                    grid%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)
                  case(2)
                    grid%connectinfo(j)%istart_donor(1) = ipntsdonor(1)
                    grid%connectinfo(j)%istart_donor(2) = ipntsdonor(2)+1
                    grid%connectinfo(j)%istart_donor(3) = ipntsdonor(3)
                    grid%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)
                    grid%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)+1
                    grid%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)
                  case(3)
                    grid%connectinfo(j)%istart_donor(1) = ipntsdonor(1)
                    grid%connectinfo(j)%istart_donor(2) = ipntsdonor(2)
                    grid%connectinfo(j)%istart_donor(3) = ipntsdonor(3)+1
                    grid%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)
                    grid%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)
                    grid%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)+1
                  end select
                end if
              end if
            end do

            do i=1,nzone
              if(trim(donorname).eq.trim(zonename(i))) then
                grid%connectinfo(j)%donor = i-1
              end if
            end do
          end do  
          
          call cg_close_f(ifile,ier)
        end if
        call mpi_barrier(mpi_comm_world,ier)
      end do
            
      if(allocated(zonename)) deallocate(zonename)
      if(allocated(famname)) deallocate(famname)
      if(allocated(fambcname)) deallocate(fambcname)
      if(allocated(xx)) deallocate(xx)
      if(allocated(yy)) deallocate(yy)
      if(allocated(zz)) deallocate(zz)

      
      select case(config%getiturb())
      case(-3,-2)
        grid%ngrd = 4
      case(-1,0)
        grid%ngrd = 5
      end select
      
      allocate(grid%grd(grid%ngrd,grid%imax+1,grid%jmax+1,grid%kmax+1))
      
      call grid%calnormal() !! must be called before calvolume !!
      call grid%calvolume() 
      if(config%getiturb().gt.-2) call grid%calydns() ! in case of turbulent, ydns must be calculated !!
      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(grid)
      implicit none
      class(t_grid), intent(inout) :: grid

      if(allocated(grid%zoneinfo))    deallocate(grid%zoneinfo)
      if(allocated(grid%x))           deallocate(grid%x)
      if(allocated(grid%cx))          deallocate(grid%cx)
      if(allocated(grid%ex))          deallocate(grid%ex)
      if(allocated(grid%tx))          deallocate(grid%tx)
      if(allocated(grid%grd))         deallocate(grid%grd)
      if(allocated(grid%bcinfo))      deallocate(grid%bcinfo)
      if(allocated(grid%connectinfo)) deallocate(grid%connectinfo)
      
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calydns(grid)
      implicit none
      class(t_grid), intent(inout) :: grid
      integer :: m,n,i,j,k,l,ii,jj,kk
      integer :: sendnum,recvnum
      integer :: ier
      integer :: status(mpi_status_size)
      integer :: recvcount(0:grid%size-1),disp(0:grid%size-1)
      integer, parameter :: dim = 3
      real(8), dimension(:), allocatable :: sendxb,sendyb,sendzb,recvxb,recvyb,recvzb
      type(t_mpitemp), dimension(:), allocatable :: mpitemp
      integer :: request_s(grid%ncon),request_r(grid%ncon),request_sa(grid%ncon),request_ra(grid%ncon)
      real(8) :: xc,yc,zc,ymin,ytemp
     
      sendnum = 0
     
      do n = 1,grid%nbc
        if(trim(grid%bcinfo(n)%bcname).eq.'BCWall') then
          if(grid%bcinfo(n)%istart(1).eq.grid%bcinfo(n)%iend(1)) then
            sendnum = sendnum + (grid%bcinfo(n)%iend(2)-grid%bcinfo(n)%istart(2)+1)*(grid%bcinfo(n)%iend(3)-grid%bcinfo(n)%istart(3)+1)
          end if
          if(grid%bcinfo(n)%istart(2).eq.grid%bcinfo(n)%iend(2)) then
            sendnum = sendnum + (grid%bcinfo(n)%iend(1)-grid%bcinfo(n)%istart(1)+1)*(grid%bcinfo(n)%iend(3)-grid%bcinfo(n)%istart(3)+1)
          end if
          if(grid%bcinfo(n)%istart(3).eq.grid%bcinfo(n)%iend(3)) then
            sendnum = sendnum + (grid%bcinfo(n)%iend(1)-grid%bcinfo(n)%istart(1)+1)*(grid%bcinfo(n)%iend(2)-grid%bcinfo(n)%istart(2)+1)
          end if
        end if
      end do
      
      allocate(sendxb(sendnum),sendyb(sendnum),sendzb(sendnum))
      m=0
      do n = 1,grid%nbc
        if(trim(grid%bcinfo(n)%bcname).eq.'BCWall') then
          if(grid%bcinfo(n)%istart(1).eq.grid%bcinfo(n)%iend(1)) then
            if(grid%bcinfo(n)%istart(1).eq.1) then
              do k = grid%bcinfo(n)%istart(3),grid%bcinfo(n)%iend(3)
                do j = grid%bcinfo(n)%istart(2),grid%bcinfo(n)%iend(2)
                  m = m + 1
                  sendxb(m) = 0.25d0*(grid%x(1,2,j,k)+grid%x(1,2,j+1,k)+grid%x(1,2,j,k+1)+grid%x(1,2,j+1,k+1))
                  sendyb(m) = 0.25d0*(grid%x(2,2,j,k)+grid%x(2,2,j+1,k)+grid%x(2,2,j,k+1)+grid%x(2,2,j+1,k+1))
                  sendzb(m) = 0.25d0*(grid%x(3,2,j,k)+grid%x(3,2,j+1,k)+grid%x(3,2,j,k+1)+grid%x(3,2,j+1,k+1))
                end do
              end do
            else
              do k = grid%bcinfo(n)%istart(3),grid%bcinfo(n)%iend(3)
                do j = grid%bcinfo(n)%istart(2),grid%bcinfo(n)%iend(2)
                  m = m + 1
                  sendxb(m) = 0.25d0*(grid%x(1,grid%imax+1,j,k)+grid%x(1,grid%imax+1,j+1,k)+grid%x(1,grid%imax+1,j,k+1)+grid%x(1,grid%imax+1,j+1,k+1))
                  sendyb(m) = 0.25d0*(grid%x(2,grid%imax+1,j,k)+grid%x(2,grid%imax+1,j+1,k)+grid%x(2,grid%imax+1,j,k+1)+grid%x(2,grid%imax+1,j+1,k+1))
                  sendzb(m) = 0.25d0*(grid%x(3,grid%imax+1,j,k)+grid%x(3,grid%imax+1,j+1,k)+grid%x(3,grid%imax+1,j,k+1)+grid%x(3,grid%imax+1,j+1,k+1))
                end do
              end do           
            end if
          end if
          if(grid%bcinfo(n)%istart(2).eq.grid%bcinfo(n)%iend(2)) then
            if(grid%bcinfo(n)%istart(2).eq.1) then
              do k = grid%bcinfo(n)%istart(3),grid%bcinfo(n)%iend(3)
                do i = grid%bcinfo(n)%istart(1),grid%bcinfo(n)%iend(1)
                  m = m + 1
                  sendxb(m) = 0.25d0*(grid%x(1,i,2,k)+grid%x(1,i+1,2,k)+grid%x(1,i,2,k+1)+grid%x(1,i+1,2,k+1))
                  sendyb(m) = 0.25d0*(grid%x(2,i,2,k)+grid%x(2,i+1,2,k)+grid%x(2,i,2,k+1)+grid%x(2,i+1,2,k+1))
                  sendzb(m) = 0.25d0*(grid%x(3,i,2,k)+grid%x(3,i+1,2,k)+grid%x(3,i,2,k+1)+grid%x(3,i+1,2,k+1))
                end do
              end do            
            else
              do k = grid%bcinfo(n)%istart(3),grid%bcinfo(n)%iend(3)
                do i = grid%bcinfo(n)%istart(1),grid%bcinfo(n)%iend(1)
                  m = m + 1
                  sendxb(m) = 0.25d0*(grid%x(1,i,grid%jmax+1,k)+grid%x(1,i+1,grid%jmax+1,k)+grid%x(1,i,grid%jmax+1,k+1)+grid%x(1,i+1,grid%jmax+1,k+1))
                  sendyb(m) = 0.25d0*(grid%x(2,i,grid%jmax+1,k)+grid%x(2,i+1,grid%jmax+1,k)+grid%x(2,i,grid%jmax+1,k+1)+grid%x(2,i+1,grid%jmax+1,k+1))
                  sendzb(m) = 0.25d0*(grid%x(3,i,grid%jmax+1,k)+grid%x(3,i+1,grid%jmax+1,k)+grid%x(3,i,grid%jmax+1,k+1)+grid%x(3,i+1,grid%jmax+1,k+1))
                end do
              end do          
            end if
          end if
          if(grid%bcinfo(n)%istart(3).eq.grid%bcinfo(n)%iend(3)) then
            if(grid%bcinfo(n)%istart(3).eq.1) then
              do j = grid%bcinfo(n)%istart(2),grid%bcinfo(n)%iend(2)
                do i = grid%bcinfo(n)%istart(1),grid%bcinfo(n)%iend(1)
                  m = m + 1
                  sendxb(m) = 0.25d0*(grid%x(1,i,j,2)+grid%x(1,i+1,j,2)+grid%x(1,i,j+1,2)+grid%x(1,i+1,j+1,2))
                  sendyb(m) = 0.25d0*(grid%x(2,i,j,2)+grid%x(2,i+1,j,2)+grid%x(2,i,j+1,2)+grid%x(2,i+1,j+1,2))
                  sendzb(m) = 0.25d0*(grid%x(3,i,j,2)+grid%x(3,i+1,j,2)+grid%x(3,i,j+1,2)+grid%x(3,i+1,j+1,2))
                end do
              end do            
            else
              do j = grid%bcinfo(n)%istart(2),grid%bcinfo(n)%iend(2)
                do i = grid%bcinfo(n)%istart(1),grid%bcinfo(n)%iend(1)
                  m = m + 1
                  sendxb(m) = 0.25d0*(grid%x(1,i,j,grid%kmax+1)+grid%x(1,i+1,j,grid%kmax+1)+grid%x(1,i,j+1,grid%kmax+1)+grid%x(1,i+1,j+1,grid%kmax+1))
                  sendyb(m) = 0.25d0*(grid%x(2,i,j,grid%kmax+1)+grid%x(2,i+1,j,grid%kmax+1)+grid%x(2,i,j+1,grid%kmax+1)+grid%x(2,i+1,j+1,grid%kmax+1))
                  sendzb(m) = 0.25d0*(grid%x(3,i,j,grid%kmax+1)+grid%x(3,i+1,j,grid%kmax+1)+grid%x(3,i,j+1,grid%kmax+1)+grid%x(3,i+1,j+1,grid%kmax+1))
                end do
              end do          
            end if
          end if      
        end if
      end do
      
      call mpi_allgather(sendnum,1,mpi_integer4,recvcount,1,mpi_integer4,mpi_comm_world,ier)
      
      disp(0) = 0
      do m=1,grid%size-1
        disp(m) = disp(m-1) + recvcount(m-1)
      end do
      
      recvnum = 0
      do m=0,grid%size-1
        recvnum = recvnum + recvcount(m)
      end do
      
      allocate(recvxb(recvnum),recvyb(recvnum),recvzb(recvnum))

      call mpi_allgatherv(sendxb,sendnum,mpi_real8,recvxb,recvcount,disp,mpi_real8,mpi_comm_world,ier)
      call mpi_allgatherv(sendyb,sendnum,mpi_real8,recvyb,recvcount,disp,mpi_real8,mpi_comm_world,ier)
      call mpi_allgatherv(sendzb,sendnum,mpi_real8,recvzb,recvcount,disp,mpi_real8,mpi_comm_world,ier)
      
      do k=2,grid%kmax
        do j=2,grid%jmax
          do i=2,grid%imax
            xc = 0.125d0*(grid%x(1,i,j,k)+grid%x(1,i+1,j,k)+grid%x(1,i,j+1,k)+grid%x(1,i,j,k+1)+grid%x(1,i+1,j+1,k)+grid%x(1,i,j+1,k+1)+grid%x(1,i+1,j,k+1)+grid%x(1,i+1,j+1,k+1))
            yc = 0.125d0*(grid%x(2,i,j,k)+grid%x(2,i+1,j,k)+grid%x(2,i,j+1,k)+grid%x(2,i,j,k+1)+grid%x(2,i+1,j+1,k)+grid%x(2,i,j+1,k+1)+grid%x(2,i+1,j,k+1)+grid%x(2,i+1,j+1,k+1))
            zc = 0.125d0*(grid%x(3,i,j,k)+grid%x(3,i+1,j,k)+grid%x(3,i,j+1,k)+grid%x(3,i,j,k+1)+grid%x(3,i+1,j+1,k)+grid%x(3,i,j+1,k+1)+grid%x(3,i+1,j,k+1)+grid%x(3,i+1,j+1,k+1))
            ymin = 1.d6
            do m =1,recvnum
              ytemp = dsqrt((xc-recvxb(m))**2+(yc-recvyb(m))**2+(zc-recvzb(m))**2)
              if(ytemp.lt.ymin) then
                ymin = ytemp
              end if
            end do
            grid%grd(5,i,j,k) = ymin
          end do
        end do
      end do
      
      deallocate(sendxb,sendyb,sendzb,recvxb,recvyb,recvzb)
      
      do k=2,grid%kmax
        do i=2,grid%imax
          grid%grd(5,i,1,k)      = grid%grd(5,i,2,k)
          grid%grd(5,i,grid%jmax+1,k) = grid%grd(5,i,grid%jmax,k)
        end do
      end do
      
      do k=2,grid%kmax
        do j=1,grid%jmax+1
          grid%grd(5,1,j,k)      = grid%grd(5,2,j,k)
          grid%grd(5,grid%imax+1,j,k) = grid%grd(5,grid%imax,j,k)
        end do
      end do

      do j=1,grid%jmax+1
        do i=1,grid%imax+1
          grid%grd(5,i,j,1)      = grid%grd(5,i,j,2)
          grid%grd(5,i,j,grid%kmax+1) = grid%grd(5,i,j,grid%kmax)
        end do
      end do
      
      allocate(mpitemp(grid%ncon))

      do n=1,grid%ncon
        mpitemp(n)%num = (abs(grid%connectinfo(n)%istart(3)-grid%connectinfo(n)%iend(3))+1) &
                        *(abs(grid%connectinfo(n)%istart(2)-grid%connectinfo(n)%iend(2))+1) &
                        *(abs(grid%connectinfo(n)%istart(1)-grid%connectinfo(n)%iend(1))+1)
        allocate(mpitemp(n)%sendbuf(mpitemp(n)%num),mpitemp(n)%recvbuf(mpitemp(n)%num))
      end do
      
      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then
          call mpi_irecv(mpitemp(n)%recvadress,21,mpi_integer,grid%connectinfo(n)%donor,grid%connectinfo(n)%donor+grid%size,mpi_comm_world,request_ra(n),ier)
          call mpi_irecv(mpitemp(n)%recvbuf,mpitemp(n)%num,mpi_real8,grid%connectinfo(n)%donor, grid%connectinfo(n)%donor,mpi_comm_world,request_r(n),ier)
        end if
      end do

      do n=1,grid%ncon
        l = 0
        do k=grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3)
          do j=grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2)
            do i=grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1)
              l = l + 1
              mpitemp(n)%sendbuf(l) = grid%grd(5,i,j,k)
            end do
          end do
        end do
        if(grid%rank.ne.grid%connectinfo(n)%donor) then 

          mpitemp(n)%sendadress = (/grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1), &
                                    grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2), &
                                    grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3), &
                                    grid%connectinfo(n)%istart_donor(1),grid%connectinfo(n)%iend_donor(1), &
                                    grid%connectinfo(n)%istart_donor(2),grid%connectinfo(n)%iend_donor(2), &
                                    grid%connectinfo(n)%istart_donor(3),grid%connectinfo(n)%iend_donor(3), &  
                                    grid%connectinfo(n)%transmat(1,1),grid%connectinfo(n)%transmat(1,2),grid%connectinfo(n)%transmat(1,3), &
                                    grid%connectinfo(n)%transmat(2,1),grid%connectinfo(n)%transmat(2,2),grid%connectinfo(n)%transmat(2,3), &
                                    grid%connectinfo(n)%transmat(3,1),grid%connectinfo(n)%transmat(3,2),grid%connectinfo(n)%transmat(3,3)/)
                     
          call mpi_isend(mpitemp(n)%sendadress,21,mpi_integer,grid%connectinfo(n)%donor,grid%rank+grid%size,mpi_comm_world,request_sa(n),ier)
          call mpi_isend(mpitemp(n)%sendbuf,mpitemp(n)%num,mpi_real8,grid%connectinfo(n)%donor,grid%rank,mpi_comm_world,request_s(n),ier)
        else ! no mpi boundary
          l = 0
          do k=grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3)
            do j=grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2)
              do i=grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1)
                ii = grid%connectinfo(n)%transmat(1,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(1,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(1,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(1)
                jj = grid%connectinfo(n)%transmat(2,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(2,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(2,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(2)
                kk = grid%connectinfo(n)%transmat(3,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(3,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(3,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(3)
                l = l + 1
                grid%grd(5,ii,jj,kk) = mpitemp(n)%sendbuf(l)
              end do
            end do
          end do
        end if    
      end do
      
      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then
          call mpi_wait(request_s(n),status,ier)
          call mpi_wait(request_r(n),status,ier)
          call mpi_wait(request_sa(n),status,ier)
          call mpi_wait(request_ra(n),status,ier)
        end if
      end do

      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then 
          l = 0
          do k=mpitemp(n)%recvadress(5),mpitemp(n)%recvadress(6)
            do j=mpitemp(n)%recvadress(3),mpitemp(n)%recvadress(4)
              do i=mpitemp(n)%recvadress(1),mpitemp(n)%recvadress(2)
                ii = mpitemp(n)%recvadress(13)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(14)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(15)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(7)
                jj = mpitemp(n)%recvadress(16)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(17)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(18)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(9)
                kk = mpitemp(n)%recvadress(19)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(20)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(21)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(11)
                l = l + 1
                grid%grd(5,ii,jj,kk) = mpitemp(n)%recvbuf(l)
              end do
            end do
          end do
        end if
      end do
      
      do n=1,grid%ncon
        deallocate(mpitemp(n)%sendbuf,mpitemp(n)%recvbuf)
      end do
      deallocate(mpitemp)
    end subroutine calydns
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calnormal(grid)
      implicit none
      class(t_grid), intent(inout) :: grid
      integer :: i,j,k,n,m,l,ii,jj,kk
      integer :: sign1,sign2
      integer :: ier,request_s(grid%ncon),request_r(grid%ncon),request_sa(grid%ncon),request_ra(grid%ncon)
      integer :: status(mpi_status_size)
      type(t_mpitemp), dimension(:), allocatable :: mpitemp
      real(8) :: xc,yc,zc,xe,ye,ze,xs,ys,zs 
      integer, parameter :: dim = 3
      
      allocate(grid%cx(3,1:grid%imax,1:grid%jmax+1,1:grid%kmax+1))
      allocate(grid%ex(3,1:grid%imax+1,1:grid%jmax,1:grid%kmax+1))
      allocate(grid%tx(3,1:grid%imax+1,1:grid%jmax+1,1:grid%kmax))
      
      do k=2,grid%kmax
        do j=2,grid%jmax
          do i=1,grid%imax
            xe = grid%x(1,i+1,j+1,k)   - grid%x(1,i+1,j,k+1)
            ye = grid%x(2,i+1,j+1,k)   - grid%x(2,i+1,j,k+1)
            ze = grid%x(3,i+1,j+1,k)   - grid%x(3,i+1,j,k+1)
            xs = grid%x(1,i+1,j+1,k+1) - grid%x(1,i+1,j,k)
            ys = grid%x(2,i+1,j+1,k+1) - grid%x(2,i+1,j,k)
            zs = grid%x(3,i+1,j+1,k+1) - grid%x(3,i+1,j,k)
            grid%cx(1,i,j,k) = 0.5d0*(ye*zs - ys*ze)
            grid%cx(2,i,j,k) = 0.5d0*(xs*ze - xe*zs)
            grid%cx(3,i,j,k) = 0.5d0*(xe*ys - xs*ye)
          end do
        end do
      end do

      do i=2,grid%imax
        do k=2,grid%kmax
          do j=1,grid%jmax
            xc = grid%x(1,i,j+1,k)   - grid%x(1,i+1,j+1,k+1)
            yc = grid%x(2,i,j+1,k)   - grid%x(2,i+1,j+1,k+1)
            zc = grid%x(3,i,j+1,k)   - grid%x(3,i+1,j+1,k+1)    
            xs = grid%x(1,i,j+1,k+1) - grid%x(1,i+1,j+1,k)
            ys = grid%x(2,i,j+1,k+1) - grid%x(2,i+1,j+1,k)
            zs = grid%x(3,i,j+1,k+1) - grid%x(3,i+1,j+1,k)
            grid%ex(1,i,j,k) = 0.5d0*(yc*zs - ys*zc)
            grid%ex(2,i,j,k) = 0.5d0*(xs*zc - xc*zs)
            grid%ex(3,i,j,k) = 0.5d0*(xc*ys - xs*yc)   
          end do
        end do
      end do

      do j=2,grid%jmax
        do i=2,grid%imax
          do k=1,grid%kmax
            xc = grid%x(1,i+1,j,k+1)   - grid%x(1,i,j+1,k+1)
            yc = grid%x(2,i+1,j,k+1)   - grid%x(2,i,j+1,k+1)
            zc = grid%x(3,i+1,j,k+1)   - grid%x(3,i,j+1,k+1)
            xe = grid%x(1,i+1,j+1,k+1) - grid%x(1,i,j,k+1)
            ye = grid%x(2,i+1,j+1,k+1) - grid%x(2,i,j,k+1)
            ze = grid%x(3,i+1,j+1,k+1) - grid%x(3,i,j,k+1)
            grid%tx(1,i,j,k) = 0.5d0*(yc*ze - ye*zc)
            grid%tx(2,i,j,k) = 0.5d0*(xe*zc - xc*ze)
            grid%tx(3,i,j,k) = 0.5d0*(xc*ye - xe*yc) 
          end do
        end do
      end do
 
      do n = 1,grid%nbc
        if((trim(grid%bcinfo(n)%bcname).eq.'BCDegeneratePoint').or.(trim(grid%bcinfo(n)%bcname).eq.'BCDegenerateLine')) then
          if(grid%bcinfo(n)%istart(1).eq.grid%bcinfo(n)%iend(1)) then
            if(grid%bcinfo(n)%istart(1).eq.1) then
              grid%cx(1,1,:,:) = 0.d0
              grid%cx(2,1,:,:) = 0.d0
              grid%cx(3,1,:,:) = 0.d0
            else
              grid%cx(1,grid%imax,:,:) = 0.d0
              grid%cx(2,grid%imax,:,:) = 0.d0
              grid%cx(3,grid%imax,:,:) = 0.d0
            end if
          end if
          if(grid%bcinfo(n)%istart(2).eq.grid%bcinfo(n)%iend(2)) then
            if(grid%bcinfo(n)%istart(2).eq.1) then
              grid%ex(1,:,1,:) = 0.d0
              grid%ex(2,:,1,:) = 0.d0
              grid%ex(3,:,1,:) = 0.d0
            else
              grid%ex(1,:,grid%jmax,:) = 0.d0
              grid%ex(2,:,grid%jmax,:) = 0.d0
              grid%ex(3,:,grid%jmax,:) = 0.d0
            end if
          end if
          if(grid%bcinfo(n)%istart(3).eq.grid%bcinfo(n)%iend(3)) then
            if(grid%bcinfo(n)%istart(3).eq.1) then
              grid%tx(1,:,:,1) = 0.d0
              grid%tx(2,:,:,1) = 0.d0
              grid%tx(3,:,:,1) = 0.d0
            else
              grid%tx(1,:,:,grid%kmax) = 0.d0
              grid%tx(2,:,:,grid%kmax) = 0.d0
              grid%tx(3,:,:,grid%kmax) = 0.d0
            end if
          end if 
        end if
      end do
      
      grid%cx(:,:,1,:) = grid%cx(:,:,2,:)
      grid%cx(:,:,grid%jmax+1,:) = grid%cx(:,:,grid%jmax,:)
      grid%cx(:,:,:,1) = grid%cx(:,:,:,2)
      grid%cx(:,:,:,grid%kmax+1) = grid%cx(:,:,:,grid%kmax)
      
      grid%ex(:,1,:,:) = grid%ex(:,2,:,:)
      grid%ex(:,grid%imax+1,:,:) = grid%ex(:,grid%imax,:,:)
      grid%ex(:,:,:,1) = grid%ex(:,:,:,2)
      grid%ex(:,:,:,grid%kmax+1) = grid%ex(:,:,:,grid%kmax)
      
      grid%tx(:,1,:,:) = grid%tx(:,2,:,:)
      grid%tx(:,grid%imax+1,:,:) = grid%tx(:,grid%imax,:,:)
      grid%tx(:,:,1,:) = grid%tx(:,:,2,:)
      grid%tx(:,:,grid%jmax+1,:) = grid%tx(:,:,grid%jmax,:)
      
      allocate(mpitemp(grid%ncon))
      
      do n=1,grid%ncon
        mpitemp(n)%num = 2*dim*(abs(grid%connectinfo(n)%istart(3)-grid%connectinfo(n)%iend(3))+1) &
                              *(abs(grid%connectinfo(n)%istart(2)-grid%connectinfo(n)%iend(2))+1) &
                              *(abs(grid%connectinfo(n)%istart(1)-grid%connectinfo(n)%iend(1))+1)
        allocate(mpitemp(n)%sendbuf(mpitemp(n)%num),mpitemp(n)%recvbuf(mpitemp(n)%num))
      end do
      
      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then
          call mpi_irecv(mpitemp(n)%recvadress,21,mpi_integer,grid%connectinfo(n)%donor,grid%connectinfo(n)%donor+grid%size,mpi_comm_world,request_ra(n),ier)
          call mpi_irecv(mpitemp(n)%recvbuf,mpitemp(n)%num,mpi_real8,grid%connectinfo(n)%donor,grid%connectinfo(n)%donor,mpi_comm_world,request_r(n),ier)
        end if
      end do

      do n=1,grid%ncon  
        l = 0
        do k=grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3)
          do j=grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2)
            do i=grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1)
              do m=1,dim
                if(grid%connectinfo(n)%istart(1).eq.grid%connectinfo(n)%iend(1)) then
                  l = l + 1
                  mpitemp(n)%sendbuf(l) = grid%ex(m,i,j,k)
                  l = l + 1
                  mpitemp(n)%sendbuf(l) = grid%tx(m,i,j,k)
                else if(grid%connectinfo(n)%istart(2).eq.grid%connectinfo(n)%iend(2)) then
                  l = l + 1
                  mpitemp(n)%sendbuf(l) = grid%tx(m,i,j,k)
                  l = l + 1
                  mpitemp(n)%sendbuf(l) = grid%cx(m,i,j,k)
                else if(grid%connectinfo(n)%istart(3).eq.grid%connectinfo(n)%iend(3)) then
                  l = l + 1
                  mpitemp(n)%sendbuf(l) = grid%cx(m,i,j,k)
                  l = l + 1
                  mpitemp(n)%sendbuf(l) = grid%ex(m,i,j,k)
                end if
              end do
            end do
          end do
        end do
        if(grid%rank.ne.grid%connectinfo(n)%donor) then 

          mpitemp(n)%sendadress = (/grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1), &
                                    grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2), &
                                    grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3), &
                                    grid%connectinfo(n)%istart_donor(1),grid%connectinfo(n)%iend_donor(1), &
                                    grid%connectinfo(n)%istart_donor(2),grid%connectinfo(n)%iend_donor(2), &
                                    grid%connectinfo(n)%istart_donor(3),grid%connectinfo(n)%iend_donor(3), &  
                                    grid%connectinfo(n)%transmat(1,1),grid%connectinfo(n)%transmat(1,2),grid%connectinfo(n)%transmat(1,3), &
                                    grid%connectinfo(n)%transmat(2,1),grid%connectinfo(n)%transmat(2,2),grid%connectinfo(n)%transmat(2,3), &
                                    grid%connectinfo(n)%transmat(3,1),grid%connectinfo(n)%transmat(3,2),grid%connectinfo(n)%transmat(3,3)/)
                     
          call mpi_isend(mpitemp(n)%sendadress,21,mpi_integer,grid%connectinfo(n)%donor,grid%rank+grid%size,mpi_comm_world,request_sa(n),ier)
          call mpi_isend(mpitemp(n)%sendbuf,mpitemp(n)%num,mpi_real8,grid%connectinfo(n)%donor,grid%rank,mpi_comm_world,request_s(n),ier)
        else ! no mpi boundary
          
          sign1 = 1
          sign2 = 1
          
          if(grid%connectinfo(n)%istart_donor(1).eq.grid%connectinfo(n)%iend_donor(1)) then
            if(grid%connectinfo(n)%istart_donor(2).gt.grid%connectinfo(n)%iend_donor(2)) then
              sign1 = -1*sign1
            end if          
            if(grid%connectinfo(n)%istart_donor(3).gt.grid%connectinfo(n)%iend_donor(3)) then
              sign2 = -1*sign2
            end if  
          else if(grid%connectinfo(n)%istart_donor(2).eq.grid%connectinfo(n)%iend_donor(2)) then
            if(grid%connectinfo(n)%istart_donor(3).gt.grid%connectinfo(n)%iend_donor(3)) then
              sign1 = -1*sign1
            end if 
            if(grid%connectinfo(n)%istart_donor(1).gt.grid%connectinfo(n)%iend_donor(1)) then
              sign2 = -1*sign2
            end if
          else if(grid%connectinfo(n)%istart_donor(3).eq.grid%connectinfo(n)%iend_donor(3)) then
            if(grid%connectinfo(n)%istart_donor(1).gt.grid%connectinfo(n)%iend_donor(1)) then
              sign1 = -1*sign1
            end if
            if(grid%connectinfo(n)%istart_donor(2).gt.grid%connectinfo(n)%iend_donor(2)) then
              sign2 = -1*sign2
            end if
          end if
          
          l = 0
          do k=grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3)
            do j=grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2)
              do i=grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1)
                ii = grid%connectinfo(n)%transmat(1,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(1,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(1,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(1)
                jj = grid%connectinfo(n)%transmat(2,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(2,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(2,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(2)
                kk = grid%connectinfo(n)%transmat(3,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(3,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(3,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(3)
                do m=1,dim
                  if(grid%connectinfo(n)%istart_donor(1).eq.grid%connectinfo(n)%iend_donor(1)) then
                    l = l + 1
                    grid%ex(m,ii,jj,kk) = sign1*mpitemp(n)%sendbuf(l)
                    l = l + 1
                    grid%tx(m,ii,jj,kk) = sign2*mpitemp(n)%sendbuf(l)
                  else if(grid%connectinfo(n)%istart_donor(2).eq.grid%connectinfo(n)%iend_donor(2)) then
                    l = l + 1
                    grid%tx(m,ii,jj,kk) = sign1*mpitemp(n)%sendbuf(l)
                    l = l + 1
                    grid%cx(m,ii,jj,kk) = sign2*mpitemp(n)%sendbuf(l)
                  else if(grid%connectinfo(n)%istart_donor(3).eq.grid%connectinfo(n)%iend_donor(3)) then
                    l = l + 1
                    grid%cx(m,ii,jj,kk) = sign1*mpitemp(n)%sendbuf(l)
                    l = l + 1
                    grid%ex(m,ii,jj,kk) = sign2*mpitemp(n)%sendbuf(l)
                  end if
                end do
              end do
            end do
          end do
        end if    
      end do

      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then
          call mpi_wait(request_s(n),status,ier)
          call mpi_wait(request_r(n),status,ier)
          call mpi_wait(request_sa(n),status,ier)
          call mpi_wait(request_ra(n),status,ier)
        end if
      end do

      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then 
          sign1 = 1
          sign2 = 1
          
          if(mpitemp(n)%recvadress(7).eq.mpitemp(n)%recvadress(8)) then
            if(mpitemp(n)%recvadress(9).gt.mpitemp(n)%recvadress(10)) then
              sign1 = -1*sign1
            end if          
            if(mpitemp(n)%recvadress(11).gt.mpitemp(n)%recvadress(12)) then
              sign2 = -1*sign2
            end if  
          else if(mpitemp(n)%recvadress(9).eq.mpitemp(n)%recvadress(10)) then
            if(mpitemp(n)%recvadress(11).gt.mpitemp(n)%recvadress(12)) then
              sign1 = -1*sign1
            end if 
            if(mpitemp(n)%recvadress(7).gt.mpitemp(n)%recvadress(8)) then
              sign2 = -1*sign2
            end if
          else if(mpitemp(n)%recvadress(11).eq.mpitemp(n)%recvadress(12)) then
            if(mpitemp(n)%recvadress(7).gt.mpitemp(n)%recvadress(8)) then
              sign1 = -1*sign1
            end if
            if(mpitemp(n)%recvadress(9).gt.mpitemp(n)%recvadress(10)) then
              sign2 = -1*sign2
            end if
          end if
          
          l = 0
          do k=mpitemp(n)%recvadress(5),mpitemp(n)%recvadress(6)
            do j=mpitemp(n)%recvadress(3),mpitemp(n)%recvadress(4)
              do i=mpitemp(n)%recvadress(1),mpitemp(n)%recvadress(2)
                ii = mpitemp(n)%recvadress(13)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(14)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(15)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(7)
                jj = mpitemp(n)%recvadress(16)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(17)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(18)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(9)
                kk = mpitemp(n)%recvadress(19)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(20)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(21)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(11)
                do m=1,dim
                  if(mpitemp(n)%recvadress(7).eq.mpitemp(n)%recvadress(8)) then
                    l = l + 1
                    grid%ex(m,ii,jj,kk) = sign1*mpitemp(n)%recvbuf(l)
                    l = l + 1
                    grid%tx(m,ii,jj,kk) = sign2*mpitemp(n)%recvbuf(l)
                  else if(mpitemp(n)%recvadress(9).eq.mpitemp(n)%recvadress(10)) then
                    l = l + 1
                    grid%tx(m,ii,jj,kk) = sign1*mpitemp(n)%recvbuf(l)
                    l = l + 1
                    grid%cx(m,ii,jj,kk) = sign2*mpitemp(n)%recvbuf(l)
                  else if(mpitemp(n)%recvadress(11).eq.mpitemp(n)%recvadress(12)) then
                    l = l + 1
                    grid%cx(m,ii,jj,kk) = sign1*mpitemp(n)%recvbuf(l)
                    l = l + 1
                    grid%ex(m,ii,jj,kk) = sign2*mpitemp(n)%recvbuf(l)
                  end if
                end do
              end do
            end do
          end do
        end if
      end do
      
      do n=1,grid%ncon
        deallocate(mpitemp(n)%sendbuf,mpitemp(n)%recvbuf)
      end do
      deallocate(mpitemp)

      ! reodering
      do n=1,grid%ncon
        do m=1,dim
          if(grid%connectinfo(n)%istart(m).ne.grid%connectinfo(n)%iend(m) ) then
            grid%connectinfo(n)%istart(m) = grid%connectinfo(n)%istart(m) + 1            
          end if

          if(grid%connectinfo(n)%istart_donor(m).ne.grid%connectinfo(n)%iend_donor(m) ) then
            if(grid%connectinfo(n)%istart_donor(m).lt.grid%connectinfo(n)%iend_donor(m) ) then
              grid%connectinfo(n)%istart_donor(m) = grid%connectinfo(n)%istart_donor(m) + 1
            else if(grid%connectinfo(n)%istart_donor(m).gt.grid%connectinfo(n)%iend_donor(m) ) then
              grid%connectinfo(n)%iend_donor(m) = grid%connectinfo(n)%iend_donor(m) + 1            
            end if
          end if
        end do
      end do
      
    end subroutine calnormal
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calvolume(grid)
      implicit none
      class(t_grid), intent(inout) :: grid
      integer :: i,j,k,n,l,ii,jj,kk
      real(8) :: pa(3),pb(3),pc(3),pd(3)
      real(8) :: volp1,volp2,volp3,volp4,volp5,volp6
      integer :: ier,request_s(grid%ncon),request_r(grid%ncon),request_sa(grid%ncon),request_ra(grid%ncon)
      integer :: status(mpi_status_size)
      type(t_mpitemp), dimension(:), allocatable :: mpitemp

      do k=2,grid%kmax
        do j=2,grid%jmax
          do i=2,grid%imax
            grid%grd(2,i,j,k) = 0.125d0*(grid%x(1,i,j,k)+grid%x(1,i+1,j,k)+grid%x(1,i,j+1,k)+grid%x(1,i,j,k+1) &
                              + grid%x(1,i+1,j+1,k)+grid%x(1,i,j+1,k+1)+grid%x(1,i+1,j,k+1)+grid%x(1,i+1,j+1,k+1))
            grid%grd(3,i,j,k) = 0.125d0*(grid%x(2,i,j,k)+grid%x(2,i+1,j,k)+grid%x(2,i,j+1,k)+grid%x(2,i,j,k+1) &
                              + grid%x(2,i+1,j+1,k)+grid%x(2,i,j+1,k+1)+grid%x(2,i+1,j,k+1)+grid%x(2,i+1,j+1,k+1))
            grid%grd(4,i,j,k) = 0.125d0*(grid%x(3,i,j,k)+grid%x(3,i+1,j,k)+grid%x(3,i,j+1,k)+grid%x(3,i,j,k+1) &
                              + grid%x(3,i+1,j+1,k)+grid%x(3,i,j+1,k+1)+grid%x(3,i+1,j,k+1)+grid%x(3,i+1,j+1,k+1))
                              
            pa(1) = grid%x(1,i,j,k) 
            pb(1) = grid%x(1,i+1,j,k)
            pc(1) = grid%x(1,i+1,j+1,k)
            pd(1) = grid%x(1,i,j+1,k)

            pa(2) = grid%x(2,i,j,k)
            pb(2) = grid%x(2,i+1,j,k)
            pc(2) = grid%x(2,i+1,j+1,k)
            pd(2) = grid%x(2,i,j+1,k)

            pa(3) = grid%x(3,i,j,k)
            pb(3) = grid%x(3,i+1,j,k)
            pc(3) = grid%x(3,i+1,j+1,k)
            pd(3) = grid%x(3,i,j+1,k)
            
            volp1 = pyramid(pa,pb,pc,pd,grid%grd(2,i,j,k),grid%grd(3,i,j,k),grid%grd(4,i,j,k))
            
            pa(1) = grid%x(1,i,j,k+1) 
            pb(1) = grid%x(1,i+1,j,k+1)
            pc(1) = grid%x(1,i+1,j+1,k+1)
            pd(1) = grid%x(1,i,j+1,k+1)

            pa(2) = grid%x(2,i,j,k+1)
            pb(2) = grid%x(2,i+1,j,k+1)
            pc(2) = grid%x(2,i+1,j+1,k+1)
            pd(2) = grid%x(2,i,j+1,k+1)

            pa(3) = grid%x(3,i,j,k+1)
            pb(3) = grid%x(3,i+1,j,k+1)
            pc(3) = grid%x(3,i+1,j+1,k+1)
            pd(3) = grid%x(3,i,j+1,k+1)
            
            volp2 = pyramid(pa,pb,pc,pd,grid%grd(2,i,j,k),grid%grd(3,i,j,k),grid%grd(4,i,j,k))
            
            pa(1) = grid%x(1,i,j,k) 
            pb(1) = grid%x(1,i+1,j,k)
            pc(1) = grid%x(1,i+1,j,k+1)
            pd(1) = grid%x(1,i,j,k+1)

            pa(2) = grid%x(2,i,j,k)
            pb(2) = grid%x(2,i+1,j,k)
            pc(2) = grid%x(2,i+1,j,k+1)
            pd(2) = grid%x(2,i,j,k+1)

            pa(3) = grid%x(3,i,j,k)
            pb(3) = grid%x(3,i+1,j,k)
            pc(3) = grid%x(3,i+1,j,k+1)
            pd(3) = grid%x(3,i,j,k+1)
            
            volp3 = pyramid(pa,pb,pc,pd,grid%grd(2,i,j,k),grid%grd(3,i,j,k),grid%grd(4,i,j,k))
            
            pa(1) = grid%x(1,i,j+1,k) 
            pb(1) = grid%x(1,i+1,j+1,k)
            pc(1) = grid%x(1,i+1,j+1,k+1)
            pd(1) = grid%x(1,i,j+1,k+1)

            pa(2) = grid%x(2,i,j+1,k)
            pb(2) = grid%x(2,i+1,j+1,k)
            pc(2) = grid%x(2,i+1,j+1,k+1)
            pd(2) = grid%x(2,i,j+1,k+1)

            pa(3) = grid%x(3,i,j+1,k)
            pb(3) = grid%x(3,i+1,j+1,k)
            pc(3) = grid%x(3,i+1,j+1,k+1)
            pd(3) = grid%x(3,i,j+1,k+1)
            
            volp4 = pyramid(pa,pb,pc,pd,grid%grd(2,i,j,k),grid%grd(3,i,j,k),grid%grd(4,i,j,k))
            
            pa(1) = grid%x(1,i,j,k) 
            pb(1) = grid%x(1,i,j+1,k)
            pc(1) = grid%x(1,i,j+1,k+1)
            pd(1) = grid%x(1,i,j,k+1)

            pa(2) = grid%x(2,i,j,k)
            pb(2) = grid%x(2,i,j+1,k)
            pc(2) = grid%x(2,i,j+1,k+1)
            pd(2) = grid%x(2,i,j,k+1)

            pa(3) = grid%x(3,i,j,k)
            pb(3) = grid%x(3,i,j+1,k)
            pc(3) = grid%x(3,i,j+1,k+1)
            pd(3) = grid%x(3,i,j,k+1)
            
            volp5 = pyramid(pa,pb,pc,pd,grid%grd(2,i,j,k),grid%grd(3,i,j,k),grid%grd(4,i,j,k))
            
            pa(1) = grid%x(1,i+1,j,k) 
            pb(1) = grid%x(1,i+1,j+1,k)
            pc(1) = grid%x(1,i+1,j+1,k+1)
            pd(1) = grid%x(1,i+1,j,k+1)

            pa(2) = grid%x(2,i+1,j,k)
            pb(2) = grid%x(2,i+1,j+1,k)
            pc(2) = grid%x(2,i+1,j+1,k+1)
            pd(2) = grid%x(2,i+1,j,k+1)

            pa(3) = grid%x(3,i+1,j,k)
            pb(3) = grid%x(3,i+1,j+1,k)
            pc(3) = grid%x(3,i+1,j+1,k+1)
            pd(3) = grid%x(3,i+1,j,k+1)
            
            volp6 = pyramid(pa,pb,pc,pd,grid%grd(2,i,j,k),grid%grd(3,i,j,k),grid%grd(4,i,j,k))
            
            grid%grd(1,i,j,k) = (volp1+volp2+volp3+volp4+volp5+volp6)/6.d0
            if(grid%grd(1,i,j,k).lt.0.d0) write(*,*) 'negative volian',i,j,k,grid%grd(1,i,j,k)

          end do
        end do
      end do

      do k=2,grid%kmax
        do i=2,grid%imax
          grid%grd(1,i,1,k) = grid%grd(1,i,2,k)
          grid%grd(2,i,1,k) = 2.d0*grid%grd(2,i,2,k) - grid%grd(2,i,3,k)
          grid%grd(3,i,1,k) = 2.d0*grid%grd(3,i,2,k) - grid%grd(3,i,3,k)
          grid%grd(4,i,1,k) = 2.d0*grid%grd(4,i,2,k) - grid%grd(4,i,3,k)
          
          grid%grd(1,i,grid%jmax+1,k) = grid%grd(1,i,grid%jmax,k)
          grid%grd(2,i,grid%jmax+1,k) = 2.d0*grid%grd(2,i,grid%jmax,k) - grid%grd(2,i,grid%jmax-1,k)
          grid%grd(3,i,grid%jmax+1,k) = 2.d0*grid%grd(3,i,grid%jmax,k) - grid%grd(3,i,grid%jmax-1,k)
          grid%grd(4,i,grid%jmax+1,k) = 2.d0*grid%grd(4,i,grid%jmax,k) - grid%grd(4,i,grid%jmax-1,k)
        end do
      end do
        
      do k=2,grid%kmax
        do j=1,grid%jmax+1
          grid%grd(1,1,j,k) = grid%grd(1,2,j,k)
          grid%grd(2,1,j,k) = 2.d0*grid%grd(2,2,j,k) - grid%grd(2,3,j,k)
          grid%grd(3,1,j,k) = 2.d0*grid%grd(3,2,j,k) - grid%grd(3,3,j,k)
          grid%grd(4,1,j,k) = 2.d0*grid%grd(4,2,j,k) - grid%grd(4,3,j,k)
          
          grid%grd(1,grid%imax+1,j,k) = grid%grd(1,grid%imax,j,k)
          grid%grd(2,grid%imax+1,j,k) = 2.d0*grid%grd(2,grid%imax,j,k) - grid%grd(2,grid%imax-1,j,k)
          grid%grd(3,grid%imax+1,j,k) = 2.d0*grid%grd(3,grid%imax,j,k) - grid%grd(3,grid%imax-1,j,k)
          grid%grd(4,grid%imax+1,j,k) = 2.d0*grid%grd(4,grid%imax,j,k) - grid%grd(4,grid%imax-1,j,k)
        end do
      end do

      do j=1,grid%jmax+1
        do i=1,grid%imax+1
          grid%grd(1,i,j,1) = grid%grd(1,i,j,2)
          grid%grd(2,i,j,1) = 2.d0*grid%grd(2,i,j,2) - grid%grd(2,3,j,3)
          grid%grd(3,i,j,1) = 2.d0*grid%grd(3,i,j,2) - grid%grd(3,3,j,3)
          grid%grd(4,i,j,1) = 2.d0*grid%grd(4,i,j,2) - grid%grd(4,3,j,3)
          
          grid%grd(1,i,j,grid%kmax+1) = grid%grd(1,i,j,grid%kmax)
          grid%grd(2,i,j,grid%kmax+1) = 2.d0*grid%grd(2,i,j,grid%kmax) - grid%grd(2,i,j,grid%kmax-1)
          grid%grd(3,i,j,grid%kmax+1) = 2.d0*grid%grd(3,i,j,grid%kmax) - grid%grd(3,i,j,grid%kmax-1)
          grid%grd(4,i,j,grid%kmax+1) = 2.d0*grid%grd(4,i,j,grid%kmax) - grid%grd(4,i,j,grid%kmax-1)
        end do
      end do
      
      allocate(mpitemp(grid%ncon))
      
      do n=1,grid%ncon
        mpitemp(n)%num = (abs(grid%connectinfo(n)%istart(3)-grid%connectinfo(n)%iend(3))+1) &
                        *(abs(grid%connectinfo(n)%istart(2)-grid%connectinfo(n)%iend(2))+1) &
                        *(abs(grid%connectinfo(n)%istart(1)-grid%connectinfo(n)%iend(1))+1)
        allocate(mpitemp(n)%sendbuf(mpitemp(n)%num),mpitemp(n)%recvbuf(mpitemp(n)%num))
      end do

      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then
          call mpi_irecv(mpitemp(n)%recvadress,21,mpi_integer,grid%connectinfo(n)%donor,grid%connectinfo(n)%donor+grid%size,mpi_comm_world,request_ra(n),ier)
          call mpi_irecv(mpitemp(n)%recvbuf,mpitemp(n)%num,mpi_real8,grid%connectinfo(n)%donor,grid%connectinfo(n)%donor,mpi_comm_world,request_r(n),ier)
        end if
      end do

      do n=1,grid%ncon
        l = 0
        do k=grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3)
          do j=grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2)
            do i=grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1)
              l = l + 1
              mpitemp(n)%sendbuf(l) = grid%grd(1,i,j,k)
            end do
          end do
        end do
        if(grid%rank.ne.grid%connectinfo(n)%donor) then 

          mpitemp(n)%sendadress = (/grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1), &
                                    grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2), &
                                    grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3), &
                                    grid%connectinfo(n)%istart_donor(1),grid%connectinfo(n)%iend_donor(1), &
                                    grid%connectinfo(n)%istart_donor(2),grid%connectinfo(n)%iend_donor(2), &
                                    grid%connectinfo(n)%istart_donor(3),grid%connectinfo(n)%iend_donor(3), &  
                                    grid%connectinfo(n)%transmat(1,1),grid%connectinfo(n)%transmat(1,2),grid%connectinfo(n)%transmat(1,3), &
                                    grid%connectinfo(n)%transmat(2,1),grid%connectinfo(n)%transmat(2,2),grid%connectinfo(n)%transmat(2,3), &
                                    grid%connectinfo(n)%transmat(3,1),grid%connectinfo(n)%transmat(3,2),grid%connectinfo(n)%transmat(3,3)/)
                     
          call mpi_isend(mpitemp(n)%sendadress,21,mpi_integer,grid%connectinfo(n)%donor,grid%rank+grid%size,mpi_comm_world,request_sa(n),ier)
          call mpi_isend(mpitemp(n)%sendbuf,mpitemp(n)%num,mpi_real8,grid%connectinfo(n)%donor,grid%rank,mpi_comm_world,request_s(n),ier)
        else ! no mpi boundary
          l = 0
          do k=grid%connectinfo(n)%istart(3),grid%connectinfo(n)%iend(3)
            do j=grid%connectinfo(n)%istart(2),grid%connectinfo(n)%iend(2)
              do i=grid%connectinfo(n)%istart(1),grid%connectinfo(n)%iend(1)
                ii = grid%connectinfo(n)%transmat(1,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(1,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(1,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(1)
                jj = grid%connectinfo(n)%transmat(2,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(2,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(2,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(2)
                kk = grid%connectinfo(n)%transmat(3,1)*(i-grid%connectinfo(n)%istart(1)) &
                   + grid%connectinfo(n)%transmat(3,2)*(j-grid%connectinfo(n)%istart(2)) &
                   + grid%connectinfo(n)%transmat(3,3)*(k-grid%connectinfo(n)%istart(3)) + grid%connectinfo(n)%istart_donor(3)
                l = l + 1
                grid%grd(1,ii,jj,kk) = mpitemp(n)%sendbuf(l)
              end do
            end do
          end do
        end if    
      end do
      
      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then
          call mpi_wait(request_s(n),status,ier)
          call mpi_wait(request_r(n),status,ier)
          call mpi_wait(request_sa(n),status,ier)
          call mpi_wait(request_ra(n),status,ier)
        end if
      end do

      do n=1,grid%ncon
        if(grid%rank.ne.grid%connectinfo(n)%donor) then 
          l = 0
          do k=mpitemp(n)%recvadress(5),mpitemp(n)%recvadress(6)
            do j=mpitemp(n)%recvadress(3),mpitemp(n)%recvadress(4)
              do i=mpitemp(n)%recvadress(1),mpitemp(n)%recvadress(2)
                ii = mpitemp(n)%recvadress(13)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(14)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(15)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(7)
                jj = mpitemp(n)%recvadress(16)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(17)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(18)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(9)
                kk = mpitemp(n)%recvadress(19)*(i-mpitemp(n)%recvadress(1)) &
                   + mpitemp(n)%recvadress(20)*(j-mpitemp(n)%recvadress(3)) &
                   + mpitemp(n)%recvadress(21)*(k-mpitemp(n)%recvadress(5)) + mpitemp(n)%recvadress(11)
                l = l + 1
                grid%grd(1,ii,jj,kk) = mpitemp(n)%recvbuf(l)
              end do
            end do
          end do
        end if
      end do
      
      do n=1,grid%ncon
        deallocate(mpitemp(n)%sendbuf,mpitemp(n)%recvbuf)
      end do
      deallocate(mpitemp)

    end subroutine calvolume
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function pyramid(pa,pb,pc,pd,xc,yc,zc)
      implicit none
      real(8), intent(in) :: pa(3) ,pb(3) ,pc(3) ,pd(3) 
      real(8), intent(in) :: xc,yc,zc
      real(8) :: pyramid
      real(8) :: csx,csy,csz,vax,vay,vaz,vbx,vby,vbz,volpx,volpy,volpz
      
      csx = 0.25d0*( pa(1)+pb(1)+pc(1)+pd(1) )
      csy = 0.25d0*( pa(2)+pb(2)+pc(2)+pd(2) )
      csz = 0.25d0*( pa(3)+pb(3)+pc(3)+pd(3) )
      
      vax = pc(1) - pa(1)
      vay = pc(2) - pa(2)
      vaz = pc(3) - pa(3)
      
      vbx = pd(1) - pb(1)
      vby = pd(2) - pb(2)
      vbz = pd(3) - pb(3)
      
      volpx = (xc - csx)*(vay*vbz - vaz*vby)
      volpy = (yc - csy)*(vaz*vbx - vax*vbz)
      volpz = (zc - csz)*(vax*vby - vay*vbx)
      
      !  make right hand rule unnecessary by take absolute value
      pyramid = dabs( volpx + volpy + volpz )
    end function pyramid
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getzonename(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      character(32) :: getzonename

      getzonename = grid%zonename
    end function getzonename
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getimax(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getimax
      
      getimax = grid%imax
    end function getimax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getjmax(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getjmax
      
      getjmax = grid%jmax
    end function getjmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getkmax(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getkmax
      
      getkmax = grid%kmax
    end function getkmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getimax_zone(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getimax_zone

      getimax_zone = grid%zoneinfo(i)%imax
    end function getimax_zone
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getjmax_zone(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getjmax_zone

      getjmax_zone = grid%zoneinfo(i)%jmax
    end function getjmax_zone
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getkmax_zone(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getkmax_zone

      getkmax_zone = grid%zoneinfo(i)%kmax
    end function getkmax_zone
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getnbc(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getnbc

      getnbc = grid%nbc
    end function getnbc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
    pure function getncon(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getncon
      
      getncon = grid%ncon
    end function getncon
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getx(grid,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getx(3)
      integer, intent(in) :: i,j,k
      
      getx = grid%x(:,i,j,k)
    end function getx
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getcx(grid,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getcx(3)
      integer, intent(in) :: i,j,k
      
      getcx = grid%cx(:,i,j,k)
    end function getcx
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getex(grid,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getex(3)
      integer, intent(in) :: i,j,k
      
      getex = grid%ex(:,i,j,k)
    end function getex
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function gettx(grid,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: gettx(3)
      integer, intent(in) :: i,j,k
      
      gettx = grid%tx(:,i,j,k)
    end function gettx
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getngrd(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getngrd
      
      getngrd = grid%ngrd
    end function getngrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getgrd(grid,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getgrd(grid%ngrd)
      integer, intent(in) :: i,j,k
      
      getgrd = grid%grd(:,i,j,k)
    end function getgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getbcname(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      character(32) :: getbcname
      
      getbcname = grid%bcinfo(i)%bcname
    end function getbcname
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getfamname(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      character(32) :: getfamname
      
      getfamname = grid%bcinfo(i)%famname
    end function getfamname
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getbcistart(grid,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i,j
      integer :: getbcistart
      
      getbcistart = grid%bcinfo(i)%istart(j)
    end function getbcistart
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getbciend(grid,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i,j
      integer :: getbciend
      
      getbciend = grid%bcinfo(i)%iend(j)
    end function getbciend
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getconnectdonor(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getconnectdonor
      integer, intent(in) :: i
      
      getconnectdonor = grid%connectinfo(i)%donor
    end function getconnectdonor
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnecttransmat(grid,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getconnecttransmat
      integer, intent(in) :: i,j,k
      
      getconnecttransmat = grid%connectinfo(i)%transmat(j,k)
    end function getconnecttransmat
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectistart(grid,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i,j
      integer :: getconnectistart
      
      getconnectistart = grid%connectinfo(i)%istart(j)
    end function getconnectistart
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectiend(grid,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i,j
      integer :: getconnectiend
      
      getconnectiend = grid%connectinfo(i)%iend(j)
    end function getconnectiend
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectistart_donor(grid,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i,j
      integer :: getconnectistart_donor
      
      getconnectistart_donor = grid%connectinfo(i)%istart_donor(j)
    end function getconnectistart_donor
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectiend_donor(grid,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i,j
      integer :: getconnectiend_donor
      
      getconnectiend_donor = grid%connectinfo(i)%iend_donor(j)
    end function getconnectiend_donor
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module grid_module
