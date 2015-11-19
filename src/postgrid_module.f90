module postgrid_module
  use config_module
  implicit none
#include <cgnslib_f.h>
#include <cgnstypes_f.h>
  private
  public :: t_grid

  type t_bcinfo
    integer :: istart(3),iend(3)
    character(32) :: bcname,famname
  end type t_bcinfo
  
  type t_connectinfo
    integer :: donor,transmat(3,3)
    integer :: istart(3),iend(3),istart_donor(3),iend_donor(3)
  end type t_connectinfo

  type t_zone
    integer :: imax,jmax,kmax,nbc,ncon
    character(32) :: zonename
    type(t_bcinfo), dimension(:), allocatable :: bcinfo
    type(t_connectinfo), dimension(:), allocatable ::connectinfo  
    real(8), dimension(:,:,:,:), allocatable :: x
    real(8), dimension(:,:,:,:), allocatable :: cx,ex,tx
    real(8), dimension(:,:,:,:), allocatable :: grd ! vol,xcen,ycen
  end type t_zone

  type t_grid
    private
    integer :: nzone,ngrd
    type(t_zone), dimension(:), allocatable :: zone
    contains
      procedure :: construct
      procedure :: destruct
      procedure, private :: calnormal
      procedure, private :: calvolume
      procedure :: getnzone
      procedure :: getimax
      procedure :: getjmax
      procedure :: getkmax
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
      character(32), dimension(:), allocatable :: famname,fambcname
      integer :: i,j,n,m
      integer :: ifile,ier
      integer :: nfam
      cgsize_t :: nmax(3),nstart(3)
      cgsize_t :: isize(9),ipnts(6),ipntsdonor(6)
      integer :: transvec(3) 
      integer :: nfambc,ngeo,ibctype
      integer :: normallist
      integer, parameter :: dim = 3
      real(8), dimension(:,:,:), allocatable :: xx,yy,zz
      
      grid%ngrd=4
      call cg_open_f('./'//trim(config%getname())//'.cgns',cg_mode_read,ifile,ier)
      if(ier.ne.cg_ok) call cg_error_exit_f  
      call cg_nzones_f(ifile,1,grid%nzone,ier)
      allocate(grid%zone(grid%nzone))

      do n= 1,grid%nzone
        call cg_zone_read_f(ifile,1,n,grid%zone(n)%zonename,isize,ier)
        grid%zone(n)%imax = isize(1)
        grid%zone(n)%jmax = isize(2)
        grid%zone(n)%kmax = isize(3)
      end do

      do n=1,grid%nzone
        allocate(grid%zone(n)%x(dim,2:grid%zone(n)%imax+1,2:grid%zone(n)%jmax+1,2:grid%zone(n)%kmax+1))
        allocate(xx(2:grid%zone(n)%imax+1,2:grid%zone(n)%jmax+1,2:grid%zone(n)%kmax+1) &
                ,yy(2:grid%zone(n)%imax+1,2:grid%zone(n)%jmax+1,2:grid%zone(n)%kmax+1) &
                ,zz(2:grid%zone(n)%imax+1,2:grid%zone(n)%jmax+1,2:grid%zone(n)%kmax+1))
          
        nmax(1) = grid%zone(n)%imax; nmax(2) = grid%zone(n)%jmax; nmax(3) = grid%zone(n)%kmax
        nstart  = 1
          
        call cg_coord_read_f(ifile,1,n,'CoordinateX',realdouble,nstart,nmax,xx,ier)
        call cg_coord_read_f(ifile,1,n,'CoordinateY',realdouble,nstart,nmax,yy,ier)
        call cg_coord_read_f(ifile,1,n,'CoordinateZ',realdouble,nstart,nmax,zz,ier)
          
        grid%zone(n)%x(1,:,:,:) = xx*config%getscale()
        grid%zone(n)%x(2,:,:,:) = yy*config%getscale()
        grid%zone(n)%x(3,:,:,:) = zz*config%getscale()
          
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
          
        call cg_nbocos_f(ifile,1,n,grid%zone(n)%nbc,ier)
        allocate(grid%zone(n)%bcinfo(grid%zone(n)%nbc))
          
        do i=1,grid%zone(n)%nbc
          call cg_goto_f(ifile,1,ier,'Zone_t',n,'ZoneBC_t',1,'BC_t',i,'end')
          call cg_famname_read_f(nodename,ier)
          call cg_boco_read_f(ifile,1,n,i,ipnts,normallist,ier)
            
          do j=1,dim
            if(ipnts(j).eq.ipnts(dim+j)) then
              if(ipnts(j).eq.1) then
                select case(j)
                case(1)
                  grid%zone(n)%bcinfo(i)%istart(1) = ipnts(1)
                  grid%zone(n)%bcinfo(i)%istart(2) = ipnts(2)+1
                  grid%zone(n)%bcinfo(i)%istart(3) = ipnts(3)+1
                  grid%zone(n)%bcinfo(i)%iend(1)   = ipnts(4)
                  grid%zone(n)%bcinfo(i)%iend(2)   = ipnts(5)
                  grid%zone(n)%bcinfo(i)%iend(3)   = ipnts(6)
                case(2)
                  grid%zone(n)%bcinfo(i)%istart(1) = ipnts(1)+1
                  grid%zone(n)%bcinfo(i)%istart(2) = ipnts(2)
                  grid%zone(n)%bcinfo(i)%istart(3) = ipnts(3)+1
                  grid%zone(n)%bcinfo(i)%iend(1)   = ipnts(4)
                  grid%zone(n)%bcinfo(i)%iend(2)   = ipnts(5)
                  grid%zone(n)%bcinfo(i)%iend(3)   = ipnts(6)
                case(3)
                  grid%zone(n)%bcinfo(i)%istart(1) = ipnts(1)+1
                  grid%zone(n)%bcinfo(i)%istart(2) = ipnts(2)+1
                  grid%zone(n)%bcinfo(i)%istart(3) = ipnts(3)
                  grid%zone(n)%bcinfo(i)%iend(1)   = ipnts(4)
                  grid%zone(n)%bcinfo(i)%iend(2)   = ipnts(5)
                  grid%zone(n)%bcinfo(i)%iend(3)   = ipnts(6)
                end select
              else
                select case(j)
                case(1)
                  grid%zone(n)%bcinfo(i)%istart(1) = ipnts(1)+1
                  grid%zone(n)%bcinfo(i)%istart(2) = ipnts(2)+1
                  grid%zone(n)%bcinfo(i)%istart(3) = ipnts(3)+1
                  grid%zone(n)%bcinfo(i)%iend(1)   = ipnts(4)+1
                  grid%zone(n)%bcinfo(i)%iend(2)   = ipnts(5)
                  grid%zone(n)%bcinfo(i)%iend(3)   = ipnts(6)
                case(2)
                  grid%zone(n)%bcinfo(i)%istart(1) = ipnts(1)+1
                  grid%zone(n)%bcinfo(i)%istart(2) = ipnts(2)+1
                  grid%zone(n)%bcinfo(i)%istart(3) = ipnts(3)+1
                  grid%zone(n)%bcinfo(i)%iend(1)   = ipnts(4)
                  grid%zone(n)%bcinfo(i)%iend(2)   = ipnts(5)+1
                  grid%zone(n)%bcinfo(i)%iend(3)   = ipnts(6)
                case(3)
                  grid%zone(n)%bcinfo(i)%istart(1) = ipnts(1)+1
                  grid%zone(n)%bcinfo(i)%istart(2) = ipnts(2)+1
                  grid%zone(n)%bcinfo(i)%istart(3) = ipnts(3)+1
                  grid%zone(n)%bcinfo(i)%iend(1)   = ipnts(4)
                  grid%zone(n)%bcinfo(i)%iend(2)   = ipnts(5)
                  grid%zone(n)%bcinfo(i)%iend(3)   = ipnts(6)+1
                end select
              end if
            end if
          end do
         
          do j=1,nfam
            if(trim(nodename).eq.trim(famname(j))) then
              grid%zone(n)%bcinfo(i)%bcname = fambcname(j)
              if(trim(grid%zone(n)%bcinfo(i)%bcname).eq.'UserDefined') then
                grid%zone(n)%bcinfo(i)%bcname = famname(j)
              end if
              grid%zone(n)%bcinfo(i)%famname = famname(j)
            end if
          end do
        end do

        call cg_n1to1_f(ifile,1,n,grid%zone(n)%ncon,ier)
        allocate(grid%zone(n)%connectinfo(grid%zone(n)%ncon))
          
        do j=1,grid%zone(n)%ncon
          call cg_1to1_read_f(ifile,1,n,j,connectname,donorname,ipnts,ipntsdonor,transvec,ier)
         
          grid%zone(n)%connectinfo(j)%transmat = 0
          grid%zone(n)%connectinfo(j)%transmat(abs(transvec(1)),1)=transvec(1)/abs(transvec(1))
          grid%zone(n)%connectinfo(j)%transmat(abs(transvec(2)),2)=transvec(2)/abs(transvec(2))
          grid%zone(n)%connectinfo(j)%transmat(abs(transvec(3)),3)=transvec(3)/abs(transvec(3))
            
          do m=1,dim
            if(ipnts(m).eq.ipnts(dim+m)) then
              if(ipnts(m).eq.1) then
                select case(m)
                case(1)
                  grid%zone(n)%connectinfo(j)%istart(1) = ipnts(1)+1
                  grid%zone(n)%connectinfo(j)%istart(2) = ipnts(2)
                  grid%zone(n)%connectinfo(j)%istart(3) = ipnts(3)
                  grid%zone(n)%connectinfo(j)%iend(1)   = ipnts(4)+1
                  grid%zone(n)%connectinfo(j)%iend(2)   = ipnts(5)
                  grid%zone(n)%connectinfo(j)%iend(3)   = ipnts(6)
                case(2)
                  grid%zone(n)%connectinfo(j)%istart(1) = ipnts(1)
                  grid%zone(n)%connectinfo(j)%istart(2) = ipnts(2)+1
                  grid%zone(n)%connectinfo(j)%istart(3) = ipnts(3)
                  grid%zone(n)%connectinfo(j)%iend(1)   = ipnts(4)
                  grid%zone(n)%connectinfo(j)%iend(2)   = ipnts(5)+1
                  grid%zone(n)%connectinfo(j)%iend(3)   = ipnts(6)
                case(3)
                  grid%zone(n)%connectinfo(j)%istart(1) = ipnts(1)
                  grid%zone(n)%connectinfo(j)%istart(2) = ipnts(2)
                  grid%zone(n)%connectinfo(j)%istart(3) = ipnts(3)+1
                  grid%zone(n)%connectinfo(j)%iend(1)   = ipnts(4)
                  grid%zone(n)%connectinfo(j)%iend(2)   = ipnts(5)
                  grid%zone(n)%connectinfo(j)%iend(3)   = ipnts(6)+1
                end select
              else
                grid%zone(n)%connectinfo(j)%istart(1) = ipnts(1)
                grid%zone(n)%connectinfo(j)%istart(2) = ipnts(2)
                grid%zone(n)%connectinfo(j)%istart(3) = ipnts(3)
                grid%zone(n)%connectinfo(j)%iend(1)   = ipnts(4)
                grid%zone(n)%connectinfo(j)%iend(2)   = ipnts(5)
                grid%zone(n)%connectinfo(j)%iend(3)   = ipnts(6)
              end if
            end if
              
            if(ipntsdonor(m).eq.ipntsdonor(dim+m)) then
              if(ipntsdonor(m).eq.1) then
                grid%zone(n)%connectinfo(j)%istart_donor(1) = ipntsdonor(1)
                grid%zone(n)%connectinfo(j)%istart_donor(2) = ipntsdonor(2)
                grid%zone(n)%connectinfo(j)%istart_donor(3) = ipntsdonor(3)
                grid%zone(n)%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)
                grid%zone(n)%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)
                grid%zone(n)%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)
              else
                select case(m)
                case(1)
                  grid%zone(n)%connectinfo(j)%istart_donor(1) = ipntsdonor(1)+1
                  grid%zone(n)%connectinfo(j)%istart_donor(2) = ipntsdonor(2)
                  grid%zone(n)%connectinfo(j)%istart_donor(3) = ipntsdonor(3)
                  grid%zone(n)%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)+1
                  grid%zone(n)%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)
                  grid%zone(n)%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)
                case(2)
                  grid%zone(n)%connectinfo(j)%istart_donor(1) = ipntsdonor(1)
                  grid%zone(n)%connectinfo(j)%istart_donor(2) = ipntsdonor(2)+1
                  grid%zone(n)%connectinfo(j)%istart_donor(3) = ipntsdonor(3)
                  grid%zone(n)%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)
                  grid%zone(n)%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)+1
                  grid%zone(n)%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)
                case(3)
                  grid%zone(n)%connectinfo(j)%istart_donor(1) = ipntsdonor(1)
                  grid%zone(n)%connectinfo(j)%istart_donor(2) = ipntsdonor(2)
                  grid%zone(n)%connectinfo(j)%istart_donor(3) = ipntsdonor(3)+1
                  grid%zone(n)%connectinfo(j)%iend_donor(1)   = ipntsdonor(4)
                  grid%zone(n)%connectinfo(j)%iend_donor(2)   = ipntsdonor(5)
                  grid%zone(n)%connectinfo(j)%iend_donor(3)   = ipntsdonor(6)+1
                end select
              end if
            end if
          end do

          do i=1,grid%nzone
            if(trim(donorname).eq.trim(grid%zone(n)%zonename)) then
              grid%zone(n)%connectinfo(j)%donor = i-1
            end if
          end do  
        end do      
        if(allocated(famname)) deallocate(famname)
        if(allocated(fambcname)) deallocate(fambcname)
        if(allocated(xx)) deallocate(xx)
        if(allocated(yy)) deallocate(yy)
        if(allocated(zz)) deallocate(zz)
        allocate(grid%zone(n)%grd(grid%ngrd,grid%zone(n)%imax+1,grid%zone(n)%jmax+1,grid%zone(n)%kmax+1))
      end do
      call cg_close_f(ifile,ier)
      
      call grid%calnormal() !! must be called before calvolume !!
      call grid%calvolume() 
      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(grid)
      implicit none
      class(t_grid), intent(inout) :: grid
      integer :: n

      do n = 1,grid%nzone
        if(allocated(grid%zone(n)%x))           deallocate(grid%zone(n)%x)
        if(allocated(grid%zone(n)%cx))          deallocate(grid%zone(n)%cx)
        if(allocated(grid%zone(n)%ex))          deallocate(grid%zone(n)%ex)
        if(allocated(grid%zone(n)%tx))          deallocate(grid%zone(n)%tx)
        if(allocated(grid%zone(n)%grd))         deallocate(grid%zone(n)%grd)
        if(allocated(grid%zone(n)%bcinfo))      deallocate(grid%zone(n)%bcinfo)
        if(allocated(grid%zone(n)%connectinfo)) deallocate(grid%zone(n)%connectinfo)
      end do
     
      if(allocated(grid%zone)) deallocate(grid%zone)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calnormal(grid)
      implicit none
      class(t_grid), intent(inout) :: grid
      integer :: i,j,k,n,m,l
      real(8) :: xc,yc,zc,xe,ye,ze,xs,ys,zs 
      integer, parameter :: dim = 3
      
      do n=1,grid%nzone
        allocate(grid%zone(n)%cx(3,1:grid%zone(n)%imax,1:grid%zone(n)%jmax+1,1:grid%zone(n)%kmax+1))
        allocate(grid%zone(n)%ex(3,1:grid%zone(n)%imax+1,1:grid%zone(n)%jmax,1:grid%zone(n)%kmax+1))
        allocate(grid%zone(n)%tx(3,1:grid%zone(n)%imax+1,1:grid%zone(n)%jmax+1,1:grid%zone(n)%kmax))
      
        do k=2,grid%zone(n)%kmax
          do j=2,grid%zone(n)%jmax
            do i=1,grid%zone(n)%imax
              xe = grid%zone(n)%x(1,i+1,j+1,k)   - grid%zone(n)%x(1,i+1,j,k+1)
              ye = grid%zone(n)%x(2,i+1,j+1,k)   - grid%zone(n)%x(2,i+1,j,k+1)
              ze = grid%zone(n)%x(3,i+1,j+1,k)   - grid%zone(n)%x(3,i+1,j,k+1)
              xs = grid%zone(n)%x(1,i+1,j+1,k+1) - grid%zone(n)%x(1,i+1,j,k)
              ys = grid%zone(n)%x(2,i+1,j+1,k+1) - grid%zone(n)%x(2,i+1,j,k)
              zs = grid%zone(n)%x(3,i+1,j+1,k+1) - grid%zone(n)%x(3,i+1,j,k)
              grid%zone(n)%cx(1,i,j,k) = 0.5d0*(ye*zs - ys*ze)
              grid%zone(n)%cx(2,i,j,k) = 0.5d0*(xs*ze - xe*zs)
              grid%zone(n)%cx(3,i,j,k) = 0.5d0*(xe*ys - xs*ye)
            end do
          end do
        end do

        do i=2,grid%zone(n)%imax
          do k=2,grid%zone(n)%kmax
            do j=1,grid%zone(n)%jmax
              xc = grid%zone(n)%x(1,i,j+1,k)   - grid%zone(n)%x(1,i+1,j+1,k+1)
              yc = grid%zone(n)%x(2,i,j+1,k)   - grid%zone(n)%x(2,i+1,j+1,k+1)
              zc = grid%zone(n)%x(3,i,j+1,k)   - grid%zone(n)%x(3,i+1,j+1,k+1)    
              xs = grid%zone(n)%x(1,i,j+1,k+1) - grid%zone(n)%x(1,i+1,j+1,k)
              ys = grid%zone(n)%x(2,i,j+1,k+1) - grid%zone(n)%x(2,i+1,j+1,k)
              zs = grid%zone(n)%x(3,i,j+1,k+1) - grid%zone(n)%x(3,i+1,j+1,k)
              grid%zone(n)%ex(1,i,j,k) = 0.5d0*(yc*zs - ys*zc)
              grid%zone(n)%ex(2,i,j,k) = 0.5d0*(xs*zc - xc*zs)
              grid%zone(n)%ex(3,i,j,k) = 0.5d0*(xc*ys - xs*yc)   
            end do
          end do
        end do

        do j=2,grid%zone(n)%jmax
          do i=2,grid%zone(n)%imax
            do k=1,grid%zone(n)%kmax
              xc = grid%zone(n)%x(1,i+1,j,k+1)   - grid%zone(n)%x(1,i,j+1,k+1)
              yc = grid%zone(n)%x(2,i+1,j,k+1)   - grid%zone(n)%x(2,i,j+1,k+1)
              zc = grid%zone(n)%x(3,i+1,j,k+1)   - grid%zone(n)%x(3,i,j+1,k+1)
              xe = grid%zone(n)%x(1,i+1,j+1,k+1) - grid%zone(n)%x(1,i,j,k+1)
              ye = grid%zone(n)%x(2,i+1,j+1,k+1) - grid%zone(n)%x(2,i,j,k+1)
              ze = grid%zone(n)%x(3,i+1,j+1,k+1) - grid%zone(n)%x(3,i,j,k+1)
              grid%zone(n)%tx(1,i,j,k) = 0.5d0*(yc*ze - ye*zc)
              grid%zone(n)%tx(2,i,j,k) = 0.5d0*(xe*zc - xc*ze)
              grid%zone(n)%tx(3,i,j,k) = 0.5d0*(xc*ye - xe*yc) 
            end do
          end do
        end do
      
        do m = 1,grid%zone(n)%nbc
          if((trim(grid%zone(n)%bcinfo(m)%bcname).eq.'BCDegeneratePoint').or.(trim(grid%zone(n)%bcinfo(m)%bcname).eq.'BCDegenerateLine')) then
            if(grid%zone(n)%bcinfo(m)%istart(1).eq.grid%zone(n)%bcinfo(m)%iend(1)) then
              if(grid%zone(n)%bcinfo(m)%istart(1).eq.1) then
                grid%zone(n)%cx(1,1,:,:) = 0.d0
                grid%zone(n)%cx(2,1,:,:) = 0.d0
                grid%zone(n)%cx(3,1,:,:) = 0.d0
              else
                grid%zone(n)%cx(1,grid%zone(n)%imax,:,:) = 0.d0
                grid%zone(n)%cx(2,grid%zone(n)%imax,:,:) = 0.d0
                grid%zone(n)%cx(3,grid%zone(n)%imax,:,:) = 0.d0
              end if
            end if
            if(grid%zone(n)%bcinfo(m)%istart(2).eq.grid%zone(n)%bcinfo(m)%iend(2)) then
              if(grid%zone(n)%bcinfo(m)%istart(2).eq.1) then
                grid%zone(n)%ex(1,:,1,:) = 0.d0
                grid%zone(n)%ex(2,:,1,:) = 0.d0
                grid%zone(n)%ex(3,:,1,:) = 0.d0
              else
                grid%zone(n)%ex(1,:,grid%zone(n)%jmax,:) = 0.d0
                grid%zone(n)%ex(2,:,grid%zone(n)%jmax,:) = 0.d0
                grid%zone(n)%ex(3,:,grid%zone(n)%jmax,:) = 0.d0
              end if
            end if
            if(grid%zone(n)%bcinfo(m)%istart(3).eq.grid%zone(n)%bcinfo(m)%iend(3)) then
              if(grid%zone(n)%bcinfo(m)%istart(3).eq.1) then
                grid%zone(n)%tx(1,:,:,1) = 0.d0
                grid%zone(n)%tx(2,:,:,1) = 0.d0
                grid%zone(n)%tx(3,:,:,1) = 0.d0
              else
                grid%zone(n)%tx(1,:,:,grid%zone(n)%kmax) = 0.d0
                grid%zone(n)%tx(2,:,:,grid%zone(n)%kmax) = 0.d0
                grid%zone(n)%tx(3,:,:,grid%zone(n)%kmax) = 0.d0
              end if
            end if 
          end if
        end do
      
        grid%zone(n)%cx(:,:,1,:) = grid%zone(n)%cx(:,:,2,:)
        grid%zone(n)%cx(:,:,grid%zone(n)%jmax+1,:) = grid%zone(n)%cx(:,:,grid%zone(n)%jmax,:)
        grid%zone(n)%cx(:,:,:,1) = grid%zone(n)%cx(:,:,:,2)
        grid%zone(n)%cx(:,:,:,grid%zone(n)%kmax+1) = grid%zone(n)%cx(:,:,:,grid%zone(n)%kmax)
      
        grid%zone(n)%ex(:,1,:,:) = grid%zone(n)%ex(:,2,:,:)
        grid%zone(n)%ex(:,grid%zone(n)%imax+1,:,:) = grid%zone(n)%ex(:,grid%zone(n)%imax,:,:)
        grid%zone(n)%ex(:,:,:,1) = grid%zone(n)%ex(:,:,:,2)
        grid%zone(n)%ex(:,:,:,grid%zone(n)%kmax+1) = grid%zone(n)%ex(:,:,:,grid%zone(n)%kmax)
      
        grid%zone(n)%tx(:,1,:,:) = grid%zone(n)%tx(:,2,:,:)
        grid%zone(n)%tx(:,grid%zone(n)%imax+1,:,:) = grid%zone(n)%tx(:,grid%zone(n)%imax,:,:)
        grid%zone(n)%tx(:,:,1,:) = grid%zone(n)%tx(:,:,2,:)
        grid%zone(n)%tx(:,:,grid%zone(n)%jmax+1,:) = grid%zone(n)%tx(:,:,grid%zone(n)%jmax,:)

      ! reodering
        do l=1,grid%zone(n)%ncon
          do m=1,dim
            if(grid%zone(n)%connectinfo(l)%istart(m).ne.grid%zone(n)%connectinfo(l)%iend(m) ) then
              grid%zone(n)%connectinfo(l)%istart(m) = grid%zone(n)%connectinfo(l)%istart(m) + 1            
            end if
  
            if(grid%zone(n)%connectinfo(l)%istart_donor(m).ne.grid%zone(n)%connectinfo(l)%iend_donor(m) ) then
              if(grid%zone(n)%connectinfo(l)%istart_donor(m).lt.grid%zone(n)%connectinfo(l)%iend_donor(m) ) then
                grid%zone(n)%connectinfo(l)%istart_donor(m) = grid%zone(n)%connectinfo(l)%istart_donor(m) + 1
              else if(grid%zone(n)%connectinfo(l)%istart_donor(m).gt.grid%zone(n)%connectinfo(l)%iend_donor(m) ) then
                grid%zone(n)%connectinfo(l)%iend_donor(m) = grid%zone(n)%connectinfo(l)%iend_donor(m) + 1            
              end if
            end if
          end do
        end do
      end do
    end subroutine calnormal
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calvolume(grid)
      implicit none
      include 'mpif.h'
      class(t_grid), intent(inout) :: grid
      integer :: i,j,k,n
      real(8) :: pa(3),pb(3),pc(3),pd(3)
      real(8) :: volp1,volp2,volp3,volp4,volp5,volp6
      
      do n=1,grid%nzone
        do k=2,grid%zone(n)%kmax
          do j=2,grid%zone(n)%jmax
            do i=2,grid%zone(n)%imax
              grid%zone(n)%grd(2,i,j,k) = 0.125d0*(grid%zone(n)%x(1,i,j,k)+grid%zone(n)%x(1,i+1,j,k)+grid%zone(n)%x(1,i,j+1,k)+grid%zone(n)%x(1,i,j,k+1) &
                                + grid%zone(n)%x(1,i+1,j+1,k)+grid%zone(n)%x(1,i,j+1,k+1)+grid%zone(n)%x(1,i+1,j,k+1)+grid%zone(n)%x(1,i+1,j+1,k+1))
              grid%zone(n)%grd(3,i,j,k) = 0.125d0*(grid%zone(n)%x(2,i,j,k)+grid%zone(n)%x(2,i+1,j,k)+grid%zone(n)%x(2,i,j+1,k)+grid%zone(n)%x(2,i,j,k+1) &
                                + grid%zone(n)%x(2,i+1,j+1,k)+grid%zone(n)%x(2,i,j+1,k+1)+grid%zone(n)%x(2,i+1,j,k+1)+grid%zone(n)%x(2,i+1,j+1,k+1))
              grid%zone(n)%grd(4,i,j,k) = 0.125d0*(grid%zone(n)%x(3,i,j,k)+grid%zone(n)%x(3,i+1,j,k)+grid%zone(n)%x(3,i,j+1,k)+grid%zone(n)%x(3,i,j,k+1) &
                                + grid%zone(n)%x(3,i+1,j+1,k)+grid%zone(n)%x(3,i,j+1,k+1)+grid%zone(n)%x(3,i+1,j,k+1)+grid%zone(n)%x(3,i+1,j+1,k+1))
                              
              pa(1) = grid%zone(n)%x(1,i,j,k) 
              pb(1) = grid%zone(n)%x(1,i+1,j,k)
              pc(1) = grid%zone(n)%x(1,i+1,j+1,k)
              pd(1) = grid%zone(n)%x(1,i,j+1,k)

              pa(2) = grid%zone(n)%x(2,i,j,k)
              pb(2) = grid%zone(n)%x(2,i+1,j,k)
              pc(2) = grid%zone(n)%x(2,i+1,j+1,k)
              pd(2) = grid%zone(n)%x(2,i,j+1,k)

              pa(3) = grid%zone(n)%x(3,i,j,k)
              pb(3) = grid%zone(n)%x(3,i+1,j,k)
              pc(3) = grid%zone(n)%x(3,i+1,j+1,k)
              pd(3) = grid%zone(n)%x(3,i,j+1,k)
            
              volp1 = pyramid(pa,pb,pc,pd,grid%zone(n)%grd(2,i,j,k),grid%zone(n)%grd(3,i,j,k),grid%zone(n)%grd(4,i,j,k))
            
              pa(1) = grid%zone(n)%x(1,i,j,k+1) 
              pb(1) = grid%zone(n)%x(1,i+1,j,k+1)
              pc(1) = grid%zone(n)%x(1,i+1,j+1,k+1)
              pd(1) = grid%zone(n)%x(1,i,j+1,k+1)

              pa(2) = grid%zone(n)%x(2,i,j,k+1)
              pb(2) = grid%zone(n)%x(2,i+1,j,k+1)
              pc(2) = grid%zone(n)%x(2,i+1,j+1,k+1)
              pd(2) = grid%zone(n)%x(2,i,j+1,k+1)

              pa(3) = grid%zone(n)%x(3,i,j,k+1)
              pb(3) = grid%zone(n)%x(3,i+1,j,k+1)
              pc(3) = grid%zone(n)%x(3,i+1,j+1,k+1)
              pd(3) = grid%zone(n)%x(3,i,j+1,k+1)
            
              volp2 = pyramid(pa,pb,pc,pd,grid%zone(n)%grd(2,i,j,k),grid%zone(n)%grd(3,i,j,k),grid%zone(n)%grd(4,i,j,k))
            
              pa(1) = grid%zone(n)%x(1,i,j,k) 
              pb(1) = grid%zone(n)%x(1,i+1,j,k)
              pc(1) = grid%zone(n)%x(1,i+1,j,k+1)
              pd(1) = grid%zone(n)%x(1,i,j,k+1)

              pa(2) = grid%zone(n)%x(2,i,j,k)
              pb(2) = grid%zone(n)%x(2,i+1,j,k)
              pc(2) = grid%zone(n)%x(2,i+1,j,k+1)
              pd(2) = grid%zone(n)%x(2,i,j,k+1)

              pa(3) = grid%zone(n)%x(3,i,j,k)
              pb(3) = grid%zone(n)%x(3,i+1,j,k)
              pc(3) = grid%zone(n)%x(3,i+1,j,k+1)
              pd(3) = grid%zone(n)%x(3,i,j,k+1)
            
              volp3 = pyramid(pa,pb,pc,pd,grid%zone(n)%grd(2,i,j,k),grid%zone(n)%grd(3,i,j,k),grid%zone(n)%grd(4,i,j,k))
            
              pa(1) = grid%zone(n)%x(1,i,j+1,k) 
              pb(1) = grid%zone(n)%x(1,i+1,j+1,k)
              pc(1) = grid%zone(n)%x(1,i+1,j+1,k+1)
              pd(1) = grid%zone(n)%x(1,i,j+1,k+1)

              pa(2) = grid%zone(n)%x(2,i,j+1,k)
              pb(2) = grid%zone(n)%x(2,i+1,j+1,k)
              pc(2) = grid%zone(n)%x(2,i+1,j+1,k+1)
              pd(2) = grid%zone(n)%x(2,i,j+1,k+1)

              pa(3) = grid%zone(n)%x(3,i,j+1,k)
              pb(3) = grid%zone(n)%x(3,i+1,j+1,k)
              pc(3) = grid%zone(n)%x(3,i+1,j+1,k+1)
              pd(3) = grid%zone(n)%x(3,i,j+1,k+1)
            
              volp4 = pyramid(pa,pb,pc,pd,grid%zone(n)%grd(2,i,j,k),grid%zone(n)%grd(3,i,j,k),grid%zone(n)%grd(4,i,j,k))
            
              pa(1) = grid%zone(n)%x(1,i,j,k) 
              pb(1) = grid%zone(n)%x(1,i,j+1,k)
              pc(1) = grid%zone(n)%x(1,i,j+1,k+1)
              pd(1) = grid%zone(n)%x(1,i,j,k+1)

              pa(2) = grid%zone(n)%x(2,i,j,k)
              pb(2) = grid%zone(n)%x(2,i,j+1,k)
              pc(2) = grid%zone(n)%x(2,i,j+1,k+1)
              pd(2) = grid%zone(n)%x(2,i,j,k+1)

              pa(3) = grid%zone(n)%x(3,i,j,k)
              pb(3) = grid%zone(n)%x(3,i,j+1,k)
              pc(3) = grid%zone(n)%x(3,i,j+1,k+1)
              pd(3) = grid%zone(n)%x(3,i,j,k+1)
            
              volp5 = pyramid(pa,pb,pc,pd,grid%zone(n)%grd(2,i,j,k),grid%zone(n)%grd(3,i,j,k),grid%zone(n)%grd(4,i,j,k))
            
              pa(1) = grid%zone(n)%x(1,i+1,j,k) 
              pb(1) = grid%zone(n)%x(1,i+1,j+1,k)
              pc(1) = grid%zone(n)%x(1,i+1,j+1,k+1)
              pd(1) = grid%zone(n)%x(1,i+1,j,k+1)

              pa(2) = grid%zone(n)%x(2,i+1,j,k)
              pb(2) = grid%zone(n)%x(2,i+1,j+1,k)
              pc(2) = grid%zone(n)%x(2,i+1,j+1,k+1)
              pd(2) = grid%zone(n)%x(2,i+1,j,k+1)

              pa(3) = grid%zone(n)%x(3,i+1,j,k)
              pb(3) = grid%zone(n)%x(3,i+1,j+1,k)
              pc(3) = grid%zone(n)%x(3,i+1,j+1,k+1)
              pd(3) = grid%zone(n)%x(3,i+1,j,k+1)
            
              volp6 = pyramid(pa,pb,pc,pd,grid%zone(n)%grd(2,i,j,k),grid%zone(n)%grd(3,i,j,k),grid%zone(n)%grd(4,i,j,k))
            
              grid%zone(n)%grd(1,i,j,k) = (volp1+volp2+volp3+volp4+volp5+volp6)/6.d0
              if(grid%zone(n)%grd(1,i,j,k).lt.0.d0) write(*,*) 'negative volian',i,j,k,grid%zone(n)%grd(1,i,j,k)

            end do
          end do
        end do

        do k=2,grid%zone(n)%kmax
          do i=2,grid%zone(n)%imax
            grid%zone(n)%grd(1,i,1,k) = grid%zone(n)%grd(1,i,2,k)
            grid%zone(n)%grd(2,i,1,k) = 2.d0*grid%zone(n)%grd(2,i,2,k) - grid%zone(n)%grd(2,i,3,k)
            grid%zone(n)%grd(3,i,1,k) = 2.d0*grid%zone(n)%grd(3,i,2,k) - grid%zone(n)%grd(3,i,3,k)
            grid%zone(n)%grd(4,i,1,k) = 2.d0*grid%zone(n)%grd(4,i,2,k) - grid%zone(n)%grd(4,i,3,k)
          
            grid%zone(n)%grd(1,i,grid%zone(n)%jmax+1,k) = grid%zone(n)%grd(1,i,grid%zone(n)%jmax,k)
            grid%zone(n)%grd(2,i,grid%zone(n)%jmax+1,k) = 2.d0*grid%zone(n)%grd(2,i,grid%zone(n)%jmax,k) - grid%zone(n)%grd(2,i,grid%zone(n)%jmax-1,k)
            grid%zone(n)%grd(3,i,grid%zone(n)%jmax+1,k) = 2.d0*grid%zone(n)%grd(3,i,grid%zone(n)%jmax,k) - grid%zone(n)%grd(3,i,grid%zone(n)%jmax-1,k)
            grid%zone(n)%grd(4,i,grid%zone(n)%jmax+1,k) = 2.d0*grid%zone(n)%grd(4,i,grid%zone(n)%jmax,k) - grid%zone(n)%grd(4,i,grid%zone(n)%jmax-1,k)
          end do
        end do
        
        do k=2,grid%zone(n)%kmax
          do j=1,grid%zone(n)%jmax+1
            grid%zone(n)%grd(1,1,j,k) = grid%zone(n)%grd(1,2,j,k)
            grid%zone(n)%grd(2,1,j,k) = 2.d0*grid%zone(n)%grd(2,2,j,k) - grid%zone(n)%grd(2,3,j,k)
            grid%zone(n)%grd(3,1,j,k) = 2.d0*grid%zone(n)%grd(3,2,j,k) - grid%zone(n)%grd(3,3,j,k)
            grid%zone(n)%grd(4,1,j,k) = 2.d0*grid%zone(n)%grd(4,2,j,k) - grid%zone(n)%grd(4,3,j,k)
          
            grid%zone(n)%grd(1,grid%zone(n)%imax+1,j,k) = grid%zone(n)%grd(1,grid%zone(n)%imax,j,k)
            grid%zone(n)%grd(2,grid%zone(n)%imax+1,j,k) = 2.d0*grid%zone(n)%grd(2,grid%zone(n)%imax,j,k) - grid%zone(n)%grd(2,grid%zone(n)%imax-1,j,k)
            grid%zone(n)%grd(3,grid%zone(n)%imax+1,j,k) = 2.d0*grid%zone(n)%grd(3,grid%zone(n)%imax,j,k) - grid%zone(n)%grd(3,grid%zone(n)%imax-1,j,k)
            grid%zone(n)%grd(4,grid%zone(n)%imax+1,j,k) = 2.d0*grid%zone(n)%grd(4,grid%zone(n)%imax,j,k) - grid%zone(n)%grd(4,grid%zone(n)%imax-1,j,k)
          end do
        end do

        do j=1,grid%zone(n)%jmax+1
          do i=1,grid%zone(n)%imax+1
            grid%zone(n)%grd(1,i,j,1) = grid%zone(n)%grd(1,i,j,2)
            grid%zone(n)%grd(2,i,j,1) = 2.d0*grid%zone(n)%grd(2,i,j,2) - grid%zone(n)%grd(2,3,j,3)
            grid%zone(n)%grd(3,i,j,1) = 2.d0*grid%zone(n)%grd(3,i,j,2) - grid%zone(n)%grd(3,3,j,3)
            grid%zone(n)%grd(4,i,j,1) = 2.d0*grid%zone(n)%grd(4,i,j,2) - grid%zone(n)%grd(4,3,j,3)
          
            grid%zone(n)%grd(1,i,j,grid%zone(n)%kmax+1) = grid%zone(n)%grd(1,i,j,grid%zone(n)%kmax)
            grid%zone(n)%grd(2,i,j,grid%zone(n)%kmax+1) = 2.d0*grid%zone(n)%grd(2,i,j,grid%zone(n)%kmax) - grid%zone(n)%grd(2,i,j,grid%zone(n)%kmax-1)
            grid%zone(n)%grd(3,i,j,grid%zone(n)%kmax+1) = 2.d0*grid%zone(n)%grd(3,i,j,grid%zone(n)%kmax) - grid%zone(n)%grd(3,i,j,grid%zone(n)%kmax-1)
            grid%zone(n)%grd(4,i,j,grid%zone(n)%kmax+1) = 2.d0*grid%zone(n)%grd(4,i,j,grid%zone(n)%kmax) - grid%zone(n)%grd(4,i,j,grid%zone(n)%kmax-1)
          end do
        end do
      end do
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
    pure function getnzone(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getnzone
      
      getnzone = grid%nzone
    end function getnzone
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getimax(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getimax
      
      getimax = grid%zone(i)%imax
    end function getimax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getjmax(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getjmax
      
      getjmax = grid%zone(i)%jmax
    end function getjmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getkmax(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getkmax
      
      getkmax = grid%zone(i)%kmax
    end function getkmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getnbc(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getnbc
      
      getnbc = grid%zone(i)%nbc
    end function getnbc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
    pure function getncon(grid,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: i
      integer :: getncon
      
      getncon = grid%zone(i)%ncon
    end function getncon
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getx(grid,n,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getx(3)
      integer, intent(in) :: n,i,j,k
      
      getx = grid%zone(n)%x(:,i,j,k)
    end function getx
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getcx(grid,n,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getcx(3)
      integer, intent(in) :: n,i,j,k
      
      getcx = grid%zone(n)%cx(:,i,j,k)
    end function getcx
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getex(grid,n,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getex(3)
      integer, intent(in) :: n,i,j,k
      
      getex = grid%zone(n)%ex(:,i,j,k)
    end function getex
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function gettx(grid,n,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: gettx(3)
      integer, intent(in) :: n,i,j,k
      
      gettx = grid%zone(n)%tx(:,i,j,k)
    end function gettx
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getngrd(grid)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getngrd
      
      getngrd = grid%ngrd
    end function getngrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getgrd(grid,n,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      real(8) :: getgrd(grid%ngrd)
      integer, intent(in) :: n,i,j,k
      
      getgrd = grid%zone(n)%grd(:,i,j,k)
    end function getgrd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getbcname(grid,n,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i
      character(32) :: getbcname
      
      getbcname = grid%zone(n)%bcinfo(i)%bcname
    end function getbcname
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getfamname(grid,n,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i
      character(32) :: getfamname
      
      getfamname = grid%zone(n)%bcinfo(i)%famname
    end function getfamname
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getbcistart(grid,n,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i,j
      integer :: getbcistart
      
      getbcistart = grid%zone(n)%bcinfo(i)%istart(j)
    end function getbcistart
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getbciend(grid,n,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i,j
      integer :: getbciend
      
      getbciend = grid%zone(n)%bcinfo(i)%iend(j)
    end function getbciend
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
    pure function getconnectdonor(grid,n,i)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getconnectdonor
      integer, intent(in) :: n,i
      
      getconnectdonor = grid%zone(n)%connectinfo(i)%donor
    end function getconnectdonor
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnecttransmat(grid,n,i,j,k)
      implicit none
      class(t_grid), intent(in) :: grid
      integer :: getconnecttransmat
      integer, intent(in) :: n,i,j,k
      
      getconnecttransmat = grid%zone(n)%connectinfo(i)%transmat(j,k)
    end function getconnecttransmat
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectistart(grid,n,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i,j
      integer :: getconnectistart
      
      getconnectistart = grid%zone(n)%connectinfo(i)%istart(j)
    end function getconnectistart
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectiend(grid,n,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i,j
      integer :: getconnectiend
      
      getconnectiend = grid%zone(n)%connectinfo(i)%iend(j)
    end function getconnectiend
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectistart_donor(grid,n,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i,j
      integer :: getconnectistart_donor
      
      getconnectistart_donor = grid%zone(n)%connectinfo(i)%istart_donor(j)
    end function getconnectistart_donor
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getconnectiend_donor(grid,n,i,j)
      implicit none
      class(t_grid), intent(in) :: grid
      integer, intent(in) :: n,i,j
      integer :: getconnectiend_donor
      
      getconnectiend_donor = grid%zone(n)%connectinfo(i)%iend_donor(j)
    end function getconnectiend_donor
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module postgrid_module
