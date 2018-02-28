module variable_module
  use mpi
  use config_module
  use grid_module
  implicit none
  private
  public :: t_variable

  type t_variable
    private
    integer :: npv,ndv,ntv,nqq
    integer :: size,rank,nsteady,ngrd,imax,jmax,kmax
    integer :: intsize,realsize
    logical :: l_csf
    integer(kind=mpi_offset_kind) :: disp
    real(8), dimension(:,:,:,:), allocatable :: pv      ! p,u,v,w,t,y1,y2,k,o
    real(8), dimension(:,:,:,:), allocatable :: dv      ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
    real(8), dimension(:,:,:,:), allocatable :: tv      ! vis,cond,emut
    real(8), dimension(:,:,:,:,:), allocatable :: qq    ! unsteady n-1,n-2 conservative variable
    real(8), dimension(:,:,:,:), allocatable :: csf,vfg ! surface tension force, volume fraction gradient
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: getpv
      procedure :: getdv
      procedure :: gettv
      procedure :: getqq
      procedure :: getvfg
      procedure :: setpv
      procedure :: setdv
      procedure :: settv
      procedure :: setqq
      procedure :: setcsf
      procedure :: setvfg
      procedure :: export_variable
      procedure :: calvfg
  end type t_variable

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(variable,config,grid)
      implicit none
      class(t_variable), intent(out) :: variable
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      integer :: ier,i

      call mpi_type_size(mpi_integer,variable%intsize,ier)
      call mpi_type_size(mpi_real8,variable%realsize,ier)

      variable%npv = config%getnpv()
      variable%ndv = config%getndv()
      variable%ntv = config%getntv()
      variable%nqq = config%getnqq()
      variable%size = config%getsize()
      variable%rank = config%getrank()
      variable%nsteady = config%getnsteady()
      variable%ngrd = grid%getngrd()
      variable%imax = grid%getimax()
      variable%jmax = grid%getjmax()
      variable%kmax = grid%getkmax()

      variable%l_csf = .false.; if(config%getcsf().ne.0) variable%l_csf = .true.

      allocate(variable%pv(variable%npv,-1:variable%imax+3,-1:variable%jmax+3,-1:variable%kmax+3))
      allocate(variable%dv(variable%ndv,-1:variable%imax+3,-1:variable%jmax+3,-1:variable%kmax+3))
      allocate(variable%tv(variable%ntv,-1:variable%imax+3,-1:variable%jmax+3,-1:variable%kmax+3))
      allocate(variable%qq(variable%npv,variable%nqq,2:variable%imax,2:variable%jmax,2:variable%kmax))
      if(variable%l_csf) then
        allocate(variable%csf(3,2:variable%imax,2:variable%jmax,2:variable%kmax))
        allocate(variable%vfg(3,1:variable%imax+1,1:variable%jmax+1,1:variable%kmax+1))
      end if

      if(variable%rank.eq.0) then
        variable%disp = 0
      else
        variable%disp = 0
        do i=0,variable%rank-1
          variable%disp = variable%disp + variable%intsize*2 + variable%realsize &
                        + variable%realsize*variable%npv*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5) &
                        + variable%realsize*variable%ntv*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5) &
                        + variable%realsize*variable%nqq*variable%npv*(grid%getimax_zone(i)-1)*(grid%getjmax_zone(i)-1)*(grid%getkmax_zone(i)-1)
        end do
      end if
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(variable)
      implicit none
      class(t_variable), intent(inout) :: variable

      if(allocated(variable%pv)) deallocate(variable%pv)
      if(allocated(variable%dv)) deallocate(variable%dv)
      if(allocated(variable%tv)) deallocate(variable%tv)
      if(allocated(variable%qq)) deallocate(variable%qq)
      if(allocated(variable%csf)) deallocate(variable%csf)
      if(allocated(variable%vfg)) deallocate(variable%vfg)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine export_variable(variable,nt_phy,nt,time)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: nt_phy,nt
      real(8), intent(in) :: time
      integer :: ier,io,num
      integer(kind=mpi_offset_kind) :: disp
      character(8) :: iter_tag

      if(variable%nsteady.ge.1) then
        write(iter_tag,'(i8.8)') nt_phy
      else
        write(iter_tag,'(i8.8)') nt
      end if

      disp = variable%disp
      call mpi_file_open(mpi_comm_world,"./out_"//trim(iter_tag)//".dat",mpi_mode_wronly+mpi_mode_create,mpi_info_null,io,ier)

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_write_all(io,nt_phy,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + variable%intsize

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_write_all(io,nt,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + variable%intsize

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      call mpi_file_write_all(io,time,1,mpi_real8,mpi_status_ignore,ier)
      disp = disp + variable%realsize

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = variable%npv*(variable%imax+5)*(variable%jmax+5)*(variable%kmax+5)
      call mpi_file_write_all(io,variable%pv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + variable%realsize*num


      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = variable%ntv*(variable%imax+5)*(variable%jmax+5)*(variable%kmax+5)
      call mpi_file_write_all(io,variable%tv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + variable%realsize*num

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = variable%nqq*variable%npv*(variable%imax-1)*(variable%jmax-1)*(variable%kmax-1)
      call mpi_file_write_all(io,variable%qq,num,mpi_real8,mpi_status_ignore,ier)

      call mpi_file_close(io,ier)

    end subroutine export_variable
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getpv(variable,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k
      real(8) :: getpv(variable%npv)

      getpv = variable%pv(:,i,j,k)

    end function getpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getdv(variable,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k
      real(8) :: getdv(variable%ndv)

      getdv = variable%dv(:,i,j,k)

    end function getdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function gettv(variable,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k
      real(8) :: gettv(variable%ntv)

      gettv = variable%tv(:,i,j,k)

    end function gettv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getqq(variable,n,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k,n
      real(8) :: getqq(variable%npv)

      getqq = variable%qq(:,n,i,j,k)

    end function getqq
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getvfg(variable,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: i,j,k
      real(8) :: getvfg(3)

      getvfg = variable%vfg(:,i,j,k)

    end function getvfg
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var

      variable%pv(n,i,j,k) = var

    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setdv(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var

      variable%dv(n,i,j,k) = var

    end subroutine setdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine settv(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var

      variable%tv(n,i,j,k) = var

    end subroutine settv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setqq(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var(variable%npv)

      variable%qq(:,n,i,j,k) = var

    end subroutine setqq
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setcsf(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var

      variable%csf(n,i,j,k) = var

    end subroutine setcsf
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setvfg(variable,n,i,j,k,var)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer, intent(in) :: n,i,j,k
      real(8), intent(in) :: var

      variable%vfg(n,i,j,k) = var

    end subroutine setvfg
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calvfg(variable,grid)
      implicit none
      class(t_variable), intent(inout) :: variable
      type(t_grid), intent(in) :: grid
      integer :: i,j,k,ll,ii,jj,kk
      integer :: m(2)
      real(8) :: vf(0:variable%imax+2,0:variable%jmax+2,0:variable%kmax+2)
      real(8) :: x_vf(7),nx(3,6),grd(variable%ngrd)

      do k=0,variable%kmax+2
        do j=0,variable%jmax+2
          do i=0,variable%imax+2
            vf(i,j,k) = variable%dv(1,i,j,k)*(variable%pv(6,i,j,k)/variable%dv(5,i,j,k) &
                                             +variable%pv(7,i,j,k)/variable%dv(6,i,j,k))
          end do
        end do
      end do

      m = (/-1,1/)
      do k=1,variable%kmax+1
        do j=1,variable%jmax+1
          do i=1,variable%imax+1
            ll = 0
             do ii = 1,2
               ll = ll+1
               x_vf(ll) = vf(i+m(ii),j,k)
             end do
             do jj = 1,2
               ll = ll+1
               x_vf(ll) = vf(i,j+m(jj),k)
             end do
             do kk = 1,2
               ll = ll+1
               x_vf(ll) = vf(i,j,k+m(kk))
             end do
             x_vf(7) = vf(i,j,k)
             nx(:,1) = -grid%getcx(i-1,j,k)
             nx(:,2) = grid%getcx(i,j,k)
             nx(:,3) = -grid%getex(i,j-1,k)
             nx(:,4) = grid%getex(i,j,k)
             nx(:,5) = -grid%gettx(i,j,k-1)
             nx(:,6) = grid%gettx(i,j,k)
             grd = grid%getgrd(i,j,k)
             variable%vfg(:,i,j,k) = calgrad_gg(variable,x_vf,nx,grd)
           end do
         end do
       end do

    end subroutine calvfg
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function calgrad_gg(variable,x,nx,grd)
      implicit none
      class(t_variable), intent(in) :: variable
      real(8), intent(in) :: x(7),nx(3,6),grd(variable%ngrd)
      real(8) :: calgrad_gg(3)
      real(8) :: xm(6)
      integer :: l

      ! surface tension force only exerted between liquid and its vapor phases
      !                       between different gases is ignored
      ! sigma should be averaged in csfsource


      !! Green-Gauss Gradient
      calgrad_gg = 0.d0

      do l=1,6
        xm(l) = 0.5d0*(x(7)+x(l))
        calgrad_gg(1) = calgrad_gg(1) + xm(l)*nx(1,l)
        calgrad_gg(2) = calgrad_gg(2) + xm(l)*nx(2,l)
        calgrad_gg(3) = calgrad_gg(3) + xm(l)*nx(3,l)
      end do

      calgrad_gg = calgrad_gg/grd(1)
    end function calgrad_gg
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module variable_module
