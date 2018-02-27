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
    integer :: size,rank,nsteady,imax,jmax,kmax
    integer :: intsize,realsize
    integer(kind=mpi_offset_kind) :: disp
    real(8), dimension(:,:,:,:), allocatable :: pv   ! p,u,v,w,t,y1,y2,k,o
    real(8), dimension(:,:,:,:), allocatable :: dv   ! rho,h,rhol,rhov,rhog,snd2,drdp,drdt,drdy1,drdy2,dhdp,dhdt,dhdy1,dhdy2,drdpv,drdtv,drdpl,drdtl
    real(8), dimension(:,:,:,:), allocatable :: tv   ! vis,cond,emut
    real(8), dimension(:,:,:,:,:), allocatable :: qq ! unsteady n-1,n-2 conservative variable
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: getpv
      procedure :: getdv
      procedure :: gettv
      procedure :: getqq
      procedure :: setpv
      procedure :: setdv
      procedure :: settv
      procedure :: setqq
      procedure :: export_variable
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
      variable%imax = grid%getimax()
      variable%jmax = grid%getjmax()
      variable%kmax = grid%getkmax()

      allocate(variable%pv(variable%npv,-1:variable%imax+3,-1:variable%jmax+3,-1:variable%kmax+3))
      allocate(variable%dv(variable%ndv,-1:variable%imax+3,-1:variable%jmax+3,-1:variable%kmax+3))
      allocate(variable%tv(variable%ntv,-1:variable%imax+3,-1:variable%jmax+3,-1:variable%kmax+3))
      allocate(variable%qq(variable%npv,variable%nqq,2:variable%imax,2:variable%jmax,2:variable%kmax))

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
end module variable_module
