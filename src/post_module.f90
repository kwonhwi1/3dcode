module post_module
  use config_module
  use grid_module
  use variable_module
  implicit none
  private
  public :: t_post
  
  type t_post
    private
    real(8), dimension(:,:,:,:), allocatable :: pv,tv
    real(8), dimension(:,:,:,:,:), allocatable :: qq
    real(8) :: pref
    integer :: size,rank,nsteady
    integer :: imax,jmax,kmax
    integer :: npv,ntv,nqq
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: export_variable
  end type t_post
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(post,config,grid,variable)
      implicit none
      class(t_post), intent(out) :: post
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      
      post%size = config%getsize()
      post%rank = config%getrank()
      post%nsteady = config%getnsteady()
      post%pref = config%getpref()
      post%imax = grid%getimax()
      post%jmax = grid%getjmax()
      post%kmax = grid%getkmax()
      post%npv = variable%getnpv()
      post%ntv = variable%getntv()
      post%nqq = variable%getnqq() 

      allocate(post%pv(post%npv,post%imax,post%jmax,post%kmax))
      allocate(post%tv(post%ntv,post%imax,post%jmax,post%kmax))
      allocate(post%qq(post%nqq,variable%getnpv(),post%imax,post%jmax,post%kmax))
      
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(post)
      implicit none
      class(t_post), intent(inout) :: post
      
      deallocate(post%pv,post%tv,post%qq)
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine export_variable(post,variable,nt_phy,nt)
      implicit none
      include 'mpif.h'
      class(t_post), intent(inout) :: post
      type(t_variable), intent(in) :: variable
      integer, intent(in) :: nt_phy,nt
      integer :: i,j,k,l,n,m,ier,io
      character(7) :: iter_tag
    
      if(post%nsteady.eq.1) then
        write(iter_tag,'(i4.4)') nt_phy
      else
        write(iter_tag,'(i7.7)') nt
      end if
      
      do m=0,post%size-1
        if(m.eq.post%rank) then
          if(m.eq.0) then
            open(newunit=io,file='./out_'//trim(iter_tag)//'.dat',status='unknown',action='write',form='unformatted')
            write(io) post%size,nt_phy,nt,post%nqq
          else
            open(newunit=io,file='./out_'//trim(iter_tag)//'.dat',status='old',action='write',position='append',form='unformatted')          
          end if
          write(io) post%rank,post%imax,post%jmax,post%kmax
          do k=2,post%kmax
            do j=2,post%jmax
              do i=2,post%imax
                post%pv(:,i,j,k) = variable%getpv(i,j,k)
                post%pv(1,i,j,k) = post%pv(1,i,j,k)+post%pref 
                post%tv(:,i,j,k) = variable%gettv(i,j,k)
                do n=1,post%nqq
                  post%qq(1,:,i,j,k) = variable%getqq(1,i,j,k)
                  post%qq(2,:,i,j,k) = variable%getqq(2,i,j,k)
                end do
              end do
            end do
          end do
          write(io) ((((post%pv(n,i,j,k),n=1,post%npv),i=2,post%imax),j=2,post%jmax),k=2,post%kmax)
          write(io) ((((post%tv(n,i,j,k),n=1,post%ntv),i=2,post%imax),j=2,post%jmax),k=2,post%kmax)
          write(io) (((((post%qq(l,n,i,j,k),l=1,post%nqq),n=1,post%npv),i=2,post%imax),j=2,post%jmax),k=2,post%kmax)
          close(io)
        end if
        call mpi_barrier(mpi_comm_world,ier)
      end do
      
    end subroutine export_variable
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module post_module
