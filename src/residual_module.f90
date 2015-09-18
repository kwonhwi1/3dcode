module residual_module
  use mpi
  use config_module
  use grid_module
  use variable_module
  implicit none
  private
  public :: t_resi
  
  type t_resi
    private
    logical :: l_converge
    integer :: npv,imax,jmax,kmax
    integer :: ntmax,npmax,nprint
    integer :: rank,io
    real(8) :: bond
    real(8), dimension(:,:,:,:), allocatable :: qres
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setqres
      procedure :: residual
      procedure :: getconverge
  end type t_resi
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(resi,config,grid,variable)
      implicit none
      class(t_resi), intent(out) :: resi
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      
      resi%rank = config%getrank()
      resi%npmax = config%getnpmax()      
      resi%ntmax = config%getntmax()
      resi%imax = grid%getimax()
      resi%jmax = grid%getjmax()
      resi%kmax = grid%getkmax()
      resi%npv = variable%getnpv()
      resi%nprint = config%getnprint()
      resi%bond = 10.d0**(-config%getbond())
      
      allocate(resi%qres(resi%npv,resi%imax,resi%jmax,resi%kmax))
      
      if(config%getiread().eq.0) then
        if(resi%rank.eq.0) then
          open(newunit=resi%io,file='./res.plt',status='unknown',action='write')
          select case(config%getiturb())
          case(-1,0)
            write(resi%io,*) 'variables="nt","resp","resu","resv","resw","rest","resy1","resy2","resk","reso"' 
          case(-2,-3)
            write(resi%io,*) 'variables="nt","resp","resu","resv","resw","rest","resy1","resy2"' 
          end select
        end if
      else !iread=1 restart
        if(resi%rank.eq.0) then
          open(newunit=resi%io,file='./res.plt',status='old',action='write',position='append')
        end if
      end if
      resi%l_converge = .false.
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(resi)
      implicit none
      class(t_resi), intent(inout) :: resi
      
      if(allocated(resi%qres)) deallocate(resi%qres)
      
      if(resi%rank.eq.0) then
        close(resi%io)
      end if
      
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setqres(resi,variable)
      implicit none
      class(t_resi), intent(inout) :: resi
      type(t_variable), intent(in) :: variable
      integer :: i,j,k
      
      do k=2,resi%kmax
        do j=2,resi%jmax
          do i=2,resi%imax
            resi%qres(:,i,j,k) = variable%getpv(i,j,k)
          end do
        end do
      end do
      
    end subroutine setqres
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine residual(resi,variable,nt_phy,nt)
      implicit none
      class(t_resi), intent(inout) :: resi
      type(t_variable), intent(in) :: variable
      integer, intent(in) :: nt_phy,nt
      integer :: i,j,k,n,ierr,iter
      integer :: iresm,jresm,kresm
      real(8) :: restest,resmax
      real(8), dimension(:) :: res(resi%npv)
      real(8), dimension(:) :: mpi_res(resi%npv),mpi_res_sum(resi%npv)

      if((resi%rank.eq.0).and.(nt.eq.1).and.(nt_phy.eq.1)) then
        write(resi%io,*) 'zone t =" "'
      end if
      
      resmax = 1.d-12
      iresm = 0
      jresm = 0
      kresm = 0
      res = 0.d0
      
      do k=2,resi%kmax
        do j=2,resi%jmax
          do i=2,resi%imax
            resi%qres(:,i,j,k) = dabs(resi%qres(:,i,j,k) - variable%getpv(i,j,k))
            do n=1,resi%npv
              res(n) = res(n) + resi%qres(n,i,j,k)**2
            end do
            restest = resi%qres(1,i,j,k)
            if(restest.gt.resmax) then
              resmax = restest
              iresm=i
              jresm=j
              kresm=k
            end if
          end do
        end do
      end do
      
      mpi_res = res
      call mpi_reduce(mpi_res,mpi_res_sum,resi%npv,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
      call mpi_bcast(mpi_res_sum,resi%npv,mpi_real8,0,mpi_comm_world,ierr)
      
      do n=1,resi%npv
        mpi_res_sum(n) = dsqrt(mpi_res_sum(n))
      end do
      iter = resi%ntmax*(nt_phy-1)+nt
      if(mod(iter,resi%nprint).eq.0) then
        if(resi%rank.eq.0) then
          if(resi%npv.eq.7) then
            write(*,77) iter, mpi_res_sum(1),mpi_res_sum(2),mpi_res_sum(3),mpi_res_sum(4),mpi_res_sum(5),mpi_res_sum(6),mpi_res_sum(7)
          else
            write(*,99) iter, mpi_res_sum(1),mpi_res_sum(2),mpi_res_sum(3),mpi_res_sum(4),mpi_res_sum(5),mpi_res_sum(6),mpi_res_sum(7),mpi_res_sum(8),mpi_res_sum(9)
          end if
          write(*,*) resi%rank,'max err',iresm,jresm,kresm,resmax
        end if
        call mpi_barrier(mpi_comm_world,ierr)
        if(resi%rank.ne.0) then
          write(*,*) resi%rank,'max err',iresm,jresm,kresm,resmax        
        end if
      end if
      
      if(resi%rank.eq.0) then
        if(resi%npv.eq.7) then
          write(resi%io,66) iter, mpi_res_sum(1),mpi_res_sum(2),mpi_res_sum(3),mpi_res_sum(4),mpi_res_sum(5),mpi_res_sum(6),mpi_res_sum(7)
        else
          write(resi%io,88) iter, mpi_res_sum(1),mpi_res_sum(2),mpi_res_sum(3),mpi_res_sum(4),mpi_res_sum(5),mpi_res_sum(6),mpi_res_sum(7),mpi_res_sum(8),mpi_res_sum(9)
        end if
      end if

      if(resi%rank.eq.0) then
        if(mpi_res_sum(1).lt.resi%bond) resi%l_converge = .true.
      end if
      
      call mpi_bcast(resi%l_converge,1,mpi_logical,0,mpi_comm_world,ierr)
      
66  format(i10,7(f30.12))
88  format(i10,9(f30.12))
77  format(i10,2(f30.12)/10x,2(f30.12)/10x,2(f30.12)/10x,2(f30.12))
99  format(i10,2(f30.12)/10x,2(f30.12)/10x,2(f30.12)/10x,2(f30.12)/10x,2(f30.12))
    end subroutine residual
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getconverge(resi)
      implicit none
      class(t_resi), intent(in) :: resi
      logical :: getconverge
      
      getconverge = resi%l_converge
    
    end function getconverge
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module residual_module
