module postvariable_module
  use config_module
  use eos_module
  implicit none
  private
  public :: t_variable
  
  type t_domain
    integer :: rank,imax,jmax,kmax
    real(8), dimension(:,:,:,:), allocatable :: pv,tv,dv
  end type t_domain
  
  type t_solution
    integer :: size,nps,nts
    type(t_domain), dimension(:), allocatable :: domain
  end type t_solution
  
  type t_variable
    private
    integer :: npv,ntv,ndv,nsolution
    type(t_solution), dimension(:), allocatable :: solution
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: getnps
      procedure :: getnts
      procedure :: getnpv
      procedure :: getntv
      procedure :: getndv
      procedure :: getnsolution
      procedure :: getpv
      procedure :: gettv
      procedure :: getdv
      procedure :: getsize
      procedure :: getrank
      procedure :: getimax
      procedure :: getjmax
      procedure :: getkmax
  end type t_variable

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(variable,config,eos,istart,iend,nsolution)
      implicit none
      class(t_variable), intent(out) :: variable
      type(t_config), intent(in) :: config
      type(t_eos), intent(in) :: eos
      integer, intent(in) :: istart,iend,nsolution
      character(7) :: iter_tag
      real(8) :: var
      integer :: i,j,k,n,m,l,o,io,nqq,iter
      
      select case(config%getiturb())
      case(-3)
        variable%npv = 7
        variable%ntv = 0  
      case(-2)
        variable%npv = 7
        variable%ntv = 2      
      case(-1,0)
        variable%npv = 9
        variable%ntv = 3      
      end select
      variable%ndv = 18

      variable%nsolution = nsolution
      
      if(variable%nsolution.eq.1) then
        if(iend.eq.istart) then
          iter = 0
        else 
          write(*,*) 'invalid nsolution'
        end if
      else if(mod((iend-istart)/config%getnexport(),variable%nsolution-1).ne.0) then
        write(*,*) 'invalid nsolution'
      else
        iter = (iend-istart)/config%getnexport()/(variable%nsolution-1)
      end if
      
      allocate(variable%solution(variable%nsolution))
      
      do l=1,variable%nsolution
        if(config%getnsteady().eq.1) then
          write(iter_tag,'(i4.4)') istart+iter*config%getnexport()*(l-1)
        else
          write(iter_tag,'(i7.7)') istart+iter*config%getnexport()*(l-1)
        end if
        
        open(newunit=io,file='./out_'//trim(iter_tag)//'.dat',status='old',action='read',form='unformatted')
        read(io) variable%solution(l)%size,variable%solution(l)%nps,variable%solution(l)%nts,nqq
        
        allocate(variable%solution(l)%domain(1:variable%solution(l)%size))
        
        do m=1,variable%solution(l)%size
          read(io) variable%solution(l)%domain(m)%rank,variable%solution(l)%domain(m)%imax,variable%solution(l)%domain(m)%jmax,variable%solution(l)%domain(m)%kmax
          
          allocate(variable%solution(l)%domain(m)%pv(variable%npv,variable%solution(l)%domain(m)%imax,variable%solution(l)%domain(m)%jmax,variable%solution(l)%domain(m)%kmax))
          allocate(variable%solution(l)%domain(m)%dv(variable%ndv,variable%solution(l)%domain(m)%imax,variable%solution(l)%domain(m)%jmax,variable%solution(l)%domain(m)%kmax))
          allocate(variable%solution(l)%domain(m)%tv(variable%ntv,variable%solution(l)%domain(m)%imax,variable%solution(l)%domain(m)%jmax,variable%solution(l)%domain(m)%kmax))
          
          read(io) ((((variable%solution(l)%domain(m)%pv(n,i,j,k),n=1,variable%npv),i=2,variable%solution(l)%domain(m)%imax),j=2,variable%solution(l)%domain(m)%jmax),k=2,variable%solution(l)%domain(m)%kmax)
          read(io) ((((variable%solution(l)%domain(m)%tv(n,i,j,k),n=1,variable%ntv),i=2,variable%solution(l)%domain(m)%imax),j=2,variable%solution(l)%domain(m)%jmax),k=2,variable%solution(l)%domain(m)%kmax)
          read(io) (((((var,o=1,nqq),n=1,variable%npv),i=2,variable%solution(l)%domain(m)%imax),j=2,variable%solution(l)%domain(m)%jmax),k=2,variable%solution(l)%domain(m)%kmax)
          
          do k=2,variable%solution(l)%domain(m)%kmax
            do j=2,variable%solution(l)%domain(m)%jmax
              do i=2,variable%solution(l)%domain(m)%imax
                call eos%deteos(variable%solution(l)%domain(m)%pv(1,i,j,k),variable%solution(l)%domain(m)%pv(5,i,j,k),variable%solution(l)%domain(m)%pv(6,i,j,k),variable%solution(l)%domain(m)%pv(7,i,j,k) &
                                ,variable%solution(l)%domain(m)%dv(:,i,j,k)) 
              end do
            end do
          end do
        end do
        
        close(io)
      end do
  
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(variable)
      implicit none
      class(t_variable), intent(inout) :: variable
      integer :: l,m
      
      do l=1,variable%nsolution
        do m=1,variable%solution(l)%size
          deallocate(variable%solution(l)%domain(m)%pv)
          deallocate(variable%solution(l)%domain(m)%dv)
          deallocate(variable%solution(l)%domain(m)%tv)
        end do
        deallocate(variable%solution(l)%domain)
      end do
      
      deallocate(variable%solution)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getnps(variable,l)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l
      integer :: getnps
      
      getnps = variable%solution(l)%nps
    end function getnps
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getnts(variable,l)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l
      integer :: getnts
      
      getnts = variable%solution(l)%nts
    end function getnts
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getnpv(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getnpv
      
      getnpv = variable%npv
      
    end function getnpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getntv(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getntv
      
      getntv = variable%ntv
      
    end function getntv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getndv(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getndv
      
      getndv = variable%ndv
      
    end function getndv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getnsolution(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getnsolution
      
      getnsolution = variable%nsolution
      
    end function getnsolution
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getpv(variable,l,m,n,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m,n,i,j,k
      real(8) :: getpv
      
      getpv = variable%solution(l)%domain(m)%pv(n,i,j,k)
    end function getpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function gettv(variable,l,m,n,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m,n,i,j,k
      real(8) :: gettv
      
      gettv = variable%solution(l)%domain(m)%tv(n,i,j,k)
    end function gettv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getdv(variable,l,m,n,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m,n,i,j,k
      real(8) :: getdv
      
      getdv = variable%solution(l)%domain(m)%dv(n,i,j,k)
    end function getdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getsize(variable,l)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l
      integer :: getsize
      getsize = variable%solution(l)%size
    end function getsize
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getrank(variable,l,m)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m
      integer :: getrank
      
      getrank = variable%solution(l)%domain(m)%rank
    end function getrank
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getimax(variable,l,m)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m
      integer :: getimax
      
      getimax = variable%solution(l)%domain(m)%imax
    end function getimax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getjmax(variable,l,m)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m
      integer :: getjmax
      
      getjmax = variable%solution(l)%domain(m)%jmax
    end function getjmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pure function getkmax(variable,l,m)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m
      integer :: getkmax
      
      getkmax = variable%solution(l)%domain(m)%kmax
    end function getkmax
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module postvariable_module
