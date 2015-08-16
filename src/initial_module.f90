module initial_module
  use config_module
  use grid_module
  use variable_module
  use eos_module
  use prop_module
  implicit none
  private
  public :: t_ini,t_ini_initial,t_ini_restart
  
  type, abstract :: t_ini
    private
    logical :: l_ini
    integer :: size,rank,imax,jmax,kmax
    integer :: iturb,nsteady,rstnum
    real(8) :: pref,uref,aoa,tref,y1ref,y2ref,kref,oref,emutref
    contains
      procedure :: construct
      procedure :: destruct
      procedure(p_initialize), deferred :: initialize    
  end type t_ini
  
  type, extends(t_ini) :: t_ini_initial
    contains
      procedure :: initialize => initial
  end type t_ini_initial

  type, extends(t_ini) :: t_ini_restart
    contains
      procedure :: initialize => restart
  end type t_ini_restart
  
  abstract interface
    subroutine p_initialize(ini,grid,variable,eos,prop,nps,nts)
      import t_ini
      import t_grid
      import t_variable
      import t_eos
      import t_prop
      implicit none
      class(t_ini), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer, intent(out) :: nps,nts
    end subroutine p_initialize
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(ini,config,grid)
      implicit none
      class(t_ini), intent(out) :: ini
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      
      ini%pref = config%getpref()
      ini%uref = config%geturef()
      ini%aoa = config%getaoa()
      ini%tref = config%gettref()
      ini%y1ref = config%gety1ref()
      ini%y2ref = config%gety2ref()
      ini%iturb = config%getiturb()
      ini%nsteady = config%getnsteady()
      ini%rstnum = config%getrstnum()
      ini%size = config%getsize()
      ini%rank = config%getrank()
      ini%kref = config%getkref()
      ini%oref = config%getoref()
      ini%emutref = config%getemutref()
      ini%imax = grid%getimax()
      ini%jmax = grid%getjmax()
      ini%kmax = grid%getkmax()
      
      ini%l_ini = .true.
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(ini)
      implicit none
      class(t_ini), intent(inout) :: ini

      ini%l_ini = .false.
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine restart(ini,grid,variable,eos,prop,nps,nts)
      implicit none
      include 'mpif.h'
      class(t_ini_restart), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer, intent(out) :: nps,nts
      integer :: i,j,k,l,n,m,io,mm,ier
      real(8) :: dv(variable%getndv())
      real(8), dimension(:), allocatable :: qq_temp
      real(8), dimension(:,:,:,:), allocatable :: pv,tv
      real(8), dimension(:,:,:,:,:), allocatable :: qq
      character(7) :: iter_tag
      integer :: size,rank,imax,jmax,kmax,nqq
      real(8) :: var1,var2,var3
      
      allocate(pv(variable%getnpv(),ini%imax,ini%jmax,ini%kmax))
      allocate(tv(variable%getntv(),ini%imax,ini%jmax,ini%kmax))
      allocate(qq_temp(variable%getnpv()))
      
      if(ini%nsteady.eq.1) then
        write(iter_tag,'(i4.4)') ini%rstnum
      else
        write(iter_tag,'(i7.7)') ini%rstnum
      end if

      do mm=0,ini%size-1
        if(mm.eq.ini%rank) then
          open(newunit=io,file='./out_'//trim(iter_tag)//'.dat',status='old',action='read',form='unformatted')
          read(io) size,nps,nts,nqq
          allocate(qq(nqq,variable%getnpv(),ini%imax,ini%jmax,ini%kmax))
          if(size.ne.ini%size) write(*,*) 'invalid size'
          
          do m=0,ini%size-1
            read(io) rank,imax,jmax,kmax
            if((rank.eq.ini%rank).and.(imax.eq.ini%imax).and.(jmax.eq.ini%jmax).and.(kmax.eq.ini%kmax)) then
              read(io) ((((pv(n,i,j,k),n=1,variable%getnpv()),i=2,imax),j=2,jmax),k=2,kmax)
              read(io) ((((tv(n,i,j,k),n=1,variable%getntv()),i=2,imax),j=2,jmax),k=2,kmax)
              read(io) (((((qq(l,n,i,j,k),l=1,nqq),n=1,variable%getnpv()),i=2,imax),j=2,jmax),k=2,kmax)    
            else       
              read(io) ((((var1,n=1,variable%getnpv()),i=2,imax),j=2,jmax),k=2,kmax)
              read(io) ((((var2,n=1,variable%getntv()),i=2,imax),j=2,jmax),k=2,kmax)
              read(io) (((((var3,l=1,nqq),n=1,variable%getnpv()),i=2,imax),j=2,jmax),k=2,kmax)
            end if
          end do
      
          close(io)
        end if
        call mpi_barrier(mpi_comm_world,ier)
      end do
      
      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            call variable%setpv(1,i,j,k,pv(1,i,j,k)-ini%pref)
            do n=2,variable%getnpv()
              call variable%setpv(n,i,j,k,pv(n,i,j,k))
            end do
        
            do n=1,variable%getntv()
              call variable%settv(n,i,j,k,tv(n,i,j,k))
            end do      
            
            call eos%deteos(pv(1,i,j,k),pv(5,i,j,k),pv(6,i,j,k),pv(7,i,j,k),dv)
            
            do n=1,variable%getndv()
              call variable%setdv(n,i,j,k,dv(n))
            end do
        
            if(nqq.ne.variable%getnqq()) then
              qq_temp(1) = dv(1)
              qq_temp(2) = dv(1)*pv(2,i,j,k)
              qq_temp(3) = dv(1)*pv(3,i,j,k)
              qq_temp(4) = dv(1)*pv(4,i,j,k)
              qq_temp(5) = dv(1)*(dv(2)+0.5d0*(pv(2,i,j,k)**2+pv(3,i,j,k)**2+pv(4,i,j,k)**2))-pv(1,i,j,k)
              do n=6,variable%getnpv()
                qq_temp(n) = dv(1)*pv(n,i,j,k)
              end do
              
              call variable%setqq(1,i,j,k,qq_temp)
              call variable%setqq(2,i,j,k,qq_temp)
            else
              do n=1,variable%getnqq()
                qq_temp = qq(n,:,i,j,k)
                call variable%setqq(n,i,j,k,qq_temp)
              end do          
            end if
          end do
        end do
      end do
      if(allocated(qq_temp)) deallocate(qq_temp)
      deallocate(pv,tv,qq)

      nts = nts + 1
      nps = nps + 1
      
      if(ini%nsteady.eq.1) then
        nts = 1
      end if
      if(nqq.eq.0) then
        if(ini%nsteady.eq.1) then
          write(*,*) 'invalid restart option : unsteady -> steady'
        else
          nps = 1
        end if
      end if
      
    end subroutine restart
#ifdef test
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initial(ini,grid,variable,eos,prop,nps,nts)
      implicit none
      class(t_ini_initial), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer, intent(out) :: nps,nts
      integer :: i,j,k,n
      real(8) :: pv(variable%getnpv()),dv(variable%getndv())
      real(8) :: tv(variable%getntv()),qq(variable%getnpv())
      
      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            call variable%setpv(1,i,j,k,0.d0)
            call variable%setpv(2,i,j,k,ini%uref*dcos(ini%aoa))
            call variable%setpv(3,i,j,k,ini%uref*dsin(ini%aoa))
            call variable%setpv(4,i,j,k,0.d0)
            call variable%setpv(5,i,j,k,ini%tref)
            call variable%setpv(6,i,j,k,ini%y1ref)
            call variable%setpv(7,i,j,k,ini%y2ref)
            
            pv = variable%getpv(i,j,k)   
        
            call eos%deteos(pv(1)+ini%pref,pv(5),pv(6),pv(7),dv)
            
            do n=1,variable%getndv()
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            if(ini%iturb.ge.-2) then
              call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            end if
            
            if(ini%iturb.ge.-1) then
              tv(3) = ini%emutref
              call variable%setpv(8,i,j,k,ini%kref)
              call variable%setpv(9,i,j,k,ini%oref)          
            end if
            
            do n=1,variable%getntv()
              call variable%settv(n,i,j,k,tv(n))
            end do
            
            if(ini%nsteady.eq.1) then
              qq(1) = dv(1)
              qq(2) = dv(1)*pv(2)
              qq(3) = dv(1)*pv(3)
              qq(4) = dv(1)*pv(4)
              qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-ini%pref
              do n=6,variable%getnpv()
                qq(n) = dv(1)*pv(n)
              end do
              
              call variable%setqq(1,i,j,k,qq)
              call variable%setqq(2,i,j,k,qq)
            end if
          end do
        end do
      end do
      nps = 1
      nts = 1
    end subroutine initial
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#elif test1
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initial(ini,grid,variable,eos,prop,nps,nts)
      implicit none
      class(t_ini_initial), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer, intent(out) :: nps,nts
      integer :: i,j,k,n
      real(8) :: pv(variable%getnpv()),dv(variable%getndv())
      real(8) :: tv(variable%getntv()),qq(variable%getnpv())
      real(8) :: x(5)
      
      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            x = grid%getgrd(i,j,k)
            if(x(2).lt.5.d0) then
           ! if(((1.d0-x(3))/(x(2)-5.d0).le.-1.d0).and.((1.d0-x(4))/(x(2)-5.d0).le.-1.d0)) then
           ! call variable%setpv(1,i,j,k,10132.5-ini%pref)
           ! call variable%setpv(2,i,j,k,0.d0)
           ! call variable%setpv(3,i,j,k,0.d0)
           ! call variable%setpv(4,i,j,k,0.d0)
           ! call variable%setpv(5,i,j,k,240.d0)
           ! call variable%setpv(6,i,j,k,ini%y1ref)
           ! call variable%setpv(7,i,j,k,ini%y2ref)  
           !else
            call variable%setpv(1,i,j,k,0.d0)
            call variable%setpv(2,i,j,k,0.d0)
            call variable%setpv(3,i,j,k,0.d0)
            call variable%setpv(4,i,j,k,0.d0)
            call variable%setpv(5,i,j,k,ini%tref)
            call variable%setpv(6,i,j,k,ini%y1ref)
            call variable%setpv(7,i,j,k,ini%y2ref)
            !end if
            else
             
            if((x(3)/(x(2)-5.d0).ge.1.d0)) then
              call variable%setpv(1,i,j,k,0.d0)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,ini%tref)
              call variable%setpv(6,i,j,k,ini%y1ref)
              call variable%setpv(7,i,j,k,ini%y2ref)
           else
            
            
            call variable%setpv(1,i,j,k,10132.5-ini%pref)
            call variable%setpv(2,i,j,k,0.d0)
            call variable%setpv(3,i,j,k,0.d0)
            call variable%setpv(4,i,j,k,0.d0)
            call variable%setpv(5,i,j,k,240.d0)
            call variable%setpv(6,i,j,k,ini%y1ref)
            call variable%setpv(7,i,j,k,ini%y2ref)    
            end if
            end if
            pv = variable%getpv(i,j,k)   
        
            call eos%deteos(pv(1)+ini%pref,pv(5),pv(6),pv(7),dv)
            
            do n=1,variable%getndv()
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            if(ini%iturb.ge.-2) then
              call prop%detprop(dv(3),dv(4),dv(5),pv(5),pv(6),pv(7),tv(1:2))
            end if
            
            if(ini%iturb.ge.-1) then
              tv(3) = ini%emutref
              call variable%setpv(8,i,j,k,ini%kref)
              call variable%setpv(9,i,j,k,ini%oref)          
            end if
            
            do n=1,variable%getntv()
              call variable%settv(n,i,j,k,tv(n))
            end do
            
            if(ini%nsteady.eq.1) then
              qq(1) = dv(1)
              qq(2) = dv(1)*pv(2)
              qq(3) = dv(1)*pv(3)
              qq(4) = dv(1)*pv(4)
              qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-ini%pref
              do n=6,variable%getnpv()
                qq(n) = dv(1)*pv(n)
              end do
              
              call variable%setqq(1,i,j,k,qq)
              call variable%setqq(2,i,j,k,qq)
            end if
          end do
        end do
      end do
      nps = 1
      nts = 1
    end subroutine initial
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif
end module initial_module

