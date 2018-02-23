module initial_module
  use mpi
  use config_module
  use grid_module
  use variable_module
  use eos_module
  implicit none
  private
  public :: t_ini,t_ini_initial,t_ini_initial_rot,t_ini_restart
  
  type, abstract :: t_ini
    private
    logical :: l_ini
    integer :: npv,ndv,ntv,nqq
    integer :: size,rank,imax,jmax,kmax
    integer :: iturb,nsteady,rstnum
    real(8) :: pref,uref,aoa,aos,tref,y1ref,y2ref,kref,oref,emutref,omega(3)
    contains
      procedure :: construct
      procedure :: destruct
      procedure(p_initialize), deferred :: initialize
  end type t_ini
  
  type, extends(t_ini) :: t_ini_initial
    contains
      procedure :: initialize => initial
  end type t_ini_initial

  type, extends(t_ini) :: t_ini_initial_rot
    contains
      procedure :: initialize => initial_rot
  end type t_ini_initial_rot

  type, extends(t_ini) :: t_ini_restart
    contains
      procedure :: initialize => restart
  end type t_ini_restart
  
  abstract interface
    subroutine p_initialize(ini,grid,variable,eos,nps,nts,timeprev)
      import t_ini
      import t_grid
      import t_variable
      import t_eos
      implicit none
      class(t_ini), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
    end subroutine p_initialize
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(ini,config,grid)
      implicit none
      class(t_ini), intent(out) :: ini
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid

      ini%npv = config%getnpv()
      ini%ndv = config%getndv()
      ini%ntv = config%getntv()
      ini%nqq = config%getnqq()
      ini%pref = config%getpref()
      ini%uref = config%geturef()
      ini%aoa = config%getaoa()
      ini%aos = config%getaos()
      ini%tref = config%gettref()
      ini%y1ref = config%gety1ref()
      ini%y2ref = config%gety2ref()
      ini%iturb = config%getiturb()
      ini%nsteady = config%getnsteady()
      ini%rstnum = config%getrstnum()
      ini%size = config%getsize()
      ini%rank = config%getrank()
      if(config%getiturb().ge.-1) then
        ini%kref = config%getkref()
        ini%oref = config%getoref()
        ini%emutref = config%getemutref()
      end if
      ini%omega = config%getomega()
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
#ifdef steadyini
    subroutine restart(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_restart), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n,io,ier,num
      integer :: intsize,realsize
      integer(kind=mpi_offset_kind) :: disp
      real(8) :: dv(ini%ndv),qq_temp(ini%npv)
      real(8), dimension(:,:,:,:), allocatable :: pv,tv
      character(8) :: iter_tag

      call mpi_type_size(mpi_integer,intsize,ier)
      call mpi_type_size(mpi_real8,realsize,ier)

      write(iter_tag,'(i8.8)') ini%rstnum

      if(ini%rank.eq.0) then
        disp = 0
      else
        disp = 0
        do i=0,ini%rank-1
          disp = disp + intsize*2 + realsize &
               + realsize*ini%npv*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5) &
               + realsize*ini%ntv*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5)
        end do
      end if

      call mpi_file_open(mpi_comm_world,"./out_"//trim(iter_tag)//".dat",mpi_mode_rdonly,mpi_info_null,io,ier)

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,nps,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + intsize

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,nts,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + intsize

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,timeprev,1,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize

      allocate(pv(ini%npv,-1:ini%imax+3,-1:ini%jmax+3,-1:ini%kmax+3))
      allocate(tv(ini%ntv,-1:ini%imax+3,-1:ini%jmax+3,-1:ini%kmax+3))

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = ini%npv*(ini%imax+5)*(ini%jmax+5)*(ini%kmax+5)
      call mpi_file_read_all(io,pv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize*num

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = ini%ntv*(ini%imax+5)*(ini%jmax+5)*(ini%kmax+5)
      call mpi_file_read_all(io,tv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize*num

      call mpi_file_close(io,ier)

      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            do n=1,ini%npv
              call variable%setpv(n,i,j,k,pv(n,i,j,k))
            end do
        
            do n=1,ini%ntv
              call variable%settv(n,i,j,k,tv(n,i,j,k))
            end do      
            
            call eos%deteos_simple(pv(1,i,j,k)+ini%pref,pv(5,i,j,k),pv(6,i,j,k),pv(7,i,j,k),dv)

            do n=1,ini%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do

            qq_temp(1) = dv(1)
            qq_temp(2) = dv(1)*pv(2,i,j,k)
            qq_temp(3) = dv(1)*pv(3,i,j,k)
            qq_temp(4) = dv(1)*pv(4,i,j,k)
            qq_temp(5) = dv(1)*(dv(2)+0.5d0*(pv(2,i,j,k)**2+pv(3,i,j,k)**2+pv(4,i,j,k)**2))-pv(1,i,j,k)-ini%pref
            do n=6,ini%npv
              qq_temp(n) = dv(1)*pv(n,i,j,k)
            end do

            call variable%setqq(1,i,j,k,qq_temp)
            call variable%setqq(2,i,j,k,qq_temp)

          end do
        end do
      end do

      if(allocated(pv)) deallocate(pv)
      if(allocated(tv)) deallocate(tv)

      nts = 1
      nps = 1
      timeprev = 0.d0

    end subroutine restart
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#else
    subroutine restart(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_restart), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n,io,ier,num
      integer :: intsize,realsize
      integer(kind=mpi_offset_kind) :: disp
      real(8) :: dv(ini%ndv),qq_temp(ini%npv)
      real(8), dimension(:,:,:,:), allocatable :: pv,tv
      real(8), dimension(:,:,:,:,:), allocatable :: qq
      character(8) :: iter_tag

      call mpi_type_size(mpi_integer,intsize,ier)
      call mpi_type_size(mpi_real8,realsize,ier)

      write(iter_tag,'(i8.8)') ini%rstnum

      if(ini%rank.eq.0) then
        disp = 0
      else
        disp = 0
        do i=0,ini%rank-1
          disp = disp + intsize*2 + realsize &
               + realsize*ini%npv*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5) &
               + realsize*ini%ntv*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5) &
               + realsize*ini%nqq*ini%npv*(grid%getimax_zone(i)-1)*(grid%getjmax_zone(i)-1)*(grid%getkmax_zone(i)-1)
        end do
      end if

      call mpi_file_open(mpi_comm_world,"./out_"//trim(iter_tag)//".dat",mpi_mode_rdonly,mpi_info_null,io,ier)

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,nps,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + intsize

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,nts,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + intsize

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,timeprev,1,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize

      allocate(pv(ini%npv,-1:ini%imax+3,-1:ini%jmax+3,-1:ini%kmax+3))
      allocate(tv(ini%ntv,-1:ini%imax+3,-1:ini%jmax+3,-1:ini%kmax+3))
      allocate(qq(ini%npv,ini%nqq,2:ini%imax,2:ini%jmax,2:ini%kmax))

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = ini%npv*(ini%imax+5)*(ini%jmax+5)*(ini%kmax+5)
      call mpi_file_read_all(io,pv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize*num

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = ini%ntv*(ini%imax+5)*(ini%jmax+5)*(ini%kmax+5)
      call mpi_file_read_all(io,tv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize*num

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = ini%nqq*ini%npv*(ini%imax-1)*(ini%jmax-1)*(ini%kmax-1)
      call mpi_file_read_all(io,qq,num,mpi_real8,mpi_status_ignore,ier)

      call mpi_file_close(io,ier)

      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            do n=1,ini%npv
              call variable%setpv(n,i,j,k,pv(n,i,j,k))
            end do
        
            do n=1,ini%ntv
              call variable%settv(n,i,j,k,tv(n,i,j,k))
            end do      
            
            call eos%deteos_simple(pv(1,i,j,k)+ini%pref,pv(5,i,j,k),pv(6,i,j,k),pv(7,i,j,k),dv)
            
            do n=1,ini%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do

            do n=1,ini%nqq
              qq_temp = qq(:,n,i,j,k)
              call variable%setqq(n,i,j,k,qq_temp)
            end do
          end do
        end do
      end do

      if(allocated(pv)) deallocate(pv)
      if(allocated(tv)) deallocate(tv)
      if(allocated(qq)) deallocate(qq)

      nts = nts + 1
      nps = nps + 1
      
      if(ini%nsteady.eq.1) then
        nts = 1
      else
        nps = 1
      end if
      
    end subroutine restart
#endif
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initial_rot(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_initial_rot), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n
      real(8) :: pv(ini%npv),dv(ini%ndv)
      real(8) :: tv(ini%ntv),qq(ini%npv)
      real(8) :: grd(grid%getngrd()),gridvel(3)

      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            grd = grid%getgrd(i,j,k)

            gridvel(1) = ini%omega(2)*grd(4)-ini%omega(3)*grd(3)
            gridvel(2) = ini%omega(3)*grd(2)-ini%omega(1)*grd(4)
            gridvel(3) = ini%omega(1)*grd(3)-ini%omega(2)*grd(2)

            call variable%setpv(1,i,j,k,0.d0)
            call variable%setpv(2,i,j,k,-gridvel(1))
            call variable%setpv(3,i,j,k,-gridvel(2))
            call variable%setpv(4,i,j,k,1.559373336d0-gridvel(3))
            call variable%setpv(5,i,j,k,ini%tref)
            call variable%setpv(6,i,j,k,ini%y1ref)
            call variable%setpv(7,i,j,k,ini%y2ref)

            pv = variable%getpv(i,j,k)

            call eos%deteos(pv(1)+ini%pref,pv(5),pv(6),pv(7),dv,tv)

            do n=1,ini%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do

            if(ini%iturb.ge.-1) then
              tv(3) = ini%emutref
              call variable%setpv(8,i,j,k,ini%kref)
              call variable%setpv(9,i,j,k,ini%oref)
            end if

            do n=1,ini%ntv
              call variable%settv(n,i,j,k,tv(n))
            end do

            if(ini%nsteady.eq.1) then
              qq(1) = dv(1)
              qq(2) = dv(1)*pv(2)
              qq(3) = dv(1)*pv(3)
              qq(4) = dv(1)*pv(4)
              qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-ini%pref
              do n=6,ini%npv
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
      timeprev = 0.d0
    end subroutine initial_rot
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef test

#elif lax3d
     subroutine initial(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_initial), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n
      real(8) :: pv(ini%npv),dv(ini%ndv)
      real(8) :: tv(ini%ntv),qq(ini%npv)
      real(8) :: x(5)
      
      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            x = grid%getgrd(i,j,k)

            if(x(2)-x(3)-x(4).lt.5.d0) then !left
              call variable%setpv(1,i,j,k,3.528d0-ini%pref)
              call variable%setpv(2,i,j,k,0.698d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,3.528d0*1.4d0/(0.4d0*1004.64d0*0.445d0))
              call variable%setpv(6,i,j,k,ini%y1ref)
              call variable%setpv(7,i,j,k,ini%y2ref)
            else !right
              call variable%setpv(1,i,j,k,0.571d0-ini%pref)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,0.571d0*1.4d0/(0.4d0*1004.64d0*0.5d0))
              call variable%setpv(6,i,j,k,ini%y1ref)
              call variable%setpv(7,i,j,k,ini%y2ref)
            end if
            pv = variable%getpv(i,j,k)   
        
            call eos%deteos(pv(1)+ini%pref,pv(5),pv(6),pv(7),dv,tv)
            
            do n=1,ini%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            if(ini%iturb.ge.-1) then
              tv(3) = ini%emutref
              call variable%setpv(8,i,j,k,ini%kref)
              call variable%setpv(9,i,j,k,ini%oref)          
            end if
            
            do n=1,ini%ntv
              call variable%settv(n,i,j,k,tv(n))
            end do
            
            if(ini%nsteady.eq.1) then
              qq(1) = dv(1)
              qq(2) = dv(1)*pv(2)
              qq(3) = dv(1)*pv(3)
              qq(4) = dv(1)*pv(4)
              qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-ini%pref
              do n=6,ini%npv
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
      timeprev = 0.d0
    end subroutine initial

#elif shocktube
    subroutine initial(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_initial), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n
      real(8) :: pv(ini%npv),dv(ini%ndv)
      real(8) :: tv(ini%ntv),qq(ini%npv)
      real(8) :: grd(grid%getngrd())
      
      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            grd = grid%getgrd(i,j,k)
            if(grd(2).lt.0.5d0) then
              call variable%setpv(1,i,j,k,120.d0/1.4d0-ini%pref)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,ini%tref)
              call variable%setpv(6,i,j,k,ini%y1ref)
              call variable%setpv(7,i,j,k,ini%y2ref)
            else
              call variable%setpv(1,i,j,k,1.2d0/1.4d0-ini%pref)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,ini%tref)
              call variable%setpv(6,i,j,k,ini%y1ref)
              call variable%setpv(7,i,j,k,ini%y2ref)
            end if
            pv = variable%getpv(i,j,k)   
        
            call eos%deteos(pv(1)+ini%pref,pv(5),pv(6),pv(7),dv,tv)
            
            do n=1,ini%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            if(ini%iturb.ge.-1) then
              tv(3) = ini%emutref
              call variable%setpv(8,i,j,k,ini%kref)
              call variable%setpv(9,i,j,k,ini%oref)          
            end if
            
            do n=1,ini%ntv
              call variable%settv(n,i,j,k,tv(n))
            end do
            
            if(ini%nsteady.eq.1) then
              qq(1) = dv(1)
              qq(2) = dv(1)*pv(2)
              qq(3) = dv(1)*pv(3)
              qq(4) = dv(1)*pv(4)
              qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-ini%pref
              do n=6,ini%npv
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
      timeprev = 0.d0
    end subroutine initial

#elif sloshing_3d
    subroutine initial(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_initial), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n
      real(8) :: pv(ini%npv),dv(ini%ndv)
      real(8) :: tv(ini%ntv),qq(ini%npv)
      real(8) :: grd(grid%getngrd())

    ! tank geometry : sivb tank
    ! test paper : Analysis of cryogenic propellant behaviour in microgravity and low thrust environments, MF Fisher, 1992 

      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            grd = grid%getgrd(i,j,k)

    ! real grid is 43.3 times smaller than test grid
          ! dimethylether
            if(dsqrt((grd(2)-0.05765d0)**2+(grd(3)+0.04878d0)**2).ge.0.19d0) then
              call variable%setpv(1,i,j,k,0.61793d6-ini%pref)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,299.85d0-ini%tref)
              call variable%setpv(6,i,j,k,0.d0)
              call variable%setpv(7,i,j,k,0.d0)
          ! air
            else if(dsqrt((grd(2)-0.05765d0)**2+(grd(3)+0.04878d0)**2).lt.0.19d0) then
              call variable%setpv(1,i,j,k,0.61793d6-ini%pref)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,299.85d0-ini%tref)
              call variable%setpv(6,i,j,k,0.d0)
              call variable%setpv(7,i,j,k,1.d0)
            end if
          end do
        end do
      end do

      call set_others(ini,variable,eos,nps,nts,timeprev)

    end subroutine initial

#elif test4tank
    subroutine initial(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_initial), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n
      real(8) :: pv(ini%npv),dv(ini%ndv)
      real(8) :: tv(ini%ntv),qq(ini%npv)
      real(8) :: grd(grid%getngrd())

    ! tank geometry : sivb tank
    ! test paper : Analysis of cryogenic propellant behaviour in microgravity and low thrust environments, MF Fisher, 1992 

      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            grd = grid%getgrd(i,j,k)

    ! real grid is 43.3 times smaller than test grid
          ! dimethylether
            if(grd(1).ge.10.d0) then
              call variable%setpv(1,i,j,k,0.61793d6-ini%pref)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,299.85d0-ini%tref)
              call variable%setpv(6,i,j,k,0.d0)
              call variable%setpv(7,i,j,k,0.d0)
          ! air
            else if(grd(1).lt.10.d0) then
              call variable%setpv(1,i,j,k,0.61793d6-ini%pref)
              call variable%setpv(2,i,j,k,0.d0)
              call variable%setpv(3,i,j,k,0.d0)
              call variable%setpv(4,i,j,k,0.d0)
              call variable%setpv(5,i,j,k,299.85d0-ini%tref)
              call variable%setpv(6,i,j,k,0.d0)
              call variable%setpv(7,i,j,k,1.d0)
            end if
          end do
        end do
      end do

      call set_others(ini,variable,eos,nps,nts,timeprev)

    end subroutine initial

#else
    subroutine initial(ini,grid,variable,eos,nps,nts,timeprev)
      implicit none
      class(t_ini_initial), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n
      real(8) :: pv(ini%npv),dv(ini%ndv)
      real(8) :: tv(ini%ntv),qq(ini%npv)
      
      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            call variable%setpv(1,i,j,k,0.d0)
            call variable%setpv(2,i,j,k,ini%uref*dcos(ini%aos)*dcos(ini%aoa))
            call variable%setpv(3,i,j,k,ini%uref*dcos(ini%aos)*dsin(ini%aoa))
            call variable%setpv(4,i,j,k,ini%uref*dsin(ini%aos))
            call variable%setpv(5,i,j,k,ini%tref)
            call variable%setpv(6,i,j,k,ini%y1ref)
            call variable%setpv(7,i,j,k,ini%y2ref)
            
            pv = variable%getpv(i,j,k)   
        
            call eos%deteos(pv(1)+ini%pref,pv(5),pv(6),pv(7),dv,tv)
            
            do n=1,ini%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            if(ini%iturb.ge.-1) then
              tv(3) = ini%emutref
              call variable%setpv(8,i,j,k,ini%kref)
              call variable%setpv(9,i,j,k,ini%oref)          
            end if
            
            do n=1,ini%ntv
              call variable%settv(n,i,j,k,tv(n))
            end do
            
            if(ini%nsteady.eq.1) then
              qq(1) = dv(1)
              qq(2) = dv(1)*pv(2)
              qq(3) = dv(1)*pv(3)
              qq(4) = dv(1)*pv(4)
              qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-ini%pref
              do n=6,ini%npv
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
      timeprev = 0.d0
    end subroutine initial
#endif
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_others(ini,variable,eos,nps,nts,timeprev)
      implicit none
      type(t_ini_initial), intent(in) :: ini
      type(t_variable), intent(inout) :: variable
      class(t_eos), intent(in) :: eos
      integer, intent(out) :: nps,nts
      real(8), intent(out) :: timeprev
      integer :: i,j,k,n
      real(8) :: pv(ini%npv),dv(ini%ndv),tv(ini%ntv),qq(ini%npv)

      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            pv = variable%getpv(i,j,k)   
write(*,*) 'p:',pv(1),'T:',pv(5)
            call eos%deteos(pv(1)+ini%pref,pv(5)+ini%tref,pv(6),pv(7),dv,tv)
            
            do n=1,ini%ndv
              call variable%setdv(n,i,j,k,dv(n))
            end do
            
            if(ini%iturb.ge.-1) then
              tv(3) = ini%emutref
              call variable%setpv(8,i,j,k,ini%kref)
              call variable%setpv(9,i,j,k,ini%oref)          
            end if
            
            do n=1,ini%ntv
              call variable%settv(n,i,j,k,tv(n))
            end do
            
            if(ini%nsteady.eq.1) then
              qq(1) = dv(1)
              qq(2) = dv(1)*pv(2)
              qq(3) = dv(1)*pv(3)
              qq(4) = dv(1)*pv(4)
              qq(5) = dv(1)*(dv(2)+0.5d0*(pv(2)**2+pv(3)**2+pv(4)**2))-pv(1)-ini%pref
              do n=6,ini%npv
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
      timeprev = 0.d0
    end subroutine set_others
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module initial_module

