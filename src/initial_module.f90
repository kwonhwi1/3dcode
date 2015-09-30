module initial_module
  use mpi
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
    real(8) :: pref,uref,aoa,aos,tref,y1ref,y2ref,kref,oref,emutref
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
    subroutine construct(ini,config,grid,variable)
      implicit none
      class(t_ini), intent(out) :: ini
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable

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
      class(t_ini_restart), intent(inout) :: ini
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(inout) :: variable
      type(t_eos), intent(in) :: eos
      type(t_prop), intent(in) :: prop
      integer, intent(out) :: nps,nts
      integer :: i,j,k,n,io,ier,num
      integer :: intsize,realsize
      integer(kind=mpi_offset_kind) :: disp
      real(8) :: dv(variable%getndv()),qq_temp(variable%getnpv())
      real(8), dimension(:,:,:,:), allocatable :: pv,tv
      real(8), dimension(:,:,:,:,:), allocatable :: qq
      character(7) :: iter_tag

      call mpi_type_size(mpi_integer,intsize,ier)
      call mpi_type_size(mpi_real8,realsize,ier)

      if(ini%nsteady.eq.1) then
        write(iter_tag,'(i4.4)') ini%rstnum
      else
        write(iter_tag,'(i7.7)') ini%rstnum
      end if

      if(ini%rank.eq.0) then
        disp = 0
      else
        disp = 0
        do i=0,ini%rank-1
          disp = disp + intsize*2 &
               + realsize*variable%getnpv()*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5) &
               + realsize*variable%getntv()*(grid%getimax_zone(i)+5)*(grid%getjmax_zone(i)+5)*(grid%getkmax_zone(i)+5) &
               + realsize*variable%getnqq()*variable%getnpv()*(grid%getimax_zone(i)-1)*(grid%getjmax_zone(i)-1)*(grid%getkmax_zone(i)-1)
        end do
      end if

      call mpi_file_open(mpi_comm_world,"./out_"//trim(iter_tag)//".dat",mpi_mode_rdonly,mpi_info_null,io,ier)

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,nps,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + intsize

      call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
      call mpi_file_read_all(io,nts,1,mpi_integer,mpi_status_ignore,ier)
      disp = disp + intsize

      allocate(pv(variable%getnpv(),-1:ini%imax+3,-1:ini%jmax+3,-1:ini%kmax+3))
      allocate(tv(variable%getntv(),-1:ini%imax+3,-1:ini%jmax+3,-1:ini%kmax+3))
      allocate(qq(variable%getnpv(),variable%getnqq(),2:ini%imax,2:ini%jmax,2:ini%kmax))

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = variable%getnpv()*(ini%imax+5)*(ini%jmax+5)*(ini%kmax+5)
      call mpi_file_read_all(io,pv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize*num

      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = variable%getntv()*(ini%imax+5)*(ini%jmax+5)*(ini%kmax+5)
      call mpi_file_read_all(io,tv,num,mpi_real8,mpi_status_ignore,ier)
      disp = disp + realsize*num

#ifdef steadyini

#else
      call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
      num = variable%getnqq()*variable%getnpv()*(ini%imax-1)*(ini%jmax-1)*(ini%kmax-1)
      call mpi_file_read_all(io,qq,num,mpi_real8,mpi_status_ignore,ier)
#endif
      call mpi_file_close(io,ier)

      do k=2,ini%kmax
        do j=2,ini%jmax
          do i=2,ini%imax
            do n=1,variable%getnpv()
              call variable%setpv(n,i,j,k,pv(n,i,j,k))
            end do
        
            do n=1,variable%getntv()
              call variable%settv(n,i,j,k,tv(n,i,j,k))
            end do      
            
            call eos%deteos(pv(1,i,j,k)+ini%pref,pv(5,i,j,k),pv(6,i,j,k),pv(7,i,j,k),dv)
            
            do n=1,variable%getndv()
              call variable%setdv(n,i,j,k,dv(n))
            end do
        
#ifdef steadyini
              qq_temp(1) = dv(1)
              qq_temp(2) = dv(1)*pv(2,i,j,k)
              qq_temp(3) = dv(1)*pv(3,i,j,k)
              qq_temp(4) = dv(1)*pv(4,i,j,k)
              qq_temp(5) = dv(1)*(dv(2)+0.5d0*(pv(2,i,j,k)**2+pv(3,i,j,k)**2+pv(4,i,j,k)**2))-pv(1,i,j,k)-ini%pref
              do n=6,variable%getnpv()
                qq_temp(n) = dv(1)*pv(n,i,j,k)
              end do

              call variable%setqq(1,i,j,k,qq_temp)
              call variable%setqq(2,i,j,k,qq_temp)
#else
            do n=1,variable%getnqq()
              qq_temp = qq(:,n,i,j,k)
              call variable%setqq(n,i,j,k,qq_temp)
            end do
#endif
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

#ifdef steadyini
      nps = 1
#endif
      
    end subroutine restart
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef test

#elif lax3d
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

#elif shocktube
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
#else
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
            call variable%setpv(2,i,j,k,ini%uref*dcos(ini%aos)*dcos(ini%aoa))
            call variable%setpv(3,i,j,k,ini%uref*dcos(ini%aos)*dsin(ini%aoa))
            call variable%setpv(4,i,j,k,ini%uref*dsin(ini%aos))
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
#endif
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module initial_module

