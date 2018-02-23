module postvariable_module
  use mpi
  use cgns
  use config_module
  use postgrid_module
  use eos_module
  implicit none
  private
  public :: t_variable

  type t_zone
    real(8), dimension(:,:,:,:), allocatable :: pv,tv,dv
  end type t_zone

  type t_solution
    integer :: size,nps,nts
    real(8) ::time
    type(t_zone), dimension(:), allocatable :: zone
  end type t_solution

  type t_variable
    private
    integer :: npv,ntv,ndv,nqq,nsolution
    type(t_solution), dimension(:), allocatable :: solution
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: getnps
      procedure :: getnts
      procedure :: getnsolution
      procedure :: getpv
      procedure :: gettv
      procedure :: getdv
      procedure :: cgnswriting
      procedure :: surface_writing
      procedure :: clcd_writing
  end type t_variable

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(variable,config,grid,eos,istart,iend,nsolution)
      implicit none
      class(t_variable), intent(out) :: variable
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      class(t_eos), intent(in) :: eos
      integer, intent(in) :: istart,iend,nsolution
      character(8) :: iter_tag
      integer :: i,j,k,m,l,io,iter,num,ier
      integer :: intsize,realsize
      integer(kind=mpi_offset_kind) :: disp

      call mpi_type_size(mpi_integer,intsize,ier)
      call mpi_type_size(mpi_real8,realsize,ier)

      variable%npv = config%getnpv()
      variable%ndv = config%getndv()
      variable%ntv = config%getntv()
      variable%nqq = config%getnqq()

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
        write(iter_tag,'(i8.8)') istart+iter*config%getnexport()*(l-1)
        allocate(variable%solution(l)%zone(grid%getnzone()))
        do m = 1,grid%getnzone()
          allocate(variable%solution(l)%zone(m)%pv(variable%npv,-1:grid%getimax(m)+3,-1:grid%getjmax(m)+3,-1:grid%getkmax(m)+3))
          allocate(variable%solution(l)%zone(m)%dv(variable%ndv,-1:grid%getimax(m)+3,-1:grid%getjmax(m)+3,-1:grid%getkmax(m)+3))
          allocate(variable%solution(l)%zone(m)%tv(variable%ntv,-1:grid%getimax(m)+3,-1:grid%getjmax(m)+3,-1:grid%getkmax(m)+3))
        end do

        call mpi_file_open(mpi_comm_world,"./out_"//trim(iter_tag)//".dat",mpi_mode_rdonly,mpi_info_null,io,ier)

        disp = 0
        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,variable%solution(l)%nps,1,mpi_integer,mpi_status_ignore,ier)
        disp = disp + intsize

        call mpi_file_set_view(io,disp,mpi_integer,mpi_integer,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,variable%solution(l)%nts,1,mpi_integer,mpi_status_ignore,ier)
        disp = disp + intsize

        call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
        call mpi_file_read_all(io,variable%solution(l)%time,1,mpi_real8,mpi_status_ignore,ier)
        disp = disp + realsize

        do m=1,grid%getnzone()
          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
          num = variable%npv*(grid%getimax(m)+5)*(grid%getjmax(m)+5)*(grid%getkmax(m)+5)
          call mpi_file_read_all(io,variable%solution(l)%zone(m)%pv,num,mpi_real8,mpi_status_ignore,ier)
          disp = disp + realsize*num

          call mpi_file_set_view(io,disp,mpi_real8,mpi_real8,'native',mpi_info_null,ier)
          num = variable%ntv*(grid%getimax(m)+5)*(grid%getjmax(m)+5)*(grid%getkmax(m)+5)
          call mpi_file_read_all(io,variable%solution(l)%zone(m)%tv,num,mpi_real8,mpi_status_ignore,ier)
          disp = disp + realsize*num
          num = variable%nqq*variable%npv*(grid%getimax(m)-1)*(grid%getjmax(m)-1)*(grid%getkmax(m)-1)
          disp = disp + realsize*num + 2*intsize + realsize
        end do

        call mpi_file_close(io,ier)

        do m=1,grid%getnzone()
          do k=2,grid%getkmax(m)
            do j=2,grid%getjmax(m)
              do i=2,grid%getimax(m)
                variable%solution(l)%zone(m)%pv(1,i,j,k) = variable%solution(l)%zone(m)%pv(1,i,j,k)+config%getpref()
                call eos%deteos_simple(variable%solution(l)%zone(m)%pv(1,i,j,k),variable%solution(l)%zone(m)%pv(5,i,j,k) &
                               ,variable%solution(l)%zone(m)%pv(6,i,j,k),variable%solution(l)%zone(m)%pv(7,i,j,k) &
                               ,variable%solution(l)%zone(m)%dv(:,i,j,k))
              end do
            end do
          end do
        end do
      end do

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(variable,grid)
      implicit none
      class(t_variable), intent(inout) :: variable
      type(t_grid), intent(in) :: grid
      integer :: l,m

      do l=1,variable%nsolution
        do m=1,grid%getnzone()
          deallocate(variable%solution(l)%zone(m)%pv)
          deallocate(variable%solution(l)%zone(m)%dv)
          deallocate(variable%solution(l)%zone(m)%tv)
        end do
        deallocate(variable%solution(l)%zone)
      end do

      deallocate(variable%solution)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine cgnswriting(variable,config,grid,nsolname)
      implicit none
      class(t_variable), intent(in) :: variable
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      integer, intent(in) :: nsolname
      integer :: ifile,ier,index_flow,index_field
      integer :: n,m
      integer(cgsize_t) :: dimensionvector
      real(8), dimension(:), allocatable :: time
      real(8), dimension(:,:,:), allocatable :: temp
      character(12), dimension(:), allocatable :: solname

      allocate(solname(variable%nsolution),time(variable%nsolution))

      if(nsolname.eq.1) then
        if(config%getnsteady().eq.1) then
          do n=1,variable%nsolution
            write(solname(n),'(i12.12)') variable%solution(n)%nps
            time(n) = dble(variable%solution(n)%nps)
          end do
        else if(config%getnsteady().eq.0) then
          do n=1,variable%nsolution
            write(solname(n),'(i12.12)') variable%solution(n)%nts
            time(n) = dble(variable%solution(n)%nts)
          end do
        end if
      else if(nsolname.eq.2) then
        do n=1,variable%nsolution
          write(solname(n),'(e12.6)') variable%solution(n)%time
          time(n) = variable%solution(n)%time
        end do
      end if

      call cg_open_f('./'//trim(config%getname())//'.cgns',cg_mode_read,ifile,ier)
      if(ier.ne.cg_ok) call cg_error_exit_f
      call cg_save_as_f(ifile,'./'//trim(config%getname())//'_result.cgns',cg_file_hdf5,0,ier)
      call cg_close_f(ifile,ier)

      call cg_open_f('./'//trim(config%getname())//'_result.cgns',cg_mode_modify,ifile,ier)
      if(ier.ne.cg_ok) call cg_error_exit_f

      do n=1,variable%nsolution
        write(*,*) 'solution=',n, 'nps=',variable%solution(n)%nps, 'nts=',variable%solution(n)%nts
        do m=1,grid%getnzone()
          write(*,*) 'writing variables to domain',m
          allocate(temp(grid%getimax(m)-1,grid%getjmax(m)-1,grid%getkmax(m)-1))
          call cg_sol_write_f(ifile,1,m,solname(n),cellcenter,index_flow,ier)
          temp = variable%solution(n)%zone(m)%pv(1,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'pressure',temp,index_field,ier)
          temp = variable%solution(n)%zone(m)%pv(2,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'u_velocity',temp,index_field,ier)
          temp = variable%solution(n)%zone(m)%pv(3,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'v_velocity',temp,index_field,ier)
          temp = variable%solution(n)%zone(m)%pv(4,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'w_velocity',temp,index_field,ier)
          temp = variable%solution(n)%zone(m)%pv(5,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'temperature',temp,index_field,ier)
          temp = variable%solution(n)%zone(m)%pv(6,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'y1',temp,index_field,ier)
          temp = variable%solution(n)%zone(m)%pv(7,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'y2',temp,index_field,ier)
          temp = variable%solution(n)%zone(m)%dv(1,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
          call cg_field_write_f(ifile,1,m,index_flow,realdouble,'density',temp,index_field,ier)
          if(config%getiturb().ge.-1) then
            temp = variable%solution(n)%zone(m)%pv(8,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
            call cg_field_write_f(ifile,1,m,index_flow,realdouble,'k',temp,index_field,ier)
            temp = variable%solution(n)%zone(m)%pv(9,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
            call cg_field_write_f(ifile,1,m,index_flow,realdouble,'omega',temp,index_field,ier)
            temp = variable%solution(n)%zone(m)%tv(3,2:grid%getimax(m),2:grid%getjmax(m),2:grid%getkmax(m))
            call cg_field_write_f(ifile,1,m,index_flow,realdouble,'emut',temp,index_field,ier)
          end if
          deallocate(temp)
        end do
        write(*,*) '-----------------------------------------------------------'
      end do

      call cg_biter_write_f(ifile,1,'TimeIterValues',variable%nsolution,ier)
      call cg_goto_f(ifile,1,ier,'BaseIterativeData_t',1,'end')
      dimensionvector = variable%nsolution
      call cg_array_write_f('TimeValues',realdouble,1,dimensionvector,time,ier)
      call cg_close_f(ifile,ier)
      deallocate(solname,time)
    end subroutine cgnswriting
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine surface_writing(variable,config,grid)
      implicit none
      class(t_variable), intent(in) :: variable
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      integer :: io,n,m,l,i,j,k
      integer :: num1,num2,zoneorder
      real(8) :: pv(variable%npv),xx(3)
      real(8), dimension(:,:,:), allocatable :: x,y,z,pressure

      open(newunit=io,file='./surface_'//trim(config%getname())//'.plt',status='unknown',action='write',form='formatted')
      write(io,*) 'variables = "x","y","z","p"'
      do n=1,variable%nsolution
        zoneorder = 0
        do l=1,grid%getnzone()
          do m=1,grid%getnbc(l)
            if((trim(grid%getfamname(l,m)).eq.'Surface')) then
              if(grid%getbcistart(l,m,1).eq.grid%getbciend(l,m,1)) then ! i-surface
                write(*,*) 'zone=',l,'bc=',m,'zoneorder=',zoneorder+1,'i-surface'
                num1 = (grid%getbciend(l,m,2)-grid%getbcistart(l,m,2)+2)
                num2 = (grid%getbciend(l,m,3)-grid%getbcistart(l,m,3)+2)
                zoneorder = zoneorder + 1
                write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                write(io,*) 'varlocation=([4]=cellcentered)'
                write(io,*) 'zonetype=ordered, datapacking=block'
                write(io,*) 'solutiontime=',n
                allocate(pressure(grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                 ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                 ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                allocate(x(grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                allocate(y(grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                allocate(z(grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)+1
                  do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)+1
                    do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                      if(grid%getbcistart(l,m,1).eq.1) then !imin
                        xx = grid%getx(l,i+1,j,k)
                      else
                        xx = grid%getx(l,i,j,k)
                      end if
                      x(i,j,k) = xx(1)
                      y(i,j,k) = xx(2)
                      z(i,j,k) = xx(3)
                    end do
                  end do
                end do
                do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                  do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                    do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                      if(grid%getbcistart(l,m,1).eq.1) then !imin
                        pv = variable%solution(n)%zone(l)%pv(:,i+1,j,k)
                      else
                        pv = variable%solution(n)%zone(l)%pv(:,i-1,j,k)
                      end if
                      pressure(i,j,k) = pv(1)
                    end do
                  end do
                end do
                write(io,*) x
                write(io,*) y
                write(io,*) z
                write(io,*) pressure
                deallocate(x,y,z,pressure)
              else if(grid%getbcistart(l,m,2).eq.grid%getbciend(l,m,2)) then
                write(*,*) 'zone=',l,'bc=',m,'zoneorder=',zoneorder+1,'j-surface'
                num1 = (grid%getbciend(l,m,1)-grid%getbcistart(l,m,1)+2)
                num2 = (grid%getbciend(l,m,3)-grid%getbcistart(l,m,3)+2)
                zoneorder = zoneorder + 1
                write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                write(io,*) 'varlocation=([4]=cellcentered)'
                write(io,*) 'zonetype=ordered, datapacking=block'
                write(io,*) 'solutiontime=',n
                allocate(pressure(grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                 ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                 ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                allocate(x(grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                allocate(y(grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                allocate(z(grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)+1
                  do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                    do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)+1
                      if(grid%getbcistart(l,m,2).eq.1) then !jmin
                        xx = grid%getx(l,i,j+1,k)
                      else
                        xx = grid%getx(l,i,j,k)
                      end if
                      x(i,j,k) = xx(1)
                      y(i,j,k) = xx(2)
                      z(i,j,k) = xx(3)
                    end do
                  end do
                end do
                do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                  do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                    do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                      if(grid%getbcistart(l,m,2).eq.1) then !jmin
                        pv = variable%solution(n)%zone(l)%pv(:,i,j+1,k)
                      else
                        pv = variable%solution(n)%zone(l)%pv(:,i,j-1,k)
                      end if
                      pressure(i,j,k) = pv(1)
                    end do
                  end do
                end do
                write(io,*) x
                write(io,*) y
                write(io,*) z
                write(io,*) pressure
                deallocate(x,y,z,pressure)
              else if(grid%getbcistart(l,m,3).eq.grid%getbciend(l,m,3)) then
                write(*,*) 'zone=',l,'bc=',m,'zoneorder=',zoneorder+1,'k-surface'
                num1 = (grid%getbciend(l,m,1)-grid%getbcistart(l,m,1)+2)
                num2 = (grid%getbciend(l,m,2)-grid%getbcistart(l,m,2)+2)
                zoneorder = zoneorder + 1
                write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                write(io,*) 'varlocation=([4]=cellcentered)'
                write(io,*) 'zonetype=ordered, datapacking=block'
                write(io,*) 'solutiontime=',n
                allocate(pressure(grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                 ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                 ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                allocate(x(grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                allocate(y(grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                allocate(z(grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                          ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                          ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                  do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)+1
                    do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)+1
                      if(grid%getbcistart(l,m,3).eq.1) then !kmin
                        xx = grid%getx(l,i,j,k+1)
                      else
                        xx = grid%getx(l,i,j,k)
                      end if
                      x(i,j,k) = xx(1)
                      y(i,j,k) = xx(2)
                      z(i,j,k) = xx(3)
                    end do
                  end do
                end do
                do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                  do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                    do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                      if(grid%getbcistart(l,m,3).eq.1) then !kmin
                        pv = variable%solution(n)%zone(l)%pv(:,i,j,k+1)
                      else
                        pv = variable%solution(n)%zone(l)%pv(:,i,j,k-1)
                      end if
                      pressure(i,j,k) = pv(1)
                    end do
                  end do
                end do
                write(io,*) x
                write(io,*) y
                write(io,*) z
                write(io,*) pressure
                deallocate(x,y,z,pressure)
              end if
            end if
          end do
        end do
      end do
      close(io)
    end subroutine surface_writing
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine clcd_writing(variable,config,grid,area)
      implicit none
      class(t_variable), intent(in) :: variable
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      real(8), intent(in) :: area
      integer :: io,n,m,l,i,j,k
      real(8) :: cx(3),x(3)
      real(8) :: pv(variable%npv),tv(variable%ntv),grd(grid%getngrd())
      real(8) :: cl,cd
      real(8) :: ref,ut,dl,vel,vis
      real(8) :: fx_p,fy_p,fz_p,fx_v,fy_v,fz_v,dfx_p,dfy_p,dfz_p,dfx_v,dfy_v,dfz_v
      real(8), parameter :: pi=4.d0*datan(1.d0)

      open(newunit=io,file='./clcd_'//trim(config%getname())//'.plt',status='unknown',action='write',form='formatted')
      write(io,*) 'variables = "time", "cl","cd"'
      write(io,*) 'zone t = " ",i=',variable%nsolution

      ref = 2.d0/config%getrhoref()/config%geturef()**2/area

      do n=1,variable%nsolution
        fx_p = 0.d0
        fy_p = 0.d0
        fz_p = 0.d0
        fx_v = 0.d0
        fy_v = 0.d0
        fz_v = 0.d0
        do l=1,grid%getnzone()
          do m=1,grid%getnbc(l)
            if(trim(grid%getfamname(l,m)).eq.'Surface') then
              do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                  do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                    if(grid%getbcistart(l,m,1).eq.grid%getbciend(l,m,1)) then
                      if(grid%getbcistart(l,m,1).eq.1) then !imin
                        cx   = grid%getcx(l,i,j,k)
                        x    = 0.25d0*(grid%getx(l,i+1,j,k)+grid%getx(l,i+1,j+1,k)+grid%getx(l,i+1,j,k+1)+grid%getx(l,i+1,j+1,k+1))
                        grd = grid%getgrd(l,i+1,j,k)
                        pv = variable%solution(n)%zone(l)%pv(:,i+1,j,k)
                        tv = variable%solution(n)%zone(l)%tv(:,i+1,j,k)
                      else ! imax
                        cx   = -grid%getcx(l,i-1,j,k)
                        x    = 0.25d0*(grid%getx(l,i,j,k)+grid%getx(l,i,j+1,k)+grid%getx(l,i,j,k+1)+grid%getx(l,i,j+1,k+1))
                        grd = grid%getgrd(l,i-1,j,k)
                        pv = variable%solution(n)%zone(l)%pv(:,i-1,j,k)
                        tv = variable%solution(n)%zone(l)%tv(:,i-1,j,k)
                      end if
                    else if(grid%getbcistart(l,m,2).eq.grid%getbciend(l,m,2)) then
                      if(grid%getbcistart(l,m,2).eq.1) then !jmin
                        cx   = grid%getex(l,i,j,k)
                        x    = 0.25d0*(grid%getx(l,i,j+1,k)+grid%getx(l,i+1,j+1,k)+grid%getx(l,i,j+1,k+1)+grid%getx(l,i+1,j+1,k+1))
                        grd = grid%getgrd(l,i,j+1,k)
                        pv = variable%solution(n)%zone(l)%pv(:,i,j+1,k)
                        tv = variable%solution(n)%zone(l)%tv(:,i,j+1,k)
                      else ! jmax
                        cx   = -grid%getex(l,i,j-1,k)
                        x    = 0.25d0*(grid%getx(l,i,j,k)+grid%getx(l,i+1,j,k)+grid%getx(l,i,j,k+1)+grid%getx(l,i+1,j,k+1))
                        grd = grid%getgrd(l,i,j-1,k)
                        pv = variable%solution(n)%zone(l)%pv(:,i,j-1,k)
                        tv = variable%solution(n)%zone(l)%tv(:,i,j-1,k)
                      end if
                    else if(grid%getbcistart(l,m,3).eq.grid%getbciend(l,m,3)) then
                      if(grid%getbcistart(l,m,3).eq.1) then !kmin
                        cx   = grid%gettx(l,i,j,k)
                        x    = 0.25d0*(grid%getx(l,i,j,k+1)+grid%getx(l,i+1,j,k+1)+grid%getx(l,i,j+1,k+1)+grid%getx(l,i+1,j+1,k+1))
                        grd = grid%getgrd(l,i,j,k+1)
                        pv = variable%solution(n)%zone(l)%pv(:,i,j,k+1)
                        tv = variable%solution(n)%zone(l)%tv(:,i,j,k+1)
                      else ! kmax
                        cx   = -grid%gettx(l,i,j,k-1)
                        x    = 0.25d0*(grid%getx(l,i,j,k)+grid%getx(l,i+1,j,k)+grid%getx(l,i,j+1,k)+grid%getx(l,i+1,j+1,k))
                        grd = grid%getgrd(l,i,j,k-1)
                        pv = variable%solution(n)%zone(l)%pv(:,i,j,k-1)
                        tv = variable%solution(n)%zone(l)%tv(:,i,j,k-1)
                      end if
                    end if

                    dfx_p = -pv(1)*cx(1)
                    dfy_p = -pv(1)*cx(2)
                    dfz_p = -pv(1)*cx(3)
                    dl = cx(1)**2+cx(2)**2+cx(3)**2
                    vel = cx(1)*pv(2)+cx(2)*pv(3)+cx(3)*pv(4)
                    select case(config%getiturb())
                    case(-1,0)
                      vis = tv(1)+tv(3)
                    case(-2)
                      vis = tv(1)
                    case default
                      vis = 0.d0
                    end select
                    ut = pv(2) - cx(1)*vel/dl
                    dfx_v = vis*ut/dsqrt((grd(2)-x(1))**2+(grd(3)-x(2))**2+(grd(4)-x(3))**2)
                    ut = pv(3) - cx(2)*vel/dl
                    dfy_v = vis*ut/dsqrt((grd(2)-x(1))**2+(grd(3)-x(2))**2+(grd(4)-x(3))**2)
                    ut = pv(4) - cx(3)*vel/dl
                    dfz_v = vis*ut/dsqrt((grd(2)-x(1))**2+(grd(3)-x(2))**2+(grd(4)-x(3))**2)

                    dfx_v = dfx_v*dsqrt(dl)
                    dfy_v = dfy_v*dsqrt(dl)
                    dfz_v = dfz_v*dsqrt(dl)

                    fx_p = fx_p + dfx_p
                    fy_p = fy_p + dfy_p
                    fz_p = fz_p + dfz_p
                    fx_v = fx_v + dfx_v
                    fy_v = fy_v + dfy_v
                    fz_v = fz_v + dfz_v
                  end do
                end do
              end do
            end if
          end do
        end do

        fx_p = fx_p*ref
        fy_p = fy_p*ref
        fz_p = fz_p*ref
        fx_v = fx_v*ref
        fy_v = fy_v*ref
        fz_v = fz_v*ref
        cl = (fy_p+fy_v)*dcos(config%getaoa()) - (fx_p+fx_v)*dsin(config%getaoa())
        cd = (fy_p+fy_v)*dsin(config%getaoa()) + (fx_p+fx_v)*dcos(config%getaoa())

        write(io,*) config%getdt_phy()*variable%solution(n)%nps, cl, cd
      end do
      close(io)
    end subroutine clcd_writing
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
    pure function getnsolution(variable)
      implicit none
      class(t_variable), intent(in) :: variable
      integer :: getnsolution

      getnsolution = variable%nsolution

    end function getnsolution
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getpv(variable,l,m,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m,i,j,k
      real(8) :: getpv(variable%npv)

      getpv = variable%solution(l)%zone(m)%pv(:,i,j,k)
    end function getpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function gettv(variable,l,m,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m,i,j,k
      real(8) :: gettv(variable%ntv)

      gettv = variable%solution(l)%zone(m)%tv(:,i,j,k)
    end function gettv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function getdv(variable,l,m,i,j,k)
      implicit none
      class(t_variable), intent(in) :: variable
      integer, intent(in) :: l,m,i,j,k
      real(8) :: getdv(variable%ndv)

      getdv = variable%solution(l)%zone(m)%dv(:,i,j,k)
    end function getdv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module postvariable_module
