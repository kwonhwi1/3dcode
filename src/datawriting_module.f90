module datawriting_module
  use config_module
  use postgrid_module
  use postvariable_module
  implicit none
#include <cgnslib_f.h>
#include <cgnstypes_f.h>
  private
  public :: t_datawriting
  
  type t_datawriting
  
    contains
      procedure :: cgnswriting
      procedure :: clcd_writing
      procedure :: surface_writing
  end type t_datawriting
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine cgnswriting(datawriting,config,variable)
      implicit none
      class(t_datawriting), intent(inout) :: datawriting
      type(t_config), intent(in) :: config
      type(t_variable), intent(in) :: variable
      integer :: ifile,ier,index_flow,index_field
      integer :: n,m,i,j,k
      real(8), dimension(:), allocatable :: time
      real(8), dimension(:,:,:,:), allocatable :: pv,dv,tv
      real(8), dimension(:,:,:), allocatable :: a,vof
      character(7), dimension(:), allocatable :: solname
 
      allocate(solname(variable%getnsolution()),time(variable%getnsolution()))
      
      if(config%getnsteady().eq.1) then
        do n=1,variable%getnsolution()
          write(solname(n),'(i4.4)') variable%getnps(n)
          time(n) = dble(variable%getnps(n))
        end do      
      else
        do n=1,variable%getnsolution()
          write(solname(n),'(i7.7)') variable%getnts(n)
          time(n) = dble(variable%getnts(n))
        end do      
      end if

      
      call cg_open_f('./'//trim(config%getname())//'.cgns',cg_mode_read,ifile,ier)
      if(ier.ne.cg_ok) call cg_error_exit_f
      call cg_save_as_f(ifile,'./'//trim(config%getname())//'_result.cgns',cg_file_hdf5,0,ier)
      call cg_close_f(ifile,ier)
      
      call cg_open_f('./'//trim(config%getname())//'_result.cgns',cg_mode_modify,ifile,ier)
      if(ier.ne.cg_ok) call cg_error_exit_f
      
      do n=1,variable%getnsolution()
        write(*,*) 'solution=',n, 'nps=',variable%getnps(n), 'nts=',variable%getnts(n)
        do m=1,variable%getsize(n)
          write(*,*) 'writing variables to domain',m
          allocate(pv(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1,variable%getnpv()))
          allocate(dv(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1,variable%getndv()))
          allocate(tv(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1,variable%getntv()))
          allocate(a(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1))
          allocate(vof(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1))
          
          do k=2,variable%getkmax(n,m)
            do j=2,variable%getjmax(n,m)
              do i=2,variable%getimax(n,m)
                pv(i-1,j-1,k-1,:) = variable%getpv(n,m,i,j,k)
                dv(i-1,j-1,k-1,:) = variable%getdv(n,m,i,j,k)
                tv(i-1,j-1,k-1,:) = variable%gettv(n,m,i,j,k)
              end do
            end do
          end do
          
          a = dsqrt(dv(:,:,:,6))
          vof = dv(:,:,:,1)*pv(:,:,:,6)/dv(:,:,:,4)
          call cg_sol_write_f(ifile,1,variable%getrank(n,m)+1,solname(n),cellcenter,index_flow,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'pressure',pv(:,:,:,1),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'uvelocity',pv(:,:,:,2),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'vvelocity',pv(:,:,:,3),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'wvelocity',pv(:,:,:,4),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'temperature',pv(:,:,:,5),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'y1',pv(:,:,:,6),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'y2',pv(:,:,:,7),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'density',dv(:,:,:,1),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'enthalpy',dv(:,:,:,2),index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'sos',a,index_field,ier)
          call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'volumefraction',vof,index_field,ier)
          if(config%getiturb().ge.-1) then
            call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'k',pv(:,:,:,8),index_field,ier)
            call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'omega',pv(:,:,:,9),index_field,ier)
            call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'emut',tv(:,:,:,3),index_field,ier)
          end if
          deallocate(pv,dv,tv,a,vof)
        end do
        write(*,*) '-----------------------------------------------------------'
      end do
      
      call cg_biter_write_f(ifile,1,'TimeIterValues',variable%getnsolution(),ier)
      call cg_goto_f(ifile,1,ier,'BaseIterativeData_t',1,'end')
      call cg_array_write_f('TimeValues',realdouble,1,variable%getnsolution(),time,ier)
      
      call cg_close_f(ifile,ier)
      deallocate(solname,time)
    end subroutine cgnswriting
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine clcd_writing(datawriting,config,grid,variable,area)
      implicit none
      class(t_datawriting), intent(inout) :: datawriting
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      real(8), intent(in) :: area
      integer :: io,n,m,l,i,j,k
      real(8) :: cx(3),x(3)
      real(8) :: pv(variable%getnpv()),tv(variable%getntv()),grd(grid%getngrd())
      real(8) :: cl,cd
      real(8) :: ref,ut,dl,vel,vis
      real(8) :: fx_p,fy_p,fz_p,fx_v,fy_v,fz_v,dfx_p,dfy_p,dfz_p,dfx_v,dfy_v,dfz_v
      real(8), parameter :: pi=4.d0*datan(1.d0)
      
      open(newunit=io,file='./clcd_'//trim(config%getname())//'.plt',status='unknown',action='write',form='formatted')
      write(io,*) 'variables = "time", "cl","cd"'
      write(io,*) 'zone t = " ",i=',variable%getnsolution()

      ref = 2.d0/config%getrhoref()/config%geturef()**2/area

      do n=1,variable%getnsolution()
        fx_p = 0.d0
        fy_p = 0.d0
        fz_p = 0.d0
        fx_v = 0.d0
        fy_v = 0.d0
        fz_v = 0.d0
        do l=1,grid%getnzone()
          do m=1,grid%getnbc(l)
            if(trim(grid%getbcname(l,m)).eq.'BCWall') then
              do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                  do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                    if(grid%getbcistart(l,m,1).eq.grid%getbciend(l,m,1)) then
                      if(grid%getbcistart(l,m,1).eq.1) then !imin
                        cx   = grid%getcx(l,i,j,k)
                        x    = 0.25d0*(grid%getx(l,i+1,j,k)+grid%getx(l,i+1,j+1,k)+grid%getx(l,i+1,j,k+1)+grid%getx(l,i+1,j+1,k+1))
                        grd = grid%getgrd(l,i+1,j,k)
                        pv = variable%getpv(n,l,i+1,j,k)
                        tv = variable%gettv(n,l,i+1,j,k)
                      else ! imax
                        cx   = grid%getcx(l,i-1,j,k)
                        x    = 0.25d0*(grid%getx(l,i,j,k)+grid%getx(l,i,j+1,k)+grid%getx(l,i,j,k+1)+grid%getx(l,i,j+1,k+1))
                        grd = grid%getgrd(l,i-1,j,k)
                        pv = variable%getpv(n,l,i-1,j,k)
                        tv = variable%gettv(n,l,i-1,j,k)
                      end if
                    else if(grid%getbcistart(l,m,2).eq.grid%getbciend(l,m,2)) then
                      if(grid%getbcistart(l,m,2).eq.1) then !jmin
                        cx   = grid%getex(l,i,j,k)
                        x    = 0.25d0*(grid%getx(l,i,j+1,k)+grid%getx(l,i+1,j+1,k)+grid%getx(l,i,j+1,k+1)+grid%getx(l,i+1,j+1,k+1))
                        grd = grid%getgrd(l,i,j+1,k)
                        pv = variable%getpv(n,l,i,j+1,k)
                        tv = variable%gettv(n,l,i,j+1,k)
                      else ! jmax
                        cx   = grid%getex(l,i,j-1,k)
                        x    = 0.25d0*(grid%getx(l,i,j,k)+grid%getx(l,i+1,j,k)+grid%getx(l,i,j,k+1)+grid%getx(l,i+1,j,k+1))
                        grd = grid%getgrd(l,i,j-1,k)
                        pv = variable%getpv(n,l,i,j-1,k)
                        tv = variable%gettv(n,l,i,j-1,k)
                      end if
                    else if(grid%getbcistart(l,m,3).eq.grid%getbciend(l,m,3)) then
                      if(grid%getbcistart(l,m,3).eq.1) then !kmin
                        cx   = grid%gettx(l,i,j,k)
                        x    = 0.25d0*(grid%getx(l,i,j,k+1)+grid%getx(l,i+1,j,k+1)+grid%getx(l,i,j+1,k+1)+grid%getx(l,i+1,j+1,k+1))
                        grd = grid%getgrd(l,i,j,k+1)
                        pv = variable%getpv(n,l,i,j,k+1)
                        tv = variable%gettv(n,l,i,j,k+1)
                      else ! kmax
                        cx   = grid%gettx(l,i,j,k-1)
                        x    = 0.25d0*(grid%getx(l,i,j,k)+grid%getx(l,i+1,j,k)+grid%getx(l,i,j+1,k)+grid%getx(l,i+1,j+1,k))
                        grd = grid%getgrd(l,i,j,k-1)
                        pv = variable%getpv(n,l,i,j,k-1)
                        tv = variable%gettv(n,l,i,j,k-1)
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
          
        write(io,*) config%getdt_phy()*variable%getnps(n), cl, cd
      end do
      close(io)
    end subroutine clcd_writing
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine surface_writing(datawriting,config,grid,variable)
      implicit none
      class(t_datawriting), intent(inout) :: datawriting
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      integer :: io,n,m,l,i,j,k
      integer :: num1,num2,zoneorder
      real(8), dimension(:,:,:,:), allocatable :: pv,x
      
      open(newunit=io,file='./surface_'//trim(config%getname())//'.plt',status='unknown',action='write',form='formatted')
      write(io,*) 'variables = "x","y","z","p"'
      do n=1,variable%getnsolution()
        zoneorder = 0
        do l=1,grid%getnzone()
          do m=1,grid%getnbc(l)
            if(trim(grid%getbcname(l,m)).eq.'BCWall') then
              if(grid%getbcistart(l,m,1).eq.grid%getbciend(l,m,1)) then ! i-surface
                if(grid%getbcistart(l,m,1).eq.1) then !imin
                  num1 = (grid%getbciend(l,m,2)-grid%getbcistart(l,m,2)+2)
                  num2 = (grid%getbciend(l,m,3)-grid%getbcistart(l,m,3)+2)
                  zoneorder = zoneorder + 1
                  write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                  write(io,*) 'varlocation=([4]=cellcentered)'
                  write(io,*) 'zonetype=ordered, datapacking=block'
                  write(io,*) 'solutiontime=',n
                  allocate(pv(variable%getnpv(),grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                               ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                               ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  allocate(x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                              ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                              ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)+1
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)+1
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        x(:,i,j,k) = grid%getx(l,i+1,j,k)
                      end do
                    end do
                  end do
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        pv(:,i,j,k)   = variable%getpv(n,l,i+1,j,k)
                      end do
                    end do
                  end do
                  write(io,*) x(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(2,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) pv(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  deallocate(x,pv)
                else
                  num1 = (grid%getbciend(l,m,2)-grid%getbcistart(l,m,2)+2)
                  num2 = (grid%getbciend(l,m,3)-grid%getbcistart(l,m,3)+2)
                  zoneorder = zoneorder + 1
                  write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                  write(io,*) 'varlocation=([4]=cellcentered)'
                  write(io,*) 'zonetype=ordered, datapacking=block'
                  write(io,*) 'solutiontime=',n
                  allocate(pv(variable%getnpv(),grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                               ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                               ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  allocate(x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                              ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                              ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)+1
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)+1
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        x(:,i,j,k) = grid%getx(l,i,j,k)
                      end do
                    end do
                  end do
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        pv(:,i,j,k)  = variable%getpv(n,l,i-1,j,k)
                      end do
                    end do
                  end do
                  write(io,*) x(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(2,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) pv(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  deallocate(x,pv)
                end if
              else if(grid%getbcistart(l,m,2).eq.grid%getbciend(l,m,2)) then
                if(grid%getbcistart(l,m,2).eq.1) then !jmin
                  num1 = (grid%getbciend(l,m,1)-grid%getbcistart(l,m,1)+2)
                  num2 = (grid%getbciend(l,m,3)-grid%getbcistart(l,m,3)+2)
                  zoneorder = zoneorder + 1
                  write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                  write(io,*) 'varlocation=([4]=cellcentered)'
                  write(io,*) 'zonetype=ordered, datapacking=block'
                  write(io,*) 'solutiontime=',n
                  allocate(pv(variable%getnpv(),grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                               ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                               ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  allocate(x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                              ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                              ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)+1
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)+1
                        x(:,i,J,k) = grid%getx(l,i,j+1,k)
                      end do
                    end do
                  end do
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        pv(:,i,J,k)   = variable%getpv(n,l,i,j+1,k)
                      end do
                    end do
                  end do
                  write(io,*) x(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(2,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) pv(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  deallocate(x,pv)
                else
                  num1 = (grid%getbciend(l,m,1)-grid%getbcistart(l,m,1)+2)
                  num2 = (grid%getbciend(l,m,3)-grid%getbcistart(l,m,3)+2)
                  zoneorder = zoneorder + 1
                  write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                  write(io,*) 'varlocation=([4]=cellcentered)'
                  write(io,*) 'zonetype=ordered, datapacking=block'
                  write(io,*) 'solutiontime=',n
                  allocate(pv(variable%getnpv(),grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                               ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                               ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  allocate(x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                              ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                              ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1))
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)+1
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)+1
                        x(:,i,J,K) = grid%getx(l,i,j,k)
                      end do
                    end do
                  end do
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        pv(:,i,J,k)  = variable%getpv(n,l,i,j-1,k)
                      end do
                    end do
                  end do
                  write(io,*) x(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(2,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3)+1)
                  write(io,*) pv(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  deallocate(x,pv)
                end if
              else if(grid%getbcistart(l,m,3).eq.grid%getbciend(l,m,3)) then
                if(grid%getbcistart(l,m,3).eq.1) then !kmin
                  num1 = (grid%getbciend(l,m,1)-grid%getbcistart(l,m,1)+2)
                  num2 = (grid%getbciend(l,m,2)-grid%getbcistart(l,m,2)+2)
                  zoneorder = zoneorder + 1
                  write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                  write(io,*) 'varlocation=([4]=cellcentered)'
                  write(io,*) 'zonetype=ordered, datapacking=block'
                  write(io,*) 'solutiontime=',n
                  allocate(pv(variable%getnpv(),grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                               ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                               ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  allocate(x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                              ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                              ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)+1
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)+1
                        x(:,i,j,k) = grid%getx(l,i,j,k+1)
                      end do
                    end do
                  end do
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        pv(:,i,j,k)   = variable%getpv(n,l,i,j,k+1)
                      end do
                    end do
                  end do
                  write(io,*) x(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  write(io,*) x(2,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  write(io,*) x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  write(io,*) pv(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  deallocate(x,pv)
                else
                  num1 = (grid%getbciend(l,m,1)-grid%getbcistart(l,m,1)+2)
                  num2 = (grid%getbciend(l,m,2)-grid%getbcistart(l,m,2)+2)
                  zoneorder = zoneorder + 1
                  write(io,*) 'zone t = "',zoneorder,'",i=',num1,',j=',num2
                  write(io,*) 'varlocation=([4]=cellcentered)'
                  write(io,*) 'zonetype=ordered, datapacking=block'
                  write(io,*) 'solutiontime=',n
                  allocate(pv(variable%getnpv(),grid%getbcistart(l,m,1):grid%getbciend(l,m,1) &
                                               ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2) &
                                               ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  allocate(x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1 &
                              ,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1 &
                              ,grid%getbcistart(l,m,3):grid%getbciend(l,m,3)))
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)+1
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)+1
                        x(:,i,j,k) = grid%getx(l,i,j,k)
                      end do
                    end do
                  end do
                  do k=grid%getbcistart(l,m,3),grid%getbciend(l,m,3)
                    do j=grid%getbcistart(l,m,2),grid%getbciend(l,m,2)
                      do i=grid%getbcistart(l,m,1),grid%getbciend(l,m,1)
                        pv(:,i,j,k)  = variable%getpv(n,l,i,j,k-1)
                      end do
                    end do
                  end do
                  write(io,*) x(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  write(io,*) x(2,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  write(io,*) x(3,grid%getbcistart(l,m,1):grid%getbciend(l,m,1)+1,grid%getbcistart(l,m,2):grid%getbciend(l,m,2)+1,grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  write(io,*) pv(1,grid%getbcistart(l,m,1):grid%getbciend(l,m,1),grid%getbcistart(l,m,2):grid%getbciend(l,m,2),grid%getbcistart(l,m,3):grid%getbciend(l,m,3))
                  deallocate(x,pv)
                end if
              end if
            end if
          end do
        end do  
      end do
      close(io)
    end subroutine surface_writing
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module datawriting_module
