module datawriting_module
  use config_module
  use postgrid_module
  use postvariable_module
  implicit none
#include 'cgnslib_f.h'
#include 'cgnstypes_f.h'
  private
  public :: t_datawriting
  
  type t_datawriting
  
    contains
      procedure :: cgnswriting
      procedure :: clcd_writing
  end type t_datawriting
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine cgnswriting(datawriting,config,variable)
      implicit none
      class(t_datawriting), intent(inout) :: datawriting
      type(t_config), intent(in) :: config
      type(t_variable), intent(in) :: variable
      integer :: ifile,ier,index_flow,index_field
      integer :: n,m,l,i,j,k
      real(8), dimension(:), allocatable :: time
      real(8), dimension(:,:,:,:), allocatable :: pv,dv,tv
      real(8), dimension(:,:,:), allocatable :: a
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
        do m=0,variable%getsize(n)-1
          write(*,*) 'writing variables to domain',m+1
          allocate(pv(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1,variable%getnpv()))
          allocate(dv(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1,variable%getndv()))
          allocate(tv(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1,variable%getntv()))
          allocate(a(variable%getimax(n,m)-1,variable%getjmax(n,m)-1,variable%getkmax(n,m)-1))
          
          do k=2,variable%getkmax(n,m)
            do j=2,variable%getjmax(n,m)
              do i=2,variable%getimax(n,m)
                do l=1,variable%getnpv()
                  pv(i-1,j-1,k-1,l) = variable%getpv(n,m,l,i,j,k)
                end do
                do l=1,variable%getndv()
                  dv(i-1,j-1,k-1,l) = variable%getdv(n,m,l,i,j,k)
                end do
                do l=1,variable%getntv()
                  tv(i-1,j-1,k-1,l) = variable%gettv(n,m,l,i,j,k)
                end do
              end do
            end do
          end do
          
          a = dsqrt(dv(:,:,:,6))
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
          if(config%getiturb().ge.-1) then
            call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'k',pv(:,:,:,8),index_field,ier)
            call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'omega',pv(:,:,:,9),index_field,ier)
            call cg_field_write_f(ifile,1,variable%getrank(n,m)+1,index_flow,realdouble,'emut',tv(:,:,:,3),index_field,ier)
          end if
          deallocate(pv,dv,tv,a)
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
    subroutine clcd_writing(datawriting,config,grid,variable)
      implicit none
      class(t_datawriting), intent(inout) :: datawriting
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid
      type(t_variable), intent(in) :: variable
      integer :: io,n,m,i,j
      real(8) :: cx(2),cl1,cl2,cl
      
      !open(newunit=io,file='./cl_'//trim(config%getname())//'.dat',status='unknown',action='write',form='formatted')
      !write(io,*) 'variables = "time", "cl"'
      !write(io,*) 'zone t = " ",i=',variable%getnsolution()
      !do n=1,variable%getnsolution()
      !  do m=1,grid%getnbc()
      !    if((trim(grid%getbcname(m)).eq.'bcwallinviscid').or.(trim(grid%getbcname(m)).eq.'bcwallviscous')) then
      !      cl1 = 0.d0
      !      cl2 = 0.d0
      !      if(grid%getbcistart(m,2).eq.grid%getbciend(m,2)) then
      !        if(grid%getbcistart(m,2).eq.1) then !jmin
      !          do j=grid%getbcistart(m,2),grid%getbciend(m,2)
      !            do i=grid%getbcistart(m,1),grid%getbciend(m,1)
      !              cx = grid%getex(i,j)
      !              cl1 = cl1 + variable%getpv(n,0,1,i,j+1)*cx(1)
      !              cl2 = cl2 + variable%getpv(n,0,1,i,j+1)*cx(2)
      !            end do
      !          end do
      !        else                                !jmax
      !          do j=grid%getbcistart(m,2),grid%getbciend(m,2)
      !            do i=grid%getbcistart(m,1),grid%getbciend(m,1)
      !              cx = - grid%getex(i,j-1)
      !              cl1 = cl1 + variable%getpv(n,0,1,i,j-1)*cx(1)
      !              cl2 = cl2 + variable%getpv(n,0,1,i,j-1)*cx(2)
      !            end do
      !          end do
      !        end if
      !      else if(grid%getbcistart(m,1).eq.grid%getbciend(m,1)) then
      !        if(grid%getbcistart(m,1).eq.1) then !imin
      !          do j=grid%getbcistart(m,2),grid%getbciend(m,2)
      !            do i=grid%getbcistart(m,1),grid%getbciend(m,1)
      !              cx = grid%getcx(i,j)
      !              cl1 = cl1 + variable%getpv(n,0,1,i+1,j)*cx(1)
      !              cl2 = cl2 + variable%getpv(n,0,1,i+1,j)*cx(2)
      !            end do
      !          end do
      !        else                                !imax
      !          do j=grid%getbcistart(m,2),grid%getbciend(m,2)
      !            do i=grid%getbcistart(m,1),grid%getbciend(m,1)
      !              cx = - grid%getcx(i-1,j)
      !              cl1 = cl1 + variable%getpv(n,0,1,i-1,j)*cx(1)
      !              cl2 = cl2 + variable%getpv(n,0,1,i-1,j)*cx(2)
      !            end do
      !          end do
      !        end if
      !      end if
      !    end if
      !  end do
      !  cl1 = cl1/(0.5d0*config%getrhoref()*config%geturef()**2)/config%getl_chord()
      !  cl2 = cl2/(0.5d0*config%getrhoref()*config%geturef()**2)/config%getl_chord()
      !  cl = dsin(config%getaoa())*cl1+dcos(config%getaoa())*cl2
      !  write(io,*) n,cl
      !end do
      !close(io)
      
    end subroutine clcd_writing
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module datawriting_module
