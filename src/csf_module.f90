module csf_module
  use config_module
  use grid_module
  use eos_module
  implicit none
  private
  public :: t_csf

  type t_csf
    private
    integer :: stencil,npv,ngrd
    real(8) :: tref
    real(8), pointer :: pv(:,:),grd(:),vfg(:,:)
    real(8), pointer :: cx1(:),cx2(:),ex1(:),ex2(:),tx1(:),tx2(:)
    contains
      procedure :: construct
      procedure :: destruct
      procedure :: setpv
      procedure :: setgrd
      procedure :: setvfg
      procedure :: setnorm
      procedure :: csfsource
  end type t_csf

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(csf,config,grid)
      implicit none
      class(t_csf), intent(out) :: csf
      type(t_config), intent(in) :: config
      type(t_grid), intent(in) :: grid

      csf%stencil = config%getstencil()
      csf%npv = config%getnpv()
      csf%ngrd = grid%getngrd()

    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(csf)
      implicit none
      class(t_csf), intent(inout) :: csf

      if(associated(csf%pv))    nullify(csf%pv)
      if(associated(csf%grd))   nullify(csf%grd)
      if(associated(csf%vfg))   nullify(csf%vfg)

      if(associated(csf%cx1))   nullify(csf%cx1)
      if(associated(csf%cx2))   nullify(csf%cx2)
      if(associated(csf%ex1))   nullify(csf%ex1)
      if(associated(csf%ex2))   nullify(csf%ex2)
      if(associated(csf%tx1))   nullify(csf%tx1)
      if(associated(csf%tx2))   nullify(csf%tx2)

   end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(csf,pv)
      implicit none
      class(t_csf), intent(inout) :: csf
      real(8), intent(in), target :: pv(csf%stencil,csf%npv)

      csf%pv => pv

    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setgrd(csf,grd)
      implicit none
      class(t_csf), intent(inout) :: csf
      real(8), intent(in), target :: grd(csf%ngrd)

      csf%grd => grd

    end subroutine setgrd
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setvfg(csf,vfg)
      implicit none
      class(t_csf), intent(inout) :: csf
      real(8), intent(in), target :: vfg(7,3)

      csf%vfg => vfg

    end subroutine setvfg
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setnorm(csf,cx1,cx2,ex1,ex2,tx1,tx2)
      implicit none
      class(t_csf), intent(inout) :: csf
      real(8), intent(in), target :: cx1(3),cx2(3),ex1(3),ex2(3),tx1(3),tx2(3)

      ! cell outward vector
      csf%cx1 => cx1
      csf%cx2 => cx2
      csf%ex1 => ex1
      csf%ex2 => ex2
      csf%tx1 => tx1
      csf%tx2 => tx2

    end subroutine setnorm
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function csfsource(csf,eos)
      implicit none
      class(t_csf), intent(in) :: csf
      class(t_eos), intent(in) :: eos
      real(8) :: csfsource(csf%npv)
      real(8) :: delta,kappa
      real(8) :: sigma(7)
      real(8) :: normvfg(7,2),vfgf(6,2),normvfgf(6,2)
      integer :: l,m
      real(8) :: dl
      real(8) :: dsdx,dsdy,dsdz

      csfsource = 0.d0
      kappa = 0.d0

      delta = dsqrt(csf%vfg(7,1)**2+csf%vfg(7,2)**2+csf%vfg(7,3)**2)
      if(delta.le.2.041) return  !! could be too restrictive !!

      do l=1,7
!        sigma(l) = eos%get_sigma(csf%pv(l,5))
        sigma = 0.071684196047d0
        dl = csf%vfg(l,1)**2+csf%vfg(l,2)**2+csf%vfg(l,3)**2
        if(dl.gt.0.d0) then
          normvfg(l,1) = csf%vfg(l,1)/dsqrt(dl)
          normvfg(l,2) = csf%vfg(l,2)/dsqrt(dl)
          normvfg(l,3) = csf%vfg(l,3)/dsqrt(dl)
        else
          normvfg(l,:) = 0.d0
        end if
      end do

!! =====  interface normal force ===== !!

!      1: use normvfgf
!      B: define normvfgf by averaging normvfg of both cells
      do m=1,6
        normvfgf(m,:) = 0.5d0*(normvfg(7,:)+normvfg(m,:))
      end do

!      b: Green-Gauss gradient
      kappa = normvfgf(1,1)*csf%cx1(1) + normvfgf(1,2)*csf%cx1(2) + normvfgf(1,3)*csf%cx1(3) &
            + normvfgf(2,1)*csf%cx2(1) + normvfgf(2,2)*csf%cx2(2) + normvfgf(2,3)*csf%cx2(3) &
            + normvfgf(3,1)*csf%ex1(1) + normvfgf(3,2)*csf%ex1(2) + normvfgf(3,3)*csf%ex1(3) &
            + normvfgf(4,1)*csf%ex2(1) + normvfgf(4,2)*csf%ex2(2) + normvfgf(4,3)*csf%ex2(3) & 
            + normvfgf(5,1)*csf%tx1(1) + normvfgf(5,2)*csf%tx1(2) + normvfgf(5,3)*csf%tx1(3) &
            + normvfgf(6,1)*csf%tx2(1) + normvfgf(6,2)*csf%tx2(2) + normvfgf(6,3)*csf%tx2(3)
      kappa = -kappa/csf%grd(1)

      csfsource(2) = sigma(7)*kappa*csf%vfg(7,1)
      csfsource(3) = sigma(7)*kappa*csf%vfg(7,2)
      csfsource(4) = sigma(7)*kappa*csf%vfg(7,3)


!! =====  interface tangentil force ===== !!

!      ! dsdx = d(sigma)/dx, dsdy = d(sigma)/dy, dsdz = d(sigma)/dz
!      finite differences using metric vectors
      dsdx = 0.5d0*( (sigma(2)-sigma(1))*(csf%cx2(1)-csf%cx1(1)) &
                    +(sigma(4)-sigma(3))*(csf%ex2(1)-csf%ex1(1)) &
                    +(sigma(6)-sigma(5))*(csf%tx2(1)-csf%tx1(1)) )/csf%grd(1)
      dsdy = 0.5d0*( (sigma(2)-sigma(1))*(csf%cx2(2)-csf%cx1(2)) &
                    +(sigma(4)-sigma(3))*(csf%ex2(2)-csf%ex1(2)) &
                    +(sigma(6)-sigma(5))*(csf%tx2(2)-csf%tx1(2)) )/csf%grd(1)
      dsdz = 0.5d0*( (sigma(2)-sigma(1))*(csf%cx2(3)-csf%cx1(3)) &
                    +(sigma(4)-sigma(3))*(csf%ex2(3)-csf%ex1(3)) &
                    +(sigma(6)-sigma(5))*(csf%tx2(3)-csf%tx1(3)) )/csf%grd(1)

      csfsource(2) = csfsource(2) + ( (1.d0-normvfg(7,1)**2)*dsdx &
                     - normvfg(7,2)*normvfg(7,1)*dsdy - normvfg(7,3)*normvfg(7,1)*dsdz )*delta
      csfsource(3) = csfsource(3) + ( (1.d0-normvfg(7,2)**2)*dsdy &
                     - normvfg(7,1)*normvfg(7,2)*dsdx - normvfg(7,3)*normvfg(7,2)*dsdz )*delta
      csfsource(4) = csfsource(4) + ( (1.d0-normvfg(7,3)**2)*dsdz &
                     - normvfg(7,1)*normvfg(7,3)*dsdx - normvfg(7,2)*normvfg(7,3)*dsdy )*delta


      csfsource(5) = csfsource(2)*csf%pv(7,2) + csfsource(3)*csf%pv(7,3) + csfsource(4)*csf%pv(7,4)
      csfsource = csfsource*csf%grd(1)
    end function csfsource
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module csf_module
