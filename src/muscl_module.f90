module muscl_module
  use config_module
  use variable_module
  implicit none
  private
  public :: t_muscl,t_tvd,t_mlp

  type, abstract :: t_muscl
    private
    integer :: stencil,npv
    real(8), pointer :: x(:,:)
    real(8) :: r(2),r1(2),r2(2),alp(2)
    procedure(p_limiter), pointer :: limiter
    contains
      procedure :: construct                                !(iturb,lim)
      procedure :: destruct
      procedure :: setpv
      procedure(p_interpolation), deferred :: interpolation !(xl(npv),xr(npv))
  end type t_muscl

  abstract interface
    subroutine p_interpolation(muscl,xl,xr)
      import t_muscl
      implicit none
      class(t_muscl), intent(inout) :: muscl
      real(8), intent(out) :: xl(muscl%npv),xr(muscl%npv)
    end subroutine p_interpolation
  end interface

  type, extends(t_muscl) :: t_tvd
    contains
      procedure :: interpolation => tvd
  end type t_tvd

  type, extends(t_muscl) :: t_mlp
    contains
      procedure :: interpolation => mlp
  end type t_mlp

  interface
    function p_limiter(muscl,n) result(phi)
      import t_muscl
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
    end function p_limiter
  end interface
  
  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine construct(muscl,config,variable)
      implicit none
      class(t_muscl), intent(out) :: muscl
      type(t_config), intent(in) :: config
      type(t_variable), intent(in) :: variable
      
      select case(config%getnlim())
      case(0)
        muscl%limiter => no_limiter
      case(1)
        muscl%limiter => minmod
      case(2)
        muscl%limiter => superbee
      case(3)
        muscl%limiter => vanleer
      case(4)
        muscl%limiter => mc
      case(5)
        muscl%limiter => vanalbada
      case(6)
        muscl%limiter => beta
      case(7)
        muscl%limiter => i_3rd
      case(8)
        muscl%limiter => i_5th
      end select
      
      muscl%stencil = config%getstencil()
      muscl%npv  = variable%getnpv()
  
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine destruct(muscl)
      implicit none
      class(t_muscl), intent(inout) :: muscl
      
      if(associated(muscl%x))       nullify(muscl%x)
      if(associated(muscl%limiter)) nullify(muscl%limiter)

    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setpv(muscl,x)
      implicit none
      class(t_muscl), intent(inout) :: muscl
      real(8), intent(in), target :: x(muscl%stencil,muscl%npv)

      muscl%x => x
    end subroutine setpv
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine tvd(muscl,xl,xr)
      implicit none
      class(t_tvd), intent(inout) :: muscl
      real(8), intent(out) :: xl(muscl%npv),xr(muscl%npv)
      integer :: k
      real(8) :: dq,dqm,dqp,dqm1,dqp1,dqmm,dqpp
      real(8), parameter :: eps = 1.d-16

      do k = 1, muscl%npv    
        dq   = muscl%x(10,k) - muscl%x(9,k)  !i+1/2
        dqm  = muscl%x(9,k)  - muscl%x(23,k) !i-1/2
        dqm1 = muscl%x(23,k) - muscl%x(37,k) !i-3/2
        dqp  = muscl%x(32,k) - muscl%x(10,k) !i+3/2
        dqp1 = muscl%x(38,k) - muscl%x(32,k) !i+5/2

        if(dabs(dqm).le.eps) then
          dqmm = 1.d0/dsign(eps,dqm)
          muscl%r(1)  = dq*dqmm
          muscl%r1(1) = dqm1*dqmm ! inverse
          muscl%r2(1) = dqp*dqmm
        else
          dqmm = 1.d0/dqm
          muscl%r(1)  = dq*dqmm
          muscl%r1(1) = dqm1*dqmm ! inverse
          muscl%r2(1) = dqp*dqmm
        end if

        if(dabs(dqp).le.eps) then
          dqpp = 1.d0/dsign(eps,dqp)
          muscl%r(2)  = dq*dqpp ! inverse
          muscl%r1(2) = dqp1*dqpp
          muscl%r2(2) = dqm*dqpp ! inverse
        else
          dqpp = 1.d0/dqp
          muscl%r(2)  = dq*dqpp ! inverse
          muscl%r1(2) = dqp1*dqpp
          muscl%r2(2) = dqm*dqpp ! inverse
        end if
      
        muscl%alp(1) = 2.d0
        muscl%alp(2) = 2.d0
      
        xl(k) = muscl%x(9,k)  + 0.5d0*muscl%limiter(1)*dqm
        xr(k) = muscl%x(10,k) - 0.5d0*muscl%limiter(2)*dqp
      end do
    end subroutine tvd
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine mlp(muscl,xl,xr)
      implicit none
      class(t_mlp), intent(inout) :: muscl
      real(8), intent(out) :: xl(muscl%npv),xr(muscl%npv)
      integer :: k
      real(8) :: dq,dqm,dqp,dqm1,dqp1,dqmm,dqpp
      real(8) :: rxy_l,rxy_r,rxz_l,rxz_r,qml,qmr,qmin,qmax
      real(8), parameter :: eps = 1.d-16, eps2 = 1.d-3

      do k = 1, muscl%npv    
        dq   = muscl%x(10,k) - muscl%x(9,k)  !i+1/2
        dqm  = muscl%x(9,k)  - muscl%x(23,k) !i-1/2
        dqm1 = muscl%x(23,k) - muscl%x(37,k) !i-3/2
        dqp  = muscl%x(32,k) - muscl%x(10,k) !i+3/2
        dqp1 = muscl%x(38,k) - muscl%x(32,k) !i+5/2

        if(dabs(dqm).le.eps) then
          dqmm = 1.d0/dsign(eps,dqm)
          muscl%r(1)  = dq*dqmm
          muscl%r1(1) = dqm1*dqmm ! inverse
          muscl%r2(1) = dqp*dqmm
        else
          dqmm = 1.d0/dqm
          muscl%r(1)  = dq*dqmm
          muscl%r1(1) = dqm1*dqmm ! inverse
          muscl%r2(1) = dqp*dqmm
        end if

        if(dabs(dqp).le.eps) then
          dqpp = 1.d0/dsign(eps,dqp)
          muscl%r(2)  = dq*dqpp ! inverse
          muscl%r1(2) = dqp1*dqpp
          muscl%r2(2) = dqm*dqpp ! inverse
        else
          dqpp = 1.d0/dqp
          muscl%r(2)  = dq*dqpp ! inverse
          muscl%r1(2) = dqp1*dqpp
          muscl%r2(2) = dqm*dqpp ! inverse
        end if

        if(dabs(muscl%x(10,k)-muscl%x(23,k)).le.eps) then
          rxy_l = dabs((muscl%x(11,k)-muscl%x(7,k))/eps)
          rxz_l = dabs((muscl%x(15,k)-muscl%x(3,k))/eps)
        else
          rxy_l = dabs((muscl%x(11,k)-muscl%x(7,k))/(muscl%x(10,k)-muscl%x(23,k)))
          rxz_l = dabs((muscl%x(15,k)-muscl%x(3,k))/(muscl%x(10,k)-muscl%x(23,k)))
        end if
        if(dabs(muscl%x(32,k)-muscl%x(9,k)).le.eps) then
          rxy_r = dabs((muscl%x(12,k)-muscl%x(8,k))/eps)
          rxz_r = dabs((muscl%x(16,k)-muscl%x(4,k))/eps)
        else      
          rxy_r = dabs((muscl%x(12,k)-muscl%x(8,k))/(muscl%x(32,k)-muscl%x(9,k)))
          rxz_r = dabs((muscl%x(16,k)-muscl%x(4,k))/(muscl%x(32,k)-muscl%x(9,k)))
        end if
        
        if((rxy_l.lt.eps2).and.(rxz_l.lt.eps2)) then
          qmin = dmin1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                       muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k), &
                       muscl%x(19,k),muscl%x(20,k),muscl%x(21,k),muscl%x(22,k),muscl%x(23,k),muscl%x(24,k),muscl%x(25,k),muscl%x(26,k),muscl%x(27,k) )
          qmax = dmax1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                       muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k), &
                       muscl%x(19,k),muscl%x(20,k),muscl%x(21,k),muscl%x(22,k),muscl%x(23,k),muscl%x(24,k),muscl%x(25,k),muscl%x(26,k),muscl%x(27,k) )
        else
          if((muscl%x(10,k)-muscl%x(23,k)).gt.0.d0) then
            qmax = dmax1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                         muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )
            qmin = dmin1(muscl%x(1,k),muscl%x(11,k),muscl%x(3,k),muscl%x(13,k),muscl%x(5,k),muscl%x(15,k),muscl%x(7,k),muscl%x(17,k),muscl%x(9,k),      &
                         muscl%x(19,k),muscl%x(20,k),muscl%x(21,k),muscl%x(22,k),muscl%x(23,k),muscl%x(24,k),muscl%x(25,k),muscl%x(26,k),muscl%x(27,k) )
          else
            qmin = dmin1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                         muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )
            qmax = dmax1(muscl%x(1,k),muscl%x(11,k),muscl%x(3,k),muscl%x(13,k),muscl%x(5,k),muscl%x(15,k),muscl%x(7,k),muscl%x(17,k),muscl%x(9,k),      &
                         muscl%x(19,k),muscl%x(20,k),muscl%x(21,k),muscl%x(22,k),muscl%x(23,k),muscl%x(24,k),muscl%x(25,k),muscl%x(26,k),muscl%x(27,k) )
          end if
        end if
        qml = dmin1(dabs(qmax-muscl%x(9,k)),dabs(qmin-muscl%x(9,k)))
        
        if((rxy_r.lt.eps2).and.(rxz_r.lt.eps2)) then
          qmin = dmin1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                       muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k), &
                       muscl%x(28,k),muscl%x(29,k),muscl%x(30,k),muscl%x(31,k),muscl%x(32,k),muscl%x(33,k),muscl%x(34,k),muscl%x(35,k),muscl%x(36,k) )
          qmax = dmax1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                       muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k), &
                       muscl%x(28,k),muscl%x(29,k),muscl%x(30,k),muscl%x(31,k),muscl%x(32,k),muscl%x(33,k),muscl%x(34,k),muscl%x(35,k),muscl%x(36,k) )  
        else
          if((muscl%x(32,k)-muscl%x(9,k)).gt.0.d0) then
            qmax = dmax1(muscl%x(10,k),muscl%x(2,k),muscl%x(12,k),muscl%x(4,k),muscl%x(14,k),muscl%x(6,k),muscl%x(16,k),muscl%x(8,k),muscl%x(18,k),     &
                         muscl%x(28,k),muscl%x(29,k),muscl%x(30,k),muscl%x(31,k),muscl%x(32,k),muscl%x(33,k),muscl%x(34,k),muscl%x(35,k),muscl%x(36,k) )
            qmin = dmin1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                         muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )                
          else
            qmin = dmin1(muscl%x(10,k),muscl%x(2,k),muscl%x(12,k),muscl%x(4,k),muscl%x(14,k),muscl%x(6,k),muscl%x(16,k),muscl%x(8,k),muscl%x(18,k),     &
                         muscl%x(28,k),muscl%x(29,k),muscl%x(30,k),muscl%x(31,k),muscl%x(32,k),muscl%x(33,k),muscl%x(34,k),muscl%x(35,k),muscl%x(36,k) )
            qmax = dmax1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
                         muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )
          end if
        end if
        qmr = dmin1(dabs(qmax-muscl%x(10,k)),dabs(qmin-muscl%x(10,k)))
        
        !if(dq.gt.0.d0) then
        !  qml = dmax1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
        !              muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )-muscl%x(9,k)
        !  qmr = dmin1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
        !              muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )-muscl%x(10,k))        
        !else
        !  qml = dmin1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
        !              muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )-muscl%x(9,k)
        !  qmr = dmax1(muscl%x(1,k),muscl%x(2,k),muscl%x(3,k),muscl%x(4,k),muscl%x(5,k),muscl%x(6,k),muscl%x(7,k),muscl%x(8,k),muscl%x(9,k),          &
        !              muscl%x(10,k),muscl%x(11,k),muscl%x(12,k),muscl%x(13,k),muscl%x(14,k),muscl%x(15,k),muscl%x(16,k),muscl%x(17,k),muscl%x(18,k) )-muscl%x(10,k))
        !end if

        if(dabs(dq).le.eps) then
          muscl%alp(1) = 2.d0*dmax1(1.d0,muscl%r(1))/eps/(1.d0+rxy_l+rxz_l)*qml
          muscl%alp(2) = 2.d0*dmax1(1.d0,muscl%r(2))/eps/(1.d0+rxy_r+rxz_r)*qmr
        else
          muscl%alp(1) = 2.d0*dmax1(1.d0,muscl%r(1))/dabs(dq)/(1.d0+rxy_l+rxz_l)*qml
          muscl%alp(2) = 2.d0*dmax1(1.d0,muscl%r(2))/dabs(dq)/(1.d0+rxy_r+rxz_r)*qmr
        end if

        muscl%alp(1) = dmax1(1.d0,dmin1(2.d0,muscl%alp(1)))
        muscl%alp(2) = dmax1(1.d0,dmin1(2.d0,muscl%alp(2)))
      
        xl(k) = muscl%x(9,k)  + 0.5d0*muscl%limiter(1)*dqm
        xr(k) = muscl%x(10,k) - 0.5d0*muscl%limiter(2)*dqp
      end do
    end subroutine mlp
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function no_limiter(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
      
      phi = 1.d0

    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function minmod(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi

      phi = dmax1(0.d0,dmin1(1.d0,muscl%r(n)))

    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function superbee(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi

      phi = dmax1(0.d0,dmin1(1.d0,muscl%alp(n)*muscl%r(n)),dmin1(muscl%alp(n),muscl%r(n))) 

    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function vanleer(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
      real(8) :: c

      c = (muscl%r(n)+dabs(muscl%r(n)))/(1.d0+dabs(muscl%r(n)))
      phi = dmax1(0.d0,dmin1(muscl%alp(n),muscl%alp(n)*muscl%r(n),c))

    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function mc(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
      real(8) :: c

      c = (1.d0+muscl%r(n))/2.d0
      phi = dmax1(0.d0,dmin1(muscl%alp(n),muscl%alp(n)*muscl%r(n),c))
        
    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function vanalbada(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
      real(8) :: c

      c =(muscl%r(n)**2+muscl%r(n))/(1.d0+muscl%r(n)**2)
      phi = dmax1(0.d0,dmin1(c*muscl%r(n),1.d0),dmin1(muscl%r(n),c))
        
    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function beta(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
      real(8) :: c

      c = dmin1(1.5d0,muscl%alp(n))
      phi = dmax1(0.d0,dmin1(c*muscl%r(n),1.d0),dmin1(muscl%r(n),c))
    
    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function i_3rd(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
      real(8) :: c

      c = (1.d0+2.d0*muscl%r(n))/3.d0
      phi = dmax1(0.d0,dmin1(muscl%alp(n),muscl%alp(n)*muscl%r(n),c))

    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function i_5th(muscl,n) result(phi)
      implicit none
      class(t_muscl), intent(in) :: muscl
      integer, intent(in) :: n
      real(8)  :: phi
      real(8) :: c
      
      c = (-2.d0*muscl%r1(n)+11.d0+24.d0*muscl%r(n)-3.d0*muscl%r2(n))/30.d0
      phi = dmax1(0.d0,dmin1(muscl%alp(n),muscl%alp(n)*muscl%r(n),c))
        
    end function
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module muscl_module
