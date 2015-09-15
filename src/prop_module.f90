module prop_module
  implicit none
  private
  public :: t_prop

  type t_iapws_coeff
    real(8), dimension(:), allocatable :: h,l
    real(8), dimension(:,:), allocatable :: hh,ll
  end type t_iapws_coeff
  
  type t_prop
    private
    procedure(p_prop_l), pointer :: prop_l
    procedure(p_prop_v), pointer :: prop_v
    procedure(p_prop_g), pointer :: prop_g
    type(t_iapws_coeff) :: prop_coeff
    contains
      procedure :: construct
      procedure, private :: set_iapws97
      procedure :: destruct
      procedure :: detprop    
  end type t_prop
 
  type t_prop2
    real(8) :: vis,cond 
  end type t_prop2

  interface
    subroutine p_prop_l(prop,rho,t,liquid)
      import t_prop
      import t_prop2
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: liquid
    end subroutine p_prop_l
    
    subroutine p_prop_v(prop,rho,t,vapor)
      import t_prop
      import t_prop2
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: vapor
    end subroutine p_prop_v
    
    subroutine p_prop_g(prop,rho,t,gas)
      import t_prop
      import t_prop2
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: gas
    end subroutine p_prop_g
  end interface

  contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine construct(prop,fluid,eostype,ngas)
      implicit none
      class(t_prop), intent(out) :: prop
      integer, intent(in) :: fluid,eostype,ngas
      
      select case(ngas)
      case(1)
        prop%prop_g => sutherland
      case(2)
        prop%prop_g => n2_v_fit
      case(3)
        prop%prop_g => o2_v_fit
      case(4)
        prop%prop_g => h2_v_fit
      case(5)
        prop%prop_g => he_v_fit
      case default
      end select
      
      select case(fluid)
      case(1) ! water
        select case(eostype)
        case(1) ! fitting
          prop%prop_l => h2o_l_fit
          prop%prop_v => h2o_v_fit
        case(2) ! iapws97
          call prop%set_iapws97()
          prop%prop_l => iapws97
          prop%prop_v => iapws97
        case(3) ! stiffened
          prop%prop_l => stiffened
          prop%prop_v => sutherland
        end select
      case(2) ! nitrogen
          prop%prop_l => n2_l_fit
          prop%prop_v => n2_v_fit
      case(3) ! oxygen
          prop%prop_l => o2_l_fit
          prop%prop_v => o2_v_fit
      case(4) ! hydrogen
          prop%prop_l => h2_l_fit
          prop%prop_v => h2_v_fit
      end select
    
    end subroutine construct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    subroutine destruct(prop)
      implicit none
      class(t_prop), intent(inout) :: prop
      
      if(associated(prop%prop_l)) nullify(prop%prop_l)
      if(associated(prop%prop_v)) nullify(prop%prop_v)
      if(associated(prop%prop_g)) nullify(prop%prop_g)

      if(allocated(prop%prop_coeff%h)) deallocate(prop%prop_coeff%h)
      if(allocated(prop%prop_coeff%l)) deallocate(prop%prop_coeff%l)
      if(allocated(prop%prop_coeff%hh)) deallocate(prop%prop_coeff%hh)
      if(allocated(prop%prop_coeff%ll)) deallocate(prop%prop_coeff%ll)
    end subroutine destruct
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine detprop(prop,rhol,rhov,rhog,t,y1,y2,tv)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rhol,rhov,rhog,t,y1,y2
      real(8), intent(out) :: tv(2)
      type(t_prop2) :: liquid,vapor,gas
      real(8) :: rho
      
      call prop%prop_l(rhol,t,liquid)
      call prop%prop_v(rhov,t,vapor)
      call prop%prop_g(rhog,t,gas)
      
      rho = 1.d0/((1.d0-y1)/rhol + y1*((1.d0-y2)/rhov + y2/rhog))      
      tv(1) = rho*(liquid%vis*(1.d0-y1)/rhol  + y1*(vapor%vis*(1.d0-y2)/rhov  &
              + gas%vis*y2/rhog))
      tv(2) = rho*(liquid%cond*(1.d0-y1)/rhol + y1*(vapor%cond*(1.d0-y2)/rhov &
              + gas%cond*y2/rhog))
      
    end subroutine detprop
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine set_iapws97(prop)
      implicit none
      class(t_prop), intent(inout) :: prop
      
      allocate(prop%prop_coeff%h(0:3),prop%prop_coeff%l(0:3))
      allocate(prop%prop_coeff%hh(0:5,0:6),prop%prop_coeff%ll(0:4,0:5))
      
      prop%prop_coeff%h(0) =   1.d0          ;    prop%prop_coeff%l(0) =   1.d0
      prop%prop_coeff%h(1) =   0.978197d0    ;    prop%prop_coeff%l(1) =   6.978267d0
      prop%prop_coeff%h(2) =   0.579829d0    ;    prop%prop_coeff%l(2) =   2.599096d0
      prop%prop_coeff%h(3) = - 0.202354d0    ;    prop%prop_coeff%l(3) = - 0.998254d0


      prop%prop_coeff%hh(0,0) =   0.5132047d0;    prop%prop_coeff%hh(1,0) =   0.3205656d0
      prop%prop_coeff%hh(0,1) =   0.2151778d0;    prop%prop_coeff%hh(1,1) =   0.7317883d0
      prop%prop_coeff%hh(0,2) = - 0.2818107d0;    prop%prop_coeff%hh(1,2) = - 1.070786d0
      prop%prop_coeff%hh(0,3) =   0.1778064d0;    prop%prop_coeff%hh(1,3) =   0.460504d0
      prop%prop_coeff%hh(0,4) = - 0.0417661d0;    prop%prop_coeff%hh(1,4) =   0.d0
      prop%prop_coeff%hh(0,5) =   0.d0       ;    prop%prop_coeff%hh(1,5) = - 0.01578386d0
      prop%prop_coeff%hh(0,6) =   0.d0       ;    prop%prop_coeff%hh(1,6) =   0.d0

      prop%prop_coeff%hh(2,0) =   0.d0       ;    prop%prop_coeff%hh(3,0) =   0.d0
      prop%prop_coeff%hh(2,1) =   1.241044d0 ;    prop%prop_coeff%hh(3,1) =   1.476783d0   
      prop%prop_coeff%hh(2,2) = - 1.263184d0 ;    prop%prop_coeff%hh(3,2) =   0.d0
      prop%prop_coeff%hh(2,3) =   0.2340379d0;    prop%prop_coeff%hh(3,3) = - 0.4924179d0
      prop%prop_coeff%hh(2,4) =   0.d0       ;    prop%prop_coeff%hh(3,4) =   0.1600435d0
      prop%prop_coeff%hh(2,5) =   0.d0       ;    prop%prop_coeff%hh(3,5) =   0.d0
      prop%prop_coeff%hh(2,6) =   0.d0       ;    prop%prop_coeff%hh(3,6) = - 0.003629481d0

      prop%prop_coeff%hh(4,0) = - 0.7782567d0;    prop%prop_coeff%hh(5,0) =   0.1885447d0
      prop%prop_coeff%hh(4,1) =   0.d0       ;    prop%prop_coeff%hh(5,1) =   0.d0
      prop%prop_coeff%hh(4,2) =   0.d0       ;    prop%prop_coeff%hh(5,2) =   0.d0
      prop%prop_coeff%hh(4,3) =   0.d0       ;    prop%prop_coeff%hh(5,3) =   0.d0
      prop%prop_coeff%hh(4,4) =   0.d0       ;    prop%prop_coeff%hh(5,4) =   0.d0
      prop%prop_coeff%hh(4,5) =   0.d0       ;    prop%prop_coeff%hh(5,5) =   0.d0
      prop%prop_coeff%hh(4,6) =   0.d0       ;    prop%prop_coeff%hh(5,6) =   0.d0

      prop%prop_coeff%ll(0,0) =   1.3293046d0  ;    prop%prop_coeff%ll(1,0) =   1.7018363d0
      prop%prop_coeff%ll(0,1) = - 0.40452437d0 ;    prop%prop_coeff%ll(1,1) = - 2.2156845d0
      prop%prop_coeff%ll(0,2) =   0.24409490d0 ;    prop%prop_coeff%ll(1,2) =   1.6511057d0
      prop%prop_coeff%ll(0,3) =   0.018660751d0;    prop%prop_coeff%ll(1,3) = - 0.76736002d0
      prop%prop_coeff%ll(0,4) = - 0.12961068d0 ;    prop%prop_coeff%ll(1,4) =   0.37283344d0
      prop%prop_coeff%ll(0,5) =   0.044809953d0;    prop%prop_coeff%ll(1,5) = - 0.11203160d0

      prop%prop_coeff%ll(2,0) =   5.2246158d0  ;    prop%prop_coeff%ll(3,0) =   8.7127675d0
      prop%prop_coeff%ll(2,1) = - 10.124111d0  ;    prop%prop_coeff%ll(3,1) = - 9.5000611d0 
      prop%prop_coeff%ll(2,2) =   4.9874687d0  ;    prop%prop_coeff%ll(3,2) =   4.3786606d0
      prop%prop_coeff%ll(2,3) = - 0.27297694d0 ;    prop%prop_coeff%ll(3,3) = - 0.91783782d0
      prop%prop_coeff%ll(2,4) = - 0.43083393d0 ;    prop%prop_coeff%ll(3,4) =   0.d0
      prop%prop_coeff%ll(2,5) =   0.13333849d0 ;    prop%prop_coeff%ll(3,5) =   0.d0

      prop%prop_coeff%ll(4,0) = - 1.8525999d0
      prop%prop_coeff%ll(4,1) =   0.93404690d0
      prop%prop_coeff%ll(4,2) =   0.d0 
      prop%prop_coeff%ll(4,3) =   0.d0  
      prop%prop_coeff%ll(4,4) =   0.d0  
      prop%prop_coeff%ll(4,5) =   0.d0      
    end subroutine set_iapws97
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine sutherland(prop,rho,t,gas)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: gas

#ifdef shocktube      
      gas%vis = 0.006d0
#else      
      gas%vis  = 1.716d-5*(t/273.11d0)**1.5d0*383.67d0/(t+110.56d0)
      !gas%cond = 0.0241d0*(t/273.d0)**1.5d0*467.d0/(t+194.d0)
#endif
      gas%cond = 1004.64d0/0.72d0*gas%vis

    end subroutine sutherland
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2o_l_fit(prop,rho,t,liquid)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: liquid
      liquid%vis = 182.29652337358428d0 - 2.138464474127563d0*t + 0.008996951839771332d0*t**2 -                                            &
            0.00001603934897568428d0*t**3 + 1.0939578759971515d-8*t**4 - 2.4506596592636383d-12*t**5 -                             &
            0.29273617960175397d0*rho+0.003279672826452764d0*t*rho - 0.000012681282250878405d0*t**2*rho +                       &
            1.8819010369235483d-8*t**3*rho - 7.174526314439514d-12*t**4*rho + 0.00011727665484754127d0*rho**2 -                 &
            1.2515260581064558d-6*t*rho**2 + 4.398491450268375d-9*t**2*rho**2 -  5.096592868004483d-12*t**3*rho**2

      liquid%cond = -2873.0920773919556d0 + 15.48896833182411d0*t + 0.023524601720546184d0*t**2 -                                          &
                0.00012515437912843214d0*t**3 + 1.2571696075836983d-7*t**4 - 3.65506256595383d-11*t**5 +                           &
                7.759665507709882d0*rho - 0.049128575886152964d0*t*rho + 1.253874327228879d-6*t**2*rho +                        &
                1.1524510867042747d-7*t**3*rho - 7.13949444766979d-11*t**4*rho - 0.0064435874488161776d0*rho**2 +               &
                0.00004365922518442514d0*t*rho**2 - 1.8922275927041228d-8*t**2*rho**2 -                                          &
                2.2527733573654497d-11*t**3*rho**2 + 1.6453197045073728d-6*rho**3 -                                              &
                1.1480302962074449d-8*t*rho**3 + 3.8673838876035004d-12*t**2*rho**3
    end subroutine h2o_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2o_v_fit(prop,rho,t,vapor)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: vapor
      
      if(t.le.347.d0) then
        
    
        vapor%vis = 0.00006818449331465641d0 - 8.671616545015754d-7*t + 4.88173425808012d-9*t**2 -                                      &
                1.3599835589740756d-11*t**3 + 1.9379102076160456d-14*t**4 - 1.1208856898634452d-17*t**5 -                          &
                4.622775716609354d-6*rho + 3.38168423217193d-8*t*rho - 8.953211921704184d-11*t**2*rho +                         &
                8.761339390837146d-14*t**3*rho - 1.1348488357421367d-17*t**4*rho - 8.226703496364258d-6*rho**2 +                &
                7.332470544669567d-8*t*rho**2 - 2.1755598606357746d-10*t**2*rho**2 +                                             &
                2.1494670148384862d-13*t**3*rho**2       
    
        vapor%cond = 0.11886575771113059d0 - 0.0015105562801725884d0*t + 8.505634478148256d-6*t**2 -                                    &
                2.366884745998674d-8*t**3 + 3.3839656098631864d-11*t**4 - 1.9620891902741127d-14*t**5 +                            &
                0.14620233232119875d0*rho - 0.0014524861710255087d0*t*rho + 5.7372705663111265d-6*t**2*rho -                    &
                1.044171880117295d-8*t**3*rho + 7.295109646462301d-12*t**4*rho + 0.02885120393186474d0*rho**2 -                 &
                0.00023252426657073257d0*t*rho**2 + 6.340942395897939d-7*t**2*rho**2 -                                           &
                5.834854288655945d-10*t**3*rho**2
    
      else if((t.gt.347.d0).and.(t.le.497.d0)) then
        
        vapor%vis = 0.0000107258705777731d0 - 7.21490280880739d-08*t + 4.06640374383612d-10*t**2 - 8.23792620042554d-13*t**3 +          &
               9.25447858609362d-16*t**4 - 4.44448240737233d-19*t**5 - 9.08279871583265d-08*rho - 5.61933547965227d-09*t*rho +   &
               3.20992571868845d-11*t**2*rho - 6.17485682146896d-14*t**3*rho + 4.06271012748734d-17*t**4*rho -                  &
               7.04507831064784d-09*rho**2 + 7.506004218597d-11*t*rho**2 - 2.03592074267208d-13*t**2*rho**2 +                   &
               1.68256110210731d-16*t**3*rho**2 - 2.8214178893372d-10*rho**3 + 1.13313589486255d-12*t*rho**3 -                  &
               1.14799089707064d-15*t**2*rho**3
        
        vapor%cond = 0.037260605103079d0 - 0.000309862515077041d0*t + 1.42237484311934d-06*t**2 - 2.72685354163864d-09*t**3 +           &
                2.81019308084187d-12*t**4 - 1.18869002456954d-15*t**5 + 0.0497554818439424d0*rho - 0.000344355978541611d0*t*rho+ &
                9.4617453819575d-07*t**2*rho - 1.20241981740753d-09*t**3*rho + 5.89956268384747d-13*t**4*rho +                  &
                0.0029213309201406d0*rho**2 - 0.0000172026461736515d0*t*rho**2 + 3.40364245522528d-08*t**2*rho**2 -             &
                2.25906653399214d-11*t**3*rho**2 + 0.0000132994306209143d0*rho**3 - 5.37175064347615d-08*t*rho**3 +             &
                5.4230558405697d-11*t**2*rho**3 + 3.30349589454178d-08*rho**4 - 6.13311964470895d-11*t*rho**4
        
      else !if((t.gt.497.d0).and.(t.le.552.d0)) then
        
        vapor%vis = 0.0000023502261679263d0 + 1.22881682908457d-08*t + 4.88348411442782d-11*t**2 - 2.72417906025321d-14*t**3 -          &
               9.05360803629602d-07*rho + 3.68189722109223d-09*t*rho - 5.35173691942562d-12*t**2*rho +                          &
               2.76825488024515d-15*t**3*rho + 1.09511339027305d-08*rho**2 - 4.86043688227258d-11*t*rho**2 +                    &
               7.67768900804663d-14*t**2*rho**2 - 4.17985519045525d-17*t**3*rho**2 - 3.14843417460143d-11*rho**3 +              &
               9.35748491847235d-14*t*rho**3 - 7.21198815398186d-17*t**2*rho**3 + 4.91736078150471d-14*rho**4 - &
               7.42494028104876d-17*t*rho**4 - 2.45596742396709d-17*rho**5
        vapor%cond =  -0.035413687557659d0 + 0.000243202319273188d0*t - 3.24154208932776d-07*t**2 + 2.43230121531957d-10*t**3 +         &
                  0.0309107955274621d0*rho - 0.000155425773741589d0*t*rho + 2.69487567888641d-07*t**2*rho -                     &
                  1.59016538638875d-10*t**3*rho + 0.000253858599042422d0*rho**2 - 9.1317389636681d-07*t*rho**2 +                &
                  8.30400668574453d-10*t**2*rho**2 + 5.64534018741425d-07*rho**3 - 1.05522723444757d-09*t*rho**3 + &
                  9.89369240215337d-10*rho**4
        
      end if
    end subroutine h2o_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine iapws97(prop,rho,t,liquid)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: liquid
      
      real(8) :: tt,rr,mu0,mu1,mu2,lam0,lam1,lam2,tt1
      integer :: i,j
      
      tt = t/647.27d0
      rr = dmax1(1.d-3,dmin1(1332.4d0,rho))/317.763d0
      tt1 = 647.27d0/t
      
      mu0 = dsqrt(tt)/(prop%prop_coeff%h(0) + prop%prop_coeff%h(1)*tt1 + prop%prop_coeff%h(2)*tt1**2 &
                     + prop%prop_coeff%h(3)*tt1**3)

      mu1 = 0.d0
      do i = 0,5
        do j = 0,6
          mu1 = mu1 + prop%prop_coeff%hh(i,j)*(tt1-1.d0)**i*(rr-1.d0)**j
        end do
      end do
      mu1 = dexp(rr*mu1)
      mu2 = 1.d0

      liquid%vis = 55.071d-6*mu0*mu1*mu2


      lam0 = dsqrt(tt)/(prop%prop_coeff%l(0) + prop%prop_coeff%l(1)*tt1 + prop%prop_coeff%l(2)*tt1**2 &
                      + prop%prop_coeff%l(3)*tt1**3)
      lam1 = 0.d0
      do i = 0,4
        do j = 0,5
          lam1 = lam1 + prop%prop_coeff%ll(i,j)*(tt1-1.d0)**i*(rr-1.d0)**j
        end do
      end do
      lam1 = dexp(rr*lam1)
      lam2 = 0.d0

      liquid%cond = 0.4945d0*(lam0*lam1 + lam2)      
      
      
    end subroutine iapws97
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine stiffened(prop,rho,t,liquid)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: liquid
      real(8) :: temp
      temp = dmax1(273.d0,t)
      liquid%vis  = 1.788d-3*dexp(-1.704d0 - 5.306d0*273.d0/temp + 7.003d0*(273.d0/temp)**2)
      liquid%cond = -0.000009438d0*t**2 + 0.007294d0*t - 0.7284d0
      
    end subroutine stiffened
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine n2_l_fit(prop,rhol,t,liquid)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rhol,t
      type(t_prop2), intent(out) :: liquid
 
      liquid%vis = -0.008430924153865833d0 - 0.000059966816988724016d0*t + 7.720708551957218d-6*t**2 -                                     &
            1.1967414335263782d-7*t**3 + 7.212269295432329d-10*t**4 - 1.5398097673111893d-12*t**5 +                                &
            0.00001859023644598734d0*rhol - 4.1209771470495197d-7*t*rhol + 2.2281552629578468d-9*t**2*rhol +                       &
            7.446355834144507d-12*t**3*rhol - 6.751810922972822d-14*t**4*rhol

        liquid%cond = -0.07034842054985035d0 - 0.005396846276546431d0*t + 0.000024029280785551317d0*t**2 +                                   &
                4.327926033229066d-7*t**3 - 1.970851600276364d-9*t**4 + 0.00030918782023753755d0*rhol +                            &
                5.749332339516528d-6*t*rhol - 4.169594059035202d-8*t**2*rhol - 1.274901660903033d-10*t**3*rhol       
    
    end subroutine n2_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine n2_v_fit(prop,rhov,t,vapor)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rhov,t
      type(t_prop2), intent(out) :: vapor

      vapor%vis = - 2.491d-7 +  7.319d-8*t + 3.483d-8*rhov + 2.299d-11*t**2 - 3.841d-10*t*rhov - 4d-13*t**3 + 2d-12*t**2*rhov

      vapor%cond = - 0.0009171d0 + 0.0001027d0*t + 5.623d-5*rhov + 3.45d-8*t**2 - 6.263d-7*t*rhov + 1.846d-6*rhov**2 -                    &
                4.202d-10*t**3+ 3.498d-9*t**2*rhov - 1.384d-8*t*rhov**2 + 1.899d-9*rhov**3
    
    end subroutine n2_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine o2_l_fit(prop,rhol,t,liquid)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rhol,t
      type(t_prop2), intent(out) :: liquid
      
      if (t.lt.70.d0) then        
        liquid%vis = 0.29727857163984944d0 - 0.0057547623155158115d0*t + 0.000023388284145741465d0*t**2 - 3.561741196149302d-8*t**3 -    &
                0.0004068453297973751d0*rhol + 6.765695641051581d-6*t*rhol - 1.2909910533677718d-8*t**2*rhol +                           &
                1.4285662220304494d-7*rhol**2 - 2.048110296670849d-9*t*rhol**2       
      else if ((t.ge.70.d0).and.(t.lt.85.d0)) then        
        liquid%vis = -0.1919477614404794d0 + 0.003934289955937033d0*t - 0.000018942400056368772d0*t**2 + 2.4109301435350475d-8*t**3 +    &
                0.00023364746564737008d0*rhol - 4.117254988282897d-6*t*rhol + 1.116979910579062d-8*t**2*rhol -                           &
                6.860174133713692d-8*rhol**2 + 9.91743209937825d-10*t*rhol**2       
      else if ((t.ge.85.d0).and.(t.lt.100.d0)) then           
        liquid%vis = 0.09025312394585623d0 - 0.0012785295990980342d0*t + 5.0094649861167464d-6*t**2 - 6.644225362343381d-9*t**3 -        &
                0.00012099184525686051d0*rhol + 1.3661658754234605d-6*t*rhol - 2.665076188321893d-9*t**2*rhol +                          &
                4.1744634881922714d-8*rhol**2 - 3.642536340672785d-10*t*rhol**2
      else        
        liquid%vis = -0.00803625174238251d0 + 0.00005922000340418791d0*t - 9.1607566002618d-8*t**2 + 8.973691940405286d-6*rhol -         &
                3.764812682010134d-8*t*rhol - 2.100299254317111d-9*rhol**2               
      end if
      
        liquid%cond = -141.61795141437668d0 + 5.216918762737852d0*t - 0.07070371464825916d0*t**2 + 0.0004343297698581477d0*t**3 -        &
                1.1844726676852275d-6*t**4 + 1.1609098514665666d-9*t**5 + 0.1833715170832602d0*rhol -                                    &
                0.006233426502002148d0*t*rhol + 0.00007382727277089705d0*t**2*rhol - 3.531707655265091d-7*t**3*rhol +                    &
                5.383705785806565d-10*t**4*rhol - 0.000059196663474975474d0*rhol**2 + 1.8465682534397517d-6*t*rhol**2 -                  &
                1.859211130116936d-8*t**2*rhol**2 + 6.107790454008584d-11*t**3*rhol**2  
      
    end subroutine o2_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine o2_v_fit(prop,rhov,t,vapor)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rhov,t
      type(t_prop2), intent(out) :: vapor
      
      vapor%vis = -2.4489559048367235d-7 + 7.430880530995061d-8*t + 2.038151724884577d-10*t**2 - 2.400679111877841d-12*t**3 +           &
            1.0710148280227404d-14*t**4 - 1.9210876205433418d-17*t**5 - 5.728768048886477d-8*rhov +                                     &
            1.4970484254660756d-9*t*rhov - 1.5629361628716688d-11*t**2*rhov + 8.490848201139843d-14*t**3*rhov -                         &
            1.859869348874556d-16*t**4*rhov - 2.879669386032337d-10*rhov**2 + 8.619466283100214d-12*t*rhov**2 -                         &
            8.481035857148986d-14*t**2*rhov**2 + 2.7560285963956266d-16*t**3*rhov**2
      
      vapor%cond = -0.0006175685891639204d0 + 0.00007775003155090338d0*t + 4.712567885412006d-7*t**2 - 4.674466049379829d-9*t**3 +      &
                2.1611405851428288d-11*t**4 - 4.1335952497308934d-14*t**5 + 0.00015494643771375157d0*rhov -                             &
                5.3917507213964886d-6*t*rhov + 8.979156786066838d-8*t**2*rhov - 6.566507202765043d-10*t**3*rhov +                       &
                1.7806526914224917d-12*t**4*rhov + 0.0000393481687054683d0*rhov**2 - 1.1614437575686187d-6*t*rhov**2 +                  &
                1.1500442064123965d-8*t**2*rhov**2 - 3.799946476827433d-11*t**3*rhov**2 + 5.266706402406586d-7*rhov**3 -                &
                9.859380026533517d-9*t*rhov**3 + 4.6781820541173406d-11*t**2*rhov**3
      
    end subroutine o2_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2_l_fit(prop,rho,t,liquid)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: liquid
      liquid = t_prop2(0.d0,0.d0)
    end subroutine h2_l_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine h2_v_fit(prop,rho,t,vapor)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: vapor
      vapor = t_prop2(0.d0,0.d0)
    end subroutine h2_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine he_v_fit(prop,rho,t,vapor)
      implicit none
      class(t_prop), intent(in) :: prop
      real(8), intent(in) :: rho,t
      type(t_prop2), intent(out) :: vapor
      vapor = t_prop2(0.d0,0.d0)
    end subroutine he_v_fit
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  end module prop_module
