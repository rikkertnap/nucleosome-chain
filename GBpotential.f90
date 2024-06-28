! Module to compute Gay-Berne implementation following Persson, jcp, 136, 2012
! LJ interactiopn between unaxial liquid crystals or uniaxial anisotropic atoms or ellipsoid
! Contains two potentials 
! (1) GBpotential: orginal GB pot( Jcp 74 1981) 
! (2) GBpotential_Persson : modification proposed by Persson

module GB_potential     

    use precision_definition
  
    implicit none

    real(dp) :: epsilonS        !  side-to-side binding energy
    real(dp) :: epsilonE        !  end-to-end binding energy 
    real(dp) :: epsilon0        !  unit of energy 
    real(dp) :: sigmaS          !  short axis sigma of ellipsoid  
    real(dp) :: sigmaE          !  long axis sigma of ellipsoid 
    real(dp) :: sigma0          !  unit of lenght 
    
    ! derived quanities    
    real(dp) :: axisratio       !  sigmaE/sigmaS     : l in Persson paper
    real(dp) :: energyratio     !  epsilonS/epsilonE : d in Persson paper 

    real(dp) :: chi             ! 
    real(dp) :: chiprime        !   

    character(len=15) :: GBtype   
    character(len=15) :: GBCOMtype 

    private :: dotproduct, shapefunc_chi, shapefunc_chiprime
    private :: shapefunc_sigma, shapefunc_epsilon
    private :: shapefunc_sigma_Persson, shapefunc_epsilon_Persson

contains  

    function shapefunc_chi(l)result(chi)

        real(dp), intent(in) :: l
        real(dp) :: chi

        chi=(l**2-1.0_dp)/(l**2+1.0_dp)
    
    end function shapefunc_chi

    function  shapefunc_chiprime(d) result(chiprime)

        real(dp), intent(in) :: d
        real(dp) :: chiprime 
         
        chiprime=(sqrt(d)-1.0_dp)/(sqrt(d)+1.0_dp)

    end function shapefunc_chiprime



    function shapefunc_sigma(u1,u2,runit,chi) result(sigma) 
    
        real(dp), intent(in) , dimension(3) :: u1, u2, runit 
        real(dp), intent(in) ::  chi
        real(dp) :: sigma

        ! local variables
        real(dp) :: a, b, c, sqraplusb, sqraminb

        a = dotproduct(runit, u1) 
        b = dotproduct(runit, u2) 
        c = dotproduct(u1 , u2) 
        sqraplusb = (a+b)**2 
        sqraminb = (a-b)**2
        sigma=1.0_dp/sqrt((1.0_dp-(chi/2.0_dp)*(sqraplusb/(1.0_dp+chi*c) +sqraminb/(1.0_dp-chi*c)))) 
   
    end function shapefunc_sigma
        
    function shapefunc_epsilon(u1,u2,runit,chi,chiprime) result(epsilon)

        real(dp), intent(in) , dimension(3) :: u1, u2, runit 
        real(dp), intent(in) :: chi, chiprime
        real(dp) :: epsilon

        ! local variables
        real(dp) ::  a, b, c, sqraplusb, sqraminb

        a = dotproduct(runit, u1) 
        b = dotproduct(runit, u2) 
        c = dotproduct(u1 , u2) 
        sqraplusb = (a+b)**2 
        sqraminb = (a-b)**2
   
        epsilon = 1.0_dp/sqrt(1.0_dp-(chi**2)*(c**2)) 
        epsilon = epsilon*(1.0_dp-(chiprime/2.0_dp)*( &
        sqraplusb/(1.0_dp+chiprime*c) +sqraminb/(1.0_dp-chiprime*c)))**2 
        
    end function shapefunc_epsilon

    ! Sets value of energy scale (epsilon0)  such that minimal energy of Gay-Berne potential 
    ! in side-to-side orientation is 1 in unit of  epsilonS ( binding side-to-side energy) 
    ! and  assumes length scale to be set to omegaS=omega0 
    
    function set_epsilon0(chi)result(epsilon0)
        real(dp), intent(in) :: chi
        real(dp) :: epsilon0

        epsilon0=sqrt(1.0-chi**2)
    
    end function    

    ! Computes energy of Gay-Berne potential  
    ! input real(dp),  u1(3), u2(3) : orientation vector of interaction ellipsoids
    !       real(dp),  rvec(3) :  center-to-center vector difference 
    ! output real(dp) VGB: value of Gay-Berne potential
    ! pre : call to init_GB_const() to set (axisratio, energyratio) chi, chiprime, epsilon0, and sigma0
    
    function GBpotential_LJ(u1,u2,rvec)result(VGB)
        
        real(dp), intent(in) , dimension(3) :: u1, u2, rvec
        real(dp) :: VGB

        ! local variables
        real(dp), dimension(3) :: runit
        real(dp) :: normr, sigma, epsilon, ratio

        normr=sqrt(dotproduct(rvec,rvec)) 
        runit=rvec/normr

        sigma=shapefunc_sigma(u1,u2,runit,chi) 
        epsilon=epsilon0 * shapefunc_epsilon(u1,u2,runit,chi,chiprime)
        ratio=1.0_dp/(normr/sigma0-sigma+1.0_dp) 

        VGB=4.0_dp*epsilon*(ratio**12-ratio**6) ! factor 4 is not in  equation in Persson 
    
    end function GBpotential_LJ


    ! Computes energy of Gay-Berne potential
    ! Van der Waals attraction only 
    ! input real(dp),  u1(3), u2(3) : orientation vector of interaction ellipsoids
    !       real(dp),  rvec(3) :  center-to-center vector difference 
    ! output real(dp) VGB: value of Gay-Berne potential
    ! pre : call to init_GB_const() to set (axisratio, energyratio) chi, chiprime, epsilon0, and sigma0
    
    function GBpotential_VdW(u1,u2,rvec)result(VGB)
        
        real(dp), intent(in) , dimension(3) :: u1, u2, rvec
        real(dp) :: VGB

        ! local variables
        real(dp), dimension(3) :: runit
        real(dp) :: normr, sigma, epsilon, ratio

        normr=sqrt(dotproduct(rvec,rvec)) 
        runit=rvec/normr

        sigma=shapefunc_sigma(u1,u2,runit,chi) 
        epsilon=epsilon0 * shapefunc_epsilon(u1,u2,runit,chi,chiprime)
        ratio=1.0_dp/(normr/sigma0-sigma+1.0_dp) 

        VGB=4.0_dp*epsilon*(-ratio**6) ! factor 4 is not in  equation in Persson 
    
    end function GBpotential_VdW

    function shapefunc_sigma_Persson(u1,u2,runit,l) result(sigma) 
    
        real(dp), intent(in) , dimension(3) :: u1, u2, runit 
        real(dp), intent(in) ::  l
        real(dp) :: sigma

        ! local variables
        real(dp) :: a, b, c, sqraplusb, sqraminb

        a = dotproduct(runit, u1) 
        b = dotproduct(runit, u2) 
       
        sigma=1.0_dp+((l-1.0_dp)/2.0_dp)*(abs(a)+abs(b)) 
   
    end function shapefunc_sigma_Persson

    function shapefunc_epsilon_Persson(u1,u2,runit,d) result(epsilon)

        real(dp), intent(in) , dimension(3) :: u1, u2, runit 
        real(dp), intent(in) :: d
        real(dp) :: epsilon

        ! local variables
        real(dp) ::  a, b

        a = dot_product(runit, u1) 
        b = dot_product(runit, u2) 
        
        epsilon = 1.0_dp+((1.0_dp/d-1.0_dp)/2.0_dp)*( abs(a)+abs(b) )
        
    end function shapefunc_epsilon_Persson

    ! Computes energy of Gay-Berne potential using Persson modifications
    ! input real(dp),  u1(3), u2(3) : orientation vector of interaction ellipsoids
    !       real(dp),  rvec(3) :  center-to-center vector difference 
    ! output real(dp) VGB: value of Gay-Berne potential
    ! pre : call to init_GB_Persson_const() to set (axisratio, energyratio) chi, chiprime, epsilon0, and sigma0

    function GBpotential_Persson_LJ(u1,u2,rvec)result(VGB)
        
        real(dp), intent(in) , dimension(3) :: u1, u2, rvec
        real(dp) :: VGB

        ! local variables
        real(dp), dimension(3) :: runit
        real(dp) :: normr,  sigma, epsilon, ratio

        normr=sqrt(dotproduct(rvec,rvec)) 
        runit=rvec/normr

        sigma=shapefunc_sigma_Persson(u1,u2,runit,axisratio) 
        epsilon=epsilon0 * shapefunc_epsilon_Persson(u1,u2,runit,energyratio)
        ratio=1.0_dp/(normr/sigma0-sigma+1.0_dp) 

        VGB=4.0_dp*epsilon*(ratio**12-ratio**6) ! factor 4 is not in  equation in Persson 
    
    end function GBpotential_Persson_LJ

    ! Computes energy of Gay-Berne potential using Persson modifications
    ! VdW attraction only 
    ! input real(dp),  u1(3), u2(3) : orientation vector of interaction ellipsoids
    !       real(dp),  rvec(3) :  center-to-center vector difference 
    ! output real(dp) VGB: value of Gay-Berne potential
    ! pre : call to init_GB_Persson_const() to set (axisratio, energyratio) chi, chiprime, epsilon0, and sigma0

    function GBpotential_Persson_VdW(u1,u2,rvec)result(VGB)
        
        real(dp), intent(in) , dimension(3) :: u1, u2, rvec
        real(dp) :: VGB

        ! local variables
        real(dp), dimension(3) :: runit
        real(dp) :: normr,  sigma, epsilon, ratio

        normr=sqrt(dotproduct(rvec,rvec)) 
        runit=rvec/normr

        sigma=shapefunc_sigma_Persson(u1,u2,runit,axisratio) 
        epsilon=epsilon0 * shapefunc_epsilon_Persson(u1,u2,runit,energyratio)
        ratio=1.0_dp/(normr/sigma0-sigma+1.0_dp) 

        VGB=4.0_dp*epsilon*(-ratio**6) ! factor 4 is not in  equation in Persson 
    
    end function GBpotential_Persson_VdW


    function GBpotential_general(u1,u2,rvec)result(VGB)
        
        real(dp), intent(in) , dimension(3) :: u1, u2, rvec
        real(dp) :: VGB

        select case (GBtype)
        case ("GayBerneLJ") 
            VGB = GBpotential_LJ(u1,u2,rvec)
        case ("GayBerneVdW") 
            VGB = GBpotential_VdW(u1,u2,rvec)
        case ("PerssonLJ")
            VGB = GBpotential_Persson_LJ(u1,u2,rvec)
        case ("PerssonVdW")
            VGB = GBpotential_Persson_VdW(u1,u2,rvec)   
        case default 
            print*,"Error in GBpotential_general"
            print*,"Wrong value GBtype : ", GBtype
            VGB = 0.0_dp
        end select

    end function GBpotential_general

    subroutine init_GB_const_defaults()

        epsilonS = 1.0_dp  
        epsilonE = 6.0_dp 

        sigmaS = 1.0_dp
        sigmaE = 0.5733672387811127_dp
        sigma0 = 10.289880194520919_dp

    end subroutine


    ! Sets Gay-Berne constants: axisratio, energyratio, sigma0, and epislon0

    subroutine init_GB_const()
    
        energyratio = epsilonS/epsilonE
        axisratio = sigmaE/sigmaS 
        
        sigmaS = sigma0
        sigmaE = axisratio * sigmaS

        chi=shapefunc_chi(axisratio)
        chiprime=shapefunc_chiprime(energyratio)

        epsilon0=set_epsilon0(chi) ! this set scale of VGB to units of epsilonS 

    end subroutine init_GB_const


    ! Sets Gay-Berne Persson constants : axisratio, energyratio, sigma0, and epislon0

    subroutine init_GB_Persson_const()
        
        energyratio = epsilonS/epsilonE
        axisratio = sigmaE/sigmaS 

        sigmaS = sigma0
        sigmaE = axisratio * sigmaS

        chi=shapefunc_chi(axisratio)
        chiprime=shapefunc_chiprime(energyratio)
        
        epsilon0=epsilonS ! this set scale of VGB_Persson to units of epsilonS 

    end subroutine init_GB_Persson_const


    subroutine init_GB_general_const()
        
        select case (GBtype)
        case ("GayBerneLJ","GayBerneVdW") 
            call  init_GB_const()
        case ("PerssonLJ","PerssonVdW")
            call init_GB_Persson_const()
        case default 
            print*,"Error in init_GB_general_const"
            print*,"Wrong value GBtype : ", GBtype
        end select 
       
    end subroutine init_GB_general_const

    subroutine print_GB_const()

        print*,"# version     = ",VERSION
        print*,"# axisratio   = ",axisratio
        print*,"# energyratio = ",energyratio
        print*,"# sigma0      = ",sigma0
        print*,"# epsilon0    = ",epsilon0
        print*,"# sigmaS      = ",sigmaS
        print*,"# sigmaE      = ",sigmaE 
    end subroutine


    subroutine print_GB_directions(u1,u2,runit)
        
        real(dp), intent(in) , dimension(3) :: u1, u2, runit

        print*,"# u1    = ",u1
        print*,"# u2    = ",u2
        print*,"# runit = ",runit

    end subroutine

    ! Inner product of vectors a and b with dimenstion 3
    ! similar to intrisic function dotproduct

    function dotproduct(a,b) result(dotprod)
       
        real(dp), dimension(3) :: a, b 
        real(dp) :: dotprod
      
        dotprod = sum(a*b)

    end function dotproduct


end module GB_potential


