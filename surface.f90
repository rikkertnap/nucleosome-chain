module surface

    use precision_definition
    use globals, only : LEFT, RIGHT
    use mathconst, only : pi
    use parameters, only : psiSR, psiSL , psiS ! surface potenital in reduced units for bc=cp

    implicit none

    real(dp) :: sigmaSurfL          ! surface density of acid on surface in nm^2
    real(dp) :: sigmaSurfR          ! surface density of acid on surface in nm^2
    real(dp) :: qsurfL, qsurfR

    real(dp), dimension(:), allocatable :: sigmaqSurfL         ! surface charge density on surface in nm^2
    real(dp), dimension(:), allocatable :: sigmaqSurfR         ! surface charge density on surface in nm^2
    real(dp), dimension(:), allocatable :: psiSurfL            ! surface potential
    real(dp), dimension(:), allocatable :: psiSurfR            ! surface potential

    !   different surface states

    real(dp) :: fdisS(6)            ! fraction of different surface states
    real(dp) :: KS(5)               ! experimemtal equilibruim constant
    real(dp) :: pKS(5)              ! experimental equilibruim constant pKS= -log[KS]
    real(dp) :: K0S(5)              ! intrinsic equilibruim constant
    real(dp) :: qS(6)               ! charge
    real(dp) :: cap                 ! capacitance

    !  taurine

    real(dp) :: fdisTaL(4),fdisTaR(4) ! fraction of different surface states
    real(dp) :: KTa(3)               ! experimemtal equilibruim constant
    real(dp) :: pKTa(3)              ! experimental equilibruim constant pKS= -log[KS]
    real(dp) :: K0Ta(3)              ! intrinsic equilibruim constant
    real(dp) :: qTa(4)               ! charge

    private  :: allocate_psiSurf_sigmaqSurf
    private  :: init_surface_quartz, init_surface_calcite, init_surface_taurine, init_surface_constcharge
    private  :: init_surface_constpotential, init_surface_clay
    private  :: KS, pKS, cap, KTa, pKTa

contains

        subroutine init_surface(bc,nsurf)

            character(len=2), intent(in) :: bc(2)
            integer, intent(in) :: nsurf

            call allocate_psiSurf_sigmaqSurf(nsurf)

            select case (bc(RIGHT))
                case ("qu")
                    call  init_surface_quartz()
                case ("cl")
                    call  init_surface_clay()
                case ("ca")
                    call init_surface_calcite()
                case ("ta")
                    call init_surface_taurine(RIGHT)
                case ("cc")
                    call init_surface_constcharge(RIGHT)
                case ("cp") 
                     call init_surface_constpotential(RIGHT)  
                case default
                    print*,"bc(RIGHT) does not match qu, cl, ca, ta, cc or cp"
            end select

            select case (bc(LEFT))
                case ("ta")
                    call init_surface_taurine(LEFT)
                case ("cc")
                    call init_surface_constcharge(LEFT)
                case ("cp") 
                    call init_surface_constpotential(LEFT)  
                case default
                    print*,"bc(LEFT) does not match ta, cc or cp"
            end select

        end subroutine init_surface

        function surface_charge(bc,psiSurf,side) result(sigmaqSurf)

            use volume,  only : nsurf

            real(dp), intent(in) :: psiSurf(:)
            character(len=2), intent(in) :: bc
            integer, intent(in) :: side

            real(dp) :: sigmaqSurf(nsurf)

            if(side==RIGHT) then
                select case (bc)
                    case ("qu")
                        sigmaqSurf = surface_charge_quartz(psiSurf)
                    case ("cl")
                        sigmaqSurf = surface_charge_clay(psiSurf)
                    case ("ca")
                        sigmaqSurf = surface_charge_calcite(psiSurf)
                    case ("ta")
                        sigmaqSurf = surface_charge_taurine(psiSurf,RIGHT)
                    case ("cc")
                        sigmaqSurf = sigmaSurfR
                    case ("cp") 
                        sigmaqSurf = surface_charge_constant_potential(RIGHT)   
                    case default
                        print*,"RIGHT: bc does not match qu, cl, ca, ta, or cc"    
                        sigmaqSurf = 0.0_dp
                end select
            elseif (side==LEFT) then
                select case (bc)
                    case ("ta")
                        sigmaqSurf = surface_charge_taurine(psiSurf,LEFT)
                    case ("cc")
                        sigmaqSurf = sigmaSurfL 
                    case ("cp") 
                        sigmaqSurf = surface_charge_constant_potential(LEFT)  
                    case default
                        print*,"LEFT: bc does not match ta, cc or cp "
                        sigmaqSurf = 0.0_dp
                end select
            else
                print*,"Error : side value is not equal to  LEFT or RIGHT"
                stop
            endif

        end function surface_charge


        subroutine allocate_psiSurf_sigmaqSurf(nsurf)

            integer, intent(in) :: nsurf

            allocate(sigmaqSurfL(nsurf))
            allocate(sigmaqSurfR(nsurf))
            allocate(psiSurfL(nsurf))
            allocate(psiSurfR(nsurf))

        end subroutine allocate_psiSurf_sigmaqSurf


        subroutine init_surface_quartz()

            use physconst, only : Na
            use parameters,  only : vsol,delta,lb

            integer :: i

            cap=1.0e-18_dp ! in F/nm^2 1F= 1C/V

            pKS(1)=   8.1_dp  !  8.1_dp !  >SOH <=> >SO- + H+ ! see Luetzenkirchen book chapter 14
            pKS(2)=  -4.4_dp  ! -4.4_dp !  >SOH2+ <=> >SOH + H+
            pKS(3)=  -1.5_dp  ! >SONa <=> >SO- + Na+
            pKS(4)=  -5.20_dp ! >SOCa+ <=> >SO- + Ca2+  Ber. Bunsenges. Phys. Chem., 98, pp 1062-1067, 1994
            pKS(5)=  -1.5_dp  ! >SOH2Cl <=>>SOH2+ Cl-

            do i=1,5
                KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid
                K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant
            enddo

            ! charges surface states
            qS(1)=-1.0_dp ! SO^-
            qS(2)=0.0_dp  ! SOH
            qS(3)=1.0_dp  ! SOH_2^+
            qS(4)=0.0_dp  ! SONa
            qS(5)=1.0_dp  ! SOCa
            qS(6)=0.0_dp  ! SOH2Cl

            ! sites density
            sigmaSurfR=8.0_dp
            sigmaSurfR = sigmaSurfR * (4.0_dp*pi*lb)*delta ! dimensionless surface charge

        end subroutine init_surface_quartz

        subroutine init_surface_calcite()

            use mathconst
            use physconst, only : Na
            use parameters,  only : vsol,delta,lb

            integer :: i

        !    cap=1.0d-18 ! in F/nm^2 1F= 1C/V
            !   1=CO3-, 2=CO3H, 3=CO3Ca+, 4=CaO^-, 5=CaOH, 6= CaOH2^+

            pKS(1)=   4.9_dp   ! >CO3H         <=> >CO3^- + H+
            pKS(2)=   2.8_dp   ! >CO3H + Ca^2+ <=> >CO3Ca^+ + H+
            pKS(3)=   12.2_dp  ! >CaOH2+       <=> >CaOH + H+
            pKS(4)=   17.0_dp  ! >CaOH         <=> >CaO- + H+

            do i=1,4
               KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid
               K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant
            enddo
            ! i=2 is different
            K0S(2)=KS(2)

            ! charges surface states
            qS(1)=-1.0_dp ! >CO3^-
            qS(2)=0.0_dp  ! >CO3H
            qS(3)=1.0_dp  ! >CO3Ca^+

            qS(4)=-1.0_dp ! >CaO^-
            qS(5)=0.0_dp  ! >CaOH
            qS(6)=1.0_dp  ! >CaOH2^+

            ! site density
            sigmaSurfR=5.0_dp
            sigmaSurfR = sigmaSurfR * (4.0_dp*pi*lb)*delta ! dimensionless surface charge

        end subroutine init_surface_calcite


        subroutine init_surface_clay()

            use physconst, only : Na
            use parameters,  only : vsol,delta,lb

            integer :: i

            ! see Luetzenkirchen book chapter 7 pages 204-206

            cap=0.90e-18_dp ! in F/nm^2 1F= 1C/V

            pKS(1)=   10.0_dp  ! 10.0  >SOH^0.5   <=> >SO^-0.5 + H+
            pKS(2)=   5.3_dp   !  5.3  >SOH_2^0.5 <=> >SOH_1.5^0 + 1/2H+
            pKS(3)=   0.20_dp  !       >SONa^0.5  <=> >SO^-0.5 + Na+
            pKS(4)=  -5.20_dp  !       >SOCa^1.5  <=> >SO^-0.5 + Ca2+
            pKS(5)=  -0.40_dp !        >SOH_2Cl^-0.5 <=> >SOH_2^0.5+Cl-

            do i=1,5
               KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid
               K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant
            enddo
            ! i=2
            !K0S(2) = KS(2)*((vsol*Na/1.0d24)**0.5) ! intrinstic equilibruim constant

            ! charges surface states
            qS(1)=-0.5_dp  ! SOH^-0.5
            qS(2)= 0.5_dp  ! SOH_2^0.5
            qS(3)= 0.0_dp  ! SOH_1.5^0
            qS(4)= 0.5_dp  ! SOHNa^0.5
            qS(5)= 1.5_dp  ! SOHCa^1.5
            qS(6)=-0.5_dp  ! SOH_2Cl^-0.5

            ! site density
            sigmaSurfR = 8.5_dp
            sigmaSurfR = sigmaSurfR * (4.0_dp*pi*lb)*delta ! dimensionless surface charge

        end subroutine init_surface_clay

        subroutine init_surface_taurine(side)

            use physconst, only : Na
            use parameters,  only : vsol,delta,lb

            integer, intent(in) :: side

            integer :: i

            pKTa(1)=   -2.0_dp    !  >SOH <=> >SO^- + H+
            pKTa(2)=   -0.42_dp   !  >SONa <=> >SO^- + Na+
            pKTa(3)=   -0.72243_dp  !  >SOCa^+ <=> >SO^- +Ca2+

            do i=1,3
                KTa(i)  = 10.0_dp**(-pKTa(i))       ! experimental equilibruim constant surface acid
                K0Ta(i) = (KTa(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant
            enddo

            ! charges surface states
            qTa(1)=-1.0_dp ! SO^-
            qTa(2)=0.0_dp  ! SOH
            qTa(3)=0.0_dp  ! SONa
            qTa(4)=1.0_dp  ! SOHCa^+


            ! site density
            if(side==RIGHT) sigmaSurfR = sigmaSurfR * 4.0_dp*pi*lb*delta ! dimensionless surface charge
            if(side==LEFT)  sigmaSurfL = sigmaSurfL * 4.0_dp*pi*lb*delta ! dimensionless surface charge

        end subroutine init_surface_taurine

        subroutine init_surface_constcharge(side)

            use parameters,  only : delta,lb

            implicit none

            integer, intent(in) :: side

            ! site density
            if(side==RIGHT) sigmaSurfR = sigmaSurfR * 4.0_dp*pi*lb*delta ! dimensionless surface charge
            if(side==LEFT)  sigmaSurfL = sigmaSurfL * 4.0_dp*pi*lb*delta ! dimensionless surface charge


        end subroutine init_surface_constcharge

        subroutine init_surface_constpotential(side)

            use globals, only : runtype
            use parameters,  only : beta_times_e  
            use volume, only : nsurf

            integer, intent(in) :: side

            integer :: s

            ! site density
            if(side==RIGHT) then 
                psiSR = psiSR * beta_times_e ! dimensionless surface potential
                do s = 1, nsurf
                    psiSurfR(s) = psiSR
                enddo
            endif     
            if(side==LEFT) then 

               if(runtype=="rangepsiL") psiSL = psiS%val  ! transfer value psiS to psiSL

                psiSL=  psiSL * beta_times_e ! dimensionless surface potential
                do s = 1, nsurf
                    psiSurfL(s) = psiSL
                enddo
            endif

        end subroutine init_surface_constpotential

         subroutine init_surface_constpotential_rangepsiL

            use globals, only : runtype
            use parameters,  only : beta_times_e 
            use volume, only : nsurf

            integer :: s 

            if(runtype=="rangepsiL") then 
                psiSL = psiS%val  ! transfer value psiS to psiSL  
                psiSL=  psiSL * beta_times_e ! dimensionless surface potential

                do s = 1, nsurf
                    psiSurfL(s) = psiSL
                enddo
            endif    

        end subroutine init_surface_constpotential_rangepsiL


        function surface_charge_quartz(psiS) result(surface_charge)

            use parameters, only : xbulk, vNa, vCl, vCa
            use volume, only : nsurf

            real(dp), intent(in) :: psiS(:)

            real(dp) :: surface_charge(nsurf)

            ! .. local variables

            real(dp) :: xS(6)
            real(dp) :: A,avfdis
            integer :: i,s

            do s=1,nsurf
                xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*exp(-psiS(s)) ! SOH/SO-
                xS(2)= ((xbulk%Hplus/xbulk%sol)/K0S(2))*exp(-psiS(s)) ! SOH2+/SOH
                xS(3)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0S(3))*exp(-psiS(s)) ! SONa/SO-
                xS(4)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0S(4))*exp(-2.0_dp*psiS(s)) ! SOCa+/SO-
                xS(5)= (((xbulk%Cl/vCl)/(xbulk%sol**(vCl)))/K0S(5))*exp(psiS(s)) ! SOH2Cl/SOH2+

                A = xS(1)*(1.0_dp+xS(2))+xS(3)+xS(4)+xS(1)*xS(2)*xS(5)
                fdisS(1)  = 1.0_dp/(1.0_dp +A) ! SO-
                fdisS(2)  = fdisS(1)*xS(1) ! SOH
                fdisS(3)  = fdisS(2)*xS(2) ! SOH2+
                fdisS(4)  = fdisS(1)*xS(3) ! SONa
                fdisS(5)  = fdisS(1)*xS(4) ! SOCa+
                fdisS(6)  = fdisS(1)*xS(5)*xS(2)*xS(1)  ! SOH2Cl

                avfdis=0.0_dp
                do i=1,6
                    avfdis=avfdis +qS(i)*fdisS(i)
                enddo

                surface_charge(s)=sigmaSurfR*avfdis
            enddo

        end function surface_charge_quartz


        function surface_charge_calcite(psiS) result(surface_charge)

            use parameters, only : xbulk, vCa, vsol
            use volume, only : nsurf

            real(dp), intent(in) :: psiS(:)
            real(dp) :: surface_charge(nsurf)

            ! .. local variables
            real(dp) :: xS(6)
            real(dp) :: A,avfdis
            integer :: i,s

            do s=1,nsurf

                xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*exp(-psiS(s)) ! CO3H/CO3- =fdisS(2)/fdisS(1)

                xS(2)= (xbulk%sol/xbulk%Hplus)*((xbulk%Ca*vsol/vCa)/(xbulk%sol**(vCa/vsol)))*K0S(2)*exp(-psiS(s))   ! CO3Ca+/CO3H = fdisS(3)/fdisS(2)
                xS(3)= ((xbulk%Hplus/xbulk%sol)/K0S(4))*exp(-psiS(s)) ! CaOH/CaO- =fdisS(5)/fdisS(4)
                xS(4)= ((xbulk%Hplus/xbulk%sol)/K0S(3))*exp(-psiS(s)) ! CaOH2^+/CaOH = fdisS(6)/fdisS(5)

                A = xS(1)*(1.0_dp+xS(2))
                fdisS(1)  = 1.0_dp/(1.0_dp +A) ! CO3^-
                fdisS(2)  = fdisS(1)*xS(1)   ! CO3H
                fdisS(3)  = fdisS(2)*xS(2)   ! CO3Ca+

                A=xS(3)*(1.0_dp+xS(4))
                fdisS(4)  = 1.0_dp/(1.0_dp +A) ! CaO^-
                fdisS(5)  = fdisS(4)*xS(3)   ! CaOH
                fdisS(6)  = fdisS(5)*xS(4)   ! CaOH2^+
                avfdis=0.0_dp
                do i=1,6
                    avfdis=avfdis +qS(i)*fdisS(i)
                enddo

                surface_charge(s)=sigmaSurfR*avfdis
            enddo

        end function surface_charge_calcite


        function surface_charge_clay(psiS) result(surface_charge)

            use parameters, only : xbulk, vNa, vCl, vCa
            use volume, only : nsurf

            real(dp), intent(in) :: psiS(:)
            real(dp) :: surface_charge(nsurf)

            ! .. local variables

            real(dp) :: xS(6)
            real(dp) :: A,avfdis
            integer :: i,s

            do s=1,nsurf

                xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*exp(-psiS(s)) ! SOH_2^0.5/SOH^-0.5
                xS(2)= (((xbulk%Hplus/xbulk%sol)**0.5_dp)/K0S(2))*exp(-0.5_dp*psiS(s)) ! SOH_1.5^0/SOH_2^0.5
                xS(3)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0S(3))*exp(-psiS(s)) ! SOHNa^0.5/SOH^-0.5
                xS(4)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0S(4))*exp(-2.0_dp*psiS(s)) ! SOHCa^1.5/SOH^-0.5
                xS(5)= (((xbulk%Cl/vCl)/(xbulk%sol**(vCl)))/K0S(5))*exp(psiS(s)) ! SOH2Cl^-0.5/SOH_2^0.5

                xS(2)=0.0 ! make zero

                A = xS(1)*(1.0_dp+xS(2))+xS(3)+xS(4)+xS(1)*xS(5)
                fdisS(1)  = 1.0_dp/(1.0_dp +A) ! SOH^-0.5
                fdisS(2)  = fdisS(1)*xS(1) ! SOH_2^0.5
                fdisS(3)  = fdisS(2)*xS(2) ! SOH_1.5^0
                fdisS(4)  = fdisS(1)*xS(3) ! SOHNa^0.5
                fdisS(5)  = fdisS(1)*xS(4) ! SOHCa^1.5
                fdisS(6)  = fdisS(1)*xS(1)*xS(5)  ! SOH_2Cl^-0.5

                avfdis=0.0_dp
                do i=1,6
                    avfdis=avfdis +qS(i)*fdisS(i)
                enddo

                surface_charge(s)=sigmaSurfR*avfdis

            enddo

        end function surface_charge_clay

        function surface_charge_taurine(psiS,side) result(surface_charge)

            use parameters, only : xbulk, vNa,  vCa
            use volume, only : nsurf

            real(dp), intent(in) :: psiS(:)
            integer, intent(in) :: side
            real(dp) :: surface_charge(nsurf)

            ! .. local variables

            real(dp) :: xS(3)
            real(dp) :: A,avfdis
            integer :: i,s

            do s=1,nsurf

                xS(1)= ((xbulk%Hplus/xbulk%sol)/K0Ta(1))*exp(-psiS(s)) ! SOH/SO^-
                xS(2)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0Ta(2))*exp(-psiS(s)) ! SONa/SO^-
                xS(3)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0Ta(3))*exp(-2.0_dp*psiS(s)) ! SOCa^+/SO^-

                A = xS(1)+xS(2)+xS(3)

                if(side==LEFT) then
                    fdisTaL(1)  = 1.0_dp/(1.0_dp + A) ! SO^-
                    fdisTaL(2)  = fdisTal(1)*xS(1) ! SOH
                    fdisTal(3)  = fdisTaL(1)*xS(2) ! SONa
                    fdisTaL(4)  = fdisTaL(1)*xS(3) ! SOCa^+

                    avfdis=0.0_dp
                    do i=1,4
                        avfdis=avfdis +qTa(i)*fdisTaL(i)
                    enddo

                    surface_charge(s)=sigmaSurfL*avfdis

                else if(side==RIGHT) then
                    fdisTaR(1)  = 1.0_dp/(1.0_dp + A) ! SO^-
                    fdisTaR(2)  = fdisTaR(1)*xS(1) ! SOH
                    fdisTaR(3)  = fdisTaR(1)*xS(2) ! SONa
                    fdisTaR(4)  = fdisTaR(1)*xS(3) ! SOCa^+

                    avfdis=0.0_dp
                    do i=1,4
                        avfdis=avfdis +qTa(i)*fdisTaR(i)
                    enddo
                    surface_charge(s)=sigmaSurfR*avfdis
                endif
            enddo

        end function surface_charge_taurine


        ! surface charge for constant surface  potential surface charge in dimensionless units

        function surface_charge_constant_potential(side) result(surface_charge)

            use globals, only : nsize
            use volume, only : nsurf, nz, nx, ny, coordtoindex
            use field, only : psi
          
            integer, intent(in) :: side
            real(dp) :: surface_charge(nsurf)

            ! .. local variables
            integer :: ix,iy, idxL, idxR, idxR2D , nshift
            
            nshift= (nsize-nx*ny)

             if(side==LEFT) then
                do ix=1,nx
                    do iy=1,ny
                        idxL = coordtoindex(ix, iy ,1)
                        surface_charge(idxL) = -2.0_dp * ( psi(idxL) - psiSL )  
                    enddo
                enddo    
            else if(side==RIGHT) then
                do ix=1,nx
                    do iy=1,ny
                        idxR = coordtoindex(ix, iy ,nz)
                        idxR2D = idxR-nshift
                        surface_charge(idxR2D) = -2.0_dp * ( psi(idxR) -psiSR ) 
                    enddo
                enddo 
            else
                print*,"Error: side value is not equal to  LEFT or RIGHT"
                stop
            endif   

        end function surface_charge_constant_potential


end module surface
