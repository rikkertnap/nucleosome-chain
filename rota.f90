! ------------------------------------------------------------|
! rotates a given chains conformation                         |  
! pre: xend = input chain                                     |
! post: xendr= rotated chain                                  |
!       return= is_postive                                    |
! ------------------------------------------------------------|

module chain_rotation

    use precision_definition
  
    implicit none

contains  

!Euler rotation Z1X2Z3

function rotation(xend,xendr,nseg)result(is_positive_z)
  
    use random
    use mathconst
    
    implicit none

    integer,  intent(in)    :: nseg  
    real(dp), intent(in)    :: xend(:,:)
    real(dp), intent(inout) :: xendr(:,:)  
    logical                 :: is_positive_z  ! output 
  
    ! ..local argument 

    integer :: i
    real(dp)  :: fac,fac1,fac2,sbe,cbe,sal,cal,sga,dist
    real(dp)  :: alfa,gama,cga,a,b,c
    
    fac=rands(seed)
    fac1=rands(seed)
    fac2=rands(seed)
    alfa=fac*2.0_dp*pi
    cbe=fac1*2.0_dp-1.0_dp
    gama=fac2*2.0_dp*pi
  
    sbe=(1.0_dp-cbe**2)**0.5_dp
    cal=cos(alfa)
    sal=sin(alfa)
    cga=cos(gama)
    sga=sin(gama)
  
    do i=1,nseg+1     ! rotation of xend result stored in xendr
        a=xend(1,i)
        b=xend(2,i)
        c=xend(3,i)
        xendr(1,i)=a*(-cbe*sal*sga+cal*cga) -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
        xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
        xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe 
    enddo
  
    is_positive_z=.true.

    i=2 ! i=1 skipped because (0,0,0)

    do while((i.le.(nseg+1)).and.(is_positive_z)) 
        if(xendr(1,i)<=0.0_dp) is_positive_z=.false.
        i=i+1
    enddo

end function rotation


! rotation only around x-axis
! Euler rotatation X1Z2X3 

function rotationXaxis(xend,xendr,nseg)result(is_positive_z)
  
    use random
    use mathconst
    
    implicit none

    integer,  intent(in)    :: nseg  
    real(dp), intent(in)    :: xend(:,:)
    real(dp), intent(inout) :: xendr(:,:)  
    logical                 :: is_positive_z  ! output 
  
    ! ..local argument 

    integer :: i
    real(dp)  :: fac,fac1,fac2,sbe,cbe,sal,cal,sga
    real(dp)  :: alfa,gama,cga,a,b,c
    
    fac = rands(seed)
    fac1= rands(seed)
    fac2= rands(seed)

    alfa=fac*2.0_dp*pi
    cbe=1.0_dp ! fac1*2.0_dp-1.0_dp
    gama=fac2*2.0_dp*pi
  
    sbe=(1.0_dp-cbe**2)**0.5_dp
    cal=cos(alfa)
    sal=sin(alfa)
    cga=cos(gama)
    sga=sin(gama)
  
    do i=1,nseg+1     ! rotation of xend result stored in xendr
        a=xend(1,i)
        b=xend(2,i)
        c=xend(3,i)

        xendr(1,i)=a ! a*cbe + b*(-cga*sbe)+c*sbe*sga
        xendr(2,i)=a*cal*sbe +b*(cal*cbe*cga-sal*sga) +c*(-cga*sal-cal*cbe*sga)
        xendr(3,i)=a*sbe*sga +b*(cal*sga+cbe*cga*sal) +c*( cal*cga-cbe*sal*sga)
    enddo
  
    is_positive_z=.true.
    i=2 
    do while((i.le.(nseg+1)).and.(is_positive_z)) 
        if(xendr(3,i)<=0.0_dp) is_positive_z=.false.
        i=i+1
    enddo

end function rotationXaxis


! returns Rmat=Rz(alpha)Rx(beta)Rz(gamma)
! exterior rotation over axis z-x-z over (gamma, beta, alpha) 
! equivalently to interior Euler rotatation over z-x'-z'' by (alpha, beta, gamma)

function rotation_matrix(alpha,beta,gamma)result(Rmat)

    real(dp), intent(in)    :: alpha,beta,gamma
    real(dp) :: Rmat(3,3)

    ! local variables
    real(dp) :: cal, sal, cbe, sbe, cga, sga

    cal=cos(alpha)
    sal=sin(alpha)
    cbe=cos(beta)    
    sbe=sin(beta)
    cga=cos(gamma)
    sga=sin(gamma)

    ! Rmat=Rz(alpha)Rx(beta)Rz(gamma)


    Rmat(1,1)= cal*cga-sal*sga*cbe
    Rmat(2,1)= sal*cga+cal*sga*cbe
    Rmat(3,1)= sbe*sga

    Rmat(1,2)= -cal*sga-sal*cbe*cga
    Rmat(2,2)= cal*cbe*cga -sal*sga
    Rmat(3,2)= sbe*cga

    Rmat(1,3)= sal*sbe
    Rmat(2,3)=-cal*sbe
    Rmat(3,3)= cbe


end function rotation_matrix

! Pre: chain(s1,:)={0,0,0)
! Post: Rotate nuclesome conformation chain onto chain_rot 
! such that the three segments (sgraftpts=(s1,s2,s3)
! are located at postion chain_rot(:,s1)=(0,0,0) and 
! chains_rot(:,s2) =(0,0,zval) ( on z-axis) and chain_rot(s3)=(0,yval,uval) ( in z-y  plane) 


subroutine rotate_nucl_chain(chain,chain_rot,sgraftpts,nseg)

    use vectornorm

    real(dp), intent(in)    :: chain(:,:)
    real(dp), intent(inout) :: chain_rot(:,:)  
    integer , intent(in)    :: sgraftpts(3) 
    integer , intent(in)    :: nseg   

    ! .. local variable
    real(dp) :: Rmat(3,3) ,vec2(3), vec3(3),vec3prime(3)
    real(dp) :: alpha, beta, gamma
    integer :: s

    ! .. determine rotation angles 
    
    vec2(:) = chain(:,sgraftpts(2))
    vec3(:) = chain(:,sgraftpts(3))

    ! atan2(Y, X) angle of X+iY between  -pi< angle < pi

    gamma = atan2(vec2(1),vec2(2)) 
    beta  = atan2(sqrt(vec2(1)**2+vec2(2)**2),vec2(3))

    Rmat = rotation_matrix(0.0_dp,beta,gamma) 

    vec3prime = Matmul(Rmat,vec3) ! new coordinate third vector

    alpha=atan2(vec3prime(1),vec3prime(2)) ! rotation such that new third vector is in z-y plane
    
    ! .. make rotation matrix

    Rmat = rotation_matrix(alpha,beta,gamma)

    ! .. apply Rmat onto chain

    do s=1,nseg
        chain_rot(:,s)=Matmul(Rmat,chain(:,s))
    enddo  


end subroutine rotate_nucl_chain

subroutine unit_test_rotation_matrix(info)

    integer, intent(out) :: info

    real(dp) :: Rmat(3,3),Rmattest(3,3),Rdiffsqr(3,3)
    real(dp) :: alpha, beta, gamma
    integer :: i,j
    real(dp) :: chain(3,4),chain_rot(3,4)
    real(dp) :: x(3), Rdiffmax
    real(dp), parameter :: epsRdiffmax= 0.00001_dp


    !  Euler matrix from Mathematica for alpha,beta, gamma
    alpha=1.20_dp
    beta=0.714_dp
    gamma=0.5_dp

    Rmat=transpose(reshape((/-0.0197023_dp, -0.791881_dp, 0.610357_dp,&
        0.949233_dp, -0.206516_dp, -0.237294_dp,&
        0.313958_dp, 0.574695_dp, 0.755748_dp /),shape(Rmat)))

    ! Rotatation
    Rmattest=rotation_matrix(alpha,beta,gamma)
    
    ! difference
    Rdiffsqr=(Rmattest-Rmat)**2
    Rdiffmax=sqrt(sum(Rdiffsqr))


    if(Rdiffmax<epsRdiffmax) then
        info=0
    else
        print*,"Test rotation matrix failed: epsRdiffmax> Rdiffmax: Rdiffmax=",Rdiffmax
        info=1
    endif       
   
end subroutine unit_test_rotation_matrix
 
subroutine unit_test_rotate_nucl_chain(info)

    use vectornorm

    integer, intent(out) :: info

    ! .. local variables

    real(dp) :: chain(3,5),chain_rot(3,5)
    integer :: nseg,sgraftpts(3),i
    real(dp) :: radius
    real(dp), parameter :: epsradius = 0.000001_dp

    ! example of a 'chain conformation'
    chain=reshape((/0.0_dp,0.0_dp,0.0_dp,&
        -1.0_dp, 0.1230_dp,-1.0_dp,&
        -0.0197023_dp, -0.791881_dp, 0.610357_dp,&
        0.949233_dp, -0.206516_dp, -0.237294_dp,&
        0.313958_dp, 0.574695_dp, 0.755748_dp /),shape(chain))

    nseg=5
    sgraftpts=(/1,2,3/)

    call rotate_nucl_chain(chain,chain_rot,sgraftpts,nseg)

        
    ! test
    info=0

    radius=L2norm(chain_rot(:,1),3)
    if(radius> epsradius) then
        print*,"First chain conformation not in origin" 
        info=1
    endif  
    
    radius=sqrt(chain_rot(1,2)**2+chain_rot(2,2)**2)
    if( radius>epsradius) then 
        print*,"Second chain conformation not on z-axis" 
        info=1
    endif

    radius=sqrt(chain_rot(1,3)**2)
    if(radius> epsradius) then
        print*,"Third chain conformation not in z-y plane" 
        info=1
    endif
           
end subroutine unit_test_rotate_nucl_chain

! identical to rotate_nucl_chain plus test of rotated conformation

subroutine rotate_nucl_chain_test(chain,chain_rot,sgraftpts,nseg,write_rotations)

    use vectornorm

    real(dp), intent(in)    :: chain(:,:)
    real(dp), intent(inout) :: chain_rot(:,:)  
    integer , intent(in)    :: sgraftpts(3) 
    integer , intent(in)    :: nseg  
    logical , intent(in)    :: write_rotations
    
    ! .. local variable
    real(dp) :: Rmat(3,3) ,vec2(3), vec3(3),vec3prime(3)
    real(dp) :: alpha, beta, gamma
    integer :: s,i

    integer   :: info 
    real(dp) :: radius
    real(dp), parameter :: epsradius = 0.000001_dp

    ! .. determine rotation angles 
    
    vec2(:) = chain(:,sgraftpts(2))
    vec3(:) = chain(:,sgraftpts(3))

    ! atan2(Y, X) angle of X+iY between  -pi< angle < pi

    gamma = atan2(vec2(1),vec2(2)) 
    beta  = atan2(sqrt(vec2(1)**2+vec2(2)**2),vec2(3))

    Rmat = rotation_matrix(0.0_dp,beta,gamma) 

    vec3prime = Matmul(Rmat,vec3) ! new coordinate third vector

    alpha=atan2(vec3prime(1),vec3prime(2)) ! rotation such that new third vector is in z-y plane
    
    ! .. make rotation matrix

    Rmat = rotation_matrix(alpha,beta,gamma)

    ! .. apply Rmat onto chain

    do s=1,nseg
        chain_rot(:,s)=Matmul(Rmat,chain(:,s))
    enddo   


    ! test
    info=0

    radius=L2norm(chain_rot(:,sgraftpts(1)),3)
    if(radius> epsradius) then
        print*,"First chain conformation not in origin" 
        info=1
    endif  
    
    radius=sqrt(chain_rot(1,sgraftpts(1))**2+chain_rot(2,sgraftpts(2))**2)
    if( radius>epsradius) then 
        print*,"Second chain conformation not on z-axis" 
        info=2
    endif

    radius=sqrt(chain_rot(1,sgraftpts(3))**2)
    if(radius> epsradius) then
        print*,"Third chain conformation not in z-y plane" 
        info=3
    endif
    
    if(info/=0.or.write_rotations) then
        print*,"alpha=",alpha
        print*,"beta=",beta
        print*,"gamma=",gamma
        do s=1,3
            do i=1,3 
                print*,s,i,chain(i,sgraftpts(s)),chain_rot(i,sgraftpts(s))
            enddo  
            print*,""
        enddo    
    endif
             
end subroutine rotate_nucl_chain_test

end module chain_rotation
    
