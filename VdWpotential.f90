! Module to compute VdW/LJ interaction energy for nucleosome conformation

module VdW_potential     

    use precision_definition
  
    implicit none

    integer :: conf_write_com=1 ! counter for write_chain_com_unitvec_lammps_trj

    private :: dotproduct

contains 

    ! Inner product of vectors a and b with dimenstion 3
    ! similar to intrisic function dot_product

    function dotproduct(a,b) result(dotprod)
       
        real(dp), dimension(3) :: a, b 
        real(dp) :: dotprod
      
        dotprod = sum(a*b)

    end function dotproduct


    ! Inits constant varialbke used VGBenergyeffective: makes array segunitvector 
    ! Inits constant parameter of module GB_potential
    ! input : integer segcm
    !         integer nmer
    ! output: integer segunitvector(nmer) 
    !         and all variable of Gay-Bern Potential in module GB_potential

    subroutine init_VGBenergyeffective(segcm,nmer,segunitvector)
        
        use globals, only : nseg, nsegtypes
        use GB_potential, only : init_GB_const

        integer, intent(in) :: segcm(:)
        integer, intent(in) :: nmer
        integer, intent(inout) :: segunitvector(:) 

        integer :: i, deltaseg

        deltaseg=68 ! number from pca python program 18/07/24

        do i=1,nmer
            segunitvector(i)=segcm(i)+deltaseg
        enddo

        call init_GB_const()
 
    end subroutine init_VGBenergyeffective

    
    ! pre : chain conformation and segcm and segunitvector
    ! post: effective VGB energy conformation using location COM of nucleosome
    ! VLJ(u1,u2,r) = Gay-Berne potential 
    ! NB This function is equivalent to VGBenergyeffective, but latter one nicier !!!

    function  VGBenergyeffective_comb(chain,nmer,segcm,segunitvector)result(Energy)

        use globals, only : nseg, nsegtypes
        use chains, only : type_of_monomer
        use GB_potential, only : GBpotential

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        integer, intent(in) :: segcm(:)
        integer, intent(in) :: segunitvector(:) 

        real(dp) :: Energy
        
        real(dp) :: Ene
        integer ::  i,j
        integer :: isegcm, isegunitvector, jsegcm, jsegunitvector
        real(dp) :: xi, xj, yi, yj, zi, zj
        real(dp) :: xiunit, xjunit, yiunit, yjunit, ziunit,zjunit
        real(dp) :: u1(3), u2(3), norm_u1, norm_u2 ! unitvector and its norm
        real(dp) :: rvec(3)  ! distance between com

        Ene=0.0_dp 
       
        do i=1,nmer

            isegcm=segcm(i) ! atom_id of AA that is (closest to) COM
            xi = chain(1,isegcm)
            yi = chain(2,isegcm)
            zi = chain(3,isegcm)

            isegunitvector=segunitvector(i) ! atom_id of AA with isegcm to make unit vector
            xiunit = chain(1,isegunitvector)
            yiunit = chain(2,isegunitvector)
            ziunit = chain(3,isegunitvector)

            u1=(/xiunit-xi,yiunit-yi,ziunit-zi/) ! unit vector

            norm_u1=sqrt(dotproduct(u1,u1)) 
            u1=u1/norm_u1

            do j=i+1,nmer

                jsegcm=segcm(j)
                xj = chain(1,jsegcm)
                yj = chain(2,jsegcm)
                zj = chain(3,jsegcm) 

                jsegunitvector=segunitvector(j)

                xjunit = chain(1,jsegunitvector)
                yjunit = chain(2,jsegunitvector)
                zjunit = chain(3,jsegunitvector)

                u2=(/xjunit-xj,yjunit-yj,zjunit-zj/)
                norm_u2=sqrt(dotproduct(u2,u2)) 
                u2=u2/norm_u2

                rvec=(/xj-xi,yj-yi,zj-zi/)
               
                Ene=Ene + GBpotential(u1,u2,rvec)
            
            enddo
        enddo     

        Energy = Ene          

    end function VGBenergyeffective_comb

    ! Calculate sum value Gay-Berne potential for nucleosome conformation 
    ! input :  real(dp :: chainchain(:,:) conformation 
    !          integer :: nmer : nnucl
    !          integer :: segcm(:) : atom_id closest to com 
    !          integer :: segunitvector(:) atom_id used to make unit vector
    ! output/return:  real(dp) :: energy 

    function  VGBenergyeffective(chain,nmer,segcm,segunitvector)result(Energy)

        use GB_potential, only : GBpotential

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        integer, intent(in) :: segcm(:)
        integer, intent(in) :: segunitvector(:) 

        real(dp) :: Energy

        real(dp) :: rcom(3,nmer)
        real(dp) :: uvector(3,nmer)


        call make_com_unit_vector_nucl(chain,nmer,segcm,segunitvector,rcom,uvector)

        Energy = VGBenergyeff(nmer,rcom,uvector)

    end function VGBenergyeffective


    !  Calculate sum of value Gay-Berne potential for given set of com and unitvectors
    !  input : real(dp) : rcom(:,:) : ceom of nucleosome 
    !          real(dp) : uvector(:,:) : orientation vector ( normalized) of each nucleosome 
    !  output/return : real(dp)  energy : effective VGB energy = VLJ(u1,u2,r)  

    function  VGBenergyeff(nmer,rcom,uvector)result(Energy)

        use GB_potential, only : GBpotential

        integer, intent(in)  :: nmer
        real(dp), intent(in) :: rcom(:,:)
        real(dp), intent(in)  :: uvector(:,:) 

        real(dp) :: Energy
        
        real(dp) :: Ene
        integer ::  i,j
        real(dp) :: uveci(3), uvecj(3) ! unitvector
        real(dp) :: rvec(3), rveci(3), rvecj(3) ! distance between com's

        Ene=0.0_dp 
       
        do i=1,nmer

            rveci = rcom(:,i)
            uveci = uvector(:,i)

            do j=i+1,nmer

                rvecj =rcom(:,j)
                uvecj =uvector(:,j) 
                rvec=rvecj-rveci

                Ene = Ene + GBpotential(uveci,uvecj,rvec)
            
            enddo
        enddo     

        Energy = Ene            

    end function VGBenergyeff


    ! Determine com and unit vector for a given nucloesome conformation
    ! input :  real(dp :: chainchain(:,:) conformation 
    !          integer :: nmer : nnucl
    !          integer :: segcm(:) : atom_id closest to com 
    !          integer :: segunitvector(:) atom_id used to make unit vector
    ! output:  real(dp) :: rcom(:,:)    : com of each nucleosome 
    !          real(dp) :: uvector(:,:) : unit vector of each nucleosome

    subroutine make_com_unit_vector_nucl(chain,nmer,segcm,segunitvector,rcom,uvector)

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        integer, intent(in) :: segcm(:)
        integer, intent(in) :: segunitvector(:) 
        real(dp), intent(inout) :: rcom(:,:)
        real(dp), intent(inout) :: uvector(:,:)
      

        integer :: i, isegcm, isegunitvec
        real(dp) :: rvec(3), rvecunitvec(3), uvec(3), norm_uvec

        do i=1,nmer

            isegcm=segcm(i) ! atom_id of AA that is (closest to) COM
            rvec = chain(:,isegcm)
        
            isegunitvec=segunitvector(i) ! atom_id of AA with isegcm to make unit vector
            rvecunitvec = chain(:,isegunitvec)

            uvec=rvecunitvec-rvec   ! unit vector
            norm_uvec=sqrt(dot_product(uvec,uvec))  
            uvec=uvec/norm_uvec  ! normalize

            rcom(:,i) = rvec
            uvector(:,i) = uvec

        enddo     

    end subroutine make_com_unit_vector_nucl

    !Computes value of parameteric function s and t of plane that is spanned by 
    !v=vec[2]-vec[1] and w=vec[3]-vec[1] at origin: x= s*v + t*w +vec[1]=0  
    !input : real(dp) :: vec(3,3) three xyz-coordinate triplets/vectors 
    !output: integer ::  s, t:  value of parameteric function plane as origin
    !        real(dp) :: vec_origin(3)

    subroutine parameter_value_origin(vec,s,t,vec_origin)

        real(dp), intent(in)   :: vec(3,3) 
        integer, intent(inout) ::  s, t
        real(dp), intent(inout) :: vec_origin(3) 

        real(dp) :: v1(3),v2(3), u1(3), u2(3),  sqr_norm_u1,  sqr_norm_u2

        v1=vec(:,2)-vec(:,1)
        v2=vec(:,3)-vec(:,1)
        
        ! orthogonal vectors
        u1=v1
        u2=v2-(dot_product(v2,u1)/dot_product(u1,u1))*u1
    
        sqr_norm_u1=dot_product(u1,u1)
        sqr_norm_u2=sqrt(dot_product(u2,u2))
    
        ! value of 
        s = -dot_product(u1,vec(:,1))/sqr_norm_u1
        t = -dot_product(u2,vec(:,1))/sqr_norm_u2
    
        ! compute (relative) origin/check
   
        vec_origin= s*u1 + t*u2 +v1

    end subroutine parameter_value_origin     
    
   
    subroutine make_com_nucl_rotation(chain,nmer,unitvector_triplets,rcom,uvector)

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        integer, intent(in) :: unitvector_triplets(:,:)  ! dimension : nmer x 3  
        real(dp), intent(inout) :: rcom(:,:)
        real(dp), intent(inout) :: uvector(:,:)
      

        integer :: i, k, s0, t0
        real(dp) :: vec(3,3), vec_origin(3)
         
        ! determine com using unitvectortriplets
        ! determine s0 and t0 parameter of reference vectors

        do i=1,nmer
            ! get vector spanning origin-plane
            do k=1,3
                vec(:,k)=chain(:,unitvector_triplets(i,k))
            enddo
            call parameter_value_origin(vec,s0,t0,vec_origin)
            rcom(:,i) = vec_origin
        enddo     

    end subroutine make_com_nucl_rotation


    ! Writes the coordinates of com of nucleosomes orientation vector associated with com to a lammps trajectory file.
    ! Unit length in Angstrom
    ! pre : file opened with open_com_unitvec_lammps_trj
    ! input: integer un_trj : unit number of file to write to.
    !        real(dp) segcom(:) : coordinates of com of nucleosomes
    !        real(dp) segunitvector : orientation vector associated with com 
    !    
    ! return: a lammpstrj file

    subroutine write_chain_com_unitvec_lammps_trj(un_trj,segcom,segunitvector)

        use globals, only : nnucl
        use volume, only : delta, nx, ny, nz 
        use quaternions, only :  rot_axis_angle, rot_axis_angle_to_quat
        use GB_potential, only : sigmaS,sigmaE  

        real(dp), dimension(:,:), intent(in) :: segcom
        real(dp), dimension(:,:), intent(in) :: segunitvector
        integer, intent(in) :: un_trj
        
        real(dp) :: x, y, z
        integer ::  n
        real(dp) :: xbox0, xbox1(3)
        integer :: idatom, item, conf
        real(dp) :: shapex,shapey,shapez
        real(dp) :: rcom(3), uvec(3), uref(3), qu(4), rotaxis(3), theta

        xbox0=0.0_dp 
        xbox1(1)=nx*delta*10.0_dp ! conversion form nm -> Angstrom
        xbox1(2)=ny*delta*10.0_dp
        xbox1(3)=nz*delta*10.0_dp

        item=0
        conf_write_com=conf_write_com+1
        
        ! axis length of Gay-Berne potential
        shapex=sigmaS
        shapey=sigmaS
        shapez=sigmaE
        ! reference unit vector
        uref=(/0.0_dp,0.0_dp,1.0_dp/) ! not sure ...


        ! write preamble 
        write(un_trj,'(14A)')'ITEM: TIMESTEP' 
        write(un_trj,*)conf_write_com 
        write(un_trj,'(21A)')'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nnucl
        write(un_trj,'(25A)')'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1(1)
        write(un_trj,*)xbox0,xbox1(2)
        write(un_trj,*)xbox0,xbox1(3)
        write(un_trj,'(70A)')'ITEM: ATOMS id type x y z quatw quati quatj quatk shapex shapey shapez'
        
        do n=1,nnucl
            
            rcom = segcom(:,n)
            uvec = segunitvector(:,n)

            call rot_axis_angle(uref, uvec, rotaxis, theta) 
            qu = rot_axis_angle_to_quat(rotaxis, theta) ! ???

            item=item+1  
            write(un_trj,*)item,idatom,rcom(1),rcom(2),rcom(3),qu(1),qu(2),qu(3),qu(4),shapex,shapey,shapez

        enddo    
        
      
    end subroutine write_chain_com_unitvec_lammps_trj


    ! Creates and opens file for lammps trajectory film
    ! Writing done by write_chain_com_unitvec_lammps_trj
    ! input:        integer : info : value 0 no errors
    ! output/return integer : un_trj : unit nubmer of file

    function open_chain_com_unitvec_lammps_trj(info)result(un_trj)

        use mpivars, only : rank
        use myutils, only : newunit, lenText
        use myio, only : myio_err_chainsfile

        integer, intent(inout) :: info

        integer :: un_trj

        ! local
        character(len=lenText) :: istr
        character(len=35) :: fname
        integer :: ios
        logical :: exist
        
        info = 0    

        write(istr,'(I4)')rank
        fname='trajcomunit'//trim(adjustl(istr))//'.lammpstrj'
        inquire(file=fname,exist=exist)
        if(.not.exist) then
            open(unit=newunit(un_trj),file=fname,status='new',action="write",iostat=ios)
            if(ios > 0 ) then
                print*, 'Error opening : ',fname,' file : iostat =', ios
                info = myio_err_chainsfile
                return
            endif
        else
            open(unit=newunit(un_trj),file=fname,position="append",status='old',action="write", iostat=ios)
            if(ios > 0 ) then
                print*, 'Error opening : ',fname,' file : iostat =', ios
                info = myio_err_chainsfile
                return
            endif
        endif   

    end function open_chain_com_unitvec_lammps_trj


end module VdW_potential    