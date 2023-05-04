! --------------------------------------------------------------|
!                                                               | 
! chainsgenerator.f90:                                          |       
! generator chains in center of box                             |
! --------------------------------------------------------------|


module chaingenerator 

    use precision_definition   
    implicit none

   ! integer, parameter :: lenfname=40
    integer :: conf_write
    real(dp) :: xgraftloop(3,2)

    real(dp), parameter :: eps_equilat=1.0e-8_dp

    private 

    public :: make_chains, make_chains_mc,read_chains_XYZ
    public :: make_charge_table, make_segcom, make_sequence_chain, make_type_of_charge_table
    public :: set_mapping_num_to_char, set_properties_chain, write_chain_struct
    
    private ::  eps_equilat

contains


subroutine make_chains(chainmethod,systype)

    use myutils, only : print_to_log, LogUnit, lenText, error_handler
    use myio, only : myio_err_chainmethod

    character(len=15), intent(in) :: chainmethod
    character(len=15), intent(in) :: systype 

    integer :: i, info
    character(len=lenText) :: text, istr

    info=0

    select case (chainmethod)
    case ("MC")
        call make_chains_mc()
    case ("FILE_XYZ")      
        call read_chains_xyz(systype,info)       
    case default
        text="chainmethod not equal to MC or FILE_XYZ"
        call print_to_log(LogUnit,text)
        print*,text
        info=myio_err_chainmethod
    end select

    if(info/=0) call error_handler(info,"make_chains")

end subroutine make_chains


!  init of cuantas polymer configurations of polymer chain anchored onto a flat surface 

subroutine make_chains_mc()
  
    use mpivars
    use globals
    use chains
    use random
    use parameters, only : geometry, lseg, write_mc_chains, isVdW, isVdWintEne
    use parameters, only : maxnchainsrotations, maxnchainsrotationsxy
    use volume, only : nx, ny, nz, delta
    use volume, only : coordinateFromLinearIndex, linearIndexFromCoordinate
    use volume, only : ut, vt
    use myutils
    use cadenas_linear
    use cadenas_sequence

    !     .. variable and constant declaractions      

    integer :: i,j,k,s,g,gn      ! dummy indices
    integer :: idx               ! index label
    integer :: ix,iy,idxtmp,ntheta
    integer :: nchains           ! number of rotations
    integer :: maxnchains        ! number of rotations
    integer :: maxntheta         ! maximum number of rotation in xy-plane
    integer :: conf              ! counts number of conformations
    integer :: allowedconf
    real(dp) :: chain(3,nseg,200) ! chain(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp) :: chain_rot(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg) ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates
    real(dp) :: xpp(nseg), ypp(nseg), zpp(nseg)  
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    real(dp) :: theta, theta_angle
    character(len=lenText) :: text, istr
    integer  :: xi,yi,zi ,un_trj, un_ene, segcenter
    real(dp) :: energy   
    logical :: saw     
   
    !  .. executable statements
    !  .. initializations of variables     
       
    conf = 1                 ! counter for conformations
    seed = 435672*(rank+1)   ! seed for random number generator  different on each node
    maxnchains = maxnchainsrotations
    maxntheta = maxnchainsrotationsxy         ! maximum number of rotation in xy-plane  
    !theta_angle = 2.0_dp*pi/maxntheta
    
    Lz = nz*delta            ! maximum height box 
    Lx = nx*delta            ! maximum width box 
    Ly = ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp           ! center box
    ycm= Ly/2.0_dp
    zcm= LZ/2.0_dp
   
    energy=0.0_dp
            
    if(write_mc_chains) then 
        conf_write=0
        un_trj= open_chain_lammps_trj()
        un_ene= open_chain_energy()
    endif   

    if(isHomopolymer.eqv..FALSE.) then 
        allocate(lsegseq(nseg))
        call make_lsegseq(lsegseq,nseg)
    endif    
  
    do while (conf<=max_confor)

        nchains= 0      ! init zero 
        if(isHomopolymer) then 
            call make_linear_chains(chain,nchains,maxnchains,nseg,lseg) ! chain generator f90
        else
            call make_linear_seq_chains(chain,nchains,maxnchains,nseg) 
        endif  

        if(write_mc_chains) then 
            energy=0.0_dp
            write(un_ene,*)energy  
            call write_chain_lammps_trj(un_trj,chain,nchains)  
        endif    
   
            
        if(geometry=="cubic") then

            do j=1,nchains   
            
                do s=1,nseg                          !  transforming form real- to lattice coordinates
                    z(s) = chain(1,s,j)+zcm
                    x(s) = chain(2,s,j)+xcm
                    y(s) = chain(3,s,j)+ycm
                enddo
                
                do s=1,nseg 

                    ! .. periodic boundary conditions in x-direction and y-direction 
                    chain_rot(2,s) = pbc(x(s),Lx)
                    chain_rot(3,s) = pbc(y(s),Ly)
                    chain_rot(1,s) = pbc(z(s),Lz)        

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_rot(2,s)/delta)+1
                    yi = int(chain_rot(3,s)/delta)+1
                    zi = int(chain_rot(1,s)/delta)+1
                    
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                    indexchain(s,conf) = idx
                    if(idx<=0) then
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                    endif
                enddo            
        
                conf = conf +1 

            
            enddo             ! end loop over nchains

        else if(geometry=="prism") then
            
            do j=1,nchains   
            
                do s=1,nseg                     !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)+zcm
                    xp(s) = chain(2,s,j)+xcm
                    yp(s) = chain(3,s,j)+ycm
                enddo

                do s=1,nseg
                        
                    ! .. prism coordinates

                    x(s) = ut(xp(s),yp(s))
                    y(s) = vt(xp(s),yp(s))
                    
                    ! .. periodic boundary conditions in x-direction and y-direction an z-direction 
                    chain_rot(2,s) = pbc(x(s),Lx)
                    chain_rot(3,s) = pbc(y(s),Ly)
                    chain_rot(1,s) = pbc(zp(s),Lz)       

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_rot(2,s)/delta)+1
                    yi = int(chain_rot(3,s)/delta)+1
                    zi = int(chain_rot(1,s)/delta)+1

                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                    indexchain(s,conf) = idx
                    if(idx<=0) then
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                    endif
                    
                enddo         ! end loop over graft points

                conf = conf +1 
             
            enddo                 ! end loop over nchains

        else 
           print*,"Error: in make_chains_mc: geometry not cubic or prism"
           print*,"stopping program"
           stop
        endif
     
    enddo                     ! end while loop
    
    !  .. end chains generation 
      
    write(istr,'(I4)')rank
    text='AB Chains generated on node '//istr
    call print_to_log(LogUnit,text)

  
    if(isHomopolymer.eqv..FALSE.) deallocate(lsegseq)

    if(write_mc_chains) then
        close(un_trj)
        close(un_ene)
    endif   

    energychain=0.0_dp ! no internal energy 

end subroutine make_chains_mc


subroutine read_chains_xyz(systype,info)

    
    ! .. argument
    character(len=15), intent(in) :: systype
    integer, intent(out) :: info

    if(systype=="nucl_neutral_sv".or.systype=="nucl_ionbin_sv") then 
        call read_chains_xyz_nucl_volume(info)
    else
        call read_chains_xyz_nucl(info)
    endif    

end subroutine


! Reads conformations from a file called traj.<rank>.xyz
! Reads energy of the conformation from a file called energy.<rank>.ene
! Reading of energy file only if isChainEnergyFile==.true.
! Format conformation file : lammps trajectory file xyz file
! number of ATOMS much equal nseg  
! 
! returns : indexchain and weightchain
!         : structural quantities:  
!           Rgsqr   = radius of gyration,
!           Rendsqr = end_to_end distance,
!           bond_angle,
!           dihedral_angle,
!           nucl_spacing = nucleosomal_spacing

subroutine read_chains_xyz_nucl(info)

    !     .. variable and constant declaractions  
    use mpivars, only : rank,size                                                                                   
    use globals
    use chains
    use random
    use parameters
    use volume, only :  nx, ny,nz, delta
    use chain_rotation, only : rotate_nucl_chain, rotate_nucl_chain_test
    use myio, only : myio_err_chainsfile, myio_err_energyfile, myio_err_index
    use myio, only : myio_err_conf, myio_err_nseg, myio_err_geometry, myio_err_equilat
    use myutils,  only :  print_to_log, LogUnit, lenText, newunit
    use Eqtriangle

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    integer :: i,j,s,rot,g,gn      ! dummy indices
    integer :: idx                 ! index label
    integer :: ix,iy,iz,idxtmp,ntheta
    integer :: nchains              ! number of rotations
    integer :: maxnchains           ! number of rotations
    integer :: maxntheta            ! maximum number of rotation in xy-plane
    integer :: conf,conffile        ! counts number of conformations  
    integer :: nsegfile             ! nseg in chain file      
    integer :: cuantasfile          ! cuantas in chain file                                              
    real(dp) :: chain(3,nseg),chain_rot(3,nseg),chain_pbc(3,nseg)  ! chains(x,i)= coordinate x of segement i
    real(dp) :: xseg(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg)    ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: xpp(nseg),ypp(nseg)
    integer  :: xi,yi,zi
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm ! sizes box and center of mass box
    real(dp) :: xpt,ypt              ! coordinates
    real(sp) :: xc,yc,zc               
    real(dp) :: energy                                             
    character(len=25) :: fname
    integer :: ios, rankfile, iosene
    character(len=30) :: str
    real(dp) :: scalefactor
    integer :: un,unw,un_ene ! unit number
    logical :: exist
    character(len=lenText) :: text,istr
    real(dp) :: d_type_num, d_atom_num
    integer :: i_type_num, i_atom_num
    logical :: isReadGood
    character(len=80), parameter  :: fmt3reals = "(5F25.16)"
    real(dp) :: equilat, equilat_rot

    ! .. executable statements   

    info=0

    ! .. open file   

    rankfile=rank                                                                                     
    
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.xyz'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,' traj file :',fname,' does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    if(isChainEnergyFile) then
        write(istr,'(I4)')rankfile
        fname='energy.'//trim(adjustl(istr))//'.ene'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un_ene),file=fname,status='old',iostat=ios)
        else
            text='energy.'//trim(adjustl(istr))//'.ene file does not exit'
            print*,text
            info = myio_err_energyfile
            return
        endif
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            info = myio_err_energyfile
            return
        endif
    endif    

    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file    
    ios=0
    scalefactor=unit_conv
    energy=0.0_dp
    seed=435672               ! seed for random number generator                                                                               
    
    ios=0

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center box
    ycm= Ly/2.0_dp
    zcm= Lz/2.0_dp

    isReadGood=.true. 

    do while ((conf<=max_confor).and.isReadGood)
    

        if(conf.ne.1) then ! skip lines
            read(un,*,iostat=ios)
            read(un,*,iostat=ios) 
        else               ! read preamble
            read(un,*,iostat=ios)nsegfile
            if(ios/=0) isReadGood=.false.
            read(un,*,iostat=ios) ! skip line
            if(ios/=0) isReadGood=.false.   
        
            if(nsegfile.ne.nseg) then 
                text="nseg chain file not equal internal nseg : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_nseg 
                return
            endif    
        endif
        
        do s=1,nseg              ! .. read form  trajectory file
    
            read(un,*,iostat=ios)xc,yc,zc
            if(ios/=0) isReadGood=.false. 
            
            xseg(1,s) = xc*scalefactor 
            xseg(2,s) = yc*scalefactor  
            xseg(3,s) = zc*scalefactor 
            
        enddo

        if(isChainEnergyFile) read(un_ene,*,iostat=ios)energy

        if(isReadGood) then ! read was succesfull  

            conffile=conffile +1 
           
            do s=1,nseg        
                chain(1,s) = xseg(1,s)-xseg(1,sgraftpts(1)) 
                chain(2,s) = xseg(2,s)-xseg(2,sgraftpts(1)) 
                chain(3,s) = xseg(3,s)-xseg(3,sgraftpts(1)) 
            enddo


          
            ! rotate chain 
            
            ! call rotate_nucl_chain(chain,chain_rot,sgraftpts,nseg)
            call rotate_nucl_chain_test(chain,chain_rot,sgraftpts,nseg,write_rotations)                

            ! extra test rotation
            equilat=equilateralness_vectors(chain(:,sgraftpts(1)),chain(:,sgraftpts(2)),chain(:,sgraftpts(3)))
            equilat_rot=equilateralness_vectors(chain_rot(:,sgraftpts(1)),chain_rot(:,sgraftpts(2)),chain_rot(:,sgraftpts(3)))
 
            if( abs(equilat-equilat_rot)>=eps_equilat )  then
                text="Equilateralness test failed : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_equilat 
                return
            endif    
            

            select case (geometry)
            case ("cubic")

                do s=1,nseg

                    chain_pbc(1,s) = pbc(chain_rot(1,s)+xcm,Lx) ! periodic boundary conditions in x and y and z direcxtion 
                    chain_pbc(2,s) = pbc(chain_rot(2,s)+ycm,Ly)
                    chain_pbc(3,s) = pbc(chain_rot(3,s)+zcm,Lz)  

                    ! transforming form real- to lattice coordinates                 
                    xi = int(chain_pbc(1,s)/delta)+1
                    yi = int(chain_pbc(2,s)/delta)+1
                    zi = int(chain_pbc(3,s)/delta)+1
                        
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                    indexchain(s,conf) = idx

                    if(idx<=0.or.idx>nsize) then   

                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text                          
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return
                    endif
                    
                enddo

                energychain(conf)      = energy

                Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                nucl_spacing(:,conf)   = nucleosomal_spacing(chain_pbc,nnucl,segcm)
                
                conf=conf+1   
                                    
            case("prism") 
                    
                do s=1,nseg
                    
                    xp(s) = chain_rot(1,s)+xcm
                    yp(s) = chain_rot(2,s)+ycm
                    zp(s) = chain_rot(3,s)+zcm

                    ! .. transformation to prism coordinates 
                    xpp(s) = ut(xp(s),yp(s))
                    ypp(s) = vt(xp(s)+ycm,yp(s))

                    ! .. periodic boundary conditions in u and v ands z direction
                    chain_pbc(1,s) = pbc(xpp(s),Lx) 
                    chain_pbc(2,s) = pbc(ypp(s),Ly)
                    chain_pbc(3,s) = pbc(zp(s),Lz)        

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_pbc(1,s)/delta)+1
                    yi = int(chain_pbc(2,s)/delta)+1
                    zi = int(chain_pbc(3,s)/delta)+1
                        
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                    indexchain(s,conf) = idx

                    if(idx<=0.or.idx>nsize) then    
                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text                          
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return
                    endif
                    
                enddo
                
                energychain(conf)      = energy

                Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                nucl_spacing(:,conf)   = nucleosomal_spacing(chain_pbc,nnucl,segcm)

                conf=conf+1   

            case default
                text="Error: in make_chains_XYZ_nucl geometry not cubic or prism: stopping program"
                call print_to_log(LogUnit,text)
                info = myio_err_geometry
                return 
                    
            end select

        endif   ! read 

    enddo       ! end while loop                                                                                                          
    ! end chains generation    
    
    conf=conf-1  ! lower by one  

    if(conf<max_confor) then
        print*,"subroutine make_chains_XYZ_nucl :" 
        print*,"conf     = ",conf," less then imposed max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas=conf   
        info=myio_err_conf        
    else
        text="Chains generated: subroutine make_chains_XYZ_nucl"
        call print_to_log(LogUnit,text)
        readinchains=conffile
        info=0
    endif


    write(istr,'(L2)')isVdWintEne
    text="isVdWintEne = "//trim(adjustl(istr))
    call print_to_log(LogUnit,text)

    if(.not.(isChainEnergyFile)) energychain=0.0_dp

    close(un) 
    if(isChainEnergyFile) close(un_ene)


    call normed_weightchains()     

    deallocate(energychain) ! free unused variables 

    
end subroutine read_chains_XYZ_nucl


! Reads conformations from a file called traj.<rank>.xyz
! Reads energy of the conformation from a file called energy.<rank>.ene
! Reading of energy file only if isChainEnergyFile==.true.
! Format conformation file : lammps trajectory file xyz file
! number of ATOMS much equal nseg  
! 
! returns : indexconf and weightchain for fcnnucl_neutral_sv 
!         : structural quantities:  
!           Rgsqr   = radius of gyration,
!           Rendsqr = end_to_end distance,
!           bond_angle,
!           dihedral_angle,
!           nucl_spacing = nucleosomal_spacing

subroutine read_chains_xyz_nucl_volume(info)

    !     .. variable and constant declaractions  
    use mpivars, only  : rank,size                                                                                   
    use globals
    use chains
    !use random
    use parameters
    use volume, only   :  nx, ny,nz, delta
    use chain_rotation, only : rotate_nucl_chain, rotate_nucl_chain_test 
    use chain_rotation, only : orientation_vector_ref, orientation_vector, rotate_chain_elem
    use myio, only     : myio_err_chainsfile, myio_err_energyfile, myio_err_index
    use myio, only     : myio_err_conf, myio_err_nseg, myio_err_geometry, myio_err_equilat
    use myutils, only  : print_to_log, LogUnit, lenText, newunit
    use Eqtriangle

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    integer :: i,j,s,rot,g,gn,k     ! dummy indices
    integer :: idx                  ! index label
    integer :: ix,iy,iz,idxtmp !,ntheta
    !integer :: nchains              ! number of rotations
    !integer :: maxnchains           ! number of rotations
    !integer :: maxntheta            ! maximum number of rotation in xy-plane
    integer :: conf,conffile        ! counts number of conformations  
    integer :: nsegfile             ! nseg in chain file      
    integer :: cuantasfile          ! cuantas in chain file                                              
    real(dp) :: chain(3,nseg),chain_rot(3,nseg),chain_pbc(3,nseg),chain_pbc_tmp(3)    ! chains(x,i)= coordinate x of segement i
    real(dp) :: xseg(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg)    ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: xpp(nseg),ypp(nseg)
    integer  :: xi,yi,zi, ri(3)
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm, Lr(3), rcm(3) ! sizes box and center of mass box
    real(dp) :: xpt,ypt              ! coordinates
    real(sp) :: xc,yc,zc            ! coordinates in single precision            
    real(dp) :: energy                                             
    character(len=25) :: fname
    character(lenText):: fname2
    integer :: ios, rankfile, iosene
    character(len=30) :: str
    real(dp) :: scalefactor
    integer :: un,unw,un_ene ! unit number
    logical :: exist
    character(len=lenText) :: text,istr
    real(dp) :: d_type_num, d_atom_num
    integer :: i_type_num, i_atom_num
    logical :: isReadGood
    integer :: nrotpts
    character(len=80), parameter  :: fmt3reals = "(5F25.16)"
    real(dp) :: equilat,equilat_rot
    integer :: nelem2(3),nsegAA2 
    integer  :: segnumAAstart(nnucl), segnumAAend(nnucl) ! segment numbers first/last AAs 
    integer  :: orient_triplet_ref(3)
    real(dp) :: orient_vector_ref(3)
    real(dp) :: orient_vectors(3,nnucl)

    type(var_darray), dimension(:,:), allocatable :: chain_elem, chain_elem_index
    type(var_darray), dimension(:,:,:), allocatable :: chain_elem_rot

    ! .. executable statements   

    info=0

    ! .. open file   

    rankfile=rank                                                                                     
    
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.xyz'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,' traj file :',fname,' does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    if(isChainEnergyFile) then
        write(istr,'(I4)')rankfile
        fname='energy.'//trim(adjustl(istr))//'.ene'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un_ene),file=fname,status='old',iostat=ios)
        else
            text='energy.'//trim(adjustl(istr))//'.ene file does not exit'
            print*,text
            info = myio_err_energyfile
            return
        endif
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            info = myio_err_energyfile
            return
        endif
    endif    

    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file    
    ios=0
    scalefactor=unit_conv
    energy=0.0_dp
    !seed=435672              ! seed for random number generator                                                                               
    
    ios=0

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center box
    ycm= Ly/2.0_dp
    zcm= Lz/2.0_dp

    rcm=(/xcm,ycm,zcm/)
    Lr=(/Lx,Ly,Lz/)

    nrotpts=sgraftpts(1)  ! nucleosome id /segment number around which to rotate whole conformation
    !nrotpts=rotation_triplets(1)
    
    ! return position (chain_elem) and number (nelem) of elements of every AA segment
                        
    call read_nucl_elements(mtpdbfname,nsegAA,nelemAA,chain_elem,typeAA,vnucl,elem_charge,info)
    call print_nucl_elements(nsegAA,nelemAA,chain_elem)
    call read_nucl_orient_triplets(orientfname,nnucl,orientation_triplets,info)

    call compute_segnumAAstart(nseg,nsegtypes,nnucl,segnumAAstart)
    call compute_segnumAAend(nnucl,nsegAA,segnumAAstart,segnumAAend)

    do i=1,nnucl
        print*," sAAstart ",segnumAAstart(i)," sAAend ",segnumAAend(i)
    enddo    

    ! compare ORDER of sgraftpts and orientation triplets !!!!!!!
    do i=1,3
        orient_triplet_ref(i)=orientation_triplets(1,i)-segnumAAstart(1)+1
        print*,"orient_triplet_ref(",i,")=",orient_triplet_ref(i)," original ",orientation_triplets(1,i)
    enddo

    ! make chain_elem relative to CA of AA of CM rotation i.e., sgraftpts
    call shift_nucl_elements(nsegAA,nelemAA,orient_triplet_ref(1),chain_elem)
    call print_nucl_elements(nsegAA,nelemAA,chain_elem)

    call make_nelem(nseg,nsegAA,nnucl,segnumAAstart,nelemAA,nelem)
  
    call allocate_indexconf(cuantas,nseg,nelem)
    call allocate_nucl_chain_elements(nnucl,nsegAA,nelemAA,chain_elem_rot) 
    call allocate_chain_elements(nseg,nelem,chain_elem_index)
    call orientation_vector_ref(chain_elem,orient_triplet_ref,orient_vector_ref)

    isReadGood=.true. 

    do while ((conf<=max_confor).and.isReadGood)
    

        if(conf.ne.1) then ! skip lines
            read(un,*,iostat=ios)
            read(un,*,iostat=ios) 
        else               ! read preamble
            read(un,*,iostat=ios)nsegfile
            if(ios/=0) isReadGood=.false.
            read(un,*,iostat=ios) ! skip line
            if(ios/=0) isReadGood=.false.   
        
            if(nsegfile.ne.nseg) then 
                text="nseg chain file not equal internal nseg : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_nseg 
                return
            endif    
        endif
        
        do s=1,nseg              ! .. read form  trajectory file
    
            read(un,*,iostat=ios)xc,yc,zc
            if(ios/=0) isReadGood=.false. 
            
            xseg(1,s) = xc*scalefactor 
            xseg(2,s) = yc*scalefactor  
            xseg(3,s) = zc*scalefactor 
            
        enddo

        if(isChainEnergyFile) read(un_ene,*,iostat=ios)energy

        if(isReadGood) then ! read was succesfull  

            conffile=conffile +1 
           
            ! -1. translation such that chain(:,nrotpts) in origin

            do s=1,nseg        
                chain(1,s) = xseg(1,s)-xseg(1,nrotpts) 
                chain(2,s) = xseg(2,s)-xseg(2,nrotpts) 
                chain(3,s) = xseg(3,s)-xseg(3,nrotpts) 
            enddo

            ! 0. rotate chain 
            
            call rotate_nucl_chain_test(chain,chain_rot,sgraftpts,nseg,write_rotations) 
           
            ! extra test rotation
            equilat=equilateralness_vectors(chain(:,sgraftpts(1)),chain(:,sgraftpts(2)),chain(:,sgraftpts(3)))
            equilat_rot=equilateralness_vectors(chain_rot(:,sgraftpts(1)),chain_rot(:,sgraftpts(2)),chain_rot(:,sgraftpts(3)))
 
            if( abs(equilat-equilat_rot)>=eps_equilat )  then
                text="Equilateralness test failed : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_equilat 
                return
            endif  

            ! 1. determine orientation chain_rot : n 
            ! 2. get orientation vector of all nnucl nucleosomes in chain_rot and of reference in chain_elem
                
            call orientation_vector(chain_rot,orientation_triplets,orient_vectors)

            ! 3. make rotation axis/angle for all nnuclnucleosome
            ! 4. apply rotation for all nucleosome to chain_elem
                
            call rotate_chain_elem(orient_vector_ref,orient_vectors,segnumAAstart,nelemAA,chain_elem,chain_elem_rot)
        
            ! 5. add chain+chain_elem_rot together
            
            call add_chain_rot_and_chain_elem_rot(nseg,nsegAA,nnucl,segnumAAstart,segnumAAend,nelemAA,&
                chain_rot,chain_elem_rot,chain_elem_index) 

            ! 6. make indexconfig i.e. place onto lattice

            select case (geometry)
            case ("cubic")

                do s=1,nseg

                    
                    chain_pbc(1,s) = pbc(chain_rot(1,s)+xcm,Lx) ! periodic boundary conditions in x and y and z direcxtion  
                    chain_pbc(2,s) = pbc(chain_rot(2,s)+ycm,Ly) 
                    chain_pbc(3,s) = pbc(chain_rot(3,s)+zcm,Lz) 
                    
                    xi  = int(chain_pbc(1,s)/delta)+1
                    yi  = int(chain_pbc(2,s)/delta)+1
                    zi  = int(chain_pbc(3,s)/delta)+1

                     ! transforming form real- to lattice coordinates   
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                    indexconf(s,conf)%elem(1) = idx ! CA element

                    if(idx<=0.or.idx>nsize) then   

                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text                          
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return

                    endif

                    ! apply elements other than CA
                    do j=2,nelem(s) 

                        chain_pbc_tmp(1) = pbc(chain_elem_index(1,s)%elem(j)+xcm,Lx) ! periodic boundary conditions in x and y and z direcxtion 
                        chain_pbc_tmp(2) = pbc(chain_elem_index(2,s)%elem(j)+ycm,Ly)
                        chain_pbc_tmp(3) = pbc(chain_elem_index(3,s)%elem(j)+zcm,Lz)  

                        ! transforming form real- to lattice coordinates                 
                        xi = int(chain_pbc_tmp(1)/delta)+1
                        yi = int(chain_pbc_tmp(2)/delta)+1
                        zi = int(chain_pbc_tmp(3)/delta)+1
                        
                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                        indexconf(s,conf)%elem(j) = idx ! all other element

                        if(idx<=0.or.idx>nsize) then   

                            text="Conformation outside box:"
                            call print_to_log(LogUnit,text)  
                            print*,text                          
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            info= myio_err_index
                            return
                        endif

                    enddo
                    
                enddo

                energychain(conf)      = energy

                Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                nucl_spacing(:,conf)   = nucleosomal_spacing(chain_pbc,nnucl,segcm)
                
                conf=conf+1   
                                    
            case("prism") 
                    
                do s=1,nseg
                    
                    xp(s) = chain_rot(1,s)+xcm
                    yp(s) = chain_rot(2,s)+ycm
                    zp(s) = chain_rot(3,s)+zcm

                    ! .. transformation to prism coordinates 
                    xpp(s) = ut(xp(s),yp(s))
                    ypp(s) = vt(xp(s)+ycm,yp(s))

                    ! .. periodic boundary conditions in u and v ands z direction
                    chain_pbc(1,s) = pbc(xpp(s),Lx) 
                    chain_pbc(2,s) = pbc(ypp(s),Ly)
                    chain_pbc(3,s) = pbc(zp(s),Lz)        

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_pbc(1,s)/delta)+1
                    yi = int(chain_pbc(2,s)/delta)+1
                    zi = int(chain_pbc(3,s)/delta)+1
                        
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                    indexchain(s,conf) = idx

                    if(idx<=0.or.idx>nsize) then    
                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text                          
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return
                    endif
                    
                enddo
                
                energychain(conf)      = energy

                Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                nucl_spacing(:,conf)   = nucleosomal_spacing(chain_pbc,nnucl,segcm)

                conf=conf+1   

            case default
                text="Error: in make_chains_XYZ_nucl geometry not cubic or prism: stopping program"
                call print_to_log(LogUnit,text)
                info = myio_err_geometry
                return 
                    
            end select

        endif   ! read 

    enddo       ! end while loop                                                                                                          
    ! end chains generation    
    
    conf=conf-1  ! lower by one  

    if(conf<max_confor) then
        print*,"subroutine make_chains_xyz_nucl_volume :" 
        print*,"conf     = ",conf," less then imposed max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas=conf   
        info=myio_err_conf        
    else
        text="Chains generated: subroutine make_chains_xyz_nucl_volume"
        call print_to_log(LogUnit,text)
        readinchains=conffile
        info=0
    endif


    write(istr,'(L2)')isVdWintEne
    text="isVdWintEne = "//trim(adjustl(istr))
    call print_to_log(LogUnit,text)

    if(.not.(isChainEnergyFile)) energychain=0.0_dp

    close(un) 
    if(isChainEnergyFile) close(un_ene)


    call normed_weightchains()     

    deallocate(energychain) ! free unused variables 
    
end subroutine read_chains_xyz_nucl_volume

subroutine read_graftpts_xyz_nucl(info)

    use mpivars, only : rank
    use parameters, only : unit_conv
    use myio, only : myio_err_chainsfile, myio_err_graft
    use myutils,  only : newunit
    use chains, only : sgraftpts

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    character(len=25) :: fname
    integer :: ios 
    real(dp) :: xc,yc,zc          
    integer :: ix,iy,iz,un,s, i,t
    integer :: rankfile
    integer :: item,moltype,nsegfile,idatom
    character(len=30) :: istr,str
    real(dp) :: xbox0,xbox1,scalefactor
    logical :: exist, isGraftItem

    ! .. executable statements 
    info=0

    !. . open file    
    rankfile=rank                                                                                           
    write(istr,'(I4)')rankfile
    fname='traj-graft.'//trim(adjustl(istr))//'.xyz'
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj-graft.rank.xyz file does not exit'
        info = myio_err_chainsfile
        return
    end if
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    end if

    ! read preamble/header
   
    scalefactor=unit_conv
    isGraftItem=.false.
    do s=1,2
        read(un,*,iostat=ios)item,xc,yc,zc
        if(item==sgraftpts(1)) then
            t=1
            isGraftItem=.true.
        else
            t=2
        end if    
        xgraftloop(1,t)=xc*scalefactor
        xgraftloop(2,t)=yc*scalefactor
        xgraftloop(3,t)=zc*scalefactor    
    end do
    
    close(un)

    if(.not.isGraftItem) then
        print*,'Error read_graftpts_xyz_loop: sgraft not in traj file.'
        info = myio_err_graft
    end if
        
end subroutine


! post:  make isAmonomer, type_of_monomer and  type_of_monomer_char 
! pre: nseg, nsegtypes, freq=chainperiod and chaintype are set

subroutine make_sequence_chain(freq,chaintype)
  
    use globals, only : nseg, nsegtypes
    use chains, only :  isAmonomer,type_of_monomer_char, type_of_monomer,ismonomer_of_type
    use parameters, only : typesfname

    integer, intent(in) :: freq
    character(len=8),intent(in)  :: chaintype

    !     .. local variables
    integer :: info
    integer :: s

    select case (chaintype ) 
    case ('altA')
        do s=1,nseg
            if(mod(s,freq).ne.0) then ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo
    case('altB') 
        do s=1,nseg
            if(mod(s,freq).eq.0) then ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo
    case('diblockA')
        do s=1,nseg
            if(s<=freq) then   ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo
    case('diblockB') 
        do s=1,nseg
            if(s>=freq) then   ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo  
    case('copolyAB')

        call read_sequence_copoly_from_file(info)

         do s=1,nseg
           if(isAmonomer(s)) then 
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo

    case('multi')
    
        call read_type_of_monomer(type_of_monomer,type_of_monomer_char,typesfname, nseg) 
        call make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)    
        
        do s=1,nseg
            isAmonomer(s)=(type_of_monomer_char(s)=="AA")
        end do

    case default
        print*,"Wrong chaintype: aborting program"
        stop
    end select

end subroutine make_sequence_chain


subroutine read_sequence_copoly_from_file(info)

    use globals, only : nseg
    use chains,  only : isAmonomer
    use myutils, only : newunit

    integer, intent(out),optional :: info

    character(len=11) :: fname
    integer :: ios, un_seq, s
    character :: char

    if (present(info)) info = 0

    write(fname,'(A11)')'sequence.in'
    open(unit=newunit(un_seq),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file sequence.in : iostat =', ios
        stop
    end if
        
    s=0
    ios=0

    do while (s<nseg.and.ios==0)

        read(un_seq,*,iostat=ios)char

        if (ios==0) then
            s=s+1 
            if(char=="A") then   ! A segment
                isAmonomer(s)=.TRUE.
            else
                isAmonomer(s)=.FALSE.
            end if
        end if    
    end do 
    if(s/=nseg) then 
        print*,"reached end of file before all elements read"
        info = 1
        stop "read sequence file failed"
    end if

end subroutine read_sequence_copoly_from_file


! Compute weight chain w=e^E/Tre^E and normalize
! layout conformations : there are  nset of conformations thus total nodes = nset 
! Each set hold different conformations

subroutine normed_weightchains()

    use mpivars
    use globals, only : cuantas
    use chains, only : energychain, logweightchain
   
    integer :: un, c, k
    real(dp) :: localsum, totalsum, logtotalsum

        
    localsum=0.0_dp    

    do c=1,cuantas
        localsum=localsum+exp(energychain(c))   
    enddo    
  
    call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
    call MPI_ALLREDUCE(localsum, totalsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr) 
    
    ! normalize
    logtotalsum=log(totalsum)
    do c=1,cuantas
        logweightchain(c)=energychain(c)-logtotalsum
    enddo    

end subroutine



function minimum_chainenergy() result(min_chainenergy)

    use  mpivars
    use  globals, only :cuantas
    use  chains, only : energychain

    real(dp) :: min_chainenergy

    integer :: conf
    
    min_chainenergy =0.0_dp
    do conf=1,cuantas      
        if(min_chainenergy > energychain(conf)) min_chainenergy=energychain(conf)
    enddo
   
end function


subroutine global_minimum_chainenergy()

    use  mpivars
    use  globals, only : cuantas
    use  chains, only : energychain, energychain_min
    use  parameters, only: isEnergyShift

    real(dp) :: localmin(2), globalmin(2)
    integer :: i

    localmin(1)=minimum_chainenergy()
    localmin(2)=rank   
    print*,"rank ",rank ," has local lowest value of", localmin(1)
    
    call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
    call MPI_ALLREDUCE(localmin, globalmin, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
    
    if (rank == 0) then  
        print*,"Rank ",globalmin(2)," has lowest value of", globalmin(1)  
    endif
    
    energychain_min =globalmin(1)

    if(isEnergyShift) then        
        call shift_energy_chain(globalmin(1))
        print*,"Shift minimum by :",energychain_min 
    endif
end subroutine



! post: substract minimum internal energy from all internal chain energies
subroutine shift_energy_chain(min_chainenergy)
    
    use  globals, only :cuantas
    use  chains, only : energychain
    
    real(dp), intent(in) :: min_chainenergy

    ! local
    integer :: conf

    do conf=1,cuantas        
        energychain(conf)=energychain(conf)-min_chainenergy
    enddo

end subroutine shift_energy_chain


subroutine set_properties_chain(freq,chaintype)
    
    use globals
    use chains
    use parameters, only : lseg

    integer, intent(in) :: freq
    character(len=8), intent(in)   :: chaintype

    call set_isHomopolymer(freq,chaintype)
    call set_lsegAA()

    if(isHomopolymer) lseg=set_lseg_Homopolymer()

end subroutine set_properties_chain

! set logical isHomoolymer 
! pre : type_of_monomer need to be initialized

subroutine set_isHomopolymer(freq,chaintype)

    use globals, only : nseg
    use chains, only : type_of_monomer, isHomopolymer

    integer, intent(in) :: freq
    character(len=8), intent(in)   :: chaintype
    
    integer :: s,  type_number

    if(chaintype=="multi") then
        type_number=type_of_monomer(1)
        do s=2,nseg
            isHomopolymer=(type_of_monomer(s)==type_number)
        enddo
    else if(freq>nseg) then 
        isHomopolymer=.TRUE.
    else 
        if(freq==0.and.(chaintype.eq."diblockA".or. chaintype=="diblockB")) then 
            isHomopolymer=.TRUE.
        else
            isHomopolymer=.FALSE.
        endif      
    endif    

end subroutine set_isHomopolymer


!  Assigns specific values to lsegAA initial set in init_lseg in parameters

subroutine set_lsegAA

    use globals
    use chains
    use parameters, only : lseg, lsegAA,lsegPAA, lsegPAMPS, lsegPEG
    use parameters, only : chainmethod 

    if(chainmethod=='MC') then  ! chain are not read in from file 
        
        select case (systype)
        case ("elect","electVdWAB","electdouble") ! diblock copolymer lseg determined in cadenas_sequence      
            ! assume A is present and type_of_monomer =1 
            lsegAA(1) = lsegPAA
            lsegAA(2) = lsegPAMPS        
        case ("electA")                 ! homopolymer weak polyacid VdW     
            lsegAA = lsegPAA
        case ("neutral","neutralnoVdW")                ! homopolymer neutral
            !lsegAA = lsegPEG
        case ("brush_mul","brush_mulnoVdW","brush","brush_neq","brushvarelec","brushborn","brushdna")
            ! lsegAA = lsegPAA !0.36287_dp
        case default
            print*,"Error: in set_lsegAA, systype=",systype
            print*,"stopping program"
            stop
        end select  
    endif   

end subroutine set_lsegAA

     
function set_lseg_Homopolymer()result(lengthseg)

    use parameters, only : lsegAA

    real(dp) :: lengthseg
    
    lengthseg=lsegAA(1)
    
end function set_lseg_Homopolymer



! .. compute periodic boundary condition in z direction
! .. maps z onto interval [0,deltaz]

function pbc(z,deltaz) result(zpbc)
    implicit none 
    real(dp), intent(in) :: deltaz
    real(dp), intent(in) :: z
    real(dp) :: zpbc 
    if(z>0) then 
        zpbc=z-int(z/deltaz)*deltaz
    else
        zpbc=z-(int(z/deltaz)-1)*deltaz
    endif   
end function


! function ipbc(ival,imax) result(intpbc)
!     implicit none 
!     integer, intent(in) :: ival
!     integer, intent(in) :: imax
!     integer :: intpbc

!     if(ival>0) then
!         intpbc=ival-int((ival-1)/imax)*imax
!     else
!         intpbc=ival-(int((ival-1)/imax)-1)*imax
!     endif

! end function


!     .. returns number of A monomer
!     .. this works only for binary AB polymers

function number_Amonomers(nseg)result(npolA)
  
    use chains, only : isAmonomer

    integer , intent(in) :: nseg
    integer :: npolA

    !     .. local variables
    integer :: s

    npolA=0
    do s=1,nseg
        if(isAmonomer(s).eqv..true.) npolA=npolA+1
    enddo

end function number_Amonomers



!  .. assign type_of_monomer and type_of_monomer_char from values in file named filename

subroutine read_type_of_monomer(type_of_monomer, type_of_monomer_char,filename, nseg)

    use mpivars
    use myutils, only : newunit
    use parameters, only : lenfname
    
    !     .. arguments 
    integer, intent(inout) :: type_of_monomer(:)
    character(len=3), intent(inout) :: type_of_monomer_char(:)
    character(lenfname), intent(in) :: filename
    integer,  intent(in) :: nseg

    !      .. local variables
    integer :: ios, un  ! un = unit number
    integer :: s
    character(80) :: istr,str,letter

    !     .. reading in of variables from file
    open(unit=newunit(un),file=filename,iostat=ios,status='old')
    if(ios/=0 ) then
        write(istr,'(I2)')ios
        str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
        print*,str
        call MPI_FINALIZE(ierr)
        stop
    endif
    
    s=0
    ios=0

    do while (s<nseg.and.ios==0)
        s=s+1
        read(un,*,iostat=ios)type_of_monomer(s),type_of_monomer_char(s)
    enddo

    if(s/=nseg.or.ios/=0) then 
        str="reached end of file before all elements read or ios error"
        print*,str
        str="read file "//trim(adjustl(filename))//" failed"
        print*,str
        call MPI_FINALIZE(ierr)
        stop
    endif

    close(un)

end subroutine read_type_of_monomer



! ismonomer_of_type is a table which row index is the segment  and column index correspond to the segment type  
! pre: type_of_monomer needs to be initialized
! post: table list of logicals indication is segement s is of type t

subroutine make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)
    
    logical, intent(out) :: ismonomer_of_type(:,:)
    integer, intent(in) :: type_of_monomer(:) 
    integer, intent(in) :: nseg
    integer, intent(in) :: nsegtypes

    ! local variable 
    integer :: s,t 

    do s=1,nseg 
        do t=1,nsegtypes   
            ismonomer_of_type(s,t)=.false.
        enddo
        ismonomer_of_type(s,type_of_monomer(s))=.true.    
    enddo

end subroutine make_type_table


! routine determines is segment type t is chargeable
! pre: zpol needs to be initialized
! post: ismonomer_chargeable list of logicals

subroutine make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
 
    logical, intent(out) :: ismonomer_chargeable(:)
    integer, intent(in) :: zpol(:,:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t 

    do t=1,nsegtypes
        if(zpol(t,1)==0.and.zpol(t,2)==0) then 
            ismonomer_chargeable(t)=.false.
        else
            ismonomer_chargeable(t)=.true.
        endif
    enddo

end subroutine make_charge_table


! routine determines is type of charge of segment type t 
! pre: zpol needs to be initialized
! post: type_of_charge list of characters

subroutine make_type_of_charge_table(type_of_charge,zpol,nsegtypes)
 
    character(len=1), intent(out) :: type_of_charge(:)
    integer, intent(in) :: zpol(:,:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t 

    do t=1,nsegtypes
        if(zpol(t,1)==0.and.zpol(t,2)==0) then 
            type_of_charge(t)="N"
        else if(zpol(t,1)==0.and.zpol(t,2)==-1) then
             type_of_charge(t)="A"
        else if(zpol(t,1)==1.and.zpol(t,2)==0) then
             type_of_charge(t)="B"  
        endif
    enddo

end subroutine make_type_of_charge_table



logical function is_polymer_neutral(ismonomer_chargeable, nsegtypes)
 
    logical, intent(in) :: ismonomer_chargeable(:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t,flag 

    flag=0
    do t=1,nsegtypes
        if(.not.ismonomer_chargeable(t)) flag=flag+1
    enddo

    if(flag==nsegtypes) then 
        is_polymer_neutral=.true. 
    else
        is_polymer_neutral=.false. 
    endif

end function is_polymer_neutral


! make a lammps trajectory file based on  indexchain

subroutine write_indexchain_lammps_trj(info)

    use mpivars, only : rank
    use globals, only : nseg, cuantas
    use myutils, only : newunit, lenText
    use volume, only : delta, nz, coordinateFromLinearIndex
    use chains, only : indexchain, type_of_monomer
    use myio, only : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    character(len=lenText) :: text, istr
    character(len=25) :: fname
    integer :: ios, un_trj 
    real(dp):: xs, ys, zs
    integer :: ix, iy, iz, i, j, k, idx
    real(dp) :: xbox0, xbox1
    integer :: idatom, item, moltype, conf

    ! open file 

    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.lammpstrj'
    open(unit=newunit(un_trj),file=fname,status='new',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

    ix=0
    iy=0
    iz=0
    xbox0=0.0_dp
    xbox1=nz*delta
    
    do conf=1,cuantas
        ! write preamble 
        write(un_trj,*)'ITEM: TIMESTEP' 
        write(un_trj,*)conf 
        write(un_trj,*)'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nseg 
        write(un_trj,*)'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)'ITEM: ATOMS id mol type xs ys zs ix iy iz'
        !  determine xs, ys, zs,
        do item=1,nseg
            idx=indexchain(item,conf)
            call coordinateFromLinearIndex(idx, i, j, k)
            xs=(i-0.5_dp)*delta/xbox1
            ys=(j-0.5_dp)*delta/xbox1
            zs=(k-0.5_dp)*delta/xbox1
            idatom=1
            moltype=type_of_monomer(item) 
            write(un_trj,*)item,idatom,moltype,xs,ys,zs,ix,iy,iz
        enddo    
    enddo
        
    close(un_trj)    


end subroutine write_indexchain_lammps_trj


function open_chain_lammps_trj(info)result(un_trj)

    use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    integer :: un_trj
    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios
    logical :: exist
    
    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.lammpstrj'
    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un_trj),file=fname,status='new',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un_trj),file=fname,status='old',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    endif    
end function


subroutine write_chain_lammps_trj(un_trj,chain,nchains)

    use globals, only : nseg
    use myutils, only : newunit, lenText
    use volume, only : delta, nz, coordinateFromLinearIndex
    use chains, only :  type_of_monomer

    real(dp), intent(in) :: chain(:,:,:)
    integer, intent(in) :: nchains
    integer , intent(in) :: un_trj
    
    character(len=lenText) :: istr
    real(dp) :: xs, ys, zs
    integer :: ix, iy, iz, i, j, k
    real(dp) :: xbox0, xbox1
    integer :: idatom, item, moltype, conf

    ix=0
    iy=0
    iz=0
    xbox0=0.0_dp
    xbox1=nz*delta
    do j=1,nchains   
        
        conf_write=conf_write+1 
        ! write preamble 
        write(un_trj,*)'ITEM: TIMESTEP' 
        write(un_trj,*)conf_write 
        write(un_trj,*)'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nseg 
        write(un_trj,*)'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)'ITEM: ATOMS id mol type xs ys zs ix iy iz'
        do item=1,nseg
            zs = chain(1,item,j)/xbox1
            xs = chain(2,item,j)/xbox1
            ys = chain(3,item,j)/xbox1 
            idatom=1
            moltype=type_of_monomer(item) 
            write(un_trj,*)item,idatom,moltype,xs,ys,zs,ix,iy,iz
        enddo    
    enddo
        
end subroutine write_chain_lammps_trj


function open_chain_energy(info)result(un_ene)

    use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    integer :: un_ene
    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios

    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.ene'
    open(unit=newunit(un_ene),file=fname,status='replace',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

end function

function VdWpotentialenergy_MC(chain,nchains)result(Energy)

    use globals, only : nseg
    use myutils, only : newunit, lenText
    use volume, only : delta, nx, ny, nz, coordinateFromLinearIndex
    use chains, only : type_of_monomer
    use parameters, only :  lsegAA,VdWeps

    real(dp), intent(in) :: chain(:,:,:)
    integer, intent(in) :: nchains

    real(dp) :: Energy(nchains)
    real(dp) :: Ene,sqrlseg,sqrdist
    integer :: k,i,j,s,t
    real(dp) :: xi,xj,yi,yj,zi,zj
    real(dp) :: Lz,Ly,Lx

    Lz= nz*delta            ! maximum height box s
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 

    do k=1,nchains
        Ene=0.0_dp      
        do i=1,nseg
            do j=i+1,nseg

                s=type_of_monomer(i)
                t=type_of_monomer(j)
               
                zi = pbc(chain(1,i,k),Lz)
                xi = pbc(chain(2,i,k),Lx)
                yi = pbc(chain(3,i,k),Ly) 
                zj = pbc(chain(1,j,k),Lz)
                xj = pbc(chain(2,j,k),Lx)
                yj = pbc(chain(3,j,k),Ly) 

                sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

                Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3

            enddo
        enddo 
        Energy(k) = Ene            
    enddo

end function VdWpotentialenergy_MC
 
! pre : chain conformation 
! post: VdW of conformation and if conformation is saw or not

subroutine VdWpotentialenergySaw(chain,Energy,saw)

    use globals, only : nseg, nsegtypes
    use chains, only : type_of_monomer
    use parameters, only :  lsegAA,VdWeps

    real(dp), intent(in)  :: chain(:,:)
    real(dp), intent(out) :: Energy
    logical, intent(out) :: saw

    real(dp) :: Ene,sqrlseg,sqrdist
    integer :: i,j,s,t
    real(dp) :: xi,xj,yi,yj,zi,zj
   
    Ene=0.0_dp      
    saw=.true.

    do i=1,nseg
        do j=i+1,nseg

            s=type_of_monomer(i)
            t=type_of_monomer(j)

            zi = chain(1,i)
            xi = chain(2,i)
            yi = chain(3,i)

            zj = chain(1,j)
            xj = chain(2,j)
            yj = chain(3,j) 

            sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
            sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

            Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3
            
            if(j/=(i+1).and.sqrdist<sqrlseg) then 
                saw=.false.
             !   print*,"overlap occured for i= ",i," and j= ",j 
            endif    
        enddo
    enddo 

    !print*,"total ",Ene," saw ",saw,' ',(lsegAA(t),t=1,nsegtypes),VdWeps 
    Energy = Ene            

end subroutine VdWpotentialenergySaw
   

 
! pre : chain conformation 
! post: VdW of conformation 
! Warning: need to to add VdWcutoff

subroutine VdWpotentialenergy(chain,Energy)

    use globals, only : nseg, nsegtypes
    use chains, only : type_of_monomer
    use parameters, only :  lsegAA,VdWeps

    real(dp), intent(in)  :: chain(:,:)
    real(dp), intent(out) :: Energy
    
    real(dp) :: Ene,sqrlseg,sqrdist !,maxdist
    integer ::  i,j,s,t
    real(dp) :: xi,xj,yi,yj,zi,zj
    
    Ene=0.0_dp 

    !maxdist=VdWcutoff*delta 
   
    do i=1,nseg
        do j=i+1,nseg

            s=type_of_monomer(i)
            t=type_of_monomer(j)

            zi = chain(1,i)
            xi = chain(2,i)
            yi = chain(3,i)

            zj = chain(1,j)
            xj = chain(2,j)
            yj = chain(3,j) 

            sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
            sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

            if(sqrdist<sqrlseg) sqrdist=sqrlseg

            Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3
            
        enddo
    enddo 

    Energy = Ene            

end subroutine VdWpotentialenergy


! .. commuter radius of gyration of a sub chain conformation
! .. sub chain confomation defined by sequence segcm = segment number 
! .. denoting center mass of  the nmer histone of the nucleosome nmer-array

function radius_gyration_com(chain,nmer,segcm) result(Rgsqr)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcm(:) 

    real(dp) :: Rgsqr
    integer :: i,j,k,isegcm,jsegcm

    Rgsqr=0.0_dp
    
    do i=1,nmer
        isegcm=segcm(i)
        do j=1,nmer
            jsegcm=segcm(j)
            do k=1,3
                Rgsqr=Rgsqr+(chain(k,isegcm)-chain(k,jsegcm))**2
            enddo
        enddo
    enddo     

    Rgsqr=Rgsqr/(2.0_dp*nmer*nmer)

end function radius_gyration_com


function radius_gyration(chain,nseg) result(Rgsqr)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nseg
    real(dp) :: Rgsqr
    integer :: i,j,k

    Rgsqr=0.0_dp
    
    do i=1,nseg
        do j=1,nseg
            do k=1,3
                Rgsqr=Rgsqr+(chain(k,i)-chain(k,j))**2
            enddo
        enddo
    enddo     

    Rgsqr=Rgsqr/(2.0_dp*nseg*nseg)

end function radius_gyration


function end_to_end_distance_com(chain,nmer,segcm) result(Rendsqr)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcm(:)  

    real(dp) :: Rendsqr
    integer :: isegcm, jsegcm, k

    Rendsqr=0.0_dp
       
    isegcm=segcm(1)
    jsegcm=segcm(nmer)

    do k=1,3
        Rendsqr=Rendsqr+(chain(k,isegcm)-chain(k,jsegcm))**2
    enddo     

end function end_to_end_distance_com

function nucleosomal_spacing(chain,nmer,segcom) result(spacing)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcom(:)  

    real(dp) :: spacing(nmer-1)
    integer :: i,k
    real(dp) :: spacingsqr
        
    spacing=0.0_dp
    if(nmer>=2) then  ! need at least 2 unit/nucleosomes
        do i=1,nmer-1
            spacingsqr=0.0_dp
            do k=1,3
                spacingsqr=spacingsqr+(chain(k,segcom(i+1))-chain(k,segcom(i)))**2
            enddo
            spacing(i)=sqrt(spacingsqr)
        enddo  
    endif           

end function nucleosomal_spacing

function bond_angles_com(chain,nmer,segcom) result(bondangle)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcom(:)
    real(dp) :: bondangle(nmer-2)

    ! .. local variables
    real(dp) :: u1(3), u2(3), absu1, absu2
    integer :: i,j,k ,isegcom, ipls1segcom


    bondangle=0.0_dp    

    if(nmer>=3) then  ! need at least 3 unit/nucleosomes
        isegcom=segcom(1)
        ipls1segcom=segcom(2)
        do k=1,3
            u1(k)=(chain(k,ipls1segcom)-chain(k,isegcom))
        enddo

        absu1=sqrt(dot_product(u1,u1))

        do i=3,nmer
            isegcom=segcom(i-1)
            ipls1segcom=segcom(i)
            do k=1,3
                u2(k)=(chain(k,ipls1segcom)-chain(k,isegcom))
            enddo
            absu2=sqrt(dot_product(u2,u2))
            bondangle(i-2)=acos(dot_product(u1,u2)/(absu2*absu1))
            u1=u2
            absu1=absu2
        enddo
    endif         

end function bond_angles_com


function dihedral_angles_com(chain,nmer,segcm) result(dihedral)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcm(:) 
    
    real(dp) :: dihedral(nmer-3)
    
    ! .. local variables
    integer  :: i,j,k 
    real(dp) :: u1(3),u2(3),u3(3),n123(3),n234(3)
    real(dp) :: absn123, absn234, absu2
    integer  :: isegcm,ipls1segcm,ipls2segcm,ipls3segcm
    real(dp) :: theta(nmer-3), costheta, x, y, sintheta, sintheta2 ,theta2,sintheta3, sintheta4
    


    if(nmer>=4) then  ! need at least 4 unit/nucleosomes

        isegcm    =segcm(1)
        ipls1segcm=segcm(2)
        ipls2segcm=segcm(3)
        ipls3segcm=segcm(4)

        do k=1,3
            u1(k)=(chain(k,ipls1segcm)-chain(k,isegcm    ))
            u2(k)=(chain(k,ipls2segcm)-chain(k,ipls1segcm))
            u3(k)=(chain(k,ipls3segcm)-chain(k,ipls2segcm))
        enddo

        n123=crossproduct(u1,u2)
        n234=crossproduct(u2,u3)

        absn123=sqrt(dot_product(n123,n123))  ! or sqrt(sum(n123*n123)) 
        absn234=sqrt(dot_product(n234,n234))
        absu2=sqrt(dot_product(u2,u2))
   
        y = dot_product(u2,crossproduct(n123,n234)) 
        x = absu2*dot_product(n123,n234)   
            
        theta(1)=atan2(y,x)  !formula See e.g. Bondel and Karplus j. comp. chem. 17 1132

        costheta= dot_product(n123,n234)/(absn123*absn234) ! cos first dihedral angle
        sintheta= dot_product(u2,crossproduct(n123,n234))/(absu2*absn123*absn234)

       ! following dihedral angles 

        do i=5,nmer 
        
            ! advance one unit
            n123=n234 
            absn123=absn234 
            u1=u2
            u2=u3 

            ipls2segcm=ipls3segcm
            ipls3segcm=segcm(i) 

            do k=1,3
                u3(k)=(chain(k,ipls3segcm)-chain(k,ipls2segcm))
            enddo

            n234=crossproduct(u2,u3)
            absn234=sqrt(dot_product(n234,n234))
            absu2=sqrt(dot_product(u2,u2))
           
            y = dot_product(u2,crossproduct(n123,n234)) 
            x = absu2*dot_product(n123,n234)
            
            theta(i-3)=atan2(y,x) ! i-3 dihedral angle

        enddo    
   
    endif   

    dihedral=theta    

end function dihedral_angles_com



! inner product of vectors a and b with dimenstion 3
! similar to intrisic function dotproduct

function dotproduct(a,b) result(dotprod)
   
    real(dp), dimension(3) :: a, b 
    real(dp) :: dotprod
  
    dotprod = sum(a*b)

end function dotproduct

function crossproduct(a,b)result(crossprod)
  
  real(dp), intent(in), dimension(3) :: a, b
  real(dp), dimension(3) :: crossprod
  
  crossprod(1) = a(2)*b(3) - a(3)*b(2)
  crossprod(2) =-a(1)*b(3) + a(3)*b(1)
  crossprod(3) = a(1)*b(2) - a(2)*b(1)

end function crossproduct



!  .. read segment id number to which we assign center of mass of one histone from file
!  .. segment id numbers stored in array segcm 
!  .. nucleosome has nnucl histones 
!  .. post segcom assigned

subroutine make_segcom(segcom,nnucl,filename)

    use mpivars
    use myutils
    use parameters, only : lenfname
    
    !     .. arguments 
    integer, intent(inout) :: segcom(:)
    character(lenfname), intent(in) :: filename
    integer,  intent(in) :: nnucl

    !      .. local variables
    integer :: ios, un  ! un = unit number
    integer :: s
    character(80) :: istr,str,letter

    !     .. reading in of variables from file
    open(unit=newunit(un),file=filename,iostat=ios,status='old')
    if(ios/=0 ) then
        write(istr,'(I2)')ios
        str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
        print*,str
        call MPI_FINALIZE(ierr)
        stop
    endif
    
    s=0
    ios=0

    do while (s<nnucl.and.ios==0)
        s=s+1
        read(un,*,iostat=ios)segcom(s)
    enddo

    if(s/=nnucl.or.ios/=0) then 
        str="reached end of file before all elements read or ios error"
        print*,str
        str="read file "//trim(adjustl(filename))//" failed"
        print*,str
        call MPI_FINALIZE(ierr)
        stop
    endif

    close(un)

end subroutine make_segcom



subroutine write_chain_struct(write_struct,info)

    use globals, only : cuantas,nnucl
    use myutils, only : lenText
    use chains, only : Rgsqr,Rendsqr,bond_angle,dihedral_angle,nucl_spacing

    implicit none 

    logical, intent(in) :: write_struct
    integer, intent(inout) :: info
 
    ! .. local
    character(len=lenText) :: filename
    integer :: un_dihedral,un_bond,un_Rg,un_Rend, un_dist
    integer :: c,s,nbonds,nangles,ndihedrals

    info=0

    if(write_struct) then
    
        filename="dihedral."
        un_dihedral=open_chain_struct_file(filename,info)
        filename="bond."
        un_bond=open_chain_struct_file(filename,info)
        filename="Rgalpha."
        un_Rg=open_chain_struct_file(filename,info)
        filename="Rend."
        un_Rend=open_chain_struct_file(filename,info)
        filename="spacing."
        un_dist=open_chain_struct_file(filename,info)
        
        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
    
        do c=1,cuantas
            write(un_Rg,*)Rgsqr(c)
            write(un_Rend,*)Rendsqr(c)
            write(un_bond,*)(bond_angle(s,c),s=1,nangles)
            write(un_dihedral,*)(dihedral_angle(s,c),s=1,ndihedrals)
            write(un_dist,*)(nucl_spacing(s,c),s=1,nbonds)
        enddo 

        close(un_dihedral)
        close(un_bond)
        close(un_Rg)
        close(un_Rend)  
        close(un_dist)
      
    endif

        
end subroutine write_chain_struct


function open_chain_struct_file(filename,info)result(un)

    use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    character(len=lenText) :: filename
    integer, optional, intent(inout) :: info

    integer :: un

    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios
    logical :: exist
    
    if (present(info))  info=0 ! init
    
    write(istr,'(I4)')rank
    fname=trim(adjustl(filename))//trim(adjustl(istr))//'.dat'

    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un),file=fname,status='new',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    endif    
end function open_chain_struct_file


subroutine set_mapping_num_to_char(type_of_monomer_num_to_char)

    use globals, only : nsegtypes,nseg
    use chains, only : type_of_monomer, type_of_monomer_char

    character(len=3), intent(inout) :: type_of_monomer_num_to_char(:)  
        
    integer :: t, k
    logical :: notfound

    do t=1,nsegtypes
        k=1
        notfound=.true.
        do while (k<=nseg .and. notfound)
            if(type_of_monomer(k)==t) then
                type_of_monomer_num_to_char(t)=type_of_monomer_char(k)
                notfound=.false. ! stops
                k=1 ! reset
            endif 
            k=k+1
        enddo        
    enddo

end subroutine set_mapping_num_to_char

! reads in MT pdb file to assign 
! input  : fname   = char(*)  filename of MTpdb file ) 
!        : nsegAA  = integer = number of AA residues
!        : typesAA = array of integer:  the s element gives number of AA type
!        : nelemAA = array of integer:  teh s element  gives number of element the sth AA consits of
! output : chain_elem(3,s)%elem(j) = position of jthe element of AA number s. Enumarate the AA in order
!        : vnucl(j,t) = real(dp) volume of AA  element j of AA type number t
!        : nsegAA   = integer = number of AAs
!        : nelem(s) = array of integer number of elements of AA number s
!        : info     = integer val=0 assignment  succesfull val/=0 error              

subroutine read_nucl_elements(fname,nsegAA,nelemAA,chain_elem,typeAA,vnucl,elem_charge,info)

    use globals, only : DEBUG, nsegtypes,nsegtypesAA
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile
    use parameters, only : lenfname, vpol, vsol
    use chains, only : var_darray, ismonomer_chargeable, mapping_num_to_char

    character(lenfname),   intent(in)                            :: fname
    integer,               intent(in)                            :: nsegAA
    integer, dimension(:), intent(in)                            :: typeAA
    integer, dimension(:), intent(inout)                         :: nelemAA
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem
    real(dp), dimension(:,:) , intent(inout)                     :: vnucl
    integer, dimension(:), allocatable, intent(inout)            :: elem_charge
    integer, intent(out),optional                                :: info 

    ! local arguments

    integer   :: ios, un, s, j, k, ttype ,t
    character(len=3) :: elem_type(20), voltype
    real(dp)  :: x(3)
    logical   :: isReadGood
    integer   :: num_elem_CA
    character(len=3) :: list_type_char_DNA(6)  ! list of DNA elements in char
    integer          :: list_type_int_DNA(6)   ! list of DNA elements in integer 
    integer          ::  AAid, DNAid
    character(len=3) :: type_char  

    if (present(info)) info = 0

    open(unit=newunit(un),file=fname,status='old',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

    allocate(chain_elem(3,nsegAA))

    do s=1,nsegAA

        read(un,*,iostat=ios) ! comment
        if(ios/=0) isReadGood=.false.
        read(un,*,iostat=ios) AAid ! type number
        if(ios/=0) isReadGood=.false.
        read(un,*,iostat=ios) nelemAA(s)
        if(ios/=0) isReadGood=.false.
        
        do k=1,3
            allocate(chain_elem(k,s)%elem(nelemAA(s)))
        enddo     
        
        elem_type=""
        do j=1,nelemAA(s)    
            
            read(un,*,iostat=ios)elem_type(j),x(1),x(2),x(3) 
            
            do k=1,3
                chain_elem(k,s)%elem(j)=x(k)
            enddo  
            
            vnucl(j,AAid)=find_vol_elem(elem_type(j))  
            ! assignment doubles since multiple s (AAs) corresponds to same type ( AAid) 
            
            if(ios/=0) isReadGood=.false.

        enddo 

        num_elem_CA=find_CA_elem(elem_type,nelemAA(s))

    
        print*,"AAId=", AAid," chargeable ",ismonomer_chargeable(AAid)

        if(ismonomer_chargeable(AAid )) then ! ismonomer_chargeable has type number as input !! 
            elem_charge(s)=find_chargeable_elem(elem_type,nelemAA(s))
        else 
            elem_charge(s)=-1
        endif      
       
        ! place CA element in to first position j=1
       
        ! print*,"elem_type=",elem_type
        ! print*,"num_elem_CA=",num_elem_CA

        if(num_elem_CA/=1) then 
            ! swapping
            call swap_char3(elem_type(num_elem_CA),elem_type(1))
            !call swap_real(vnucl(num_elem_CA,s), vnucl(1,s))
            call swap_real(vnucl(num_elem_CA,AAid), vnucl(1,AAid))

            do k=1,3
                call swap_real(chain_elem(k,s)%elem(num_elem_CA),chain_elem(k,s)%elem(1))
            enddo
        endif  

        num_elem_CA=find_CA_elem(elem_type,nelemAA(s))

        if(DEBUG) then 
            print*,"elem_type = ",elem_type
            print*,"num_elem_CA = ",num_elem_CA
        endif

    enddo
    
    ! .. assign volume to vnucl that is from  DNA elements 
    ! .. above loops ssign vnucl only includes AA that are in the MT pdb file  
    
    list_type_char_DNA = (/"P  ","S  ","A  ","C  ","G  ","T  "/)
    list_type_int_DNA  = 0

    ! .. make list_type_int_DNA"

    do k=1,6                                ! 6 = number DNA elements 
        type_char= list_type_char_DNA(k)
        do t=1,nsegtypes
            if( mapping_num_to_char(t) == type_char)  then 
                list_type_int_DNA(k)=t
            endif
        enddo 
    enddo

    do k=1,6
        DNAid=list_type_int_DNA(k)
        print*,"DNAid=",DNAid
        if(DNAid/=0) then  ! in case the conformation is incomplete adn does not 
            vnucl(1,DNAid) = vpol(DNAid)*vsol
        endif
    enddo 

    ! .. end assignment DNA volume 

    if(DEBUG) then 
        print*,"DEBUG: output of read_nucl_elements"
        do s=1,nsegAA
            do k=1,3
                do j=1,nelemAA(s)
                    print*,"s=",s,"k=",k,"j=",j," ",chain_elem(k,s)%elem(j)  
                enddo 
            enddo 
        enddo
    endif     
        
    if(DEBUG) then 
        print*,"DEBUG: output of read_nucl_elements volumes"
        do t=1,nsegtypes
            do j=1,nsegtypesAA
                print*,"j=",j," t=",t," ",vnucl(j,t)
            enddo 
        enddo
    endif       



end subroutine read_nucl_elements

! Shift chain_elem(k,s)%elem(j) by position of the CA of AA that is the origin of the "CM" rotation triplet
! translation vector is chain_CA(k)=chain_elem(k,segnumcm)%elem(1)

subroutine shift_nucl_elements(nsegAA,nelemAA,segnumcm,chain_elem)
    
    use chains, only : var_darray ! type def  

    integer, intent(in)                  :: nsegAA
    integer, dimension(:), intent(in)    :: nelemAA
    integer, intent(in)                  :: segnumcm
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem

    integer :: s,k,j,i
    real(dp) :: chain_CA(3)

    print*,"segnumcm=",segnumcm
 
    ! translation vector
    do i=1,3 
        chain_CA(i)=chain_elem(i,segnumcm)%elem(1) ! 1= CA
    enddo     
    ! substract 
    do s=1,nsegAA
        do j=1,nelemAA(s)
            do k=1,3
                chain_elem(k,s)%elem(j)=chain_elem(k,s)%elem(j) - chain_CA(k)  
            enddo
        enddo        
   enddo     

end subroutine shift_nucl_elements


function find_CA_elem(elem_type,nelem)result(elem)

    character(len=3), intent(in) :: elem_type(:)
    integer, intent(in) :: nelem

    integer :: elem, j 
    logical :: IsNotFound

    IsNotFound=.true.
    j=0
    do while (IsNotFound.and.j<=nelem)
        j=j+1
        if(elem_type(j)=="CA") IsNotFound=.false.
    enddo 

    elem=j 
    if(j==0) print*,"error No CA found "

end function find_CA_elem



! Assign/finds volume that matches with elem_type
! return : vol_elem
!          if no match found vol_elem is set to 0   

function find_vol_elem(elem_type)result(vol_elem)

    use parameters, only : vnucl_type, vnucl_type_char

    character(len=3), intent(in) :: elem_type

    real(dp):: vol_elem

    integer :: elem, j , nelem_types
    logical :: IsNotFound

    nelem_types=size(vnucl_type)
    IsNotFound=.true.
    j=0

    do while (IsNotFound.and.j<=nelem_types)
        j=j+1
        ! print*,j,vnucl_type_char(j),elem_type,vnucl_type_char(j)==elem_type
        if(vnucl_type_char(j)==elem_type) IsNotFound=.false.
    enddo 

    if(j/=0) then
        vol_elem=vnucl_type(j)
    else 
        vol_elem=0.0_dp
        print*,"error vol nucl not found "
    endif    

end function find_vol_elem

function find_chargeable_elem(elem_type,nelem)result(elem)

    use parameters, only : vnucl_type, vnucl_type_char,vnucl_type_isChargeable

    character(len=3), intent(in) :: elem_type(:)
    integer, intent(in) :: nelem

    integer :: elem, j, k, ncharge, num_elem_charge, nelem_types 
    logical :: IsNotFound

    nelem_types=size(vnucl_type)
    ncharge=0
    do k=1,nelem
        IsNotFound=.true.
        j=0
        do while (IsNotFound.and.j<=nelem_types)
            j=j+1
            if(elem_type(k)==vnucl_type_char(j)) IsNotFound=.false.
        enddo 
        if(vnucl_type_isChargeable(j)) then  ! j element is charged 
            ncharge=ncharge+1  
            num_elem_charge=j 
       endif
    enddo
    if(ncharge==1) then ! oke should be only one charge
        elem = num_elem_charge
    else if (ncharge >1 ) then 
        elem =  -1
        print*,"Error found more then one element AA charged"
    else 
        elem = -1
        print*,"Error found no element AA charged"
    endif    

end function find_chargeable_elem



subroutine swap_int(a,b)
    
    integer, intent(inout) :: a, b
    integer :: tmp 

    tmp = a
    a = b 
    b = tmp 

end subroutine  swap_int   


subroutine swap_real(a,b)
    
    real(dp), intent(inout) :: a, b
    real(dp) :: tmp 

    tmp = a
    a = b 
    b = tmp 
    
end subroutine 

subroutine swap_char3(a,b)
    
    character(len=3), intent(inout) :: a, b
    character(len=3) :: tmp 

    tmp = a
    a = b 
    b = tmp 

end subroutine  swap_char3



subroutine read_nucl_orient_triplets(fname,nnucl,nuc_orient_triplet,info)

    use myutils, only : newunit, lenText
    use parameters, only : lenfname
    use myio, only : myio_err_chainsfile, myio_err_readfile 

    character(lenfname),    intent(in)      :: fname
    integer, intent(in)                    :: nnucl
    integer, dimension(:,:), intent(inout) :: nuc_orient_triplet
    integer, intent(out), optional         :: info

    ! local arguments

    integer  :: ios, un, s,  i
    integer :: ix(3)
    logical :: isReadGood

    if (present(info)) info = 0
    isReadGood=.true.

    open(unit=newunit(un),file=fname,status='old',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_readfile
        return
    endif

    do i=1,nnucl
        read(un,*,iostat=ios)ix(1),ix(2),ix(3) 
        if(ios/=0) isReadGood=.false.
        nuc_orient_triplet(i,:)=ix     
    enddo
 
    if(isReadGood.eqv..False. ) then
        print*, 'Error reading : ',fname,' file : iostat =', ios
        if (present(info)) info= myio_err_readfile 
    endif    

end subroutine read_nucl_orient_triplets

! pre : nelem( s) has been assigned
! post : chain_elem(k,s)%elem(s) allocated

subroutine allocate_chain_elements(nsegAA,nelemAA,chain_elem)

    use chains, only : var_darray

    integer, intent(in)               :: nsegAA
    integer, dimension(:), intent(in) :: nelemAA
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem

    ! local arguments

    integer   :: s, k, j
    
    allocate(chain_elem(3,nsegAA))

    do s=1,nsegAA
        do k=1,3
            allocate(chain_elem(k,s)%elem(nelemAA(s)))
        enddo     
    enddo 

    !call print_nucl_elements(nsegAA,nelemAA,chain_elem)

end subroutine allocate_chain_elements


subroutine allocate_nucl_chain_elements(nnucl,nsegAA,nelemAA,chain_elem_rot)

    use chains, only : var_darray

    integer, intent(in)               :: nnucl
    integer, intent(in)               :: nsegAA
    integer, dimension(:), intent(in) :: nelemAA
    type(var_darray), dimension(:,:,:), allocatable, intent(inout) :: chain_elem_rot

    ! local arguments

    integer   :: s, k, n
    
    allocate(chain_elem_rot(3,nsegAA,nnucl))

    do n=1,nnucl
        do s=1,nsegAA
            do k=1,3
                allocate(chain_elem_rot(k,s,n)%elem(nelemAA(s)))
            enddo     
        enddo
    enddo     


end subroutine allocate_nucl_chain_elements

subroutine print_nucl_elements(nsegAA,nelemAA,chain_elem)

    use chains, only : var_darray

    integer, intent(in)               :: nsegAA
    integer, dimension(:), intent(in) :: nelemAA
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem

    integer   :: s, k, j

    do s=1,nsegAA
        do k=1,3
            do j=1,nelemAA(s)
                print*,"s=",s,"k=",k,"j=",j," ",chain_elem(k,s)%elem(j)
            enddo     
        enddo
    enddo  

end subroutine print_nucl_elements

! Computes segment numbers that correspond to first AA of nth nucleosome.
! inputs : nseg      = number of segments
!          nsegtypes = number of segment types
!          nnucl     = number of nuclesomes
! output : segAAstart = array of integer = nth val= segment number that is first AA of nth nucleosome

subroutine compute_segnumAAstart(nseg,nsegtypes,nnucl,segnumAAstart)

    use chains, only : type_of_monomer, type_of_monomer_char, mapping_num_to_char  

    integer, intent(in) :: nseg 
    integer, intent(in) :: nsegtypes
    integer, intent(in) :: nnucl 
    integer, dimension(:), intent(inout) :: segnumAAstart

    ! .. local arguments

    character(len=3) :: list_type_char_DNA(6)  ! list of DNA elements in char
    integer          :: list_type_int_DNA(6)   ! list of DNA elements in integer 
    integer          :: i, s, t, k, nnucl_counter, type_int
    character(len=3) :: type_char  
    logical          :: isMonomerDNA, isPreviousMonomerDNA

    ! .. init list DNA segments  
    
    list_type_char_DNA = (/"P  ","S  ","A  ","C  ","G  ","T  "/)
    list_type_int_DNA  = 0

    ! .. make list_type_init_DNA"

    do k=1,6                                ! 6 = number DNA elements 
        type_char= list_type_char_DNA(k)
        do t=1,nsegtypes
            if( mapping_num_to_char(t) == type_char)  then 
                list_type_int_DNA(k)=t
            endif
        enddo 
    enddo

    ! locate segment numbers that is correspond to first AA of nth nucl 

    nnucl_counter=1
    isPreviousMonomerDNA=.true.
    
    do s=1,nseg 
        type_int=type_of_monomer(s)
        isMonomerDNA=.false.
        do k=1,6
            if(type_int == list_type_int_DNA(k)) isMonomerDNA=.true.        
        enddo     
        if(isPreviousMonomerDNA .and. .not. isMonomerDNA) then 
            segnumAAstart(nnucl_counter)=s 
            isPreviousMonomerDNA =.false.
            nnucl_counter=nnucl_counter+1
        endif     
        if(isMonomerDNA) isPreviousMonomerDNA =.true.
    enddo 
                
end subroutine compute_segnumAAstart


! Computes segment numbers that correspond to last AA of nth nucleosome.
! inputs : nnucl     = number of nuclesomes
!          nsegAA    = number of AA segments per nucleosome
!          segnumAAstart = segment numbers that correspond to first AA of nth nucleosome.
! output : segnumAAsend = array of integer = nth val= segment number that is last AA of nth nucleosome

subroutine compute_segnumAAend(nnucl,nsegAA,segnumAAstart,segnumAAend)

    integer, intent(in) :: nnucl
    integer, intent(in) :: nsegAA
    integer, dimension(:), intent(in) :: segnumAAstart
    integer, dimension(:), intent(inout) :: segnumAAend 

    integer :: n 

    do n=1,nnucl 
        segnumAAend(n)=segnumAAstart(n)+nsegAA-1
    enddo 

end subroutine compute_segnumAAend 

subroutine make_nelem(nseg,nsegAA,nnucl,segnumAAstart,nelemAA,nelem)

    integer,               intent(in) :: nseg
    integer,               intent(in) :: nsegAA
    integer,               intent(in) :: nnucl
    integer, dimension(:), intent(in) :: segnumAAstart
    integer, dimension(:), intent(in) :: nelemAA
    integer, dimension(:), intent(inout) :: nelem

    integer :: n , sn, s

    nelem=1

    do n=1,nnucl 
        do sn=1,nsegAA
            s=sn+segnumAAstart(n)-1
            nelem(s)=nelemAA(sn)
        enddo 
    enddo 
    
end subroutine make_nelem


! Computes chain_index(k,s)%elem(j) coordinate (component k) of segment s and element j
! by adding chain_elem_rot to chain_rot(k,s)
! Note chain_elem_rot distance relative to cm such be relative to location CA of every AA  
! inputs : nseg      = number of segments
!          nsegAA    = number of AA segments per nucleosome
!          nnucl     = number of nucleosomes
!          segnumAAstart = segment numbers that corresponds to first AA of nth nucleosome
!          segnumAAsend  = segment numbers that corresponds to last AA of nth nucleosome
!          nelemAA   = array list number elements for each segments
!          chain_rot = array of coordinates of chain for AA only the CA elements
!          chain_elem_rot => chain_elem_rot(k,sn,n)%elem(j) orienteted AA elements for n nucleosome
! output : chain_elem_index => chain_index(k,s)%elem(j) = coordinate (component k) of segment s and element j


subroutine add_chain_rot_and_chain_elem_rot(nseg,nsegAA,nnucl,segnumAAstart,segnumAAend,nelemAA,&
    chain_rot,chain_elem_rot,chain_elem_index) 

    use chains, only : var_darray

    integer, intent(in) :: nseg
    integer, intent(in) :: nsegAA
    integer, intent(in) :: nnucl
    integer, dimension(:), intent(in) :: segnumAAstart
    integer, dimension(:), intent(in) :: segnumAAend 
    integer, dimension(:), intent(in) :: nelemAA
    real(dp), dimension(:,:), intent(in) :: chain_rot 
    type(var_darray), dimension(:,:,:), allocatable, intent(in) :: chain_elem_rot
    type(var_darray), dimension(:,:), allocatable, intent(inout) ::chain_elem_index

    integer :: n, s, j, k, sn

    ! chain_elem_rot distance relative to cm such be relative to location CA of every AA  

    ! AA coordinates
    do n=1,nnucl                    ! loop over number of nucleosomes
        do sn=1,nsegAA              ! loop over number of AA 
            s=sn+segnumAAstart(n)-1 ! segment number 
            do j=1,nelemAA(sn)      ! loop over number of chain elements for segment s 
                do k=1,3
                    chain_elem_index(k,s)%elem(j)= chain_rot(k,s)+(chain_elem_rot(k,sn,n)%elem(j)-chain_elem_rot(k,sn,n)%elem(1))
                enddo
            enddo 
        enddo 
    enddo  

    ! DNA coordinates

    do s=1,segnumAAstart(1)-1 
        do k=1,3
            chain_elem_index(k,s)%elem(1)= chain_rot(k,s)
        enddo  
    enddo        
    do n=2,nnucl-1 
        do s=segnumAAend(n)+1,segnumAAstart(n+1)-1 
            do k=1,3
                chain_elem_index(k,s)%elem(1)= chain_rot(k,s)
            enddo 
        enddo 
    enddo 
    do s=segnumAAend(nnucl)+1,nseg 
        do k=1,3
            chain_elem_index(k,s)%elem(1)= chain_rot(k,s)
        enddo  
    enddo           

end subroutine add_chain_rot_and_chain_elem_rot

end module chaingenerator
