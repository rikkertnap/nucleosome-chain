!     .. module file of chains variables


module chains
  
    use globals
    implicit none

    type var_darray
        real(dp), allocatable :: elem(:)
    end type var_darray

    type var_iarray
        integer, allocatable :: elem(:)
    end type var_iarray
    
    type var_char3array
        character(len=3), allocatable :: elem(:)
    end type var_char3array

    type(var_iarray), allocatable               :: indexconfpair(:,:)       ! indexconfpair(s,alpha)%elem(j) = layer number of conf alpha and 
                                                                            ! segment number s and neighbor j used for distributed volume 
    integer, dimension(:,:), allocatable        :: nneigh                   ! number of neigbors or pairs of segment s in conf alpha used only phosphates  

    type(var_iarray), allocatable               :: indexconf(:,:)           ! indexconf(s,alpha)%elem(j) = layer number of conf alpha and segment number s and element j
                                                                            ! used for distributed volume 
    integer, dimension(:,:), allocatable        :: indexchain               ! indexchain(s,alpha) = layer number of conf alpha and segment number s
    logical, dimension(:), allocatable          :: isAmonomer               ! isAmonomer(s) =.true. if s is a "A" monomoer  
    integer, dimension(:), allocatable          :: index_phos               ! index_phos() = list of all layer number that contain phophates
    integer, dimension(:), allocatable          :: inverse_index_phos       ! inverse_index_phos() = invser of index_phos
    integer, dimension(:), allocatable          :: type_of_monomer          ! type of monomer represented as a number
    character(len=3), dimension(:), allocatable :: type_of_monomer_char     ! type of monomer represented as one-three letters
    character(len=3), dimension(:), allocatable :: mapping_num_to_char      ! mapping of type of monomer as a number to one-three letters
    
    logical, dimension(:,:), allocatable        :: ismonomer_of_type        ! ismomomer_of_type(s,t)= true if segment number "s" is of type "t" otherwise false 
    logical, dimension(:), allocatable          :: ismonomer_chargeable     ! ismonomer_chargeabl(s)=true if segment number type "t" is acid or base  
    character(len=1), dimension(:), allocatable :: type_of_charge           ! either "A" =Acid, "B"=base or "N"=neutral for segment number type "t" 
    real(dp), dimension(:), allocatable         :: energychain              ! energy chain   
    real(dp)                                    :: energychain_min          ! mimimum energy chain
    real(dp), dimension(:),   allocatable       :: logweightchain    
    real(dp), dimension(:),   allocatable       :: energychainLJ            ! Lennard Jones internal VdW energy 
    logical                                     :: isHomopolymer
    double precision, dimension(:),allocatable  :: lsegseq                  ! segment length only needed for copolymer
    
    ! sgraftpts used to be in volume.f90
    integer                                     :: sgraftpts(3)             ! triplet of unit number of histone that is rotated into fixed orientation
    integer, dimension(:,:), allocatable        :: orientation_triplets     ! triplet of unit number for all nnucl histone
    integer, dimension(:,:), allocatable        :: unitvector_triplets      ! triplet of for norm vectos and com for all nnucl histone
    integer, dimension(:), allocatable          :: nelem                    ! number of elements of every segment
    integer, dimension(:), allocatable          :: nelemAA                  ! number of elements of every AA segment
    integer, dimension(:), allocatable          :: typeAA                   ! type of number of elements of every AA segment
    integer, dimension(:), allocatable          :: elem_charge              ! element number that is chargeable 
    character(len=3), dimension(:,:), allocatable :: nucl_elem_type         ! type of nuclesome element in character 

    ! chain stuctural quantities

    integer, dimension(:), allocatable          :: segcm                    ! monomer or units of chain closest to cm of histone  
    integer, dimension(:), allocatable          :: segunitvector            ! monomer or units of chain such that with segcm make an unitvecotr  
    real(dp), dimension(:), allocatable         :: Rgsqr                    ! radius of gyration 
    real(dp), dimension(:), allocatable         :: Rendsqr                  ! end-to-end distance
    real(dp), dimension(:,:), allocatable       :: bond_angle               ! bond angle
    real(dp), dimension(:,:), allocatable       :: dihedral_angle           ! dihedralangle
    real(dp), dimension(:,:), allocatable       :: nucl_spacing             ! spacing or distance between Nuclesome  
    real(dp), dimension(:), allocatable         :: avbond_angle             ! average bond angle
    real(dp), dimension(:), allocatable         :: avdihedral_angle         ! average dihedral angle
    real(dp), dimension(:), allocatable         :: avnucl_spacing           ! average spacing or distance between Nuclesome
    real(dp)                                    :: avRgsqr                  ! radius of gyration 
    real(dp)                                    :: avRendsqr                ! end-to-end distance


    ! .. pairing parameters 

    real(dp) :: distphoscutoff ! distance allow between two phosphate to be a pair
    integer  :: maxneigh       ! maximum of neigbors 
    
    ! 
    integer  :: len_index_phos ! length of array index_phos
 
contains 

    subroutine allocate_chains(cuantas,nnucl,nseg,nsegAA,nsegtypes,nsegtypesAA,maxnchains,maxnchainsxy)

        integer, intent(in) :: cuantas,nnucl,nseg,nsegAA,nsegtypes,nsegtypesAA
        integer, intent(in) :: maxnchains,maxnchainsxy

        integer :: maxcuantas
    
        maxcuantas=cuantas+maxnchains*maxnchainsxy     ! .. extra  because of  nchain rotations
       
        allocate(indexchain(nseg,maxcuantas))
        allocate(energychain(maxcuantas))
        allocate(logweightchain(maxcuantas))
        allocate(energychainLJ(maxcuantas))
        allocate(isAmonomer(nseg)) 
        allocate(type_of_monomer(nseg)) 
        allocate(type_of_monomer_char(nseg))
        allocate(mapping_num_to_char(nsegtypes)) 
        allocate(ismonomer_of_type(nseg,nsegtypes)) 
        allocate(ismonomer_chargeable(nsegtypes))
        allocate(type_of_charge(nsegtypes))

        ! chain stuctural quantities
        allocate(segcm(nnucl))
        allocate(Rgsqr(maxcuantas))
        allocate(Rendsqr(maxcuantas))
        allocate(bond_angle(nnucl-2,maxcuantas))
        allocate(dihedral_angle(nnucl-3,maxcuantas))
        allocate(nucl_spacing(nnucl-1,maxcuantas)) 
        allocate(avnucl_spacing(nnucl-1)) 
        allocate(avbond_angle(nnucl-2))
        allocate(avdihedral_angle(nnucl-3))  

        ! rotational and orientational segments 
        allocate(orientation_triplets(nnucl,3))
        allocate(nelem(nseg))
        allocate(nelemAA(nsegAA))
        allocate(typeAA(nsegtypesAA))
        allocate(elem_charge(nsegtypes))

          
    end subroutine allocate_chains
 

    ! Allocates indexconf(s,alpha)%elem(j) = layer number of conf alpha and segment number s and element j
    ! used for distrubuted volume
    ! When used indexchain is not needed and can be deallocated
    ! inputs: dimension of indexconf: cuantas, nseg and  nelem(:)   

    subroutine allocate_indexconf(cuantas,nseg,nelem)

        integer, intent(in) :: cuantas,nseg
        integer, intent(in) :: nelem(:)
    
        integer :: c, s

        allocate(indexconf(nseg,cuantas))
        do c=1,cuantas
            do s=1,nseg
                allocate(indexconf(s,c)%elem(nelem(s)))
            enddo
        enddo        
        
    end subroutine allocate_indexconf

   
    ! Allocates indexconfpair(s,alpha)%elem(j) = layer number of conf alpha and segment number s and neigbor pair j
    ! used for distributed volume
    ! When used indexchain is not needed and can be deallocated
    ! inputs: dimension of indexconfpair: cuantas, nseg and  nelem(:) 

    subroutine allocate_indexconfpair(cuantas,nseg)

        integer, intent(in) :: cuantas,nseg

        allocate(indexconfpair(nseg,cuantas))  
        
    end subroutine allocate_indexconfpair

   
    ! Allocates neigh : neigbors that segment number s in conf alpha has 
    ! used only for segment s that is a phosphate

    subroutine allocate_nneighbor(cuantas,nseg)

        integer, intent(in) :: cuantas,nseg

        allocate(nneigh(nseg,cuantas))
                  
    end subroutine allocate_nneighbor
   

end module chains
