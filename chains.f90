!     .. module file of chains variables
module chains
  
    use globals
    implicit none

    type var_array
        integer, allocatable :: elem(:)
    end type var_array
    
    type(var_array), allocatable                :: indexconf(:,:)           ! indexconf(s,alpha)%elem(j) = layer number of conf alpha and segment number s and element j
                                                                            ! used for distrubuted volume 
    integer, dimension(:,:), allocatable        :: indexchain               ! indexchain(s,alpha) = layer number of conf alpha and segment number s
    logical, dimension(:), allocatable          :: isAmonomer               ! isAmonomer(s) =.true. if s is a "A" monomoer  
    integer, dimension(:), allocatable          :: type_of_monomer          ! type of monomer represented as a number
    character(len=3), dimension(:), allocatable :: type_of_monomer_char     ! type of monomer represented as one-three letters
    character(len=3), dimension(:), allocatable :: mapping_num_to_char      ! mapping of type of monomer as a number to one-three letters
    
    logical, dimension(:,:), allocatable        :: ismonomer_of_type        ! ismomomer_of_type(s,t)= true if segment number "s" is of type "t" otherwise false 
    logical, dimension(:), allocatable          :: ismonomer_chargeable     ! ismonomer_chargeabl(s)=true if segment number type "t" is acid or base  
    character(len=1), dimension(:), allocatable :: type_of_charge           ! either "A" =Acid, "B"=base or "N"=neutral for segment number type "t" 
    real(dp), dimension(:), allocatable         :: energychain              ! energy chain   
    real(dp)                                    :: energychain_min          ! mimimum energy chain
    real(dp), dimension(:),   allocatable       :: logweightchain           !  
    logical                                     :: isHomopolymer
    double precision, dimension(:),allocatable  :: lsegseq                  ! segment length only needed for copolymer


    ! chain stuctural quantities

    integer, dimension(:), allocatable          :: segcm                    ! monomerw or units of chain closes to cm of histone  
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

contains


    subroutine allocate_chains(cuantas,nnucl,nseg,nsegtypes,maxnchains,maxnchainsxy)

        integer, intent(in) :: cuantas,nnucl,nseg,nsegtypes
        integer, intent(in) :: maxnchains,maxnchainsxy

        integer :: maxcuantas
    
        maxcuantas=cuantas+maxnchains*maxnchainsxy     ! .. extra  because of  nchain rotations
        
        allocate(indexchain(nseg,maxcuantas))
        allocate(energychain(maxcuantas))
        allocate(logweightchain(maxcuantas))
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
          
    end subroutine allocate_chains
 

    ! Allocates indexconf(s,alpha)%elem(j) = layer number of conf alpha and segment number s and element j
    ! used for distrubuted volume
    ! When used indexchain not need and can be deallocated
    ! inputs: dimension of indexconf: nseg, cuantas, nelem  

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

   


end module chains
