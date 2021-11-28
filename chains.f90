!     .. module file of chains variables
module chains
  
    use globals
    implicit none
  
    integer, dimension(:,:), allocatable :: indexchain               ! index(alpha,s)= layer number of conf alpha and segment number s
    integer, dimension(:,:), allocatable :: indexchain_init 
    logical, dimension(:), allocatable      :: isAmonomer               ! isAmonomer(s) =.true. if s is a "A" monomoer  
    integer, dimension(:), allocatable      :: type_of_monomer          ! type of monomer represented as a number
    character(len=3), dimension(:), allocatable :: type_of_monomer_char ! type of monomer represented as two letters
    logical, dimension(:,:), allocatable    :: ismonomer_of_type        ! ismomomer_of_type(s,t)= true if segment number "s" is of type "t" otherwise false 
    logical, dimension(:), allocatable      :: ismonomer_chargeable     ! ismonomer_chargeabl(s)=true if segment number type "t" is acid or base  
    real(dp), dimension(:), allocatable     :: energychain              ! energy chain   
    real(dp), dimension(:), allocatable     :: energychain_init         ! energy chain   
    real(dp) :: energychain_min              ! mimimum energy chain
    real(dp), dimension(:),   allocatable    :: logweightchain
    logical :: isHomopolymer
    double precision, dimension(:),allocatable :: lsegseq  ! only needed for copolymer


contains


    subroutine allocate_chains(cuantas,nseg,nsegtypes,maxnchains,maxnchainsxy)

        integer, intent(in) :: cuantas,nseg,nsegtypes
        integer, intent(in) :: maxnchains,maxnchainsxy

        integer :: maxcuantas
    
        maxcuantas=cuantas+maxnchains*maxnchainsxy     ! .. extra  because of  nchain rotations
        allocate(indexchain(nseg,maxcuantas))
        allocate(indexchain_init(nseg,maxcuantas))
        allocate(energychain(maxcuantas))
        allocate(energychain_init(maxcuantas))
        allocate(logweightchain(maxcuantas))
        allocate(isAmonomer(nseg)) 
        allocate(type_of_monomer(nseg)) 
        allocate(type_of_monomer_char(nseg))
        allocate(ismonomer_of_type(nseg,nsegtypes)) 
        allocate(ismonomer_chargeable(nsegtypes))

    end subroutine allocate_chains
  
end module chains
