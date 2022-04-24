module MConfig

!     The origin is located at the center of the simulation box.

use MDimension
use MUtil

implicit none

! Signature for identifying configuration files:
character(4), parameter :: FourCC = 'DEMC'

! Version of this module (that should be impressed in any output file):
integer(1), parameter, private :: Version = 2_1

! Number of bytes for data in compact configuration format:
integer, parameter, private :: ncc = 2

! Tags for complete or compact format (please do not change them):
integer(1), parameter :: Complete = 72_1, &
                         Compact  = 28_1

! Scaling factors of the system:
real(rb) :: LengthScale, &
            TimeScale,   &
            MassScale

! Dimensions of the simulation box:
real(rb) :: Lbox(dim), InvLbox(dim)

! Derived type concerning component properties:
type TComponent
  integer(ib) :: Number              ! Current number of particles
  real(rb)    :: Radius,          &  ! Radius of the component
                 Density,         &  ! Density of the component
                 Stiffness,       &  ! Stiffness of a component
                 ShearStiffness,  &  ! Shear stiffness of a component
                 Mass,            &  ! Mass of a component
                 Inertia             ! Moment of inertia of the component
end type TComponent

! Module variables used to handle wall shapes:
integer(ib), parameter :: NS = 2_ib ! Number of shapes
integer(ib), parameter :: Rectangle = 1_ib,  &
                          Cone      = 2_ib

! Module variables used to optimize overlap calculations:
integer(ib), private :: NWS(NS)

type, private :: TConeComp
  real(rb) :: MaxAlphaU, MinAlphaWSq, MaxAlphaWSq
end type TConeComp
type(TConeComp), allocatable, private :: ConeComp(:,:)

type, private :: TCone
  real(rb) :: Rap, SinTheta, CosTheta
end type TCone
type(TCone), allocatable, private :: ConeProp(:)

! Global variables:
real(rb) :: Time = 0.0_rb

integer(ib) :: NC = 0_ib, &  ! Number of components
               NP = 0_ib, &  ! Number of particles
               NW = 0_ib     ! Number of walls

type (TComponent), allocatable :: Component(:)
type (TParticle),  allocatable :: Particle(:)
type (TDirection), allocatable :: Direction(:)
type (TWall),      allocatable :: Wall(:)

real(rb) :: FluidDensity

! Generic prodedures:
interface SaveConfiguration
  module procedure SaveConfToFile
  module procedure SaveConfToUnit
end interface SaveConfiguration

interface ReadConfiguration
  module procedure ReadConfFromFile
  module procedure ReadConfFromUnit
end interface ReadConfiguration

contains

!===============================================================================

  subroutine SaveConfToFile( File, Form, Error )

    character(*), intent(in) :: File
    integer(1), intent(in)   :: Form
    integer, intent(out), optional :: Error

    integer :: Unit

    Unit = BinOpen( File, 'write' )
    if (present(Error)) then
      call SaveConfToUnit( Unit, Form, Error )
    else
      call SaveConfToUnit( Unit, Form )
    end if
    call BinClose(Unit)

  end subroutine SaveConfToFile

!===============================================================================

  subroutine SaveConfToUnit( Unit, Form, Error )

    integer,    intent(in) :: Unit
    integer(1), intent(in) :: Form
    integer, intent(out), optional :: Error

    !     This subroutine saves the current configuration into an open binary
    ! file, i.e., a file that was opened using the following command:
    !
    ! Unit = BinOpen( file = <file name>, action = "write" )
    !
    ! Two file formats are supported:
    !
    ! 1) Complete format -> All the kinematic properties of the comfiguration
    !                       are saved as real(rb) variables.
    !
    ! 2) Compact format ->  Only the geometrical properties of the configuration
    !                       are saved, being the radii of the components, the
    !                       positions of the particles, and the direction vectors
    !                       saved as integer(ncc) and the others as real(rb)
    !                       variables.
    !
    ! Sintax:
    !      call SaveConfiguration( Unit, Complete )
    ! or
    !      call SaveConfiguration( Unit, Compact )

    integer      :: IOStat = 0
    integer(ncc) :: R(dim), phi(dim-1)
    integer(ib)  :: i, comp
    real(rb)     :: Factor, A, B, Highest
    logical(1)   :: Full, Vectors

    Full = (Form == Complete)

    Highest = real(huge(1_ncc)-1,rb)
    Factor = 2.0_rb*Highest/maxval(Lbox)
    A = 0.8_rb*Highest/acos(-1.0_rb)
    B = -0.6_rb*Highest

    call BinWrite( Unit, FourCC, IOStat )
    call BinWrite( Unit, Version, IOStat )
    call BinWrite( Unit, int(dim,1), IOStat )
    call BinWrite( Unit, Form, IOStat )

    call BinWrite( Unit, LengthScale, IOStat )
    call BinWrite( Unit, TimeScale, IOStat )
    if (Full) call BinWrite( Unit, MassScale, IOStat )
    call BinWrite( Unit, Time, IOStat )
    call BinWrite( Unit, Lbox, dim, IOStat )
    call BinWrite( Unit, NC, IOStat )

    do i = 1_ib , NC
      call BinWrite( Unit, Component(i)%Number, IOStat )
      call BinWrite( Unit, Component(i)%Radius, IOStat )
      if (Full) then
        call BinWrite( Unit, Component(i)%Density, IOStat )
        call BinWrite( Unit, Component(i)%Stiffness, IOStat )
        call BinWrite( Unit, Component(i)%ShearStiffness, IOStat )
        call BinWrite( Unit, Component(i)%Mass, IOStat )
        call BinWrite( Unit, Component(i)%Inertia, IOStat )
      end if
    end do

    call BinWrite( Unit, NW, IOStat )

    do i = 1_ib , NW

      call BinWrite( Unit, Wall(i)%Shape, IOStat )
      call BinWrite( Unit, Wall(i)%L1, IOStat )
      call BinWrite( Unit, Wall(i)%L2, IOStat )
      call BinWrite( Unit, Wall(i)%L3, IOStat )
      if (Full) then
        call BinWrite( Unit, Wall(i)%U, dim, IOStat )
        call BinWrite( Unit, Wall(i)%W, dim, IOStat )
        call BinWrite( Unit, Wall(i)%N, dim, IOStat )
      else
        phi = int(A*CartToSpheric(Wall(i)%U) + B,ncc)
        call BinWrite( Unit, phi, dim-1, IOStat )
        phi = int(A*CartToSpheric(Wall(i)%W) + B,ncc)
        call BinWrite( Unit, phi, dim-1, IOStat )
      end if
      call BinWrite( Unit, Wall(i)%Position, dim, IOStat )
      if (Full) then
        call BinWrite( Unit, Wall(i)%Velocity, dim, IOStat )
        call WriteAxialVector( Unit, Wall(i)%Omega, IOStat )
        call BinWrite( Unit, Wall(i)%Stiffness, IOStat )
        call BinWrite( Unit, Wall(i)%ShearStiffness, IOStat )
      end if

    end do

    Vectors = allocated(Direction)

    call BinWrite( Unit, Vectors, IOStat )

    if (Full) then

      do i = 1_ib , NP
        call BinWrite( Unit, Particle(i)%Index, IOStat )
        call BinWrite( Unit, Particle(i)%Position, dim, IOStat )
        call BinWrite( Unit, Particle(i)%Velocity, dim, IOStat )
        call WriteAxialVector( Unit, Particle(i)%Omega, IOStat )
        if (Vectors) then
          call BinWrite( Unit, Direction(i)%U, dim, IOStat )
          call BinWrite( Unit, Direction(i)%W, dim, IOStat )
          call BinWrite( Unit, Direction(i)%N, dim, IOStat )
        end if
      end do

    else

      do comp = 1_ib, NC
        do i = 1_ib, NP
          if ( Particle(i)%Index == comp ) then
            R = int(Factor*Particle(i)%Position,ncc)
            call BinWrite( Unit, R, dim, IOStat )
            if (Vectors) then
              phi = int(A*CartToSpheric(Direction(i)%U) + B,ncc)
              call BinWrite( Unit, phi, dim-1, IOStat )
              phi = int(A*CartToSpheric(Direction(i)%W) + B,ncc)
              call BinWrite( Unit, phi, dim-1, IOStat )
            end if
          end if
        end do
      end do

    end if

    if (present(Error)) then
       Error = IOStat
    else if (IOStat /= 0) then
      stop
    end if

  end subroutine SaveConfToUnit

!===============================================================================

  subroutine ReadConfFromFile( File, Error )

    character(*), intent(in) :: File
    integer, intent(out), optional :: Error

    integer :: Unit

    Unit = BinOpen( File, 'read' )
    if (present(Error)) then
      call ReadConfFromUnit( Unit, Error )
    else
      call ReadConfFromUnit( Unit )
    end if
    call BinClose(Unit)

  end subroutine ReadConfFromFile

!===============================================================================

  subroutine ReadConfFromUnit( Unit, Error )

    integer, intent(in) :: Unit
    integer, intent(out), optional :: Error

    !     This subroutine reads a configuration from an open binary file, i.e.,
    ! a file that was opened using the following command:
    !
    ! Unit = BinOpen( file = <file name>, action = "write" )
    !
    ! Two file formats are supported:
    !
    ! 1) Complete format -> All the kinematic properties of the comfiguration
    !                       were saved as real(rb) variables.
    !
    ! 2) Compact format ->  Only the geometrical properties of the configuration
    !                       were saved, being the radii of the components, the
    !                       positions of the particles, and the direction vectors
    !                       as integer(ncc) and the others as real(rb) variables.
    !
    ! Sintax:
    !      call ReadConfiguration( Unit )
    ! or
    !      call ReadConfiguration( Unit )

    integer      :: i, j, IOStat = 0
    integer(1)   :: Form, FileVersion, FileDim
    integer(ncc) :: R(dim), phi(dim-1)
    integer(ib)  :: comp
    real(rb)     :: Factor, A, B, Pi, Highest
    logical(1)   :: Full, Vectors
    character(4) :: Signature

    Highest = real(huge(1_ncc)-1,rb)
    Pi = acos(-1.0_rb)
    A = 1.25_rb*Pi/Highest
    B = 0.75_rb*Pi

    call BinRead( Unit, Signature, IOStat )

    if (Signature /= FourCC) then
      if (present(Error)) then
        Error = 10
        return
      else
        write(*,'("Error: Not a configuration file.")')
        stop
      end if
    end if

    call BinRead( Unit, FileVersion, IOStat )
    call BinRead( Unit, FileDim, IOStat )
    if (FileDim /= dim) then
      if (present(Error)) then
        Error = 10
        return
      else
        write(*,'("Error: Invalid dimension")')
      end if
    end if
    call BinRead( Unit, Form, IOStat )

    Full = Form == Complete

    call BinRead( Unit, LengthScale, IOStat )
    call BinRead( Unit, TimeScale, IOStat )
    if (Full) call BinRead( Unit, MassScale, IOStat )
    call BinRead( Unit, Time, IOStat )
    call BinRead( Unit, Lbox, dim, IOStat )

    Factor = 0.5_rb*maxval(Lbox)/Highest

    call BinRead( Unit, NC, IOStat )

    if (allocated(Component)) deallocate(Component)
    allocate( Component(NC) )
    do i = 1_ib , NC
      call BinRead( Unit, Component(i)%Number, IOStat )
      call BinRead( Unit, Component(i)%Radius, IOStat )
      if (Full) then
        call BinRead( Unit, Component(i)%Density, IOStat )
        call BinRead( Unit, Component(i)%Stiffness, IOStat )
        if (FileVersion == 2_1) then
          call BinRead( Unit, Component(i)%ShearStiffness, IOStat )
        else
          Component(i)%ShearStiffness = 2.0_rb/7.0_rb*Component(i)%Stiffness
        end if
        call BinRead( Unit, Component(i)%Mass, IOStat )
        call BinRead( Unit, Component(i)%Inertia, IOStat )
      end if
    end do

    call BinRead( Unit, NW )
    if (allocated(Wall)) deallocate(Wall)
    allocate( Wall(NW) )

    do i = 1_ib , NW

      call BinRead( Unit, Wall(i)%Shape, IOStat )
      call BinRead( Unit, Wall(i)%L1, IOStat )
      call BinRead( Unit, Wall(i)%L2, IOStat )
      call BinRead( Unit, Wall(i)%L3, IOStat )
      if (Full) then
        call BinRead( Unit, Wall(i)%U, dim, IOStat )
        call BinRead( Unit, Wall(i)%W, dim, IOStat )
        call BinRead( Unit, Wall(i)%N, dim, IOStat )
      else
        call BinRead( Unit, phi, dim-1, IOStat )
        Wall(i)%U = SphericToCart( A*real(phi,rb) + B )
        call BinRead( Unit, phi, dim-1, IOStat )
        Wall(i)%W = SphericToCart( A*real(phi,rb) + B )
        call GramSchmidt( Wall(i)%U, Wall(i)%W )
        Wall(i)%N = Wall(i)%U .x. Wall(i)%W
      end if
      call BinRead( Unit, Wall(i)%Position, dim, IOStat )
      if (Full) then
        call BinRead( Unit, Wall(i)%Velocity, dim, IOStat )
        call ReadAxialVector( Unit, Wall(i)%Omega, IOStat )
        if (FileVersion == 2_1) then
          call BinRead( Unit, Wall(i)%Stiffness, IOStat )
          call BinRead( Unit, Wall(i)%ShearStiffness, IOStat )
        end if
      end if

    end do

    !call ConfigureWalls

    NP = sum(Component%Number)
    if (allocated(Particle)) deallocate(Particle)
    allocate( Particle(NP) )

    call BinRead( Unit, Vectors, IOStat )
    if (allocated(Direction)) deallocate(Direction)
    if (Vectors) then
      allocate(Direction(NP))
    end if

    if (Full) then

      do i = 1_ib , NP
        call BinRead( Unit, Particle(i)%Index, IOStat )
        call BinRead( Unit, Particle(i)%Position, dim, IOStat )
        call BinRead( Unit, Particle(i)%Velocity, dim, IOStat )
        call ReadAxialVector( Unit, Particle(i)%Omega, IOStat )
        if (Vectors) then
          call BinRead( Unit, Direction(i)%U, dim, IOStat )
          call BinRead( Unit, Direction(i)%W, dim, IOStat )
          call BinRead( Unit, Direction(i)%N, dim, IOStat )
        end if
      end do

    else

      i = 0
      do comp = 1, NC
        do j = 1, Component(comp)%Number
          i = i + 1
          call BinRead( Unit, R, dim, IOStat )
          Particle(i)%Index = comp
          Particle(i)%Position = Factor*real(R,rb)
          if (Vectors) then
            call BinRead( Unit, phi, dim-1, IOStat )
            Direction(i)%U = SphericToCart( A*real(phi,rb) + B )
            call BinRead( Unit, phi, dim-1, IOStat )
            Direction(i)%W = SphericToCart( A*real(phi,rb) + B )
            call GramSchmidt( Direction(i)%U, Direction(i)%W )
            Direction(i)%N = Direction(i)%U .x. Direction(i)%W
          end if
        end do
      end do

    end if

    if (present(Error)) then
      Error = IOStat
    else if ( IOStat /= 0 ) then
      write(*,'("Error: BinFiles error number ",I2)') IOStat
      stop
    end if

  end subroutine ReadConfFromUnit

!===============================================================================

  subroutine ReadParticlePositions( ConfigFile )

    character(*), intent(in) :: ConfigFile

    integer(ib) :: i, SavedNC, SavedNW
    real(rb) :: SavedLbox(dim)
    real(rb) :: SavedLengthScale, SavedTimeScale, SavedMassScale
    real(rb) :: SavedTime
    type (TComponent), allocatable :: SavedComponent(:)
    type (TWall), allocatable :: SavedWall(:)

    SavedLbox = Lbox
    SavedTime = Time
    SavedLengthScale = LengthScale
    SavedTimeScale = TimeScale
    SavedMassScale = MassScale
    SavedNC = NC
    allocate( SavedComponent(NC) )
    SavedComponent = Component
    deallocate( Component )
    SavedNW = NW
    if (NW /= 0_ib) then
      allocate( SavedWall(NW) )
      SavedWall = Wall
      deallocate(Wall)
    end if

    call ReadConfiguration( ConfigFile )

    if (LengthScale /= SavedLengthScale) then
      forall (i=1_ib:NP)
        Particle(i)%Position = LengthScale/SavedLengthScale*Particle(i)%Position
      end forall
    end if

    Lbox = SavedLbox
    Time = SavedTime
    LengthScale = SavedLengthScale
    TimeScale = SavedTimeScale
    MassScale = SavedMassScale
    NC = SavedNC
    deallocate(Component)
    allocate(Component(NC))
    Component = SavedComponent
    deallocate(SavedComponent)
    do i = 1_ib, NC
      Component(i)%Number = count(Particle%Index == i)
    end do
    NW = SavedNW
    if (NW /= 0_ib) then
      if (allocated(Wall)) deallocate(Wall)
      allocate(Wall(NW))
      Wall = SavedWall
      deallocate(SavedWall)
    end if

  end subroutine ReadParticlePositions

!===============================================================================

  function RandomConfiguration( Dummy ) result ( Fit )

    integer(4), intent(inout) :: Dummy
    logical :: Fit

    integer, parameter :: MaxAttempts = 1000

    Fit = allocated(Component).and.(NC > 0_ib).and.all(Lbox > 0.0_rb)
    if (.not.Fit) return

    NP = sum(Component%Number)
    if (allocated(Particle)) deallocate(Particle)
    allocate( Particle(NP) )

  end function RandomConfiguration

!===============================================================================

  function SimpleLatticeConfiguration( Dummy, Diameter ) result ( Fit )

    integer(4), intent(inout) :: Dummy
    real(rb), intent(in), optional :: Diameter

    logical :: Fit

    real(rb)    :: MinLCell, Lcell(dim), Ri, Dij, Nij(dim), Qji(dim)
    integer(ib) :: iCell, i, Last, j, m, n, MCell(dim), ic(dim), NCell
    integer(ib), allocatable :: ListC(:), ListP(:)
    logical :: Possible, Cylinder, PBC(dim) = .false.
    type(TParticle) :: P

    Fit = allocated(Component).and.(NC > 0_ib).and.all(Lbox > 0.0_rb)
    if (.not.Fit) return

    NP = sum(Component%Number)
    if (allocated(Particle)) deallocate(Particle)
    allocate( Particle(NP) )

    MinLcell = 2.00001_rb*maxval(Component%Radius)

    MCell = int(Lbox/MinLcell,ib)
    where (MCell < 1)
      MCell = 1
    end where

    LCell = Lbox/real(MCell,rb)

    NCell = product(MCell)

    Fit = Ncell >= Np

    if (Fit) then

      allocate( ListC(NCell), ListP(NP) )

      Cylinder = present(Diameter)

!      ListC = RandomList(NCell,Dummy) - 1_ib
      ListC = (/(m,m=0,NCell-1)/)
      ListP = RandomList(NP,Dummy)

      m = 1_ib
      j = 1_ib
      P%Index = 0_ib
      Last = 0_ib
      do while ( (m <= NCell).and.(j <= Np) )

        if (j > Last) then
          P%Index = P%Index + 1_ib
          Last = Last + Component(P%Index)%Number
          Ri = Component(P%Index)%Radius
        end if

        i = ListP(j)

        iCell = ListC(m)

        ic = LatticeCoordinates( iCell, MCell )
        where (MCell == 1)
          P%Position = 0.0_rb
        elsewhere
          P%Position = (Lbox-Lcell)*(real(ic,rb)/real(Mcell-1_ib,rb)-0.5_rb)
        end where

        if (Cylinder) then
          Possible = P%Position(1)**2 + P%Position(2)**2 <= (0.5_rb*Diameter - Ri)**2
        else
          Possible = .true.
        end if

        n = 0_ib
        do while ( Possible .and. (n < NW) )
          n = n + 1_ib
          Possible = .not.Overlap( P, Wall(n), Ri, PBC, Dij, Nij, Qji )
        end do

        if (Possible) then
          Particle(i) = P
          j = j + 1_ib
        end if

        m = m + 1_ib

      end do

      Fit = j > Np

      if (Fit) then
        forall(i=1:dim)
          Particle%Velocity(i) = 0.0_rb
          Particle%Omega(i) = 0.0_rb
          Particle%Acceleration(i) = 0.0_rb
          Particle%Alpha(i) = 0.0_rb
        end forall
      end if

      deallocate( ListC, ListP )

    end if

  end function SimpleLatticeConfiguration

!===============================================================================

  subroutine ConfigureWalls

    integer(ib) :: j, Ord, Ind, Shape
    real(rb) :: Ri
    real(rb) :: L_2, L, Ra, Rp, Aux ! Variables for conic walls

    ! Define wall order:
    NWS = 0_ib
    do j = 1_ib, NW
      Shape = Wall(j)%Shape
      NWS(Shape) = NWS(Shape) + 1
      Wall(j)%Order = NWS(Shape)
    end do

    ! Allocate variable for conic walls:
    if (.not.allocated(ConeComp)) then
      allocate( ConeComp(NWS(Cone),NC), ConeProp(NWS(Cone)) )
    end if

    ! Define the auxiliary variables for every wall:
    do j = 1_ib, NW

      Ord = Wall(j)%Order

      select case (Wall(j)%Shape)
        case (Rectangle)
        case (Cone)

          L_2 = Wall(j)%L1
          L  = 2.0_rb*L_2
          Ra = Wall(j)%L2
          Rp = Wall(j)%L3

          Aux = sqrt(L**2 + (Rp-Ra)**2)
          ConeProp(Ord)%Rap = Aux
          ConeProp(Ord)%SinTheta = L/Aux
          ConeProp(Ord)%CosTheta = (Rp-Ra)/Aux

          do Ind = 1_ib, NC

            Ri = Component(Ind)%Radius

            ConeComp(Ord,Ind)%MaxAlphaU = L_2 + Ri
            ConeComp(Ord,Ind)%MinAlphaWSq = max(0.0_rb, min(Ra,Rp) - Ri)**2
            ConeComp(Ord,Ind)%MaxAlphaWSq = (max(Ra,Rp) + Ri)**2

          end do

      end select

    end do

  end subroutine ConfigureWalls

!===============================================================================

  logical function Overlap( iPart, jWall, Ri, PBC, Dij, Nij, Qji )

    type(TParticle), intent(in) :: iPart
    type(TWall), intent(in) :: jWall
    real(rb), intent(in) :: Ri
    real(rb), intent(out) :: Dij, Nij(dim), Qji(dim)
    logical, intent(in) :: PBC(dim)

    integer :: k
    real(rb) :: Rij(dim), AlphaU, AlphaW, AlphaN,   &
                absAlphaU, absAlphaW, absAlphaN,    &
                AlphaWSq, DijSq, RijW(dim), Lambda, &
                SignedDij, Lby2_AlphaU, Ra_AlphaW

    type(TConeComp) :: MyConeComp
    type(TCone)     :: MyCone

    Rij = jWall%Position - iPart%Position

    forall (k=1:dim,PBC(k))
      Rij(k) = Rij(k) - Lbox(k)*nint( Rij(k)*InvLbox(k) )
    end forall

    select case (jWall%Shape)

      case (Rectangle)

        AlphaN = Rij .dot. jWall%N
        absAlphaN = abs(AlphaN)
        Overlap = absAlphaN < Ri

        if (Overlap) then
          AlphaU = Rij .dot. jWall%U
          AlphaW = Rij .dot. jWall%W
          absAlphaU = abs(alphaU)
          absAlphaW = abs(alphaW)
          Overlap = (absAlphaU < jWall%L1 + Ri) .and. &
                    (absAlphaW < jWall%L2 + Ri)
          if (Overlap) then
            if ( (absAlphaU < jWall%L1) .and. &
                 (absAlphaW < jWall%L2) ) then
              Nij = sign(1.0_rb,AlphaN)*jWall%N
              Dij = absAlphaN
              Qji = -(AlphaU*jWall%U + AlphaW*jWall%W)
            else
              AlphaU = max(-jWall%L1,min(AlphaU,jWall%L1))
              AlphaW = max(-jWall%L2,min(AlphaW,jWall%L2))
              Qji = -(AlphaU*jWall%U + AlphaW*jWall%W)
              Nij = Rij + Qji
              DijSq = sum(Nij*Nij)
              Overlap = DijSq < Ri*Ri
              if (Overlap) then
                Dij = sqrt(DijSq)
                Nij = Nij/Dij
              end if
            end if
          end if
        end if

      case (Cone)

        AlphaU = Rij .dot. jWall%U
        MyConeComp = ConeComp(jWall%Order,iPart%Index)
        Overlap = abs(AlphaU) < MyConeComp%MaxAlphaU
        if (Overlap) then
          RijW = Rij - AlphaU*jWall%U
          AlphaWSq = sum(RijW*RijW)
          Overlap = (AlphaWSq >= MyConeComp%MinAlphaWSq) .and. &
                    (AlphaWSq <  MyConeComp%MaxAlphaWSq)
          if (Overlap) then
            AlphaW = sqrt(AlphaWSq)
            Lby2_AlphaU = jWall%L1 - AlphaU
            Ra_AlphaW   = jWall%L2 - AlphaW
            MyCone = ConeProp(jWall%Order)
            Lambda = Lby2_AlphaU*MyCone%SinTheta - &
                     Ra_AlphaW*MyCone%CosTheta
            if ( (Lambda >= 0.0_rb).and.(Lambda <= MyCone%Rap) ) then
              SignedDij = -(Lby2_AlphaU*MyCone%CosTheta + &
                            Ra_AlphaW*MyCone%SinTheta)
              Dij = abs(SignedDij)
              Overlap = Dij < Ri
              if (Overlap) then
                Nij = MyCone%SinTheta*Wij(RijW,AlphaW,jWall%U) + &
                      MyCone%CosTheta*jWall%U
                Qji = SignedDij*Nij - Rij
                Nij = (Qji + Rij)/Dij
              end if
            else
              Lambda = max(0.0_rb, min(Lambda, MyCone%Rap))
              Nij = (MyCone%SinTheta*Lambda - Lby2_AlphaU)*jWall%U - &
                    (MyCone%CosTheta*Lambda + Ra_AlphaW)*Wij(RijW,AlphaW,jWall%U)
              DijSq = sum(Nij*Nij)
              Overlap = DijSq < Ri*Ri
              if (Overlap) then
                Dij = sqrt(DijSq)
                Qji = Nij - Rij
                Nij = Nij/Dij
              end if
            end if
          end if
        end if

    end select

    contains

      function Wij( RijW, AlphaW, U ) result ( W )

        real(rb), intent(in) :: RijW(dim), AlphaW, U(dim)
        real(rb) :: W(dim)

        real(rb), parameter :: Tol = 1.e-5_rb

        logical :: Loop = .true.
        real(rb) :: WdotU

        if (AlphaW == 0.0_rb) then
          do while (Loop)
            W = RandomUnitVector( Dummy )
            WdotU = sum(W*U)
            Loop = abs( abs(WdotU) - 1.0_rb ) < Tol
          end do
          W = W - WdotU*U
          W = W/sqrt(sum(W*W))
        else
          W = RijW/AlphaW
        end if

      end function Wij

  end function Overlap

!===============================================================================

end module MConfig