module MDEM

use MForces
use MUtil

implicit none

! Variables denoting the total simulation time and the integration step:
real(rb) :: TotalTime, TimeStep
integer(ib) :: NumberOfEulerIntervals

! Auxiliary module variables used for the Adams-Bashforth-Moulton
! predictor-corrector integration method:
type, private :: TPrevious
  real(rb) :: Velocity(dim),     &
              Acceleration(dim), &
              Alpha(dim)
end type TPrevious
integer :: MethodOrder, NPrev
type(TParticle), allocatable, private :: Predicted(:)
type(TPrevious), allocatable, private :: Previous(:,:)
integer(ib), allocatable, private :: ipc(:)
real(rb), private :: AB0, AM0
real(rb), allocatable, private :: AB(:), AM(:)

! Auxiliary module variables used for enhancing efficiency:
real(rb), private ::TimeStepBy2, TimeStepSqBy2,  &
                    pHalfLbox(dim), mHalfLbox(dim)

! Module variables used for controling the particle number:
integer(ib) :: MaxNP

! Module variables for controling the particle number:
integer(1), parameter :: Rect   = 1_1, &
                         Circle = 2_1

! Module variables for controling the movement of the walls:
logical :: WallsInMovement = .false.
type TWallMovement
  real(rb) :: ConstantFactor(dim),    &
              Amplitude(dim),         &
              VelocityAmplitude(dim), &
              TwoPiFrequency(dim),    &
              PhaseAngle(dim)
end type
type (TWallMovement), allocatable :: WM(:)

contains

!===============================================================================

  subroutine SaveIntegrationVariables( Unit )

    integer, intent(in) :: Unit

    integer(ib) :: i, j

    ! Save the variables of the integration method:
    call BinWrite( Unit, ipc, NPrev )
    do j = 1_ib, NPrev
      do i = 1_ib, NP
        call BinWrite( Unit, Previous(j,i)%Velocity, dim )
        call BinWrite( Unit, Previous(j,i)%Acceleration, dim )
        call BinWrite( Unit, Previous(j,i)%Alpha, dim )
      end do
    end do

  end subroutine SaveIntegrationVariables

!===============================================================================

  subroutine ReadIntegrationVariables( Unit )

    integer, intent(in) :: Unit

    integer(ib) :: i, j

    allocate( ipc(NPrev), Previous(NPrev,NP), Predicted(NP) )
    call BinRead( Unit, ipc, NPrev )
    do j = 1_ib, NPrev
      do i = 1_ib, NP
        call BinRead( Unit, Previous(j,i)%Velocity, dim )
        call BinRead( Unit, Previous(j,i)%Acceleration, dim )
        call BinRead( Unit, Previous(j,i)%Alpha, dim )
      end do
    end do
    Predicted%Index = Particle(1:NP)%Index

  end subroutine ReadIntegrationVariables

!===============================================================================

  subroutine InitializeIntegrationVariables

    integer :: i

    allocate( ipc(NPrev), Previous(NPrev,NP), Predicted(NP) )
    Predicted%Index = Particle(1:NP)%Index
    forall (i=1:dim)
      Previous%Velocity(i) = 0.0_rb
      Previous%Acceleration(i) = 0.0_rb
      Previous%Alpha(i) = 0.0_rb
    end forall
    forall (i=1:NPrev) ipc(i) = i

  end subroutine InitializeIntegrationVariables

!===============================================================================

  subroutine InitializeAuxiliaries

    real(rb) :: Aux

    ! Save some frequently used expressions:
    TimeStepBy2 = 0.5_rb*TimeStep
    TimeStepSqBy2 = 0.5_rb*TimeStep**2

    pHalfLbox = 0.5_rb*Lbox
    mHalfLbox = -pHalfLbox

    ! Initialize some variables that concern to the integration method:
    NPrev = MethodOrder - 1_ib
    allocate( AB(NPrev), AM(NPrev) )

    ! Define the coefficients of the Adams-Bashforth and Adams-Moulton
    ! formulae:
    select case (MethodOrder)

    case (4)

      Aux = TimeStep/24.0_rb
      AB0   =  55.0_rb*Aux ; AM0   =  9.0_rb*Aux
      AB(1) = -59.0_rb*Aux ; AM(1) = 19.0_rb*Aux
      AB(2) =  37.0_rb*Aux ; AM(2) = -5.0_rb*Aux
      AB(3) =  -9.0_rb*Aux ; AM(3) =  1.0_rb*Aux

    case (3)

      Aux = TimeStep/12.0_rb

      AB0   =  23.0_rb*Aux ; AM0   =  5.0_rb*Aux
      AB(1) = -16.0_rb*Aux ; AM(1) =  8.0_rb*Aux
      AB(2) =   5.0_rb*Aux ; AM(2) = -1.0_rb*Aux

    case (2)

      AB0   =  1.5_rb*TimeStep ; AM0   = 0.5_rb*TimeStep
      AB(1) = -0.5_rb*TimeStep ; AM(1) = 0.5_rb*TimeStep

    case (1)

      AB0 = TimeStep ; AM0 = TimeStep

    end select

    ! Initialize cell list:
    call InitializeCellList
 
    ! Initialize auxiliaries for the calculation of forces:
    call InitializeForce

  end subroutine InitializeAuxiliaries

!===============================================================================

  subroutine FinalizeDEM

    deallocate( Kn, Knw, Yn, Ynw )

    deallocate( MinDij, MinDijSq )

    deallocate( CellIndex )
    call DestroyCellStructure

    deallocate( Previous, Predicted )

    call DestroyContactList

  end subroutine FinalizeDEM

!===============================================================================

  function PotentialEnergy() result ( Ep )

    real(rb) :: Ep

    integer :: i

    Ep = 0.0_rb
    do i = 1, NP
      Ep = Ep - Component(Particle(i)%Index)%Mass *  &
                (Particle(i)%Position .dot. Gravity)
    end do

  end function PotentialEnergy 

!===============================================================================

  function KineticEnergy() result ( Ek )

    real(rb) :: Ek

    integer :: i, Indx

    Ek = 0.0_rb
    do i = 1 , NP
      Indx = Particle(i)%Index
      Ek = Ek + Component(Indx)%Mass*sum(Particle(i)%Velocity**2) + &
                Component(Indx)%Inertia*sum(Particle(i)%Omega**2)
    end do
    Ek = 0.5_rb*Ek

  end function KineticEnergy

!===============================================================================

  function Energy() result ( E )

    real(rb) :: E

    E = PotentialEnergy() + KineticEnergy()

  end function Energy

!===============================================================================

  subroutine Euler( TimeStep )

    real(rb), intent(in) :: TimeStep

    logical, parameter :: Update = .true.

    integer(ib) :: i
    real(rb)    :: TimeStepSqBy2

    TimeStepSqBy2 = 0.5_rb*TimeStep**2

    call ComputeAccelerations( Particle, NP, TimeStep, Update )

    do i = 1_ib, NP

      Particle(i)%Position = Particle(i)%Position +  &
                             Particle(i)%Velocity*TimeStep +  &
                             Particle(i)%Acceleration*TimeStepSqBy2

      Particle(i)%Velocity = Particle(i)%Velocity +  &
                             Particle(i)%Acceleration*TimeStep

      Particle(i)%Omega = Particle(i)%Omega +  &
                          Particle(i)%Alpha*TimeStep

    end do

  end subroutine Euler

!===============================================================================

  subroutine MoveParticlesEuler

    integer(ib) :: i

    real(rb) :: NewTimeStep

    NewTimeStep = TimeStep/real(NumberOfEulerIntervals,rb)

    ipc = cshift( ipc, shift = -1 )
    forall (i=1_ib:NP)
      Previous(ipc(1),i)%Velocity = Particle(i)%Velocity
      Previous(ipc(1),i)%Acceleration = Particle(i)%Acceleration
      Previous(ipc(1),i)%Alpha = Particle(i)%Alpha
    end forall

    do i = 1_ib, NumberOfEulerIntervals
      call Euler( NewTimeStep )
      Time = Time + NewTimeStep
      if (WallsInMovement) call MoveWalls( Time )
    end do

    call ApplyPeriodicBoundaryConditions( Particle, NP )

  end subroutine MoveParticlesEuler

!===============================================================================

  subroutine MoveParticles( step )

    integer(ib), intent(in) :: step

    logical, parameter :: Update = .true.,  &
                          DoNotUpdate = .false.
    integer(ib) :: i, k
    real(rb) :: ABk, AMk
    type(TPrevious) :: Sum, Prev

    ! Prediction step (4th order Adams-Bashforth):
    call ComputeAccelerations( Particle, NP, TimeStep, Update )
    do i = 1_ib, NP

      Sum%Velocity = 0.0_rb
      Sum%Acceleration = 0.0_rb
      Sum%Alpha = 0.0_rb
      do k = 1_ib, NPrev
        ABk = AB(k)
        Prev = Previous(ipc(k),i)
        Sum%Velocity = Sum%Velocity + ABk*Prev%Velocity
        Sum%Acceleration = Sum%Acceleration + ABk*Prev%Acceleration
        Sum%Alpha = Sum%Alpha + ABk*Prev%Alpha
      end do

      Predicted(i)%Position = Particle(i)%Position     + &
                              AB0*Particle(i)%Velocity + &
                              Sum%Velocity

      Predicted(i)%Velocity = Particle(i)%Velocity         + &
                              AB0*Particle(i)%Acceleration + &
                              Sum%Acceleration

      Predicted(i)%Omega = Particle(i)%Omega     + &
                           AB0*Particle(i)%Alpha + &
                           Sum%Alpha

    end do

    call ApplyPeriodicBoundaryConditions( Predicted, NP )

    ipc = cshift( ipc, shift = -1 )
    forall (i=1_ib:NP,NPrev>0)
      Previous(ipc(1),i)%Velocity = Particle(i)%Velocity
      Previous(ipc(1),i)%Acceleration = Particle(i)%Acceleration
      Previous(ipc(1),i)%Alpha = Particle(i)%Alpha
    end forall

    ! Advance in time:
    Time = real(step,rb)*TimeStep

    ! Calculate movement of walls:
    if (WallsInMovement) call MoveWalls( Time )

    ! Correction step (4th order Adams-Moulton):
    call ComputeAccelerations( Predicted, NP, TimeStep, DoNotUpdate )
    do i = 1_ib, NP

      Sum%Velocity = 0.0_rb
      Sum%Acceleration = 0.0_rb
      Sum%Alpha = 0.0_rb
      do k = 1_ib, NPrev
        AMk = AM(k)
        Prev = Previous(ipc(k),i)
        Sum%Velocity = Sum%Velocity + AMk*Prev%Velocity
        Sum%Acceleration = Sum%Acceleration + AMk*Prev%Acceleration
        Sum%Alpha = Sum%Alpha + AMk*Prev%Alpha
      end do

      Particle(i)%Position = Particle(i)%Position      + &
                             AM0*Predicted(i)%Velocity + &
                             Sum%Velocity

      Particle(i)%Velocity = Particle(i)%Velocity          + &
                             AM0*Predicted(i)%Acceleration + &
                             Sum%Acceleration

      Particle(i)%Omega = Particle(i)%Omega      + &
                          AM0*Predicted(i)%Alpha + &
                          Sum%Alpha
    end do

    call ApplyPeriodicBoundaryConditions( Particle, NP )

  end subroutine MoveParticles

!===============================================================================

  subroutine ApplyPeriodicBoundaryConditions( Particle, NP )

    integer(ib), intent(in) :: NP
    type(TParticle), intent(inout) :: Particle(NP)

    integer(ib) :: i
    real(rb) :: Pos(NP)

    do i = 1_ib, dim
      if ( PBC(i) ) then
        Pos = Particle%Position(i)
        Particle%Position(i) = Pos - Lbox(i)*nint( Pos*InvLbox(i) )
      end if
    end do

  end subroutine ApplyPeriodicBoundaryConditions

!===============================================================================

  subroutine MoveWalls( Time )

    real(rb), intent(in) :: Time

    integer(ib) :: j
    real(rb) :: Phase(dim)
    type(TWallMovement) :: jWM

    do j = 1_ib, NW
      jWM = WM(j)
      Phase = jWM%TwoPiFrequency*Time + jWM%PhaseAngle
      Phase = Phase - TwoPi*real(int(Phase/TwoPi,ib),rb)
      Wall(j)%Position = jWM%ConstantFactor + jWM%Amplitude*cos(Phase)
      Wall(j)%Velocity = jWM%VelocityAmplitude*sin(Phase)
    end do

  end subroutine MoveWalls

!===============================================================================

  subroutine PourParticle( TargetNumber, PourShape, PourSize, MaxAttempts )

    integer(ib), intent(in) :: TargetNumber(:), MaxAttempts
    integer(1), intent(in)  :: PourShape
    real(rb), intent(in)    :: PourSize(:)

    integer     :: attempt, i
    integer(ib) :: NNR, iCell, NPSC, j, Rij(dim),    &
                   NotReached(size(TargetNumber)),   &
                   Aux(size(TargetNumber))
    real(rb)    :: Rad, Angle, Di, Ri, Dij, Nij(dim), Qji(dim)
    logical     :: Possible
    type(TParticle) :: iPart, jPart

!<test>
integer(ib) :: iBig(1)
!</test>

    NotReached = 0_ib
    forall (i=1:NC) Aux(i) = i
    NotReached = pack(Aux, &
                      mask = Component%Number < TargetNumber,  &
                      vector = NotReached )
    NNR = count( NotReached /= 0_ib )

    if ( NNR > 0_ib ) then

      ! Choose a component to insert:
      j = int(real(NNR,rb)*ranf(Dummy),ib) + 1_ib
      iPart%Index = NotReached(j)

      ! Save the component radius and diameter:
!      Ri = Component(iPart%Index)%Radius
!<test>
      Ri = maxval(Component%Radius)
!</test>

      Di = 2.0_rb*Ri

      ! Define the height of the particle insertion:
      iPart%Position(dim) = pHalfLbox(dim) - Ri

      ! Try to insert a particle:
      Possible = .false.
      attempt = 0
      do while ( (.not.Possible).and.(attempt < MaxAttempts) )
        Possible = .true.
        attempt = attempt + 1

        ! Choose the position for inserting the particle:
        select case (PourShape)
          case (Rect)
            do i = 1, dim-1
              iPart%Position(i) = (PourSize(i) - Di)*(ranf(Dummy) - 0.5_rb)
            end do
          case (Circle)
            Rad = 0.5_rb*ranf(Dummy)*(PourSize(1) - Di)
            Angle = TwoPi*ranf(Dummy)
            iPart%Position(1) = Rad*sin(Angle)
            iPart%Position(2) = Rad*cos(Angle)
        end select

        ! Verify the cell in which the particle are:
        iCell = CurrentCell3D( iPart%Position(1), &
                               iPart%Position(2), &
                               iPart%Position(3)  )

        ! Create a Monte Carlo-like list of neighbour particles:
        call HinderNeighbourCellsList( iCell )
        call NeighbourCells( iCell, MonteCarlo, PBC(1), PBC(2), PBC(3) )
        call AllNeighbourParticles( iCell, NPSC, Neighbour )
        call RecoverNeighbourCellsList( iCell )


!<test>
iBig = maxloc(Component%Radius)
!</test>

        ! Search for particle-particle overlap:
        j = 0_ib
        do while ( Possible.and.(j < NPSC) )
          j = j + 1_ib
          jPart = Particle(Neighbour(j))
          Rij   = jPart%Position - iPart%Position
!<test>
!          Possible = sum(Rij*Rij) > MinDijSq(iPart%Index,jPart%Index)
          Possible = sum(Rij*Rij) > MinDijSq(iBig(1),jPart%Index)
!</test>
        end do

        ! Search for particle-wall overlap:
        j = 0_ib
        do while ( Possible.and.(j < NW) )
          j = j + 1_ib
!<test>
          Possible = .not.Overlap( iPart, Wall(j), Component(iBig(1))%Radius, PBC, Dij, Nij, Qji )
!          Possible = .not.Overlap( iPart, Wall(j), Ri, PBC, Dij, Nij, Qji )
!</test>
        end do

      end do

      ! Insert the particle, if possible:
      if (Possible) then

        ! Verify if it is necessary to expand some variables in at least 5%:
        if (NP == MaxNP) then
          call ResizeVariables( max( MaxNP + 1_ib, int(1.05_rb*real(MaxNP,rb),ib) ) )
        end if
  
        ! Carry out insersion:
        iPart%Velocity = (/0.0_rb,0.0_rb,-1.e-4_rb*LengthScale/TimeScale/)
        iPart%Omega = 0.0_rb
        call InsertParticle( iPart )

      end if

    end if

  end subroutine PourParticle

!===============================================================================

  subroutine InsertParticle( Part )

    type(TParticle), intent(in) :: Part

    integer :: j

    ! Update the number of particles:
    NP = NP + 1_ib

    ! Add one particle to component "Index":
    Component(Part%Index)%Number = Component(Part%Index)%Number + 1_ib

    ! Update the number of particles:
    Particle(NP) = Part
    Predicted(NP)%Index = Part%Index
    forall(j=1:NPrev)
      Previous(j,NP)%Velocity = 0.0_rb
      Previous(j,NP)%Acceleration = 0.0_rb
      Previous(j,NP)%Alpha = 0.0_rb
    end forall

    ! Calculate the cell index of particle NP and add it to the lattice:
    CellIndex(NP) = CurrentCell3D( Part%Position(1), &
                                   Part%Position(2), &
                                   Part%Position(3)  )
    call AddParticle( NP, CellIndex(NP) )

  end subroutine InsertParticle

!===============================================================================

  subroutine ResizeVariables( NewMaxNP )

    integer(ib), intent(in) :: NewMaxNP

    type(TParticle), allocatable :: Part(:)
    type(TPrevious), allocatable :: Prev(:,:)

    ! Redefine the maximum number of particles:
    MaxNP = NewMaxNP

    ! Allocate auxiliary variables:
    allocate( Part(NP), Prev(NPrev,NP) )

    ! Resize variable "Particle":
    Part = Particle
    deallocate( Particle )
    allocate( Particle(MaxNP) )
    Particle(1:NP) = Part

    ! Resize variable "Predicted":
    deallocate( Predicted )
    allocate( Predicted(MaxNP) )
    Predicted(1:NP)%Index = Particle(1:NP)%Index

    ! Resize variable "Previous":
    Prev = Previous
    deallocate( Previous )
    allocate( Previous(NPrev,MaxNP) )
    Previous(:,1:NP) = Prev

    ! Resize variable "CellIndex":
    Part%Index = CellIndex
    deallocate( CellIndex )
    allocate( CellIndex(MaxNP) )
    CellIndex(1:NP) = Part%Index

    ! Deallocate auxiliary variable:
    deallocate( Part, Prev )

    ! Resize the contact list:
    call ResizeContactList( NP, MaxNP )

  end subroutine ResizeVariables

!===============================================================================

  subroutine UpdateParticleList

    integer(ib) :: i, iCell, NPSC
    type(TParticle) :: iPart

    ! Examine all the particles:
    do i = NP, 1_ib, -1_ib

      ! Save the properties of particle "i":
      iPart = Particle(i)

      ! Verify if particle "i" has ran out from the box:
      if ( any( (iPart%Position > pHalfLbox).or.  &
                (iPart%Position < mHalfLbox) ) ) then

        ! Update the number of particles of the corresponding components:
        Component(iPart%Index)%Number = Component(iPart%Index)%Number - 1_ib

        ! Create a Monte Carlo-like list of neighbour particles:
        iCell = CellIndex(i)
        call HinderNeighbourCellsList( iCell )
        call NeighbourCells( iCell, MonteCarlo, PBC(1), PBC(2), PBC(3) )
        call AllNeighbourParticles( iCell, NPSC, Neighbour )
        call RecoverNeighbourCellsList( iCell )

        ! Remove particle "i" from the cell structure:
        call RemoveParticleFromCell( i, iCell )

        ! Remove particle "i" from the contact list:
        call RemoveParticleFromContactList( i , NPSC, Neighbour )

        ! If necessary, put particle "NP" in the place of particle "i":
        if (i /= NP) then

          Particle(i) = Particle(NP)
          Previous(:,i) = Previous(:,NP)

          ! Do the replacement in the cell structure:
          call ReplaceParticleIndexInCell( NP, CellIndex(NP), i )
          CellIndex(i) = CellIndex(NP)

          ! Do the replacement in the contact list:
          call ChangeIndexInContactList( NP, i , NPSC, Neighbour )

        end if

        ! Update the number of components:
        NP = NP - 1_ib

        ! Contract some variables if MaxNP is at least 10% greater than NP:
        if ( MaxNP >= int(1.1_rb*real(NP,rb),ib) ) then
          call ResizeVariables( int(1.05_rb*real(NP,rb),ib) )
        end if

      end if

    end do


  end subroutine UpdateParticleList

!===============================================================================

end module MDEM