module MForces

use MTypes
use MDimension
use MCells
use MContacts
use MConfig
use MFluid


real(rb) :: Gravity(dim)

real(rb), allocatable :: En(:,:), Enw(:,:), &
                         Et(:,:), Etw(:,:), &
                         Kn(:,:), Knw(:,:), &
                         Kt(:,:), Ktw(:,:), &
                         Yn(:,:), Ynw(:,:), &
                         Yt(:,:), Ytw(:,:), &
                         Mu(:,:), Muw(:,:)

real(rb), allocatable :: MinDij(:,:), MinDijSq(:,:),  &
                         InvMass(:), InvInertia(:),   &
                         EffectiveMass(:),            &
                         MuSq(:,:), MuwSq(:,:)

! Module variables used for controling neighbours lists:
integer, parameter, private :: MaxNeighbours = 50
integer(ib), parameter, private :: MaxNCell = 30000_ib

logical :: PBC(dim)
integer(ib), allocatable :: CellIndex(:)
integer(ib) :: MCellx, MCelly, Mcellz, NCell, &
               Neighbour(MaxNeighbours)

! Module variables for the calculation of forces:
real(rb), private :: Ri, Rj, Rij(dim), DijSq, Dij, Nij(dim), Qji(dim),   &
                     Vij(dim), VijC(dim), iForce(dim),                   &
                     ijForceN(dim),  ijForceT(dim),                      &
                     iTorque(dim), Nij_x_ijForceT(dim), ijDispN(dim),    &
                     ijDispT(dim), ijDispT0(dim), VijT(dim),             &
                     ijForceT0(dim), ijDeltaDT(dim), ijDeltaVT(dim),     &
                     VijT0(dim), DeltaDijSq, HalfTimeStep

logical, private :: Contacting

!<test>
integer(ib) :: Ncol = 0_ib
!</test>

contains

!===============================================================================

  subroutine InitializeForce

    integer(ib) :: i, j

    ! Save the distances for which particles of two components are in contact
    ! and the squares of such distances as well:
    allocate( MinDij(NC,NC), MinDijSq(NC,NC) )
    do i = 1_ib, NC
      MinDij(i,i) = 2.0_rb*Component(i)%Radius
      MinDijSq(i,i) = MinDij(i,i)**2
      do j = i + 1_ib, NC
        MinDij(i,j) = Component(i)%Radius + Component(j)%Radius
        MinDij(j,i) = MinDij(i,j)
        MinDijSq(i,j) = MinDij(i,j)**2
        MinDijSq(j,i) = MinDijSq(i,j)
      end do
    end do

    ! Save the inverse of the mass of each component:
    allocate( InvMass(NC), InvInertia(NC), EffectiveMass(NC) )
    InvMass = 1.0_rb/Component%Mass
    InvInertia = 1.0_rb/Component%Inertia
    EffectiveMass = Component%Mass* &
                   (Component%Density - FluidDensity)/Component%Density

    ! Save thee square of each coefficient of friction:
    allocate( MuSq(NC,NC), MuwSq(NC,NW) )
    MuSq  = Mu**2
    MuwSq = Muw**2

    ! Initialize variables for fluid calculations:
    if (FluidFlow) call InitializeFluid
 
  end subroutine InitializeForce

!===============================================================================

  subroutine InitializeCellList

    real(rb) :: CutOff, Lx, Ly, Lz

    ! Prepare for using the cell structure to maintain neighbour lists:
    CutOff = 2.0_rb*maxval(Component%Radius)
    Lx = Lbox(1)
    Ly = Lbox(2)
    Lz = Lbox(3)

    call InitCellStructure3D( MCellx, MCelly, Mcellz, NCell, &
                              Lx, Ly, Lz, CutOff, MaxNCell   )

    call FullNeighbourCells3D( MolecularDynamics, PBC(1), PBC(2), PBC(3) )

  end subroutine InitializeCellList

!===============================================================================

  subroutine InitializeCellDistribution

    allocate( CellIndex(NP) )
    call DistributeParticles3D( CellIndex, Particle(1:NP)%Position(1), &
                                           Particle(1:NP)%Position(2), &
                                           Particle(1:NP)%Position(3)  )

  end subroutine InitializeCellDistribution

!===============================================================================

  subroutine ReadCellDistribution( Unit )

    integer, intent(in) :: Unit

    allocate( CellIndex(NP) )
    CellIndex(1:NP) = CurrentCell3D( Particle(1:NP)%Position(1), &
                                     Particle(1:NP)%Position(2), &
                                     Particle(1:NP)%Position(3)  )
    call ReadCellStructure( Unit )

  end subroutine ReadCellDistribution

!===============================================================================

  subroutine SaveCellDistribution( Unit )

    integer, intent(in) :: Unit

    call SaveCellStructure( Unit )

  end subroutine SaveCellDistribution

!===============================================================================

  subroutine ComputeAccelerations( Particle, NP, TimeStep, Update )

    integer(ib), intent(in) :: NP
    type(TParticle), intent(inout) :: Particle(NP)
    real(rb), intent(in) :: TimeStep
    logical, intent(in) :: Update

    integer(ib) :: i, j, k, iCell, ii, jj, NPCC, NPSC

    real(rb) :: Force(NP,dim), Torque(NP,dim)

    type(TParticle) :: iPart, jPart
    type(TWall)     :: jWall

    HalfTimeStep = 0.5_rb*TimeStep

    ! Redistribute the particles throughout the lattice:
    call RedistributeParticles3D( CellIndex, Particle%Position(1),  &
                                             Particle%Position(2),  &
                                             Particle%Position(3) )

    ! Zero the resultant forces and torques:
    Force = 0.0_rb
    Torque = 0.0_rb

    if (FluidFlow) call UpdateVelocityField( Particle )

    ! Examine every cell:
    do iCell = 1_ib, NCell
    !do iCell = Ncell, 1_ib, -1_ib

      ! Create a list of particles in the current and surrounding cells:
      call NeighbourParticles( iCell, NPCC, NPSC, Neighbour )

      ! Consider every particle in the current cell:
      do ii = 1_ib, NPCC

        ! Particle "i" is the ii-th neighbour of the list:
        i = Neighbour(ii)

        ! Save the properties of particle "i" in scalar variables:
        iPart  = Particle(i)
        Ri = Component(iPart%Index)%Radius

        ! Set up the gravitational and buoyance forces:
        iForce = EffectiveMass(iPart%Index)*Gravity

!write(*,'("Fg, Fa = ", ES,"  ")',advance="no") iForce(3)*MassScale*LengthScale/TimeScale**2
!stop
        ! Zero the accumulators for the torque on particle "i":
        iTorque = 0.0_rb

        ! Examine the interaction of Particle "i" with every particle in
        ! the current and surrounding cells that was still not considered:
        do jj = ii + 1_ib, NPSC

          ! Particle "j" is the jj-th neighbour of the list:
          j = Neighbour(jj)

          ! Save the properties of particle "i" in scalar variables:
          jPart = Particle(j)
          Rj = Component(jPart%Index)%Radius

          ! Calculate the relative position of Particles "i" and "j":
          Rij = jPart%Position - iPart%Position

          ! Apply periodic boundary conditions:
          forall (k=1_ib:dim,PBC(k))
            Rij(k) = Rij(k) - Lbox(k)*nint( Rij(k)*InvLbox(k) )
          end forall

          ! Calculate the square of their center-center distance:
          DijSq = sum(Rij*Rij)

          ! Verify if "i" and "j" were contacting in the previous step:
          call VerifyPairContact( i, j, Contacting, ijForceT0, VijT0 )

          ! Verify if the particles are overlapping:
          DeltaDijSq = MinDijSq(iPart%Index,jPart%Index) - DijSq
          if ( DeltaDijSq > 0.0_rb ) then

            !<test>
            if (.not.Contacting) Ncol = Ncol + 1_ib
            !</test>

            ! Calculate the center-center distance between "i" and "j":
            Dij = sqrt(DijSq)

            ! Calculate the unit vector that points from the center of
            ! particle "i" to the center of particle "j":
            Nij = Rij/Dij

            ! Calculate the normal displacement:
            ijDispN = (Dij - MinDij(iPart%Index,jPart%Index))*Nij

            ! Calculate the center-center relative velocity of "i" and "j":
            Vij = jPart%Velocity - iPart%Velocity

            ! Calculate the relative velocity at the contact:
            VijC = Vij - ( (Ri*iPart%Omega + Rj*jPart%Omega) .x. Nij )

            ! Perform calculation:
            call CommonCalculations( Rij,                          &
                                     Kn(iPart%Index,jPart%Index),  &
                                     Yn(iPart%Index,jPart%Index),  &
                                     Kt(iPart%Index,jPart%Index),  &
                                     Yt(iPart%Index,jPart%Index),  &
                                     MuSq(iPart%Index,jPart%Index) )

            ! Add the reciprocal of the contact force to the resultant
            ! force on particle "j":
            Force(j,:) = Force(j,:) - ijForceN - ijForceT

            ! Calculate the torque on particle "j":
            Torque(j,:) = Torque(j,:) + Rj*Nij_x_ijForceT

            ! Update the contact status, if necessary:
            if (Update) call UpdatePairContact( ijForceT, VijT, Contacting )

          else

            ! If "i" and "j" were in contact in the previous step, it is
            ! necessary to remove the pair from the contact list:
            if (Contacting .and. Update) call RemovePairContact

          end if

        end do ! jj

        ! Verify the interaction of Particle "i" with every wall:
        do j = 1_ib, NW

          ! Save the properties of wall "j" in a scalar variable:
          jWall = Wall(j)

          ! Verify if "i" and "j" were in contact in the previous step:
          call VerifyWallContact( i, j, Contacting, ijForceT0, VijT0 )

          ! Verify if Particle "i" and Wall "j" are overlapping:
          if ( Overlap(iPart, jWall, Ri, PBC, Dij, Nij, Qji) ) then

            DeltaDijSq = Ri*Ri - Dij*Dij

            ! Calculate the normal displacement:
            ijDispN = (Dij - Ri)*Nij

            ! Calculate the center-contact relative velocity of "i" and "j":
            Vij = jWall%Velocity + (jWall%Omega .x. Qji) - iPart%Velocity

            ! Calculate the relative velocity at the contact:
            VijC = Vij - (Ri*iPart%Omega .x. Nij)

            ! Calculate forces:
            call CommonCalculations( jWall%Position + Qji - iPart%Position,  &
                                     Knw(iPart%Index,j), Ynw(iPart%Index,j), &
                                     Ktw(iPart%Index,j), Ytw(iPart%Index,j), &
                                     MuwSq(iPart%Index,j)                    )

            ! Update the contact status, if necessary:
            if (Update) call UpdateWallContact( ijForceT, VijT, Contacting )

          else

            ! If "i" and "j" were in contact in the previous step, it is
            ! necessary to remove the pair from the contact list:
            if (Contacting .and. Update) call RemoveWallContact

          end if

        end do ! j

        if (FluidFlow) iForce = iForce + DragForce(iPart)

        ! Update the resultant force and torque on Particle "i":
        Force(i,:) = Force(i,:) + iForce
        Torque(i,:) = Torque(i,:) + iTorque

      end do ! ii

    end do ! iCell

!print*, Force(1,:)*MassScale*LengthScale/TimeScale**2
!print*, Particle(1)%Velocity*LengthScale/TimeScale

    ! Calculate the accelerations of every particle:
    do i = 1_ib, NP
      Particle(i)%Acceleration = Force(i,:)*InvMass(Particle(i)%Index)
      Particle(i)%Alpha = Torque(i,:)*InvInertia(Particle(i)%Index)
    end do

  end subroutine ComputeAccelerations

!===============================================================================

    subroutine CommonCalculations( Rij, Kn, Yn, Kt, Yt, MuSq )

      real(rb), intent(in) :: Rij(dim), Kn, Yn, Kt, Yt, MuSq

      real(rb) :: VijN(dim), MuNormForceNSq, NormForceTSq,  &
                  RijVij, VijSq, PartialDt

      ! Calculate the components of the contact relative velocity:
      VijN = (VijC .dot. Nij)*Nij
      VijT = VijC - VijN

      ! Calculate the normal contact force between "i" and "j":
      ijForceN = Kn*ijDispN + Yn*VijN

      ! Calculate the tangential contact force between "i" and "j":
      if (Contacting) then
        VijT0 = ChangePlane(VijT0,Nij)
        ijDeltaDT = HalfTimeStep*(VijT0 + VijT)
        ijDeltaVT = VijT - VijT0
        ijForceT = ChangePlane(ijForceT0,Nij) + Kt*ijDeltaDT + Yt*ijDeltaVT
      else
        VijSq = sum(Vij*Vij)
        if (VijSq > 0.0_rb) then
          RijVij = Rij .dot. Vij
          PartialDt = min( 2.0_rb*HalfTimeStep,  &
                      (RijVij + sqrt(RijVij*RijVij + VijSq*DeltaDijSq))/VijSq )
        else
          PartialDt = 0.0_rb
        end if
        ijForceT = (Kt*PartialDt + Yt)*VijT
      end if

      ! Calculate the norms of the normal and tangential forces:
      MuNormForceNSq = MuSq*sum(ijForceN**2)
      NormForceTSq = sum(ijForceT**2)

      ! Verify if sliding occurs:
      if (NormForceTSq > MuNormForceNSq) then
        ijForceT = sqrt(MuNormForceNSq/NormForceTSq)*ijForceT
      end if

      ! Add the contact force to the accumulator for the force on "i":
      iForce = iForce + ijForceN + ijForceT

      ! Calculate the torque caused by the tangential force:
      Nij_x_ijForceT = Nij .x. ijForceT
      iTorque = iTorque + (Ri*Nij .x. ijForceT)

    end subroutine CommonCalculations

!===============================================================================

  function ChangePlane( Vector, Normal ) result ( NewVector )

    real(rb), intent(in) :: Vector(dim), Normal(dim)
    real(rb) :: NewVector(dim)

    real(rb) :: Projection(dim), NormSq

    Projection = Vector - (Vector .dot. Normal)*Normal
    NormSq = sum(Projection**2)
    if ( NormSq /= 0.e0_rb ) then
      NewVector = sqrt(sum(Vector**2)/NormSq)*Projection
    else
      NewVector = 0.e0_rb
    end if

  end function

!===============================================================================

end module MForces
