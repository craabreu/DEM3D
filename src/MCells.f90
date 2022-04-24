module MCells

use MTypes
use BinFiles

implicit none

!      This module contains subroutines and functions for handling lists of
! neighbour particles based on the cell structure algorithm (Allen & Tildesley,
! p.149) for Monte Carlo and Molecular Dynamics simulations.
!
! Author...: Charlles R. A. Abreu
! Location.: Escola de Química - Universidade Federal do Rio de Janeiro - Brazil
! E-mail...: rubber@eq.ufrj.br
! Date.....: May 22, 2003

! Defined types:

type, private :: TCellParticle
  integer(ib) :: Index
  type (TCellParticle), pointer :: Next => null()
end type

type, private :: TCell
  type (TCellParticle), pointer :: FirstParticle  => null()
  type (TCellParticle), pointer :: FirstNeighbour => null()
end type

! Module constants:

integer(ib), parameter :: MonteCarlo        = 1_ib,  &
                          MolecularDynamics = 2_ib

! Module variables:

integer(ib), private :: Nx, Ny, Nz, Nt, Nxm1, Nym1, Nzm1

real(rb), private    :: NxByLx, NyByLy, NzByLz,   &
                        NxBy2,  NyBy2,  NzBy2

!integer(ib), private :: MaxNeighbours

type (TCell), allocatable, private :: Cell(:)

type (TCellParticle), pointer, private :: Protected

! Interfaces for generic procedures:

interface InitCellStructure
  module procedure InitCellStructure2D,  &
                   InitCellStructure3D
end interface InitCellStructure

interface FullNeighbourCells
  module procedure FullNeighbourCells2D,  &
                   FullNeighbourCells3D
end interface FullNeighbourCells

interface NeighbourCells
  module procedure NeighbourCells2D,  &
                   NeighbourCells3D
end interface NeighbourCells

interface CurrentCell
  module procedure CurrentCell2D,  &
                   CurrentCell3D
end interface CurrentCell

interface DistributeParticles
  module procedure DistributeParticles2D,  &
                   DistributeParticles3D
end interface DistributeParticles

interface RedistributeParticles
  module procedure RedistributeParticles2D,  &
                   RedistributeParticles3D
end interface RedistributeParticles

! Local procedures:

private :: VerifyPBC, ApplyPBC

contains

!==============================================================================

  subroutine InitCellStructure2D( MCellx, MCelly, NCell, Lx, Ly,  &
                                  CutOff, MaxNCell )

    integer(ib), intent(out) :: MCellx, MCelly, NCell
    real(rb), intent(in)     :: Lx, Ly, CutOff
    integer(ib), intent(in)  :: MaxNCell

    ! This subroutine must be called before any other procedure in this module.
    !
    ! Lx       - Dimension of the simulation box in the direction x
    ! Ly       - Dimension of the simulation box in the direction y
    ! CutOff   - Maximum cut-off distance of the particles
    ! MaxNCell - Maximum allowed number of cells

    real(rb)    :: LCell
    integer(ib) :: Nx4, Ny4, N4

    ! Predefine the number of cells in every direction:
    Nx4 = int(Lx/CutOff,ib)
    Ny4 = int(Ly/CutOff,ib)

    ! Predefine the total number of cells:
    N4 = Nx4*Ny4

    ! Verify if the maximum allowable number of cells was exceeded:
    if ( N4 > MaxNCell ) then

      ! Define the minimum side length allowed for the cells:
      LCell = sqrt( Lx*Ly/real(MaxNCell,rb) )

      ! Redefine the number of cells in every direction:
      Nx = int(Lx/LCell,ib)
      Ny = int(Ly/LCell,ib)

      ! Predefine the total number of cells:
      Nt = Nx*Ny

    else

      ! Confirm the numbers of cells:
      Nx = Nx4
      Ny = Ny4
      Nt = N4

    end if

    ! Define auxiliary variables for forthcoming calculations:
    NxByLx = real(Nx,rb)/Lx
    NyByLy = real(Ny,rb)/Ly

    NxBy2  = real(Nx,rb)*0.5_rb
    NyBy2  = real(Ny,rb)*0.5_rb

    Nxm1   = Nx - 1_ib
    Nym1   = Ny - 1_ib

    ! Initialize variables for handling cell lists:
    call AllocateCellList( Nt )

    ! Setup output variables:
    Mcellx = Nx
    MCelly = Ny
    NCell  = Nt

  end subroutine InitCellStructure2D

!-------------------------------------------------------------------------------

  subroutine InitCellStructure3D( MCellx, MCelly, Mcellz, NCell,  &
                                  Lx, Ly, Lz, CutOff, MaxNCell    )

    integer(ib), intent(out) :: MCellx, MCelly, Mcellz, NCell
    real(rb), intent(in)    :: Lx, Ly, Lz, CutOff
    integer(ib), intent(in) :: MaxNCell

    ! This subroutine must be called before any other procedure in this module.
    !
    ! Lx       - Dimension of the simulation box in the direction x
    ! Ly       - Dimension of the simulation box in the direction y
    ! Lz       - Dimension of the simulation box in the direction y
    ! CutOff   - Maximum cut-off distance of the particles
    ! MaxNCell - Maximum allowed number of cells

    real(rb)    :: LCell
    integer(ib) :: Nx4, Ny4, Nz4, N4

    ! Predefine the number of cells in every direction:
    Nx4 = int(Lx/CutOff,ib)
    Ny4 = int(Ly/CutOff,ib)
    Nz4 = int(Lz/CutOff,ib)

    ! Predefine the total number of cells:
    N4 = Nx4*Ny4*Nz4

    ! Verify if the maximum allowable number of cells was exceeded:
    if ( N4 > MaxNCell ) then

      ! Define the minimum side length allowed for the cells:
      LCell = ( Lx*Ly*Lz/real(MaxNCell,rb) )**(1.0_rb/3.0_rb)

      ! Redefine the number of cells in every direction:
      Nx = int(Lx/LCell,ib)
      Ny = int(Ly/LCell,ib)
      Nz = int(Lz/LCell,ib)

      ! Predefine the total number of cells:
      Nt = Nx*Ny*Nz

    else

      ! Confirm the numbers of cells:
      Nx = Nx4
      Ny = Ny4
      Nz = Nz4
      Nt = N4

    end if

    ! Define auxiliary variables for forthcoming calculations:
    NxByLx = real(Nx,rb)/Lx
    NyByLy = real(Ny,rb)/Ly
    NzByLz = real(Nz,rb)/Lz

    NxBy2  = real(Nx,rb)*0.5_rb
    NyBy2  = real(Ny,rb)*0.5_rb
    NzBy2  = real(Nz,rb)*0.5_rb

    Nxm1   = Nx - 1_ib
    Nym1   = Ny - 1_ib
    Nzm1   = Nz - 1_ib

    ! Initialize variables for handling cell lists:
    call AllocateCellList( Nt )

    ! Setup output variables:
    Mcellx = Nx
    MCelly = Ny
    Mcellz = Nz
    NCell  = Nt

  end subroutine InitCellStructure3D

!==============================================================================

  subroutine FullNeighbourCells2D( Method, PBCx, PBCy )

    integer(ib), intent(in) :: Method
    logical, intent(in) :: PBCx, PBCy

    !     This routine creates a list of the neighbour cells of every cell.
    !  The simulation method, as well as the presence of rigid walls and/or
    !  periodic boundary conditions are considered. For Molecular Dynamics
    !  calculations, only one cell of each pair of opposed neighbours is
    !  accounted for, because Newton's Third Law can be used to improve
    !  efficiency.
    !
    !  Method - MonteCarlo or MolecularDynamics
    !  PBCx   - .true. if periodic boundary conditions apply inthe direction x
    !  PBCy   - .true. if periodic boundary conditions apply inthe direction y

    integer(ib) :: ix, iy, iCell, xi, xf, yi, yf
    logical     :: MC

    MC = Method == MonteCarlo

    do iy = 0_ib , Nym1

      ! Verify if periodic boundary conditions apply to the direction y:
      call VerifyPBC( yi, yf, iy, Ny, PBCy )

      do ix = 0_ib , Nxm1

        ! Verify if periodic boundary conditions apply to the direction x:
        call VerifyPBC( xi, xf, ix, Nx, PBCx )

        ! Calculate the index of the current central cell:
        iCell = 1_ib + ix + iy*Nx

        ! Add all the neighbour cells to the list of iCell:
        call AddAllNeighbours2D( iCell, ix, iy, PBCx, PBCy, xi, xf, yi, yf, MC )

      end do ! ix

    end do ! iy

  end subroutine FullNeighbourCells2D

!-------------------------------------------------------------------------------

  subroutine FullNeighbourCells3D( Method, PBCx, PBCy, PBCz )

    integer(ib), intent(in) :: Method
    logical :: PBCx, PBCy, PBCz

    !     This routine creates a list of the neighbour cells of every cell.
    !  The simulation method, as well as the presence of rigid walls and/or
    !  periodic boundary conditions are considered. For Molecular Dynamics
    !  calculations, only one cell of each pair of opposed neighbours is
    !  accounted for, because Newton's Third Law can be used to improve
    !  efficiency.
    !
    !  Method - MonteCarlo or MolecularDynamics
    !  PBCx   - .true. if periodic boundary conditions apply inthe direction x
    !  PBCy   - .true. if periodic boundary conditions apply inthe direction y
    !  PBCz   - .true. if periodic boundary conditions apply inthe direction z

    integer(ib) :: ix, iy, iz, iCell, xi, xf, yi, yf, zi, zf
    logical     :: MC

    MC = Method == MonteCarlo

    ! Pass through all the cells:
    do iz = 0_ib , Nzm1

      ! Verify if periodic boundary conditions apply to the direction z:
      call VerifyPBC( zi, zf, iz, Nz, PBCz )

      do iy = 0_ib , Nym1

        ! Verify if periodic boundary conditions apply to the direction y:
        call VerifyPBC( yi, yf, iy, Ny, PBCy )

        do ix = 0_ib , Nxm1

          ! Verify if periodic boundary conditions apply to the direction x:
          call VerifyPBC( xi, xf, ix, Nx, PBCx )

          ! Calculate the index of the current central cell:
          iCell = 1_ib + ix + Nx*(iy + Ny*iz)

          ! Add all the neighbour cells to the list of iCell:
          call AddAllNeighbours3D( iCell, ix, iy, iz, PBCx, PBCy, PBCz,  &
                                   xi, xf, yi, yf, zi, zf, MC )

        end do ! ix

      end do ! iy

    end do ! iz

  end subroutine FullNeighbourCells3D

! ==============================================================================

  subroutine NeighbourCells2D( iCell, Method, PBCx, PBCy )

    integer(ib), intent(in) :: iCell, Method
    logical, intent(in) :: PBCx, PBCy

    integer(ib) :: ix, iy, xi, xf, yi, yf

    iy = (iCell - 1_ib)/Nx
    ix = (iCell - 1_ib) - iy*Nx

    ! Verify if periodic boundary conditions apply:
    call VerifyPBC( xi, xf, ix, Nx, PBCx )
    call VerifyPBC( yi, yf, iy, Ny, PBCy )

    ! Add all the neighbour cells to the list of iCell:
    call AddAllNeighbours2D( iCell, ix, iy, PBCx, PBCy, xi, xf, yi, yf,  &
                             Method == MonteCarlo )

  end subroutine NeighbourCells2D

!-------------------------------------------------------------------------------

  subroutine NeighbourCells3D( iCell, Method, PBCx, PBCy, PBCz )

    integer(ib), intent(in) :: iCell, Method
    logical, intent(in) :: PBCx, PBCy, PBCz

    integer(ib) :: ix, iy, iz, xi, xf, yi, yf, zi, zf, aux1, aux2

    aux1 = iCell - 1_ib
    iz   = aux1/(Nx*Ny)
    aux2 = iz*Ny
    iy   = aux1/Nx - aux2
    ix   = aux1 - (iy + aux2)*Nx

    ! Verify if periodic boundary conditions apply:
    call VerifyPBC( xi, xf, ix, Nx, PBCx )
    call VerifyPBC( yi, yf, iy, Ny, PBCy )
    call VerifyPBC( zi, zf, iz, Nz, PBCz )

    ! Add all the neighbour cells to the list of iCell:
    call AddAllNeighbours3D( iCell, ix, iy, iz, PBCx, PBCy, PBCz,  &
                             xi, xf, yi, yf, zi, zf, Method == MonteCarlo )

  end subroutine NeighbourCells3D

!===============================================================================

  subroutine AddAllNeighbours2D( iCell, ix, iy, PBCx, PBCy, xi, xf, yi, yf, MC )

    integer(ib), intent(in) :: iCell, ix, iy, xi, xf, yi, yf
    logical, intent(in) :: PBCx, PBCy
    logical, intent(in)     :: MC

    integer(ib) :: jx, jy, kx, ky, jCell, kCell

    ! Define variations of y around the central cell:
    do ky = yi, yf

      ! Define the coordinate y of the current neighbour cell:
      jy = iy + ky

      ! Apply or not the periodic bounds in the direction y:
      call ApplyPBC( jy, Ny, PBCy )

      ! Define variations of x around the central cell:
      do kx = xi, xf

        kCell = kx + ky*Nx

        if ( (MC.and.(kCell /= 0_ib)).or.(kCell > 0_ib) ) then

          ! Define the coordinate x of the current neighbour cell:
          jx = ix + kx

          ! Apply or not the periodic bounds in the direction x:
          call ApplyPBC( jx, Nx, PBCx )

          ! Calculate the index of the current neighbour cell:
          jCell = 1_ib + jx + jy*Nx

          ! Add jCell to the list of neighbours of iCell:
          call AddNeighbour( jCell, iCell )

        end if

      end do ! kx
    end do ! ky

  end subroutine AddAllNeighbours2D

!-------------------------------------------------------------------------------

  subroutine AddAllNeighbours3D( iCell, ix, iy, iz, PBCx, PBCy, PBCz,  &
                                 xi, xf, yi, yf, zi, zf, MC )

    integer(ib), intent(in) :: iCell, ix, iy, iz, xi, xf, yi, yf, zi, zf
    logical, intent(in) :: PBCx, PBCy, PBCz
    logical, intent(in)     :: MC

    integer(ib) :: jx, jy, jz, kx, ky, kz, jCell, kCell

    ! Define variations of z around the central cell:
    do kz = zi, zf

      ! Define the coordinate z of the current neighbour cell:
      jz = iz + kz

      ! Apply or not the periodic bounds in the direction z:
      call ApplyPBC( jz, Nz, PBCz )

      ! Define variations of y around the central cell:
      do ky = yi, yf

        ! Define the coordinate y of the current neighbour cell:
        jy = iy + ky

        ! Apply or not the periodic bounds in the direction y:
        call ApplyPBC( jy, Ny, PBCy )

        ! Define variations of x around the central cell:
        do kx = xi, xf

          kCell = kx + Nx*(ky + Ny*kz)

          if ( (MC.and.(kCell /= 0_ib)).or.(kCell > 0_ib) ) then

            ! Define the coordinate x of the current neighbour cell:
            jx = ix + kx

            ! Apply or not the periodic bounds in the direction x:
            call ApplyPBC( jx, Nx, PBCx )

            ! Calculate the index of the current neighbour cell:
            jCell = 1_ib + jx + Nx*(jy + Ny*jz)

            ! Add jCell to the list of neighbours of iCell:
            call AddNeighbour( jCell, iCell )

          end if

        end do ! kx
      end do ! ky
    end do ! kz

  end subroutine AddAllNeighbours3D

! ==============================================================================

  elemental function CurrentCell2D( x, y ) result ( iCell )

    real(rb), intent(in) :: x, y
    integer(ib) :: iCell

    ! This function returns the index of the cell in which the point (x,y) is. 

    integer(ib) :: ix, iy

    ix = int(NxByLx*x + NxBy2,ib)
    iy = int(NyByLy*y + NyBy2,ib)

    iCell = 1_ib + ix + iy*Nx

  end function CurrentCell2D

!-------------------------------------------------------------------------------

  elemental function CurrentCell3D( x, y, z ) result ( iCell )

    real(rb), intent(in) :: x, y, z
    integer(ib) :: iCell

    ! This function returns the index of the cell in which the point (x,y,z) is.

    integer(ib) :: ix, iy, iz

    ix = max( 0_ib, min( int(NxByLx*x + NxBy2,ib), Nxm1 ) )
    iy = max( 0_ib, min( int(NyByLy*y + NyBy2,ib), Nym1 ) )
    iz = max( 0_ib, min( int(NzByLz*z + NzBy2,ib), Nzm1 ) )

    iCell = 1_ib + ix + (iy + iz*Ny)*Nx

  end function CurrentCell3D

!===============================================================================

  subroutine DistributeParticles2D( CellIndex, Rx, Ry )

    integer(ib), intent(out) :: CellIndex(:)
    real(rb), intent(in)     :: Rx(:), Ry(:)

    !    This subroutine distribute all the particles throughout the cell struc-
    ! ture, based on the location of the center of mass of each particle.
    !
    ! Rx - Vector of the coordinates x of the centers of mass of the particles
    ! Ry - Vector of the coordinates y of the centers of mass of the particles

    integer(ib) :: i, Np

    Np = size(Rx)

    CellIndex(1:Np) = CurrentCell2D( Rx, Ry )

    do i = 1 , Np
      call AddParticle( i, CellIndex(i) )
    end do

  end subroutine DistributeParticles2D

!-------------------------------------------------------------------------------

  subroutine DistributeParticles3D( CellIndex, Rx, Ry, Rz )

    integer(ib), intent(out) :: CellIndex(:)
    real(rb), intent(in)     :: Rx(:), Ry(:), Rz(:)

    !    This subroutine distribute all the particles throughout the cell struc-
    ! ture, based on the location of the center of mass of each particle.
    !
    ! Rx - Vector of the coordinates x of the centers of mass of the particles
    ! Ry - Vector of the coordinates y of the centers of mass of the particles
    ! Rz - Vector of the coordinates z of the centers of mass of the particles

    integer(ib) :: i, Np

    Np = size(Rx)

    CellIndex(1:Np) = CurrentCell3D( Rx, Ry, Rz )

    do i = 1 , Np
      call AddParticle( i, CellIndex(i) )
    end do

  end subroutine DistributeParticles3D

!===============================================================================

  subroutine RedistributeParticles2D( CellIndex, Rx, Ry )

    integer(ib), intent(inout) :: CellIndex(:)
    real(rb), intent(in)       :: Rx(:), Ry(:)

    !    This subroutine distribute all the particles throughout the cell struc-
    ! ture, based on the location of the center of mass of each particle.
    !
    ! Rx - Vector of the coordinates x of the centers of mass of the particles
    ! Ry - Vector of the coordinates y of the centers of mass of the particles

    integer(ib) :: i, Old, New

    do i = 1 , size(Rx)
      Old = CellIndex(i)
      New = CurrentCell2D( Rx(i), Ry(i) )
      if ( New /= Old ) then
        call TransferParticle( i, Old, New )
        CellIndex(i) = New
      end if
    end do

!    integer(ib) :: i, NewCellIndex(size(Rx)), Old, New
!
!    NewCellIndex = CurrentCell2D( Rx, Ry )
!
!    do i = 1 , size(Rx)
!      Old = CellIndex(i)
!      New = NewCellIndex(i)
!      if ( New /= Old ) then
!        call TransferParticle( i, Old, New )
!        CellIndex(i) = New
!      end if
!    end do

  end subroutine RedistributeParticles2D

!-------------------------------------------------------------------------------

  subroutine RedistributeParticles3D( CellIndex, Rx, Ry, Rz )

    integer(ib), intent(inout) :: CellIndex(:)
    real(rb), intent(in)       :: Rx(:), Ry(:), Rz(:)

    !    This subroutine distribute all the particles throughout the cell struc-
    ! ture, based on the location of the center of mass of each particle.
    !
    ! Rx - Vector of the coordinates x of the centers of mass of the particles
    ! Ry - Vector of the coordinates y of the centers of mass of the particles
    ! Rz - Vector of the coordinates z of the centers of mass of the particles

    integer(ib) :: i, Old, New

    do i = 1 , size(Rx)
      Old = CellIndex(i)
      New = CurrentCell3D( Rx(i), Ry(i), Rz(i) )
      if ( New /= Old ) then
        call TransferParticle( i, Old, New )
        CellIndex(i) = New
      end if
    end do

!    integer(ib) :: i, NewCellIndex(size(Rx)), Old, New
!
!    NewCellIndex = CurrentCell3D( Rx, Ry, Rz )
!
!    do i = 1 , size(Rx)
!      Old = CellIndex(i)
!      New = NewCellIndex(i)
!      if ( New /= Old ) then
!        call TransferParticle( i, Old, New )
!        CellIndex(i) = New
!      end if
!    end do

  end subroutine RedistributeParticles3D

!===============================================================================

  subroutine DiscardNeighbourCells( iCell )

    integer(ib), intent(in) :: iCell

    type (TCellParticle), pointer :: Aux, Current

    ! Destruct Lists of neighbour cells of iCell:
    Current => Cell(iCell)%FirstNeighbour
    do while ( associated(Current) )
      Aux => Current
      Current => Current%Next
      deallocate( Aux )
    end do
    Cell(iCell)%FirstNeighbour => null()

  end subroutine DiscardNeighbourCells

! ==============================================================================

  subroutine AddNeighbour( jCell, iCell )

    integer(ib), intent(in) :: jCell, iCell

    !    This subroutine adds a cell index "jCell" to the list of neighbour
    !  cells of another cell "iCell".

    type (TCellParticle), pointer :: Current

    allocate( Current )
    Current%Index = jCell
    Current%Next  => Cell(iCell)%FirstNeighbour
    Cell(iCell)%FirstNeighbour => Current

  end subroutine AddNeighbour

! ==============================================================================

  subroutine AddParticle( i, iCell )

    integer(ib), intent(in) :: i, iCell

    !    This subroutine adds a particle index "i" to the list of particles
    ! that lay in the cell "iCell".

    type (TCellParticle), pointer :: Current

    allocate( Current )
    Current%Index = i
    Current%Next  => Cell(iCell)%FirstParticle
    Cell(iCell)%FirstParticle => Current

  end subroutine AddParticle

! ==============================================================================

  subroutine RemoveParticleFromCell( i, iCell )

    integer(ib), intent(in) :: i, iCell

    !    This subroutine removes a particle index "i" from the list of particles
    ! that lay in the cell "iCell".

    type (TCellParticle), pointer :: Aux, Current

    Current => Cell(iCell)%FirstParticle
    if ( Current%Index == i ) then
      Cell(iCell)%FirstParticle => Current%Next
      deallocate( Current )
    else
      do while ( Current%Next%Index /= i )
        Current => Current%Next
      end do
      Aux => Current%Next
      Current%Next => Aux%Next
      deallocate( Aux )
    end if

  end subroutine RemoveParticleFromCell

! ==============================================================================

  subroutine ReplaceParticleIndexInCell( i, iCell, j )

    integer(ib), intent(in) :: i, iCell, j

    !   This subroutine replaces the index "i" in the cell "iCell"
    ! by the index "j".

    type (TCellParticle), pointer :: Current

    Current => Cell(iCell)%FirstParticle
    do while (Current%Index /= i)
      Current => Current%Next
    end do
    Current%Index = j

  end subroutine ReplaceParticleIndexInCell

! ==============================================================================

  subroutine TransferParticle( i, iCell, jCell )

    integer(ib), intent(in) :: i, iCell, jCell

    !    This subroutine transfers a particle "i" from the list of particles of a
    ! cell "iCell" to the list of particles of another cell "jCell".

    type (TCellParticle), pointer :: Aux, Current

    ! Remove particle "i" from iCell:
    Current => Cell(iCell)%FirstParticle

    if ( Current%Index == i ) then
      Cell(iCell)%FirstParticle => Current%Next
      Aux => Current
    else
      do while ( Current%Next%Index /= i )
        Current => Current%Next
      end do
      Aux => Current%Next
      Current%Next => Aux%Next
    end if

    ! Add particle "i" to jCell:
    Aux%Next => Cell(jCell)%FirstParticle
    Cell(jCell)%FirstParticle => Aux

  end subroutine TransferParticle

! ==============================================================================

  subroutine AllocateCellList( NCell )

    integer(ib), intent(in) :: NCell

    ! Allocate the array of pointers Cell:
    if (allocated(Cell)) then
      write(*,'("Warning: Cell structure already existed.")')
      deallocate(Cell)
    end if
    allocate( Cell(NCell) )

  end subroutine AllocateCellList

! ==============================================================================

  subroutine NeighbourParticles( iCell, NumberInCell, &
                                 NumberOfNeighbours,  &
                                 Neighbour            )

    integer(ib), intent(in)  :: iCell
    integer(ib), intent(out) :: NumberInCell, NumberOfNeighbours, Neighbour(:)

    type (TCellParticle), pointer :: CurrentParticle, CurrentCell

    NumberInCell = 0_ib

    CurrentParticle => Cell(iCell)%FirstParticle

    do while ( associated(CurrentParticle) )
      NumberInCell = NumberInCell + 1_ib
      Neighbour(NumberInCell) = CurrentParticle%Index
      CurrentParticle => CurrentParticle%Next
    end do

    if (NumberInCell /= 0_ib) then

      NumberOfNeighbours = NumberInCell

      CurrentCell => Cell(iCell)%FirstNeighbour
      do while ( associated(CurrentCell) )
        CurrentParticle => Cell(CurrentCell%Index)%FirstParticle
        do while ( associated(CurrentParticle) )
          NumberOfNeighbours = NumberOfNeighbours + 1_ib
          Neighbour(NumberOfNeighbours) = CurrentParticle%Index
          CurrentParticle => CurrentParticle%Next
        end do
        CurrentCell => CurrentCell%Next
      end do

    end if

  end subroutine NeighbourParticles

! ==============================================================================

  subroutine AllNeighbourParticles( iCell, NumberOfNeighbours, Neighbour )

    integer(ib), intent(in)  :: iCell
    integer(ib), intent(out) :: NumberOfNeighbours, Neighbour(:)

    type (TCellParticle), pointer :: Particle, CurrentCell

    NumberOfNeighbours = 0_ib

    Particle => Cell(iCell)%FirstParticle
    do while ( associated(Particle) )
      NumberOfNeighbours = NumberOfNeighbours + 1_ib
      Neighbour(NumberOfNeighbours) = Particle%Index
      Particle => Particle%Next
    end do

    CurrentCell => Cell(iCell)%FirstNeighbour
    do while ( associated(CurrentCell) )
      Particle => Cell(CurrentCell%Index)%FirstParticle
      do while ( associated(Particle) )
        NumberOfNeighbours = NumberOfNeighbours + 1_ib
        Neighbour(NumberOfNeighbours) = Particle%Index
        Particle => Particle%Next
      end do
      CurrentCell => CurrentCell%Next
    end do

  end subroutine AllNeighbourParticles

! ==============================================================================

  subroutine DestroyCellStructure

    integer(ib) :: iCell
    type (TCellParticle), pointer :: Current, Aux

    do iCell = 1 , Nt

      ! Destruct Lists of neighbour cells:
      Current => Cell(iCell)%FirstNeighbour
      do while ( associated(Current) )
        Aux => Current
        Current => Current%Next
        deallocate( Aux )
      end do
      Cell(iCell)%FirstNeighbour => null()

      ! Destruct lists of particles:
      Current => Cell(iCell)%FirstParticle
      do while ( associated(Current) )
        Aux => Current
        Current => Current%Next
        deallocate( Aux )
      end do
      Cell(iCell)%FirstParticle => null()

    end do ! iCell

    ! Deallocate cell list:
    deallocate( Cell )

  end subroutine DestroyCellStructure

!==============================================================================-

  subroutine PrintCellStructure( Unit, Np )

    integer, intent(in)     :: Unit
    integer(ib), intent(in) :: Np

    integer(ib) :: iCell, nsc, nsp
    character   :: sc, sp
    type (TCellParticle), pointer :: Current

    nsc = 1 + int(log10(real(Nt,rb)),ib)
    write(sc,'(I1)') nsc

    nsp = 1 + int(log10(real(Np,rb)),ib)
    write(sp,'(I1)') nsp

    write(Unit,*)
    write(Unit,'(" ** Structure of cells **")')
    write(Unit,*)
    write(Unit,'("Number of cells in the direction x: ",I'//sc//')') Nx
    write(Unit,'("Number of cells in the direction y: ",I'//sc//')') Ny
    write(Unit,'("Number of cells in the direction z: ",I'//sc//')') Nz
    write(Unit,*)
    write(Unit,'("Total number of cells: ",I'//sc//')') Nt
    write(Unit,*)

    do iCell = 1 , Nt

      write(Unit,'("Cell [",I'//sc//',"]: ")') iCell

      ! Print neighbour cells:
      write(Unit,'("   Neighbours: ")',advance='no')
      Current => Cell(iCell)%FirstNeighbour
      if (associated(Current)) then
        do while ( associated(Current%Next) )
          write(Unit,'(I'//sc//',", ")',advance='no') Current%Index
          Current => Current%Next
        end do
        write(Unit,'(I'//sc//')',advance='no') Current%Index
      end if
      write(Unit,*)

      ! Print particles:
      write(Unit,'("   Particles : ")',advance='no')
      Current => Cell(iCell)%FirstParticle
      if (associated(Current)) then
        do while ( associated(Current%Next) )
          write(Unit,'(I'//sp//',", ")',advance='no') Current%Index
          Current => Current%Next
        end do
        write(Unit,'(I'//sp//')',advance='no') Current%Index
      end if
      write(Unit,*)
      write(Unit,*)

    end do

  end subroutine PrintCellStructure

!===============================================================================

  subroutine SaveCellStructure( Unit )

    integer, intent(in) :: Unit

    integer(ib) :: iCell, NumberOfNeighbours, jCell
    type (TCellParticle), pointer :: Current

    call BinWrite( Unit, Nt )

    do iCell = 1 , Nt

      Current => Cell(iCell)%FirstParticle
      NumberOfNeighbours = 0
      do while (associated(Current))
        Current => Current%Next
        NumberOfNeighbours = NumberOfNeighbours + 1_ib
      end do
      call BinWrite( Unit, NumberOfNeighbours )
      Current => Cell(iCell)%FirstParticle
      do jCell = 1, NumberOfNeighbours
        call BinWrite( Unit, Current%Index )
        Current => Current%Next
      end do

    end do

  end subroutine SaveCellStructure

!===============================================================================

  subroutine ReadCellStructure( Unit )

    integer, intent(in) :: Unit

    integer(ib) :: MyNt, iCell, NumberOfNeighbours, jCell
    type (TCellParticle), pointer :: Current

    call BinRead( Unit, MyNt )
    if (MyNt /= Nt) then
      print*, 'My Nt = ', MyNt
      print*, 'Nt = ', Nt
      write(*,*) 'Error reading cell structure: numbers of cells do not match.'
      stop
    end if

    do iCell = 1 , Nt

      call BinRead( Unit, NumberOfNeighbours )

      if (NumberOfNeighbours > 0_ib) then

        allocate(Current)
        Cell(iCell)%FirstParticle => Current
        do jCell = 1_ib, NumberOfNeighbours - 1_ib
          call BinRead( Unit, Current%Index )
          allocate(Current%Next)
          Current => Current%Next
        end do
        call BinRead( Unit, Current%Index )

      else

        Cell(iCell)%FirstParticle => null()

      end if

    end do

  end subroutine ReadCellStructure

!===============================================================================

  subroutine VerifyPBC( ini, fin, i, MCell, PBC )
    integer(ib), intent(out) :: ini, fin
    integer(ib), intent(in)  :: i, MCell
    logical, intent(in)  :: PBC
    ini = -1_ib
    fin = +1_ib
    if (.not.PBC) then
      if ( i == 0_ib ) ini = 0_ib
      if ( i == MCell - 1_ib ) fin = 0_ib
    end if
  end subroutine VerifyPBC

!===============================================================================

  subroutine ApplyPBC( j, MCell, PBC )
    integer(ib), intent(inout) :: j
    integer(ib), intent(in) :: MCell
    logical, intent(in) :: PBC
    if ( PBC ) then
      if ( j == MCell ) j = 0_ib
      if ( j == -1_ib ) j = MCell - 1_ib
    end if
  end subroutine ApplyPBC

!===============================================================================

  subroutine HinderNeighbourCellsList( iCell )

    integer(ib), intent(in) :: iCell

    Protected => Cell(iCell)%FirstNeighbour
    Cell(iCell)%FirstNeighbour => null()

  end subroutine HinderNeighbourCellsList

!===============================================================================

  subroutine RecoverNeighbourCellsList( iCell )

    integer(ib), intent(in) :: iCell

    call DiscardNeighbourCells( iCell )

    Cell(iCell)%FirstNeighbour => Protected

  end subroutine RecoverNeighbourCellsList

!===============================================================================

end module MCells