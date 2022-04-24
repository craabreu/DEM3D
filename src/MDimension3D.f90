module MDimension

! Flat Walls: Shape = Rectangle
!    L1 - Half width
!    L2 - Half height
!    U  - Direction of the width
!    W  - Direction of the height
!
! Cylindrical Wall: Shape = Cylinder
!    L1 - Half height
!    L2 - Radius
!    U  - Direction of the main axis
!
! Annular Wall: Shape = Annulus
!    L1 - External radius
!    L2 - Internal radius
!    N  - Normal direction

use MTypes
use BinFiles

implicit none

! Dimension of the systems:
integer, parameter :: dim = 3

! Derived type concernig the degrees of freedom of particles:
type TParticle
  integer(ib) :: Index                ! Index for particle type
  real(rb)    :: Position(dim),     & ! Position of the center of mass
                 Velocity(dim),     & ! Velocity of the center of mass
                 Acceleration(dim), & ! Acceleration of the center of mass
                 Omega(dim),        & ! Angular velocity around the center of mass
                 Alpha(dim)           ! Angular Acceleration around the center of mass
end type TParticle

! Derived type concernig direction vectors of the particles:
type TDirection
  real(rb) :: U(dim), &
              W(dim), &
              N(dim)
end type TDirection

! Derived type concernig the degrees of freedom of walls:
type TWall
  integer(ib) :: Shape,          & ! Index for the type of the wall
                 Order             ! Order in the set of equal shaped walls
  real(rb)    :: L1,             & ! First length parameter
                 L2,             & ! Second length parameter
                 L3,             & ! Third length parameter
                 Stiffness,      & ! Wall stiffness
                 ShearStiffness, & ! Wall shear stiffness 
                 U(dim),         & ! Directional unit vector
                 W(dim),         & ! Directional unit vector
                 N(dim),         & ! Directional unit vector (N = U x W)
                 Position(dim),  & ! Position
                 Velocity(dim),  & ! Velocity
                 Omega(dim)        ! Angular velocity
end type TWall

! Definition of dot (inner) product:
interface operator (.dot.)
  module procedure DotProduct
end interface operator (.dot.)

! Definition of cross (vector) product:
interface operator (.x.)
  module procedure CrossProduct
end interface operator (.x.)

contains

!===============================================================================

  function DotProduct( a, b ) result( c )

    real(rb), intent(in) :: a(dim), b(dim)
    real(rb) :: c

    c = dot_product(a,b)

  end function DotProduct

!===============================================================================

  function CrossProduct( a, b ) result( c )

    real(rb), intent(in) :: a(dim), b(dim)
    real(rb) :: c(dim)

    c = (/ a(2)*b(3)-a(3)*b(2), &
           a(3)*b(1)-a(1)*b(3), &
           a(1)*b(2)-a(2)*b(1) /)

  end function CrossProduct

!===============================================================================

  subroutine WriteAxialVector( Unit, a, IOStat )

    integer, intent(in)  :: Unit
    real(rb), intent(in) :: a(dim)
    integer, intent(out) :: IOStat

    call BinWrite( Unit, a, dim, IOStat )

  end subroutine WriteAxialVector

!===============================================================================

  subroutine ReadAxialVector( Unit, a, IOStat )

    integer, intent(in)   :: Unit
    real(rb), intent(out) :: a(dim)
    integer, intent(out)  :: IOStat

    call BinRead( Unit, a, dim, IOStat )

  end subroutine ReadAxialVector

!===============================================================================

  function LatticeIndex( ic, Mcell ) result ( iCell )

    integer(ib), intent(in) :: ic(dim), MCell(dim)
    integer(ib) :: iCell

    iCell = ic(1) + (ic(2) + ic(3)*MCell(2))*MCell(1)

  end function LatticeIndex

!===============================================================================

  function LatticeCoordinates( iCell, MCell ) result ( ic )

    integer(ib), intent(in) :: iCell, MCell(dim)
    integer(ib) :: ic(dim)

    integer(ib) :: Aux

    Aux = MCell(1)*MCell(2)

    ic(3) = iCell/Aux
    ic(2) = (iCell-ic(3)*Aux)/Mcell(1)
    ic(1) = iCell - ic(3)*Aux - ic(2)*MCell(1)

  end function LatticeCoordinates

!===============================================================================

  function CartToSpheric( U ) result ( phi )

    real(rb), intent(in) :: U(dim)
    real(rb) :: phi(dim-1)

    ! U(1) = cos(phi(1))*cos(phi(2))
    ! U(2) = sin(phi(1))*cos(phi(2))
    ! U(3) = sin(phi(2))

    real(rb), parameter :: PiBy2 = 1.57079632679489662_rb
    real(rb) :: cos_phi2, aux

    phi(2) = asin( U(3) )

    cos_phi2 = cos(phi(2))

    if (cos_phi2 == 0.0_rb) then

      phi(1) = 0.0_rb

    else

      aux = 1.0_rb/cos_phi2

      if ( aux*U(1) >= 0.0_rb ) then  !  1st or 4th Quadrant

        phi(1) = asin( aux*U(2) )

      else ! 2nd or 3rd Quadrant

        phi(1) = acos( aux*U(2) ) + PiBy2

      end if

    end if

  end function CartToSpheric

!===============================================================================

  function SphericToCart( phi ) result ( U )

    real(rb), intent(in) :: phi(dim-1)
    real(rb) :: U(dim)

    ! U(1) = cos(phi(1))*cos(phi(2))
    ! U(2) = sin(phi(1))*cos(phi(2))
    ! U(3) = sin(phi(2))

    real(rb) :: cos_phi2

    cos_phi2 = cos(phi(2))
    U(1) = cos(phi(1))*cos_phi2
    U(2) = sin(phi(1))*cos_phi2
    U(3) = sin(phi(2))

  end function SphericToCart

!===============================================================================

  function Rotate( V, Theta, U ) result ( R )

    real(rb), intent(in) :: V(dim), Theta, U(dim)
    real(rb) :: R(dim)

    !    Executes a rotation of a vector V through and angle Theta about the
    ! axis given by the unitary vector U.
    !
    ! Reference: Rutherford Aris, Vectors, Tensors, and the Basic Equations of
    !            Fluid Mechanics, Dover Pub. Inc., New York, 1962, p. 263

    real(rb) :: SinTheta, CosTheta, Aux

    SinTheta = sin(Theta)
    CosTheta = cos(Theta)
    Aux = sum(U*V)*(1.0_rb - CosTheta)

    R = V*CosTheta + U*Aux + (/U(2)*V(3)-V(2)*U(3), &
                               U(3)*V(1)-V(3)*U(1), &
                               U(1)*V(2)-V(1)*U(2)/)*SinTheta

  end function Rotate

!===============================================================================

end module MDimension