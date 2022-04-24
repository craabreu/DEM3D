module mFluid

use MTypes
use MDimension
use MConfig
use MUtil

implicit none

logical,  public  :: FluidFlow
integer,  private :: NL
real(rb), private :: FluidFlowRate,    &
                     FluidViscosity,   &
                     Delta_Z,          &
                     Area,             &
                     Vf0,              &
                     H_2,              &
                     mH_2

real(rb) :: GlobalPorosity

real(rb), parameter :: Pi_2 = Pi/2.0_rb

real(rb), allocatable, private :: Vf(:), Ef(:)

contains

!===============================================================================

  subroutine ReadFluidSpecifications( In, Error )

    integer, intent(in) :: In
    character(*), intent(in) :: Error

    character(50) :: String

    read(In,*)
    read(In,*) String
    String = UpperCase(adjustl(String))
    if ( (trim(String) /= 'YES').and.(trim(String) /= 'NO') ) then
      write(*,'(A37)') Error
      write(*,'("The answer for the question about simulation of '// &
              'fluid flow must be ""Yes"" or ""No"".")')
      stop
    end if
    FluidFlow = trim(String) == 'YES'

    if (FluidFlow) then

      read(In,*)
      read(In,*) FluidViscosity
      read(In,*)
      read(In,*) FluidFlowRate
      read(In,*)
      read(In,*) NL

    end if

  end subroutine ReadFluidSpecifications

!===============================================================================

  subroutine ScaleFluidProperties

    FluidViscosity = FluidViscosity*LengthScale*TimeScale/MassScale
    FluidFlowRate = FluidFlowRate*TimeScale/LengthScale**3

  end subroutine ScaleFluidProperties

!===============================================================================

  subroutine InitializeFluid

    integer(ib) :: i
    real(rb)    :: FluidVolume

    Area = Lbox(1)*Lbox(2)   ! <--- VALID ONLY FOR SQUARE CONTAINER
    Delta_Z = Lbox(3)/NL
    Vf0 = FluidFlowRate/Area

    H_2 = Lbox(3)/2.0_rb
    mH_2 = -H_2

    allocate( Vf(0:NL), Ef(0:NL) )

    FluidVolume = product(Lbox)
    do i = 1, NC
      FluidVolume = FluidVolume - Component(i)%Number*Component(i)%Radius**3*(4.0_rb*Pi/3.0_rb)
    end do
    GlobalPorosity = FluidVolume/product(Lbox)

  end subroutine InitializeFluid

!===============================================================================

  subroutine UpdateVelocityField( Particle )

    type(TParticle), intent(in) :: Particle(:)

    integer(ib) :: k, i, Ind
    real(rb)    :: Zk, Sum_A, Sum_B, Sum_C, Ai, Diff2, R2, Vs, Es, Dens_Ai

    do k = 0, NL

      Zk = mH_2 + real(k,rb)*Delta_Z

      Sum_A = 0.0_rb
      Sum_B = 0.0_rb
      Sum_C = 0.0_rb
      do i = 1, NP
        Ind = Particle(i)%Index
        Diff2 = (Zk - Particle(i)%Position(3))**2
        R2 = Component(Ind)%Radius**2
        if ( Diff2 < R2 ) then
          Ai = Pi*(R2 - Diff2)
          Sum_A = Sum_A + Ai
          Dens_Ai = Ai*Component(Ind)%Density
          Sum_B = Sum_B + Dens_Ai*Particle(i)%Velocity(3)
          Sum_C = Sum_C + Dens_Ai
        end if
      end do

      if (Sum_A == 0.0_rb) then
        Vf(k) = Vf0
        Ef(k) = 1.0_rb
      else
        Es = Sum_A/Area
        Ef(k) = 1.0_rb - Es
        Vs = Sum_B/Sum_C
        Vf(k) = (Vf0 - Es*Vs)/Ef(k)
      end if

    end do

  end subroutine UpdateVelocityField

!===============================================================================

  function DragForce( iPart ) result ( FD )

    type(TParticle), intent(in) :: iPart
    real(rb) :: FD(dim)

    integer(ib) :: j
    real(rb)    :: Zi, jr, Vfz, Efz, Vr(dim), Ri, Re, Vra, CD, mBeta, n, n0, x, Cd1, fi

    Zi  = iPart%Position(3)

    jr  = (Zi + H_2)/Delta_Z
    j   = int(jr)
    Vfz = Vf(j)
    Vfz = Vfz + (jr - j)*(Vf(j+1) - Vfz)
    Efz = Ef(j)
    Efz = Efz + (jr - j)*(Ef(j+1) - Efz)

    Vr = -iPart%Velocity
    Vr(3) = Vfz + Vr(3)

    Vra = sqrt(sum(Vr*Vr))

    if (Vra > 0.0_rb) then

      Ri = Component(iPart%Index)%Radius

      Re = 2.0_rb*Efz*Ri*Vra*FluidDensity/FluidViscosity

	  n = 3.6_rb
      n0 = n + 1
      do while (abs(n-n0) > 1.e-8_rb)
        n0 = n
        x = Re/(Efz**n)
        n = (4.8_rb + 0.42_rb*x**0.75_rb)/(1.0_rb + 0.175_rb*x**0.75_rb)
      end do
      Cd1 = (0.63_rb + 4.8_rb/sqrt(x))**2_ib
      CD = (0.63_rb + 4.8_rb/sqrt(Re))**2_ib
      fi = (Cd1/CD)*Efz**(2.0_rb*(1.0_rb - n))

      FD = CD*Pi_2*Ri*Ri*FluidDensity*Vra*Vr*fi

!print*, "Re = ", Re
!print*, "FluidDensity = ", FluidDensity
!print*, "mBeta = ", mBeta
!print*, "CD = ", CD
!print*, "Vra = ", Vra*LengthScale/TimeScale
!print*, GlobalPorosity
!print*, "FD = ", FD*MassScale*LengthScale/TimeScale**2
!stop

!write(*,'(ES)') FD(3)*MassScale*LengthScale/TimeScale**2

    else

      FD = 0.0_rb

    end if

  end function DragForce

!===============================================================================

  function DragForce_OLD( iPart ) result ( FD )

    type(TParticle), intent(in) :: iPart
    real(rb) :: FD(dim)

    integer(ib) :: j
    real(rb)    :: Zi, jr, Vfz, Efz, Vr(dim), Ri, Re, Vra, CD, mBeta

!    Zi  = iPart%Position(3)

!    jr  = (Zi + H_2)/Delta_Z
!    j   = int(jr)
!    Vfz = Vf(j)
!    Vfz = Vfz + (jr - j)*(Vf(j+1) - Vfz)
!    Efz = Ef(j)
!    Efz = Efz + (jr - j)*(Ef(j+1) - Efz)

    Vr = -iPart%Velocity
!    Vr(3) = Vfz + Vr(3)

    Vra = sqrt(sum(Vr*Vr))

    if (Vra > 0.0_rb) then

      Ri = Component(iPart%Index)%Radius

!Vr = (/0.0_rb,0.0_rb,-0.0027615445048414_rb/)*TimeScale/LengthScale
!Vra = 0.0027615445048414*TimeScale/LengthScale

      Re = 2.0_rb*Ri*Vra*FluidDensity/FluidViscosity

      CD = (0.63_rb + 4.8_rb/sqrt(Re))**2

      mBeta = 0.65_rb*exp(-0.5_rb*(1.5_rb - log10(Re))**2) - 4.7_rb

!      FD = Efz**mBeta*CD*Pi_2*Ri*Ri*FluidDensity*Vra*Vr
      FD = GlobalPorosity**mBeta*CD*Pi_2*Ri*Ri*FluidDensity*Vra*Vr


!print*, "Re = ", Re
!print*, "FluidDensity = ", FluidDensity
!print*, "mBeta = ", mBeta
!print*, "CD = ", CD
!print*, "Vra = ", Vra*LengthScale/TimeScale
!print*, GlobalPorosity
!print*, "FD = ", FD*MassScale*LengthScale/TimeScale**2
!stop

!write(*,'(ES)') FD(3)*MassScale*LengthScale/TimeScale**2

    else

      FD = 0.0_rb

    end if

  end function DragForce_OLD

!===============================================================================

end module mFluid
