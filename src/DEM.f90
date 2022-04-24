program DEM

! To do:
! - Subroutine SaveState: Save in a different file and then rename it to StateFile

use MDEM
use MUtil
use MTime
use MFluid

implicit none

integer, parameter :: Outlet = 25

character(100) :: BaseName
character(104) :: ConfigFile, &
                  StateFile,  &
                  StateFile2, &
                  ReportFile, &
                  OutletFile

character(4), parameter :: ConfigExtension = '.cfg',  &
                           StateExtension  = '.stt',  &
                           StateExtension2 = '.st2',  &
                           ReportExtension = '.rpt',  &
                           OutletExtension = '.csv'

real(rb) :: MinStepsPerCollision, &
            MinCollisionTime,     &
            StateInterval,        &
            SnapshotInterval,     &
            ShowInterval

integer(ib) :: NumberOfSteps, &
               InitialStep,   &
               ISaveState,    &
               ISaveSnapshot, &
               IShow

logical, allocatable :: Smooth(:)

logical(lb) :: ReadInitialConfig
character(104) :: InitialConfigFile

logical(lb) :: Pouring
integer(1) :: PouringShape
real(rb) :: PouringInterval
integer(ib) :: MaxAttempts, IPouring
integer(ib), allocatable :: TargetNumber(:)
real(rb), allocatable :: PouringSize(:)

logical(lb) :: Resume, Finished
type(RealTime) :: Startup, Duration
integer(ib) :: NumberOfInterruptions
character(28) :: StartTime, FinishTime
character(28), allocatable :: Interruption(:), Retaking(:)

integer(ib) :: Step, Conf, Nconf, NDigits, MyNP, iTable = 0_ib

!============================== Executable  code ===============================

write(*,'("+",50("-"),"+")')
write(*,'("|",11(" "),"Dynamics of Granular Systems",11(" "),"|")')
write(*,'("|",13(" "),"Discrete Element Method",14(" "),"|")')
write(*,'("+",50("-"),"+")')

call ReadSpecifications( 5 )
call CalculateProperties( NC, NW )
write(*,'("Total number of time steps: ",I10)') NumberOfSteps
call ScaleVariables
call InitializeAuxiliaries

if (Resume) then
  if (Exists(StateFile)) then
    call ReadState( StateFile )
  else if (Exists(StateFile2)) then
    call ReadState( StateFile2 )
  else
    call InitializeState
  end if
else
  call InitializeState
end if

if (.not.Finished) then

  if (InitialStep == 1) then
    do Step = 1, MethodOrder - 1
      call MoveParticlesEuler
    end do
    InitialStep = MethodOrder
  end if

  do Step = InitialStep, NumberOfSteps
    call PouringAttempt( Step )
    call MoveParticles( Step )
    call UpdateParticleList
    call PerformPeriodicTasks( Step )
  end do

  call FinalizeDEM

  call SaveConfiguration( ConfigFile, Complete )

  Duration = Startup + ElapsedTime()

end if

call GenerateReport

contains

!============================= Internal Procedures =============================

  subroutine PouringAttempt( Step )

    integer(ib), intent(in) :: Step

    if ( Pouring .and. (mod(Step,IPouring) == 0) ) then
      call PourParticle( TargetNumber, PouringShape, PouringSize, MaxAttempts )
    end if

  end subroutine PouringAttempt

!===============================================================================

  subroutine PerformPeriodicTasks( Step )

    integer(ib), intent(in) :: Step

    integer(ib)  :: i
    real(rb)     :: V(NC)
    character(3) :: CC

    Finished = Step == NumberOfSteps

    if ( mod(Step,ISaveSnapshot) == 0 ) then
      Conf = Conf + 1_ib
      call SaveConfiguration( FileSeries(ConfigFile,Conf,NDigits), Compact )
    end if

    write(CC,'(I3)') NC

    if ( (mod(Step,IShow) == 0) .or. Finished ) then

      call ShowTable( Step )

      open( unit=Outlet, file=OutletFile, status='old', position='append' )
      do i = 1, NC
        V(i) = LengthScale/TimeScale*sum(Particle%Velocity(3),Particle%Index == i)/count(Particle%Index == i)
      end do
      write(Outlet,'(ES15.7,'//CC//'(",",ES15.7),",",I10)') TimeScale*Time, V, Ncol
      close(Outlet)

    end if

!    if (NP /= MyNP) then
!      open( unit=Outlet, file=OutletFile, status='old', position='append' )
!      write(Outlet,*) TimeScale*Time, Component%Number, NP
!      close(Outlet)
!      MyNP = NP
!    end if

    if ( (mod(Step,ISaveState) == 0) .or. Finished ) then
      if ( (mod(Step,2*ISaveState) == 0) .or. Finished ) then
        call SaveState( StateFile2, Step )
      else
        call SaveState( StateFile, Step )
      end if
    end if

  end subroutine PerformPeriodicTasks

!===============================================================================

  subroutine ShowTable( Step )

    integer(ib), intent(in) :: Step

    integer(ib), parameter :: RestartTable = 21
    real(rb) :: Percentage, TimeSI
    character(14) :: Elapsed

    Percentage = 100.0_rb*real(Step,rb)/real(NumberOfSteps,rb)
    TimeSI = TimeScale*Time
    Elapsed = TimeToStr(Startup + ElapsedTime())
    Elapsed = adjustr(Elapsed)

    if ( mod(iTable,RestartTable) == 0 ) then
      write(*,'("+",15("-"),"+",15("-"),"+",15("-"),"+",17("-"),"+")')
      write(*,'("|     Time      |   Particles   '// &
                '|   Concluded   |    Real Time    |")')
      write(*,'("+",15("-"),"+",15("-"),"+",15("-"),"+",17("-"),"+")')
    end if
    write(*,'("|",F12.5," s |",I14," |",F12.2," % |",A14," h |")') &
          TimeSI, NP, Percentage, Elapsed
    if ( Finished ) then
      write(*,'("+",15("-"),"+",15("-"),"+",15("-"),"+",17("-"),"+")')
    end if
    iTable = iTable + 1_ib

  end subroutine ShowTable

!===============================================================================

  subroutine GenerateReport

    integer, parameter :: Out = 18

    integer(ib) :: i

    write(*,'("")')
    write(*,'(A)') "Simulation time: "//trim(TimeToStr(Duration))
    write(*,'("")')
    write(*,'("Saving report in <",A,">")') trim(ReportFile)

    open( unit=Out, file=ReportFile, status = 'replace' )

    write(Out,'("Simulation started on ",A28)') StartTime
    write(Out,'("")')
    do i = 1_ib, NumberOfInterruptions
      write(Out,'("Simulation interrupted on ",A28)') Interruption(i)
      write(Out,'("Simulation restarted on ",A28)') Retaking(i)
      write(Out,'("")')
    end do
    write(Out,'("Simulation concluded on ",A28)') FinishTime
    write(Out,'("")')
    write(Out,'(A)') "Time spent in the simulation: "// &
                     trim(TimeToStr(Duration))//" h"

    write(Out,'("")')
    write(Out,'("Snapshots:")')
    write(Out,'("First File: <",A,">")') FileSeries(ConfigFile,0_ib,NDigits)
    write(Out,'("Last File:  <",A,">")') FileSeries(ConfigFile,NConf,NDigits)

    write(Out,'("")')
    write(Out,'("Number of integration steps: ",I10)') NumberOfSteps

    close(Out)

  end subroutine GenerateReport

!===============================================================================

  subroutine CalculateProperties( NC, NW )

    integer(ib), intent(in) :: NC, NW

    real(rb), parameter :: Pi = 3.141592653589793238_rb

    integer(ib) :: i, j
    real(rb)    :: iMass, iStiffness, iShearStiffness,  &
                   Meff(NC,NC),  Meffw(NC,NW),          &
                   LogE(NC,NC), LogEw(NC,NW),           &
                   Gamma(NC,NC), Gammaw(NC,NW)

    write(*,'("Calculating properties...")')

    ! Calculate the mass of all the components:
    Component%Mass = 4.0_rb*Pi/3.0_rb*Component%Radius**3*Component%Density

    ! Calculate the moment of inertia of all the components:
    Component%Inertia = 2.0_rb/5.0_rb*Component%Mass*Component%Radius**2

    ! Calculate the effective mass and the spring contants of the
    ! component-component pairs:
    allocate( Kn(NC,NC), Kt(NC,NC) )
    do i = 1_ib, NC
      iMass = Component(i)%Mass
      iStiffness = Component(i)%Stiffness
      iShearStiffness = Component(i)%ShearStiffness
      Meff(i,i) = 0.5_rb*iMass
      Kn(i,i) = 0.5_rb*iStiffness
      Kt(i,i) = 0.5_rb*iShearStiffness
      do j = i + 1_ib, NC
        Meff(i,j) = iMass*Component(j)%Mass/(iMass + Component(j)%Mass)
        Meff(j,i) = Meff(i,j)
        Kn(i,j) = iStiffness*Component(j)%Stiffness /  &
                  (iStiffness + Component(j)%Stiffness)
        Kn(j,i) = Kn(i,j)
        Kt(i,j) = iShearStiffness*Component(j)%ShearStiffness /  &
                  (iShearStiffness + Component(j)%ShearStiffness)
        Kt(j,i) = Kt(i,j)
      end do
    end do

    ! Calculate the damping coefficients of the component-component pairs:
    allocate( Yn(NC,NC) )
    LogE = log(En)
    Gamma = -2.0_rb*logE/sqrt(logE**2 + Pi**2)
    Yn    = Gamma*sqrt(Meff*Kn)

    ! Calculate the tangential damping coefficients of the
    ! component-component pairs:
    allocate( Yt(NC,NC) )
    LogE = log(Et)
    Gamma = -2.0_rb*logE/sqrt(logE**2 + Pi**2)
    Yt = Gamma*sqrt(Meff*Kt)

    ! Calculate the effective mass and the spring contants of the
    ! component-wall pairs:
    allocate( Knw(NC,NW), Ktw(NC,NW) )
    do i = 1_ib, NC
      iMass = Component(i)%Mass
      iStiffness = Component(i)%Stiffness
      iShearStiffness = Component(i)%ShearStiffness
      do j = 1_ib, NW
        Meffw(i,j) = iMass
        Knw(i,j) = iStiffness*Wall(j)%Stiffness /  &
                   (iStiffness + Wall(j)%Stiffness)
        Ktw(i,j) = iShearStiffness*Wall(j)%ShearStiffness /  &
                   (iShearStiffness + Wall(j)%ShearStiffness)
      end do
    end do

    ! Calculate the damping coefficients of the component-wall pairs:
    allocate( Ynw(NC,NW) )
    LogEw = log(Enw)
    Gammaw = -2.0_rb*logEw/sqrt(logEw**2 + Pi**2)
    Ynw    = Gammaw*sqrt(Meffw*Knw)

    ! Calculate the tangential damping coefficients of the
    ! component-wall pairs:
    allocate( Ytw(NC,NW) )
    LogEw = log(Etw)
    Gammaw = -2.0_rb*logEw/sqrt(logEw**2 + Pi**2)
    Ytw    = Gammaw*sqrt(Meffw*Ktw)

    ! Calculate the minimum characteristic collision time:
    MinCollisionTime = TwoPi*min(minval(Meff/sqrt(4.0_rb*Meff*Kn - Yn**2)),  &
                                 minval(Meffw/sqrt(4.0_rb*Meffw*Knw - Ynw**2)))

    ! Calculate the time Step and the number of steps:
    TimeStep = MinCollisionTime/MinStepsPerCollision
    NumberOfSteps = nint(TotalTime/TimeStep,ib)
    TimeStep = TotalTime/real(NumberOfSteps,rb)

    ! Calculate the number of steps for saving the state of the
    ! simulation, for saving snapshots, and for showing intermediary
    ! results:
    ISaveState = max(1_ib,nint(StateInterval/TimeStep,ib))
    ISaveSnapshot = max(1_ib,nint(SnapshotInterval/TimeStep,ib))
    IShow = max(1_ib,nint(ShowInterval/TimeStep,ib))

    ! Calculate the number of steps for attempting particle pouring:
    IPouring = max(1_ib,nint(PouringInterval/TimeStep,ib))

    ! Calculate the total number of configurations to be saved and the
    ! number of digits necessary to enumerate the files:
    Nconf = NumberOfSteps/ISaveSnapshot
    NDigits = int(log10(real(Nconf,rb)),ib) + 1_ib

    ! Deallocate variables that will no longer be used:
    deallocate( En, Enw )

  end subroutine CalculateProperties

!===============================================================================

  subroutine ReadSpecifications( In )

    integer, intent(in) :: In

    !   This procedure is used to read the specifications of a simulation from
    ! an open external file whose unit is passed through the argument "In".
    ! All the physical properties must be expressed in the MKS system.

    character(37), parameter :: Error = 'Error reading the specification file:'
    integer(ib) :: i
    real(rb) :: MaxRadius, Diameter, Frequency(dim), Phi(dim)
    character(50) :: String
    type(TComponent) :: iComp
    type(TWall) :: iWall

    write(*,'("Reading specifications...")')

    ! Read the dimensions of the simulation box:
    read(In,*)
    read(In,*) Lbox
    if (any(Lbox < 0.0_rb)) then
      write(*,'(A37)') Error
      write(*,'("At least one of the specified box dimensions is negative.")')
      stop
    end if

    ! Read the number of components:
    read(In,*)
    read(In,*) NC
    if (NC <= 0_ib) then
      write(*,'(A37)') Error
      write(*,'("The number of components cannot be zero or negative.")')
      stop
    end if

    ! Read the properties of the components (Number, Radius, Density,
    ! and Stiffness):
    read(In,*)
    allocate( Component(NC) )
    do i = 1_ib, NC
      read(In,*) iComp%Number, Diameter,  &
                 iComp%Density, iComp%Stiffness, iComp%ShearStiffness
      iComp%Radius = 0.5_rb*Diameter
      if ( iComp%Number < 0_ib ) then
        write(*,'(A37)') Error
        write(*,'("The initial number of particles for '//  &
                'component ",I3," cannot be negative.")') i
        stop
      end if
      if ( (iComp%Radius    <= 0.0_rb) .or. &
           (iComp%Density   <= 0.0_rb) .or. &
           (iComp%StiffNess <= 0.0_rb) ) then
        write(*,'(A37)') Error
        write(*,'("The radius, density, and stiffness of component "'//  &
                ',I3," cannot be negative.")') i
        stop
      end if
      if ( iComp%ShearStiffNess < 0.0_rb ) then
        write(*,'(A37)') Error
        write(*,'("The shear stiffness of component "'//  &
                ',I3," must be positive.")') i
        stop
      end if
      Component(i) = iComp
    end do

    ! Read the control of initial configuration:
    read(In,*)
    read(In,*) String
    ReadInitialConfig = UpperCase(trim(adjustl(String))) /= 'NO'
    if (ReadInitialConfig) then
      if ( Exists(trim(String)) ) then
        InitialConfigFile = trim(String)
      else
        write(*,'(A37)') Error
        write(*,'("The specified file with the initial configuration, <'// &
                trim(String)//'>, does not exist.")')
        stop
      end if
    end if

    ! Read the pouring control:
    read(In,*)
    read(In,*) String
    Pouring = UpperCase(adjustl(String)) /= 'NO'
    if (Pouring) then
      if ( Exists(trim(String)) ) then
        call ReadPouringSpecifications( trim(String) )
      else
        write(*,'(A37)') Error
        write(*,'("The specified file with the pouring specifications, <'// &
                trim(String)//'>, does not exist.")')
        stop
      end if
    end if

    ! Read the number of walls:
    read(In,*)
    read(In,*) NW

    if (NW < 0_ib) then
      write(*,'(A37)') Error
      write(*,'("The specified number of walls must be positive or zero.")')
      stop
    else if (NW >= 0_ib) then
      allocate( Wall(NW), Smooth(NW) )
    end if

    ! Read the properties of the walls (Shape, L1, L2, L3, Stiffness,
    ! Position, U, and W):
    read(In,*)
    MaxRadius = maxval(Component%Radius)
    do i = 1_ib, NW

      read(In,*) iWall%Shape, iWall%L1, iWall%L2, iWall%L3, iWall%Stiffness, &
                 iWall%ShearStiffness, iWall%Position, iWall%U, iWall%W

      iWall%L1 = 0.5_rb*iWall%L1
      iWall%L2 = 0.5_rb*iWall%L2
      iWall%L3 = 0.5_rb*iWall%L3

      if ( (iWall%Shape < 1_ib).or.(iWall%Shape > NS) ) then
        write(*,'(A37)') Error
        write(*,'("The shape specified for wall ",I3," is invalid.")') i
        stop
      else if (iWall%Stiffness <= 0.0_rb) then
        write(*,'(A37)') Error
        write(*,'("The stiffness specified for wall ",I3," is'//  &
                ' must be positive.")') i
        stop
      else if ( (iWall%L1 < 0.0_rb) .or. &
                (iWall%L2 < 0.0_rb) .or. &
                (iWall%L3 < 0.0_rb)      ) then
        write(*,'(A37)') Error
        write(*,'("At least one of the dimensions specified for '// &
                'wall ",I3," is negative.")') i
        stop
      else if ( sum((iWall%U .x. iWall%W)**2) == 0.0_rb ) then
        write(*,'(A37)') Error
        write(*,'("The specified direction vectors for wall ",I3,'// &
                '" must be linearly independent.")') i
        stop
      end if

      call GramSchmidt( iWall%U, iWall%W)
      iWall%N = iWall%U .x. iWall%W

      Wall(i) = iWall

    end do

    ! Read the normal coefficients of restitution of the component-component pairs:
    read(In,*)
    allocate( En(NC,NC) )
    do i = 1_ib, NC
      read(In,*) En(i,:)
    end do
    if ( any(En <= 0.0_rb).or.any(En > 1.0_rb) ) then
      write(*,'(A37)') Error
      write(*,'("All the component-component normal coefficients of restitution'// &
              ' must lie in the interval (0,1].")')
      stop
    end if

    ! Read the tangential coefficients of restitution of the component-component pairs:
    read(In,*)
    allocate( Et(NC,NC) )
    do i = 1_ib, NC
      read(In,*) Et(i,:)
    end do
    if ( any(Et <= 0.0_rb).or.any(Et > 1.0_rb) ) then
      write(*,'(A37)') Error
      write(*,'("All the component-component tangential coefficients of restitution'// &
              ' must lie in the interval (0,1].")')
      stop
    end if

    ! Read the coefficients of friction of the component-component pairs:
    read(In,*)
    allocate( Mu(NC,NC) )
    do i = 1_ib, NC
      read(In,*) Mu(i,:)
    end do
    if ( any(Mu <= 0.0_rb) ) then
      write(*,'("All the component-component coefficients'//  &
              ' of friction must be positive.")')
      stop
    end if

    ! Read the normal coefficients of restitution of the component-wall pairs:
    read(In,*)
    allocate( Enw(NC,NW) )
    do i = 1_ib, NC
      read(In,*) Enw(i,:)
    end do
    if ( any(Enw <= 0.0_rb).or.any(Enw > 1.0_rb) ) then
      write(*,'("All the component-wall normal coefficients of restitution'//  &
              ' must lie in the interval (0,1].")')
      stop
    end if

    ! Read the tangential coefficients of restitution of the component-wall pairs:
    read(In,*)
    allocate( Etw(NC,NW) )
    do i = 1_ib, NC
      read(In,*) Etw(i,:)
    end do
    if ( any(Etw <= 0.0_rb).or.any(Etw > 1.0_rb) ) then
      write(*,'("All the component-wall tangential coefficients of restitution'//  &
              ' must lie in the interval (0,1].")')
      stop
    end if

    ! Read the coefficients of friction of the component-wall pairs:
    read(In,*)
    allocate( Muw(NC,NW) )
    do i = 1_ib, NC
      read(In,*) Muw(i,:)
    end do
    if ( any(Muw <= 0.0_rb) ) then
      write(*,'(A37)') Error
      write(*,'("All the component-wall coefficients'//  &
              ' of friction must be positive.")')
      stop
    end if

    ! Read the periodic boundary conditions:
    read(In,*)
    read(In,*) PBC

    ! Read the gravity acceleration vector:
    read(In,*)
    read(In,*) Gravity

    ! Read the fluid density:
    read(In,*)
    read(In,*) FluidDensity

    ! Read a seed for random number generation:
    read(In,*)
    read(In,*) Dummy
    if (Dummy >= 0) then
      write(*,'(A37)') Error
      write(*,'("The seed for random number generation must be negative.")')
    end if

    ! Read the base for file names:
    read(In,*)
    read(In,*) BaseName
    i = scan(BaseName,'.')
    if ( (i /= 0).and.(i /= len_trim(BaseName)) ) then
      write(*,'(A37)') Error
      write(*,'("A base for file names should not contain an extension.")')
      stop
    end if

    ! Define the names of the configuration file and state file:
    ConfigFile  = trim(BaseName)//ConfigExtension
    StateFile   = trim(BaseName)//StateExtension
    StateFile2  = trim(BaseName)//StateExtension2
    ReportFile  = trim(BaseName)//ReportExtension
    OutletFile  = trim(BaseName)//OutletExtension

    ! Read the total time to be simulated:
    read(In,*)
    read(In,*) TotalTime
    if (TotalTime <= 0.0_rb) then
      write(*,'(A37)') Error
      write(*,'("The specified time for the simulation must be positive.")')
      stop
    end if

    ! Read the minimum number of time steps per collision:
    read(In,*)
    read(In,*) MinStepsPerCollision
    if (MinStepsPerCollision <= 0.0_rb) then
      write(*,'(A37)') Error
      write(*,'("The minimum number of steps per collision'//  &
              ' must be positive.")')
      stop
    end if

    ! Read the order of the predictor-corrector integration method:
    read(In,*)
    read(In,*) MethodOrder
    if ( (MethodOrder < 1).or.(MethodOrder > 4) ) then
      write(*,'(A37)') Error
      write(*,'("The specified order for the integration method must'//  &
              ' be 1, 2, 3, or 4.")')
      stop
    end if

    ! Read the number of intervals in each initial euler time step:
    read(In,*)
    read(In,*) NumberOfEulerIntervals

    if ( NumberOfEulerIntervals < 1 ) then
      write(*,'(A37)') Error
      write(*,'("The specified number of intervals in each initial Euler'//  &
              ' time step must be greater than zero.")')
      stop
    end if

    ! Read the time interval for saving the state of the simulation:
    read(In,*)
    read(In,*) StateInterval
    if ( (StateInterval <= 0.0_rb) .or. (StateInterval > TotalTime) ) then
      write(*,'(A37)') Error
      write(*,'("The specified interval for saving the state of the '// &
              'simulation is invalid.")')
      stop
    end if

    ! Read the time interval for saving snapshots of the simulation:
    read(In,*)
    read(In,*) SnapshotInterval
    if ((SnapshotInterval <= 0.0_rb).or.(SnapshotInterval > TotalTime)) then
      write(*,'(A37)') Error
      write(*,'("The specified interval for saving snapshots of the '// &
              'simulation is invalid.")')
      stop
    end if

    ! Read the time interval for showing intermediary results:
    read(In,*)
    read(In,*) ShowInterval
    if ((ShowInterval <= 0.0_rb).or.(ShowInterval > TotalTime)) then
      write(*,'(A37)') Error
      write(*,'("The specified interval for reporting the progress of the '// &
              'simulation is invalid.")')
      stop
    end if

    ! Read the resuming control:
    read(In,*)
    read(In,*) String
    String = UpperCase(adjustl(String))
    if ( (trim(String) /= 'YES').and.(trim(String) /= 'NO') ) then
      write(*,'(A37)') Error
      write(*,'("The answer for the question about simulation remuming '// &
              'must be ""Yes"" or ""No"".")')
      stop
    end if
    Resume = trim(String) == 'YES'

    ! Read the wall movement control:
    ! Read the resuming control:
    read(In,*)
    read(In,*) String
    String = UpperCase(adjustl(String))
    if ( (trim(String) /= 'YES').and.(trim(String) /= 'NO') ) then
      write(*,'(A37)') Error
      write(*,'("The answer for the question about wall movement '// &
              'must be ""Yes"" or ""No"".")')
      stop
    end if
    WallsInMovement = trim(String) == 'YES'

    read(In,*)
    if (WallsInMovement) then
      allocate( WM(NW) )
      do i = 1_ib, NW
        read(In,*) WM(i)%Amplitude, Frequency, Phi
        if (any(Frequency < 0.0_rb)) then
          write(*,'(A37)') Error
          write(*,'("The frequency vector for the movement of wall ",I3,'// &
                  '" cannot have any negative element.")') i
          stop
        end if
        WM(i)%TwoPiFrequency = TwoPi*Frequency
        WM(i)%PhaseAngle = TwoPi*Phi
        WM(i)%VelocityAmplitude = -WM(i)%TwoPiFrequency*WM(i)%Amplitude
        WM(i)%ConstantFactor = Wall(i)%Position - WM(i)%Amplitude*cos(WM(i)%PhaseAngle)
      end do
    end if

    ! Read the fluid specifications:
    call ReadFluidSpecifications( In, Error )

  end subroutine ReadSpecifications

!===============================================================================

  subroutine ReadPouringSpecifications( FileName )

    character(*), intent(in) :: FileName

    integer, parameter :: In = 73
    character(31), parameter :: Error = 'Error reading the pouring file:'

    integer(ib) :: i, Npour
    character(50) :: String

    ! Allocate the target number:
    allocate( TargetNumber(NC) )

    ! Open the pouring file:
    open( unit=In, file=FileName, status='old')

    ! Read the target number:
    read(In,*)
    do i = 1_ib, NC
      read(In,*) TargetNumber(i)
      if ( TargetNumber(i) < 0_ib ) then
        write(*,'(A31)') Error
        write(*,'("The desired number of particles for '//  &
                'component ",I3," cannot be negative.")') i
        stop
      end if
    end do

    ! Read the shape of the pouring origin surface:
    read(In,*)
    read(In,*) String
    String = UpperCase(adjustl(String))
    if ( trim(String) == 'RECTANGLE' ) then
      PouringShape = Rect
      Npour = 2
    else if ( trim(String) == 'CIRCLE' ) then
      PouringShape = Circle
      Npour = 1
    else
      write(*,'(A31)') Error
      write(*,'("The shape specified for the pouring origin surface'//  &
              ' is invalid.")')
      stop
    end if

    ! Read the size of the pouring origin surface:
    allocate( PouringSize(NPour) )
    read(In,*)
    read(In,*) PouringSize(1:NPour)

    ! Read the maximum number of attempts in each time step:
    read(In,*)
    read(In,*) MaxAttempts

    ! Read the interval for pouring attempts:
    read(In,*)
    read(In,*) PouringInterval

    close(In)

  end subroutine ReadPouringSpecifications

!===============================================================================

  subroutine ScaleVariables

    integer(ib) :: i

    write(*,'("Scaling variables...")')

    ! Determine the scale factors:
    LengthScale = 2.0_rb*minval(Component%Radius) ! The minimum diameter
    MassScale   = minval(Component%Mass)          ! The minimum mass
    TimeScale   = MinCollisionTime                ! The minimum collision time

    print*, 'Length Scale = ', LengthScale, ' m.'
    print*, 'Time Scale = ', TimeScale, ' s.'

    ! Calculate the time Step in reduced unit:

    ! Scale all the dimensional variables:
    Lbox = Lbox/LengthScale
    InvLbox = 1.0_rb/Lbox
    Component%Radius = Component%Radius/LengthScale
    Component%Mass = Component%Mass/MassScale
    Component%Density = Component%Density*LengthScale**3/MassScale
    Component%Inertia = Component%Inertia/(LengthScale**2*MassScale)
    Component%Stiffness = Component%Stiffness*TimeScale**2/MassScale
    Component%ShearStiffness = Component%ShearStiffness*TimeScale**2/MassScale
    Wall%L1 = Wall%L1/LengthScale
    Wall%L2 = Wall%L2/LengthScale
    Wall%L3 = Wall%L3/LengthScale
    Wall%Stiffness = Wall%Stiffness*TimeScale**2/MassScale
    Wall%ShearStiffness = Wall%ShearStiffness*TimeScale**2/MassScale
    forall(i=1:NW)
      Wall(i)%Position = Wall(i)%Position/LengthScale
      Wall(i)%Velocity = Wall(i)%Velocity*TimeScale/LengthScale
    end forall
    if (WallsInMovement) then
      forall(i=1:NW)
        WM(i)%ConstantFactor = WM(i)%ConstantFactor/LengthScale
        WM(i)%Amplitude = WM(i)%Amplitude/LengthScale
        WM(i)%VelocityAmplitude = WM(i)%VelocityAmplitude*TimeScale/LengthScale
        WM(i)%TwoPiFrequency = WM(i)%TwoPiFrequency*TimeScale
      end forall
    end if
    Kn  = Kn*TimeScale**2/MassScale
    Yn  = Yn*TimeScale/MassScale
    Knw = Knw*TimeScale**2/MassScale
    Ynw = Ynw*TimeScale/MassScale
    Kt  = Kt*TimeScale**2/MassScale
    Yt  = Yt*TimeScale/MassScale
    Ktw = Ktw*TimeScale**2/MassScale
    Ytw = Ytw*TimeScale/MassScale
    Gravity = Gravity*TimeScale**2/LengthScale
    FluidDensity = FluidDensity*LengthScale**3/MassScale
    TotalTime = TotalTime/TimeScale
    TimeStep = TimeStep/TimeScale
    PouringSize = PouringSize/LengthScale
    PouringInterval = PouringInterval/TimeScale

    if (FluidFlow) call ScaleFluidProperties

  end subroutine ScaleVariables

!===============================================================================

  subroutine InitializeState

    write(*,'("* New simulation *")')

    ! Assert that the simulation is not finished:
    Finished = .false._lb

    ! Set up the startup date and time of the simulation:
    StartTime = StrDateAndTime()

    ! Initialize number of interruptions:
    NumberOfInterruptions = 0

    ! Initialize start time as zero:
    Startup = RealTime(0,0,0,0)

    ! Configure wall variables:
    call ConfigureWalls

    ! Read or create an initial configuration:
    if (ReadInitialConfig) then

      write(*,'("Reading initial configuration...")')
      call ReadParticlePositions( InitialConfigFile )

    else

      if (any(Component%Number > 0_ib)) then
        write(*,'("Creating initial configuration...")')
        if (.not.SimpleLatticeConfiguration(Dummy)) then
          write(*,'("It was not possible to construct an initial configuration.")')
          stop
        end if
      else
        allocate( Particle(0) )
      end if

    end if

    ! Initialize the time:
    Time = 0.0_rb

    ! Initialize maximum number of particles:
    MaxNP = NP

    ! Distribute particles throughout the cells:
    call InitializeCellDistribution

    ! Initialize data related to numerical integration:
    call InitializeIntegrationVariables

    ! Initialize the data structure for storage of contacting pairs:
    call InitializeContacts( NP )

    ! Set the initial Step as 1:
    InitialStep = 1_ib

    ! Save the initial configuration:
    Conf = 0_ib
    call SaveConfiguration( FileSeries(ConfigFile,Conf,NDigits), Compact )

    ! Open the outlet file:
    open( unit = Outlet, file = OutletFile, status = 'replace' )
    write(Outlet,'("'//trim(OutletTitles(NC))//'")')
    close(Outlet)
    MyNP = NP

    ! Start the time counter:
    call StartTimeCounter

    ! Start showing the table of simulation progress:
    call ShowTable( 0_ib )

  end subroutine InitializeState

!===============================================================================

  subroutine SaveState( StateFile, Step )

    character(*), intent(in) :: StateFile
    integer(ib), intent(in) :: Step

    integer :: Out, i

    ! Open a binary file for writing:
    Out = BinOpen( file = StateFile, action = 'write' )

    ! Save the variable that accounts for the end of the simulation:
    call BinWrite( Out, Finished )

    ! Save the startup moment of the simulation:
    call BinWrite( Out, StartTime )

    ! Save the number of interruptions:
    call BinWrite( Out, NumberOfInterruptions )

    ! Save all the interruptions and restarts:
    do i = 1, NumberOfInterruptions
      call BinWrite( Out, Interruption(i) )
      call BinWrite( Out, Retaking(i) )
    end do

    ! Save current date and time:
    FinishTime = StrDateAndTime()
    call BinWrite( Out, FinishTime )

    ! Save the simulation time:
    Duration = Startup + ElapsedTime()
    call BinWrite( Out, Duration%Hours    )
    call BinWrite( Out, Duration%Minutes  )
    call BinWrite( Out, Duration%Seconds  )
    call BinWrite( Out, Duration%MiliSecs )

    ! Save the configuration:
    call SaveConfiguration( Out, Complete )

    ! Save the particle distribution:
    call SaveCellDistribution( Out )

    ! Save the variables that concern to the integration method:
    call SaveIntegrationVariables( Out )

    ! Save the contact list:
    call SaveContactList( Out )

    ! Save variables that concert to random number generation:
    call BinWrite( Out, aam )
    call BinWrite( Out, Dummy )
    call BinWrite( Out, iix )
    call BinWrite( Out, iiy )
    call BinWrite( Out, iik )

    ! Save the index of the current configuration:
    call BinWrite( Out, Conf )

    ! Save the next initial Step:
    call BinWrite( Out, Step + 1_ib )

    ! Save the number of colisions:
    call BinWrite( Out, Ncol )

    ! Close the binary file:
    call BinClose(Out)

  end subroutine SaveState

!===============================================================================

  subroutine ReadState( StateFile )

    character(*), intent(in) :: StateFile

    integer :: In, i, N
    type(RealTime) :: Time

    ! Open a binary file for reading:
    In = BinOpen( file = StateFile, action = 'read' )

    ! Save the variable that accounts for the end of the simulation:
    call BinRead( In, Finished )

    if (Finished) then
      write(*,'("* Simulation already finished *")')
    else
      write(*,'("* Previously interrupted simulation *")')
      write(*,'("Reading the last state of the simulation...")')
    end if

    ! Read the startup moment of the simulation:
    call BinRead( In, StartTime )

    ! Read the number of interruptions:
    call BinRead( In, NumberOfInterruptions )

    ! Update the number of interruptions and allocate the variables:
    N = NumberOfInterruptions

    if (.not.Finished) NumberOfInterruptions = NumberOfInterruptions + 1_ib
    allocate( Interruption(NumberOfInterruptions),  &
              Retaking(NumberOfInterruptions) )

    ! Read all the previous interruptions and restarts:
    do i = 1_ib, N
      call BinRead( In, Interruption(i) )
      call BinRead( In, Retaking(i) )
    end do

    ! Read date and time at the last interruption or at the simulation's end:
    if (Finished) then
      call BinRead( In, FinishTime )
    else
      call BinRead( In, Interruption(NumberOfInterruptions) )
      Retaking(NumberOfInterruptions) = StrDateAndTime()
    end if

    ! Read the simulation time:
    call BinRead( In, Time%Hours    )
    call BinRead( In, Time%Minutes  )
    call BinRead( In, Time%Seconds  )
    call BinRead( In, Time%MiliSecs )

    if (Finished) then
      Duration = Time
    else
      Startup = Time
    end if

    ! Read the configuration:
    call ReadConfiguration( In )

    ! Initialize maximum number of particles:
    MaxNP = NP

    ! Configure wall variables:
    call ConfigureWalls

    ! Read the particle distribution:
    call ReadCellDistribution( In )

    ! Read the variables that concern to the integration method:
    call ReadIntegrationVariables( In )

    ! Read the contact list:
    call ReadContactList( In )

    ! Read variables that concert to random number generation:
    call BinRead( In, aam )
    call BinRead( In, Dummy )
    call BinRead( In, iix )
    call BinRead( In, iiy )
    call BinRead( In, iik )

    ! Read the index of the current configuration:
    call BinRead( In, Conf )

    ! Read the initial Step:
    call BinRead( In, InitialStep )

    ! Read the number of collisions:
    call BinRead( In, Ncol )

    ! Close the binary file:
    call BinClose(In)

    if (.not.Finished) then
      write(*,'("Resuming it.")')
    end if

    ! Open the outlet file:
    if (.not.Exists(OutletFile)) then
      write(*,*) 'Warning: The outlet file was not found. Creating a new one.'
      open( unit = Outlet, file = OutletFile, status = 'replace' )
      write(Outlet,'("'//trim(OutletTitles(NC))//'")')
      !write(Outlet,'("Time(s), <Vt1>(m/s), <Vt2>(m/s), Ncol")')
      close(Outlet)
    end if
    MyNP = NP

    ! Start the time counter:
    call StartTimeCounter

    ! Start showing the table of simulation progress:
    call ShowTable( InitialStep )

  end subroutine ReadState

!===============================================================================

  function OutletTitles( NC ) result ( Titles )

    integer(ib), intent(in) :: NC
    character(256) :: Titles

    integer(ib) :: Index
    character(5) :: StrIndex

    Titles = "Time(s)"
    do Index = 1_ib, NC
      write(StrIndex,'(I5)') Index
      Titles = trim(Titles)//", V"//trim(adjustl(StrIndex))
    end do
    Titles = trim(Titles)//", Ncol"

  end function OutletTitles

!===============================================================================

end program DEM
