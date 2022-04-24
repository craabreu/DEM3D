module MTime

! Module variables for the time counter:
type RealTime
  integer(4) :: Hours
  integer(1) :: Minutes,  &
                Seconds
  integer(2) :: MiliSecs
end type RealTime

interface operator (+)
  module procedure AddTimes
end interface operator (+)

interface operator (-)
  module procedure SubtractTimes
end interface operator (-)

interface assignment(=)
  module procedure AssignTime
end interface assignment(=)

type(RealTime), private :: InitialTime

contains

!===============================================================================

  subroutine StartTimeCounter

    InitialTime = GetRelativeTime()

  end subroutine StartTimeCounter

!===============================================================================

  function ElapsedTime() result ( Time )

    type(RealTime) :: Time

     Time = GetRelativeTime() - InitialTime

  end function ElapsedTime

!===============================================================================

  subroutine AssignTime( Time1, Time2 )

    type(RealTime), intent(out) :: Time1
    type(RealTime), intent(in)  :: Time2

    Time1%Hours    = Time2%Hours
    Time1%Minutes  = Time2%Minutes
    Time1%Seconds  = Time2%Seconds
    Time1%MiliSecs = Time2%MiliSecs

  end subroutine AssignTime

!===============================================================================

  function AddTimes( Time1, Time2 ) result( ResultTime )

    type(RealTime), intent(in) :: Time1, Time2
    type(RealTime) :: ResultTime

    integer(4) :: H, M, S, MS

    H = Time1%Hours + Time2%Hours
    M = Time1%Minutes + Time2%Minutes
    S = Time1%Seconds + Time2%Seconds
    MS = Time1%MiliSecs + Time2%MiliSecs

    if (MS >= 1000_4) then
      MS = MS - 1000_4
      S = S + 1_4
    end if

    if (S >= 60_4) then
      S = S - 60_4
      M = M + 1_4
    end if

    if (M >= 60_4) then
      M = M - 60_4
      H = H + 1_4
    end if

    ResultTime = RealTime( H, M, S, MS )

  end function AddTimes

!===============================================================================

  function SubtractTimes( Time1, Time2 ) result( ResultTime )

    type(RealTime), intent(in) :: Time1, Time2
    type(RealTime) :: ResultTime

    integer(4) :: H, M, S, MS

    H  = Time1%Hours - Time2%Hours
    M  = Time1%Minutes  - Time2%Minutes
    S  = Time1%Seconds  - Time2%Seconds
    MS = Time1%MiliSecs - Time2%MiliSecs

    if (MS < 0_4) then
      MS = 1000_4 + MS
      S = S - 1_4
    end if

    if (S < 0_4) then
      S = 60_4 + S
      M = M - 1_4
    end if

    if (M < 0_4) then
      M = 60_4 + M
      H = H - 1_4
    end if

    ResultTime = RealTime( H, M, S, MS )

  end function SubtractTimes

!===============================================================================

  function GetRelativeTime() result ( Time )

    type(RealTime) :: Time

    integer, parameter :: DM(12) = (/  0,  31,  59,  90, 120, 151, &
                                     181, 212, 243, 273, 304, 334 /)

    integer :: DateTime(8)
    integer(4) :: Year, Month, Day, Hour, NY, NCD, NCH
    integer(1) :: Minute, Second
    integer(2) :: MiliSec

    call date_and_time( values = DateTime )

    Year    = DateTime(1)
    Month   = DateTime(2)
    Day     = DateTime(3)
    Hour    = DateTime(5)
    Minute  = DateTime(6)
    Second  = DateTime(7)
    MiliSec = DateTime(8)

    NY = Year - 1901_4                               ! Number of concluded years
    NCD = 365_4*NY + NY/4_4 + DM(Month) + Day - 1_4   ! Number of concluded days
    if (mod(Year,4_4) == 0_4) NCD = NCD + 1_4            ! Bissextile correction
    NCH = 24_4*NCD + Hour                            ! Number of concluded hours

    Time = RealTime( NCH, Minute, Second, MiliSec )

  end function GetRelativeTime

!===============================================================================

  function TimeToStr( Time ) result ( StrTime )

    type(RealTime), intent(in) :: Time
    character(20) :: StrTime

    character(2) :: CM, CS
    character(3) :: CMS, CNCH
    integer      :: NCH

    write(CM, '(I2)') Time%Minutes
    write(CS, '(I2)') Time%Seconds
    write(CMS,'(I3)') Time%MiliSecs

    select case (Time%Hours)
      case (:-1)
        NCH = int(log10(real(-Time%Hours))) + 2
      case (0)
        NCH = 1
      case (1:)
        NCH = int(log10(real(Time%Hours))) + 1
    end select
    write(CNCH,'(I3)') NCH

    if (CM(1:1)  == ' ') CM(1:1)  = '0'
    if (CS(1:1)  == ' ') CS(1:1)  = '0'
    if (CMS(1:1) == ' ') CMS(1:1) = '0'
    if (CMS(2:2) == ' ') CMS(2:2) = '0'

    write(StrTime,fmt='(I'//CNCH//',":",A2,":",A2,".",A3)') &
              Time%Hours, CM, CS, CMS

    StrTime = adjustl(StrTime)

  end function TimeToStr

!===============================================================================

  function StrDateAndTime() result ( String )

    character(28) :: String

    character(3), parameter :: StrMonth(12) = (/'Jan', 'Feb', 'Mar', 'Apr', &
                                              'May', 'Jun', 'Jul', 'Aug', &
                                              'Sep', 'Oct', 'Nov', 'Dec' /)

    character(2) :: CD, CH, CM, CS, Ord = 'th'
    character(4) :: CY
    integer      :: DateTime(8)
    integer(4)   :: Year, Month, Day, Hour
    integer(1)   :: Minute, Second

    call date_and_time( values = DateTime )

    Year    = DateTime(1)
    Month   = DateTime(2)
    Day     = DateTime(3)
    Hour    = DateTime(5)
    Minute  = DateTime(6)
    Second  = DateTime(7)

    select case (mod(Day,10))
      case (1)
        if (Day /= 11) Ord = 'st'
      case (2)
        if (Day /= 12) Ord = 'nd'
      case (3)
        if (Day /= 13) Ord = 'rd'
    end select

    write(CD,'(I2)') Day
    write(CY,'(I4)') Year
    write(CH,'(I2)') Hour
    write(CM,'(I2)') Minute
    write(CS,'(I2)') Second


    if (CM(1:1) == " ") CM(1:1) = "0"
    if (CS(1:1) == " ") CS(1:1) = "0"


    String = StrMonth(Month)//" "//CD//Ord//", "//CY// &
             " at "//CH//":"//CM//":"//CS//" h"

  end function StrDateAndTime

!===============================================================================

end module MTime