module BinFiles

!       ## Handling of Sequential Binary Files ##
!
!           VERSION FOR COMPAQ VISUAL FORTRAN
!
!   Author: Charlles R. A. Abreu (Rio de Janeiro, Brazil)
!   Emails: rubber@eq.ufrj.br / craabreu@ig.com.br
!
!   Version 1.0: July 17th, 2002
!
!
!   Indexes of error messages:
!
!     1 - All the available units for binary file operations are already in use.
!     2 - Incorrect file name.
!     3 - The argument <Action> must be 'read' or 'write'.
!     4 - Cannot open file.
!     5 - The specified unit is not available for binary file operations.
!     6 - The passed unit does not correspond to a binary file.
!     7 - End of file exceeded.
!     8 - The file was not opened for reading.
!     9 - The file was not opened for writing.

! Define the allowable range of unit indexes:
integer, parameter, private :: MinUnit = 90,   &
                               MaxUnit = 99

! Define the size (number of bytes) of the buffers (it must be a power of 2):
integer, parameter, private :: BufferExponent = 10
integer, parameter, private :: BufferSize = 2**BufferExponent

! Declare a logical array for assuring the hendling of binary files:
logical, private :: IsBinary(MinUnit:MaxUnit) = .false.

! Declare a logical array for verifying the I/O status of a file:
logical, private :: Writing(MinUnit:MaxUnit)

! Declare generic names for the routines:

interface BinRead ! Generic routine for reading from a binary file

  module procedure ReadString,       &   ! Single character variable
                   ReadInteger_1,    &   ! Single integer(1) variable
                   ReadInteger_2,    &   ! Single integer(2) variable
                   ReadInteger_4,    &   ! Single integer(4) variable
                   ReadReal_4,       &   ! Single real(4) variable
                   ReadReal_8,       &   ! Single real(8) variable
                   ReadLogical_1,    &   ! Single logical(1) variable
                   ReadLogical_2,    &   ! Single logical(2) variable
                   ReadLogical_4,    &   ! Single logical(4) variable
                   ReadComplex_4,    &   ! Single complex(4) variable
                   ReadComplex_8,    &   ! Single complex(8) variable
                   ReadSeqString,    &   ! Sequence of character variables
                   ReadSeqInteger_1, &   ! Sequence of integer(1) variables
                   ReadSeqInteger_2, &   ! Sequence of integer(2) variables
                   ReadSeqInteger_4, &   ! Sequence of integer(4) variables
                   ReadSeqReal_4,    &   ! Sequence of real(4) variables
                   ReadSeqReal_8,    &   ! Sequence of real(8) variables
                   ReadSeqLogical_1, &   ! Sequence of logical(1) variables
                   ReadSeqLogical_2, &   ! Sequence of logical(2) variables
                   ReadSeqLogical_4, &   ! Sequence of logical(4) variables
                   ReadSeqComplex_4, &   ! Sequence of complex(4) variables
                   ReadSeqComplex_8      ! Sequence of complex(8) variables

end interface BinRead

interface BinWrite ! Generic routine for writing to a binary file

  module procedure WriteString,       &   ! Single character variable
                   WriteInteger_1,    &   ! Single integer(1) variable
                   WriteInteger_2,    &   ! Single integer(2) variable
                   WriteInteger_4,    &   ! Single integer(4) variable
                   WriteReal_4,       &   ! Single real(4) variable
                   WriteReal_8,       &   ! Single real(8) variable
                   WriteLogical_1,    &   ! Single logical(1) variable
                   WriteLogical_2,    &   ! Single logical(2) variable
                   WriteLogical_4,    &   ! Single logical(4) variable
                   WriteComplex_4,    &   ! Single complex(4) variable
                   WriteComplex_8,    &   ! Single complex(8) variable
                   WriteSeqString,    &   ! Sequence of character variables
                   WriteSeqInteger_1, &   ! Sequence of integer(1) variables
                   WriteSeqInteger_2, &   ! Sequence of integer(2) variables
                   WriteSeqInteger_4, &   ! Sequence of integer(4) variables
                   WriteSeqReal_4,    &   ! Sequence of real(4) variables
                   WriteSeqReal_8,    &   ! Sequence of real(8) variables
                   WriteSeqLogical_1, &   ! Sequence of logical(1) variables
                   WriteSeqLogical_2, &   ! Sequence of logical(2) variables
                   WriteSeqLogical_4, &   ! Sequence of logical(4) variables
                   WriteSeqComplex_4, &   ! Sequence of complex(4) variables
                   WriteSeqComplex_8      ! Sequence of complex(8) variables

end interface BinWrite

contains

!==============================================================================

  integer function BinOpen( File, Action, IOStat, Append ) result ( Unit )

    character(*), intent(in)       :: File, Action
    integer, intent(out), optional :: IOStat
    logical, intent(in), optional  :: Append

    !
    !  ## This function opens a binary file for input or output operations and
    !     return the associated unit index ##
    !
    !  Arguments:
    !    Unit   - An index to be associated to the file
    !             (MinUnit <= Unit <= MaxUnit)
    !    File   - Name of the file to be opened
    !    Action - 'read' if the file is to be opened for input or
    !             'write' if it is to be opened for output.
    !

    integer Error
    character(7) Status, Position
    logical      Exists, Opened

    ! Search for a valid unit index to be opened:
    Unit = MinUnit - 1
    Opened = .true.
    do while ( Opened .and. (Unit < MaxUnit) )
      Unit = Unit + 1
      inquire( unit = Unit , opened = Opened )
    end do
    if (Opened) then
      if (present(IOStat)) then
        IOStat = 1
        return
      else
        write(*,'("Error handling binary file: all the available units '// &
                  'for binary file operations are already in use.")')
        stop
      end if
    end if

    ! Select the status of the file:
    select case ( trim(Action) )

      case ('read') ! The file must be opened for reading:

        ! Verify if the file already exists:
        inquire( file = File, exist = Exists )
        if (.not.Exists) then
          if (present(IOStat)) then
            IOStat = 2
            return
          else
            stop "Error handling binary file: incorrect file name."
          end if
        end if

        ! Define reading status:
        Status = 'old'
        Position = 'asis'
        Writing(Unit) = .false.

      case ('write') ! The file must be opened for writing:

        ! Define writing status:
        if (present(Append)) then
          if (Append) then
            Status = 'old'
            Position = 'append'
          else
            Status = 'replace'
            Position = 'asis'
          end if
        else
          Status = 'replace'
          Position = 'asis'
        end if

        Writing(Unit) = .true.

      case default   ! Wrong value of the argument <Action>:

        if (present(IOstat)) then
          IOStat = 3
          return
        else
          write(*,'("Error handling binary file: the argument <Action> '// &
                    'must be ''read'' or ''write''.")')
          stop
        end if

    end select

    ! Attempt to open the file:
    open( unit = Unit, file = File, status = Status, action = Action,  &
          form = 'binary', position = Position, iostat = Error )

    ! Verify if an error occurs during the attempt:
    if ( Error /= 0 ) then
      if (present(IOStat)) then
        IOStat = 4
        return
      else
        stop "Error handling binary file: cannot open file."
      endif
    end if

    ! Set the considered unit as a binary one:
    IsBinary(Unit) = .true.

  end function BinOpen

!==============================================================================

  subroutine BinClose( Unit, IOStat )

    integer, intent(in)            :: Unit
    integer, intent(out), optional :: IOStat

    !
    !  ## This subroutine closes a binary file ##
    !
    !  Argument:
    !    Unit - The index of the file to be closed.
    !

    ! Verify if the passed unit lies inside the correct range:
    if ( (Unit < MinUnit).or.(Unit > MaxUnit) ) then
      if (present(IOStat)) then
        IOStat = 5
        return
      else
        write(*,'("Error handling binary files: the specified unit is '// &
                  'not available for binary file operations.")')
        stop
      end if
    end if

    ! Verify if the file is already opened as a binary one:
    if ( IsBinary(Unit) ) then

      ! Close the file:
      close ( Unit )

      ! Erase the binary status:
      IsBinary(Unit) = .false.

    else

      if (present(IOStat)) then
        IOStat = 6
        return
      else
        write(*,'("Error handling binary files: the passed unit does '// &
                  'not correspond to a binary file.")')
        stop
      end if

    end if

  end subroutine BinClose

!==============================================================================

  subroutine ReadString( Unit, String, IOStat )

    integer, intent(in)            :: Unit
    character(*), intent(out)      :: String
    integer, intent(out), optional :: IOStat

    !
    !  ## This subroutine reads a character variable from a file buffer ##
    !
    !  Arguments:
    !    Unit   - The index of the considered file.
    !    String - A character variable of any length.
    !    IOStat - An optional output integer variable aiming at error handling
    !             (IOStat == 0 if the task is successful)
    !

    integer Status

    ! Verify if the file was not opened for reading:
    if ( Writing(Unit) ) then
      if (present(IOStat)) then
        IOStat = 8
        return
      else
        stop "Error handling binary file: the file was not opened for reading."
      end if
    end if

    read(Unit,iostat=Status) String

    if (present(IOStat)) then
      IOStat = Status
    else
      if (Status /= 0) stop "Error handling binary files: end of file exceeded."
    end if

  end subroutine ReadString

!==============================================================================

  subroutine WriteString( Unit, String, IOStat )

    integer, intent(IN)            :: Unit
    character(*), intent(IN)       :: String
    integer, intent(out), optional :: IOStat

    !
    !  ## This subroutine writes a character variable to a file buffer
    !     and flushes the buffer to an external file, if necessary ##
    !
    !  Arguments:
    !    Unit   - The index of the considered file.
    !    String - A character variable of any length.
    !    IOStat - An optional output integer variable aiming at error handling
    !             (IOStat == 0 if the task is successful)
    !

    ! Verify if the file was not opened for writing:
    if (.not.Writing(Unit)) then
      if (present(IOStat)) then
        IOStat = 9
        return
      else
        stop "Error handling binary file: the file was not opened for writing."
      end if
    end if

    ! Write the String to the file buffer:
    write(Unit) String

  end subroutine WriteString

!==============================================================================

  subroutine ReadInteger_1( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    integer(1), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(1) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadInteger_1

!==============================================================================

  subroutine WriteInteger_1( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    integer(1), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(1) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteInteger_1

!==============================================================================

  subroutine ReadInteger_2( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    integer(2), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(2) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadInteger_2

!==============================================================================

  subroutine WriteInteger_2( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    integer(2), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(2) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteInteger_2

!==============================================================================

  subroutine ReadInteger_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    integer(4), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(4) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadInteger_4

!==============================================================================

  subroutine WriteInteger_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    integer(4), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(4) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteInteger_4

!==============================================================================

  subroutine ReadReal_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    real(4), intent(out)           :: Value
    integer, intent(out), optional :: IOStat

    character(4) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadReal_4

!==============================================================================

  subroutine WriteReal_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    real(4), intent(in)            :: Value
    integer, intent(out), optional :: IOStat

    character(4) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteReal_4

!==============================================================================

  subroutine ReadReal_8( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    real(8), intent(out)           :: Value
    integer, intent(out), optional :: IOStat

    character(8) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadReal_8

!==============================================================================

  subroutine WriteReal_8( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    real(8), intent(in)            :: Value
    integer, intent(out), optional :: IOStat

    character(8) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteReal_8

!==============================================================================

  subroutine ReadLogical_1( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    logical(1), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(1) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadLogical_1

!==============================================================================

  subroutine WriteLogical_1( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    logical(1), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(1) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteLogical_1

!==============================================================================

  subroutine ReadLogical_2( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    logical(2), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(2) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadLogical_2

!==============================================================================

  subroutine WriteLogical_2( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    logical(2), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(2) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteLogical_2

!==============================================================================

  subroutine ReadLogical_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    logical(4), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(4) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadLogical_4

!==============================================================================

  subroutine WriteLogical_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    logical(4), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(4) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteLogical_4

!==============================================================================

  subroutine ReadComplex_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    complex(4), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(8) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadComplex_4

!==============================================================================

  subroutine WriteComplex_4( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    complex(4), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(8) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteComplex_4

!==============================================================================

  subroutine ReadComplex_8( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    complex(8), intent(out)        :: Value
    integer, intent(out), optional :: IOStat

    character(16) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value = transfer( String, Value )

  end subroutine ReadComplex_8

!==============================================================================

  subroutine WriteComplex_8( Unit, Value, IOStat )

    integer, intent(in)            :: Unit
    complex(8), intent(in)         :: Value
    integer, intent(out), optional :: IOStat

    character(16) String

    String = transfer( Value, String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteComplex_8

!==============================================================================

  subroutine ReadSeqString( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    character(*), intent(inout)    :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(Number*len(Value(1))) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqString

!==============================================================================

  subroutine WriteSeqString( Unit, Value, Number, IOStat )

    integer, intent(in)                 :: Unit
    character(*), intent(inout)         :: Value(*)
    integer, intent(in)                 :: Number
    integer, intent(out), optional :: IOStat

    character(Number*len(Value(1))) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqString

!==============================================================================

  subroutine ReadSeqInteger_1( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    integer(1), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqInteger_1

!==============================================================================

  subroutine WriteSeqInteger_1( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    integer(1), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqInteger_1

!==============================================================================

  subroutine ReadSeqInteger_2( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    integer(2), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(2*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqInteger_2

!==============================================================================

  subroutine WriteSeqInteger_2( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    integer(2), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(2*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqInteger_2

!==============================================================================

  subroutine ReadSeqInteger_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    integer(4), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(4*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqInteger_4

!==============================================================================

  subroutine WriteSeqInteger_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    integer(4), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(4*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqInteger_4

!==============================================================================

  subroutine ReadSeqReal_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    real(4), intent(out)           :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(4*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqReal_4

!==============================================================================

  subroutine WriteSeqReal_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    real(4), intent(in)            :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(4*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqReal_4

!==============================================================================

  subroutine ReadSeqReal_8( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    real(8), intent(out)           :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(8*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqReal_8

!==============================================================================

  subroutine WriteSeqReal_8( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    real(8), intent(in)            :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(8*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqReal_8

!==============================================================================

  subroutine ReadSeqLogical_1( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    logical(1), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqLogical_1

!==============================================================================

  subroutine WriteSeqLogical_1( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    logical(1), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqLogical_1

!==============================================================================

  subroutine ReadSeqLogical_2( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    logical(2), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(2*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqLogical_2

!==============================================================================

  subroutine WriteSeqLogical_2( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    logical(2), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(2*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqLogical_2

!==============================================================================

  subroutine ReadSeqLogical_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    logical(4), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(4*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqLogical_4

!==============================================================================

  subroutine WriteSeqLogical_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    logical(4), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(4*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqLogical_4

!==============================================================================

  subroutine ReadSeqComplex_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    complex(4), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(8*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqComplex_4

!==============================================================================

  subroutine WriteSeqComplex_4( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    complex(4), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(8*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqComplex_4

!==============================================================================

  subroutine ReadSeqComplex_8( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    complex(8), intent(out)        :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(16*Number) String

    if ( present(IOStat) ) then
      call ReadString( Unit, String, IOStat )
    else
      call ReadString( Unit, String )
    end if
    Value(1:Number) = transfer( String, Value(1:Number) )

  end subroutine ReadSeqComplex_8

!==============================================================================

  subroutine WriteSeqComplex_8( Unit, Value, Number, IOStat )

    integer, intent(in)            :: Unit
    complex(8), intent(in)         :: Value(*)
    integer, intent(in)            :: Number
    integer, intent(out), optional :: IOStat

    character(16*Number) String

    String = transfer( Value(1:Number), String )
    if (present(IOStat)) then
      call WriteString( Unit, String, IOStat )
    else
      call WriteString( Unit, String )
    end if

  end subroutine WriteSeqComplex_8

!==============================================================================

end module BinFiles