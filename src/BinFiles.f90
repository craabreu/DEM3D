module BinFiles

!
!
!       ## Handling of Sequential Binary Files ##
!
!   Author: Charlles R. A. Abreu (Rio de Janeiro, Brazil)
!   Emails: rubber@eq.ufrj.br / craabreu@ig.com.br
!
!   Version 1.0: July 17th, 2002
!
!   Note: This module was only tested with the compiler
!         Absoft Pro Fortran for Linux.
!
!
!   Numbers of error messages:
!
!     1 - All the available units for binary operations are already in use.
!     2 - Incorrect file name.
!     3 - The argument <Action> must be 'read' or 'write'.
!     4 - Cannot open file.
!     5 - The specified unit is not available for binary operations.
!     6 - The passed unit does not correspond to a binary file.
!     7 - End of file exceeded.
!     8 - The file was not opened for reading.
!     9 - The file was not opened for writing.
!

implicit none

! Define the allowable range of unit indexes:
integer, parameter, private :: MinUnit = 90,   &
                               MaxUnit = 99

! Define the size (number of bytes) of the buffers (it must be a power of 2):
integer, parameter, private :: BufferExponent = 10
integer, parameter, private :: BufferSize = 2**BufferExponent

! Define a structured type for containing parameters of a binary file:
type TBinFile
  sequence
  integer               FilePos        ! Current position (record) in a file
  character(BufferSize) Buffer         ! Buffer for containing input/output data
  integer               BufferPos      ! Current position in the buffer
  integer               LastBufferPos  ! Last position of an unfulfilled buffer
  logical               Writing        ! Tag for identifying writing/reading status
end type TBinFile

! Declare an array for containing various binary files:
type(TBinFile), private :: BinFile(MinUnit:MaxUnit)

! Declare a logical array for assuring the hendling of binary files:
logical, private :: IsBinary(MinUnit:MaxUnit) = .false.

! Declare generic names for the routines:

interface BinRead ! Generic routine for reading from a binary file

  module procedure ReadString,       &   ! Reads a single character variable
                   ReadInteger_1,    &   ! Reads a single integer(1) variable
                   ReadInteger_2,    &   ! Reads a single integer(2) variable
                   ReadInteger_4,    &   ! Reads a single integer(4) variable
                   ReadReal_4,       &   ! Reads a single real(4) variable
                   ReadReal_8,       &   ! Reads a single real(8) variable
                   ReadLogical_1,    &   ! Reads a single logical(1) variable
                   ReadLogical_2,    &   ! Reads a single logical(2) variable
                   ReadLogical_4,    &   ! Reads a single logical(4) variable
                   ReadComplex_4,    &   ! Reads a single complex(4) variable
                   ReadComplex_8,    &   ! Reads a single complex(8) variable
                   ReadSeqString,    &   ! Reads a sequence of character variables
                   ReadSeqInteger_1, &   ! Reads a sequence of integer(1) variables
                   ReadSeqInteger_2, &   ! Reads a sequence of integer(2) variables
                   ReadSeqInteger_4, &   ! Reads a sequence of integer(4) variables
                   ReadSeqReal_4,    &   ! Reads a sequence of real(4) variables
                   ReadSeqReal_8,    &   ! Reads a sequence of real(8) variables
                   ReadSeqLogical_1, &   ! Reads a sequence of logical(1) variables
                   ReadSeqLogical_2, &   ! Reads a sequence of logical(2) variables
                   ReadSeqLogical_4, &   ! Reads a sequence of logical(4) variables
                   ReadSeqComplex_4, &   ! Reads a sequence of complex(4) variables
                   ReadSeqComplex_8      ! Reads a sequence of complex(8) variables

end interface

interface BinWrite ! Generic routine for writing to a binary file

  module procedure WriteString,       &   ! Writes a single character variable
                   WriteInteger_1,    &   ! Writes a single integer(1) variable
                   WriteInteger_2,    &   ! Writes a single integer(2) variable
                   WriteInteger_4,    &   ! Writes a single integer(4) variable
                   WriteReal_4,       &   ! Writes a single real(4) variable
                   WriteReal_8,       &   ! Writes a single real(8) variable
                   WriteLogical_1,    &   ! Writes a single logical(1) variable
                   WriteLogical_2,    &   ! Writes a single logical(2) variable
                   WriteLogical_4,    &   ! Writes a single logical(4) variable
                   WriteComplex_4,    &   ! Writes a single complex(4) variable
                   WriteComplex_8,    &   ! Writes a single complex(8) variable
                   WriteSeqString,    &   ! Writes a sequence of character variables
                   WriteSeqInteger_1, &   ! Writes a sequence of integer(1) variables
                   WriteSeqInteger_2, &   ! Writes a sequence of integer(2) variables
                   WriteSeqInteger_4, &   ! Writes a sequence of integer(4) variables
                   WriteSeqReal_4,    &   ! Writes a sequence of real(4) variables
                   WriteSeqReal_8,    &   ! Writes a sequence of real(8) variables
                   WriteSeqLogical_1, &   ! Writes a sequence of logical(1) variables
                   WriteSeqLogical_2, &   ! Writes a sequence of logical(2) variables
                   WriteSeqLogical_4, &   ! Writes a sequence of logical(4) variables
                   WriteSeqComplex_4, &   ! Writes a sequence of complex(4) variables
                   WriteSeqComplex_8      ! Writes a sequence of complex(8) variables

end interface

contains

!==============================================================================

  integer function BinOpen( File, Action, IOStat ) result ( Unit )

    character(*), intent(in)       :: File, Action
    integer, intent(out), optional :: IOStat

    !
    !  ## This function opens a binary file for input or output operations and
    !     return the associated unit index ##
    !
    !  Arguments:
    !    Unit   - An index to be associated to the file (MinUnit <= Unit <= MaxUnit)
    !    File   - Name of the file to be opened
    !    Action - 'read' if the file is to be opened for input or
    !             'write' if it is to be opened for output.
    !

    integer Error
    character(7) Status
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
        stop "Error handling binary file: all the available units for binary operations are already in use."
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
        BinFile(Unit)%Writing = .false.

      case ('write') ! The file must be opened for writing:

        ! Define writing status:
        Status = 'replace'
        BinFile(Unit)%Writing = .true.

      case default   ! Wrong value of the argument <Action>:

        if (present(IOstat)) then
          IOStat = 3
          return
        else
          stop "Error handling binary file: the argument <Action> must be 'read' or 'write'."
        end if

    end select

    ! Initialize the file and buffer positions:
    BinFile(Unit)%FilePos   = 0
    BinFile(Unit)%BufferPos = 0

    ! Initialize the last buffer position:
    BinFile(Unit)%LastBufferPos = BufferSize

    ! Attempt to open the file:
    open( unit = Unit, file = File, status = Status, action = Action,  &
          form = 'unformatted', access='direct', recl = BufferSize,    &
          iostat = Error )

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

    type (TBinFile) File

    ! Verify if the passed unit lies inside the correct range:
    if ( (Unit < MinUnit).or.(Unit > MaxUnit) ) then
      if (present(IOStat)) then
        IOStat = 5
        return
      else
        write(*,'("Error handling binary files: the specified unit is not available for binary operations.")')
        stop
      end if
    end if

    ! Verify if the file is already opened as a binary one:
    if ( IsBinary(Unit) ) then

      ! Save the file parameters in a single variable:
      File = BinFile(Unit)

      ! Flush the last buffer, if it is necessary:
      if ( File%Writing .and. (File%BufferPos /= 0) ) then
        call FlushBuffer( Unit, File )
      end if

      ! Close the file:
      close ( Unit )

      ! Erase the binary status:
      IsBinary(Unit) = .false.

    else

      if (present(IOStat)) then
        IOStat = 6
        return
      else
        stop "Error handling binary files: the passed unit does not correspond to a binary file."
      end if

    end if

  end subroutine BinClose

!==============================================================================

  subroutine FlushBuffer( Unit, File )

    integer, intent(in)            :: Unit
    type (TBinFile), intent(inout) :: File

    integer Size, Pos, LastRec, Error
    character(256) FileName

    ! Inquire the name of the file associated to Unit:
    inquire( Unit, Name = FileName )

    ! Reduce the record length to a half and point out the beggining of the buffer:
    Size = BufferSize / 2
    LastRec = 2*File%FilePos
    Pos = 0
    do while ( Size > 0 )

      if ( File%BufferPos - Pos >= Size ) then

        ! Close and reopen the file with the new record length:
        close(Unit)
        open( unit = Unit, file = FileName, status = 'old', action = 'write', &
              form = 'unformatted', access = 'direct', recl = Size )

        ! Flush a slice:
        LastRec = LastRec + 1
        write(Unit,rec=LastRec) File%Buffer(Pos+1:Pos+Size)
        Pos  = Pos + Size

      end if

      ! Reduce the record length to a half:
      Size = Size / 2
      LastRec = 2*LastRec

    end do

  end subroutine FlushBuffer

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
    if ( BinFile(Unit)%Writing ) then
      if (present(IOStat)) then
        IOStat = 8
        return
      else
        stop "Error handling binary file: the file was not opened for reading."
      end if
    end if

    call ReadFromBuffer( Unit, String, BinFile(Unit), Status )

    if (present(IOStat)) then
      IOStat = Status
    else
      if (Status /= 0) stop "Error handling binary files: end of file exceeded."
    end if

    contains

    !--------------------------------------------------------------------------

      recursive subroutine ReadFromBuffer( Unit, String, File, Status )

        integer, intent(in)            :: Unit
        character(*), intent(out)      :: String
        type (TBinFile), intent(inout) :: File
        integer, intent(out)           :: Status

        integer Pos, Length, StringEnd, LastPos, Slice

        ! Save the current position in the buffer:
        Pos = File%BufferPos

        ! Verify if it is necessary to read the buffer from the external file:
        if (Pos == 0) then
          call ReadBufferFromFile( Unit, File )
          Pos = 1
        end if

        ! Save some parameters:
        Length    = len(String)         ! Length of the argument String
        StringEnd = Pos + Length - 1    ! Possible last position after placing
        LastPos   = File%LastBufferPos  ! Last available position of the buffer

        ! Verify if it is possible to read the whole string from the buffer:
        if ( StringEnd <= LastPos ) then

          ! Read the string from the buffer and update the current position:
          String = File%Buffer(Pos:StringEnd)
          File%BufferPos = mod(StringEnd + 1, BufferSize + 1)

          ! Set the variable representing the error status (no errors):
          Status = 0

        else if ( LastPos < BufferSize ) then

          ! Take a slice from the buffer and place it in the beginning of the string:
          Slice = LastPos - Pos + 1
          String(1:Slice) = File%Buffer(Pos:LastPos)

          ! Set the variable representing the error status (error occured):
          Status = 7

        else

          ! Take a slice from the end of the buffer and place it in the beggining of the string:
          Slice = BufferSize - Pos + 1
          String(1:Slice) = File%Buffer(Pos:BufferSize)

          ! Zero the current position and read the remaining of the string from the buffer:
          File%BufferPos = 0
          call ReadFromBuffer( Unit, String(Slice+1:Length), File, Status )

        end if

      end subroutine ReadFromBuffer

    !--------------------------------------------------------------------------

      subroutine ReadBufferFromFile( Unit, File )

        integer, intent(in)            :: Unit
        type (TBinFile), intent(inout) :: File

        integer Error, Rec, Size, Pos
        character(256) FileName

        ! Attempt to read the next record of the file:
        Rec = File%FilePos + 1
        read(Unit,rec=Rec,iostat=Error) File%Buffer

        ! An error indicates that the last record was reached:
        if ( Error /= 0 ) then

          ! Inquire the name of the file associated to Unit:
          inquire( Unit, Name = FileName )

          ! Reduce the record length to a half and point out the beggining of the buffer:
          Size = BufferSize / 2
          Pos  = 0
          do while ( Size > 0 )

            ! Close and reopen the file with the new record length:
            close(Unit)
            open( unit = Unit, file = FileName, status = 'old',   &
                  form = 'unformatted', access = 'direct', recl = Size )

            ! Attempt to read the new record:
            Rec = 2*Rec - 1
            read(Unit,rec=Rec,iostat=Error) File%Buffer(Pos+1:Pos+Size)

            ! If the attempt was successful, update the file and buffer positions:
            if (Error == 0) then
              Pos  = Pos + Size
              Rec  = Rec + 1
            end if

            ! Reduce the record length to a half:
            Size = Size / 2

          end do

          ! Set up the last buffer position:
          File%LastBufferPos = Pos

          ! Set up the current buffer position:
          File%BufferPos = 1

        end if

        ! Update the file position:
        File%FilePos = Rec

      end subroutine ReadBufferFromFile

    !--------------------------------------------------------------------------

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
    if (.not.BinFile(Unit)%Writing) then
      if (present(IOStat)) then
        IOStat = 9
        return
      else
        stop "Error handling binary file: the file was not opened for writing."
      end if
    end if

    ! Write the String to the file buffer:
    call WriteToBuffer( Unit, String, BinFile(Unit) )

    contains

    !--------------------------------------------------------------------------

      recursive subroutine WriteToBuffer( Unit, String, File )

        integer, intent(in)            :: Unit
        character(*), intent(in)       :: String
        type (TBinFile), intent(inout) :: File

        integer Pos, Length, LastPos, Slice

        ! Save some parameters:
        Pos     = File%BufferPos     ! Current position in the buffer
        Length  = len(String)        ! Length of the argument String
        LastPos = Pos + Length       ! Possible last position after placing

        ! Verify if it is possible to write the whole string to the buffer:
        if ( LastPos <= BufferSize ) then

          ! Write the string to the buffer and save the new current position:
          File%Buffer(Pos+1:LastPos) = String
          File%BufferPos = LastPos

        else

          ! Take a slice of the string and write it until the end of the buffer:
          Slice = BufferSize - Pos
          File%Buffer(Pos+1:BufferSize) = String(1:Slice)

          ! Flush the buffer to the external file:
          File%FilePos = File%FilePos + 1
          write(Unit,rec=File%FilePos) File%Buffer

          ! Zero the current position and write the remaining of the string:
          File%BufferPos = 0
          call WriteToBuffer( Unit, String(Slice+1:Length), File )

        end if

      end subroutine WriteToBuffer

    !--------------------------------------------------------------------------

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