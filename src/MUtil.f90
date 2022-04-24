module MUtil

use MTypes
use MDimension

real(rb), parameter :: Pi = 3.1415926535897932384626433832795_rb, &
                       TwoPi = 2.0_rb*Pi

! Global variables for the random number generator:
real(8)    :: aam
integer(4) :: Dummy, iix = -1_4, iiy = -1_4, iik

contains

!===============================================================================

  real(8) function ranf( idum )

    implicit none
    integer(4), intent(inout) :: idum

    ! "Minimal" random number generator of Park and Miller combined with
    ! a Marsaglia shift sequence. Returns a uniform random deviate between
    ! 0.0 and 1.0 (exclusive of the endpoint values). This fully portable,
    ! generator has the "traditional" (not Fortran 90) calling sequence with
    ! a scalar random deviate as the returned function value: call with idum
    ! a negative integer to initialize; thereafter, do not alter idum except
    ! to reinitialize.
    ! The period of this generator is about 3.1E+18

    integer(4), parameter :: IA = 16807, IM = 2147483647,  &
                             IQ = 127773, IR = 2836

    if ( (idum <= 0) .or. (iiy < 0) ) then   !  Initialize.
      aam  = nearest(1.0,-1.0)/IM
      iiy  = ior(ieor(888889999,abs(idum)),1)
      iix  = ieor(777755555,abs(idum))
      idum = abs(idum)+1                     ! Set idum positive.
    end if
    iix = ieor(iix,ishft(iix,13))         ! Marsaglia shift sequence with
    iix = ieor(iix,ishft(iix,-17))        ! period 2^32-1.
    iix = ieor(iix,ishft(iix,5))
    iik = iiy/IQ                          ! Park-Miller sequence by Schrage's
    iiy = IA*(iiy-iik*IQ)-IR*iik          ! method,period 2^31-2.
    if (iiy < 0) iiy = iiy+IM
    ranf = aam*ior(iand(IM,ieor(iix,iiy)),1) ! Combine the two generators with
                                             ! masking to ensure nonzero value.
  end function ranf

!===============================================================================

  function RandomUnitVector( Dummy ) result ( U )

    integer(4), intent(inout) :: Dummy
    real(rb) :: U(3)

    !
    ! This subroutine generate a random unit vector uniformly distributed 
    ! according to the Marsaglia's algorithm (Allen & Tildesley, p.349)
    !

    real(rb) :: Dot, Sq, X, Y

    Dot = 1.1_rb
    do while ( Dot > 1.0_rb )
      X = 1.0_rb - 2.0_rb*Ranf(Dummy)
      Y = 1.0_rb - 2.0_rb*Ranf(Dummy)
      Dot = X*X + Y*Y
    end do
    Sq = 2.0_rb*sqrt(1.0_rb - Dot)
    U = (/ X*Sq, Y*Sq, 1.0_rb - 2.0_rb*Dot /)

  end function RandomUnitVector

!===============================================================================

  function RandomList( N, idum ) result ( List )

    integer(ib), intent(in)   :: N
    integer(4), intent(inout) :: idum
    integer(ib) :: List(N)

    type TLocalList
      integer(ib) :: I
      real(rb)    :: R
    end type TLocalList

    type (TLocalList) :: LocalList(N), Aux
    integer :: i, j, M
    logical :: Replaced

    do i = 1 , N
      LocalList(i)%R = ranf(idum)
    end do

    forall (i=1:N) LocalList(i)%I = i

    call Sort( LocalList%R, LocalList%I, N)

    Replaced = .true.
    M = N - 1
    do while (Replaced)
      Replaced = .false.
      do i = 1 , M
        j = i + 1
        if ( LocalList(i)%R > LocalList(j)%R ) then
          Aux = LocalList(i)
          LocalList(i) = LocalList(j)
          LocalList(j) = Aux
          Replaced = .true.
        end if
      end do
      M = M - 1
    end do

    List = LocalList%I

  end function RandomList

!===============================================================================

  function FileSeries( FileName, Index, NumberOfDigits ) result( Name )

    integer(ib), intent(in)      :: Index, NumberOfDigits
    character(*), intent(in) :: FileName

    character(len(FileName)) :: Name

    integer(ib) :: i, j, D, Digit(0:NumberOfDigits-1), Ext, Last
    character   :: Str
    character(10) :: Extension

    do i = NumberOfDigits-1, 0, -1
      D = Index
      do j = NumberOfDigits-1, i + 1, -1
        D = D - 10**j*Digit(j)
      end do
      D = D / 10**i
      Digit(i) = D
    end do

    Last = len_trim( FileName )
    Ext  = scan( FileName, '.' )

    if ( Ext /= 0 ) then
      Extension  = FileName(Ext:Last)
      Name = FileName(1:Ext-1)//'_'
    else
      Extension = ''
      Name = FileName(1:Last)//'_'
    end if

    do i = NumberOfDigits-1, 0, -1
      write(Str,'(i1)') Digit(i)
      Name = trim(Name)//Str
    end do
    Name = trim(Name)//trim(Extension)

  end function FileSeries

!===============================================================================

  subroutine GramSchmidt( U, W )

    real(rb), intent(inout) :: U(dim), W(dim)

    U = U/sqrt(sum(U*U))
    W = W - (U .dot. W)*U
    W = W/sqrt(sum(W*W))

  end subroutine GramSchmidt

!===============================================================================

  subroutine Sort( X, Ind, N )

    integer(ib), intent(in)  :: N
    integer(ib), intent(out) :: Ind(N)
    real(rb), intent(inout)  :: X(N)

    integer(ib) :: i

    ! A recursive sorting routine Adapted from The design of Well-Structured
    ! and Correct Programs, S. Alagic, Springer-Verlag, 1978.

    forall (i=1:N) Ind(i) = i
    call Qsort( X, Ind, 1, N, N )

    contains

      recursive subroutine Qsort( X, Ind, M, N, NN )

        integer(ib), intent(in) :: M, N, NN
        integer(ib), intent(inout) :: Ind(NN)
        real(rb), intent(inout) :: X(NN)

        integer :: I, J

        if (M < N) then
          call Partit(X, I, J, M, N, Ind, NN) ! divide in two
          call Qsort(X, Ind, M, J, NN)        ! Sort left part
          call Qsort(X, Ind, I, N, NN)        ! Sort right part
        end if

      end subroutine Qsort

      subroutine Partit( A, I, J, Left, Right, Ind, NN )

        integer(ib), intent(in) :: Left, Right, NN
        integer(ib), intent(inout) :: I, J
        real(rb), intent(inout) :: A(NN)
        integer(ib), intent(inout) :: Ind(NN)

        integer(ib) :: IAux
        real(rb) :: Pivot, Aux

        Pivot = A( (Left + Right)/2 )
        I = Left
        J = right
        do while (I <= J)
          do while (A(I) < Pivot)
            I = I + 1
          end do
          do while (Pivot < A(J))
            J = J - 1
          end do
          if (I <= J) then
            Aux  = A(I)
            A(I) = A(J)
            A(J) = Aux
            IAux = Ind(I)
            Ind(I) = Ind(J)
            Ind(J) = IAux
            I = I + 1
            J = J - 1
          end if
        end do

      end subroutine Partit

  end subroutine Sort

!===============================================================================

  function Exists( FileName ) result( Status )

    character(*) FileName

    logical :: Status

    inquire( file = FileName, exist = Status )

  end function Exists

!===============================================================================

  function ArrayToStr( a ) result ( s )
    character, intent(in) :: a(:)
    character(size(a)) :: s
    integer :: i
    do i = 1 , size(a)
      s(i:i) = a(i)
    end do
  end function ArrayToStr

!===============================================================================

  function UpperCase( String ) result ( UpCase )

    character(*), intent(in) :: String
    character(len(String)) :: UpCase

    integer :: i, Code
    character :: C

    do i = 1, len(String)
      C = String(i:i)
      Code = ichar(C)
      if ( (Code >= 97) .and. (Code <= 122) ) then
        UpCase(i:i) = char(Code-32)
      else
        UpCase(i:i) = C
      end if
    end do

  end function UpperCase

!===============================================================================

end module MUtil