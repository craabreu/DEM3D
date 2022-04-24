module MContacts

use MTypes
use MDimension
use BinFiles

implicit none

!         This module contais subroutines that handle the contacts
! between pairs  of particles and  between particles  and walls in
! DEM simulations. For each pair, one needs to know the tangential
! displacement since to contact time.
!
! Author        : Charlles R. A. Abreu
! E-mail        : rubber@eq.ufrj.br
! First Version : January 26th, 2002
! Revision      : August 25th, 2003
! Revision      : September 23th, 2003
!

type, private :: TPartner
  integer(ib) :: Index
  real(rb)    :: Vector1(dim)
  real(rb)    :: Vector2(dim)
  type(TPartner), pointer :: Next => null()
end type TPartner

type, private :: TMainParticle
  type(TPartner), pointer :: FirstPartner => null()
  type(TPartner), pointer :: FirstWall => null()
end type TMainParticle

type(TMainParticle), allocatable, private :: MainParticle(:)

type(TPartner), pointer, private :: CurrentPair, PreviousPair

integer(ib), private :: imin, imax, iwall

contains

!===============================================================================

  subroutine InitializeContacts( NP )

    integer(ib), intent(in) :: NP

    allocate( MainParticle(NP) )

  end subroutine InitializeContacts

!===============================================================================

  subroutine ResizeContactList( OldNP, NewNP )

    integer(ib), intent(in) :: OldNP, NewNP

    type (TMainParticle) :: Aux(NewNP)

    if ( NewNP > OldNP ) then
      Aux(1:OldNP) = MainParticle
    else
      Aux = MainParticle(1:NewNP)
    end if

    deallocate( MainParticle )

    allocate( MainParticle(NewNP) )

    MainParticle = Aux

  end subroutine ResizeContactList

!===============================================================================

  subroutine SaveContactList( Unit )

    integer, intent(in) :: Unit

    integer(ib) :: NP, i, j, NPartners
    type(TPartner), pointer :: First, Current

    NP = size(MainParticle)
    call BinWrite( Unit, NP )
    do i = 1_ib, NP

      do j = 1_ib, 2_ib

        select case (j)
          case (1_ib)
            First => MainParticle(i)%FirstPartner
          case (2_ib)
            First => MainParticle(i)%FirstWall
        end select

        ! Count partners:
        Current => First
        NPartners = 0_ib
        do while (associated(Current))
          NPartners = NPartners + 1_ib
          Current => Current%Next
        end do

        call BinWrite(Unit,NPartners)

        Current => First
        do while (associated(Current))
          call BinWrite(Unit,Current%Index)
          call BinWrite(Unit,Current%Vector1,dim)
          call BinWrite(Unit,Current%Vector2,dim)
          Current => Current%Next
        end do

      end do ! j

    end do ! i

  end subroutine SaveContactList

!===============================================================================

  subroutine ReadContactList( Unit )

    integer, intent(in) :: Unit

    integer(ib) :: NP, i, j, k, NPartners
    type(TPartner), pointer :: First, Current

    call BinRead( Unit, NP )
    allocate( MainParticle(NP) )

    do i = 1_ib, NP

      do j = 1_ib, 2_ib

        call BinRead(Unit,NPartners)

        if (NPartners > 0_ib) then

          allocate(Current)
          First => Current
          do k = 1_ib, NPartners - 1_ib
            call BinRead(Unit,Current%Index)
            call BinRead(Unit,Current%Vector1,dim)
            call BinRead(Unit,Current%Vector2,dim)
            allocate(Current%Next)
            Current => Current%Next
          end do ! k
          call BinRead(Unit,Current%Index)
          call BinRead(Unit,Current%Vector1,dim)
          call BinRead(Unit,Current%Vector2,dim)

        else

          First => null()

        end if

        select case (j)
          case (1_ib)
            MainParticle(i)%FirstPartner => First
          case (2_ib)
            MainParticle(i)%FirstWall => First
        end select

      end do ! j

    end do ! i

  end subroutine ReadContactList

!===============================================================================

  subroutine ChangeIndexInContactList( NP, i, NPSC, Neighbour )

    integer(ib), intent(in) :: NP, i, NPSC, Neighbour(NPSC)

    ! This subroutine changes the index "NP" in the contact list by another
    ! index "iNew", making the convenient list manipulations. For this routine
    ! to function properly,  it is necessary to call beforehand the subroutine
    ! RemoveParticleFromContactList( i ).

    type(TPartner), pointer :: Current, Previous, Point, First, Last
    integer(ib) :: jj, j

    ! Search for "NP" as secondary partner in the list and,
    ! where it appears, replace it by "i":
    First => null()
    do jj = 1_ib, NPSC

      j = Neighbour(jj)
      Current => MainParticle(j)%FirstPartner
      if (associated(Current)) then

        if (j < i) then

          Previous => null()
          Point => null()
          do while (associated(Current%Next))
            Previous => Current
            if (Current%Index < i) Point => Current
            Current => Current%Next
          end do
          if (Current%Index == NP) then
            if (associated(Previous)) then
              Previous%Next => Current%Next
              Current%Index = i
              if (associated(Point)) then
                Current%Next => Point%Next
                Point%Next => Current
              else
                Current%Next => MainParticle(j)%FirstPartner
                MainParticle(j)%FirstPartner => Current
              end if
            else
              Current%Index = i
            end if
          end if

        else

          Previous => null()
          do while (associated(Current%Next))
            Previous => Current
            Current => Current%Next
          end do
          if (Current%Index == NP) then
            Current%Index = j
            if (associated(Previous)) then
              Previous%Next => null()
            else
              MainParticle(j)%FirstPartner => null()
            end if
            if (associated(First)) then
              Last%Next => Current
              Last => Last%Next
            else
              First => Current
              Last => First
            end if
          end if

        end if

      end if

    end do ! j

    MainParticle(i)%FirstPartner => First
    MainParticle(i)%FirstWall => MainParticle(NP)%FirstWall
    MainParticle(NP)%FirstWall => null()

  end subroutine ChangeIndexInContactList

!===============================================================================

  subroutine PrintList( Unit, i )

    integer, intent(in) :: Unit
    integer(ib), intent(in) :: i

    type(TPartner), pointer :: Current

    write(Unit,'("List of particle ",I3,":")') i
    Current => MainParticle(i)%FirstPartner
    write(Unit,fmt='("Particles: ")',advance='no')
    if (associated(Current)) then
      do while (associated(Current))
        write(Unit,fmt='(I3,"(",3D15.8,")")',advance='no') Current%Index,   &
                                                           Current%Vector1, &
                                                           Current%Vector2
        Current => Current%Next
      end do
      write(Unit,*)
    else
      write(Unit,*) 'Empty list'
    end if

    Current => MainParticle(i)%FirstWall
    write(Unit,fmt='("Walls: ")',advance='no')
    if (associated(Current)) then
      do while (associated(Current))
        write(Unit,fmt='(I3,"(",3D15.8,")")',advance='no') Current%Index,   &
                                                           Current%Vector1, &
                                                           Current%Vector2
        Current => Current%Next
      end do
      write(Unit,*)
    else
      write(Unit,*) 'Empty list'
    end if
    write(Unit,*)

  end subroutine PrintList

!===============================================================================

  logical function ValidList()

    integer(ib) :: i, j, k, m
    type(TPartner), pointer :: Current

    ValidList = .true.
    do i = 1, size(MainParticle)
      do m = 1, 2
        if (m==1) then
          Current => MainParticle(i)%FirstPartner
        else
          Current => MainParticle(i)%FirstWall
        end if
        if (associated(Current)) then
          j = Current%Index
          do while (associated(Current%Next))
            Current => Current%Next
            k = Current%Index
            if (k < j) then
              ValidList = .false.
              return
            end if
            j = k
          end do ! while
        end if
      end do ! m
    end do ! i

  end function ValidList

!===============================================================================

  subroutine VerifyPairContact( i, j, Contacting, Vector1, Vector2 )

    integer(ib), intent(in)  :: i, j
    logical, intent(out) :: Contacting
    real(rb),    intent(out) :: Vector1(dim), Vector2(dim)

    ! Verify the minimum and maximum indexes:
    imin = min(i,j)
    imax = max(i,j)

    ! Point to the first partner of the main particle (if it exists):
    CurrentPair  => MainParticle(imin)%FirstPartner

    ! Verify if the particle "imin" has at least one partner:
    Contacting = associated(CurrentPair)

    if (Contacting) then

      ! Search for the partner "imax" or the end of the list:
      PreviousPair => CurrentPair
      do while ( (CurrentPair%Index < imax).and.associated(CurrentPair%Next) )
        PreviousPair => CurrentPair
        CurrentPair  => CurrentPair%Next
      end do ! while

      ! Verify if the partner "imax" was found:
      Contacting = ( CurrentPair%Index == imax )

      ! Return the vectors if there is a contact:
      if (Contacting) then
        Vector1 = CurrentPair%Vector1
        Vector2 = CurrentPair%Vector2
      else
        Vector1 = 0.0_rb
        Vector2 = 0.0_rb
      end if

    else

      Vector1 = 0.0_rb
      Vector2 = 0.0_rb

    end if

  end subroutine VerifyPairContact

!===============================================================================

  subroutine UpdatePairContact( Vector1, Vector2, Contacting )

    real(rb), intent(in)    :: Vector1(dim), Vector2(dim)
    logical, intent(in) :: Contacting

    type (TPartner), pointer :: Aux

    if (Contacting) then

      CurrentPair%Vector1 = Vector1
      CurrentPair%Vector2 = Vector2

    else

      ! Allocate the auxiliary pointer and set up its target:
      allocate(Aux)
      Aux%Index   =  imax
      Aux%Vector1 =  Vector1
      Aux%Vector2 =  Vector2
      Aux%Next    => null()
  
      if ( associated( CurrentPair ) ) then
  
        if ( CurrentPair%Index < imax ) then
  
          CurrentPair%Next => Aux
  
        else
  
          Aux%Next => CurrentPair
          if ( CurrentPair%Index == PreviousPair%Index ) then
            MainParticle(imin)%FirstPartner => Aux
          else
            PreviousPair%Next => Aux
          end if
  
        end if
  
      else
  
        MainParticle(imin)%FirstPartner => Aux
  
      end if

    end if

  end subroutine UpdatePairContact

!===============================================================================

  subroutine RemovePairContact

    if ( PreviousPair%Index == imax ) then
      MainParticle(imin)%FirstPartner => CurrentPair%Next
    else
      PreviousPair%Next => CurrentPair%Next
    end if

    deallocate( CurrentPair )

  end subroutine RemovePairContact

!===============================================================================

  subroutine VerifyWallContact( i, w, Contacting, Vector1, Vector2 )

    integer(ib), intent(in)  :: i, w
    logical, intent(out) :: Contacting
    real(rb),    intent(out) :: Vector1(dim), Vector2(dim)

    ! Verify the minimum and maximum indexes:
    imin  = i
    iwall = w

    ! Point to the first Wall of the main particle (if it exists):
    CurrentPair  => MainParticle(imin)%FirstWall

    ! Verify if the particle "imin" has at least one Wall:
    Contacting = associated(CurrentPair)

    if (Contacting) then

      ! Search for the Wall "iwall" or the end of the list:
      PreviousPair => CurrentPair
      do while ( (CurrentPair%Index < iwall).and.associated(CurrentPair%Next) )
        PreviousPair => CurrentPair
        CurrentPair  => CurrentPair%Next
      end do ! while

      ! Verify if the Wall "iwall" was found:
      Contacting = ( CurrentPair%Index == iwall )

      ! Return the vectors if there is a contact:
      if (Contacting) then
        Vector1 = CurrentPair%Vector1
        Vector2 = CurrentPair%Vector2
      else
        Vector1 = 0.0_rb
        Vector2 = 0.0_rb
      end if

    else

      Vector1 = 0.0_rb
      Vector2 = 0.0_rb

    end if

  end subroutine VerifyWallContact

!===============================================================================

  subroutine UpdateWallContact( Vector1, Vector2, Contacting )

    real(rb), intent(in)    :: Vector1(dim), Vector2(dim)
    logical, intent(in) :: Contacting

    type (TPartner), pointer :: Aux

    if (Contacting) then

      CurrentPair%Vector1 = Vector1
      CurrentPair%Vector2 = Vector2

    else

      ! Allocate the auxiliary pointer and set up its target:
      allocate(Aux)
      Aux%Index   =  iwall
      Aux%Vector1 =  Vector1
      Aux%Vector2 =  Vector2
      Aux%Next    => null()

      if ( associated( CurrentPair ) ) then

        if ( CurrentPair%Index < iwall ) then

          CurrentPair%Next => Aux

        else

          Aux%Next => CurrentPair
          if ( CurrentPair%Index == PreviousPair%Index ) then
            MainParticle(imin)%FirstWall => Aux
          else
            PreviousPair%Next => Aux
          end if

        end if

      else

        MainParticle(imin)%FirstWall => Aux

      end if

    end if

  end subroutine UpdateWallContact

!===============================================================================

  subroutine RemoveWallContact

    if ( PreviousPair%Index == iwall ) then
      MainParticle(imin)%FirstWall => CurrentPair%Next
    else
      PreviousPair%Next => CurrentPair%Next
    end if

    deallocate( CurrentPair )

  end subroutine RemoveWallContact

!===============================================================================

  subroutine DestroyContactList

    integer(ib) :: i, NP
    type(TPartner), pointer :: Current, Aux

    NP = size(MainParticle)
    do i = 1_ib, NP
      Current => MainParticle(i)%FirstPartner
      do while (associated(Current))
        Aux => Current
        Current => Current%Next
        deallocate(Aux)
      end do
      Current => MainParticle(i)%FirstWall
      do while (associated(Current))
        Aux => Current
        Current => Current%Next
        deallocate(Aux)
      end do
    end do ! i
    deallocate(MainParticle)

  end subroutine DestroyContactList

!===============================================================================

  subroutine RemoveParticleFromContactList( i, NPSC, Neighbour )

    integer(ib), intent(in) :: i, NPSC, Neighbour(NPSC)

    integer(ib) :: jj, j
    type (TPartner), pointer :: Aux, Current, Previous

    ! Remove particle "i" where it appears as a secondary partner:
    do jj = 1_ib, NPSC

      j = Neighbour(jj)

      ! Point to the first partner of the particle "j" (if it exists):
      Current  => MainParticle(j)%FirstPartner

      if (associated(Current)) then

        ! Search for partner "i" or for the end of the list:
        Previous => Current
        do while ( (Current%Index < i).and.associated(Current%Next) )
          Previous => Current
          Current  => Current%Next
        end do ! while

        ! If found, remove it from the list:
        if (Current%Index == i) then
          if ( Previous%Index == i ) then
            MainParticle(j)%FirstPartner => Current%Next
          else
            Previous%Next => Current%Next
          end if
          deallocate( Current )
        end if

      end if

    end do

    ! Remove all the secondary partners of particle "i":
    Current => MainParticle(i)%FirstPartner
    do while (associated(Current))
      Aux => Current
      Current => Current%Next
      deallocate( Aux )
    end do
    Current => MainParticle(i)%FirstWall
    do while (associated(Current))
      Aux => Current
      Current => Current%Next
      deallocate( Aux )
    end do
    MainParticle(i)%FirstPartner => null()
    MainParticle(i)%FirstWall => null()

  end subroutine RemoveParticleFromContactList

!===============================================================================

end module MContacts