!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    21/02/2020
!! Project: Umbrella
!! File:    umbrella.f90 
!! Copyright: © 2020 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!*******************************************************************************************
!!
!! Analysis tool to calculate the free energy of exciton (or other carrier)
!! transfer from the WC_POT.xyz output files provided by the modified CPMD 
!! Code according to https://doi.org/10.1039/C9CP06419B
!!
!! The program expects <nwind> files to be present in the directory, named
!! <filnam>_<i>.xyz, with <i> in [1,<nwind>].
!! Usage:
!!   umbrella <nwind> <temp> <xi_start> <xi_end> <orb. pos.> <filename> [#integration points]
!! arguments in <> are necessary, arguments in [] are optional.
!! <nwind> is the number of sampling windows, <temp> is the simulation temperature,
!! <xi_start> and <xi_end> are start and end point of the reaction coordinate,
!! <orb.pos.> is the position of the orbital of interest in the xyz files, 
!! <filnam> is the filename (see above),
!! [nintp] is the number of integration points for the free energy, default is 1000.
!!
!!*******************************************************************************************
program Umbrella
   implicit none
   !! Precision
   integer, parameter         :: DP = Selected_real_kind(15, 307)
   !! Input: #windows,position of orbital in xyz file, #integration points
   integer                    :: nwind, pospart, niint = 1000, sym = 0
   !! Input: start and end point of reaction coordinate, temperature
   real(kind=DP)              :: xstart, xend, T
   !! Input: file name
   character(20)              :: filnam
   !! Internal: number of timesteps per window
   integer, allocatable       :: nsteps(:)
   !! Internal: loop indices for windows and atoms/orbitals and integration, read error, #particles
   integer                    :: iw, k, iint, npart, ierror
   !! Internal: arrays for first and second moment of xi, variance, reference xi, force constants
   real(kind=DP), allocatable :: sumxi(:), sumxi2(:), sigma2(:), xref(:), Kon(:)
   !! Internal: name, coordinates of orbital, distance and energy,RuKu stuff, stepsize
   real(kind=DP)              :: xi, x, y, z, dist, energy, F, dF, k1, k2, k3, k4, dx
   !! Internal: character for total filename,scratch for element name
   character(30)              :: filnamiw, rchar
   !! Internal: Hartree energy
   real(kind=DP), parameter   :: EH = 27.21138505_DP
   !! 

   ! read input
   call Input( nwind, T, xstart, xend, pospart, filnam, niint, sym )
   
   ! allocate arrays
   allocate( sumxi(1:nwind), sumxi2(1:nwind), nsteps(1:nwind) )
   allocate( xref(1:nwind), Kon(1:nwind), sigma2(1:nwind) )

   ! read xi and calculate average for each window
   do iw = 1, nwind
      write(filnamiw,'(I20)') iw
      filnamiw = Adjustl( filnamiw )
      filnamiw = Trim( filnam ) // '_' // Trim( filnamiw ) // '.xyz'
      filnamiw = Trim( filnamiw )
      open( unit = 42, file = filnamiw, status = 'old', action = 'read', iostat = ierror)

      sumxi(iw) = 0.0_DP
      sumxi2(iw) = 0.0_DP
      nsteps(iw) = 0
      
      !if( ierror /= 0 ) exit

      do
        read(42,*,iostat=ierror) npart
        read(42,*,iostat=ierror) 
        if( ierror /= 0 ) exit
        ! skip COMs and unwanted orbitals
        do k = 1, pospart - 1
           read(42,*)
        end do
        read(42,*,iostat=ierror) rchar, x, y, z, dist, energy, xi
        ! skip rest
        do k = pospart + 1, npart
           read(42,*)
        end do

        ! calculate <xi> and <xi^2> and nsteps
        sumxi(iw) = sumxi(iw) + xi
        sumxi2(iw) = sumxi2(iw) + xi ** 2
        nsteps(iw) = nsteps(iw) + 1
     end do

     Kon(iw) = 2.0_DP * energy / dist ** 2
     xref(iw) = xstart + Real( iw - 1, kind = DP ) / Real( nwind - 1, kind = DP ) * ( xend - xstart )
     sumxi(iw) = sumxi(iw) / Real( nsteps(iw), kind = DP )
     sumxi2(iw) = sumxi2(iw) / Real( nsteps(iw), kind = DP )
     sigma2(iw) = sumxi2(iw) - sumxi(iw) ** 2
      
     write(*,*) iw, sumxi(iw), sigma2(iw), nsteps(iw), Kon(iw), xref(iw)

     close( unit = 42 )
  end do

  ! symmetrize free energy surface
  ! middle could be cheated, simulating twice as long might be better
  if( sym == 1 ) then
     do iw = 1, nwind / 2 + 1
        sumxi(iw) = ( sumxi(iw) + ( xend - xstart ) - sumxi(nwind + 1 - iw) ) * 0.5_DP
        sumxi(nwind+1-iw)=(xend-xstart)-sumxi(iw)

        sigma2(iw) = ( sigma2(iw) + sigma2(nwind + 1 - iw) ) * 0.5_DP
        sigma2(nwind + 1 - iw) = sigma2(iw)

        nsteps(iw) = nsteps(iw) + nsteps(nwind + 1 - iw)
        nsteps(nwind + 1 - iw) = nsteps(iw)
     end do
     do iw = 1, nwind
        write(*,*) iw, sumxi(iw), sigma2(iw), nsteps(iw), Kon(iw), xref(iw)
     end do
  end if

  ! output for free energy
  open( unit = 43, file = "Fenergy.dat", status = 'replace', action = 'write', iostat = ierror)
  
  ! integration
  F = 0.0_DP
  dx = ( xend - xstart ) / Real( niint, kind = DP )
  x = xstart

  do iint = 1, niint
  
     ! Runge-Kutta 4
     call Force( nwind, sigma2, nsteps, sumxi, Kon, xref, T, x, k1 )
     call Force( nwind, sigma2, nsteps, sumxi, Kon, xref, T, x + 0.5_DP * dx, k2 )
     call Force( nwind, sigma2, nsteps, sumxi, Kon, xref, T, x + 0.5_DP * dx, k3 )
     call Force( nwind, sigma2, nsteps, sumxi, Kon, xref, T, x + dx, k4 )
     
     dF = ( k1 + 2.0_DP * k2 + 2.0_DP * k3 + k4 ) / 6.0_DP
     ! alternative: simple
     !dF=k4

     ! update free energy
     F = F + dF * dx
     
     ! update position
     x = x + dx

     write(43,*) x, F * EH, dF * EH
  enddo

  close( unit = 43 )
 
end program Umbrella


!! This subroutine reads input from command line
!! returns with error unless 7 or 8 arguments are given
subroutine Input( nwind, T, xstart, xend, pospart, filnam, niint, sym)
   implicit none
   !! Precision
   integer, parameter           :: DP = Selected_real_kind(15, 307)
   !! Input: #windows,position of orbital in xyz file,#integration points
   integer, intent(inout)       :: nwind, pospart, niint, sym
   !! Input: start and end point of reaction coordinate, temperature
   real(kind=DP), intent(inout) :: xstart, xend, T
   !! Input: file name
   character(20), intent(inout) :: filnam
   !! Internal: character for reading input
   character(20)                :: rchar
   character(8)                 :: intchar
   !! 
   
   ! read input, need 7 to 9 parameters
   if( Iargc() == 6 ) then
      call Getarg( 1, intchar )
      read(intchar,'(I8)') nwind
      call Getarg( 2, rchar )
      read(rchar,'(F20.8)') T
      call Getarg( 3, rchar )
      read(rchar,'(F20.8)') xstart
      call Getarg( 4, rchar )
      read(rchar,'(F20.8)') xend
      call Getarg( 5, intchar )
      read(intchar,'(I8)') pospart
      call Getarg( 6, filnam )
   elseif( Iargc() == 7 ) then
      call Getarg( 1, intchar )
      read(intchar,'(I8)') nwind
      call Getarg( 2, rchar )
      read(rchar,'(F20.8)') T
      call Getarg( 3, rchar )
      read(rchar,'(F20.8)') xstart
      call Getarg( 4, rchar )
      read(rchar,'(F20.8)') xend
      call Getarg( 5, intchar )
      read(intchar,'(I8)') pospart
      call Getarg( 6, filnam )
      call Getarg( 7, intchar )
      read(intchar,'(I8)') niint
   elseif( Iargc() == 8 ) then
      call Getarg( 1, intchar )
      read(intchar,'(I8)') nwind
      call Getarg( 2, rchar )
      read(rchar,'(F20.8)') T
      call Getarg( 3, rchar )
      read(rchar,'(F20.8)') xstart
      call Getarg( 4, rchar )
      read(rchar,'(F20.8)') xend
      call Getarg( 5, intchar )
      read(intchar,'(I8)') pospart
      call Getarg( 6, filnam )
      call Getarg( 7, intchar )
      read(intchar,'(I8)') niint
      call Getarg( 8, intchar )
      read(intchar,'(I8)') sym
   else
      write(*,*) 'Usage: umbrella <#windows> <temperature> <xi_start> <xi_end> <orb. pos.> <filename> [#integration points]'
      write(*,*) 'The xyz files names have to be filename_i.xyz, with an integer 0 < i < #windows+1'
      error stop
   endif
   
end subroutine Input

! This subroutine calculates the total mean force for a given x at a given temperature
pure subroutine Force( nwind, sigma2, nsteps, sumxi, Kon, xref, T, x, k )
   implicit none
   !! Precision
   integer, parameter           :: DP = Selected_real_kind(15, 307)
   !! Input: #windows
   integer, intent(in)          :: nwind
   !1 Input: variance,average of xi, force constants, reference xi
   real(kind=DP), intent(in)    :: sigma2(1:nwind), sumxi(1:nwind), Kon(1:nwind), xref(1:nwind)
   !! Input: temperature,position
   real(kind=DP), intent(in)    :: T, x
   !! Input: #timesteps per window
   integer, intent(in)          :: nsteps(1:nwind)
   !! Output: Runge-Kutta variable k
   real(kind=DP), intent(out)   :: k
   !! Internal: variables for umbrella integration
   real(kind=DP)                :: suma
   real(kind=DP)                :: a(1:nwind), dFu(1:nwind)
   !! Internal: window index
   integer                      :: iw
   !! Internal: Pi and Boltzmann
   real(kind=DP), parameter     :: PI = 3.1415926535897932_DP
   real(kind=DP), parameter     :: KB = 3.1668114e-6_DP !8.6173303d-5
   !! 
   
   suma = 0.0_DP
   do iw = 1, nwind
      ! calculate weights and their sum
      a(iw) = 1.0_DP / ( DSqrt( 2.0_DP * PI * sigma2(iw) ) ) * DExp( -0.5_DP * ( x - sumxi(iw) ) ** 2 / sigma2(iw) ) 
      a(iw) = Real( nsteps(iw), kind = DP ) * a(iw)
      suma = suma + a(iw)

      ! calculate mean force per window
      dFu(iw) = KB * T * ( x - sumxi(iw) ) / sigma2(iw) - Kon(iw) * ( x - xref(iw) )
   enddo

   ! total mean force is average of all mean forces
   k = 0.0_DP
   do iw = 1, nwind
      k = k + a(iw) / suma * dFu(iw)
   enddo
   
end subroutine Force
