!-----------------------------------------------------------------------
! PSCF - Polymer Self-Consistent Field Theory
! Copyright (2002-2016) Regents of the University of Minnesota
! contact: David Morse, morse012@umn.edu
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation. A copy of this license is included in
! the LICENSE file in the top-level PSCF directory. 
!-----------------------------------------------------------------------
!****m scf/fft_mod
! PURPOSE
!    Derived type fft_plan
!    Wrapper subroutines for fftw-3
! COMMENTS
!    Consider use even/odd fftw, and r-to-r fftw in the future
! SOURCE
!-----------------------------------------------------------------------
module fft_mod
   use const_mod
   use omp_lib
   implicit none

   PRIVATE 
   PUBLIC :: fft_plan
   PUBLIC :: fft_plan_many  

   !Generic interface
   PUBLIC :: create_fft_plan    ! initialize an fft_plan
   PUBLIC :: create_fft_plan_many

   !Generic interface 
   PUBLIC :: fftc               ! complex FFT for 1, 2, or 3D

   !Generic interface
   PUBLIC :: fft                ! Forward FFT for 1, 2, or 3D
   PUBLIC :: fft_many           ! Many Forward FFT for 1, 2, or 3D

   !Generic interface
   PUBLIC :: ifft               ! Inverse FFT for 1, 2, or 3D
   PUBLIC :: ifft_many          ! Many Inverse FFT for 1, 2, or 3D

   !***
   ! Parameters required by fftw3 supplied in fftw3.f
   integer, parameter :: FFTW_ESTIMATE=64
   integer, parameter :: FFTW_FORWARD=-1 ,FFTW_BACKWARD=1

   !-------------------------------------------------------------------
   !****t fft_mod/fft_plan
   ! TYPE
   !    fft_plan 
   ! PURPOSE
   !    Contains grid dimensions for FFT grid and integer pointers to
   !    the "plan" structures used by the FFTW package
   ! SOURCE
   !-------------------------------------------------------------------
   type fft_plan
      integer    ::  n(3)   ! grid dimensions, 0 for unused dimensions
      integer*8  ::  f      ! fftw plan object for forward transform
      integer*8  ::  r      ! fftw plan object for inverse transform
   end type fft_plan
   !***

   !-------------------------------------------------------------------
   !****t fft_mod/fft_plan_many
   ! TYPE
   !    fft_plan_many 
   ! PURPOSE
   !    Contains grid dimensions for FFT grid and integer pointers to
   !    the "plan" structures used by the FFTW package
   ! SOURCE
   !-------------------------------------------------------------------
   type fft_plan_many
      integer    ::  n(3) 
      integer    ::  howmany_l
      integer    ::  howmany_m 
      integer*8,allocatable  ::  f(:,:)      ! fftw plan object for forward transform
      integer*8,allocatable  ::  r(:,:)      ! fftw plan object for inverse transform
   end type fft_plan_many
   !***

   interface create_fft_plan 
      module procedure create_fft_plan
      module procedure create_fft_plan_many
   end interface
   interface fft 
      module procedure fft 
      module procedure fft_many
   end interface
   interface ifft 
      module procedure ifft 
      module procedure ifft_many 
   end interface

contains

   !-------------------------------------------------------------------
   !****p fft_mod/create_fft_plan
   ! SUBROUTINE
   !    create_fft_plan
   ! PURPOSE
   !    Creates an fft_plan object for grids with dimensions 
   !    ngrid(1),..,ngrid(dim)
   !-------------------------------------------------------------------
   subroutine create_fft_plan(ngrid,plan,fft_c2c)
   integer,intent(IN)             :: ngrid(3) ! dimensions of grid
   type(fft_plan),intent(OUT)     :: plan
   logical, optional, intent(IN)  :: fft_c2c
   !***
   integer                        :: void 

   plan%n=ngrid
   end subroutine create_fft_plan
   !===================================================================

   !-------------------------------------------------------------------
   !****p fft_mod/create_fft_plan_many
   ! SUBROUTINE
   !    create_fft_plan_many
   ! PURPOSE
   !    Creates an fft_plan object for many transoformation of 
   !    grids with dimensions ngrid(1),..,ngrid(dim)
   !-------------------------------------------------------------------
   subroutine create_fft_plan_many(ngrid,howmany_l,howmany_m,plan_many)
   integer,intent(IN)             :: ngrid(3)    ! dimensions of grid
   integer,intent(IN)             :: howmany_l   !
   integer,intent(IN)             :: howmany_m   !

   type(fft_plan_many),intent(OUT):: plan_many
   !***
   integer                        :: error,void

   call omp_set_num_threads(omp_get_num_procs()/2) 

   print *, omp_get_num_procs()/2, 'threads are used.'
   plan_many%n       = ngrid
   plan_many%howmany_l = howmany_l
   plan_many%howmany_m = howmany_m

   ALLOCATE(plan_many%f(0:howmany_l-1,0:howmany_m-1), STAT= error) 
   if (error /= 0) stop "plan_many%f(howmany) allocation error!" 
   ALLOCATE(plan_many%r(0:howmany_l-1,0:howmany_m-1), STAT= error) 
   if (error /= 0) stop "plan_many%r(howmany) allocation error!" 

   end subroutine create_fft_plan_many
   !===================================================================


   !-------------------------------------------------------------------
   !****p fft_mod/fft
   ! SUBROUTINE
   !     fft(plan,in,out)
   ! PURPOSE
   !    Calculates forward fft of in, returns result in out.
   !    Wrapper for 1, 2, & 3 dimensional real -> complex transforms
   ! ARGUMENTS
   !    plan    - fft plan object
   !    in, out - real(long) 3D arrays 
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine fft(plan,in,out)
   type(fft_plan),intent(IN)   :: plan
   real(long), intent(IN)      :: in(0:,0:,0:)
   complex(long), intent(OUT)  :: out(0:,0:,0:)
   !***

   call dfftw_plan_dft_r2c_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
                              in,out,FFTW_ESTIMATE)
   call dfftw_execute(plan%f)
   call dfftw_destroy_plan(plan%f)
   end subroutine fft
   !===================================================================

   !-------------------------------------------------------------------
   !****p fft_mod/fft_many
   ! SUBROUTINE
   !     fft_many(plan_many,in,out)
   ! PURPOSE
   !    Calculates forward fft of in(s), returns result in out(s).
   !    Wrapper for 1, 2, & 3 dimensional real -> complex transforms
   ! ARGUMENTS
   !    plan_many    - fft plan_many object
   !    in, out      - real(long) 3D arrays 
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine fft_many(plan_many,in,out)
   type(fft_plan_many),intent(IN)   :: plan_many
   real(long),intent(in)            :: in(0:,0:,0:,0:,0:)
   complex(long), intent(OUT)       :: out(0:,0:,0:,0:,0:)
   !***
   integer :: l,m  ! looping variables


   do l = 0, plan_many%howmany_l-1
   do m = 0, plan_many%howmany_m-1
      call dfftw_plan_dft_r2c_3d(plan_many%f(l,m),plan_many%n(1), plan_many%n(2),plan_many%n(3),&
                              in(:,:,:,l,m), out(:,:,:,l,m), FFTW_ESTIMATE) 
   enddo 
   enddo

   !$OMP PARALLEL DO 
   do l = 0, plan_many%howmany_l-1
   do m = 0, plan_many%howmany_m-1 
      call dfftw_execute(plan_many%f(l,m))
   enddo
   enddo 
   !$OMP END PARALLEL DO  

   do l = 0, plan_many%howmany_l-1
   do m = 0, plan_many%howmany_m-1
      call dfftw_destroy_plan(plan_many%f(l,m))
   enddo
   enddo 

   end subroutine fft_many



   !-------------------------------------------------------------------
   !****p fft_mod/ifft
   ! SUBROUTINE
   !     ifft(plan,in,out)
   ! PURPOSE
   !    Calculates inverse fft of real array in, returns in out.
   !    Wrapper for 1, 2, & 3 dimensional complex -> real transforms
   ! ARGUMENTS
   !    plan - fft plan object
   !    in   - complex(long) 3D input array
   !    out  - real(long) 3D input array
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine ifft(plan,in,out)
   type(fft_plan),intent(IN)   :: plan
   complex(long), intent(IN)   :: in(0:,0:,0:)
   real(long), intent(OUT)     :: out(0:,0:,0:)
   !***

   call dfftw_plan_dft_c2r_3d(plan%r,plan%n(1),plan%n(2),plan%n(3),&!
                              in,out,FFTW_ESTIMATE)
   call dfftw_execute(plan%r)
   call dfftw_destroy_plan(plan%r)
   end subroutine ifft
   !===================================================================

   !-------------------------------------------------------------------
   !****p fft_mod/ifft_many
   ! SUBROUTINE
   !     ifft_many(plan_many,in,out)
   ! PURPOSE
   !    Calculates a series of inverse fft of real array ins, returns in outs.
   !    Wrapper for 1, 2, & 3 dimensional complex -> real transforms
   ! ARGUMENTS
   !    plan - fft plan object
   !    in   - complex(long) 3D input array
   !    out  - real(long) 3D input array
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------

   subroutine ifft_many(plan_many,in,out)
   type(fft_plan_many),intent(IN)   :: plan_many
   complex(long),intent(in)            :: in(0:,0:,0:,0:,0:)
   real(long), intent(OUT)       :: out(0:,0:,0:,0:,0:)
   !***
   integer :: l,m  ! looping variables

   do l = 0, plan_many%howmany_l-1
   do m = 0, plan_many%howmany_m-1
      call dfftw_plan_dft_c2r_3d(plan_many%r(l,m),plan_many%n(1), plan_many%n(2),plan_many%n(3),&
                              in(:,:,:,l,m), out(:,:,:,l,m), FFTW_ESTIMATE) 
   enddo 
   enddo 

   !$OMP PARALLEL DO 
   do l = 0, plan_many%howmany_l-1
   do m = 0, plan_many%howmany_m-1 
      call dfftw_execute(plan_many%r(l,m))
   enddo
   enddo 
   !$OMP END PARALLEL DO  

   do l = 0, plan_many%howmany_l-1
   do m = 0, plan_many%howmany_m-1
      call dfftw_destroy_plan(plan_many%r(l,m))
   enddo
   enddo 

   end subroutine ifft_many

   !-------------------------------------------------------------------
   !****p fft_mod/fftc
   ! SUBROUTINE
   !     fftc(direction,plan,in,out)
   ! PURPOSE
   !    Calculates forward fft of in, returns result in out.
   !    Wrapper for 1, 2, & 3 dimensional real -> complex transforms
   ! ARGUMENTS
   !    plan    - fft plan object
   !    in, out - real(long) 3D arrays 
   ! COMMENT
   !    in and out are dimensioned 0:ngrid(i)-1 for all i <= dim, 
   !    and 0:0 for any unused dimensions with dim < i <= 3
   !-------------------------------------------------------------------
   subroutine fftc(direction,plan,in,out)
   integer,intent(IN)          :: direction
   type(fft_plan),intent(IN)   :: plan
   complex(long), intent(IN)   :: in(0:,0:,0:)
   complex(long), intent(OUT)  :: out(0:,0:,0:)
   !***
   if (direction == 1) then
      call dfftw_plan_dft_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
             in,out,FFTW_FORWARD,FFTW_ESTIMATE)
   else
      call dfftw_plan_dft_3d(plan%f,plan%n(1),plan%n(2),plan%n(3),&!
             in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
   end if
   call dfftw_execute(plan%f)
   call dfftw_destroy_plan(plan%f)

   if (direction == +1) out = out/dcmplx( dble(plan%n(1)*plan%n(2)*plan%n(3)) , 0.0_long)
   end subroutine fftc
   !===================================================================


end module fft_mod
