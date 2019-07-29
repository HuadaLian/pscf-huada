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
!****m  scf/grid_mod
! PURPOSE
!    Declare dimensions ngrid(3) of fft grid
!    Declare global data structures that are defined on an fft grid
!    Routines to allocate, deallocate and manipulate grid data
!    Functions G_to_fft, G_to_bz, norm to manipulate wavevectors
! SOURCE
!-----------------------------------------------------------------------
module grid_mod
   use const_mod
   use io_mod
   use chemistry_mod
   use kronrod_mod 
   !use string_mod
   implicit none

   private

   ! Public data structures
   public :: ngrid          ! dimensions ngrid(:)=(N1,N2,N3) of grid
   public :: chain_step     ! step size chain_step(:Nblock, :Nchain) 
   public :: lbar           ! = max level of spherical harmonics l
   public :: N_sph          ! # of spherical harmonics 
   public :: angularf_grid  ! position and weights of orientation dimension (forward)                
   public :: angularr_grid  ! position and weights of orientation dimension (reversed)              
   public :: Rjj            ! Auxilliary maxtrix for wormlike block 
   public :: Yjf            ! Auxilliary maxtrix for wormlike block 
   public :: Yjr            ! Auxilliary maxtrix for wormlike block 
   public :: Yjfinv            ! Auxilliary maxtrix for wormlike block 
   public :: Yjrinv            ! Auxilliary maxtrix for wormlike block 

   public :: rho_grid       ! rho on r-grid   (r,N_monomer)
   public :: omega_grid     ! omega on r-grid (r,N_monomer)
   public :: ksq_grid       ! k**2 on k-grid 

   ! Public procedures
   public :: input_grid
   public :: input_chainstep
   public :: output_chainstep
   public :: allocate_grid
   public :: deallocate_grid
   public :: make_ksq
   public :: max_Gabs
   public :: G_to_fft
   public :: G_to_bz
   public :: Greal_to_bz
   public :: norm

   integer                  :: ngrid(3)
   real(long), ALLOCATABLE  :: chain_step(:,:) 
   integer                  :: lbar 
   integer                  :: N_sph 
   real(long), ALLOCATABLE  :: angularf_grid(:,:)     ! angularf_grid(1:3(theta,phi,weight),2n+1=(lbar+1)**2)
   real(long), ALLOCATABLE  :: angularr_grid(:,:)     ! angularr_grid(1:3(theta,phi,weight),2n+1=(lbar+1)**2)

   real(long), ALLOCATABLE  :: Rjj(:,:,:)            ! Rjj(1:3,0:N_sph-1,0:N_sph-1) 
   real(long), ALLOCATABLE  :: Yjf(:,:)            ! Yj(0:N_sph-1,0:N_sph-1) 
   real(long), ALLOCATABLE  :: Yjr(:,:)            ! Yj(0:N_sph-1,0:N_sph-1) 
   real(long), ALLOCATABLE  :: Yjfinv(:,:)            ! Yj(0:N_sph-1,0:N_sph-1) 
   real(long), ALLOCATABLE  :: Yjrinv(:,:)            ! Yj(0:N_sph-1,0:N_sph-1) 

   real(long), ALLOCATABLE  :: rho_grid(:,:,:,:)
   real(long), ALLOCATABLE  :: omega_grid(:,:,:,:)
   real(long), ALLOCATABLE  :: ksq_grid(:,:,:)
   !***

   !---------------------------------------------------------------
   !****v grid_mod/ngrid
   ! VARIABLE
   !   integer ngrid(3) :: Number of grid points in each direction
   !***
   !---------------------------------------------------------------
   !****v grid_mod/rho_grid
   ! VARIABLE
   !   real(long) rho_grid(0:,0:,0:,N_monomer)
   !              density of each monomer type
   !***
   !---------------------------------------------------------------
   !****v grid_mod/omega_grid
   ! VARIABLE
   !   real(long) omega_grid(0:,0:,0:,N_monomer)
   !              potential for each monomer type
   !***
   !---------------------------------------------------------------
   !****v grid_mod/ksq_grid
   ! VARIABLE
   !   real(long) ksq_grid(0:,0:,0:)
   !              wave vector magnitude squared, for FFT use
   !***
   !---------------------------------------------------------------

contains

   !--------------------------------------------------------------   
   !****p grid_mod/input_grid
   ! SUBROUTINE
   !   input_grid
   ! PURPOSE
   !   Input integer array ngrid(1:dim) from input script
   ! SOURCE
   !---------------------------------------------------------------
   subroutine input_grid
   !***
   use io_mod, only : input

   call input(ngrid(:),dim,f='A')
   if (dim < 3) then
      ngrid(dim+1:3) = 1
   end if

   call input(lbar,'lbar')

   N_sph = (lbar+1)**2
   call set_angular_grid(lbar/2)

   end subroutine input_grid
   !==============================================================

   ! create the disretization of orientation dimension (theta,phi)  
   subroutine set_angular_grid(n)
      implicit none 
      integer, intent(IN) :: n
      real(long) :: eps=0.000001D+00
      real(long) :: x(n+1)    ! Kronrod abscissas half of(-1,1) 
      real(long) :: w1(n+1)   ! Gauss weights 
      real(long) :: w2(n+1)   ! Kronrod weights half of (-1,1)

      integer :: j            ! looping variable 
      integer :: k, k2, l, l2 ! looping variable and local index 
      integer :: s                ! sign +1 or -1 
      integer :: error 
      real(long) :: twopi     ! const 2pi
      !*** 

      twopi = 4.0_long*acos(0.0_long) 

      if (lbar - 2*n .GT. 0) stop 'lbar had better to be an even number.'

      ALLOCATE(angularf_grid(1:3, 0:N_sph-1), STAT=error) 
      if (error/=0) stop 'Allocating angularf_grid is in trouble.'
      ALLOCATE(angularr_grid(1:3, 0:N_sph-1), STAT=error) 
      if (error/=0) stop 'Allocating angularr_grid is in trouble.'

      ! calculate the kronrod abscissas on (-1,1), only the (1,0) part is stored. 
      ! x  => abscissas 
      ! w1 => weights 
      call kronrod(n, eps, x, w1, w2) 

      ! translate the abscissas on (-1,1) to (0,pi) and (0,2pi) respectively 
      ! the weights should be multiplied by pi/2.0 and 2*pi/2.0 respectively 
      ! weights also contains the sin(theta) part for the S2 integral 
      angularf_grid = 1.0_long 
      angularr_grid = 1.0_long 
      j = 0  
      do k= 1, 2*n+1 
         do l = 1, 2*n+1 
            s = 1 
            if (k .LE. n+1) then 
               s = (-1) * s 
               k2 = k 
            elseif (k .GT. n+1) then
               s =  1 * s 
               k2 = 2*(n+1) - k 
            endif 
            angularf_grid(1,j) = (s * x(k2) + 1.0_long)/2.0_long*(twopi/2.0_long)
            angularr_grid(1,j) = (twopi/2.0_long)-(s * x(k2) + 1.0_long)/2.0_long * (twopi/2.0_long)

            angularf_grid(3,j) = angularf_grid(3,j)*sin(angularf_grid(1,j))*w1(k2) /2.0_long * (twopi/2.0_long)
            angularr_grid(3,j) = angularr_grid(3,j)*sin(angularr_grid(1,j))*w1(k2) /2.0_long * (twopi/2.0_long)

            s = 1 
            if (l .LE. n+1) then 
               s = (-1) * s 
               l2 = l 
            elseif (l .GT. n+1) then
               s =  1 * s 
               l2 = 2*(n+1) - l 
            endif 
            angularf_grid(2,j) = (s * x(l2) + 1.0_long)/2.0_long * twopi
            angularr_grid(2,j) = (twopi/2.0_long)+(s * x(l2) + 1.0_long)/2.0_long * (twopi)
            if (angularr_grid(2,j) > (4.0_long*acos(0.0_long))) then 
               angularr_grid(2,j) = angularr_grid(2,j) - twopi
            endif

            angularf_grid(3,j) = angularf_grid(3,j)*w1(l2) /2.0_long * (twopi) 
            angularr_grid(3,j) = angularr_grid(3,j)*w1(l2) /2.0_long * (twopi) 

            j = j + 1 
         enddo 
      enddo 
       
   end subroutine set_angular_grid



   !--------------------------------------------------------------   
   !****p grid_mod/input_chainstep
   ! SUBROUTINE
   !   input_chainstep
   ! PURPOSE
   !   Input real array chain_step( :,: ) from input script
   ! SOURCE
   !---------------------------------------------------------------
   subroutine input_chainstep(i_unit,fmt_flag)
   !***

   integer,     intent(IN) :: i_unit      ! input file unit #
   character(1),intent(IN) :: fmt_flag    ! F = formatted, U = unformatted
   !***

   integer           :: error
   integer           :: j                 ! looping variables
   character(len = 100) :: comment_line   

   select case(fmt_flag)
   case('F')

      ! Set io defaults 
      call set_io_units(i=i_unit,o=6)
      call set_echo(1)                  !echo inputs
      call set_com_style('A','A','A')   ! comment above data 
      call set_com_use('R')             ! replace comment in echo 

      ALLOCATE(chain_step(1:N_blk_max,1:N_chain),STAT=error)
      if(error /= 0) stop "chain_step allocation error!"
   
      read(i_unit,*) comment_line 
      call output(trim(comment_line),f='N',j='L')
      do j = 1, N_chain
         call input(chain_step(:,j),N_block(j),f='N')
      enddo

   case('U')
      ALLOCATE(chain_step(1:N_blk_max,1:N_chain),STAT=error)
      if(error /= 0) stop "chain_step allocation error!"
   
      do j=1, N_chain
         read(i_unit) chain_step(1:N_block(j),j)
      enddo
   case default
      write(6,*) 'Error: Illegal fmt_flag in input_chainstep'
      stop
   end select
       
   end subroutine input_chainstep
   !==============================================================


   !--------------------------------------------------------------------
   !****p* grid_mod/output_chainstep
   ! SUBROUTINE
   !    output_chainstep(o_unit,fmt_flag)
   ! PURPOSE
   !    Output information about linear block copolymers and homopolymers
   !    Allocate related arrays
   ! SOURCE
   !--------------------------------------------------------------------
   subroutine output_chainstep(o_unit,fmt_flag)
   use io_mod
   use chemistry_mod, only : N_block, N_chain

   integer,      intent(IN) :: o_unit   ! output file unit #
   character(1), intent(IN) :: fmt_flag ! F = formatted, U = unformatted
   !***

   integer              :: i, j         ! loop indices

   select case(fmt_flag)
   case('F')

      ! Set io defaults
      call set_io_units(o=o_unit)
      call set_com_style('A','A','A')   ! comment above data

      call output("chain_step",f='N',j='L')
      do j = 1,N_chain
         call output(chain_step(:,j),N_block(j),f='N')
      enddo 
   case('U')

      ! Polymer properties
      if (N_chain > 0) then
         do j=1, N_chain
            write(o_unit) chain_step(1:N_block(j),j)
         enddo
      endif

   case default

      write(6,*) 'Error: Illegal fmt_flag in output_chainstep'
      stop

   end select
   
   end subroutine output_chainstep
   !===================================================================



   !--------------------------------------------------------------   
   !****p grid_mod/allocate_grid
   ! SUBROUTINE
   !   allocate_grid
   ! PURPOSE
   !   Allocate memory needed by the grid module data structures
   !   rho_grid, omega_grid, ksq_grid
   ! SOURCE
   !---------------------------------------------------------------
   subroutine allocate_grid(N_monomer)
   implicit none
   integer, intent(IN) :: N_monomer
   !***

   integer  :: nx, ny, nz
   integer  :: error

   nx=ngrid(1)-1
   ny=ngrid(2)-1
   nz=ngrid(3)-1

   ALLOCATE(rho_grid(0:nx,0:ny,0:nz,N_monomer),STAT=error)
   if (error /= 0) stop "rho_grid allocation error!"

   ALLOCATE(omega_grid(0:nx,0:ny,0:nz,N_monomer),STAT=error)
   if (error /= 0) stop "rho_grid allocation error!"

   ALLOCATE(ksq_grid(0:(nx+1)/2,0:ny,0:nz),STAT=error)
   if (error /= 0) stop "rho_grid allocation error!"

   end subroutine allocate_grid
   !==============================================================


   !--------------------------------------------------------------   
   !****p grid_mod/make_ksq
   ! SUBROUTINE
   !   make_ksq
   ! PURPOSE
   !   Calculate values of k**2 and store in ksq_grid(:,:,:)
   ! COMMENT
   !   calculates norm(G_to_bz(G)) for each wavevector G
   ! SOURCE
   !--------------------------------------------------------------   
   subroutine make_ksq(G_basis)
   implicit none
   real(long)  :: G_basis(:,:)
   !***
   integer :: G(3), Gbz(3)
   integer :: i1, i2, i3

   do i1=0, ngrid(1)/2
      G(1)=i1
      do i2=0, ngrid(2)-1
         G(2)=i2
         do i3=0, ngrid(3)-1
            G(3)=i3
            Gbz=G_to_bz(G)
            ksq_grid(i1,i2,i3)=norm(Gbz,G_basis)
         enddo
      enddo
   enddo

   end subroutine make_ksq
   !==============================================================


   !--------------------------------------------------------------   
   !****p  grid_mod/deallocate_grid
   ! SUBROUTINE
   !    destroy_grid
   ! PURPOSE
   !    deallocate(return) the memory used by grid_data
   ! SOURCE
   !--------------------------------------------------------------   
   subroutine deallocate_grid
   implicit none
   !***
   integer     :: error      ! deallocation error index

   DEALLOCATE(ksq_grid,STAT=error)
   if (error /= 0) stop "rho_grid deallocation error!"

   DEALLOCATE(rho_grid,STAT=error)
   if (error /= 0) stop "rho_grid deallocation error!"

   DEALLOCATE(omega_grid,STAT=error)
   if (error /= 0) stop "rho_grid deallocation error!"

   end subroutine deallocate_grid
   !==============================================================


   !--------------------------------------------------------------   
   !****p  grid_mod/G_to_fft
   ! SUBROUTINE
   !   G_to_fft(G)
   ! PURPOSE
   !   Shift any reciprocal wave vector G = [G(1),..,G(dim)] to 
   !   an equivalent wavevector with indices 0 <= G(i) < ngrid(i) 
   ! SOURCE
   !--------------------------------------------------------------   
   function G_to_fft(G)
   implicit none
   integer             :: G_to_fft(3)
   integer, intent(IN) :: G(:)
   !***
   integer             :: i
   G_to_fft=0
   do i=1,dim
      G_to_fft(i)=modulo( G(i), ngrid(i) )
   end do
   end function G_to_fft
   !==============================================================


   !------------------------------------------------------------------
   !****p  grid_mod/G_to_bz
   ! SUBROUTINE
   !    G_to_bz
   ! PURPOSE
   !   Shift any reciprocal wave vector of the crystal to the first 
   !   Brillouin zone of the grid unit cell, i.e., to the equivalent
   !   wavevector with the smallest absolute magnitude. 
   ! ARGUMENTS
   !   integer G(3) - array containing integer indices of wavevector
   ! SOURCE
   !------------------------------------------------------------------
   function G_to_bz(G)
   use unit_cell_mod, only : G_basis
   implicit none
   integer             :: G_to_bz(3)
   integer, intent(IN) :: G(:)         ! wave vector index
   !***

   ! local variables
   integer, parameter :: mup=1, mdn=-1
   integer            :: G_min(3),G_try(3)
   integer            :: i1,i2,i3
   real(long)         :: Gsq_min

   G_try   = 0
   G_min   = 0
   Gsq_min = 1.0E+10
   select case(dim)
   case(1)
      do i1=mup,mdn,-1
         G_try(1)=G(1)+i1*ngrid(1)
         call choose
      enddo
   case(2)
      do i1=mup,mdn,-1
         G_try(1)=G(1)+i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=G(2)+i2*ngrid(2) 
            call choose
         enddo
      enddo
   case(3)
      do i1=mup,mdn,-1
         G_try(1)=G(1)+i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=G(2)+i2*ngrid(2)
            do i3=mup,mdn,-1
               G_try(3)=G(3)+i3*ngrid(3)
               call choose
            enddo
         enddo
      enddo
   end select
   G_to_bz=G_min

   contains ! internal subroutine choose

     !-----------------------------------------
      subroutine choose
      real(long), parameter :: delta=1.0E-8
      real(long)            :: Gsq
      Gsq=norm(G_try,G_basis)
      if (Gsq < Gsq_min - delta) then
         G_min(1:dim) = G_try(1:dim)
         Gsq_min      = Gsq 
      endif
      end subroutine choose
     !-----------------------------------------

   end function G_to_bz
   !==============================================================

   function Greal_to_bz(Greal)
   use unit_cell_mod, only : G_basis
   implicit none
   real(long)                :: Greal_to_bz(3)
   real(long),intent(IN)     :: Greal(:)         ! wave vector index
!   real(long),intent(IN)  :: G_basis(:,:)

!  local variables
   integer, parameter :: mup=1, mdn=-1
   real(long)            :: G_min(3),G_try(3)
   integer            :: i1,i2,i3
   real(long)         :: Gsq_min

   G_try   = 0
   G_min   = 0
   Gsq_min = 1.0E+10

   select case(dim)
   case(1)
      do i1=mup,mdn,-1
         G_try(1)=Greal(1) + i1 * ngrid(1)
         call choose
      enddo
   case(2)
      do i1=mup,mdn,-1
         G_try(1)=Greal(1)+ i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=Greal(2)+i2*ngrid(2) 
            call choose
         enddo
      enddo
   case(3)
      do i1=mup,mdn,-1
         G_try(1)=Greal(1)+i1*ngrid(1)
         do i2=mup,mdn,-1
            G_try(2)=Greal(2)+i2*ngrid(2)
            do i3=mup,mdn,-1
               G_try(3)=Greal(3)+i3*ngrid(3)
               call choose
            enddo
         enddo
      enddo
   end select
   Greal_to_bz=G_min

   contains      ! internal subroutine
     !-----------------------------------------
      subroutine choose
        use group_mod,only: operator(.dot.)
      real(long), parameter :: delta=1.0E-8
      real(long)            :: Gsq
      real(long)            :: v(3)
      integer               ::i
      
      v = 0.0_long
      do i = 1,dim
         v(i) = G_try(:) .dot. G_basis(:,i)
      end do
      Gsq = (v .dot. v)
!      Gsq=norm(G_try,G_basis)
      if (Gsq < Gsq_min - delta) then
         G_min(1:dim) = G_try(1:dim)
         Gsq_min      = Gsq 
      endif
      end subroutine choose
     !-----------------------------------------

    end function Greal_to_bz
! ==============================================================




   !--------------------------------------------------------------   
   !****p  grid_mod/max_Gabs
   ! SUBROUTINE
   !   max_Gabs
   ! PURPOSE
   !   calculate the max magniture of vectors in FFT grid
   ! SOURCE
   !--------------------------------------------------------------   
   function max_Gabs(G_basis)
   implicit none
   real(long)            :: max_Gabs
   real(long),intent(IN) :: G_basis(:,:)
   !***
   real(long)            :: length,tmp
   real(long),parameter  :: delta=1.0E-10
   integer               :: G(3),Gbz(3)
   integer               :: i1,i2,i3
   length=1.0E-15
   do i3=0,ngrid(3)-1
      G(3)=i3
      do i2=0,ngrid(2)-1
         G(2)=i2
         do i1=0,ngrid(1)-1
            G(1)=i1
            Gbz=G_to_bz(G)
            tmp=norm(Gbz,G_basis)
            if(tmp > length+delta) length=tmp
         end do
      end do
   end do
   max_Gabs=sqrt(length)
   end function max_Gabs
   !==============================================================


   !--------------------------------------------------------------
   !****p  grid_mod/norm
   ! SUBROUTINE
   !   norm(G,G_basis)
   ! PURPOSE
   !   Calculate the squared magnitude of vector G
   ! SOURCE
   !--------------------------------------------------------------
   function norm(G,G_basis)
   use group_mod,     only : operator(.dot.)
   implicit none
   real(long)             :: norm
   integer,intent(IN)     :: G(:)
   real(long),intent(IN)  :: G_basis(:,:)
   !***
   real(long)             :: Gvec(3)

   Gvec=G.dot.G_basis
   norm=Gvec.dot.Gvec
   end function norm
   !=================================================================


end module grid_mod
