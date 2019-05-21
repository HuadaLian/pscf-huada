!-----------------------------------------------------------------------
! PSCF - Polymer Self-Consistent Field Theory
! Copyright (2002-2016) Regents of the University of Minnesota
! contact: David Morse, morse012@umn.edu
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation. A copy of this license is included in
! the LICENSE file in the top-level PSCF directory. 
!****m scf/chain_mod ---------------------------------------------------
! MODULE
!   chain_mod
! PURPOSE
!   Define derived type chain_grid_type, and allocate arrays of this type.
!   Subroutines to make and destroy objects of chain_type.
! SOURCE
!-----------------------------------------------------------------------
module chain_mod
   use const_mod
   use lebedev_mod
   use chemistry_mod
   use grid_mod
   use fft_mod, only : fft_plan
   implicit none

   private
   public :: chain_grid_type
   public :: null_chain_grid
   public :: make_chain_grid
   public :: destroy_chain_grid
   !***

   public :: allocate_q, allocate_qw    

   !****t chain_mod/chain_grid_type  -----------------------------------
   !  TYPE 
   !     chain_grid_type
   !  PURPOSE
   !     Data structures defining discretization of s for a chain.
   !     Pointers to qf(r,s) and qr(r,s) functions for a chain
   !     Pointers to del_qf and del_qr functions for perturbation theory
   !  SOURCE
   !-------------------------------------------------------------------
   type chain_grid_type
      logical                 :: block_exist(1:2)   ! 1:exist_q, 2:exist_qw
      integer,    pointer     :: block_bgn_lst(:,:)   ! 1st element and last elemetn of block
      real(long), pointer     :: block_ds(:)        ! step size for block
      real(long), pointer     :: qf(:,:,:,:)        ! function qf(x,y,z,s) 
      real(long), pointer     :: qr(:,:,:,:)        ! function qr(x,y,z,s) 
      real(long), pointer     :: rho(:,:,:,:)       ! rho(x,y,z,s)
      real(long), pointer     :: qwf(:,:,:,:,:)   ! qwf(x,y,z,g,s)
      real(long), pointer     :: qwr(:,:,:,:,:)   ! qwr(x,y,z,g,s)
      real(long), pointer     :: qwj(:,:,:,:,:)     ! function qwj(x,y,z,j,s)
      type(fft_plan)          :: plan               ! fft plan, see fft_mod
      real(long)              :: bigQ               ! chain partition func.
      complex(long), pointer  :: del_qf(:,:,:,:)    ! perturbation in qf
      complex(long), pointer  :: del_qr(:,:,:,:)    ! perturbation in qr
      complex(long)           :: delQ               ! perturbation in bigQ
   end type
   !***

   logical       :: allocate_q = .FALSE.    ! By default 
   logical       :: allocate_qw= .FALSE.    ! By default  

contains

   !------------------------------------------------------------------
   !****p chain_mod/make_chain_grid
   ! PURPOSE
   !    Nullify the pointers to ensure they have the "dissociated"
   !    status.
   ! SOURCE
   !------------------------------------------------------------------
   subroutine null_chain_grid(chain)
   implicit none

   type(chain_grid_type), intent(INOUT)  :: chain
   !***

   nullify( chain%block_bgn_lst, &
            chain%block_ds,  &
            chain%qf,        &
            chain%qr,        &
            chain%qwf,       &
            chain%qwr,       &
            chain%qwj,       &
            chain%rho,       &
            chain%del_qf,    &
            chain%del_qr)

   end subroutine null_chain_grid
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !****p chain_mod/make_chain_grid
   ! SUBROUTINE
   !    make_chains(chain,plan,nblk,blk_length,ds)
   ! PURPOSE
   !    allocate memory for a single chain_grid_type variable
   !    initiate chain(:)%block_bgn_lst, chain(:)%block_ds
   ! ARGUMENTS
   !    plan       = grid dimensions and FFT plan, see fft_mod
   !    nblk       = # of blocks of a single chain
   !    blk_length = block lengths
   !    blk_type   = block types 
   !    ds         = segment length used to discretize the block
   ! COMMENTS
   !    The # of segments for each block need to be even, because
   !    of the use of Simpson's rule in density/stress calculation
   ! SOURCE
   !------------------------------------------------------------------
   subroutine &
       make_chain_grid(chain,plan,nblk,blk_length,blk_type,ds,perturb,order)
   implicit none
   type(chain_grid_type), intent(INOUT)  :: chain
   type(fft_plan),        intent(IN)     :: plan
   integer,               intent(IN)     :: nblk
   real(long),            intent(IN)     :: blk_length(:)
   character(20),         intent(IN)     :: blk_type(:)
   real(long),            intent(IN)     :: ds(:)
   logical,optional,      intent(IN)     :: perturb
   integer,optional,      intent(IN)     :: order
   !***
   real(long)            :: ds_iblk            ! step size of the ith block according to type              
   integer               :: iblk               ! index to block
   integer               :: bgnsGaus           ! index of beginers for Gaussian block
   integer               :: bgnsWorm           ! index of beginers for Wormlike block
   integer               :: nx,ny,nz,nt,i   ! loop indices

   integer               :: error              ! allocation error-message

   chain%block_exist = .FALSE. 
   do i=1, nblk
      if (blk_type(i)=='Gaussian') then
         chain%block_exist(1) = .TRUE.
         allocate_q =.TRUE. 
      elseif (blk_type(i)=='Wormlike') then
         chain%block_exist(2) = .TRUE.
         allocate_qw = .TRUE. 
      else
         stop 'Invalid type of block'
      endif
   enddo

   nx=plan%n(1)-1
   ny=plan%n(2)-1
   nz=plan%n(3)-1
   if (chain%block_exist(2)) then
      nt = lebedev_order -1
   endif
   chain%plan=plan

   if ( .not. associated(chain%block_bgn_lst) ) then
      allocate(chain%block_bgn_lst(1:2,nblk),STAT=error)
      if (error /= 0) stop 'chain%block_bgn_lst allocation error'
   end if

   if ( .not. (associated(chain%block_ds)) ) then
      allocate(chain%block_ds(nblk),STAT=error)
      if (error /= 0) stop 'chain%block_ds allocation error'
   end if

   if ( .not. associated(chain%rho) ) then
      allocate(chain%rho(0:nx,0:ny,0:nz,nblk),STAT=error)
      if (error /= 0) stop 'chain%rho allocation error!'
   end if

   bgnsGaus=1
   bgnsWorm=1

   do i=1, nblk    ! loop over blocks

      ds_iblk = ds(i)
      iblk=int(blk_length(i)/ds_iblk/2.0_long+0.5_long)

      if (iblk == 0) then
         iblk = 1
         chain%block_ds(i)=blk_length(i)/2.0_long
      else
         chain%block_ds(i)=blk_length(i)/dble(iblk)/2.0_long
      endif

      if ( present(order) ) then
        chain%block_ds(i) = chain%block_ds(i) / 2.0_long**order
        iblk = iblk * 2**order
      end if

      ! Notice: By convention, the end of each block is always stored
      !         in the first element of next block. 
      ! propagator for gaus and worm are saved in two array. 
      select case(blk_type(i))
      case ("Gaussian")
         chain%block_bgn_lst(1,i) = bgnsGaus
         bgnsGaus = bgnsGaus + iblk * 2
         chain%block_bgn_lst(2,i) = bgnsGaus-1
      case ('Wormlike')
         chain%block_bgn_lst(1,i) = bgnsWorm
         bgnsWorm = bgnsWorm + iblk * 2
         chain%block_bgn_lst(2,i) = bgnsWorm-1
      case default
         write(6,*) 'Error: Invalid block_type'
         exit 
      end select
   end do

   if ( chain%block_exist(1) ) then
      allocate(chain%qf(0:nx,0:ny,0:nz,bgnsGaus),STAT=error)
      if (error /= 0) stop "chain%qf allocation error!"

      allocate(chain%qr(0:nx,0:ny,0:nz,bgnsGaus),STAT=error)
      if (error /= 0) stop "chain%qr allocation error!"
   end if

   if ( chain%block_exist(2) ) then
      allocate(chain%qwf(0:nx,0:ny,0:nz,0:nt,bgnsWorm),STAT=error)
      if (error /= 0) stop "chain%qwf allocation error!"

      allocate(chain%qwr(0:nx,0:ny,0:nz,0:nt,bgnsWorm),STAT=error)
      if (error /= 0) stop "chain%qwr allocation error!"

      allocate(chain%qwj(0:nx,0:ny,0:nz,(lbar+1)**2,bgnsWorm),STAT=error)
      if (error /= 0) stop "chain%qwr allocation error!"
   end if

   if ( (present(perturb)) .and. perturb ) then
      allocate(chain%del_qf(0:nx,0:ny,0:nz,bgnsGaus),STAT=error)
      if (error /= 0) stop "chain%qf allocation error!"
      
      allocate(chain%del_qr(0:nx,0:ny,0:nz,bgnsGaus),STAT=error)
      if (error /= 0) stop "chain%qr allocation error!"
    end if

   end subroutine make_chain_grid
   !==================================================================


   !-------------------------------------------------------------------
   !****p chain_mod/destroy_chain_grid
   ! SUBROUTINE
   !    destroy_chain_grid(chain)
   ! PURPOSE
   !    Deallocate memories use by chain%...
   ! ARGUMENTS
   !    chain = the chain_grid_type to be deallocated
   ! SOURCE
   !-------------------------------------------------------------------
   subroutine destroy_chain_grid(chain)
   implicit none
   type(chain_grid_type) :: chain
   !***
   integer  :: error     ! deallocation error msg

   if ( associated(chain%block_bgn_lst) ) then
      deallocate(chain%block_bgn_lst,STAT=error)
      if (error /= 0) stop "chain%block_bgn_lst deallocation error!"
   endif

   if ( associated(chain%block_ds) ) then
      deallocate(chain%block_ds,STAT=error)
      if (error /= 0) stop "chain%block_ds deallocation error!"
   endif

   if ( associated(chain%rho) ) then
      deallocate(chain%rho,STAT=error)
      if (error /= 0) stop "chain%rho deallocation error!"
   endif

   if ( associated(chain%qf) ) then
      deallocate(chain%qf,STAT=error)
      if (error /= 0) stop "chain%qf deallocation error!"
   endif

   if ( associated(chain%qr) ) then
      deallocate(chain%qr,STAT=error)
      if (error /= 0) stop "chain%qr deallocation error!"
   endif

   if ( associated(chain%qwf) ) then
      deallocate(chain%qf,STAT=error)
      if (error /= 0) stop "chain%qf deallocation error!"
   endif

   if ( associated(chain%qwr) ) then
      deallocate(chain%qr,STAT=error)
      if (error /= 0) stop "chain%qr deallocation error!"
   endif

   if ( associated(chain%qwj) ) then
      deallocate(chain%qr,STAT=error)
      if (error /= 0) stop "chain%qr deallocation error!"
   endif


   if ( associated(chain%del_qf) ) then
      deallocate(chain%del_qf,STAT=error)
      if (error /= 0) stop "chain%qf deallocation error!"
   endif

   if ( associated(chain%del_qr) ) then
      deallocate(chain%del_qr,STAT=error)
      if (error /= 0) stop "chain%qr deallocation error!"
   endif

   end subroutine destroy_chain_grid
   !====================================================================

end module chain_mod
