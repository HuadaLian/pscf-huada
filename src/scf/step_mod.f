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
!****m  scf/step_mod
! MODULE
!   step_mod
! PURPOSE
!   Implements pseudo-spectral algorithm for integration of the modified
!   diffusion equation. The algorithm combines the operator-splitting 
!   method of Rasmussen and Kaloskas with Richardson extrapolation to 
!   obtain an algorithm with errors of O(ds**4). 
!
!   Subroutine init_step must be called once to allocate the FFT arrays
!   used by the module. Subroutine make_propg must be called once at the 
!   beginning of each block of a block copolymer, to set variables that
!   are used throughout that block. Subroutine step implements one
!   'time' step of the integration algorithm. 
! SOURCE
!-----------------------------------------------------------------------
module step_mod
   use omp_lib 
   use SHTOOLS
   use const_mod
   use fft_mod
   use grid_mod
   use chemistry_mod 
   use chain_mod 
   use group_mod 
   implicit none

   private

   public :: init_step      ! allocate array needed by module
   public :: make_propg     ! evaluate exp(-ds*b^2*nabla/6), exp(-ds*omega/2)
   !***

   public :: step_gaussian
   public :: step_wormlike_euler
   public :: step_wormlike_bdf3
   
   ! Generic interfaces 
   public :: set_initial              ! set initial condition of q(r,s+1) or q(r,u,s+1)          
   !Auxilliary array for Gaussian chain
   real(long), allocatable    :: exp_omega1(:,:,:)
   real(long), allocatable    :: exp_omega2(:,:,:)
   real(long), allocatable    :: exp_ksq1(:,:,:)
   real(long), allocatable    :: exp_ksq2(:,:,:)

   real(long), allocatable    :: q1(:,:,:)
   real(long), allocatable    :: q2(:,:,:)
   real(long), allocatable    :: qr(:,:,:)
   complex(long), allocatable :: qk(:,:,:)

   !Auxilliary array for Wormlike chain
   real(long), allocatable :: Yj(:,:)              ! 2D - real spherical harmonics on orientation grid u 
   complex, allocatable :: Gk(:,:,:,:,:,:,:) ! 2D - vector function Gk(n1,n2,n3,j,j',i_blk,i_chain)    
   !**

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Generic Interfaces
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   !------------------------------------------------------------------
   !****p step_mod/set_initial
   ! FUNCTION
   !    set_initial(qin, qout) 
   !
   ! COMMENT
   !    The initial condition may involve the change of rank of arrays. 
   !    The rank of input qin and qout would choose appropriate method 
   !    for setting initial condition. 
   !     
   !    worm_to_gauss : qin(:,:,:,:,:), qout(:,:,:) 
   !    gauss_to_worm : qin(:,:,:), qout(:,:,:,:,:)  
   ! SOURCE
   !------------------------------------------------------------------
   interface set_initial 
      module procedure worm_to_gauss        ! from wormlike block to gaussian block       
      module procedure gauss_to_worm        ! from gaussian block to wormlike block     
   end interface 
   !***

contains

   !----------------------------------------------------------------
   !****p step_mod/init_step
   ! SUBROUTINE
   !   init_step(n)
   ! PURPOSE
   !   allocate memory to module arrays
   ! ARGUMENTS
   !   integer n(3)  - grid dimensions
   ! SOURCE
   !----------------------------------------------------------------
   subroutine init_step(n, chains)
   implicit none

   integer, intent(IN) :: n(3)
   type(chain_grid_type), intent(IN)       :: chains(:)
   !***

   integer :: error

   integer     :: N_sph                         ! # of spherical harmonics

   integer     :: index_worm                    !  
   integer     :: i,k,l,m                        ! looping variables

   integer,allocatable   :: j_index(:,:) ! mapping j index to (l,m) index of harmonics 

   ! local variables needed by spherical harmoncis 
   real*8,allocatable    :: p(:)
   integer               :: index
   real*8                :: z 
   real(long)            :: twopi 
   real(long)            :: mphi
    

   if (allocate_q) then 
      ALLOCATE(exp_omega1(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "exp_omega1 allocation error"

      ALLOCATE(exp_omega2(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "exp_omega2 allocation error"

      ALLOCATE(exp_ksq1(0:n(1)/2,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "exp_ksq1 allocation error"

      ALLOCATE(exp_ksq2(0:n(1)/2,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "exp_ksq2 allocation error"

      ALLOCATE(q1(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "q1 allocation error"

      ALLOCATE(q2(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "q2 allocation error"

      ALLOCATE(qr(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "qr allocation error"

      ALLOCATE(qk(0:n(1)/2,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "qk allocation error"
   endif

   ! if wormlike chain exists, allocate Gjj(k) and initialize it. 
   if (allocate_qw) then
      N_sph = (lbar+1)**2 
 
      ! mapping j -> (l,m)
      ALLOCATE(j_index(1:2, 1:N_sph), STAT=error)
      if (error /= 0) stop "j_index allocation error"
      i = 1 
      do l = 0, lbar
         do m = -l, l 
            j_index(:,i) = (/ l, m /)
            i = i + 1
         enddo
      enddo 
 
      ALLOCATE(Yj(0:lebedev_order-1,1:N_sph), STAT=error)
      if (error /= 0) stop "Yj allocation error"
      ! initialization Yj(u,j) in accordance to lebedev_grid 
      ALLOCATE(p(1:(lbar+1)*(lbar+2)/2), STAT=error)
      if (error /= 0) stop "Legendre polynomials allocation error"
      
      twopi = 4.0_long*acos(0.0_long)
      !$OMP PARALLEL 
      do i = 0, lebedev_order-1 
         z       = lebedev_tp_grid(1,i)
         call PlmBar(p, lbar, cos(z)) 
         do k = 1, N_sph
            l = j_index(1,k)
            m = j_index(2,k)
            if (m .GE. 0) then 
               index=PlmIndex(l,  m) 
               mphi = cos(m*lebedev_tp_grid(2,i))/(sqrt(2.0_long*twopi)) 
            else
               index=PlmIndex(l, -m) 
               mphi = sin(-m*lebedev_tp_grid(2,i))/(sqrt(2.0_long*twopi))
            endif 

            Yj(i,k) = p(index)*mphi/(sqrt(2.0_long*twopi))
         enddo
      enddo
      !$OMP END PARALLEL

      ALLOCATE(Gk(0:n(1)-1,0:n(2)-1,0:n(3)-1,1:N_sph,1:N_sph,N_worm_blk_max,N_chain), STAT=error)

      if (error /= 0) stop "Gk allocation error"

      ! initialization Gjj(k) for each wormlike block in all species 
      do i = 1, N_chain
         do k = 1, N_worm_block(i)
            index_worm = Index_worm_block(k,i) 
            if (block_type(index_worm,i)=='Wormlike') then 
               call make_Gk_matrix(k, index_worm, i)   
            endif
         enddo
      enddo

   endif
      contains 
         !----------------------------------------------------------------
         !****p init_step/make_Gk_matrix
         ! SUBROUTINE
         !   make_Gk_matrix(k, index_worm,i)
         ! PURPOSE
         !   evaluate Gk(n1,n2,n3,j,j',index_worm,ichain)
         ! ARGUMENTS
         !   integer  k    - kth wormlike block  
         !   integer  index_worm - index of wormlike block in the chain
         !   integer  i    - ith chain 
         ! SOURCE
         !----------------------------------------------------------------
         subroutine make_Gk_matrix(k, index_worm,i)
            integer,intent(IN) :: index_worm, k, i  
            !*** 

            integer       :: l,m                      ! (l,m) index of spherical harmonics           
            integer       :: j,jp                     ! index of spherical harmonics             
            integer       :: n1,n2,n3                 ! looping variables on k-grid
            integer       :: G(3), Gbz(3)             ! waves in FBZ
            real(long)    :: rjj_element(3)                   ! element of Rjj 
            real(long)    :: deltajj                  ! deltajj=1 when j=j'
            real(long)    :: h                        ! stepsize 
            real(long)    :: blk_length               ! length of kth block in ith chain
            real(long)    :: monomer_size             ! kuhn monomer size
            real(long)    :: twopi                    ! constant 2pi 
            complex(long)       :: Gjjinv(1:N_sph,1:N_sph)  ! local variable to store inverse of Gjj(k)                    

            twopi = 4.0_long * acos(0.0_long)  
            blk_length = block_length(index_worm,i)
            monomer_size =kuhn(block_monomer(index_worm,i)) 
            h = chains(i)%block_ds(index_worm) 

            do n1=0,n(1)-1
               G(1) = n1 
               do n2=0,n(2)-1
                  G(2) = n2 
                  do n3=0,n(3)-1
                     G(3) = n3
                     !$OMP PARALLEL DO collapse(2) private(deltajj,Gbz,l,rjj_element)               
                     do j = 1, N_sph
                        do jp =1, N_sph
                           if (j==jp) then
                              deltajj=1.0_long
                           else
                              deltajj=0.0_long
                           endif
                           Gbz =G_to_bz(G)
                           l = j_index(1,j)
                           rjj_element= Rjj(:,j,jp)

                           Gk(n1,n2,n3,j,jp,k,i) = cmplx(           &
                              11.0_long/6.0_long*deltajj+h*l*(l+1.0_long)*(blk_length/monomer_size)*deltajj ,&  !real part 
                              h*blk_length*twopi*dot_product(rjj_element,Gbz)               )  !imaginary part  

                        enddo 
                     enddo 
                     !$OMP END PARALLEL DO

                     Gjjinv = inverse(Gk(n1,n2,n3,1:N_sph,1:N_sph,k,i),N_sph) 
                     Gk(n1,n2,n3,1:N_sph,1:N_sph,k,i) = Gjjinv  
                 enddo 
              enddo 
            enddo 
            end subroutine make_Gk_matrix 
                 
   end subroutine init_step
   !================================================================



   !----------------------------------------------------------------
   !****p step_mod/make_propg
   ! SUBROUTINE
   !   make_propg(blk_type, ds,b,omega)
   ! PURPOSE
   !   For Gaussian: evaluate exp_ksq, exp_omega's
   !   For Wormlike:  
   ! ARGUMENTS
   !   real(long) ds     - step size
   !   real(long) b      - statistical segment length
   !   real(long) ksq    - k^2 on grid
   !   real(long) omega  - omega field on grid
   ! SOURCE
   !----------------------------------------------------------------
   subroutine make_propg(blk_type,ds,b,omega)
   implicit none

   character(20), intent(IN) :: blk_type 
   real(long)   , intent(IN) :: ds
   real(long)   , intent(IN) :: b
   real(long)   , intent(IN) :: omega(0:,0:,0:)
   !***

   !Auxilliary coef for Gaussian chain 
   real(long)  :: lap_coeff, pot_coeff

   select case(blk_type)
   case ("Gaussian")
      lap_coeff = ds * b**2 /  6.0_long
      pot_coeff = ds / 2.0_long

      exp_ksq1 = exp( - lap_coeff * ksq_grid )
      exp_ksq2 = exp( - lap_coeff / 2.0_long * ksq_grid )

      exp_omega1 = exp( - pot_coeff * omega )
      exp_omega2 = exp( - pot_coeff / 2.0_long * omega )

   case ("Wormlike")

   case default
      write(6,*) 'Error: Invalid blk_type in make_propg'
      call exit(0)
   end select

   end subroutine make_propg
   !================================================================

   !----------------------------------------------------------------
   !****p step_mod/step_wormlike_euler
   ! SUBROUTINE
   !   step_wormlike_euler
   ! PURPOSE
   !   Calculate first few steps in euler method.
   ! ARGUMENTS
   !   real q_in      -  input  q(r,s)
   !   real q_out     -  output q(r,s+-ds)
   !   fft_plan plan  -  see fft_mod for details
   ! SOURCE
   !----------------------------------------------------------------
   subroutine step_wormlike_euler(q_in, q_out, plan) 
   implicit none

   real(long), intent(IN)     :: q_in(0:,0:,0:)
   real(long), intent(OUT)    :: q_out(0:,0:,0:)
   type(fft_plan), intent(IN) :: plan
   !***

   end subroutine step_wormlike_euler


   !----------------------------------------------------------------
   !****p step_mod/step_wormlike_bdf3
   ! SUBROUTINE
   !   step_wormlike_bdf3
   ! PURPOSE
   !   Calculate one step in 
   !   third order backward differentiation formula (BDF3) method for 
   !   integrating the modified diffusion equation of wormlike chain.
   ! ARGUMENTS
   !   real q_in      -  input  q(r,s)
   !   real q_out     -  output q(r,s+-ds)
   !   fft_plan plan  -  see fft_mod for details
   ! SOURCE
   !----------------------------------------------------------------
   subroutine step_wormlike_bdf3(q_in, q_out, plan) 
   implicit none

   real(long), intent(IN)     :: q_in(0:,0:,0:)
   real(long), intent(OUT)    :: q_out(0:,0:,0:)
   type(fft_plan), intent(IN) :: plan
   !***

   end subroutine step_wormlike_bdf3

   !----------------------------------------------------------------
   !****p step_mod/step_gaussian
   ! SUBROUTINE
   !   step_gaussian
   ! PURPOSE
   !   Calculate one step in pseudo-spectral algorithm for 
   !   integrating the modified diffusion equation.
   ! ARGUMENTS
   !   real q_in      -  input  q(r,s)
   !   real q_out     -  output q(r,s+-ds)
   !   fft_plan plan  -  see fft_mod for details
   ! SOURCE
   !----------------------------------------------------------------
   subroutine step_gaussian(q_in, q_out, plan) 
   implicit none

   real(long), intent(IN)     :: q_in(0:,0:,0:)
   real(long), intent(OUT)    :: q_out(0:,0:,0:)
   type(fft_plan), intent(IN) :: plan
   !***

   ! local variables
   real(long)     :: r_npoints

   r_npoints = dble(plan%n(1)*plan%n(2)*plan%n(3))

   qr=exp_omega1*q_in            ! qr(r) = exp{-omega(r)*ds/2)*q_in(r)
   call fft(plan,qr,qk)          ! qk    = F[qr]
   qk=exp_ksq1*qk                ! qk(k) = exp(-ds*(k*b)**2/6)*qk(k)
   call ifft(plan,qk,qr)         ! qr    = F^{-1}[qk]
   q1=exp_omega1*qr/r_npoints    ! q1    = exp^{-omega(r)*ds/2)*qr(r)

   qr=exp_omega2*q_in            ! qr(r) = exp^{-omega(r)*ds/4)*q_in(r)
   call fft(plan,qr,qk)          ! qk    = F[qr]
   qk=exp_ksq2*qk                ! qk(k) = exp(-ds*(k*b)**2/12)*qk(k)
   call ifft(plan,qk,qr)         ! qr    = F^{-1}[qk]
   qr=exp_omega1*qr/r_npoints    ! q2    = exp^{-omega(r)*ds/2)*qr(r)
   call fft(plan,qr,qk)          ! qk    = F[qr]
   qk=exp_ksq2*qk                ! qk(k) = exp(-ds*(k*b)**2/12)*qk(k)
   call ifft(plan,qk,qr)         ! qr    = F^{-1}[qk]
   q2=exp_omega2*qr/r_npoints    ! q2    = exp^{-omega(r)*ds/4)*qr(r)

   q_out = ( 4.0_long * q2 - q1 ) / 3.0_long

   end subroutine step_gaussian
   !================================================================

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Definitions of function set_initial
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !----------------------------------------------------------------
   !****ip step_mod/worm_to_gauss 
   ! SUBROUTINE
   !   worm_to_gauss 
   ! PURPOSE
   !   Set the initial condition of gaussian block
   !   given the previous block is a wormlike block  
   !   
   ! ARGUMENTS
   !   real q_in      -  input  q(r,u,s)
   !   real q_out     -  output q(r,s+ds)
   ! SOURCE
   !----------------------------------------------------------------
   subroutine worm_to_gauss(q_in, q_out) 
   implicit none

   real(long), intent(IN)     :: q_in(0:,0:,0:,0:,0:)
   real(long), intent(INOUT)    :: q_out(0:,0:,0:)
   !***

   end subroutine worm_to_gauss 
   !================================================================


   !----------------------------------------------------------------
   !****ip step_mod/gauss_to_worm 
   ! SUBROUTINE
   !   gauss_to_worm 
   ! PURPOSE
   !   Set the initial condition of wormlike block
   !   given the previous block is a gaussian block  
   !   
   ! ARGUMENTS
   !   real q_in      -  input  q(r,s)
   !   real q_out     -  output q(r,u,s+ds)
   ! SOURCE
   !----------------------------------------------------------------
   subroutine gauss_to_worm(q_in, q_out) 
   implicit none

   real(long), intent(IN)     :: q_in(0:,0:,0:)
   real(long), intent(INOUT)    :: q_out(0:,0:,0:,0:,0:)
   !***

   end subroutine gauss_to_worm 
   !================================================================

end module step_mod
