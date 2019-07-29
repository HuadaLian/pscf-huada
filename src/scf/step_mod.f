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
   use unit_cell_mod, only : G_basis 
   implicit none

   private

   public :: init_step      ! allocate array needed by module
   !***

   public :: make_propg     ! Precalculated arrays for stepping
   public :: step_gaussian
   public :: step_wormlike
   
   !Auxilliary arrays 
   real(long), allocatable    :: exp_omega1(:,:,:)
   real(long), allocatable    :: exp_omega2(:,:,:)
   real(long), allocatable    :: exp_ksq1(:,:,:)
   real(long), allocatable    :: exp_ksq2(:,:,:)

   real(long)   , allocatable    :: q1(:,:,:)
   real(long)   , allocatable    :: q2(:,:,:)
   real(long)   , allocatable    :: qr(:,:,:)
   complex(long), allocatable    :: qk(:,:,:)

   real(long), allocatable    :: expw_omega1(:,:,:,:,:)
   real(long), allocatable    :: expw_omega2(:,:,:,:,:)
   complex(long), allocatable    :: uk(:,:,:,:,:)            ! (kx,ky,kz,theta,phi) 
   complex(long), allocatable    :: cos_isin_uk1(:,:,:,:,:)   ! (kx,ky,kz,theta,phi) 
   complex(long), allocatable    :: cos_isin_uk2(:,:,:,:,:)   ! (kx,ky,kz,theta,phi) 
   complex(long), allocatable    :: cos_isin_ukr1(:,:,:,:,:)   ! (kx,ky,kz,theta,phi) 
   complex(long), allocatable    :: cos_isin_ukr2(:,:,:,:,:)   ! (kx,ky,kz,theta,phi) 

   real(long)   , allocatable    :: exp_cilm(:,:,:) 
   real(long)   , allocatable    :: exp_dif1(:,:,:)             ! (l,m) 
   real(long)   , allocatable    :: exp_dif2(:,:,:)             ! (l,m) 

   real(long)   , allocatable    :: qu1(:,:,:,:,:)
   real(long)   , allocatable    :: qu2(:,:,:,:,:)
   real(long)   , allocatable    :: qru(:,:,:,:,:)
   real(long)   , allocatable    :: qrucilm(:,:,:,:,:,:)
   real(long)   , allocatable    :: qrucilm_inv(:,:,:,:,:,:)
   complex(long), allocatable    :: qku(:,:,:,:,:) 
   !**
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

   integer :: i,j,k,l,m       ! looping variables 
   integer :: error
   real(long) :: theta, phi   ! unit: rad 
   real(long) :: u(3) 
   integer    :: n1, n2, n3
   integer    :: G(3), Gbz(3) 
   real(long) :: G_basis_local(3,3) 
   real(long) :: G_bz_wave(3) 

   ! intensive work begins.
   call omp_set_num_threads(omp_get_num_procs()/2)

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

      ALLOCATE(expw_omega1(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "expw_omega1 allocation error"
      ALLOCATE(expw_omega2(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "expw_omega2 allocation error"

      ALLOCATE(qu1(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "qu1 allocation error"
      ALLOCATE(qu2(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "qu2 allocation error"

      ALLOCATE(qru(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "qru allocation error"
      ALLOCATE(qku(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "qku allocation error"

      ALLOCATE(qrucilm(0:n(1)-1,0:n(2)-1,0:n(3)-1,1:2,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "qrucilm allocation error"
      ALLOCATE(qrucilm_inv(0:n(1)-1,0:n(2)-1,0:n(3)-1,1:2,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "qrucilm2 allocation error"

      ALLOCATE(uk(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "uk allocation error"
      ALLOCATE(cos_isin_uk1(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "cos_isin_uk allocation error"
      ALLOCATE(cos_isin_uk2(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "cos_isin_uk allocation error"
      ALLOCATE(cos_isin_ukr1(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "cos_isin_uk allocation error"
      ALLOCATE(cos_isin_ukr2(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "cos_isin_uk allocation error"

      ALLOCATE(exp_dif1(1:2,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "exp_dif1 allocation error"
      ALLOCATE(exp_dif2(1:2,0:lmax,0:2*lmax), STAT=error)
      if (error /= 0) stop "exp_dif2 allocation error"

      G_basis_local = 0.0_long 
      select case (dim)
      case(1)
         G_basis_local(1,1) = G_basis(1,1)         
      case(2)                                      
         do n1=1,2                                 
         do n2=1,2                                 
            G_basis_local(n1,n2) = G_basis(n1,n2)  
         enddo                                     
         enddo                                     
      case(3)                                      
         do n1=1,3                                 
         do n2=1,3                                 
            G_basis_local(n1,n2) = G_basis(n1,n2)  
         enddo
         enddo                                     
      end select                                   

      ! uk 
      do i= 0, n(1)/2
      G(1) = i 
      do j= 0, n(2)-1
      G(2) = j 
      do k= 0, n(3)-1
      G(3) = k 

      Gbz = G_to_bz(G) 
      G_bz_wave = 1.0_long*(Gbz.dot.G_basis_local)

         do l= 0, lmax 
         do m= 0, 2*lmax 
            theta =   acos(0.0_long)*(zero(lmax-l)+1.0_long) 
            phi   = m*acos(0.0_long)*4.0_long / (2*lmax+1) 
            u     = (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
            uk(i,j,k,l,m) = dot_product( G_bz_wave,u)
         enddo 
         enddo 
      enddo
      enddo
      enddo

   endif

   end subroutine init_step
   !================================================================

   !----------------------------------------------------------------
   !****p step_mod/step_wormlike
   ! SUBROUTINE
   !   step_wormlike
   ! PURPOSE
   !   Integrate the propagator by pseudo-spectral method  
   ! ARGUMENTS
   !   real q_in      -  input  q(r,s)
   !   real q_out     -  output q(r,s+ds)
   !   fft_plan plan  -  see fft_mod for details
   ! SOURCE
   !----------------------------------------------------------------
   subroutine step_wormlike(qwf_in,qwf_out,plan_many,dir) 
   implicit none

   real(long), intent(IN)     :: qwf_in(0:,0:,0:,0:,0:) 
   real(long), intent(OUT)    :: qwf_out(0:,0:,0:,0:,0:) 
   type(fft_plan_many), intent(IN) :: plan_many
   integer,    intent(IN)     :: dir 

   integer        :: i,j,k      
   integer        :: l,p            ! looping variable 

   real(long)     :: r_npoints
   !***

   r_npoints = dble(plan_many%n(1)*plan_many%n(2)*plan_many%n(3))

   ! Full step 
   qru = expw_omega1*qwf_in 
   call fft_many(plan_many, qru, qku) 
   if (dir==-1) then 
      qku = cos_isin_ukr1*qku
   elseif (dir==1) then
   qku = cos_isin_uk1 *qku 
   endif
   call ifft_many(plan_many, qku, qru) 
   qru = qru/r_npoints 
   do i=0,ngrid(1)-1 
   do j=0,ngrid(2)-1 
   do k=0,ngrid(3)-1 
      call SHExpandGLQ(qrucilm(i,j,k,:,:,:), lmax, qru(i,j,k,:,:), weight, plx, csphase=1)  
      qrucilm(i,j,k,:,:,:) = exp_dif1* qrucilm(i,j,k,:,:,:) 
      call MakeGridGLQ(qru(i,j,k,:,:),qrucilm(i,j,k,:,:,:),lmax,plx,csphase=1)
   enddo
   enddo
   enddo
   call fft_many(plan_many, qru, qku) 
   if (dir==-1) then 
      qku = cos_isin_ukr1*qku
   elseif (dir==1) then
   qku = cos_isin_uk1 *qku 
   endif
   call ifft_many(plan_many, qku, qu1) 
   qu1 = qu1 / r_npoints 
   qu1 = expw_omega1*qu1 

   ! Half step h/2 
   qru = expw_omega2*qwf_in  
   call fft_many(plan_many, qru, qku) 
   if (dir==-1) then 
      qku = cos_isin_ukr2*qku
   elseif (dir==1) then
   qku = cos_isin_uk2 *qku 
   endif

   call ifft_many(plan_many, qku, qru) 
   qru = qru / r_npoints 
   do i=0,ngrid(1)-1 
   do j=0,ngrid(2)-1 
   do k=0,ngrid(3)-1 
      call SHExpandGLQ(qrucilm(i,j,k,:,:,:), lmax, qru(i,j,k,:,:), weight, plx, csphase=1)  
      qrucilm(i,j,k,:,:,:) = exp_dif2* qrucilm(i,j,k,:,:,:) 
      call MakeGridGLQ(qru(i,j,k,:,:),qrucilm(i,j,k,:,:,:),lmax,plx,csphase=1)
   enddo
   enddo
   enddo
   call fft_many(plan_many, qru, qku) 
   if (dir==-1) then 
      qku = cos_isin_ukr2*qku
   elseif (dir==1) then
   qku = cos_isin_uk2 *qku 
   endif
   call ifft_many(plan_many, qku, qru) 
   qru = qru / r_npoints 
   qru = expw_omega2*qru  

   qru = expw_omega2*qru  
   call fft_many(plan_many, qru, qku) 
   if (dir==-1) then 
      qku = cos_isin_ukr2*qku
   elseif (dir==1) then
   qku = cos_isin_uk2 *qku 
   endif
   call ifft_many(plan_many, qku, qru) 
   qru = qru / r_npoints 
   do i=0,ngrid(1)-1 
   do j=0,ngrid(2)-1 
   do k=0,ngrid(3)-1 
      call SHExpandGLQ(qrucilm(i,j,k,:,:,:), lmax, qru(i,j,k,:,:), weight, plx, csphase=1)  
      qrucilm(i,j,k,:,:,:) = exp_dif2* qrucilm(i,j,k,:,:,:) 
      call MakeGridGLQ(qru(i,j,k,:,:),qrucilm(i,j,k,:,:,:),lmax,plx,csphase=1)
   enddo
   enddo
   enddo
   call fft_many(plan_many, qru, qku) 
   if (dir==-1) then 
      qku = cos_isin_ukr2*qku
   elseif (dir==1) then
   qku = cos_isin_uk2 *qku 
   endif
   call ifft_many(plan_many, qku, qu2) 
   qu2 = qu2 / r_npoints 
   qu2 = expw_omega2*qu2  

   qwf_out = (4.0_long * qu2 - qu1) / 3.0_long  

   end subroutine step_wormlike

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

   !----------------------------------------------------------------
   !****p step_mod/make_propg
   ! SUBROUTINE
   !   make_propg(ds,b,omega)
   ! PURPOSE
   !   For Gaussian: evaluate exp_ksq, exp_omega's
   ! ARGUMENTS
   !   character(20) blk_type - type of block 
   !   real(long)    ds       - step size
   !   real(long)    b        - statistical segment length
   !   real(long)    ksq      - k^2 on grid
   !   real(long)    omega    - omega field on grid
   ! SOURCE
   !----------------------------------------------------------------
   subroutine make_propg(ds,b,omega)
   implicit none

   real(long)   , intent(IN) :: ds
   real(long)   , intent(IN) :: b
   real(long)   , intent(IN) :: omega(0:,0:,0:)
   !***

   !Auxilliary coef for Gaussian chain 
   real(long)  :: pot_coeff, lap_coeff, adv_coeff, dif_coeff 

   !looping variables 
   integer     :: l,m 
        
   if (allocate_q) then 
      ! Laplace operator 
      lap_coeff = ds * b**2 /  6.0_long
      exp_ksq1 = exp( - lap_coeff * ksq_grid )
      exp_ksq2 = exp( - lap_coeff / 2.0_long * ksq_grid )

      ! Potential term 
      pot_coeff = ds / 2.0_long
      exp_omega1 = exp( - pot_coeff * omega )
      exp_omega2 = exp( - pot_coeff / 2.0_long * omega )
   endif 

   if (allocate_qw) then 
      ! Potential term 
      pot_coeff = ds / 2.0_long

      !$OMP PARALLEL DO COLLAPSE(2) 
      do l=0,lmax 
      do m=0,2*lmax 
         expw_omega1(:,:,:,l,m) = exp( - pot_coeff * omega )
         expw_omega2(:,:,:,l,m) = exp( - pot_coeff / 2.0_long * omega )
      enddo
      enddo
      !$OMP END PARALLEL DO 

      ! Advection term 
      adv_coeff   = ds / 2.0_long 
      cos_isin_uk1 = cos( adv_coeff * uk)+ cmplx(0.0_long,1.0_long) * sin( adv_coeff * uk)
      cos_isin_ukr1 = cos(-adv_coeff * uk)+ cmplx(0.0_long,1.0_long) * sin(-adv_coeff * uk)

      adv_coeff   = (ds / 2.0_long) / 2.0_long  
      cos_isin_uk2 = cos( adv_coeff * uk)+ cmplx(0.0_long,1.0_long) * sin( adv_coeff * uk)
      cos_isin_ukr2 = cos(-adv_coeff * uk)+ cmplx(0.0_long,1.0_long) * sin(-adv_coeff * uk)

      ! Diffusion term 
      dif_coeff   = ds/b  
      do l=0,lmax
         exp_dif1(:,l,:)  = exp( -dif_coeff * l * (l+1) ) 
         exp_dif2(:,l,:)  = exp( -dif_coeff /2.0_long * l * (l+1) ) 
      enddo 
   endif 

   end subroutine make_propg
   !================================================================

end module step_mod
