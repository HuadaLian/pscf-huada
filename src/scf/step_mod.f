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
   !***

   ! Generic interfaces 
   public :: make_propg     ! Precalculated arrays for stepping

   public :: step_gaussian
   public :: step_wormlike_euler
   public :: step_wormlike_bdf3
   public :: qw_decompose
   public :: qw_sum 
   
   !Auxilliary array for Gaussian chain
   real(long), allocatable    :: omega_local(:,:,:)
   real(long), allocatable    :: exp_omega1(:,:,:)
   real(long), allocatable    :: exp_omega2(:,:,:)
   real(long), allocatable    :: exp_ksq1(:,:,:)
   real(long), allocatable    :: exp_ksq2(:,:,:)

   real(long)   , allocatable    :: q1(:,:,:)
   real(long)   , allocatable    :: q2(:,:,:)
   real(long)   , allocatable    :: qr(:,:,:)
   complex(long), allocatable :: qk(:,:,:)

   real(long)   , allocatable    :: qwj1(:,:,:,:) 
   real(long)   , allocatable    :: qwj2(:,:,:,:) 
   real(long)   , allocatable    :: qwjr(:,:,:,:) 

   complex(long), allocatable :: qwjk(:,:,:,:)
   real(long)   , allocatable    :: omegaqwr(:,:,:,:)
   complex(long), allocatable :: omegaqwk(:,:,:,:)
   real(long)                 :: ds_blk_local 
   real(long)   , allocatable :: Pjr(:,:,:,:)
   complex(long), allocatable :: Pjk(:,:,:,:) 
   complex(long), allocatable :: Gjjqw(:,:,:) 
   complex(long), allocatable :: Gklocal(:,:,:,:,:) 
   complex(long), allocatable :: Gkhalfslocal(:,:,:,:,:) 
   complex(long), allocatable :: Gkinvlocal(:,:,:,:,:) 


   !Auxilliary array for Wormlike chain
   !real(long), allocatable :: Yjf(:,:)              ! 2D - real spherical harmonics on orientation grid u 
   !real(long), allocatable :: Yjr(:,:)              ! 2D - real spherical harmonics on orientation grid u 

   !real(long), allocatable :: Yjfinv(:,:)           ! inverse of Yjf 
   !real(long), allocatable :: Yjrinv(:,:)           ! inverse of Yjr 

   complex, allocatable :: Gkf(:,:,:,:,:,:)       ! For euler method 
   complex, allocatable :: Gkr(:,:,:,:,:,:)       ! For euler method 

   complex, allocatable :: Gkfhalfs(:,:,:,:,:,:)  ! For richardson extrapolation in euler method
   complex, allocatable :: Gkrhalfs(:,:,:,:,:,:)  ! For richardson extrapolation in euler method

   complex, allocatable :: Gkfinv(:,:,:,:,:,:) ! 2D - vector function Gk(n1,n2,n3,j,j',i_blk,i_chain)    
   complex, allocatable :: Gkrinv(:,:,:,:,:,:) ! 2D - vector function Gk(n1,n2,n3,j,j',i_blk,i_chain)    
   !**

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Generic Interfaces
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   !------------------------------------------------------------------
   !****p step_mod/make_propg
   ! FUNCTION
   !    make_propg
   !
   ! COMMENT
   !         
   !    make_propg_gaussian : 
   !    make_propg_wormlike : 
   ! SOURCE
   !------------------------------------------------------------------
   interface make_propg 
      module procedure make_propg_gaussian 
      module procedure make_propg_wormlike 
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

   integer     :: index_worm                    !  
   integer     :: i,k,l,m                        ! looping variables

   integer,allocatable   :: j_index(:,:) ! mapping j index to (l,m) index of harmonics 

   ! local variables needed by spherical harmoncis 
   real*8,allocatable    :: p(:)
   integer               :: index
   real*8                :: z 
   real(long)            :: twopi 
   real(long)            :: mphi


   ! intensive work begins.
   call omp_set_num_threads(omp_get_num_procs()/2)

   ALLOCATE(omega_local(0:n(1)-1,0:n(2)-1,0:n(3)-1), STAT=error)
   if (error /= 0) stop "omega_local allocation error"

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

      ALLOCATE(qwj1(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "qwj1 allocation error"

      ALLOCATE(qwj2(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "qwj2 allocation error"

      ALLOCATE(qwjr(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "qwjr allocation error"

      ALLOCATE(qwjk(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "qwjk allocation error"

      ALLOCATE(omegaqwr(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "omegaqwr allocation error"

      ALLOCATE(omegaqwk(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "omegaqwk allocation error"

      ALLOCATE(Pjr(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "Pjr allocation error"

      ALLOCATE(Pjk(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "Pjk allocation error"

      ALLOCATE(Gklocal(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "omega_local allocation error"

      ALLOCATE(Gkinvlocal(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "omega_local allocation error"

      ALLOCATE(Gkhalfslocal(0:n(1)/2,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1), STAT=error)
      if (error /= 0) stop "omega_local allocation error"                

      ALLOCATE(Gjjqw(0:n(1)/2,0:n(2)-1,0:n(3)-1), STAT=error)
      if (error /= 0) stop "Gjjqw allocation error"

      ! mapping j -> (l,m)
      ALLOCATE(j_index(1:2, 0:N_sph-1), STAT=error)
      if (error /= 0) stop "j_index allocation error"
      i = 0 
      do l = 0, lbar
         do m = -l, l 
            j_index(:,i) = (/ l, m /)
            i = i + 1
         enddo
      enddo 
 
      ! initialization Yj(u,j) in accordance to gauss-kronrod grid  

      !ALLOCATE(p((lbar+1)*(lbar+2)/2), STAT=error)
      !if (error /= 0) stop "Legendre polynomials allocation error"
      !
      !twopi = 4.0_long*acos(0.0_long)
      !do i = 0, N_sph-1 
      !   z       = angularf_grid(1,i)
      !   call PlmBar(p, lbar, cos(z)) 
      !   do k = 0, N_sph-1
      !      l = j_index(1,k)
      !      m = j_index(2,k)
      !      if (m .GE. 0) then 
      !         index=PlmIndex(l,  m) 
      !         mphi = cos(m*angularf_grid(2,i))/(sqrt(2.0_long*twopi)) 
      !      else
      !         index=PlmIndex(l, -m) 
      !         mphi = sin(-m*angularf_grid(2,i))/(sqrt(2.0_long*twopi))
      !      endif 
      !      Yjf(k,i) = p(index)*mphi
      !      Yjr(k,i) =-p(index)*mphi 
      !   enddo

      !enddo

      !Yjfinv = inverse(Yjf,N_sph) 
      !Yjrinv = -1.0_long*Yjfinv  

      ALLOCATE(Gkf(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1,1:N_worm_blk_max), STAT=error)
      if (error /= 0) stop "Gk allocation error"
      ALLOCATE(Gkr(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1,1:N_worm_blk_max), STAT=error)
      if (error /= 0) stop "Gk allocation error"

      ALLOCATE(Gkfinv(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1,1:N_worm_blk_max), STAT=error)
      if (error /= 0) stop "Gkinv allocation error"
      ALLOCATE(Gkrinv(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1,1:N_worm_blk_max), STAT=error)
      if (error /= 0) stop "Gkinv allocation error"

      ALLOCATE(Gkfhalfs(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1,1:N_worm_blk_max), STAT=error)
      if (error /= 0) stop "Gkhalfs allocation error"
      ALLOCATE(Gkrhalfs(0:n(1)-1,0:n(2)-1,0:n(3)-1,0:N_sph-1,0:N_sph-1,1:N_worm_blk_max), STAT=error)
      if (error /= 0) stop "Gkhalfs allocation error"

      ! initialization Gjj(k) for each wormlike block in all species 
      call make_Gk_matrix()   

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
         subroutine make_Gk_matrix()
            use unit_cell_mod, only : G_basis
            integer   :: index_worm, k, i  
            !*** 

            integer       :: l,m                      ! (l,m) index of spherical harmonics           
            integer       :: j,jp                     ! index of spherical harmonics             
            integer       :: n1,n2,n3                 ! looping variables on k-grid
            integer       :: G(3), Gbz(3)             ! waves in FBZ
            real(long)    :: G_basis_local(3,3)       ! G_basis_copy
            real(long)    :: G_bz_wave(3)
            real(long)    :: rjj_element(3)                   ! element of Rjj 
            real(long)    :: deltajj                  ! deltajj=1 when j=j'
            real(long)    :: h                        ! stepsize 
            real(long)    :: monomer_size             ! kuhn monomer size
            complex(long) :: Gjjinv(0:N_sph-1,0:N_sph-1)  ! local variable to store inverse of Gjj(k)                    


            G_basis_local = 0.0 
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

            do i = 1,N_chain
            do k = 1,N_block(i) 
               if (block_type(k,i) =='Wormlike') then 
                  monomer_size =  kuhn(block_monomer(k,i)) 
                  h            =  chains(i)%block_ds(k) 
                  index_worm   =  Index_worm_block(k,i) 

                  do n1=0,n(1)-1
                     G(1) = n1 
                     do n2=0,n(2)-1
                        G(2) = n2 
                        do n3=0,n(3)-1
                           G(3) = n3
                           Gbz =G_to_bz(G)
                           G_bz_wave = Gbz.dot.G_basis_local

                           do j = 0, N_sph-1
                              do jp =0, N_sph-1
                                 if (j==jp) then
                                    deltajj=1.0_long
                                 else
                                    deltajj=0.0_long
                                 endif
                                 l = j_index(1,j)
                                 rjj_element= Rjj(:,j,jp)

                                 Gkf(n1,n2,n3,j,jp,index_worm)      = cmplx(           &
                                    11.0_long/6.0_long*deltajj+h*l*(l+1.0_long)*(chain_length(i)/monomer_size)*deltajj ,&  !real part 
                                    -h*chain_length(i)*dot_product(rjj_element,G_bz_wave)    )  !imaginary part  
                                 Gkr(n1,n2,n3,j,jp,index_worm)      = cmplx(           &
                                    11.0_long/6.0_long*deltajj+h*l*(l+1.0_long)*(chain_length(i)/monomer_size)*deltajj ,&  !real part 
                                    +h*chain_length(i)*dot_product(rjj_element,G_bz_wave)   )  !imaginary part  
                              enddo 
                           enddo 

                           ! Gkinv for bdf3 
                           Gkfinv(n1,n2,n3,0:N_sph-1,0:N_sph-1,index_worm) = &
                                      inverse(Gkf(n1,n2,n3,0:N_sph-1,0:N_sph-1,index_worm),N_sph) 
                           Gkrinv(n1,n2,n3,0:N_sph-1,0:N_sph-1,index_worm) = &
                                      inverse(Gkr(n1,n2,n3,0:N_sph-1,0:N_sph-1,index_worm),N_sph) 

                           ! Gk and Gkhalfs for euler method 
                           do j = 0, N_sph-1
                              do jp =0, N_sph-1
                                 if (j==jp) then
                                    deltajj=1.0_long
                                 else
                                    deltajj=0.0_long
                                 endif
                                 l = j_index(1,j)
                                 rjj_element= Rjj(:,j,jp)

                                 Gkf(n1,n2,n3,j,jp,index_worm)      = cmplx(                              &
                                    deltajj-h*l*(l+1.0_long)*chain_length(i)/monomer_size*deltajj,     &           ! real part 
                                    +h*chain_length(i)*dot_product(rjj_element,G_bz_wave) )  !imaginary part  
                                 Gkr(n1,n2,n3,j,jp,index_worm)      = cmplx(                              &
                                    deltajj-h*l*(l+1.0_long)*chain_length(i)/monomer_size*deltajj,     &           ! real part 
                                    -h*chain_length(i)*dot_product(rjj_element,G_bz_wave)   )  !imaginary part  

                                 Gkfhalfs(n1,n2,n3,j,jp,index_worm)      = cmplx(                             &
                                    deltajj-h/2.0_long*l*(l+1.0_long)*chain_length(i)/monomer_size*deltajj,&                  !real part 
                                    +h/2.0_long*chain_length(i)*dot_product(rjj_element,G_bz_wave)        )  !imaginary part

                                 Gkrhalfs(n1,n2,n3,j,jp,index_worm)      = cmplx(                             &
                                    deltajj-h/2.0_long*l*(l+1.0_long)*chain_length(i)/monomer_size*deltajj,&                  !real part 
                                    -h/2.0_long*chain_length(i)*dot_product(rjj_element,G_bz_wave)        )  !imaginary part
                              enddo 
                           enddo 

                       enddo 
                    enddo 
                  enddo 
               endif
               enddo
               enddo
            end subroutine make_Gk_matrix 
                 
   end subroutine init_step
   !================================================================



   !----------------------------------------------------------------
   !****p step_mod/qw_decompose
   ! SUBROUTINE
   !   qw_decompose
   ! PURPOSE
   !   expand the propagator of wormlike chain in spherical harmonics 
   !   q(r,u,s) -> sum ( q_j(r,s) Ylm(u) )
   ! ARGUMENTS
   !   real q_in      -  input  q(r,u,s)
   !   real q_out     -  output q(r,j,s)
   ! SOURCE
   !----------------------------------------------------------------
   subroutine qw_decompose(qw, qwj,dir) 
   implicit none

   real(long), intent(IN)     :: qw(0:,0:,0:,0:)
   real(long), intent(OUT)    :: qwj(0:,0:,0:,0:)
   integer,    intent(IN)     :: dir 

   integer                    :: i,j,k,l   !looping variables 
   !***

   if (dir == 1) then 
      !$OMP PARALLEL DO COLLAPSE(4)
      do i=0,ngrid(1)-1
      do j=0,ngrid(2)-1
      do k=0,ngrid(3)-1
      do l=0,N_sph-1
         qwj(i,j,k,l) = Yjfinv(l,:).dot.qw(i,j,k,:)
      enddo
      enddo
      enddo
      enddo 
      !$OMP END PARALLEL DO 
   elseif (dir == -1) then
      !$OMP PARALLEL DO COLLAPSE(4)
      do i=0,ngrid(1)-1
      do j=0,ngrid(2)-1
      do k=0,ngrid(3)-1
      do l=0,N_sph-1
         qwj(i,j,k,l) = Yjrinv(l,:).dot.qw(i,j,k,:)
      enddo
      enddo
      enddo
      enddo 
      !$OMP END PARALLEL DO 
   else
      stop 'wrong dir in qw_decompose'
   endif

   end subroutine qw_decompose


   !----------------------------------------------------------------
   !****p step_mod/qw_sum
   ! SUBROUTINE
   !   qw_sum
   ! PURPOSE
   !   expand the propagator of wormlike chain in spherical harmonics 
   !   q(r,u,s) -> sum ( q_j(r,s) Ylm(u) )
   ! ARGUMENTS
   !   real q_in      -  input  q(r,u,s)
   !   real q_out     -  output q(r,j,s)
   ! SOURCE
   !----------------------------------------------------------------
   subroutine qw_sum(qwj, qw,dir) 
   implicit none

   real(long), intent(IN)     :: qwj(0:,0:,0:,0:)
   real(long), intent(OUT)    :: qw(0:,0:,0:,0:)
   integer,    intent(IN)     :: dir 

   integer                    :: i,j,k,l   !looping variables 
   !***
   if (dir == 1) then 
      !$OMP PARALLEL DO COLLAPSE(4)
      do i=0,ngrid(1)-1
      do j=0,ngrid(2)-1
      do k=0,ngrid(3)-1
      do l=0,N_sph-1
         qw(i,j,k,l) = dot_product(Yjf(l,:),qwj(i,j,k,:))
      enddo
      enddo
      enddo
      enddo 
      !$OMP END PARALLEL DO 
   elseif (dir == -1) then 
      !$OMP PARALLEL DO COLLAPSE(4)
      do i=0,ngrid(1)-1
      do j=0,ngrid(2)-1
      do k=0,ngrid(3)-1
      do l=0,N_sph-1
         qw(i,j,k,l) = dot_product(Yjr(l,:),qwj(i,j,k,:))
      enddo
      enddo
      enddo
      enddo 
      !$OMP END PARALLEL DO 
   else
      stop 'wrong dir in qw_sum '
   endif

   end subroutine qw_sum

   !----------------------------------------------------------------
   !****p step_mod/step_wormlike_euler
   ! SUBROUTINE
   !   step_wormlike_euler
   ! PURPOSE
   !   Calculate first few steps in modified euler method.
   ! ARGUMENTS
   !   real q_in      -  input  q(r,s)
   !   real q_out     -  output q(r,s+-ds)
   !   fft_plan plan  -  see fft_mod for details
   ! SOURCE
   !----------------------------------------------------------------
   subroutine step_wormlike_euler(qwj_in, qwj_out, qwf_out,plan_many,dir) 
   implicit none

   real(long), intent(IN)     :: qwj_in(0:,0:,0:,0:)  
   real(long), intent(OUT)    :: qwj_out(0:,0:,0:,0:) 

   real(long), intent(OUT)    :: qwf_out(0:,0:,0:,0:) 

   type(fft_plan_many), intent(IN) :: plan_many
   integer,    intent(IN)     :: dir 

   integer        :: i,j,k      
   integer        :: l,p            ! looping variable 

   real(long)     :: r_npoints
   !***

   r_npoints = dble(plan_many%n(1)*plan_many%n(2)*plan_many%n(3))

   ! half step size 
   call fft(plan_many, qwj_in, qwjk) 
   !$OMP PARALLEL
   do l=0,N_sph-1
      omegaqwr(:,:,:,l) = -1.0_long*ds_blk_local/2.0_long * omega_local * qwj_in(:,:,:,l)        ! (- or + )ds/2*omega(r)*qwj(r,s) 
   enddo
   !$OMP END PARALLEL 
   call fft(plan_many, omegaqwr, omegaqwk) 
   !$OMP PARALLEL DO private(Gjjqw)
   do l = 0, N_sph -1 
      Gjjqw = cmplx(0.0_long,0.0_long)                         
      do p=0, N_sph-1
         Gjjqw = Gjjqw + Gkhalfslocal(:,:,:,l,p)*qwjk(:,:,:,p)      ! sum_{j,jp} G_{j,jp}(k)*qw_jp(k)
      enddo 
      omegaqwk(:,:,:,l) = omegaqwk(:,:,:,l) + Gjjqw  
   enddo 
   !$OMP END PARALLEL DO
   call ifft(plan_many, omegaqwk, qwj2) 
   qwj2 = qwj2/r_npoints             

   call fft(plan_many, qwj2, qwjk) 
   ! half step size 
   !$OMP PARALLEL
   do l=0,N_sph-1
      omegaqwr(:,:,:,l) = -1.0_long*ds_blk_local/2.0_long * omega_local * qwj2(:,:,:,l)        ! (- or + )ds/2*omega(r)*qwj(r,s) 
   enddo
   !$OMP END PARALLEL 
   call fft(plan_many, omegaqwr, omegaqwk) 
   !$OMP PARALLEL DO private(Gjjqw)
   do l = 0, N_sph -1 
      Gjjqw = cmplx(0.0_long,0.0_long)                         
      do p=0, N_sph-1
         Gjjqw = Gjjqw + Gkhalfslocal(:,:,:,l,p)*qwjk(:,:,:,p)      ! sum_{j,jp} G_{j,jp}(k)*qw_jp(k)
      enddo 
      omegaqwk(:,:,:,l) = omegaqwk(:,:,:,l) + Gjjqw  
   enddo 
   !$OMP END PARALLEL DO
   call ifft(plan_many, omegaqwk, qwj2) 
   qwj2 = qwj2/r_npoints             


   ! full step size 
   call fft(plan_many, qwj_in, qwjk) 
   !$OMP PARALLEL
   do l=0,N_sph-1
      omegaqwr(:,:,:,l) = -1.0_long*ds_blk_local * omega_local * qwj_in(:,:,:,l)        ! ds/2*omega(r)*qwj(r,s) 
   enddo
   !$OMP END PARALLEL 

   call fft(plan_many, omegaqwr, omegaqwk) 

   !$OMP PARALLEL DO private(Gjjqw)
   do l = 0, N_sph-1 
      Gjjqw = cmplx(0.0_long,0.0_long)                         
      do p=0, N_sph-1
         Gjjqw = Gjjqw + Gklocal(:,:,:,l,p)*qwjk(:,:,:,p)      ! sum_{j,jp} G_{j,jp}(k)*qw_jp(k)
      enddo 
      omegaqwk(:,:,:,l) = omegaqwk(:,:,:,l) + Gjjqw  
   enddo 
   !$OMP END PARALLEL DO

   call ifft(plan_many, omegaqwk, qwj1) 
   qwj1 = qwj1/r_npoints             

   !$OMP PARALLEL
   do l=0,N_sph-1 
      qwj_out(:,:,:,l) = (2.0_long*qwj2(:,:,:,l)-qwj1(:,:,:,l))/1.0_long
   enddo 
   !$OMP END PARALLEL

   call qw_sum(qwj_out,qwf_out,dir)

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
   subroutine step_wormlike_bdf3(qwj3h_in,qwj2h_in,qwj1h_in,qwj_out,qwf_out,plan_many,dir) 
   implicit none

   real(long), intent(IN)     :: qwj3h_in(0:,0:,0:,0:)  
   real(long), intent(IN)     :: qwj2h_in(0:,0:,0:,0:)  
   real(long), intent(IN)     :: qwj1h_in(0:,0:,0:,0:)  

   real(long), intent(OUT)    :: qwj_out(0:,0:,0:,0:) 
   real(long), intent(OUT)    :: qwf_out(0:,0:,0:,0:) 

   type(fft_plan_many), intent(IN) :: plan_many
   integer,    intent(IN)     :: dir 

   integer        :: i,j,k      
   integer        :: l,p            ! looping variable 

   real(long)     :: r_npoints
   !***

   r_npoints = dble(plan_many%n(1)*plan_many%n(2)*plan_many%n(3))

   !$OMP PARALLEL DO 
   do l=0,N_sph-1
      Pjr(:,:,:,l) = 3.0_long*qwj1h_in(:,:,:,l)&
                    -3.0_long/2.0_long*qwj2h_in(:,:,:,l) &
                    +1.0_long/3.0_long*qwj3h_in(:,:,:,l) & 
                    -ds_blk_local*omega_local*(3.0_long*qwj1h_in(:,:,:,l)&
                    -3.0_long*qwj2h_in(:,:,:,l) + qwj3h_in(:,:,:,l))
   enddo 
   !$OMP END PARALLEL DO 

   call fft(plan_many, Pjr, Pjk) 

   qwjk = cmplx(0.0_long, 0.0_long)
   !$OMP PARALLEL DO COLLAPSE(4)
   do i=0,ngrid(1)/2
   do j=0,ngrid(2)-1
   do k=0,ngrid(3)-1
   do l=0,N_sph-1
      qwjk(i,j,k,l) = Gkinvlocal(i,j,k,l,:).dot.Pjk(i,j,k,:)
   enddo
   enddo
   enddo
   enddo
   !$OMP END PARALLEL DO 

   call ifft(plan_many, qwjk, qwj_out) 
   qwj_out = qwj_out /r_npoints  

   call qw_sum(qwj_out,qwf_out,dir)

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
   ! Definitions of function make_propg 
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !----------------------------------------------------------------
   !****p step_mod/make_propg_gaussian
   ! SUBROUTINE
   !   make_propg(blk_type, ds,b,omega)
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
   subroutine make_propg_gaussian(ds,b,omega)
   implicit none

   real(long)   , intent(IN) :: ds
   real(long)   , intent(IN) :: b
   real(long)   , intent(IN) :: omega(0:,0:,0:)
   !***

   !Auxilliary coef for Gaussian chain 
   real(long)  :: lap_coeff, pot_coeff

   lap_coeff = ds * b**2 /  6.0_long
   pot_coeff = ds / 2.0_long

   exp_ksq1 = exp( - lap_coeff * ksq_grid )
   exp_ksq2 = exp( - lap_coeff / 2.0_long * ksq_grid )

   exp_omega1 = exp( - pot_coeff * omega )
   exp_omega2 = exp( - pot_coeff / 2.0_long * omega )
   end subroutine make_propg_gaussian
   !================================================================



   !----------------------------------------------------------------
   !****p step_mod/make_propg_wormlike
   ! SUBROUTINE
   !   make_propg(ds,b,omega, gk, gkhalfs)
   ! PURPOSE
   !   For Gaussian: evaluate exp_ksq, exp_omega's
   ! ARGUMENTS
   !   real(long)    ds       - step size
   !   real(long)    b        - statistical segment length
   !   real(long)    ksq      - k^2 on grid
   !   real(long)    omega    - omega field on grid
   !   complex(long) Gk       - Gjjk array                         (optional)
   !   complex(long) Gkhalfs  - Gjjk array with half step size     (optional)            
   ! SOURCE
   !----------------------------------------------------------------
   subroutine make_propg_wormlike(ds,b,omega,index_worm,dir)
   implicit none

   real(long)   , intent(IN) :: ds
   real(long)   , intent(IN) :: b
   real(long)   , intent(IN) :: omega(0:,0:,0:)
   integer      , intent(IN) :: index_worm
   integer      , intent(IN) :: dir ! direction 
   !***

   ds_blk_local = ds  
   omega_local  = omega
   if (dir==1) then 
      Gklocal     = Gkf(0:ngrid(1)/2,:,:,:,:,index_worm) 
      Gkhalfslocal= Gkfhalfs(0:ngrid(1)/2,:,:,:,:,index_worm) 
      Gkinvlocal  = Gkfinv(0:ngrid(1)/2,:,:,:,:,index_worm) 
   elseif (dir== -1) then 
      Gklocal     = Gkr(0:ngrid(1)/2,:,:,:,:,index_worm) 
      Gkhalfslocal= Gkrhalfs(0:ngrid(1)/2,:,:,:,:,index_worm) 
      Gkinvlocal  = Gkrinv(0:ngrid(1)/2,:,:,:,:,index_worm) 
   else
      stop 'Could not make propgator.'
   endif



   end subroutine make_propg_wormlike
   !================================================================

end module step_mod
