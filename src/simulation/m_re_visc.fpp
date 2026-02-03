!>
!! @file m_re_visc.f90
!! @brief Contains module m_re_visc for computing Re_visc (Reynolds number for viscosity)
!!        Supports both Newtonian and non-Newtonian fluids for 1D, 2D, and 3D cases
#:include 'macros.fpp'

!> @brief The module contains two routines that compute viscosity-related quantities:
!!
!!  1. s_compute_re_visc: returns Re_visc_per_phase(fluid, dir) = 1/mu (inverse viscosity).
!!     USE FOR: mixture Re (CFL / viscous dt in m_variables_conversion), stress terms
!!     (m_viscous), and Riemann solvers (m_riemann_solvers). Callers form mixture Re
!!     as 1/sum(alpha/Re_visc_per_phase) = 1/sum(alpha*mu).
!!
!!  2. s_compute_phase_viscosity_at_cell: returns mu_shear_phases, mu_bulk_phases = mu (viscosity).
!!     USE FOR: gradients of viscosity (dmu/dx, dmu/dy, dmu/dz) and face reconstruction
!!     in m_viscous when actual viscosity values are needed.
module m_re_visc

    use m_derived_types        !< Definitions of the derived types
    use m_global_parameters    !< Definitions of the global parameters
    use m_hb_function          !< Herschel-Bulkley non-Newtonian viscosity

    private; public :: s_compute_re_visc, &
 s_compute_phase_viscosity_at_cell

contains

    !> Computes velocity gradients at a single cell (j,k,l) using finite differences
    !! Works for 1D, 2D, and 3D cases
    !! For non-Newtonian fluids with sufficient buffer, uses 4th order central differences
    !! For Newtonian fluids or near boundaries, uses 2nd order differences
    !! Bounds checking uses idwbuff to ensure safe array access in MPI decomposed domains
    !! @param q_prim_vf Primitive variables
    !! @param j x index
    !! @param k y index
    !! @param l z index
    !! @param D_xx Output: du/dx
    !! @param D_yy Output: dv/dy
    !! @param D_zz Output: dw/dz
    !! @param D_xy Output: 0.5*(du/dy + dv/dx)
    !! @param D_xz Output: 0.5*(du/dz + dw/dx)
    !! @param D_yz Output: 0.5*(dv/dz + dw/dy)
    pure subroutine s_compute_velocity_gradients_at_cell(q_prim_vf, j, k, l, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: j, k, l
        real(wp), intent(out) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz

        logical :: use_4th_order_x, use_4th_order_y, use_4th_order_z
        logical :: j_in_bounds, k_in_bounds, l_in_bounds
        real(wp) :: coeff_4th
        integer :: j_lo, j_hi, k_lo, k_hi, l_lo, l_hi

        ! Array bounds (including ghost cells)
        j_lo = idwbuff(1)%beg
        j_hi = idwbuff(1)%end
        k_lo = idwbuff(2)%beg
        k_hi = idwbuff(2)%end
        l_lo = idwbuff(3)%beg
        l_hi = idwbuff(3)%end

        ! First check if the cell itself is within bounds
        j_in_bounds = (j >= j_lo) .and. (j <= j_hi)
        k_in_bounds = (k >= k_lo) .and. (k <= k_hi)
        l_in_bounds = (l >= l_lo) .and. (l <= l_hi)

        ! If the cell is out of bounds, return zero gradients
        if (.not. (j_in_bounds .and. k_in_bounds .and. l_in_bounds)) then
            D_xx = 0._wp
            D_yy = 0._wp
            D_zz = 0._wp
            D_xy = 0._wp
            D_xz = 0._wp
            D_yz = 0._wp
            return
        end if

        ! Check if we can use 4th order (non-Newtonian with sufficient buffer AND safe array access)
        use_4th_order_x = any_non_newtonian .and. (buff_size >= 4) .and. &
                          (j - 2 >= j_lo) .and. (j + 2 <= j_hi)
        use_4th_order_y = any_non_newtonian .and. (buff_size >= 4) .and. &
                          (k - 2 >= k_lo) .and. (k + 2 <= k_hi) .and. (n > 0)
        use_4th_order_z = any_non_newtonian .and. (buff_size >= 4) .and. &
                          (l - 2 >= l_lo) .and. (l + 2 <= l_hi) .and. (p > 0)

        ! Compute D_xx = du/dx
        if (use_4th_order_x) then
            coeff_4th = 1._wp/(x_cc(j - 2) - 8._wp*x_cc(j - 1) - x_cc(j + 2) + 8._wp*x_cc(j + 1))
            D_xx = coeff_4th*(q_prim_vf(momxb)%sf(j - 2, k, l) - 8._wp*q_prim_vf(momxb)%sf(j - 1, k, l) + &
                              8._wp*q_prim_vf(momxb)%sf(j + 1, k, l) - q_prim_vf(momxb)%sf(j + 2, k, l))
        else if (j - 1 >= j_lo .and. j + 1 <= j_hi) then
            ! 2nd order central difference
            D_xx = (q_prim_vf(momxb)%sf(j + 1, k, l) - q_prim_vf(momxb)%sf(j - 1, k, l)) / &
                   (x_cc(j + 1) - x_cc(j - 1))
        else if (j + 1 <= j_hi) then
            ! Forward difference
            D_xx = (q_prim_vf(momxb)%sf(j + 1, k, l) - q_prim_vf(momxb)%sf(j, k, l)) / &
                   (x_cc(j + 1) - x_cc(j))
        else if (j - 1 >= j_lo) then
            ! Backward difference
            D_xx = (q_prim_vf(momxb)%sf(j, k, l) - q_prim_vf(momxb)%sf(j - 1, k, l)) / &
                   (x_cc(j) - x_cc(j - 1))
        else
            D_xx = 0._wp
        end if

        ! Compute D_yy = dv/dy (2D and 3D only)
        if (n > 0) then
            if (use_4th_order_y) then
                coeff_4th = 1._wp/(y_cc(k - 2) - 8._wp*y_cc(k - 1) - y_cc(k + 2) + 8._wp*y_cc(k + 1))
                D_yy = coeff_4th*(q_prim_vf(momxb + 1)%sf(j, k - 2, l) - 8._wp*q_prim_vf(momxb + 1)%sf(j, k - 1, l) + &
                                  8._wp*q_prim_vf(momxb + 1)%sf(j, k + 1, l) - q_prim_vf(momxb + 1)%sf(j, k + 2, l))
            else if (k - 1 >= k_lo .and. k + 1 <= k_hi) then
                ! 2nd order central difference
                D_yy = (q_prim_vf(momxb + 1)%sf(j, k + 1, l) - q_prim_vf(momxb + 1)%sf(j, k - 1, l)) / &
                       (y_cc(k + 1) - y_cc(k - 1))
            else if (k + 1 <= k_hi) then
                ! Forward difference
                D_yy = (q_prim_vf(momxb + 1)%sf(j, k + 1, l) - q_prim_vf(momxb + 1)%sf(j, k, l)) / &
                       (y_cc(k + 1) - y_cc(k))
            else if (k - 1 >= k_lo) then
                ! Backward difference
                D_yy = (q_prim_vf(momxb + 1)%sf(j, k, l) - q_prim_vf(momxb + 1)%sf(j, k - 1, l)) / &
                       (y_cc(k) - y_cc(k - 1))
            else
                D_yy = 0._wp
            end if
        else
            D_yy = 0._wp
        end if

        ! Compute D_zz = dw/dz (3D only)
        if (p > 0) then
            if (use_4th_order_z) then
                coeff_4th = 1._wp/(z_cc(l - 2) - 8._wp*z_cc(l - 1) - z_cc(l + 2) + 8._wp*z_cc(l + 1))
                D_zz = coeff_4th*(q_prim_vf(momxb + 2)%sf(j, k, l - 2) - 8._wp*q_prim_vf(momxb + 2)%sf(j, k, l - 1) + &
                                  8._wp*q_prim_vf(momxb + 2)%sf(j, k, l + 1) - q_prim_vf(momxb + 2)%sf(j, k, l + 2))
            else if (l - 1 >= l_lo .and. l + 1 <= l_hi) then
                ! 2nd order central difference
                D_zz = (q_prim_vf(momxb + 2)%sf(j, k, l + 1) - q_prim_vf(momxb + 2)%sf(j, k, l - 1)) / &
                       (z_cc(l + 1) - z_cc(l - 1))
            else if (l + 1 <= l_hi) then
                ! Forward difference
                D_zz = (q_prim_vf(momxb + 2)%sf(j, k, l + 1) - q_prim_vf(momxb + 2)%sf(j, k, l)) / &
                       (z_cc(l + 1) - z_cc(l))
            else if (l - 1 >= l_lo) then
                ! Backward difference
                D_zz = (q_prim_vf(momxb + 2)%sf(j, k, l) - q_prim_vf(momxb + 2)%sf(j, k, l - 1)) / &
                       (z_cc(l) - z_cc(l - 1))
            else
                D_zz = 0._wp
            end if
        else
            D_zz = 0._wp
        end if

        ! Compute D_xy = 0.5*(du/dy + dv/dx) (2D and 3D only)
        if (n > 0) then
            if (use_4th_order_x .and. use_4th_order_y) then
                ! 4th order central difference for both terms
                coeff_4th = 1._wp/(y_cc(k - 2) - 8._wp*y_cc(k - 1) - y_cc(k + 2) + 8._wp*y_cc(k + 1))
                D_xy = 0.5_wp*coeff_4th*(q_prim_vf(momxb)%sf(j, k - 2, l) - 8._wp*q_prim_vf(momxb)%sf(j, k - 1, l) + &
                                          8._wp*q_prim_vf(momxb)%sf(j, k + 1, l) - q_prim_vf(momxb)%sf(j, k + 2, l))
                coeff_4th = 1._wp/(x_cc(j - 2) - 8._wp*x_cc(j - 1) - x_cc(j + 2) + 8._wp*x_cc(j + 1))
                D_xy = D_xy + 0.5_wp*coeff_4th*(q_prim_vf(momxb + 1)%sf(j - 2, k, l) - 8._wp*q_prim_vf(momxb + 1)%sf(j - 1, k, l) + &
                                                 8._wp*q_prim_vf(momxb + 1)%sf(j + 1, k, l) - q_prim_vf(momxb + 1)%sf(j + 2, k, l))
            else if (j - 1 >= j_lo .and. j + 1 <= j_hi .and. k - 1 >= k_lo .and. k + 1 <= k_hi) then
                ! 2nd order central difference
                D_xy = 0.5_wp*((q_prim_vf(momxb)%sf(j, k + 1, l) - q_prim_vf(momxb)%sf(j, k - 1, l)) / &
                                (y_cc(k + 1) - y_cc(k - 1)) + &
                                (q_prim_vf(momxb + 1)%sf(j + 1, k, l) - q_prim_vf(momxb + 1)%sf(j - 1, k, l)) / &
                                (x_cc(j + 1) - x_cc(j - 1)))
            else
                ! Near boundary: use available one-sided differences
                D_xy = 0._wp
                if (k - 1 >= k_lo .and. k + 1 <= k_hi) then
                    D_xy = 0.5_wp*(q_prim_vf(momxb)%sf(j, k + 1, l) - q_prim_vf(momxb)%sf(j, k - 1, l)) / &
                                   (y_cc(k + 1) - y_cc(k - 1))
                end if
                if (j - 1 >= j_lo .and. j + 1 <= j_hi) then
                    D_xy = D_xy + 0.5_wp*(q_prim_vf(momxb + 1)%sf(j + 1, k, l) - q_prim_vf(momxb + 1)%sf(j - 1, k, l)) / &
                                          (x_cc(j + 1) - x_cc(j - 1))
                end if
            end if
        else
            D_xy = 0._wp
        end if

        ! Compute D_xz = 0.5*(du/dz + dw/dx) (3D only)
        if (p > 0) then
            if (use_4th_order_x .and. use_4th_order_z) then
                ! 4th order central difference
                coeff_4th = 1._wp/(z_cc(l - 2) - 8._wp*z_cc(l - 1) - z_cc(l + 2) + 8._wp*z_cc(l + 1))
                D_xz = 0.5_wp*coeff_4th*(q_prim_vf(momxb)%sf(j, k, l - 2) - 8._wp*q_prim_vf(momxb)%sf(j, k, l - 1) + &
                                          8._wp*q_prim_vf(momxb)%sf(j, k, l + 1) - q_prim_vf(momxb)%sf(j, k, l + 2))
                coeff_4th = 1._wp/(x_cc(j - 2) - 8._wp*x_cc(j - 1) - x_cc(j + 2) + 8._wp*x_cc(j + 1))
                D_xz = D_xz + 0.5_wp*coeff_4th*(q_prim_vf(momxb + 2)%sf(j - 2, k, l) - 8._wp*q_prim_vf(momxb + 2)%sf(j - 1, k, l) + &
                                                 8._wp*q_prim_vf(momxb + 2)%sf(j + 1, k, l) - q_prim_vf(momxb + 2)%sf(j + 2, k, l))
            else if (j - 1 >= j_lo .and. j + 1 <= j_hi .and. l - 1 >= l_lo .and. l + 1 <= l_hi) then
                ! 2nd order central difference
                D_xz = 0.5_wp*((q_prim_vf(momxb)%sf(j, k, l + 1) - q_prim_vf(momxb)%sf(j, k, l - 1)) / &
                                (z_cc(l + 1) - z_cc(l - 1)) + &
                                (q_prim_vf(momxb + 2)%sf(j + 1, k, l) - q_prim_vf(momxb + 2)%sf(j - 1, k, l)) / &
                                (x_cc(j + 1) - x_cc(j - 1)))
            else
                D_xz = 0._wp
            end if
        else
            D_xz = 0._wp
        end if

        ! Compute D_yz = 0.5*(dv/dz + dw/dy) (3D only)
        if (p > 0 .and. n > 0) then
            if (use_4th_order_y .and. use_4th_order_z) then
                ! 4th order central difference
                coeff_4th = 1._wp/(z_cc(l - 2) - 8._wp*z_cc(l - 1) - z_cc(l + 2) + 8._wp*z_cc(l + 1))
                D_yz = 0.5_wp*coeff_4th*(q_prim_vf(momxb + 1)%sf(j, k, l - 2) - 8._wp*q_prim_vf(momxb + 1)%sf(j, k, l - 1) + &
                                          8._wp*q_prim_vf(momxb + 1)%sf(j, k, l + 1) - q_prim_vf(momxb + 1)%sf(j, k, l + 2))
                coeff_4th = 1._wp/(y_cc(k - 2) - 8._wp*y_cc(k - 1) - y_cc(k + 2) + 8._wp*y_cc(k + 1))
                D_yz = D_yz + 0.5_wp*coeff_4th*(q_prim_vf(momxb + 2)%sf(j, k - 2, l) - 8._wp*q_prim_vf(momxb + 2)%sf(j, k - 1, l) + &
                                                 8._wp*q_prim_vf(momxb + 2)%sf(j, k + 1, l) - q_prim_vf(momxb + 2)%sf(j, k + 2, l))
            else if (k - 1 >= k_lo .and. k + 1 <= k_hi .and. l - 1 >= l_lo .and. l + 1 <= l_hi) then
                ! 2nd order central difference
                D_yz = 0.5_wp*((q_prim_vf(momxb + 1)%sf(j, k, l + 1) - q_prim_vf(momxb + 1)%sf(j, k, l - 1)) / &
                                (z_cc(l + 1) - z_cc(l - 1)) + &
                                (q_prim_vf(momxb + 2)%sf(j, k + 1, l) - q_prim_vf(momxb + 2)%sf(j, k - 1, l)) / &
                                (y_cc(k + 1) - y_cc(k - 1)))
            else
                D_yz = 0._wp
            end if
        else
            D_yz = 0._wp
        end if

    end subroutine s_compute_velocity_gradients_at_cell

    !> Main interface: Computes Re_visc per-phase with all logic internal
    !! Works for 1D, 2D, and 3D cases
    !! Handles both Newtonian and non-Newtonian fluids automatically
    !! @param q_prim_vf Primitive variables (required for gradient computation if gradients not provided)
    !! @param alpha_visc Volume fractions (required for mixture calculation context)
    !! @param j x index
    !! @param k y index
    !! @param l z index
    !! @param grad_x_vf Optional: Pre-computed x-direction gradients (if provided, used instead of computing from q_prim_vf)
    !! @param grad_y_vf Optional: Pre-computed y-direction gradients
    !! @param grad_z_vf Optional: Pre-computed z-direction gradients
    !! @param Re_visc_per_phase Output: Re_visc_per_phase(fluid, direction) = 1/mu (inverse viscosity).
    !!                                  Use for mixture Re (CFL, stress); form mixture as 1/sum(alpha/Re_visc_per_phase).
    !!                                  Newtonian: Re input value; non-Newtonian: 1/mu from HB model.
    pure subroutine s_compute_re_visc(q_prim_vf, alpha_visc, j, k, l, Re_visc_per_phase, &
                                      grad_x_vf, grad_y_vf, grad_z_vf)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp), dimension(num_fluids), intent(in) :: alpha_visc
        integer, intent(in) :: j, k, l
        real(wp), dimension(num_fluids, 2), intent(out) :: Re_visc_per_phase
        ! Optional pre-computed velocity gradients. Passed in from callers that already
        ! have gradients available (e.g. viscous routines / Riemann solvers).
        !
        ! NOTE: Use assumed-shape to accept both (num_dims) and (num_vels) caller arrays.
        type(scalar_field), dimension(:), intent(in), optional :: grad_x_vf, grad_y_vf, grad_z_vf

        logical :: has_non_newtonian, in_bounds
        real(wp) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp) :: shear_rate, mu_fluid
        integer :: i, q, r

        ! Initialize all to default
        !$acc loop seq
        do i = 1, 2
            !$acc loop seq
            do q = 1, num_fluids
                Re_visc_per_phase(q, i) = dflt_real
            end do
        end do

        ! Check if indices are within valid array bounds
        in_bounds = (j >= idwbuff(1)%beg) .and. (j <= idwbuff(1)%end)
        if (n > 0) in_bounds = in_bounds .and. (k >= idwbuff(2)%beg) .and. (k <= idwbuff(2)%end)
        if (p > 0) in_bounds = in_bounds .and. (l >= idwbuff(3)%beg) .and. (l <= idwbuff(3)%end)

        ! If out of bounds, return default Newtonian values
        if (.not. in_bounds) then
            !$acc loop seq
            do i = 1, 2
                !$acc loop seq
                do q = 1, Re_size(i)
                    Re_visc_per_phase(Re_idx(i, q), i) = fluid_pp(Re_idx(i, q))%Re(i)
                end do
            end do
            return
        end if

        ! Check if any fluid is non-Newtonian
        has_non_newtonian = .false.
        !$acc loop seq
        do i = 1, num_fluids
            if (fluid_pp(i)%non_newtonian) then
                has_non_newtonian = .true.
                exit
            end if
        end do

        

        if (has_non_newtonian) then
            ! Non-Newtonian: need to compute velocity gradients
            if (present(grad_x_vf) .and. present(grad_y_vf) .and. present(grad_z_vf)) then
                ! Use provided gradients
                D_xx = grad_x_vf(1)%sf(j, k, l)
                if (n > 0) then
                    D_yy = grad_y_vf(2)%sf(j, k, l)
                    D_xy = 0.5_wp*(grad_y_vf(1)%sf(j, k, l) + grad_x_vf(2)%sf(j, k, l))
                else
                    D_yy = 0._wp
                    D_xy = 0._wp
                end if
                if (p > 0) then
                    D_zz = grad_z_vf(3)%sf(j, k, l)
                    D_xz = 0.5_wp*(grad_z_vf(1)%sf(j, k, l) + grad_x_vf(3)%sf(j, k, l))
                    D_yz = 0.5_wp*(grad_z_vf(2)%sf(j, k, l) + grad_y_vf(3)%sf(j, k, l))
                else
                    D_zz = 0._wp
                    D_xz = 0._wp
                    D_yz = 0._wp
                end if
            else
                ! Compute gradients from q_prim_vf
                call s_compute_velocity_gradients_at_cell(q_prim_vf, j, k, l, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
            end if

            ! Compute shear rate from gradients
            shear_rate = f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)

            ! For each phase, compute Re_visc
            ! i=1 is shear viscosity, i=2 is bulk viscosity
            !$acc loop seq
            do q = 1, num_fluids
                if (fluid_pp(q)%non_newtonian) then
                    ! Non-Newtonian: compute shear mu from HB model with Papanastasiou regularization
                    mu_fluid = f_compute_hb_viscosity( &
                        fluid_pp(q)%tau0, &
                        fluid_pp(q)%K, &
                        fluid_pp(q)%n, &
                        fluid_pp(q)%mu_min, &
                        fluid_pp(q)%mu_max, &
                        shear_rate, &
                        fluid_pp(q)%hb_m)
                    ! Shear Reynolds number = 1/mu_shear
                    Re_visc_per_phase(q, 1) = 1._wp/mu_fluid
                    ! Bulk Reynolds number = 1/mu_bulk (use user-specified mu_bulk)
                    ! If mu_bulk is 0, set Re_bulk to a very large value (no bulk viscosity effect)
                    if (fluid_pp(q)%mu_bulk > 0._wp) then
                        Re_visc_per_phase(q, 2) = 1._wp/fluid_pp(q)%mu_bulk
                    else
                        Re_visc_per_phase(q, 2) = dflt_real
                    end if
                else
                    ! Newtonian: return Re input values
                    !$acc loop seq
                    do i = 1, 2
                        if (Re_size(i) > 0) then
                            !$acc loop seq
                            do r = 1, Re_size(i)
                                if (Re_idx(i, r) == q) then
                                    Re_visc_per_phase(q, i) = fluid_pp(q)%Re(i)
                                    exit
                                end if
                            end do
                        end if
                    end do
                end if
            end do
        else
            ! All Newtonian: return Re input values per-phase
            !$acc loop seq
            do i = 1, 2
                !$acc loop seq
                do q = 1, Re_size(i)
                    Re_visc_per_phase(Re_idx(i, q), i) = fluid_pp(Re_idx(i, q))%Re(i)
                end do
            end do
        end if

    end subroutine s_compute_re_visc

    !> Computes per-phase viscosity (mu) at a single cell (j,k,l)
    !! Returns mu_shear_phases, mu_bulk_phases = mu (viscosity). Use for dmu/dx, face reconstruction in m_viscous.
    !! For non-Newtonian phases: uses Herschel-Bulkley model; for Newtonian: mu = 1/Re.
    !! @param q_prim_vf Primitive variables
    !! @param j x index
    !! @param k y index
    !! @param l z index
    !! @param mu_shear_phases Output: Shear viscosity for each phase
    !! @param mu_bulk_phases Output: Bulk viscosity for each phase
    !! @param grad_x_vf Optional: Pre-computed x-direction gradients
    !! @param grad_y_vf Optional: Pre-computed y-direction gradients
    !! @param grad_z_vf Optional: Pre-computed z-direction gradients
    pure subroutine s_compute_phase_viscosity_at_cell(q_prim_vf, j, k, l, mu_shear_phases, mu_bulk_phases, &
                                                       grad_x_vf, grad_y_vf, grad_z_vf)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in) :: j, k, l
        real(wp), dimension(num_fluids), intent(out) :: mu_shear_phases, mu_bulk_phases
        type(scalar_field), dimension(:), intent(in), optional :: grad_x_vf, grad_y_vf, grad_z_vf

        logical :: has_non_newtonian, in_bounds
        real(wp) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp) :: shear_rate
        integer :: i, q, r

        ! Initialize all phases to zero
        !$acc loop seq
        do q = 1, num_fluids
            mu_shear_phases(q) = 0._wp
            mu_bulk_phases(q) = 0._wp
        end do

        ! Check if indices are within valid array bounds
        in_bounds = (j >= idwbuff(1)%beg) .and. (j <= idwbuff(1)%end)
        if (n > 0) in_bounds = in_bounds .and. (k >= idwbuff(2)%beg) .and. (k <= idwbuff(2)%end)
        if (p > 0) in_bounds = in_bounds .and. (l >= idwbuff(3)%beg) .and. (l <= idwbuff(3)%end)

        ! If out of bounds, return Newtonian values
        if (.not. in_bounds) then
            !$acc loop seq
            do q = 1, num_fluids
                !$acc loop seq
                do r = 1, Re_size(1)
                    if (Re_idx(1, r) == q) then
                        if (fluid_pp(q)%Re(1) > sgm_eps) then
                            mu_shear_phases(q) = 1._wp/fluid_pp(q)%Re(1)
                        end if
                        exit
                    end if
                end do
                !$acc loop seq
                do r = 1, Re_size(2)
                    if (Re_idx(2, r) == q) then
                        if (fluid_pp(q)%Re(2) > sgm_eps) then
                            mu_bulk_phases(q) = 1._wp/fluid_pp(q)%Re(2)
                        end if
                        exit
                    end if
                end do
            end do
            return
        end if

        ! Check if any fluid is non-Newtonian
        has_non_newtonian = .false.
        !$acc loop seq
        do i = 1, num_fluids
            if (fluid_pp(i)%non_newtonian) then
                has_non_newtonian = .true.
                exit
            end if
        end do

        ! Compute velocity gradients (needed for non-Newtonian shear rate)
        if (has_non_newtonian) then
            if (present(grad_x_vf) .and. present(grad_y_vf) .and. present(grad_z_vf)) then
                ! Use provided gradients
                D_xx = grad_x_vf(1)%sf(j, k, l)
                if (n > 0) then
                    D_yy = grad_y_vf(2)%sf(j, k, l)
                    D_xy = 0.5_wp*(grad_y_vf(1)%sf(j, k, l) + grad_x_vf(2)%sf(j, k, l))
                else
                    D_yy = 0._wp
                    D_xy = 0._wp
                end if
                if (p > 0) then
                    D_zz = grad_z_vf(3)%sf(j, k, l)
                    D_xz = 0.5_wp*(grad_z_vf(1)%sf(j, k, l) + grad_x_vf(3)%sf(j, k, l))
                    D_yz = 0.5_wp*(grad_z_vf(2)%sf(j, k, l) + grad_y_vf(3)%sf(j, k, l))
                else
                    D_zz = 0._wp
                    D_xz = 0._wp
                    D_yz = 0._wp
                end if
            else
                ! Compute gradients from q_prim_vf
                call s_compute_velocity_gradients_at_cell(q_prim_vf, j, k, l, D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
            end if

            ! Compute shear rate from gradients
            shear_rate = f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz)
        else
            shear_rate = 0._wp
        end if

        ! Compute viscosity for each phase
        !$acc loop seq
        do q = 1, num_fluids
            if (fluid_pp(q)%non_newtonian) then
                ! Non-Newtonian: compute shear mu from HB model with Papanastasiou regularization
                mu_shear_phases(q) = f_compute_hb_viscosity( &
                    fluid_pp(q)%tau0, &
                    fluid_pp(q)%K, &
                    fluid_pp(q)%n, &
                    fluid_pp(q)%mu_min, &
                    fluid_pp(q)%mu_max, &
                    shear_rate, &
                    fluid_pp(q)%hb_m)
                ! Bulk viscosity
                if (fluid_pp(q)%mu_bulk > 0._wp) then
                    mu_bulk_phases(q) = fluid_pp(q)%mu_bulk
                else
                    mu_bulk_phases(q) = 0._wp
                end if
            else
                ! Newtonian: get mu from Re (mu = 1/Re)
                !$acc loop seq
                do r = 1, Re_size(1)
                    if (Re_idx(1, r) == q) then
                        if (fluid_pp(q)%Re(1) > sgm_eps) then
                            mu_shear_phases(q) = 1._wp/fluid_pp(q)%Re(1)
                        end if
                        exit
                    end if
                end do
                !$acc loop seq
                do r = 1, Re_size(2)
                    if (Re_idx(2, r) == q) then
                        if (fluid_pp(q)%Re(2) > sgm_eps) then
                            mu_bulk_phases(q) = 1._wp/fluid_pp(q)%Re(2)
                        end if
                        exit
                    end if
                end do
            end if
        end do

    end subroutine s_compute_phase_viscosity_at_cell

end module m_re_visc
