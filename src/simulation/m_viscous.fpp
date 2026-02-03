!>
!! @file m_viscous.f90
!! @brief Contains module m_viscous
#:include 'macros.fpp'

!> @brief The module contains the subroutines used to compute viscous terms.
module m_viscous

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_weno

    use m_helper

    use m_finite_differences

    use m_hb_function          !< Herschel-Bulkley non-Newtonian viscosity
    use m_re_visc              !< Re_visc computation module

    private; public s_get_viscous, &
 s_compute_viscous_stress_tensor, &
 s_initialize_viscous_module, &
 s_reconstruct_cell_boundary_values_visc_deriv, &
 s_compute_viscous_gradient_correction, &
 s_finalize_viscous_module

    type(int_bounds_info) :: iv
    type(int_bounds_info) :: is1_viscous, is2_viscous, is3_viscous
    !$acc declare create(is1_viscous, is2_viscous, is3_viscous, iv)

    !real(wp), allocatable, dimension(:, :) :: Res_viscous
    !$acc declare create(Res_viscous)

contains

    impure subroutine s_initialize_viscous_module

        ! Note: Re_visc is now computed dynamically using m_re_visc module
        ! No need to pre-compute Res_viscous since it depends on velocity gradients
        ! which are not available at initialization time (especially for non-Newtonian fluids)
        !
        ! Previous implementation (commented out):
        !     integer :: i, j !< generic loop iterators
        !     @:ALLOCATE(Res_viscous(1:2, 1:maxval(Re_size)))
        !     do i = 1, 2
        !         do j = 1, Re_size(i)
        !             Res_viscous(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
        !         end do
        !     end do
        !     !$acc update device(Res_viscous, Re_idx, Re_size)
        !     !$acc enter data copyin(is1_viscous, is2_viscous, is3_viscous, iv)
        
        !$acc update device(Re_idx, Re_size)
        !$acc enter data copyin(is1_viscous, is2_viscous, is3_viscous, iv)
    end subroutine s_initialize_viscous_module

    !> The purpose of this subroutine is to compute the viscous
    !      stress tensor for the cells directly next to the axis in
    !      cylindrical coordinates. This is necessary to avoid the
    !      1/r singularity that arises at the cell boundary coinciding
    !      with the axis, i.e., y_cb(-1) = 0.
    !  @param q_prim_vf Cell-average primitive variables
    !  @param grad_x_vf Cell-average primitive variable derivatives, x-dir
    !  @param grad_y_vf Cell-average primitive variable derivatives, y-dir
    !  @param grad_z_vf Cell-average primitive variable derivatives, z-dir
    subroutine s_compute_viscous_stress_tensor(q_prim_vf, grad_x_vf, grad_y_vf, grad_z_vf, &
                                               tau_Re_vf, &
                                               ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(num_dims), intent(in) :: grad_x_vf, grad_y_vf, grad_z_vf
        type(scalar_field), dimension(1:sys_size), intent(inout) :: tau_Re_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz

        real(wp) :: rho_visc, gamma_visc, pi_inf_visc, alpha_visc_sum  !< Mixture variables
        real(wp), dimension(2) :: Re_visc
        real(wp), dimension(num_fluids) :: alpha_visc, alpha_rho_visc
        real(wp), dimension(num_fluids, 2) :: Re_visc_per_phase

        real(wp), dimension(num_dims, num_dims) :: tau_Re

        integer :: i, j, k, l, q !< Generic loop iterator

        is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

        !$acc update device(is1_viscous, is2_viscous, is3_viscous)

        !$acc parallel loop collapse(3) gang vector default(present)
        do l = is3_viscous%beg, is3_viscous%end
            do k = is2_viscous%beg, is2_viscous%end
                do j = is1_viscous%beg, is1_viscous%end
                    !$acc loop seq
                    do i = momxb, E_idx
                        tau_Re_vf(i)%sf(j, k, l) = 0._wp
                    end do
                end do
            end do
        end do
        if (shear_stress) then    ! Shear stresses
            !$acc parallel loop collapse(3) gang vector default(present) private(alpha_visc, alpha_rho_visc, Re_visc, tau_Re )
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                ! Compute Re_visc per-phase using the unified module
                                ! Pass all necessary inputs: q_prim_vf, alpha_visc, indices, and pre-computed gradients
                                call s_compute_re_visc(q_prim_vf, alpha_visc, j, k, l, Re_visc_per_phase, &
                                                       grad_x_vf, grad_y_vf, grad_z_vf)

                                ! Compute mixture Re_visc from per-phase values
                                ! Re_visc(i) = 1 / sum(alpha_visc(q) / Re_visc_per_phase(q, i))
                                ! Only sum over phases that have valid Re_visc_per_phase values
                                !$acc loop seq
                                do i = 1, 2
                                    Re_visc(i) = 0._wp
                                    !$acc loop seq
                                    do q = 1, num_fluids
                                        if (Re_visc_per_phase(q, i) /= dflt_real .and. Re_visc_per_phase(q, i) > sgm_eps) then
                                            Re_visc(i) = Re_visc(i) + alpha_visc(q)/Re_visc_per_phase(q, i)
                                        end if
                                    end do
                                    
                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)
                                   
                                
                                end do
                                
                            end if
                        end if

                        tau_Re(2, 1) = (grad_y_vf(1)%sf(j, k, l) + &
                                        grad_x_vf(2)%sf(j, k, l))/ &
                                       Re_visc(1)

                        tau_Re(2, 2) = (4._wp*grad_y_vf(2)%sf(j, k, l) &
                                        - 2._wp*grad_x_vf(1)%sf(j, k, l) &
                                        - 2._wp*q_prim_vf(momxb + 1)%sf(j, k, l)/y_cc(k))/ &
                                       (3._wp*Re_visc(1))
                        !$acc loop seq
                        do i = 1, 2
                            tau_Re_vf(contxe + i)%sf(j, k, l) = &
                                tau_Re_vf(contxe + i)%sf(j, k, l) - &
                                tau_Re(2, i)

                            tau_Re_vf(E_idx)%sf(j, k, l) = &
                                tau_Re_vf(E_idx)%sf(j, k, l) - &
                                q_prim_vf(contxe + i)%sf(j, k, l)*tau_Re(2, i)
                        end do
                    end do
                end do
            end do
        end if

        if (bulk_stress) then    ! Bulk stresses
            !$acc parallel loop collapse(3) gang vector default(present) private(alpha_visc, alpha_rho_visc, Re_visc, tau_Re )
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                ! Compute Re_visc per-phase using the unified module
                                ! Pass all necessary inputs: q_prim_vf, alpha_visc, indices, and pre-computed gradients
                                call s_compute_re_visc(q_prim_vf, alpha_visc, j, k, l, Re_visc_per_phase, &
                                                       grad_x_vf, grad_y_vf, grad_z_vf)

                                ! Compute mixture Re_visc from per-phase values
                                ! Re_visc(i) = 1 / sum(alpha_visc(q) / Re_visc_per_phase(q, i))
                                ! Only sum over phases that have valid Re_visc_per_phase values
                                !$acc loop seq
                                do i = 1, 2
                                    Re_visc(i) = 0._wp
                                    !$acc loop seq
                                    do q = 1, num_fluids
                                        if (Re_visc_per_phase(q, i) /= dflt_real .and. Re_visc_per_phase(q, i) > sgm_eps) then
                                            Re_visc(i) = Re_visc(i) + alpha_visc(q)/Re_visc_per_phase(q, i)
                                        end if
                                    end do
                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)
                                end do
                            end if
                        end if

                        tau_Re(2, 2) = (grad_x_vf(1)%sf(j, k, l) + &
                                        grad_y_vf(2)%sf(j, k, l) + &
                                        q_prim_vf(momxb + 1)%sf(j, k, l)/y_cc(k))/ &
                                       Re_visc(2)

                        tau_Re_vf(momxb + 1)%sf(j, k, l) = &
                            tau_Re_vf(momxb + 1)%sf(j, k, l) - &
                            tau_Re(2, 2)

                        tau_Re_vf(E_idx)%sf(j, k, l) = &
                            tau_Re_vf(E_idx)%sf(j, k, l) - &
                            q_prim_vf(momxb + 1)%sf(j, k, l)*tau_Re(2, 2)

                    end do
                end do
            end do
        end if

        if (p == 0) return

        if (shear_stress) then    ! Shear stresses
            !$acc parallel loop collapse(3) gang vector default(present) private(alpha_visc, alpha_rho_visc, Re_visc, Re_visc_per_phase, tau_Re )
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                ! Compute Re_visc per-phase using the unified module
                                ! Pass all necessary inputs: q_prim_vf, alpha_visc, indices, and pre-computed gradients
                                call s_compute_re_visc(q_prim_vf, alpha_visc, j, k, l, Re_visc_per_phase, &
                                                       grad_x_vf, grad_y_vf, grad_z_vf)

                                ! Compute mixture Re_visc from per-phase values
                                ! Re_visc(i) = 1 / sum(alpha_visc(q) / Re_visc_per_phase(q, i))
                                ! Only sum over phases that have valid Re_visc_per_phase values
                                !$acc loop seq
                                do i = 1, 2
                                    Re_visc(i) = 0._wp
                                    !$acc loop seq
                                    do q = 1, num_fluids
                                        if (Re_visc_per_phase(q, i) /= dflt_real .and. Re_visc_per_phase(q, i) > sgm_eps) then
                                            Re_visc(i) = Re_visc(i) + alpha_visc(q)/Re_visc_per_phase(q, i)
                                        end if
                                    end do
                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)
                                end do
                            end if
                        end if

                        tau_Re(2, 2) = -(2._wp/3._wp)*grad_z_vf(3)%sf(j, k, l)/y_cc(k)/ &
                                       Re_visc(1)

                        tau_Re(2, 3) = ((grad_z_vf(2)%sf(j, k, l) - &
                                         q_prim_vf(momxe)%sf(j, k, l))/ &
                                        y_cc(k) + grad_y_vf(3)%sf(j, k, l))/ &
                                       Re_visc(1)

                        !$acc loop seq
                        do i = 2, 3
                            tau_Re_vf(contxe + i)%sf(j, k, l) = &
                                tau_Re_vf(contxe + i)%sf(j, k, l) - &
                                tau_Re(2, i)

                            tau_Re_vf(E_idx)%sf(j, k, l) = &
                                tau_Re_vf(E_idx)%sf(j, k, l) - &
                                q_prim_vf(contxe + i)%sf(j, k, l)*tau_Re(2, i)
                        end do

                    end do
                end do
            end do
        end if

        if (bulk_stress) then    ! Bulk stresses
            !$acc parallel loop collapse(3) gang vector default(present) private(alpha_visc, alpha_rho_visc, Re_visc, Re_visc_per_phase, tau_Re )
            do l = is3_viscous%beg, is3_viscous%end
                do k = -1, 1
                    do j = is1_viscous%beg, is1_viscous%end

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho_visc(i) = q_prim_vf(i)%sf(j, k, l)
                            if (bubbles_euler .and. num_fluids == 1) then
                                alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                            else
                                alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                            end if
                        end do

                        if (bubbles_euler) then
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    rho_visc = rho_visc + alpha_rho_visc(i)
                                    gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                    pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                                end do
                            else
                                rho_visc = alpha_rho_visc(1)
                                gamma_visc = gammas(1)
                                pi_inf_visc = pi_infs(1)
                            end if
                        else
                            rho_visc = 0._wp
                            gamma_visc = 0._wp
                            pi_inf_visc = 0._wp

                            alpha_visc_sum = 0._wp

                            if (mpp_lim) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_visc(i) = max(0._wp, alpha_rho_visc(i))
                                    alpha_visc(i) = min(max(0._wp, alpha_visc(i)), 1._wp)
                                    alpha_visc_sum = alpha_visc_sum + alpha_visc(i)
                                end do

                                alpha_visc = alpha_visc/max(alpha_visc_sum, sgm_eps)

                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                rho_visc = rho_visc + alpha_rho_visc(i)
                                gamma_visc = gamma_visc + alpha_visc(i)*gammas(i)
                                pi_inf_visc = pi_inf_visc + alpha_visc(i)*pi_infs(i)
                            end do

                            if (viscous) then
                                ! Compute Re_visc per-phase using the unified module
                                ! Pass all necessary inputs: q_prim_vf, alpha_visc, indices, and pre-computed gradients
                                call s_compute_re_visc(q_prim_vf, alpha_visc, j, k, l, Re_visc_per_phase, &
                                                       grad_x_vf, grad_y_vf, grad_z_vf)

                                ! Compute mixture Re_visc from per-phase values
                                ! Re_visc(i) = 1 / sum(alpha_visc(q) / Re_visc_per_phase(q, i))
                                ! Only sum over phases that have valid Re_visc_per_phase values
                                !$acc loop seq
                                do i = 1, 2
                                    Re_visc(i) = 0._wp
                                    !$acc loop seq
                                    do q = 1, num_fluids
                                        if (Re_visc_per_phase(q, i) /= dflt_real .and. Re_visc_per_phase(q, i) > sgm_eps) then
                                            Re_visc(i) = Re_visc(i) + alpha_visc(q)/Re_visc_per_phase(q, i)
                                        end if
                                    end do
                                    Re_visc(i) = 1._wp/max(Re_visc(i), sgm_eps)
                                end do
                            end if
                        end if

                        tau_Re(2, 2) = grad_z_vf(3)%sf(j, k, l)/y_cc(k)/ &
                                       Re_visc(2)

                        tau_Re_vf(momxb + 1)%sf(j, k, l) = &
                            tau_Re_vf(momxb + 1)%sf(j, k, l) - &
                            tau_Re(2, 2)

                        tau_Re_vf(E_idx)%sf(j, k, l) = &
                            tau_Re_vf(E_idx)%sf(j, k, l) - &
                            q_prim_vf(momxb + 1)%sf(j, k, l)*tau_Re(2, 2)

                    end do
                end do
            end do
        end if
    end subroutine s_compute_viscous_stress_tensor

    !>  Computes viscous terms
    !!  @param q_cons_vf Cell-averaged conservative variables
    !!  @param q_prim_vf Cell-averaged primitive variables
    !!  @param rhs_vf Cell-averaged RHS variables
    subroutine s_get_viscous(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                             dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, &
                             qL_prim, &
                             qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                             dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n, &
                             qR_prim, &
                             q_prim_qp, &
                             dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, &
                             ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), &
            intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf, &
                             qL_prim_rsy_vf, qR_prim_rsy_vf, &
                             qL_prim_rsz_vf, qR_prim_rsz_vf

        type(vector_field), dimension(num_dims), intent(inout) :: qL_prim, qR_prim

        type(vector_field), intent(in) :: q_prim_qp

        type(vector_field), dimension(1:num_dims), &
            intent(inout) :: dqL_prim_dx_n, dqR_prim_dx_n, &
                             dqL_prim_dy_n, dqR_prim_dy_n, &
                             dqL_prim_dz_n, dqR_prim_dz_n

        type(vector_field), dimension(1), intent(inout) :: dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, l

        do i = 1, num_dims

            iv%beg = mom_idx%beg; iv%end = mom_idx%end

            !$acc update device(iv)

            call s_reconstruct_cell_boundary_values_visc( &
                q_prim_qp%vf(iv%beg:iv%end), &
                qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, &
                qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, &
                i, qL_prim(i)%vf(iv%beg:iv%end), qR_prim(i)%vf(iv%beg:iv%end), &
                ix, iy, iz)
        end do

        if (weno_Re_flux) then
            ! Compute velocity gradient at cell centers using scalar
            ! divergence theorem
            do i = 1, num_dims
                if (i == 1) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dx_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dx, m, buff_size)
                elseif (i == 2) then
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dy_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dy, n, buff_size)
                else
                    call s_apply_scalar_divergence_theorem( &
                        qL_prim(i)%vf(iv%beg:iv%end), &
                        qR_prim(i)%vf(iv%beg:iv%end), &
                        dq_prim_dz_qp(1)%vf(iv%beg:iv%end), i, &
                        ix, iy, iz, iv, dz, p, buff_size)
                end if
            end do

        else ! Compute velocity gradient at cell centers using finite differences

            iv%beg = mom_idx%beg; iv%end = mom_idx%end
            !$acc update device(iv)

            is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

            !$acc update device(is1_viscous, is2_viscous, is3_viscous)

            ! For non-Newtonian fluids with sufficient buffer, use 4th order FD
            ! For Newtonian fluids, use original 1st order one-sided differences
            if (any_non_newtonian .and. buff_size >= 2) then
                ! 4th order: use 4-point stencil for face gradients
                ! dqL at face j-1/2 uses cells j-2, j-1, j, j+1
                ! Formula: (-f_{j-2} + 27*f_{j-1} - 27*f_{j} + f_{j+1}) / (24*h) for uniform grid
                ! For non-uniform, we use: (f_{j} - f_{j-1}) / (x_j - x_{j-1}) as base
                ! and add higher-order correction

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = iy%beg, iy%end
                        do j = is1_viscous%beg + 1, is1_viscous%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                if (j > is1_viscous%beg + 1 .and. j < is1_viscous%end) then
                                    ! Interior: 4th order using 4-point stencil
                                    ! Approximate derivative at face using cubic interpolation
                                    dqL_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                        (-q_prim_qp%vf(i)%sf(j - 2, k, l) + &
                                         27._wp*q_prim_qp%vf(i)%sf(j - 1, k, l) - &
                                         27._wp*q_prim_qp%vf(i)%sf(j, k, l) + &
                                         q_prim_qp%vf(i)%sf(j + 1, k, l))/ &
                                        (24._wp*(x_cc(j) - x_cc(j - 1)))
                                else
                                    ! Near boundary: fall back to 1st order
                                    dqL_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                        (q_prim_qp%vf(i)%sf(j, k, l) - &
                                         q_prim_qp%vf(i)%sf(j - 1, k, l))/ &
                                        (x_cc(j) - x_cc(j - 1))
                                end if
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg, is2_viscous%end
                        do j = is1_viscous%beg, is1_viscous%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                if (j > is1_viscous%beg .and. j < is1_viscous%end - 1) then
                                    ! Interior: 4th order
                                    dqR_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                        (-q_prim_qp%vf(i)%sf(j - 1, k, l) + &
                                         27._wp*q_prim_qp%vf(i)%sf(j, k, l) - &
                                         27._wp*q_prim_qp%vf(i)%sf(j + 1, k, l) + &
                                         q_prim_qp%vf(i)%sf(j + 2, k, l))/ &
                                        (24._wp*(x_cc(j + 1) - x_cc(j)))
                                else
                                    ! Near boundary: fall back to 1st order
                                    dqR_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                        (q_prim_qp%vf(i)%sf(j + 1, k, l) - &
                                         q_prim_qp%vf(i)%sf(j, k, l))/ &
                                        (x_cc(j + 1) - x_cc(j))
                                end if
                            end do
                        end do
                    end do
                end do
            else
                ! Newtonian: original 1st order one-sided differences
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = iy%beg, iy%end
                        do j = is1_viscous%beg + 1, is1_viscous%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                    (q_prim_qp%vf(i)%sf(j, k, l) - &
                                     q_prim_qp%vf(i)%sf(j - 1, k, l))/ &
                                    (x_cc(j) - x_cc(j - 1))
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg, is2_viscous%end
                        do j = is1_viscous%beg, is1_viscous%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dx_n(1)%vf(i)%sf(j, k, l) = &
                                    (q_prim_qp%vf(i)%sf(j + 1, k, l) - &
                                     q_prim_qp%vf(i)%sf(j, k, l))/ &
                                    (x_cc(j + 1) - x_cc(j))
                            end do
                        end do
                    end do
                end do
            end if

            if (n > 0) then

                if (any_non_newtonian .and. buff_size >= 2) then
                    ! 4th order for y-direction
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg + 1, is2_viscous%end
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end
                                    if (j > is2_viscous%beg + 1 .and. j < is2_viscous%end) then
                                        ! Interior: 4th order
                                        dqL_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                            (-q_prim_qp%vf(i)%sf(k, j - 2, l) + &
                                             27._wp*q_prim_qp%vf(i)%sf(k, j - 1, l) - &
                                             27._wp*q_prim_qp%vf(i)%sf(k, j, l) + &
                                             q_prim_qp%vf(i)%sf(k, j + 1, l))/ &
                                            (24._wp*(y_cc(j) - y_cc(j - 1)))
                                    else
                                        ! Near boundary: fall back to 1st order
                                        dqL_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                            (q_prim_qp%vf(i)%sf(k, j, l) - &
                                             q_prim_qp%vf(i)%sf(k, j - 1, l))/ &
                                            (y_cc(j) - y_cc(j - 1))
                                    end if
                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end
                                    if (j > is2_viscous%beg .and. j < is2_viscous%end - 1) then
                                        ! Interior: 4th order
                                        dqR_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                            (-q_prim_qp%vf(i)%sf(k, j - 1, l) + &
                                             27._wp*q_prim_qp%vf(i)%sf(k, j, l) - &
                                             27._wp*q_prim_qp%vf(i)%sf(k, j + 1, l) + &
                                             q_prim_qp%vf(i)%sf(k, j + 2, l))/ &
                                            (24._wp*(y_cc(j + 1) - y_cc(j)))
                                    else
                                        ! Near boundary: fall back to 1st order
                                        dqR_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                            (q_prim_qp%vf(i)%sf(k, j + 1, l) - &
                                             q_prim_qp%vf(i)%sf(k, j, l))/ &
                                            (y_cc(j + 1) - y_cc(j))
                                    end if
                                end do
                            end do
                        end do
                    end do
                else
                    ! Newtonian: original 1st order
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg + 1, is2_viscous%end
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end
                                    dqL_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                        (q_prim_qp%vf(i)%sf(k, j, l) - &
                                         q_prim_qp%vf(i)%sf(k, j - 1, l))/ &
                                        (y_cc(j) - y_cc(j - 1))
                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg, is3_viscous%end
                        do j = is2_viscous%beg, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end
                                    dqR_prim_dy_n(2)%vf(i)%sf(k, j, l) = &
                                        (q_prim_qp%vf(i)%sf(k, j + 1, l) - &
                                         q_prim_qp%vf(i)%sf(k, j, l))/ &
                                        (y_cc(j + 1) - y_cc(j))
                                end do
                            end do
                        end do
                    end do
                end if

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do j = is2_viscous%beg + 1, is2_viscous%end
                        do k = is1_viscous%beg + 1, is1_viscous%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqL_prim_dx_n(1)%vf(i)%sf(k, j - 1, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j - 1, l))

                                dqL_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                     dqL_prim_dx_n(2)%vf(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do j = is2_viscous%beg, is2_viscous%end - 1
                        do k = is1_viscous%beg + 1, is1_viscous%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = &
                                    (dqL_prim_dx_n(1)%vf(i)%sf(k, j + 1, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j + 1, l) + &
                                     dqL_prim_dx_n(1)%vf(i)%sf(k, j, l) + &
                                     dqR_prim_dx_n(1)%vf(i)%sf(k, j, l))

                                dqR_prim_dx_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                     dqR_prim_dx_n(2)%vf(i)%sf(k, j, l)

                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg + 1, is2_viscous%end - 1
                        do j = is1_viscous%beg + 1, is1_viscous%end
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqL_prim_dy_n(2)%vf(i)%sf(j - 1, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j - 1, k, l))

                                dqL_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                     dqL_prim_dy_n(1)%vf(i)%sf(j, k, l)

                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg + 1, is2_viscous%end - 1
                        do j = is1_viscous%beg, is1_viscous%end - 1
                            !$acc loop seq
                            do i = iv%beg, iv%end
                                dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = &
                                    (dqL_prim_dy_n(2)%vf(i)%sf(j + 1, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j + 1, k, l) + &
                                     dqL_prim_dy_n(2)%vf(i)%sf(j, k, l) + &
                                     dqR_prim_dy_n(2)%vf(i)%sf(j, k, l))

                                dqR_prim_dy_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                     dqR_prim_dy_n(1)%vf(i)%sf(j, k, l)

                            end do
                        end do
                    end do
                end do

                if (p > 0) then

                    if (any_non_newtonian .and. buff_size >= 2) then
                        ! 4th order for z-direction
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do j = is3_viscous%beg + 1, is3_viscous%end
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg, is1_viscous%end
                                    !$acc loop seq
                                    do i = iv%beg, iv%end
                                        if (j > is3_viscous%beg + 1 .and. j < is3_viscous%end) then
                                            ! Interior: 4th order
                                            dqL_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                                (-q_prim_qp%vf(i)%sf(k, l, j - 2) + &
                                                 27._wp*q_prim_qp%vf(i)%sf(k, l, j - 1) - &
                                                 27._wp*q_prim_qp%vf(i)%sf(k, l, j) + &
                                                 q_prim_qp%vf(i)%sf(k, l, j + 1))/ &
                                                (24._wp*(z_cc(j) - z_cc(j - 1)))
                                        else
                                            ! Near boundary: fall back to 1st order
                                            dqL_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                                (q_prim_qp%vf(i)%sf(k, l, j) - &
                                                 q_prim_qp%vf(i)%sf(k, l, j - 1))/ &
                                                (z_cc(j) - z_cc(j - 1))
                                        end if
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do j = is3_viscous%beg, is3_viscous%end - 1
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg, is1_viscous%end
                                    !$acc loop seq
                                    do i = iv%beg, iv%end
                                        if (j > is3_viscous%beg .and. j < is3_viscous%end - 1) then
                                            ! Interior: 4th order
                                            dqR_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                                (-q_prim_qp%vf(i)%sf(k, l, j - 1) + &
                                                 27._wp*q_prim_qp%vf(i)%sf(k, l, j) - &
                                                 27._wp*q_prim_qp%vf(i)%sf(k, l, j + 1) + &
                                                 q_prim_qp%vf(i)%sf(k, l, j + 2))/ &
                                                (24._wp*(z_cc(j + 1) - z_cc(j)))
                                        else
                                            ! Near boundary: fall back to 1st order
                                            dqR_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                                (q_prim_qp%vf(i)%sf(k, l, j + 1) - &
                                                 q_prim_qp%vf(i)%sf(k, l, j))/ &
                                                (z_cc(j + 1) - z_cc(j))
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    else
                        ! Newtonian: original 1st order
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do j = is3_viscous%beg + 1, is3_viscous%end
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg, is1_viscous%end
                                    !$acc loop seq
                                    do i = iv%beg, iv%end

                                        dqL_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                            (q_prim_qp%vf(i)%sf(k, l, j) - &
                                             q_prim_qp%vf(i)%sf(k, l, j - 1))/ &
                                            (z_cc(j) - z_cc(j - 1))
                                    end do
                                end do
                            end do
                        end do

                        !$acc parallel loop collapse(3) gang vector default(present)
                        do j = is3_viscous%beg, is3_viscous%end - 1
                            do l = is2_viscous%beg, is2_viscous%end
                                do k = is1_viscous%beg, is1_viscous%end
                                    !$acc loop seq
                                    do i = iv%beg, iv%end

                                        dqR_prim_dz_n(3)%vf(i)%sf(k, l, j) = &
                                            (q_prim_qp%vf(i)%sf(k, l, j + 1) - &
                                             q_prim_qp%vf(i)%sf(k, l, j))/ &
                                            (z_cc(j + 1) - z_cc(j))
                                    end do
                                end do
                            end do
                        end do
                    end if

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do k = is2_viscous%beg, is2_viscous%end
                            do j = is1_viscous%beg + 1, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(j - 1, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j - 1, k, l))

                                    dqL_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                         dqL_prim_dz_n(1)%vf(i)%sf(j, k, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do k = is2_viscous%beg, is2_viscous%end
                            do j = is1_viscous%beg, is1_viscous%end - 1
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(j + 1, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j + 1, k, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(j, k, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(j, k, l))

                                    dqR_prim_dz_n(1)%vf(i)%sf(j, k, l) = 25.e-2_wp* &
                                                                         dqR_prim_dz_n(1)%vf(i)%sf(j, k, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do j = is2_viscous%beg + 1, is2_viscous%end
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(k, j - 1, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j - 1, l))

                                    dqL_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                         dqL_prim_dz_n(2)%vf(i)%sf(k, j, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do l = is3_viscous%beg + 1, is3_viscous%end - 1
                        do j = is2_viscous%beg, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = &
                                        (dqL_prim_dz_n(3)%vf(i)%sf(k, j + 1, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j + 1, l) + &
                                         dqL_prim_dz_n(3)%vf(i)%sf(k, j, l) + &
                                         dqR_prim_dz_n(3)%vf(i)%sf(k, j, l))

                                    dqR_prim_dz_n(2)%vf(i)%sf(k, j, l) = 25.e-2_wp* &
                                                                         dqR_prim_dz_n(2)%vf(i)%sf(k, j, l)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_viscous%beg + 1, is3_viscous%end
                        do l = is2_viscous%beg + 1, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqL_prim_dy_n(2)%vf(i)%sf(k, l, j - 1) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j - 1))

                                    dqL_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqL_prim_dy_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_viscous%beg, is3_viscous%end - 1
                        do l = is2_viscous%beg + 1, is2_viscous%end - 1
                            do k = is1_viscous%beg, is1_viscous%end
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dy_n(2)%vf(i)%sf(k, l, j + 1) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j + 1) + &
                                         dqL_prim_dy_n(2)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dy_n(2)%vf(i)%sf(k, l, j))

                                    dqR_prim_dy_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqR_prim_dy_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_viscous%beg + 1, is3_viscous%end
                        do l = is2_viscous%beg, is2_viscous%end
                            do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                !$acc loop seq
                                do i = iv%beg, iv%end

                                    dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqL_prim_dx_n(1)%vf(i)%sf(k, l, j - 1) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j - 1))

                                    dqL_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqL_prim_dx_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do j = is3_viscous%beg, is3_viscous%end - 1
                        do l = is2_viscous%beg, is2_viscous%end
                            do k = is1_viscous%beg + 1, is1_viscous%end - 1
                                !$acc loop seq
                                do i = iv%beg, iv%end
                                    dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = &
                                        (dqL_prim_dx_n(1)%vf(i)%sf(k, l, j + 1) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j + 1) + &
                                         dqL_prim_dx_n(1)%vf(i)%sf(k, l, j) + &
                                         dqR_prim_dx_n(1)%vf(i)%sf(k, l, j))

                                    dqR_prim_dx_n(3)%vf(i)%sf(k, l, j) = 25.e-2_wp* &
                                                                         dqR_prim_dx_n(3)%vf(i)%sf(k, l, j)

                                end do
                            end do
                        end do
                    end do

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                                   dq_prim_dx_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i), &
                                                   dq_prim_dz_qp(1)%vf(i))
                    end do

                else

                    do i = iv%beg, iv%end
                        call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                                   dq_prim_dx_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i), &
                                                   dq_prim_dy_qp(1)%vf(i))
                    end do

                end if

            else

                do i = iv%beg, iv%end
                    call s_compute_fd_gradient(q_prim_qp%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i), &
                                               dq_prim_dx_qp(1)%vf(i))
                end do

            end if

        end if

    end subroutine s_get_viscous

    subroutine s_reconstruct_cell_boundary_values_visc(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                       norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf

        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:), intent(inout) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z
        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l

        ! Reconstruction in s1-direction

        if (norm_dir == 1) then
            is1_viscous = ix; is2_viscous = iy; is3_viscous = iz
            weno_dir = 1; is1_viscous%beg = is1_viscous%beg + weno_polyn
            is1_viscous%end = is1_viscous%end - weno_polyn

        elseif (norm_dir == 2) then
            is1_viscous = iy; is2_viscous = ix; is3_viscous = iz
            weno_dir = 2; is1_viscous%beg = is1_viscous%beg + weno_polyn
            is1_viscous%end = is1_viscous%end - weno_polyn

        else
            is1_viscous = iz; is2_viscous = iy; is3_viscous = ix
            weno_dir = 3; is1_viscous%beg = is1_viscous%beg + weno_polyn
            is1_viscous%end = is1_viscous%end - weno_polyn

        end if

        !$acc update device(is1_viscous, is2_viscous, is3_viscous, iv)

        if (n > 0) then
            if (p > 0) then
                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                            weno_dir, &
                            is1_viscous, is2_viscous, is3_viscous)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                            weno_dir, &
                            is1_viscous, is2_viscous, is3_viscous)
            end if
        else
            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                        weno_dir, &
                        is1_viscous, is2_viscous, is3_viscous)
        end if

        if (viscous) then
            if (weno_Re_flux) then
                if (norm_dir == 2) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_viscous%beg, is3_viscous%end
                            do j = is1_viscous%beg, is1_viscous%end
                                do k = is2_viscous%beg, is2_viscous%end
                                    vL_prim_vf(i)%sf(k, j, l) = vL_y(j, k, l, i)
                                    vR_prim_vf(i)%sf(k, j, l) = vR_y(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do j = is1_viscous%beg, is1_viscous%end
                            do k = is2_viscous%beg, is2_viscous%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    vL_prim_vf(i)%sf(l, k, j) = vL_z(j, k, l, i)
                                    vR_prim_vf(i)%sf(l, k, j) = vR_z(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_viscous%beg, is3_viscous%end
                            do k = is2_viscous%beg, is2_viscous%end
                                do j = is1_viscous%beg, is1_viscous%end
                                    vL_prim_vf(i)%sf(j, k, l) = vL_x(j, k, l, i)
                                    vR_prim_vf(i)%sf(j, k, l) = vR_x(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_reconstruct_cell_boundary_values_visc

    subroutine s_reconstruct_cell_boundary_values_visc_deriv(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, &
                                                             norm_dir, vL_prim_vf, vR_prim_vf, ix, iy, iz)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        real(wp), dimension(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, iv%beg:), intent(inout) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z
        type(scalar_field), dimension(iv%beg:iv%end), intent(inout) :: vL_prim_vf, vR_prim_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer, intent(IN) :: norm_dir

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l
        ! Reconstruction in s1-direction

        if (norm_dir == 1) then
            is1_viscous = ix; is2_viscous = iy; is3_viscous = iz
            weno_dir = 1; is1_viscous%beg = is1_viscous%beg + weno_polyn
            is1_viscous%end = is1_viscous%end - weno_polyn

        elseif (norm_dir == 2) then
            is1_viscous = iy; is2_viscous = ix; is3_viscous = iz
            weno_dir = 2; is1_viscous%beg = is1_viscous%beg + weno_polyn
            is1_viscous%end = is1_viscous%end - weno_polyn

        else
            is1_viscous = iz; is2_viscous = iy; is3_viscous = ix
            weno_dir = 3; is1_viscous%beg = is1_viscous%beg + weno_polyn
            is1_viscous%end = is1_viscous%end - weno_polyn

        end if

        !$acc update device(is1_viscous, is2_viscous, is3_viscous, iv)

        if (n > 0) then
            if (p > 0) then

                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                            weno_dir, &
                            is1_viscous, is2_viscous, is3_viscous)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                            weno_dir, &
                            is1_viscous, is2_viscous, is3_viscous)
            end if
        else

            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                        weno_dir, &
                        is1_viscous, is2_viscous, is3_viscous)
        end if

        if (viscous) then
            if (weno_Re_flux) then
                if (norm_dir == 2) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_viscous%beg, is3_viscous%end
                            do j = is1_viscous%beg, is1_viscous%end
                                do k = is2_viscous%beg, is2_viscous%end
                                    vL_prim_vf(i)%sf(k, j, l) = vL_y(j, k, l, i)
                                    vR_prim_vf(i)%sf(k, j, l) = vR_y(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do j = is1_viscous%beg, is1_viscous%end
                            do k = is2_viscous%beg, is2_viscous%end
                                do l = is3_viscous%beg, is3_viscous%end
                                    vL_prim_vf(i)%sf(l, k, j) = vL_z(j, k, l, i)
                                    vR_prim_vf(i)%sf(l, k, j) = vR_z(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                elseif (norm_dir == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do i = iv%beg, iv%end
                        do l = is3_viscous%beg, is3_viscous%end
                            do k = is2_viscous%beg, is2_viscous%end
                                do j = is1_viscous%beg, is1_viscous%end
                                    vL_prim_vf(i)%sf(j, k, l) = vL_x(j, k, l, i)
                                    vR_prim_vf(i)%sf(j, k, l) = vR_x(j, k, l, i)
                                end do
                            end do
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_reconstruct_cell_boundary_values_visc_deriv

    !>  The purpose of this subroutine is to employ the inputted
        !!      left and right cell-boundary integral-averaged variables
        !!      to compute the relevant cell-average first-order spatial
        !!      derivatives in the x-, y- or z-direction by means of the
        !!      scalar divergence theorem.
        !!  @param vL_vf Left cell-boundary integral averages
        !!  @param vR_vf Right cell-boundary integral averages
        !!  @param dv_ds_vf Cell-average first-order spatial derivatives
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_apply_scalar_divergence_theorem(vL_vf, vR_vf, &
                                                 dv_ds_vf, &
                                                 norm_dir, &
                                                 ix, iy, iz, iv_in, &
                                                 dL, dim, buff_size_in)

        ! arrays of cell widths
        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(in) :: vL_vf, vR_vf

        type(scalar_field), &
            dimension(iv%beg:iv%end), &
            intent(inout) :: dv_ds_vf

        integer, intent(in) :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz, iv_in
        integer, intent(in) :: dim, buff_size_in
        real(wp), dimension(-buff_size_in:dim + buff_size_in), intent(in) :: dL

        integer :: i, j, k, l !< Generic loop iterators

        is1_viscous = ix
        is2_viscous = iy
        is3_viscous = iz
        iv = iv_in

        !$acc update device(is1_viscous, is2_viscous, is3_viscous, iv)

        ! First-Order Spatial Derivatives in x-direction
        if (norm_dir == 1) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg + 1, is1_viscous%end - 1
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1._wp/((1._wp + wa_flg)*dL(j)) &
                                *(wa_flg*vL_vf(i)%sf(j + 1, k, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j - 1, k, l))
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in x-direction

            ! First-Order Spatial Derivatives in y-direction
        elseif (norm_dir == 2) then

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg + 1, is2_viscous%end - 1
                    do j = is1_viscous%beg, is1_viscous%end
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1._wp/((1._wp + wa_flg)*dL(k)) &
                                *(wa_flg*vL_vf(i)%sf(j, k + 1, l) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k - 1, l))
                        end do
                    end do
                end do
            end do

            ! END: First-Order Spatial Derivatives in y-direction

            ! First-Order Spatial Derivatives in z-direction
        else

            ! A general application of the scalar divergence theorem that
            ! utilizes the left and right cell-boundary integral-averages,
            ! inside each cell, or an arithmetic mean of these two at the
            ! cell-boundaries, to calculate the cell-averaged first-order
            ! spatial derivatives inside the cell.

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_viscous%beg + 1, is3_viscous%end - 1
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        !$acc loop seq
                        do i = iv%beg, iv%end
                            dv_ds_vf(i)%sf(j, k, l) = &
                                1._wp/((1._wp + wa_flg)*dL(l)) &
                                *(wa_flg*vL_vf(i)%sf(j, k, l + 1) &
                                  + vR_vf(i)%sf(j, k, l) &
                                  - vL_vf(i)%sf(j, k, l) &
                                  - wa_flg*vR_vf(i)%sf(j, k, l - 1))
                        end do
                    end do
                end do
            end do

        end if
        ! END: First-Order Spatial Derivatives in z-direction

    end subroutine s_apply_scalar_divergence_theorem

    !>  Computes the scalar gradient fields via finite differences
        !!  @param var Variable to compute derivative of
        !!  @param grad_x First coordinate direction component of the derivative
        !!  @param grad_y Second coordinate direction component of the derivative
        !!  @param grad_z Third coordinate direction component of the derivative
        !!  @param norm Norm of the gradient vector
        !!  Note: For non-Newtonian fluids, uses 4th order central differences in interior
        !!        For Newtonian fluids, uses 2nd order central differences (original author's method)
    subroutine s_compute_fd_gradient(var, grad_x, grad_y, grad_z)

        type(scalar_field), intent(in) :: var
        type(scalar_field), intent(inout) :: grad_x
        type(scalar_field), intent(inout) :: grad_y
        type(scalar_field), intent(inout) :: grad_z
        type(int_bounds_info) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators
        real(wp) :: fd_coeff_m2, fd_coeff_m1, fd_coeff_p1, fd_coeff_p2 !< 4th order FD coefficients

        ix%beg = 1 - buff_size; ix%end = m + buff_size - 1
        if (n > 0) then
            iy%beg = 1 - buff_size; iy%end = n + buff_size - 1
        else
            iy%beg = 0; iy%end = 0
        end if

        if (p > 0) then
            iz%beg = 1 - buff_size; iz%end = p + buff_size - 1
        else
            iz%beg = 0; iz%end = 0
        end if

        is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

        !$acc update device(is1_viscous, is2_viscous, is3_viscous)

        if (any_non_newtonian .and. buff_size >= 2) then
            ! Non-Newtonian fluids: use 4th order central differences for higher accuracy
            ! The 4th order stencil accesses j±2, so we need to restrict loop bounds
            ! to ensure we don't access outside the allocated array bounds.
            ! Array is allocated from -buff_size to m+buff_size.
            ! For 4th order, loop from max(ix%beg, 2-buff_size) to min(ix%end, m+buff_size-2)
            !$acc parallel loop collapse(3) gang vector default(present) private(fd_coeff_m2, fd_coeff_m1, fd_coeff_p1, fd_coeff_p2)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = max(is1_viscous%beg, 2 - buff_size), min(is1_viscous%end, m + buff_size - 2)
                        ! 4th order central difference: f' = (f_{-2} - 8*f_{-1} + 8*f_{+1} - f_{+2}) / (12*h)
                        ! For non-uniform grids, we compute coefficients inline
                        fd_coeff_m2 = 1._wp/(x_cc(j - 2) - 8._wp*x_cc(j - 1) - x_cc(j + 2) + 8._wp*x_cc(j + 1))
                        fd_coeff_m1 = -8._wp*fd_coeff_m2
                        fd_coeff_p1 = -fd_coeff_m1
                        fd_coeff_p2 = -fd_coeff_m2
                        grad_x%sf(j, k, l) = fd_coeff_m2*var%sf(j - 2, k, l) + &
                                             fd_coeff_m1*var%sf(j - 1, k, l) + &
                                             fd_coeff_p1*var%sf(j + 1, k, l) + &
                                             fd_coeff_p2*var%sf(j + 2, k, l)
                    end do
                end do
            end do

            if (n > 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(fd_coeff_m2, fd_coeff_m1, fd_coeff_p1, fd_coeff_p2)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = max(is2_viscous%beg, 2 - buff_size), min(is2_viscous%end, n + buff_size - 2)
                        do j = is1_viscous%beg, is1_viscous%end
                            fd_coeff_m2 = 1._wp/(y_cc(k - 2) - 8._wp*y_cc(k - 1) - y_cc(k + 2) + 8._wp*y_cc(k + 1))
                            fd_coeff_m1 = -8._wp*fd_coeff_m2
                            fd_coeff_p1 = -fd_coeff_m1
                            fd_coeff_p2 = -fd_coeff_m2
                            grad_y%sf(j, k, l) = fd_coeff_m2*var%sf(j, k - 2, l) + &
                                                 fd_coeff_m1*var%sf(j, k - 1, l) + &
                                                 fd_coeff_p1*var%sf(j, k + 1, l) + &
                                                 fd_coeff_p2*var%sf(j, k + 2, l)
                        end do
                    end do
                end do
            end if

            if (p > 0) then
                !$acc parallel loop collapse(3) gang vector default(present) private(fd_coeff_m2, fd_coeff_m1, fd_coeff_p1, fd_coeff_p2)
                do l = max(is3_viscous%beg, 2 - buff_size), min(is3_viscous%end, p + buff_size - 2)
                    do k = is2_viscous%beg, is2_viscous%end
                        do j = is1_viscous%beg, is1_viscous%end
                            fd_coeff_m2 = 1._wp/(z_cc(l - 2) - 8._wp*z_cc(l - 1) - z_cc(l + 2) + 8._wp*z_cc(l + 1))
                            fd_coeff_m1 = -8._wp*fd_coeff_m2
                            fd_coeff_p1 = -fd_coeff_m1
                            fd_coeff_p2 = -fd_coeff_m2
                            grad_z%sf(j, k, l) = fd_coeff_m2*var%sf(j, k, l - 2) + &
                                                 fd_coeff_m1*var%sf(j, k, l - 1) + &
                                                 fd_coeff_p1*var%sf(j, k, l + 1) + &
                                                 fd_coeff_p2*var%sf(j, k, l + 2)
                        end do
                    end do
                end do
            end if
        else
            ! Newtonian fluids or insufficient buffer: use 2nd order central differences (original author's method)
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3_viscous%beg, is3_viscous%end
                do k = is2_viscous%beg, is2_viscous%end
                    do j = is1_viscous%beg, is1_viscous%end
                        grad_x%sf(j, k, l) = &
                            (var%sf(j + 1, k, l) - var%sf(j - 1, k, l))/ &
                            (x_cc(j + 1) - x_cc(j - 1))
                    end do
                end do
            end do

            if (n > 0) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg, is2_viscous%end
                        do j = is1_viscous%beg, is1_viscous%end
                            grad_y%sf(j, k, l) = &
                                (var%sf(j, k + 1, l) - var%sf(j, k - 1, l))/ &
                                (y_cc(k + 1) - y_cc(k - 1))
                        end do
                    end do
                end do
            end if

            if (p > 0) then
                !$acc parallel loop collapse(3) gang vector default(present)
                do l = is3_viscous%beg, is3_viscous%end
                    do k = is2_viscous%beg, is2_viscous%end
                        do j = is1_viscous%beg, is1_viscous%end
                            grad_z%sf(j, k, l) = &
                                (var%sf(j, k, l + 1) - var%sf(j, k, l - 1))/ &
                                (z_cc(l + 1) - z_cc(l - 1))
                        end do
                    end do
                end do
            end if
        end if

        !$acc parallel loop collapse(2) gang vector default(present)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                grad_x%sf(idwbuff(1)%beg, k, l) = &
                    (-3._wp*var%sf(idwbuff(1)%beg, k, l) + 4._wp*var%sf(idwbuff(1)%beg + 1, k, l) - var%sf(idwbuff(1)%beg + 2, k, l))/ &
                    (x_cc(idwbuff(1)%beg + 2) - x_cc(idwbuff(1)%beg))
                grad_x%sf(idwbuff(1)%end, k, l) = &
                    (+3._wp*var%sf(idwbuff(1)%end, k, l) - 4._wp*var%sf(idwbuff(1)%end - 1, k, l) + var%sf(idwbuff(1)%end - 2, k, l))/ &
                    (x_cc(idwbuff(1)%end) - x_cc(idwbuff(1)%end - 2))
            end do
        end do
        if (n > 0) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    grad_y%sf(j, idwbuff(2)%beg, l) = &
                        (-3._wp*var%sf(j, idwbuff(2)%beg, l) + 4._wp*var%sf(j, idwbuff(2)%beg + 1, l) - var%sf(j, idwbuff(2)%beg + 2, l))/ &
                        (y_cc(idwbuff(2)%beg + 2) - y_cc(idwbuff(2)%beg))
                    grad_y%sf(j, idwbuff(2)%end, l) = &
                        (+3._wp*var%sf(j, idwbuff(2)%end, l) - 4._wp*var%sf(j, idwbuff(2)%end - 1, l) + var%sf(j, idwbuff(2)%end - 2, l))/ &
                        (y_cc(idwbuff(2)%end) - y_cc(idwbuff(2)%end - 2))
                end do
            end do
            if (p > 0) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_z%sf(j, k, idwbuff(3)%beg) = &
                            (-3._wp*var%sf(j, k, idwbuff(3)%beg) + 4._wp*var%sf(j, k, idwbuff(3)%beg + 1) - var%sf(j, k, idwbuff(3)%beg + 2))/ &
                            (z_cc(idwbuff(3)%beg + 2) - z_cc(is3_viscous%beg))
                        grad_z%sf(j, k, idwbuff(3)%end) = &
                            (+3._wp*var%sf(j, k, idwbuff(3)%end) - 4._wp*var%sf(j, k, idwbuff(3)%end - 1) + var%sf(j, k, idwbuff(3)%end - 2))/ &
                            (z_cc(idwbuff(3)%end) - z_cc(idwbuff(3)%end - 2))
                    end do
                end do
            end if
        end if

        if (bc_x%beg <= BC_GHOST_EXTRAP) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(0, k, l) = (-3._wp*var%sf(0, k, l) + 4._wp*var%sf(1, k, l) - var%sf(2, k, l))/ &
                                         (x_cc(2) - x_cc(0))
                end do
            end do
        end if
        if (bc_x%end <= BC_GHOST_EXTRAP) then
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    grad_x%sf(m, k, l) = (3._wp*var%sf(m, k, l) - 4._wp*var%sf(m - 1, k, l) + var%sf(m - 2, k, l))/ &
                                         (x_cc(m) - x_cc(m - 2))
                end do
            end do
        end if
        if (n > 0) then
            if (bc_y%beg <= BC_GHOST_EXTRAP .and. bc_y%beg /= BC_NULL) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, 0, l) = (-3._wp*var%sf(j, 0, l) + 4._wp*var%sf(j, 1, l) - var%sf(j, 2, l))/ &
                                             (y_cc(2) - y_cc(0))
                    end do
                end do
            end if
            if (bc_y%end <= BC_GHOST_EXTRAP) then
                !$acc parallel loop collapse(2) gang vector default(present)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        grad_y%sf(j, n, l) = (3._wp*var%sf(j, n, l) - 4._wp*var%sf(j, n - 1, l) + var%sf(j, n - 2, l))/ &
                                             (y_cc(n) - y_cc(n - 2))
                    end do
                end do
            end if
            if (p > 0) then
                if (bc_z%beg <= BC_GHOST_EXTRAP) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, 0) = &
                                (-3._wp*var%sf(j, k, 0) + 4._wp*var%sf(j, k, 1) - var%sf(j, k, 2))/ &
                                (z_cc(2) - z_cc(0))
                        end do
                    end do
                end if
                if (bc_z%end <= BC_GHOST_EXTRAP) then
                    !$acc parallel loop collapse(2) gang vector default(present)
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            grad_z%sf(j, k, p) = &
                                (3._wp*var%sf(j, k, p) - 4._wp*var%sf(j, k, p - 1) + var%sf(j, k, p - 2))/ &
                                (z_cc(p) - z_cc(p - 2))
                        end do
                    end do
                end if
            end if
        end if

    end subroutine s_compute_fd_gradient

    !>  Computes the viscosity gradient correction term for non-Newtonian fluids
    !!  This term arises from the variable viscosity in the Navier-Stokes equations:
    !!  ∇·(μS) = μ∇²u + ∇μ·S (where S is the strain rate tensor)
    !!  The standard MFC flux computation handles μ∇²u via interface viscosity,
    !!  but the ∇μ·S term is missing for variable μ. This subroutine computes
    !!  that correction term.
    !!
    !!  PHASE-BASED APPROACH: For multi-phase flows, we compute per-phase viscosity
    !!  gradients and weight them by volume fraction:
    !!  correction = Σ_q α_q * ∇μ_q · S
    !!  This is more accurate than computing the gradient of mixture viscosity.
    !!
    !!  For momentum equation i: correction = Σ_q α_q * (dμ_q/dx_j * S_ij)
    !!  For energy equation: correction = Σ_q α_q * (dμ_q/dx_j * S_ij * u_i)
    !!
    !!  @param q_prim_vf Cell-averaged primitive variables
    !!  @param dq_prim_dx_vf Cell-averaged x-derivatives of velocity
    !!  @param dq_prim_dy_vf Cell-averaged y-derivatives of velocity
    !!  @param dq_prim_dz_vf Cell-averaged z-derivatives of velocity
    !!  @param rhs_vf Cell-averaged RHS variables (to be corrected)
    !!  @param ix, iy, iz Index bounds
    subroutine s_compute_viscous_gradient_correction(q_prim_vf, &
                                                      dq_prim_dx_vf, dq_prim_dy_vf, dq_prim_dz_vf, &
                                                      rhs_vf, ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: dq_prim_dx_vf, dq_prim_dy_vf, dq_prim_dz_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! Local variables - per-phase viscosities
        real(wp), dimension(num_fluids) :: mu_shear_phases, mu_bulk_phases
        real(wp), dimension(num_fluids) :: mu_shear_phases_xp, mu_shear_phases_xm
        real(wp), dimension(num_fluids) :: mu_shear_phases_yp, mu_shear_phases_ym
        real(wp), dimension(num_fluids) :: mu_shear_phases_zp, mu_shear_phases_zm
        real(wp), dimension(num_fluids) :: mu_shear_phases_xp2, mu_shear_phases_xm2
        real(wp), dimension(num_fluids) :: mu_shear_phases_yp2, mu_shear_phases_ym2
        real(wp), dimension(num_fluids) :: mu_shear_phases_zp2, mu_shear_phases_zm2
        real(wp), dimension(num_fluids) :: dmu_dx_phases, dmu_dy_phases, dmu_dz_phases
        real(wp), dimension(num_fluids) :: alpha_visc
        real(wp) :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz
        real(wp) :: div_u
        real(wp) :: S_xx, S_yy, S_zz, S_xy, S_xz, S_yz
        real(wp) :: corr_mom_x, corr_mom_y, corr_mom_z, corr_energy
        real(wp) :: phase_corr_x, phase_corr_y, phase_corr_z
        real(wp) :: u_vel, v_vel, w_vel
        real(wp) :: coeff_4th
        logical :: use_4th_order_x, use_4th_order_y, use_4th_order_z
        integer :: j, k, l, i, q

        ! Only proceed if any fluid is non-Newtonian
        if (.not. any_non_newtonian) return

        is1_viscous = ix; is2_viscous = iy; is3_viscous = iz

        !$acc update device(is1_viscous, is2_viscous, is3_viscous)

        !$acc parallel loop collapse(3) gang vector default(present) &
        !$acc private(mu_shear_phases, mu_bulk_phases, alpha_visc) &
        !$acc private(mu_shear_phases_xp, mu_shear_phases_xm, mu_shear_phases_yp, mu_shear_phases_ym, mu_shear_phases_zp, mu_shear_phases_zm) &
        !$acc private(mu_shear_phases_xp2, mu_shear_phases_xm2, mu_shear_phases_yp2, mu_shear_phases_ym2, mu_shear_phases_zp2, mu_shear_phases_zm2) &
        !$acc private(dmu_dx_phases, dmu_dy_phases, dmu_dz_phases) &
        !$acc private(du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz) &
        !$acc private(div_u, S_xx, S_yy, S_zz, S_xy, S_xz, S_yz) &
        !$acc private(corr_mom_x, corr_mom_y, corr_mom_z, corr_energy, phase_corr_x, phase_corr_y, phase_corr_z) &
        !$acc private(u_vel, v_vel, w_vel, coeff_4th, use_4th_order_x, use_4th_order_y, use_4th_order_z)
        do l = is3_viscous%beg, is3_viscous%end
            do k = is2_viscous%beg, is2_viscous%end
                do j = is1_viscous%beg, is1_viscous%end

                    ! Get volume fractions at this cell
                    !$acc loop seq
                    do i = 1, num_fluids
                        if (bubbles_euler .and. num_fluids == 1) then
                            alpha_visc(i) = 1._wp - q_prim_vf(E_idx + i)%sf(j, k, l)
                        else
                            alpha_visc(i) = q_prim_vf(E_idx + i)%sf(j, k, l)
                        end if
                    end do

                    ! Determine if 4th order can be used for viscosity gradient
                    ! Since s_compute_phase_viscosity_at_cell also uses 4th-order FD internally
                    ! (which accesses j±2), and we call it at j±2, the total stencil width is 4 cells.
                    ! We need buff_size >= 4 and sufficient distance from boundaries.
                    use_4th_order_x = (buff_size >= 4) .and. (j > 3) .and. (j < m - 3)
                    use_4th_order_y = (buff_size >= 4) .and. (n > 0) .and. (k > 3) .and. (k < n - 3)
                    use_4th_order_z = (buff_size >= 4) .and. (p > 0) .and. (l > 3) .and. (l < p - 3)

                    ! Initialize per-phase gradients to zero
                    !$acc loop seq
                    do q = 1, num_fluids
                        dmu_dx_phases(q) = 0._wp
                        dmu_dy_phases(q) = 0._wp
                        dmu_dz_phases(q) = 0._wp
                    end do

                    ! ===== X-DIRECTION GRADIENT =====
                    ! Note: We don't pass pre-computed gradients here because the indices
                    ! (j±1, j±2, etc.) may be outside the local gradient array bounds.
                    ! Instead, let s_compute_phase_viscosity_at_cell compute gradients
                    ! directly from q_prim_vf which includes ghost cells.
                    if (use_4th_order_x) then
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j + 1, k, l, mu_shear_phases_xp, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j - 1, k, l, mu_shear_phases_xm, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j + 2, k, l, mu_shear_phases_xp2, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j - 2, k, l, mu_shear_phases_xm2, mu_bulk_phases)
                        coeff_4th = 1._wp/(x_cc(j - 2) - 8._wp*x_cc(j - 1) - x_cc(j + 2) + 8._wp*x_cc(j + 1))
                        !$acc loop seq
                        do q = 1, num_fluids
                            dmu_dx_phases(q) = coeff_4th*(mu_shear_phases_xm2(q) - 8._wp*mu_shear_phases_xm(q) + &
                                                          8._wp*mu_shear_phases_xp(q) - mu_shear_phases_xp2(q))
                        end do
                    else if (j > 0 .and. j < m) then
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j + 1, k, l, mu_shear_phases_xp, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j - 1, k, l, mu_shear_phases_xm, mu_bulk_phases)
                        !$acc loop seq
                        do q = 1, num_fluids
                            dmu_dx_phases(q) = (mu_shear_phases_xp(q) - mu_shear_phases_xm(q))/(x_cc(j + 1) - x_cc(j - 1))
                        end do
                    end if

                    ! ===== Y-DIRECTION GRADIENT =====
                    if (use_4th_order_y) then
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k + 1, l, mu_shear_phases_yp, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k - 1, l, mu_shear_phases_ym, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k + 2, l, mu_shear_phases_yp2, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k - 2, l, mu_shear_phases_ym2, mu_bulk_phases)
                        coeff_4th = 1._wp/(y_cc(k - 2) - 8._wp*y_cc(k - 1) - y_cc(k + 2) + 8._wp*y_cc(k + 1))
                        !$acc loop seq
                        do q = 1, num_fluids
                            dmu_dy_phases(q) = coeff_4th*(mu_shear_phases_ym2(q) - 8._wp*mu_shear_phases_ym(q) + &
                                                          8._wp*mu_shear_phases_yp(q) - mu_shear_phases_yp2(q))
                        end do
                    else if (n > 0 .and. k > 0 .and. k < n) then
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k + 1, l, mu_shear_phases_yp, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k - 1, l, mu_shear_phases_ym, mu_bulk_phases)
                        !$acc loop seq
                        do q = 1, num_fluids
                            dmu_dy_phases(q) = (mu_shear_phases_yp(q) - mu_shear_phases_ym(q))/(y_cc(k + 1) - y_cc(k - 1))
                        end do
                    end if

                    ! ===== Z-DIRECTION GRADIENT =====
                    if (use_4th_order_z) then
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k, l + 1, mu_shear_phases_zp, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k, l - 1, mu_shear_phases_zm, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k, l + 2, mu_shear_phases_zp2, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k, l - 2, mu_shear_phases_zm2, mu_bulk_phases)
                        coeff_4th = 1._wp/(z_cc(l - 2) - 8._wp*z_cc(l - 1) - z_cc(l + 2) + 8._wp*z_cc(l + 1))
                        !$acc loop seq
                        do q = 1, num_fluids
                            dmu_dz_phases(q) = coeff_4th*(mu_shear_phases_zm2(q) - 8._wp*mu_shear_phases_zm(q) + &
                                                          8._wp*mu_shear_phases_zp(q) - mu_shear_phases_zp2(q))
                        end do
                    else if (p > 0 .and. l > 0 .and. l < p) then
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k, l + 1, mu_shear_phases_zp, mu_bulk_phases)
                        call s_compute_phase_viscosity_at_cell(q_prim_vf, j, k, l - 1, mu_shear_phases_zm, mu_bulk_phases)
                        !$acc loop seq
                        do q = 1, num_fluids
                            dmu_dz_phases(q) = (mu_shear_phases_zp(q) - mu_shear_phases_zm(q))/(z_cc(l + 1) - z_cc(l - 1))
                        end do
                    end if

                    ! ===== GET VELOCITY GRADIENTS =====
                    du_dx = dq_prim_dx_vf(momxb)%sf(j, k, l)
                    if (n > 0) then
                        du_dy = dq_prim_dy_vf(momxb)%sf(j, k, l)
                        dv_dx = dq_prim_dx_vf(momxb + 1)%sf(j, k, l)
                        dv_dy = dq_prim_dy_vf(momxb + 1)%sf(j, k, l)
                    else
                        du_dy = 0._wp
                        dv_dx = 0._wp
                        dv_dy = 0._wp
                    end if
                    if (p > 0) then
                        du_dz = dq_prim_dz_vf(momxb)%sf(j, k, l)
                        dv_dz = dq_prim_dz_vf(momxb + 1)%sf(j, k, l)
                        dw_dx = dq_prim_dx_vf(momxb + 2)%sf(j, k, l)
                        dw_dy = dq_prim_dy_vf(momxb + 2)%sf(j, k, l)
                        dw_dz = dq_prim_dz_vf(momxb + 2)%sf(j, k, l)
                    else
                        du_dz = 0._wp
                        dv_dz = 0._wp
                        dw_dx = 0._wp
                        dw_dy = 0._wp
                        dw_dz = 0._wp
                    end if

                    ! ===== COMPUTE STRAIN RATE TENSOR =====
                    div_u = du_dx + dv_dy + dw_dz

                    ! Deviatoric strain rate: S_ij = (du_i/dx_j + du_j/dx_i)/2 - (1/3)*div(u)*delta_ij
                    S_xx = du_dx - (1._wp/3._wp)*div_u
                    S_yy = dv_dy - (1._wp/3._wp)*div_u
                    S_zz = dw_dz - (1._wp/3._wp)*div_u
                    S_xy = 0.5_wp*(du_dy + dv_dx)
                    S_xz = 0.5_wp*(du_dz + dw_dx)
                    S_yz = 0.5_wp*(dv_dz + dw_dy)

                    ! ===== COMPUTE PHASE-WEIGHTED CORRECTION =====
                    ! correction = Σ_q α_q * 2 * (∇μ_q · S)
                    corr_mom_x = 0._wp
                    corr_mom_y = 0._wp
                    corr_mom_z = 0._wp

                    !$acc loop seq
                    do q = 1, num_fluids
                        ! Per-phase correction: 2 * (dμ_q/dx_j * S_ij)
                        phase_corr_x = 2._wp*(dmu_dx_phases(q)*S_xx + dmu_dy_phases(q)*S_xy + dmu_dz_phases(q)*S_xz)
                        phase_corr_y = 2._wp*(dmu_dx_phases(q)*S_xy + dmu_dy_phases(q)*S_yy + dmu_dz_phases(q)*S_yz)
                        phase_corr_z = 2._wp*(dmu_dx_phases(q)*S_xz + dmu_dy_phases(q)*S_yz + dmu_dz_phases(q)*S_zz)

                        ! Weight by volume fraction and accumulate
                        corr_mom_x = corr_mom_x + alpha_visc(q)*phase_corr_x
                        corr_mom_y = corr_mom_y + alpha_visc(q)*phase_corr_y
                        corr_mom_z = corr_mom_z + alpha_visc(q)*phase_corr_z
                    end do

                    ! ===== ENERGY CORRECTION =====
                    u_vel = q_prim_vf(momxb)%sf(j, k, l)
                    if (n > 0) then
                        v_vel = q_prim_vf(momxb + 1)%sf(j, k, l)
                    else
                        v_vel = 0._wp
                    end if
                    if (p > 0) then
                        w_vel = q_prim_vf(momxb + 2)%sf(j, k, l)
                    else
                        w_vel = 0._wp
                    end if

                    corr_energy = u_vel*corr_mom_x + v_vel*corr_mom_y + w_vel*corr_mom_z

                    ! ===== ADD CORRECTIONS TO RHS =====
                    rhs_vf(momxb)%sf(j, k, l) = rhs_vf(momxb)%sf(j, k, l) + corr_mom_x
                    if (n > 0) then
                        rhs_vf(momxb + 1)%sf(j, k, l) = rhs_vf(momxb + 1)%sf(j, k, l) + corr_mom_y
                    end if
                    if (p > 0) then
                        rhs_vf(momxb + 2)%sf(j, k, l) = rhs_vf(momxb + 2)%sf(j, k, l) + corr_mom_z
                    end if
                    rhs_vf(E_idx)%sf(j, k, l) = rhs_vf(E_idx)%sf(j, k, l) + corr_energy

                end do
            end do
        end do

    end subroutine s_compute_viscous_gradient_correction

    impure subroutine s_finalize_viscous_module()

        ! Res_viscous was removed - Re_visc is now computed dynamically via m_re_visc

    end subroutine s_finalize_viscous_module

end module m_viscous
