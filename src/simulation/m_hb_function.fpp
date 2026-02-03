!>
!! @file m_hb_function.f90
!! @brief Contains module m_hb_function for Herschel-Bulkley non-Newtonian viscosity
!!        with Papanastasiou regularization for improved low shear rate behavior
#:include 'macros.fpp'

!> @brief The module contains functions to compute non-Newtonian viscosity
!!        using the Herschel-Bulkley model with Papanastasiou regularization
!!
!!        Standard Herschel-Bulkley: μ = τ₀/γ̇ + K·γ̇^(n-1)
!!        This has a singularity as γ̇ → 0
!!
!!        Papanastasiou regularization: μ = (τ₀/γ̇)·[1 - exp(-m·γ̇)] + K·γ̇^(n-1)
!!        As γ̇ → 0: μ → τ₀·m (finite, no singularity)
!!        As γ̇ → ∞: recovers standard HB model
!!
!!        Parameter m controls the transition sharpness:
!!        - Large m (>1000): Sharp transition, closer to ideal Bingham/HB behavior
!!        - Small m (~100): Smoother transition, better numerical stability
module m_hb_function

    use m_derived_types        !< Definitions of the derived types
    use m_global_parameters    !< Definitions of the global parameters

    private; public :: f_compute_hb_viscosity, &
                      f_compute_shear_rate_from_components

contains

    !> Computes Herschel-Bulkley viscosity with Papanastasiou regularization
    !! 
    !! The Papanastasiou model provides a smooth, regularized approximation:
    !!   μ = (τ₀/γ̇)·[1 - exp(-m·γ̇)] + K·γ̇^(n-1)
    !!
    !! This avoids the singularity at γ̇ = 0 present in the standard HB model.
    !! As m → ∞, the model approaches the ideal HB behavior.
    !!
    !! @param tau0 Yield stress
    !! @param K Consistency index
    !! @param n Flow behavior index
    !! @param mu_min Minimum viscosity limit
    !! @param mu_max Maximum viscosity limit
    !! @param shear_rate Shear rate magnitude (γ̇)
    !! @param hb_m Papanastasiou regularization parameter (m)
    !! @return Viscosity
    pure function f_compute_hb_viscosity(tau0, K, n, mu_min, mu_max, shear_rate, hb_m) result(mu)
        !$acc routine seq

        real(wp), intent(in) :: tau0, K, n, mu_min, mu_max, shear_rate, hb_m
        real(wp) :: mu
        real(wp) :: yield_term, power_law_term, exp_term

        exp_term = exp(-hb_m*shear_rate)
        yield_term = tau0*(1._wp - exp_term)/shear_rate
        power_law_term = K*(shear_rate**(n - 1._wp))

        mu = yield_term + power_law_term
        mu = min(max(mu, mu_min), mu_max)

    end function f_compute_hb_viscosity

    !> Computes shear rate from strain rate tensor components
    !! Works for 1D, 2D, and 3D cases (set D_zz, D_xz, D_yz to 0 for 2D/1D)
    !!
    !! The shear rate is computed as: γ̇ = √(2·D_ij·D_ij)
    !! where D_ij is the strain rate tensor
    !!
    !! @param D_xx Normal strain rate in x direction (du/dx)
    !! @param D_yy Normal strain rate in y direction (dv/dy)
    !! @param D_zz Normal strain rate in z direction (dw/dz)
    !! @param D_xy Shear strain rate component xy: 0.5*(du/dy + dv/dx)
    !! @param D_xz Shear strain rate component xz: 0.5*(du/dz + dw/dx)
    !! @param D_yz Shear strain rate component yz: 0.5*(dv/dz + dw/dy)
    !! @return Shear rate magnitude (γ̇)
    pure function f_compute_shear_rate_from_components(D_xx, D_yy, D_zz, D_xy, D_xz, D_yz) result(shear_rate)
        !$acc routine seq

        real(wp), intent(in) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
        real(wp) :: shear_rate

        ! Calculate shear rate from strain rate components
        ! |γ̇| = √(2·D_ij·D_ij) where D_ij is the strain rate tensor
        ! For symmetric tensor: 2·D_ij·D_ij = 2*(D_xx² + D_yy² + D_zz² + 2*(D_xy² + D_xz² + D_yz²))
        shear_rate = sqrt(2._wp*(D_xx*D_xx + D_yy*D_yy + D_zz*D_zz + &
                                 2._wp*(D_xy*D_xy + D_xz*D_xz + D_yz*D_yz)))
        
        ! Small minimum to prevent exactly zero (for numerical safety in downstream calculations)
        ! With Papanastasiou regularization, this can be very small
        shear_rate = max(shear_rate, 1.0e-8_wp)

    end function f_compute_shear_rate_from_components

end module m_hb_function
