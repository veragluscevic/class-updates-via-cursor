#!/usr/bin/env python3
"""
example_dmeff.py — Example script demonstrating the dmeff extension to CLASS.

This script shows how to:
  1. Run CLASS with dmeff parameters from Python
  2. Access dmeff-specific derived quantities
  3. Compare dmeff results against vanilla Lambda-CDM
  4. Verify the accuracy of the dmeff port

Usage (from the class_dmeff_uptodate/ directory):
    python test_dmeff/example_dmeff.py

Prerequisites:
    - Build CLASS:  make clean && make
    - Install the Python wrapper:  cd python && python setup.py install
"""

import sys
import os
import numpy as np

try:
    from classy import Class
except ImportError:
    print("ERROR: Cannot import classy. Install the Python wrapper first:")
    print("  cd python && python setup.py install")
    sys.exit(1)


def run_vanilla():
    """Run standard Lambda-CDM (no dmeff) for comparison."""
    cosmo = Class()
    cosmo.set({
        'h': 0.6732,
        'T_cmb': 2.7255,
        'omega_b': 0.022032,
        'omega_cdm': 0.12,
        'A_s': 2.1e-9,
        'n_s': 0.9660,
        'tau_reio': 0.0543,
        'output': 'tCl,mPk',
        'l_max_scalars': 2500,
        'P_k_max_1/Mpc': 5.0,
    })
    cosmo.compute()
    return cosmo


def run_dmeff_coulomb():
    """Run CLASS with Coulomb-like dmeff interaction (npow = -4)."""
    cosmo = Class()
    cosmo.set({
        'h': 0.6732,
        'T_cmb': 2.7255,
        'omega_b': 0.022032,
        'omega_cdm': 0.0,         # all dark matter is dmeff
        'omega_dmeff': 0.12,
        'A_s': 2.1e-9,
        'n_s': 0.9660,
        'tau_reio': 0.0543,
        'm_dmeff': 0.01,           # mass in GeV
        'N_dmeff': 1,              # one interaction term
        'sigma_dmeff': 1e-41,      # cross section in cm^2
        'npow_dmeff': -4,          # velocity power law index
        'dmeff_target': 'hydrogen',# scattering target
        'Vrel_dmeff': 29.0,        # initial bulk velocity in km/s
        'output': 'tCl,mPk',
        'l_max_scalars': 2500,
        'P_k_max_1/Mpc': 5.0,
    })
    cosmo.compute()
    return cosmo


def run_dmeff_constant():
    """Run CLASS with constant cross section (npow = 0)."""
    cosmo = Class()
    cosmo.set({
        'h': 0.6732,
        'T_cmb': 2.7255,
        'omega_b': 0.022032,
        'omega_cdm': 0.0,
        'omega_dmeff': 0.12,
        'A_s': 2.1e-9,
        'n_s': 0.9660,
        'tau_reio': 0.0543,
        'm_dmeff': 0.01,
        'N_dmeff': 1,
        'sigma_dmeff': 1e-30,
        'npow_dmeff': 0,
        'dmeff_target': 'baryon',
        'Vrel_dmeff': 29.0,
        'output': 'tCl,mPk',
        'l_max_scalars': 2500,
        'P_k_max_1/Mpc': 5.0,
    })
    cosmo.compute()
    return cosmo


def main():
    print("=" * 60)
    print("  dmeff Extension — Example and Verification Script")
    print("=" * 60)

    # ---------------------------------------------------------------
    # 1. Run models
    # ---------------------------------------------------------------
    print("\n[1] Running vanilla Lambda-CDM...")
    vanilla = run_vanilla()
    print("    Done. Omega_cdm = {:.6f}".format(vanilla.Omega0_cdm()))

    print("\n[2] Running dmeff with Coulomb-like interaction (npow=-4)...")
    dmeff_coulomb = run_dmeff_coulomb()
    print("    Done. Omega_dmeff = {:.6f}".format(dmeff_coulomb.Omega0_dmeff()))

    print("\n[3] Running dmeff with constant cross section (npow=0)...")
    dmeff_constant = run_dmeff_constant()
    print("    Done. Omega_dmeff = {:.6f}".format(dmeff_constant.Omega0_dmeff()))

    # ---------------------------------------------------------------
    # 2. Access derived parameters
    # ---------------------------------------------------------------
    print("\n" + "-" * 60)
    print("  Derived Parameters")
    print("-" * 60)

    derived = dmeff_coulomb.get_current_derived_parameters(
        ['z_dmeff_decoupling', 'z_rec', 'tau_rec', '100*theta_s']
    )
    print("  z_dmeff_decoupling = {:.4e}".format(derived['z_dmeff_decoupling']))
    print("  z_rec              = {:.2f}".format(derived['z_rec']))
    print("  100*theta_s        = {:.6f}".format(derived['100*theta_s']))
    print("  Omega0_dmeff       = {:.6f}".format(dmeff_coulomb.Omega0_dmeff()))

    # ---------------------------------------------------------------
    # 3. Compare CMB power spectra
    # ---------------------------------------------------------------
    print("\n" + "-" * 60)
    print("  CMB TT Power Spectrum Comparison")
    print("-" * 60)

    cl_vanilla = vanilla.lensed_cl(2500)
    cl_coulomb = dmeff_coulomb.lensed_cl(2500)
    cl_constant = dmeff_constant.lensed_cl(2500)

    ells = np.arange(2, 2501)
    factor = ells * (ells + 1) / (2 * np.pi) * 1e12  # convert to D_l in muK^2

    print("\n  l(l+1)C_l^TT / 2pi  [muK^2]")
    print("  {:>6s}  {:>12s}  {:>12s}  {:>12s}".format(
        "l", "Vanilla", "Coulomb", "Constant"))
    for ell in [10, 100, 500, 1000, 2000]:
        idx = ell - 2
        d_van = cl_vanilla['tt'][ell] * ell * (ell + 1) / (2 * np.pi) * 1e12
        d_cou = cl_coulomb['tt'][ell] * ell * (ell + 1) / (2 * np.pi) * 1e12
        d_con = cl_constant['tt'][ell] * ell * (ell + 1) / (2 * np.pi) * 1e12
        print("  {:>6d}  {:>12.4f}  {:>12.4f}  {:>12.4f}".format(
            ell, d_van, d_cou, d_con))

    # ---------------------------------------------------------------
    # 4. Compare matter power spectra
    # ---------------------------------------------------------------
    print("\n" + "-" * 60)
    print("  Matter Power Spectrum P(k) Comparison  [z=0]")
    print("-" * 60)

    k_values = [0.01, 0.05, 0.1, 0.5, 1.0]
    print("\n  {:>10s}  {:>12s}  {:>12s}  {:>12s}  {:>12s}".format(
        "k [h/Mpc]", "Vanilla", "Coulomb", "Constant", "Coul/Van"))
    for k in k_values:
        k_per_Mpc = k * vanilla.h()
        pk_van = vanilla.pk(k_per_Mpc, 0.0)
        pk_cou = dmeff_coulomb.pk(k_per_Mpc, 0.0)
        pk_con = dmeff_constant.pk(k_per_Mpc, 0.0)
        ratio = pk_cou / pk_van if pk_van > 0 else 0
        print("  {:>10.3f}  {:>12.4e}  {:>12.4e}  {:>12.4e}  {:>12.4f}".format(
            k, pk_van, pk_cou, pk_con, ratio))

    # ---------------------------------------------------------------
    # 5. Verification checks
    # ---------------------------------------------------------------
    print("\n" + "-" * 60)
    print("  Verification Checks")
    print("-" * 60)

    all_ok = True

    # Check that Omega0_dmeff is set correctly
    omega_dmeff = dmeff_coulomb.Omega0_dmeff() * dmeff_coulomb.h()**2
    expected = 0.12
    rel_err = abs(omega_dmeff - expected) / expected * 100
    ok = rel_err < 0.1
    all_ok = all_ok and ok
    print("  omega_dmeff = {:.6f} (expected {:.6f}): {:.4f}% — {}".format(
        omega_dmeff, expected, rel_err, "PASS" if ok else "FAIL"))

    # Check that vanilla has zero Omega0_dmeff
    ok = vanilla.Omega0_dmeff() == 0.0
    all_ok = all_ok and ok
    print("  Vanilla Omega0_dmeff = {:.6f}: {}".format(
        vanilla.Omega0_dmeff(), "PASS" if ok else "FAIL"))

    # Check z_dmeff_decoupling is reasonable
    z_dec = derived['z_dmeff_decoupling']
    ok = z_dec > 0
    all_ok = all_ok and ok
    print("  z_dmeff_decoupling = {:.4e}: {}".format(
        z_dec, "PASS" if ok else "FAIL"))

    # Check that dmeff suppresses small-scale power
    k_test = 1.0 * vanilla.h()
    pk_van = vanilla.pk(k_test, 0.0)
    pk_cou = dmeff_coulomb.pk(k_test, 0.0)
    suppression = pk_cou / pk_van
    ok = suppression < 1.0  # dmeff should suppress power at small scales
    all_ok = all_ok and ok
    print("  P(k=1)/P_vanilla(k=1) = {:.4f} (expect < 1.0): {}".format(
        suppression, "PASS" if ok else "FAIL"))

    print("\n" + "=" * 60)
    if all_ok:
        print("  ALL VERIFICATION CHECKS PASSED")
    else:
        print("  SOME CHECKS FAILED — see above for details")
    print("=" * 60)

    # Cleanup
    vanilla.struct_cleanup()
    dmeff_coulomb.struct_cleanup()
    dmeff_constant.struct_cleanup()


if __name__ == '__main__':
    main()
