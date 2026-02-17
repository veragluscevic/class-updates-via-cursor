#!/usr/bin/env python3
"""
Compare CLASS output files against reference baselines.

Usage:
    python compare_outputs.py TEST_NAME [--check CHECKS] [--tolerance TOL]
                                        [--compare-to OTHER_TEST]
                                        [--output-dir DIR] [--reference-dir DIR]

Examples:
    python compare_outputs.py test_coulomb --check all
    python compare_outputs.py test_coulomb --check thermodynamics
    python compare_outputs.py test_coulomb --check cl,pk
    python compare_outputs.py test_coulomb_newtonian --compare-to test_coulomb --check cl,pk --tolerance 0.01

The script compares output files in output/ (or --output-dir) against reference
files in test_dmeff/reference/TEST_NAME/ (or --reference-dir).

Checks available:
    background       - Compare rho_dmeff(z), T_dmeff(z) from background file
    thermodynamics   - Compare T_dmeff(z), T_b(z), rates from thermodynamics file
    cl               - Compare C_l^TT at representative multipoles
    pk               - Compare P(k) at representative wavenumbers
    tk               - Compare transfer functions (d_dmeff if present)
    all              - Run all applicable checks
"""

import argparse
import os
import re
import sys
import numpy as np


def load_class_output(filepath):
    """Load a CLASS output file, skipping comment lines.

    CLASS header format uses numbered columns like:
        #  1:z   2:conf. time [Mpc]  3:x_e  ...  11:T_dmeff
    We extract the N:label pairs and use N-1 as 0-indexed column.
    """
    if not os.path.exists(filepath):
        return None, None
    header_line = None
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_line = line
            else:
                break
    if header_line is None:
        return None, None

    # Extract "N:label" patterns from the header line
    # This handles multi-word column names by only taking the first word after N:
    matches = re.findall(r'(\d+):(\S+)', header_line)
    col_names = {}
    for num_str, label in matches:
        idx = int(num_str) - 1  # convert to 0-indexed
        col_names[idx] = label

    data = np.loadtxt(filepath)
    return col_names, data


def find_column(col_names, target):
    """Find column index matching target name (case-insensitive).

    col_names is a dict {col_idx: label} from load_class_output.
    First tries exact match, then substring match.
    """
    target_lower = target.lower()
    for idx, name in col_names.items():
        if name.lower() == target_lower:
            return idx
    for idx, name in col_names.items():
        if target_lower in name.lower():
            return idx
    return None


def relative_error(ref, new, floor=1e-30):
    """Compute relative error, using floor to avoid division by zero."""
    denom = np.maximum(np.abs(ref), floor)
    return np.abs(new - ref) / denom


def interpolate_to_grid(x_ref, y_ref, x_new, y_new, x_points):
    """Interpolate both reference and new data to the same x points."""
    ref_interp = np.interp(x_points, x_ref, y_ref)
    new_interp = np.interp(x_points, x_new, y_new)
    return ref_interp, new_interp


def check_background(test_name, output_dir, ref_dir, tolerance):
    """Compare background output: rho_dmeff, T_dmeff."""
    print("\n--- Background check ---")
    ref_file = os.path.join(ref_dir, f"{test_name}_background.dat")
    new_file = os.path.join(output_dir, f"{test_name}_background.dat")

    ref_names, ref_data = load_class_output(ref_file)
    new_names, new_data = load_class_output(new_file)

    if ref_data is None:
        print(f"  SKIP: Reference file not found: {ref_file}")
        return True
    if new_data is None:
        print(f"  FAIL: Output file not found: {new_file}")
        return False

    passed = True
    z_ref = ref_data[:, 0]
    z_new = new_data[:, 0]

    z_checks = [0, 100, 1000, 10000]

    # Check rho_dmeff
    col_ref = find_column(ref_names, 'rho_dmeff')
    col_new = find_column(new_names, 'rho_dmeff')
    if col_ref is not None and col_new is not None:
        for z_val in z_checks:
            r_ref = np.interp(z_val, z_ref, ref_data[:, col_ref])
            r_new = np.interp(z_val, z_new, new_data[:, col_new])
            rel_err = abs(r_new - r_ref) / max(abs(r_ref), 1e-30) * 100
            status = "PASS" if rel_err < tolerance else "FAIL"
            if status == "FAIL":
                passed = False
            print(f"  rho_dmeff(z={z_val:>6}): ref={r_ref:.6e}, new={r_new:.6e}, "
                  f"rel_err={rel_err:.4f}% [{status}]")
    else:
        if col_ref is not None:
            print(f"  FAIL: rho_dmeff column not found in new output")
            passed = False
        else:
            print(f"  SKIP: rho_dmeff column not in reference (vanilla test?)")

    # Check T_dmeff
    col_ref = find_column(ref_names, 'T_dmeff')
    col_new = find_column(new_names, 'T_dmeff')
    if col_ref is not None and col_new is not None:
        for z_val in [10, 100, 1000]:
            t_ref = np.interp(z_val, z_ref, ref_data[:, col_ref])
            t_new = np.interp(z_val, z_new, new_data[:, col_new])
            rel_err = abs(t_new - t_ref) / max(abs(t_ref), 1e-30) * 100
            status = "PASS" if rel_err < tolerance else "FAIL"
            if status == "FAIL":
                passed = False
            print(f"  T_dmeff(z={z_val:>6}):  ref={t_ref:.6e}, new={t_new:.6e}, "
                  f"rel_err={rel_err:.4f}% [{status}]")
    else:
        if col_ref is not None:
            print(f"  FAIL: T_dmeff column not found in new output")
            passed = False
        else:
            print(f"  SKIP: T_dmeff column not in reference (vanilla test?)")

    return passed


def check_thermodynamics(test_name, output_dir, ref_dir, tolerance):
    """Compare thermodynamics output: T_dmeff, T_b, rates."""
    print("\n--- Thermodynamics check ---")
    ref_file = os.path.join(ref_dir, f"{test_name}_thermodynamics.dat")
    new_file = os.path.join(output_dir, f"{test_name}_thermodynamics.dat")

    ref_names, ref_data = load_class_output(ref_file)
    new_names, new_data = load_class_output(new_file)

    if ref_data is None:
        print(f"  SKIP: Reference file not found: {ref_file}")
        return True
    if new_data is None:
        print(f"  FAIL: Output file not found: {new_file}")
        return False

    passed = True
    z_ref = ref_data[:, 0]
    z_new = new_data[:, 0]
    z_checks = [10, 100, 1000]

    # Check Tb (baryon temperature)
    col_ref = find_column(ref_names, 'Tb')
    col_new = find_column(new_names, 'Tb')
    if col_ref is not None and col_new is not None:
        for z_val in z_checks:
            t_ref = np.interp(z_val, z_ref, ref_data[:, col_ref])
            t_new = np.interp(z_val, z_new, new_data[:, col_new])
            rel_err = abs(t_new - t_ref) / max(abs(t_ref), 1e-30) * 100
            status = "PASS" if rel_err < tolerance else "FAIL"
            if status == "FAIL":
                passed = False
            print(f"  Tb(z={z_val:>6}):      ref={t_ref:.6e}, new={t_new:.6e}, "
                  f"rel_err={rel_err:.4f}% [{status}]")

    # Check T_dmeff
    col_ref = find_column(ref_names, 'T_dmeff')
    col_new = find_column(new_names, 'T_dmeff')
    if col_ref is not None and col_new is not None:
        for z_val in z_checks:
            t_ref = np.interp(z_val, z_ref, ref_data[:, col_ref])
            t_new = np.interp(z_val, z_new, new_data[:, col_new])
            rel_err = abs(t_new - t_ref) / max(abs(t_ref), 1e-30) * 100
            status = "PASS" if rel_err < tolerance else "FAIL"
            if status == "FAIL":
                passed = False
            print(f"  T_dmeff(z={z_val:>6}):  ref={t_ref:.6e}, new={t_new:.6e}, "
                  f"rel_err={rel_err:.4f}% [{status}]")
    else:
        if col_ref is not None:
            print(f"  FAIL: T_dmeff column not found in new output")
            passed = False
        else:
            print(f"  SKIP: T_dmeff column not in reference (vanilla test?)")

    # Check rate_dmeff_mom
    col_ref = find_column(ref_names, 'rate_dmeff_mom')
    col_new = find_column(new_names, 'rate_dmeff_mom')
    if col_ref is not None and col_new is not None:
        for z_val in z_checks:
            r_ref = np.interp(z_val, z_ref, ref_data[:, col_ref])
            r_new = np.interp(z_val, z_new, new_data[:, col_new])
            if abs(r_ref) > 1e-30:
                rel_err = abs(r_new - r_ref) / abs(r_ref) * 100
                status = "PASS" if rel_err < tolerance else "FAIL"
                if status == "FAIL":
                    passed = False
                print(f"  rate_mom(z={z_val:>6}): ref={r_ref:.6e}, new={r_new:.6e}, "
                      f"rel_err={rel_err:.4f}% [{status}]")
            else:
                print(f"  rate_mom(z={z_val:>6}): ref={r_ref:.6e}, new={r_new:.6e} "
                      f"(ref near zero, skip rel check)")
    else:
        if col_ref is not None:
            print(f"  FAIL: rate_dmeff_mom column not found in new output")
            passed = False

    return passed


def check_cl(test_name, output_dir, ref_dir, tolerance, compare_to=None):
    """Compare CMB C_l^TT at representative multipoles."""
    print("\n--- C_l check ---")
    ref_test = compare_to if compare_to else test_name
    ref_file = os.path.join(ref_dir, f"{ref_test}_cl.dat")
    new_file = os.path.join(output_dir, f"{test_name}_cl.dat")

    ref_names, ref_data = load_class_output(ref_file)
    new_names, new_data = load_class_output(new_file)

    if ref_data is None:
        print(f"  SKIP: Reference file not found: {ref_file}")
        return True
    if new_data is None:
        print(f"  FAIL: Output file not found: {new_file}")
        return False

    passed = True
    ell_col = find_column(ref_names, 'l')
    if ell_col is None:
        ell_col = 0
    tt_col_ref = find_column(ref_names, 'TT')
    tt_col_new = find_column(new_names, 'TT')
    if tt_col_ref is None:
        tt_col_ref = 1
    if tt_col_new is None:
        tt_col_new = 1

    ell_ref = ref_data[:, ell_col]
    ell_new = new_data[:, ell_col]
    tt_ref = ref_data[:, tt_col_ref]
    tt_new = new_data[:, tt_col_new]

    ell_checks = [10, 100, 500, 1000, 2000]
    for ell_val in ell_checks:
        c_ref = np.interp(ell_val, ell_ref, tt_ref)
        c_new = np.interp(ell_val, ell_new, tt_new)
        rel_err = abs(c_new - c_ref) / max(abs(c_ref), 1e-30) * 100
        status = "PASS" if rel_err < tolerance else "FAIL"
        if status == "FAIL":
            passed = False
        print(f"  C_l^TT(l={ell_val:>5}):  ref={c_ref:.6e}, new={c_new:.6e}, "
              f"rel_err={rel_err:.4f}% [{status}]")

    # Also compute max relative error over all ell
    common_ells = np.intersect1d(ell_ref.astype(int), ell_new.astype(int))
    if len(common_ells) > 0:
        ref_interp = np.interp(common_ells, ell_ref, tt_ref)
        new_interp = np.interp(common_ells, ell_new, tt_new)
        mask = np.abs(ref_interp) > 1e-30
        if mask.any():
            max_rel = np.max(np.abs(new_interp[mask] - ref_interp[mask]) / np.abs(ref_interp[mask])) * 100
            print(f"  Max rel error over all l: {max_rel:.4f}%")

    return passed


def check_pk(test_name, output_dir, ref_dir, tolerance, compare_to=None):
    """Compare matter power spectrum P(k) at representative wavenumbers."""
    print("\n--- P(k) check ---")
    ref_test = compare_to if compare_to else test_name
    ref_file = os.path.join(ref_dir, f"{ref_test}_pk.dat")
    new_file = os.path.join(output_dir, f"{test_name}_pk.dat")

    ref_names, ref_data = load_class_output(ref_file)
    new_names, new_data = load_class_output(new_file)

    if ref_data is None:
        print(f"  SKIP: Reference file not found: {ref_file}")
        return True
    if new_data is None:
        print(f"  FAIL: Output file not found: {new_file}")
        return False

    passed = True
    k_ref = ref_data[:, 0]
    k_new = new_data[:, 0]
    pk_ref = ref_data[:, 1]
    pk_new = new_data[:, 1]

    # Use log interpolation for P(k)
    logk_ref = np.log10(k_ref)
    logk_new = np.log10(k_new)
    logpk_ref = np.log10(pk_ref)
    logpk_new = np.log10(pk_new)

    k_checks = [0.01, 0.1, 1.0]
    for k_val in k_checks:
        logk_val = np.log10(k_val)
        if logk_val < logk_ref[0] or logk_val > logk_ref[-1]:
            print(f"  P(k={k_val:.3f}): SKIP (outside range)")
            continue
        p_ref = 10**np.interp(logk_val, logk_ref, logpk_ref)
        p_new = 10**np.interp(logk_val, logk_new, logpk_new)
        rel_err = abs(p_new - p_ref) / max(abs(p_ref), 1e-30) * 100
        status = "PASS" if rel_err < tolerance else "FAIL"
        if status == "FAIL":
            passed = False
        print(f"  P(k={k_val:.3f} h/Mpc): ref={p_ref:.6e}, new={p_new:.6e}, "
              f"rel_err={rel_err:.4f}% [{status}]")

    # Max relative error
    k_min = max(k_ref[0], k_new[0])
    k_max = min(k_ref[-1], k_new[-1])
    k_common = np.logspace(np.log10(k_min), np.log10(k_max), 200)
    ref_interp = 10**np.interp(np.log10(k_common), logk_ref, logpk_ref)
    new_interp = 10**np.interp(np.log10(k_common), logk_new, logpk_new)
    mask = ref_interp > 1e-30
    if mask.any():
        max_rel = np.max(np.abs(new_interp[mask] - ref_interp[mask]) / ref_interp[mask]) * 100
        print(f"  Max rel error over all k: {max_rel:.4f}%")

    return passed


def check_tk(test_name, output_dir, ref_dir, tolerance):
    """Compare transfer functions, including d_dmeff if present."""
    print("\n--- Transfer function check ---")
    ref_file = os.path.join(ref_dir, f"{test_name}_tk.dat")
    new_file = os.path.join(output_dir, f"{test_name}_tk.dat")

    ref_names, ref_data = load_class_output(ref_file)
    new_names, new_data = load_class_output(new_file)

    if ref_data is None:
        print(f"  SKIP: Reference file not found: {ref_file}")
        return True
    if new_data is None:
        print(f"  FAIL: Output file not found: {new_file}")
        return False

    passed = True
    k_ref = ref_data[:, 0]
    k_new = new_data[:, 0]

    # Check d_dmeff if present
    col_ref = find_column(ref_names, 'd_dmeff')
    col_new = find_column(new_names, 'd_dmeff')
    if col_ref is not None and col_new is not None:
        logk_ref = np.log10(k_ref)
        logk_new = np.log10(k_new)
        k_checks = [0.01, 0.1, 1.0]
        for k_val in k_checks:
            logk_val = np.log10(k_val)
            if logk_val < logk_ref[0] or logk_val > logk_ref[-1]:
                print(f"  d_dmeff(k={k_val:.3f}): SKIP (outside range)")
                continue
            d_ref = np.interp(logk_val, logk_ref, ref_data[:, col_ref])
            d_new = np.interp(logk_val, logk_new, new_data[:, col_new])
            if abs(d_ref) > 1e-30:
                rel_err = abs(d_new - d_ref) / abs(d_ref) * 100
                status = "PASS" if rel_err < tolerance else "FAIL"
                if status == "FAIL":
                    passed = False
                print(f"  d_dmeff(k={k_val:.3f}): ref={d_ref:.6e}, new={d_new:.6e}, "
                      f"rel_err={rel_err:.4f}% [{status}]")
            else:
                print(f"  d_dmeff(k={k_val:.3f}): ref near zero, skip rel check")
    else:
        if col_ref is not None:
            print(f"  FAIL: d_dmeff column not found in new output")
            passed = False
        else:
            print(f"  SKIP: d_dmeff column not in reference (vanilla test?)")

    # Check d_b (baryon transfer)
    col_ref = find_column(ref_names, 'd_b')
    col_new = find_column(new_names, 'd_b')
    if col_ref is not None and col_new is not None:
        logk_ref = np.log10(k_ref)
        logk_new = np.log10(k_new)
        for k_val in [0.01, 0.1]:
            logk_val = np.log10(k_val)
            d_ref = np.interp(logk_val, logk_ref, ref_data[:, col_ref])
            d_new = np.interp(logk_val, logk_new, new_data[:, col_new])
            if abs(d_ref) > 1e-30:
                rel_err = abs(d_new - d_ref) / abs(d_ref) * 100
                status = "PASS" if rel_err < tolerance else "FAIL"
                if status == "FAIL":
                    passed = False
                print(f"  d_b(k={k_val:.3f}):     ref={d_ref:.6e}, new={d_new:.6e}, "
                      f"rel_err={rel_err:.4f}% [{status}]")

    return passed


def main():
    parser = argparse.ArgumentParser(
        description="Compare CLASS output files against reference baselines.")
    parser.add_argument('test_name', help='Name of the test case (e.g., test_coulomb)')
    parser.add_argument('--check', default='all',
                        help='Comma-separated checks: background,thermodynamics,cl,pk,tk,all')
    parser.add_argument('--tolerance', type=float, default=0.1,
                        help='Tolerance in percent (default: 0.1%%)')
    parser.add_argument('--compare-to', default=None,
                        help='Compare against a different test name (for gauge checks)')
    parser.add_argument('--output-dir', default=None,
                        help='Directory containing new output files (default: output/)')
    parser.add_argument('--reference-dir', default=None,
                        help='Directory containing reference files')
    args = parser.parse_args()

    # Determine directories
    script_dir = os.path.dirname(os.path.abspath(__file__))
    class_dir = os.path.dirname(script_dir)

    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.join(class_dir, 'output')

    if args.reference_dir:
        ref_dir = args.reference_dir
    else:
        ref_test = args.compare_to if args.compare_to else args.test_name
        ref_dir = os.path.join(script_dir, 'reference', ref_test)

    print(f"Test: {args.test_name}")
    print(f"Reference dir: {ref_dir}")
    print(f"Output dir: {output_dir}")
    print(f"Tolerance: {args.tolerance}%")

    checks = args.check.lower().split(',')
    if 'all' in checks:
        checks = ['background', 'thermodynamics', 'cl', 'pk', 'tk']

    all_passed = True

    for check in checks:
        check = check.strip()
        if check == 'background':
            ok = check_background(args.test_name, output_dir, ref_dir, args.tolerance)
        elif check == 'thermodynamics':
            ok = check_thermodynamics(args.test_name, output_dir, ref_dir, args.tolerance)
        elif check == 'cl':
            ok = check_cl(args.test_name, output_dir, ref_dir, args.tolerance, args.compare_to)
        elif check == 'pk':
            ok = check_pk(args.test_name, output_dir, ref_dir, args.tolerance, args.compare_to)
        elif check == 'tk':
            ok = check_tk(args.test_name, output_dir, ref_dir, args.tolerance)
        else:
            print(f"\nUnknown check: {check}")
            ok = False

        if not ok:
            all_passed = False

    print("\n" + "=" * 50)
    if all_passed:
        print(f"RESULT: ALL CHECKS PASSED for {args.test_name}")
    else:
        print(f"RESULT: SOME CHECKS FAILED for {args.test_name}")
    print("=" * 50)

    sys.exit(0 if all_passed else 1)


if __name__ == '__main__':
    main()
