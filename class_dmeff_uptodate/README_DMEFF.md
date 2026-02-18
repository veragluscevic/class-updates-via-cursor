# CLASS with dmeff: Dark Matter–Baryon Interactions

This is **CLASS v3.3.4** extended with the **dmeff** module for modeling
velocity-dependent dark matter–baryon scattering. The dmeff implementation was
ported from an earlier CLASS v2.9.4 codebase and validated to < 0.1% accuracy
against the original.


## Overview

Standard cold dark matter (CDM) is collisionless — it interacts with ordinary
matter only through gravity. The dmeff extension replaces some or all of the CDM
with a dark matter species that scatters off baryons (hydrogen, helium, or
electrons) with a momentum-transfer cross section that follows a power law in
relative velocity:

    sigma(v) = sigma_0 * (v / c)^n

This models a broad class of dark matter interactions:

| Power law n | Physical model                        | Example          |
|-------------|---------------------------------------|------------------|
| n = -4      | Coulomb-like / electric dipole        | Millicharged DM  |
| n = -2      | Light mediator                        | Dark photon      |
| n = 0       | Contact interaction (constant sigma)  | Heavy mediator   |
| n = +2      | p-wave scattering                     |                  |
| n = +4      | d-wave scattering                     |                  |


## Building

```bash
cd class_dmeff_uptodate
make clean && make
```

The build produces the `class` executable and `libclass.a` library. No
additional dependencies beyond standard CLASS requirements are needed.

To build the Python wrapper:

```bash
cd python
python setup.py install
```


## Quick Start: Running with dmeff

### From the command line

Create an INI file with dmeff parameters (see `test_dmeff/test_coulomb.ini`
for a complete example):

```ini
# Cosmological parameters
h = 0.6732
T_cmb = 2.7255
omega_b = 0.022032
omega_cdm = 0.0          # set to 0 if all DM is dmeff
omega_dmeff = 0.12        # dmeff density parameter

# dmeff interaction parameters
m_dmeff = 0.01            # dark matter mass in GeV
N_dmeff = 1               # number of interaction terms
sigma_dmeff = 1e-41       # cross section in cm^2
npow_dmeff = -4           # velocity power law index (>= -4)
dmeff_target = hydrogen   # target: baryon, hydrogen, helium, or electron
Vrel_dmeff = 29.0         # initial bulk velocity in km/s

# Standard CLASS output
output = tCl,mPk,mTk
```

Run:

```bash
./class your_dmeff_run.ini
```

### From Python

```python
from classy import Class

cosmo = Class()
cosmo.set({
    'h': 0.6732,
    'T_cmb': 2.7255,
    'omega_b': 0.022032,
    'omega_cdm': 0.0,
    'omega_dmeff': 0.12,
    'm_dmeff': 0.01,          # GeV
    'N_dmeff': 1,
    'sigma_dmeff': 1e-41,     # cm^2
    'npow_dmeff': -4,
    'dmeff_target': 'hydrogen',
    'Vrel_dmeff': 29.0,       # km/s
    'output': 'tCl,mPk',
})
cosmo.compute()

# Access results
cl = cosmo.lensed_cl(2500)
pk = cosmo.pk(0.1, 0.0)  # P(k=0.1/Mpc, z=0)

# Access dmeff-specific derived quantities
print("Omega_dmeff:", cosmo.Omega0_dmeff())
derived = cosmo.get_current_derived_parameters(['z_dmeff_decoupling'])
print("z_dmeff_decoupling:", derived['z_dmeff_decoupling'])

cosmo.struct_cleanup()
```


## dmeff Parameters Reference

### Density specification (choose one)

| Parameter     | Description                               | Units    | Default |
|---------------|-------------------------------------------|----------|---------|
| `Omega_dmeff` | Density parameter Omega_dmeff             | —        | 0       |
| `omega_dmeff` | Physical density omega_dmeff = Omega*h^2  | —        | 0       |
| `f_dmeff`     | Fraction of CDM that is dmeff (0 to 1)    | —        | 0       |

Only one of these three may be specified. If `f_dmeff` is used, it
automatically adjusts `omega_cdm` downward to keep the total dark matter
density fixed.

### Interaction parameters

| Parameter      | Description                               | Units    | Default    |
|----------------|-------------------------------------------|----------|------------|
| `m_dmeff`      | Dark matter particle mass                 | GeV      | 1.0        |
| `N_dmeff`      | Number of interaction terms               | —        | 0          |
| `sigma_dmeff`  | Cross section prefactor(s)                | cm^2     | 0          |
| `log10sigma_dmeff` | Log10 of cross section (alternative)  | log10(cm^2) | —       |
| `npow_dmeff`   | Velocity power law index(es)              | —        | 0          |
| `dmeff_target` | Scattering target(s)                      | string   | hydrogen   |
| `Vrel_dmeff`   | Initial DM–baryon bulk velocity           | km/s     | 0          |

When `N_dmeff > 1`, the parameters `sigma_dmeff`, `npow_dmeff`, and
`dmeff_target` are comma-separated lists. Example for two channels:

```ini
N_dmeff = 2
sigma_dmeff = 1e-41,1e-35
npow_dmeff = -4,-2
dmeff_target = hydrogen,electron
```

### Constraints

- `npow_dmeff` must be >= -4 for each term
- `sigma_dmeff` and `log10sigma_dmeff` are mutually exclusive
- All array parameters must have exactly `N_dmeff` entries
- When using `f_dmeff`, `omega_cdm` must be nonzero


## Output Quantities

### Background output (`write background = yes`)

| Column          | Description                                            |
|-----------------|--------------------------------------------------------|
| `rho_dmeff`     | dmeff energy density [units of 3c^2 H0^2 / 8pi G]    |
| `T_dmeff`       | dmeff temperature [K]                                  |
| `Vrel_dmeff`    | DM–baryon relative bulk velocity [c]                   |
| `dkappa_dmeff`  | Momentum exchange rate [Mpc^-1]                        |
| `dkappaT_dmeff` | Heat exchange rate [Mpc^-1]                            |
| `cdmeff2`       | dmeff sound speed squared [c^2]                        |

### Thermodynamics output (`write thermodynamics = yes`)

| Column            | Description                                          |
|-------------------|------------------------------------------------------|
| `T_dmeff`         | dmeff temperature [K]                                |
| `rate_dmeff_mom`  | Momentum exchange rate [Mpc^-1]                      |
| `rate_dmeff_mom'` | Derivative of momentum exchange rate                 |
| `rate_dmeff_temp` | Heat exchange rate [Mpc^-1]                          |
| `c_dmeff^2`       | dmeff sound speed squared [c^2]                      |

### Transfer functions (`output = mTk`)

| Column     | Description                                    |
|------------|------------------------------------------------|
| `d_dmeff`  | dmeff density transfer function                |
| `t_dmeff`  | dmeff velocity transfer function               |

### Derived parameters (Python wrapper)

| Method / key          | Description                              |
|-----------------------|------------------------------------------|
| `Omega0_dmeff()`      | Present-day dmeff density parameter      |
| `z_dmeff_decoupling`  | Redshift where dmeff thermally decouples |


## Physics Summary

The dmeff module adds three coupled physical effects:

### 1. Temperature evolution

The dmeff temperature T_chi evolves according to:

    dT_chi/dtau = -2 a H T_chi + 2 R_heat (T_b - T_chi)

where R_heat is the heat exchange rate summed over all interaction terms. At
early times when R_heat >> a*H, the dmeff temperature is tightly coupled to
the baryon temperature. At late times when the interaction rate drops below the
Hubble rate, T_chi cools adiabatically as a^(-2).

### 2. Perturbation equations

The dmeff density and velocity perturbations (delta_chi, theta_chi) evolve as:

    delta_chi' = -(theta_chi + metric_continuity)
    theta_chi' = -a'/a * theta_chi + k^2 * c_chi^2 * delta_chi
                 + metric_euler + R_mom * (theta_b - theta_chi)

where R_mom is the momentum exchange rate and c_chi^2 is the dmeff sound
speed. The drag term R_mom * (theta_b - theta_chi) couples dmeff to baryons.

### 3. Back-reaction on baryons

The baryons experience a reciprocal drag force:

    theta_b' += (rho_chi / rho_b) * R_mom * (theta_chi - theta_b)

During tight coupling (before recombination), this modifies the TCA
denominator from 1/(1+R) to 1/(1+R+beta*R), where beta encodes the relative
strength of dmeff-baryon and photon-baryon coupling.


## Validating the Port

To verify the dmeff implementation is working correctly:

```bash
# Build CLASS
make clean && make

# Run the full validation test suite
bash test_dmeff/run_all_tests.sh
```

You should see `OVERALL: 6/6 tests PASSED`. Any individual check reporting
< 0.1% relative error is a pass. Differences of ~0.20% at l=2000 and ~0.30%
in T_b(z=10) are expected CLASS version baselines (HyRec update), not dmeff
errors — they appear identically in the vanilla test.

See `test_dmeff/README.md` for detailed verification procedures, acceptance
criteria, and interpretation guidance.


## Test Scripts

Four scripts in `test_dmeff/` support validation and reproducibility. All
are run from the `class_dmeff_uptodate/` directory unless noted otherwise.

### generate_test_data.sh — Recreate everything from scratch

```bash
# Run from the REPOSITORY ROOT (one level above class_dmeff_uptodate/)
bash class_dmeff_uptodate/test_dmeff/generate_test_data.sh
```

This is the most important script for reproducibility. It regenerates all
test INI files and reference outputs from scratch, with no dependency on
pre-existing files. It:

1. **Creates 7 test INI files** in both `class_dmeff/test_dmeff/` and
   `class_dmeff_uptodate/test_dmeff/`, with all parameter values embedded
   directly in the script (6 physics tests + 1 Newtonian gauge variant).
2. **Builds `class_dmeff/`** (the original CLASS v2.9.4 with working dmeff)
   if the executable does not already exist.
3. **Runs all 6 test cases** in `class_dmeff/` to generate ground-truth
   reference outputs.
4. **Copies the 30 reference `.dat` files** into
   `class_dmeff_uptodate/test_dmeff/reference/<test_name>/`.

**Prerequisite**: The `class_dmeff/` directory must contain the original
CLASS v2.9.4 source code, and a C compiler must be available.

**When to use**: Run this if the test INI files or reference outputs are
missing (e.g., after a fresh clone, or if untracked files were cleaned up).

### run_all_tests.sh — One-command validation

```bash
bash test_dmeff/run_all_tests.sh
```

Runs all 6 test cases end-to-end: deletes stale output files, executes
the CLASS binary, compares each output against reference baselines using
`compare_outputs.py`, and prints a color-coded summary. Exit code 0 means
all tests passed; exit code 1 means at least one failed.

### compare_outputs.py — Detailed per-test comparison

```bash
# All checks for one test
python test_dmeff/compare_outputs.py test_coulomb --check all

# Specific checks with custom tolerance
python test_dmeff/compare_outputs.py test_coulomb --check cl,pk --tolerance 0.05

# Gauge consistency (compare two different runs)
python test_dmeff/compare_outputs.py test_coulomb_newtonian \
    --compare-to test_coulomb --check cl,pk --tolerance 0.01
```

Compares CLASS output files in `output/` against reference files in
`test_dmeff/reference/`. Reports relative errors at specific check points
(z values for thermodynamics, l values for C_l, k values for P(k)) and
prints PASS/FAIL for each. Available checks: `background`,
`thermodynamics`, `cl`, `pk`, `tk`, or `all`.

### example_dmeff.py — Python wrapper demonstration

```bash
python test_dmeff/example_dmeff.py
```

Demonstrates using the dmeff extension from Python. Runs three models
(vanilla Lambda-CDM, Coulomb-like dmeff, constant cross-section dmeff),
compares CMB and matter power spectra side by side, accesses dmeff-specific
derived parameters (`Omega0_dmeff()`, `z_dmeff_decoupling`), and runs
self-verification checks. Requires the Python wrapper to be installed
(`cd python && python setup.py install`).


## Key Files Modified

The dmeff extension modifies the following CLASS source files:

| File                          | Changes                                      |
|-------------------------------|----------------------------------------------|
| `include/background.h`       | dmeff fields and background indices          |
| `include/thermodynamics.h`   | dmeff target enum, thermo indices, functions |
| `include/perturbations.h`    | dmeff perturbation and source indices        |
| `include/precisions.h`       | dmeff integration tolerance                  |
| `source/input.c`             | dmeff parameter parsing and validation       |
| `source/background.c`        | dmeff density, placeholder T_dmeff values    |
| `source/thermodynamics.c`    | Temperature ODE, rate computation, write-back|
| `source/perturbations.c`     | Evolution equations, TCA modifications       |
| `source/output.c`            | Thermodynamics output column headers         |
| `source/fourier.c`           | Halofit/HMcode warning for dmeff             |
| `python/cclassy.pxd`         | Cython struct declarations                   |
| `python/classy.pyx`          | Python wrapper methods                       |
| `explanatory.ini`            | Parameter documentation                      |


## References

The dmeff formalism is based on:

- Gluscevic & Boddy (2018), "Constraints on Scattering of keV–TeV Dark Matter
  with Protons in the Early Universe", Physical Review Letters 121, 081301
- Boddy & Gluscevic (2018), "First Cosmological Constraint on the Effective
  Theory of Dark Matter-Proton Interactions", Physical Review D 98, 083510
- Boddy, Gluscevic, et al. (2018), "Critical assessment of CMB limits on dark
  matter-baryon scattering: New treatment of the relative bulk velocity",
  Physical Review D 98, 123506
