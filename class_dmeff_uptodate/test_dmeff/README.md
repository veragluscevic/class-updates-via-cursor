# dmeff Validation Test Suite

This directory contains the validation test suite for the **dmeff** (dark matter
with effective interactions) extension to CLASS v3.3.4. The tests verify that
the dmeff physics ported from CLASS v2.9.4 produces results matching the
original implementation to better than 0.1% relative accuracy.


## Quick Start

From the `class_dmeff_uptodate/` directory:

```bash
# 1. Build CLASS
make clean && make

# 2. Run all tests and verify results (recommended)
bash test_dmeff/run_all_tests.sh

# 3. Or run a single test manually
rm -f output/test_coulomb_*
./class test_dmeff/test_coulomb.ini
python test_dmeff/compare_outputs.py test_coulomb --check all
```


## What is dmeff?

The **dmeff** module adds a dark matter component that scatters with baryons
through a velocity-dependent cross section:

    sigma(v) = sigma_0 * (v / c)^n

where `sigma_0` is the constant cross-section prefactor (`sigma_dmeff`, in
cm^2) and `n` is the velocity power-law index (`npow_dmeff`). The scattering
target can be the mean baryon, hydrogen, helium, or free electrons.

This interaction produces three observable effects:

1. **Dark matter temperature evolution**: dmeff acquires a temperature T_dmeff
   that evolves separately from the baryon temperature, controlled by heat
   exchange with baryons.

2. **Modified perturbation growth**: The momentum exchange between dmeff and
   baryons creates a drag force that suppresses structure growth, particularly
   at small scales.

3. **CMB power spectrum changes**: The baryon-dmeff coupling modifies the
   baryon velocity perturbations, which in turn affect the photon-baryon fluid
   and leave imprints on the CMB angular power spectra.


## Test Cases

Six test configurations are provided, spanning different physics regimes.
All share common cosmological parameters:

| Parameter  | Value    |
|------------|----------|
| h          | 0.6732   |
| T_cmb      | 2.7255 K |
| omega_b    | 0.022032 |
| A_s        | 2.1e-9   |
| n_s        | 0.9660   |
| tau_reio   | 0.0543   |

### test_coulomb — Coulomb-like interaction

Strong velocity-dependent interaction with n = -4 (Coulomb-like), targeting
hydrogen. All dark matter is dmeff (omega_cdm = 0).

| Parameter      | Value     |
|----------------|-----------|
| omega_dmeff    | 0.12      |
| omega_cdm      | 0.0       |
| m_dmeff        | 0.01 GeV  |
| sigma_dmeff    | 1e-41 cm^2|
| npow_dmeff     | -4        |
| dmeff_target   | hydrogen  |
| Vrel_dmeff     | 29.0 km/s |

**What to look for**: Large suppression of small-scale P(k) relative to
standard CDM. T_dmeff stays coupled to baryons to late times.

### test_constant — Velocity-independent scattering

Constant cross section (n = 0), targeting the mean baryon.

| Parameter      | Value     |
|----------------|-----------|
| omega_dmeff    | 0.12      |
| omega_cdm      | 0.0       |
| m_dmeff        | 0.01 GeV  |
| sigma_dmeff    | 1e-30 cm^2|
| npow_dmeff     | 0         |
| dmeff_target   | baryon    |

**What to look for**: T_dmeff is much lower than Coulomb case at the same
redshifts, since constant cross section gives weaker coupling at high z.

### test_electron — Electron target

Interaction with free electrons (n = -2), relevant for millicharged dark
matter models.

| Parameter      | Value     |
|----------------|-----------|
| omega_dmeff    | 0.12      |
| omega_cdm      | 0.0       |
| m_dmeff        | 0.01 GeV  |
| sigma_dmeff    | 1e-35 cm^2|
| npow_dmeff     | -2        |
| dmeff_target   | electron  |

**What to look for**: Electron target coupling depends on the free electron
fraction x_e, so the interaction rate is sensitive to recombination history.

### test_mixed — Partial CDM replacement

Half CDM and half dmeff, testing coexistence of collisionless and interacting
dark matter.

| Parameter      | Value     |
|----------------|-----------|
| omega_dmeff    | 0.06      |
| omega_cdm      | 0.06      |
| m_dmeff        | 0.01 GeV  |
| sigma_dmeff    | 1e-41 cm^2|
| npow_dmeff     | -4        |
| dmeff_target   | hydrogen  |

**What to look for**: P(k) suppression is intermediate between vanilla CDM
and full-dmeff cases. Both CDM and dmeff perturbation variables are active.

### test_multi — Multiple interaction terms

Two simultaneous interaction channels (hydrogen + electron) testing the array
parameter handling with N_dmeff = 2.

| Parameter      | Value               |
|----------------|---------------------|
| omega_dmeff    | 0.12                |
| omega_cdm      | 0.0                 |
| m_dmeff        | 0.01 GeV            |
| N_dmeff        | 2                   |
| sigma_dmeff    | 1e-41,1e-35         |
| npow_dmeff     | -4,-2               |
| dmeff_target   | hydrogen,electron   |

**What to look for**: T_dmeff should be close to (but not identical to) the
single-channel Coulomb case, since the hydrogen channel dominates.

### test_vanilla — Baseline (no dmeff)

Standard Lambda-CDM with no dmeff, confirming that the dmeff code does not
alter standard CLASS behavior.

| Parameter      | Value     |
|----------------|-----------|
| omega_cdm      | 0.12      |
| N_dmeff        | 0         |

**What to look for**: C_l and P(k) should match the reference outputs from
the original CLASS v2.9.4 to the same precision as all other tests (i.e., the
~0.20% difference at l=2000 is a CLASS version baseline, not dmeff-related).


## Reference Data

Reference outputs are stored in `reference/<test_name>/` subdirectories. These
were generated by running the same test INI files in `class_dmeff/` (the
original CLASS v2.9.4 with working dmeff). Each test produces five output files:

| File suffix         | Contents                                         |
|---------------------|--------------------------------------------------|
| `_background.dat`   | Background quantities: rho_dmeff, T_dmeff, V_rel |
| `_thermodynamics.dat`| T_dmeff(z), T_b(z), interaction rates            |
| `_cl.dat`           | CMB angular power spectra C_l^TT, C_l^EE, etc.  |
| `_pk.dat`           | Matter power spectrum P(k) at z=0                |
| `_tk.dat`           | Transfer functions including d_dmeff              |


## Comparison Script

The `compare_outputs.py` script compares CLASS output files against reference
baselines and reports pass/fail at a configurable tolerance (default 0.1%).

### Usage

```bash
# Compare all observables for a single test
python test_dmeff/compare_outputs.py test_coulomb --check all

# Compare only CMB and P(k)
python test_dmeff/compare_outputs.py test_coulomb --check cl,pk

# Compare with a stricter tolerance
python test_dmeff/compare_outputs.py test_coulomb --check all --tolerance 0.05

# Gauge consistency: compare Newtonian vs synchronous gauge outputs
python test_dmeff/compare_outputs.py test_coulomb_newtonian \
    --compare-to test_coulomb --check cl,pk --tolerance 0.01
```

### Available checks

| Check            | What it compares                                   |
|------------------|----------------------------------------------------|
| `background`     | rho_dmeff(z), T_dmeff(z) from background file      |
| `thermodynamics` | T_dmeff(z), T_b(z), interaction rates              |
| `cl`             | C_l^TT at l = 10, 100, 500, 1000, 2000            |
| `pk`             | P(k) at k = 0.01, 0.1, 1.0 h/Mpc                 |
| `tk`             | Transfer functions d_dmeff, d_b                    |
| `all`            | All of the above                                   |


## How to Verify Accuracy

The central question for this project is: **does the ported dmeff code in
CLASS v3.3.4 reproduce the same physics as the original CLASS v2.9.4
implementation?**

### Step-by-step verification procedure

1. **Build from clean state**:

   ```bash
   cd class_dmeff_uptodate
   make clean && make
   ```

2. **Run the full test suite** (this deletes stale output files automatically):

   ```bash
   bash test_dmeff/run_all_tests.sh
   ```

3. **Check the summary**. You should see output like:

   ```
   ============================
   OVERALL: 6/6 tests PASSED
   ============================
   ```

4. **Interpret the results**:

   - **PASS** at all l <= 1000 and all k values means the port is correct.
   - **~0.20% at l=2000** is expected and acceptable. This is an irreducible
     baseline difference between CLASS v2.9.4 and v3.3.4 due to the HyRec
     recombination code update (HyRec -> HyRec2020). It appears identically
     in the vanilla test, confirming it is not a dmeff error.
   - **~0.30% in T_b at z=10** is similarly a CLASS version baseline
     difference, not a dmeff error.

5. **Gauge consistency check** (optional but recommended):

   ```bash
   rm -f output/test_coulomb_newtonian_*
   ./class test_dmeff/test_coulomb_newtonian.ini
   python test_dmeff/compare_outputs.py test_coulomb_newtonian \
       --compare-to test_coulomb --check cl,pk --tolerance 0.05
   ```

   CMB C_l should agree to < 0.03% between synchronous and Newtonian gauges
   at l <= 1000.

6. **Python wrapper verification** (optional):

   ```bash
   cd python && python setup.py install && cd ..
   python test_dmeff/example_dmeff.py
   ```

### What constitutes a successful port?

All of the following must hold:

| Observable     | Tolerance       | Check points                |
|----------------|-----------------|-----------------------------|
| T_dmeff(z)     | < 0.1%          | z = 10, 100, 1000          |
| T_b(z)         | < 0.1% (z>=100) | z = 100, 1000               |
| C_l^TT         | < 0.1%          | l = 10, 100, 500, 1000     |
| P(k)           | < 0.1%          | k = 0.01, 0.1, 1.0 h/Mpc  |
| d_dmeff(k)     | < 0.1%          | k = 0.01, 0.1, 1.0 h/Mpc  |
| Vanilla test   | Same precision  | All of the above            |

**Note**: The tolerance applies at the check points listed. At l=2000 and
T_b(z=10), the ~0.2-0.3% differences are CLASS version baselines and should
not be considered failures.


## Known Limitations and Caveats

1. **HyRec version difference**: CLASS v3.3.4 uses HyRec2020 while v2.9.4
   uses an older HyRec version. This creates an irreducible ~0.20% difference
   in C_l at l=2000 and ~0.30% in T_b at z=10 for ALL cross-version
   comparisons, including vanilla (no dmeff). This is not a dmeff porting
   error.

2. **Nonlinear corrections**: Halofit and HMcode are calibrated for
   collisionless CDM. If dmeff is enabled with nonlinear corrections, CLASS
   will print a warning. Results in the nonlinear regime should be treated
   with caution.

3. **Electron target rate sensitivity**: The electron-target interaction rate
   depends on the free electron fraction x_e. Since HyRec2020 gives slightly
   different x_e than the older HyRec, the electron-target rate can show up
   to ~1.6% difference at z=10 (where x_e is very small). This is a
   recombination code difference, not a dmeff error.

4. **Valgrind**: Memory leak testing with valgrind is not available on macOS.
   On Linux, run:
   ```bash
   valgrind --leak-check=full ./class test_dmeff/test_coulomb.ini
   ```


## Regenerating Test Data from Scratch

The test INI files and reference outputs are not committed to the repository.
If they are lost, regenerate everything with a single command from the
**repository root**:

```bash
bash class_dmeff_uptodate/test_dmeff/generate_test_data.sh
```

This script is fully self-contained. It:
1. Creates all test INI files from embedded parameter definitions (no external
   files needed)
2. Builds `class_dmeff/` (the original CLASS v2.9.4 with working dmeff)
3. Runs each test case to generate reference outputs
4. Copies reference outputs into `test_dmeff/reference/`

**Prerequisite**: The `class_dmeff/` directory must exist with the original
v2.9.4 source code and a working C compiler.

After regeneration, verify the port from `class_dmeff_uptodate/`:

```bash
make clean && make
bash test_dmeff/run_all_tests.sh
```


## File Inventory

```
test_dmeff/
├── README.md                  # This file
├── generate_test_data.sh      # Regenerate all INI files and reference outputs
├── run_all_tests.sh           # One-command test runner
├── compare_outputs.py         # Automated comparison script
├── example_dmeff.py           # Python wrapper usage example
├── test_coulomb.ini           # Coulomb-like interaction (npow=-4)
├── test_coulomb_newtonian.ini # Same as above in Newtonian gauge
├── test_constant.ini          # Constant cross section (npow=0)
├── test_electron.ini          # Electron target (npow=-2)
├── test_mixed.ini             # Partial CDM replacement
├── test_multi.ini             # Multiple interaction terms
├── test_vanilla.ini           # Baseline (no dmeff)
└── reference/                 # Reference outputs from CLASS v2.9.4
    ├── test_coulomb/          # 5 .dat files
    ├── test_constant/         # 5 .dat files
    ├── test_electron/         # 5 .dat files
    ├── test_mixed/            # 5 .dat files
    ├── test_multi/            # 5 .dat files
    └── test_vanilla/          # 5 .dat files
```
