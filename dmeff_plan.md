# Port dmeff Dark Matter-Baryon Interaction Physics to CLASS v3.3.4

This ExecPlan is a living document. The sections Progress, Surprises & Discoveries, Decision Log, and Outcomes & Retrospective must be kept up to date as work proceeds.

This document must be maintained in accordance with .agent/PLANS.md at the repository root.


## Purpose / Big Picture

After implementing this plan, CLASS v3.3.4 will support a new dark matter component called "dmeff" that interacts with baryons through velocity-dependent momentum transfer. Users will be able to model scenarios where some or all of the cold dark matter undergoes scattering with baryons, hydrogen, helium, or electrons, with configurable cross sections that follow power-law velocity dependence. This enables testing observational constraints on dark matter interactions from cosmic microwave background observations and large-scale structure data.

A user can verify this works by creating an input file that specifies dmeff parameters (such as Omega_dmeff equals 0.12, m_dmeff equals 1.0 GeV, N_dmeff equals 1, sigma_dmeff equals 1e-30 square centimeters, npow_dmeff equals negative four, and dmeff_target equals hydrogen), running the CLASS executable with that input file, and observing output files containing dmeff density evolution, temperature evolution, relative velocity, interaction rates, and transfer functions for the dmeff density and velocity perturbations. The results should match reference outputs generated from the existing class_dmeff implementation based on CLASS v2.9.4 to within 0.1 percent relative accuracy.


## Progress

- [ ] Phase 0: Create baseline validation tests from class_dmeff
  - [ ] Identify five representative test cases with varying dmeff parameters
  - [ ] Run test cases in class_dmeff directory and save reference outputs
  - [ ] Document expected outputs for background, thermodynamics, perturbations, and transfer functions
  - [ ] Create test_dmeff directory in class_dmeff_uptodate with test parameter files
- [ ] Phase 1: Port background module for dmeff density evolution
  - [ ] Add dmeff parameters to background structure in include/background.h
  - [ ] Add background indices for dmeff quantities
  - [ ] Modify background_init to allocate dmeff arrays and count quantities
  - [ ] Modify background_solve to compute dmeff density evolution
  - [ ] Modify background_functions to store dmeff quantities in output vector
  - [ ] Modify background_free to deallocate dmeff arrays
  - [ ] Add dmeff columns to background output
  - [ ] Validate: compare rho_dmeff(z) against baseline (< 0.01% error)
- [ ] Phase 2: Port thermodynamics module for temperature evolution and interaction rates
  - [ ] Add thermodynamics indices for dmeff temperature and rates
  - [ ] Add enum select_dmeff_target for target particle types
  - [ ] Add z_dmeff_decoupling derived parameter
  - [ ] Implement thermodynamics_dmeff_rate function for momentum and heat exchange
  - [ ] Implement thermodynamics_dmeff_temperature for ODE integration
  - [ ] Implement thermodynamics_dmeff_derivs for temperature derivatives
  - [ ] Modify thermodynamics_init to call temperature integration
  - [ ] Modify thermodynamics_at_z to interpolate dmeff quantities
  - [ ] Add dmeff columns to thermodynamics output
  - [ ] Validate: compare T_dmeff(z), dkappa_dmeff(z) against baseline (< 0.1% error)
- [ ] Phase 3: Port perturbations module for density and velocity evolution
  - [ ] Add perturbation indices for delta_dmeff and theta_dmeff
  - [ ] Add source function indices
  - [ ] Modify perturb_vector_init to allocate dmeff perturbation space
  - [ ] Modify perturb_initial_conditions to set dmeff initial values
  - [ ] Modify perturb_einstein to include dmeff in metric source terms
  - [ ] Implement gauge transformation functions for dmeff
  - [ ] Modify perturb_derivs to add dmeff evolution equations
  - [ ] Modify source functions to output dmeff perturbations
  - [ ] Validate: compare delta_dmeff(k,tau) and theta_dmeff(k,tau) against baseline
- [ ] Phase 4: Port input module for parameter parsing
  - [ ] Add parameter parsing for Omega_dmeff, omega_dmeff, f_dmeff
  - [ ] Add parsing for m_dmeff with unit conversion from GeV to kg
  - [ ] Add parsing for N_dmeff and array parameters
  - [ ] Add parsing for sigma_dmeff and log10sigma_dmeff arrays
  - [ ] Add parsing for npow_dmeff array with validation (>= -4)
  - [ ] Add parsing for dmeff_target string array
  - [ ] Add parsing for Vrel_dmeff with unit conversion from km/s to m/s
  - [ ] Implement validation for array length consistency
  - [ ] Set has_dmeff flag based on parameters
  - [ ] Update explanatory.ini with dmeff parameter documentation
  - [ ] Validate: test all parameter input formats and error handling
- [ ] Phase 5: Port transfer module
  - [ ] Add transfer indices for dmeff density and velocity
  - [ ] Modify transfer_init to count dmeff transfer types
  - [ ] Modify transfer_functions to extract and store dmeff transfers
  - [ ] Add dmeff columns to transfer output files
  - [ ] Validate: compare T_dmeff(k) against baseline
- [ ] Phase 6: Port fourier module warning (formerly nonlinear)
  - [ ] Locate warning code in class_dmeff/source/nonlinear.c
  - [ ] Port warning to appropriate location in source/fourier.c
  - [ ] Validate: verify warning appears when dmeff and Halofit/HMcode both enabled
- [ ] Phase 7: Port Python wrapper
  - [ ] Add Omega0_dmeff to background structure declaration in python/cclassy.pxd
  - [ ] Add z_dmeff_decoupling to thermodynamics structure declaration
  - [ ] Add Omega0_dmeff method to python/classy.pyx
  - [ ] Add z_dmeff_decoupling to derived parameters in get_current_derived_parameters
  - [ ] Validate: test Python access to dmeff parameters and derived quantities
- [ ] Phase 8: Build and integration testing
  - [ ] Clean build: make clean and make
  - [ ] Build Python wrapper: cd python and python setup.py install
  - [ ] Run all five test cases in class_dmeff_uptodate
  - [ ] Compare outputs against baseline reference (< 0.1% relative error)
  - [ ] Run vanilla CLASS test (N_dmeff equals 0) and verify no changes
  - [ ] Check for memory leaks with valgrind
  - [ ] Verify both synchronous and Newtonian gauges produce consistent results
- [ ] Documentation
  - [ ] Create test_dmeff/README.md explaining test suite
  - [ ] Document dmeff physics in comments where appropriate
  - [ ] Create example notebooks demonstrating dmeff usage


## Surprises & Discoveries

This section will document unexpected behaviors, bugs, optimizations, or insights discovered during implementation. As work proceeds, observations will be recorded here with evidence.


## Decision Log

This section records every design decision made while working on the plan.


## Outcomes & Retrospective

This section will summarize outcomes, gaps, and lessons learned at major milestones and at completion.


## Context and Orientation

The repository contains three directories at the top level. The directory class_dmeff contains CLASS version 2.9.4 with dmeff physics already implemented. The directory class_dmeff_uptodate contains pristine CLASS version 3.3.4 without dmeff, representing the target for this port. The directory class_public contains the current online vanilla CLASS for reference and must not be modified.

CLASS is a Boltzmann code written in C that computes cosmological observables. It is organized into modules that each handle a stage of the calculation. The background module computes expansion history and homogeneous densities. The thermodynamics module computes ionization history, optical depth, and visibility function. The perturbations module evolves linear density and velocity perturbations for all species. The transfer module extracts transfer functions from the perturbation evolution. The primordial module sets initial conditions from inflation. The fourier module (called nonlinear in older versions) computes nonlinear corrections. The harmonic module (called spectra in older versions) computes angular power spectra. Each module has a source file in source/ and a header in include/ that defines structures and function declarations.

The dmeff implementation adds a new cosmological fluid representing dark matter that interacts with baryons. The term dmeff stands for dark matter with effective interactions. Unlike cold dark matter (CDM) which is completely collisionless, dmeff undergoes momentum transfer with baryons, hydrogen, helium, or electrons depending on configuration. The cross section for this scattering can be velocity-dependent following a power law: sigma(v) equals sigma_zero times v to the power n, where sigma_zero is sigma_dmeff and n is npow_dmeff. The interaction causes heat exchange between dmeff and baryons, giving dmeff a separate temperature that evolves in time. In the perturbation equations, dmeff velocity perturbations couple to baryon velocity perturbations through a drag term proportional to the momentum exchange rate.

The implementation in class_dmeff touches approximately 630 locations across eight C modules. The core physics is in three modules. Background adds dmeff as a matter component with density rho_dmeff equals Omega0_dmeff times (H_zero/H) squared divided by a cubed, tracking this density plus derived quantities like temperature and relative velocity. Thermodynamics integrates a differential equation for the dmeff temperature backward from today to high redshift, computing the momentum exchange rate dkappa_dmeff and heat exchange rate dkappaT_dmeff at each step based on the target particle properties and the velocity-dependent cross section. Perturbations adds two dynamical variables delta_dmeff and theta_dmeff representing density and velocity divergence perturbations, with evolution equations that couple to the baryon perturbations through the drag term dkappa_dmeff times (theta_baryon minus theta_dmeff).

Between CLASS version 2.9.4 and version 3.3.4, the CLASS developers made architectural changes. The nonlinear module was renamed to fourier with structure prefix changed from nl to fo. The spectra module was renamed to harmonic with structure prefix changed from sp to hr. The header precision_macros.h was renamed to macros_precision.h. External dependencies were reorganized into an external/ subdirectory, and HyRec was updated from an older version to HyRec2020. A new distortions module was added for spectral distortion calculations. These changes mean that porting dmeff is not a simple file copy but requires understanding what was changed in class_dmeff and reapplying those changes to the corresponding locations in class_dmeff_uptodate, accounting for the renamed modules and updated file structure.


## Plan of Work

The work proceeds in phases that each deliver a testable increment. Phase 0 establishes the reference baseline by running the existing dmeff implementation and saving outputs. Phases 1 through 7 port dmeff module by module. Phase 8 integrates and validates everything together.

Phase 0 creates validation tests. Start by examining class_dmeff/ to understand what test parameter files already exist. Create five representative test cases covering different physics regimes: one with npow_dmeff equals negative four to test Coulomb-like interactions, one with npow_dmeff equals zero for velocity-independent scattering, one targeting electrons instead of baryons, one with partial replacement of CDM by dmeff to test coexistence, and one with multiple interaction terms to test the array handling. For each test case, create an ini file in class_dmeff/ if not already present, run the CLASS executable, and save the output files as reference data. Create a directory test_dmeff/ in class_dmeff_uptodate/ to hold the test ini files and reference outputs. Document the expected outputs for each quantity: background quantities like rho_dmeff at z equals 0 and z equals 1100, thermodynamics quantities like T_dmeff and z_dmeff_decoupling, perturbation outputs for several k values, and transfer function shapes.

Phase 1 ports the background module. Open include/background.h in class_dmeff_uptodate/ and add the dmeff parameter declarations by examining class_dmeff/include/background.h to see exactly what was added. The structure needs double Omega0_dmeff for the present-day density parameter, double m_dmeff for the particle mass in kilograms, int N_dmeff for the number of interaction terms, double pointer sigma_dmeff for an array of cross sections, double pointer npow_dmeff for an array of velocity power indices, double Vrel_dmeff for the initial relative bulk velocity, int has_dmeff as a boolean flag, and enum select_dmeff_target pointer dmeff_target for an array of target particle types. Add indices as integers: index_bg_rho_dmeff, index_bg_Tdmeff, index_bg_Vrel_dmeff, index_bg_dkappa_dmeff, index_bg_dkappaT_dmeff, and index_bg_cdmeff2 for the speed of sound squared.

Open source/background.c and find background_init. In class_dmeff/source/background.c, search for has_dmeff to see where arrays are allocated. Add memory allocation for sigma_dmeff, npow_dmeff, and dmeff_target using malloc with size N_dmeff times the appropriate type size. Increment the background quantity counter pba->bg_size to account for dmeff quantities if has_dmeff is true. Find where the background indices are assigned, which is typically a sequence of lines setting index_bg_rho_cdm, index_bg_rho_baryon, and so on, and insert the dmeff index assignments in that sequence.

In background_solve, locate where densities are computed for each species. This is typically in a loop over time or scale factor. Add computation of rho_dmeff using the formula rho_dmeff equals Omega0_dmeff times pvecback[index_bg_H] squared divided by a cubed, where a is the scale factor and H is extracted from the background vector. Include rho_dmeff in the total matter density calculation by adding it to the sum that defines rho_matter.

In background_functions, find where the background vector pvecback is filled with current values. Add assignments for each dmeff index: pvecback[index_bg_rho_dmeff] equals rho_dmeff computed above, and similarly for other quantities. At this stage, some quantities like Tdmeff come from thermodynamics interpolation, but since thermodynamics is not yet ported, set them to zero or skip them for now.

In background_free, add deallocation calls using free() for sigma_dmeff, npow_dmeff, and dmeff_target if has_dmeff is true.

In background_output_titles and background_output_data, add dmeff columns. These functions build strings or write columns for the background output file. Add entries for rho_dmeff, Tdmeff, and other dmeff quantities.

To validate Phase 1, compile with make clean followed by make. If compilation succeeds, run CLASS with a minimal dmeff test case setting Omega_dmeff but N_dmeff equals zero (no interactions, just density tracking). Check that the background output file contains a rho_dmeff column and that values match expectations from the reference. This tests that the background module correctly tracks the dmeff density without yet implementing temperature or interactions.

Phase 2 ports the thermodynamics module. This is the most complex phase because it involves integrating a differential equation for the dmeff temperature. Open include/thermodynamics.h in class_dmeff_uptodate/ and add indices: index_th_Tdmeff, index_th_dkappa_dmeff, index_th_ddkappa_dmeff for the derivative of the momentum exchange rate, index_th_dkappaT_dmeff for the heat exchange rate, and index_th_cdmeff2 for the speed of sound squared. Add double z_dmeff_decoupling to the thermodynamics structure as a derived parameter. Add the enum select_dmeff_target with values baryon, hydrogen, helium, electron (these are not quoted, they are C enum identifiers).

Open source/thermodynamics.c and locate the thermodynamics_init function. This is where the thermodynamics module sets up tables and performs integrations. Find where ionization fractions and optical depth are computed, which typically involves integrating an ODE backward in time. Add a call to a new function thermodynamics_dmeff_temperature that will integrate the dmeff temperature. Define this function earlier in the file.

The function thermodynamics_dmeff_temperature sets up an ODE integrator similar to how CLASS integrates the ionization fraction. The ODE to integrate is dT_dmeff over d tau equals negative two times a times H times T_dmeff plus two times rate_heat times (T_baryon minus T_dmeff), where tau is conformal time, a is the scale factor, H is the Hubble parameter, and rate_heat is the heat exchange rate. This equation describes how dmeff temperature cools adiabatically (first term) and heats or cools through collisions with baryons (second term).

To compute rate_heat, implement a function thermodynamics_dmeff_rate. This function loops over the N_dmeff interaction terms. For each term, it identifies the target particle (baryon, hydrogen, helium, or electron) and extracts the target mass and number density from the thermodynamics tables. It computes the thermal velocity squared as v_th squared equals k_B times T_target over m_target plus (Vrel/c) squared, where k_B is Boltzmann constant, T_target is the target temperature (baryon temperature for all targets except electron which uses electron temperature if tracked separately, but typically electron temperature equals baryon temperature in the tight-coupling regime), m_target is the target mass, and Vrel is the relative bulk velocity. The normalization constant cn equals two to the power ((npow+5)/2) times Gamma(3 + npow/2) divided by (3 times sqrt(pi)), where Gamma is the gamma function from the C math library. The momentum transfer rate is rate_mom equals a times rho_target times cn times sigma/(m_dmeff + m_target) times (v_th squared / c squared) to the power ((npow+1)/2). The heat transfer rate is rate_heat equals rate_mom times m_dmeff / (m_dmeff + m_target). Sum the contributions from all N_dmeff terms to get total dkappa_dmeff (sum of rate_mom) and total dkappaT_dmeff (sum of rate_heat).

The function thermodynamics_dmeff_derivs provides the derivatives for the ODE solver. It calls thermodynamics_dmeff_rate to get dkappaT_dmeff, then returns dy equals negative two times a times H times y plus two times dkappaT_dmeff times (T_baryon minus y), where y is the current T_dmeff value.

After integrating the ODE, compute the decoupling redshift z_dmeff_decoupling. This is defined as the redshift where dkappa_dmeff times tau_zero is less than ten, where tau_zero is the conformal time today. Scan the interpolation table to find where this condition is first met. Store the result in pth->z_dmeff_decoupling.

In thermodynamics_init, after calling thermodynamics_dmeff_temperature, increment pth->th_size to account for the dmeff quantities. Assign indices in the sequence where other thermodynamics indices are assigned.

Create the function thermodynamics_at_z or modify it if it exists. This function interpolates thermodynamics quantities at a given redshift. Add interpolation for Tdmeff, dkappa_dmeff, and related quantities from the tables built during integration. Compute the speed of sound squared as c_dmeff squared equals k_B over (m_dmeff times c squared) times (T_dmeff minus dT_dmeff over (3 times a times H)), where dT_dmeff is the time derivative of temperature. This derivative can be approximated from the interpolation table or computed directly from the ODE formula.

In thermodynamics output functions, add columns for Tdmeff, dkappa_dmeff, dkappaT_dmeff, and cdmeff2.

To validate Phase 2, compile and run a test case with interactions enabled (N_dmeff equals one, nonzero sigma_dmeff). Check that the thermodynamics output file contains Tdmeff and dkappa_dmeff columns. Compare Tdmeff at several redshifts against the reference baseline from class_dmeff. Verify that z_dmeff_decoupling is computed and matches the reference value. Check that dkappa_dmeff decreases at low redshift (late times) as the interaction rate drops due to cosmic expansion diluting the baryon density.

Phase 3 ports the perturbations module. Open include/perturbations.h and add indices index_pt_delta_dmeff and index_pt_theta_dmeff for the perturbation vector. Add source function indices index_tp_delta_dmeff and index_tp_theta_dmeff for the transfer function output. These indices are integers in the perturbations structure.

Open source/perturbations.c and find perturb_vector_init. This function allocates space for the perturbation vector y which holds all the perturbation variables (photon, baryon, CDM, neutrino perturbations, etc.). Increment the vector size to accommodate delta_dmeff and theta_dmeff. Assign indices in the sequence where other perturbation indices are assigned, typically after the CDM perturbations.

In perturb_initial_conditions, set the initial values for delta_dmeff and theta_dmeff. For adiabatic initial conditions in the radiation-dominated era, all matter perturbations start with the same amplitude. Find where delta_cdm and theta_cdm are initialized and add analogous initialization for delta_dmeff and theta_dmeff. In synchronous gauge, for super-horizon modes, delta_dmeff is set to a constant proportional to the curvature perturbation, and theta_dmeff starts at zero or a small value depending on the mode. Refer to class_dmeff/source/perturbations.c for the exact formulas, which depend on k (the wavenumber) and the initial curvature.

In perturb_einstein, include dmeff in the Einstein equations. The Einstein equations relate metric perturbations (h and eta in synchronous gauge, phi and psi in Newtonian gauge) to the energy-momentum tensor sourced by all species. Find where rho_cdm times delta_cdm and similar terms are summed for the density source. Add rho_dmeff times delta_dmeff to that sum. Similarly, for the velocity source (which is rho plus p times theta), add rho_dmeff times theta_dmeff since dmeff pressure is small but nonzero due to its temperature, specifically p_dmeff equals rho_dmeff times cdmeff2.

Implement gauge transformation functions. In CLASS, perturbations can be computed in synchronous gauge and then transformed to Newtonian gauge or vice versa. Find the gauge transformation block in perturbations.c. It typically has an if statement checking ppt->gauge. In the synchronous to Newtonian transformation, apply delta_dmeff transforms to delta_dmeff minus three times (a prime over a) times alpha and theta_dmeff transforms to theta_dmeff plus k squared times alpha, where alpha is the gauge transformation variable related to the time shift. For the inverse transformation, reverse the signs. Refer to class_dmeff/source/perturbations.c lines around 5567 and 8160 for the exact code.

In perturb_derivs, add the evolution equations for dmeff. This function computes dy/dtau for all perturbation variables. Define metric_continuity and metric_euler, which are gauge-dependent metric terms. In synchronous gauge, metric_continuity equals h_prime over two and metric_euler equals zero. In Newtonian gauge, metric_continuity equals negative three times phi_prime and metric_euler equals k squared times psi. Then set dy[index_pt_delta_dmeff] equals negative (y[index_pt_theta_dmeff] plus metric_continuity) and dy[index_pt_theta_dmeff] equals negative (a_prime over a) times y[index_pt_theta_dmeff] plus k squared times cdmeff2 times y[index_pt_delta_dmeff] plus metric_euler plus dkappa_dmeff times (theta_baryon minus y[index_pt_theta_dmeff]). The last term is the drag coupling to baryons.

In the source function section of perturbations.c, find where transfer functions are prepared. Add delta_dmeff and theta_dmeff to the sources output vector so the transfer module can extract them.

To validate Phase 3, compile and run a test case. After the run, check that the perturbations output file (if any) contains delta_dmeff and theta_dmeff for various k modes and times. More importantly, the transfer module will use these, so validation will be clearer after Phase 5. For now, check that the code runs without crashing and that the perturbation vector sizes are correct.

Phase 4 ports the input module. Open source/input.c and find the input_read_parameters function. This large function reads the ini file and parses parameters. Search for where Omega_cdm or omega_cdm is read to find the pattern. Add reading for Omega_dmeff, omega_dmeff, and f_dmeff. These three are mutually exclusive ways to specify the dmeff density. If Omega_dmeff is present, use it directly. If omega_dmeff is present (lowercase omega is the physical density parameter omega equals Omega times h squared), convert to Omega_dmeff by dividing by h squared. If f_dmeff is present (fraction of total dark matter), compute Omega_dmeff as f_dmeff times the sum of Omega_cdm and Omega_dmeff_total, but since Omega_dmeff_total is what we are solving for, this requires adjusting Omega_cdm to make room. Specifically, set Omega_dmeff equals f_dmeff times (Omega_cdm_input plus Omega_dmeff_input) and Omega_cdm equals (1 minus f_dmeff) times (Omega_cdm_input plus Omega_dmeff_input), where input values are what the user provided. The logic in class_dmeff handles this; replicate it exactly.

Read m_dmeff as a double in GeV and convert to kilograms by multiplying by the electron volt mass conversion factor (approximately 1.78e-27 kg per GeV, but CLASS has a constant for this defined in common.h).

Read N_dmeff as an integer. Then read the array parameters. CLASS input.c has utilities for reading arrays from strings. Read sigma_dmeff as an array of doubles with N_dmeff entries. Also support log10sigma_dmeff where the user specifies log base ten of the cross section; if this is used, compute sigma_dmeff[i] equals pow(10, log10sigma_dmeff[i]) and convert from cm squared to m squared by multiplying by 1e-4. Read npow_dmeff as an array of doubles with N_dmeff entries. Validate that each entry is greater than or equal to negative four; if not, return an error. Read dmeff_target as an array of strings with N_dmeff entries. For each string, compare to "baryon", "hydrogen", "helium", "electron" and set the enum value accordingly. If the string does not match any of these, return an error message indicating invalid target.

Read Vrel_dmeff as a double in km/s and convert to m/s by multiplying by 1000.

After reading all parameters, set has_dmeff. If Omega_dmeff is greater than zero and N_dmeff is greater than zero, set has_dmeff to _TRUE_. Otherwise set to _FALSE_.

In the section of input.c where parameters are echoed or validated, add checks for array length consistency. If N_dmeff is greater than zero but the sigma_dmeff array length does not match, return an error. Similarly for npow_dmeff and dmeff_target.

Update class_dmeff_uptodate/explanatory.ini by copying the dmeff parameter documentation from class_dmeff/explanatory.ini. This file serves as a reference for users. The documentation should appear in the section where other matter component parameters are documented, typically after omega_cdm. Include comments explaining that Omega_dmeff, omega_dmeff, and f_dmeff are mutually exclusive, that N_dmeff controls the number of interaction terms, that sigma_dmeff can be specified directly or as log10sigma_dmeff, that npow_dmeff must be at least negative four, that dmeff_target accepts baryon, hydrogen, helium, or electron, and that Vrel_dmeff is the initial relative bulk velocity in km/s.

To validate Phase 4, compile and test that the input parser reads dmeff parameters correctly. Create an ini file with Omega_dmeff equals 0.12, run CLASS, and check that the background module receives the correct value. Test error handling by providing mismatched array lengths or an invalid target string and verify that CLASS returns a clear error message. Test all three density specification methods (Omega, omega, f) and confirm they produce the expected Omega_dmeff.

Phase 5 ports the transfer module. Open include/transfer.h and add indices index_tr_delta_dmeff and index_tr_theta_dmeff as integers in the transfer structure. Open source/transfer.c and find transfer_init. This function counts how many transfer types to compute. If has_dmeff is true, increment the transfer type counter by two (one for delta, one for theta). Assign indices in the section where index_tr_delta_cdm and similar are assigned.

In transfer_functions, locate where transfer functions are extracted from the perturbation sources. For each k mode, after the perturbation module has computed sources, the transfer module reads delta_dmeff and theta_dmeff from the sources vector and stores them in the transfer output array. Add these extractions parallel to how delta_cdm is handled.

In the transfer output section, add dmeff columns to the transfer function file. The output file typically has columns for k, T_delta_cdm, T_theta_baryon, and so on. Add T_delta_dmeff and T_theta_dmeff.

To validate Phase 5, compile and run a test case. Check that the transfer function output file (usually named something like test_tk.dat or transfer_functions.dat) contains dmeff columns. Compare the transfer functions T_delta_dmeff(k) against the reference baseline. The shape should be similar to the CDM transfer function but with suppression at small scales if the interaction is strong, reflecting the drag from baryons.

Phase 6 ports the warning about nonlinear corrections. In class_dmeff/source/nonlinear.c, search for "dmeff" to find the warning code. It likely appears in nonlinear_init and issues a warning if has_dmeff is true and the user has requested Halofit or HMcode. The warning states that nonlinear corrections may not be reliable in the presence of DM-baryon interactions because Halofit and HMcode are calibrated on standard CDM simulations.

In class_dmeff_uptodate/source/fourier.c (the renamed module), locate the corresponding init function, probably fourier_init. Find a suitable place to insert the warning, such as after reading the method parameter. If has_dmeff is true and the method is not linear, print a warning message using the CLASS warning or message system (typically fprintf(stderr,...) or a CLASS-specific warning macro). The exact text can be: "Warning: using dmeff dark matter interactions with nonlinear power spectrum method; Halofit/HMcode may be inaccurate as they are calibrated for collisionless CDM."

To validate Phase 6, compile and run a test case with output equals mPk (matter power spectrum) which triggers the nonlinear module. Verify that the warning message appears in the output or log.

Phase 7 ports the Python wrapper. Open python/cclassy.pxd in class_dmeff_uptodate. This file declares the C structures for Cython. Find the background structure declaration (it will have cdef extern from "background.h") and add "double Omega0_dmeff" as a new field. Find the thermodynamics structure and add "double z_dmeff_decoupling".

Open python/classy.pyx. Find where methods like Omega0_cdm() are defined. Add a similar method:

    def Omega0_dmeff(self):
        return self.ba.Omega0_dmeff

This allows Python users to call cosmo.Omega0_dmeff() and retrieve the dmeff density parameter.

Find the get_current_derived_parameters method. This method has a large if-elif chain for different derived parameter names. Add an elif block:

    elif name == 'z_dmeff_decoupling':
        value = self.th.z_dmeff_decoupling

This allows Python users to call cosmo.get_current_derived_parameters(['z_dmeff_decoupling']) and retrieve the decoupling redshift.

To validate Phase 7, compile the Python wrapper. In the class_dmeff_uptodate directory, run "cd python" and "python setup.py install" (or "pip install ." if using a virtual environment). Then start a Python interpreter, import classy, create a Classy instance, set parameters including dmeff parameters, call compute(), and then retrieve Omega0_dmeff() and get_current_derived_parameters(['z_dmeff_decoupling']). Verify that the returned values match expectations.

Phase 8 performs integration testing. Run "make clean" followed by "make" in the class_dmeff_uptodate directory to ensure a fresh build. Check for compiler warnings and errors. If successful, the executable "class" should exist in the main directory.

Run the five test cases created in Phase 0. For each test case, execute "./class test_dmeff/test_NAME.ini" where NAME is the test case name. CLASS will produce output files in the output directory specified in the ini file. Compare each output against the reference baseline. Use a script or manual inspection to compute relative errors for key quantities. For background, check rho_dmeff at z equals 0 and z equals 1100. For thermodynamics, check Tdmeff at several redshifts and z_dmeff_decoupling. For transfer functions, check T_delta_dmeff at several k values spanning two decades. All relative errors should be less than 0.1 percent.

Run a vanilla CLASS test. Create an ini file with no dmeff parameters (or Omega_dmeff equals 0 and N_dmeff equals 0). Run CLASS and compare the output to the same run using class_public or to a reference vanilla CLASS output. The outputs should match exactly or to machine precision, confirming that dmeff additions do not affect normal CLASS operation when disabled.

Check for memory leaks. If valgrind is available, run "valgrind --leak-check=full ./class test_dmeff/test_coulomb.ini". Inspect the output for any "definitely lost" or "possibly lost" memory blocks. If leaks are found, review the memory allocation and deallocation in background_free and other cleanup functions.

Verify gauge consistency. Run a test case in synchronous gauge (the default) and then run the same case in Newtonian gauge by adding "gauge = newtonian" to the ini file. Compare the final outputs (e.g., CMB power spectra or matter power spectra). These should match to within numerical precision, confirming that gauge transformation logic is correct.

At the completion of Phase 8, the port is functionally complete. All dmeff physics is present and validated. Write the Outcomes & Retrospective section summarizing what was achieved, any discrepancies from the plan, and lessons learned.


## Concrete Steps

The following commands assume the working directory is the repository root /Users/veragluscevic/research/repositories/class-updates-via-cursor.

Phase 0 Preparation:

    cd class_dmeff
    ./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini

(This runs a baseline cosmology. If it succeeds, output files appear in output/ directory. If there are existing dmeff test ini files, use those. Otherwise, create test cases by editing ini files.)

    mkdir -p ../class_dmeff_uptodate/test_dmeff
    # Copy or create test ini files in test_dmeff/
    # Run each test case and save outputs
    ./class test_dmeff/test_coulomb.ini
    cp output/* ../class_dmeff_uptodate/test_dmeff/reference_coulomb/

(Repeat for each test case, saving reference outputs in subdirectories.)

Phases 1-7 Editing:

For each phase, edit the specified files in class_dmeff_uptodate. After editing, test compilation:

    cd class_dmeff_uptodate
    make clean
    make

If compilation succeeds, run a simple test to confirm the module works:

    ./class test_dmeff/test_simple.ini

If the test fails or crashes, debug by examining error messages and comparing code against class_dmeff.

Phase 8 Integration Testing:

    cd class_dmeff_uptodate
    make clean
    make
    ./class test_dmeff/test_coulomb.ini
    ./class test_dmeff/test_constant.ini
    ./class test_dmeff/test_electron.ini
    ./class test_dmeff/test_mixed.ini
    ./class test_dmeff/test_multi.ini

For each test, compare output files to references. Use a script to compute relative errors. Example for background:

    python -c "
    import numpy as np
    data = np.loadtxt('output/test_coulomb_background.dat')
    ref = np.loadtxt('test_dmeff/reference_coulomb/test_coulomb_background.dat')
    reldiff = np.abs(data - ref) / (np.abs(ref) + 1e-30)
    print('Max relative error:', np.max(reldiff))
    "

(Expected output: Max relative error less than 0.001.)

Python wrapper testing:

    cd python
    python setup.py install
    python -c "
    from classy import Class
    cosmo = Class()
    cosmo.set({'Omega_dmeff': 0.12, 'm_dmeff': 1.0, 'N_dmeff': 1,
               'sigma_dmeff': 1e-30, 'npow_dmeff': -4, 'dmeff_target': 'hydrogen',
               'output': 'tCl'})
    cosmo.compute()
    print('Omega0_dmeff:', cosmo.Omega0_dmeff())
    print('z_dmeff_decoupling:', cosmo.get_current_derived_parameters(['z_dmeff_decoupling']))
    "

(Expected output: Omega0_dmeff approximately 0.12, z_dmeff_decoupling some value greater than 1000.)

Memory leak check:

    valgrind --leak-check=full ./class test_dmeff/test_coulomb.ini 2>&1 | grep "definitely lost"

(Expected output: zero bytes definitely lost.)


## Validation and Acceptance

After completing Phase 8, the implementation is accepted if the following behaviors are observed.

Run a dmeff test case by executing "./class test_dmeff/test_coulomb.ini" from the class_dmeff_uptodate directory. The command should complete without error, producing output files in the specified output directory. Open the background output file and verify that it contains columns for rho_dmeff and Tdmeff. Compute the relative difference between rho_dmeff at z equals 1100 in this output and the corresponding value in the reference output from class_dmeff. The relative difference must be less than 0.01 percent.

Open the thermodynamics output file and verify that it contains columns for Tdmeff and dkappa_dmeff. Extract Tdmeff at redshift z equals 10 from both the new output and the reference. The relative difference must be less than 0.1 percent. Check that z_dmeff_decoupling is written to the log or to a derived parameters file. The value should be a redshift greater than zero and less than the maximum redshift in the computation (typically 1e4 for standard runs). Compare this value to the reference; it should agree to within 1 percent.

Open the transfer function output file and find the columns for T_delta_dmeff and T_theta_dmeff. Plot T_delta_dmeff(k) versus k and compare visually to the reference transfer function. The curves should overlap. Compute the relative difference at k equals 0.1 per Mpc and at k equals 1.0 per Mpc. Each relative difference must be less than 0.1 percent.

Run a vanilla CLASS test by executing "./class test_vanilla.ini" where test_vanilla.ini has standard cosmological parameters and no dmeff parameters. Compare the output CMB C_l values to those from class_public using the same parameters. The values should match to within 0.01 percent, demonstrating that dmeff additions do not alter vanilla CLASS behavior.

Run the Python test script shown in Concrete Steps. The script should print Omega0_dmeff approximately 0.12 and z_dmeff_decoupling a positive value. The script should not raise exceptions. This confirms the Python wrapper exposes dmeff parameters.

Run valgrind with a test case. The valgrind output summary should report zero bytes definitely lost and zero errors. This confirms no memory leaks.

Run a test case with gauge equals newtonian in the ini file. Compare the output CMB TT power spectrum to the same case run in synchronous gauge. The two spectra should match to within 0.01 percent, confirming gauge transformation correctness.

If all the above conditions hold, the dmeff port is successful and the feature is complete.


## Idempotence and Recovery

The editing steps are idempotent in the sense that if Phase N is partially completed, it can be resumed by continuing the edits without undoing prior work. Compilation will fail until all required edits in a phase are complete, but failed compilation does not corrupt the source tree. Running "make clean" removes object files and the executable, allowing a fresh build.

If a test case produces incorrect outputs, compare the code in class_dmeff_uptodate to class_dmeff to identify missing or incorrect changes. Use "diff -u class_dmeff/source/background.c class_dmeff_uptodate/source/background.c" to see differences. If a module is incorrectly ported, re-examine the source and re-edit to match the reference.

If valgrind reports memory leaks, use valgrind's detailed output to identify which allocation is not freed. Add the corresponding free() call in the appropriate cleanup function.

If the Python wrapper fails to build, check that the Cython declarations match the C structure definitions. Ensure that field names and types are exact. Python import errors often indicate ABI mismatches; rebuild CLASS with the same compiler flags used for the Python extension.

The repository state can be restored at any time by reverting class_dmeff_uptodate to its initial state (pristine CLASS v3.3.4). The class_dmeff and class_public directories are never modified and serve as stable references.


## Artifacts and Notes

This section will be populated with transcripts, diffs, or snippets as implementation proceeds. Examples of what will be recorded: compiler output showing successful build, sample output from running a test case, valgrind summary, Python test output.


## Interfaces and Dependencies

The implementation uses the existing CLASS module structure and does not introduce new external dependencies. The key interfaces are the structures defined in the header files.

In include/background.h, the background structure must contain:

    double Omega0_dmeff;
    double m_dmeff;
    int N_dmeff;
    double * sigma_dmeff;
    double * npow_dmeff;
    double Vrel_dmeff;
    int has_dmeff;
    enum select_dmeff_target * dmeff_target;
    int index_bg_rho_dmeff;
    int index_bg_Tdmeff;
    int index_bg_Vrel_dmeff;
    int index_bg_dkappa_dmeff;
    int index_bg_dkappaT_dmeff;
    int index_bg_cdmeff2;

The enum select_dmeff_target is defined in include/thermodynamics.h:

    enum select_dmeff_target {baryon, hydrogen, helium, electron};

In include/thermodynamics.h, the thermodynamics structure must contain:

    int index_th_Tdmeff;
    int index_th_dkappa_dmeff;
    int index_th_ddkappa_dmeff;
    int index_th_dkappaT_dmeff;
    int index_th_cdmeff2;
    double z_dmeff_decoupling;

In include/perturbations.h, the perturbs structure must contain:

    int index_pt_delta_dmeff;
    int index_pt_theta_dmeff;

The transfer structure in include/transfer.h must contain:

    int index_tr_delta_dmeff;
    int index_tr_theta_dmeff;

New functions to define in source/thermodynamics.c:

    int thermodynamics_dmeff_rate(
        struct background * pba,
        struct thermodynamics * pth,
        double z,
        double * pvecback,
        double * pvecthermo,
        double * dkappa_dmeff,
        double * dkappaT_dmeff,
        double * ddkappa_dmeff
    );

    int thermodynamics_dmeff_temperature(
        struct background * pba,
        struct thermodynamics * pth
    );

    int thermodynamics_dmeff_derivs(
        double tau,
        double * y,
        double * dy,
        void * params,
        ErrorMsg error_message
    );

These functions are called from thermodynamics_init and implement the temperature evolution integration.

The physics constants needed are defined in include/common.h or include/constant.h (depending on CLASS version naming). Key constants: _c_ for speed of light in m/s, _k_B_ for Boltzmann constant in J/K, _m_H_ for hydrogen mass, _m_e_ for electron mass. The conversion from GeV to kg is approximately 1.78266192e-27, but CLASS may define this as a constant.

The ODE integrator used by CLASS is typically evolver_rk or evolver_ndf from the evolver utility. The function thermodynamics_dmeff_temperature will call the same integrator used for other thermodynamics integrations, passing thermodynamics_dmeff_derivs as the derivative function.
