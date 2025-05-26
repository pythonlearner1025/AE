import math

# ==============================================================================
# --- SCRIPT FOR ESTIMATING NOISE EQUIVALENT CURRENT (NEC) FOR -----------
# --- ACOUSTOELECTRIC IMAGING (AEI) / TRANSCRANIAL ACOUSTOELECTRIC -------
# --- BRAIN IMAGING (tABI) -------------------------------------------------
# ==============================================================================
#
# This script calculates the theoretical NEC based on simplified assumptions
# for thermal noise and AE signal generation.
#
# NEC is the neuronal source current that produces an AE signal equal to
# the system's noise level. A lower NEC indicates better sensitivity.
#
# Formulae used:
# 1. Sensitivity (S) = P0 * K * rho0 * J_L_norm * (V_focus / A_source)
# 2. Thermal Noise Voltage (V_thermal) = sqrt(4 * k_B * T * R_det * delta_f_eff)
# 3. Noise Equivalent Current (NEC) = V_thermal / S
#
# DISCLAIMER: This is a simplified model. The "Lead Field Current Density
# Magnitude" (J_L_norm) is a particularly complex parameter that is highly
# dependent on specific electrode geometry, tissue conductivities, and source
# location; the value used here is a rough order-of-magnitude estimate for
# transcranial sensing of a relatively deep source. Actual experimental NEC
# may vary due to other noise sources, model simplifications, and specific
# experimental conditions.
#
# ==============================================================================
# --- 1. DEFINE CONSTANTS AND PARAMETERS -------------------------------------
# ==============================================================================
# You can modify the values in this section to explore their impact.

# --- Physical Constants ---
K_B = 1.380649e-23  # Boltzmann's constant (J/K)

# --- Ultrasound and Material Properties ---
P0 = 0.67e6          # Peak Ultrasound Pressure at focus (Pascals, Pa)
                    # Example: 1 MPa, typical for biologically safe intensities.

K_AE = 1.0e-9       # Acoustoelectric Interaction Constant (1/Pa or Pa^-1)
                    # Example: 0.1% per MPa = 1e-3 / 1e6 Pa = 1e-9 Pa^-1.

RHO0 = 3.0          # Average Tissue Resistivity (Ohm*m)
                    # Example: ~3 Ohm*m for brain tissue/saline gel.

V_FOCUS = 8.0e-9    # Ultrasound Focal Volume (m^3)
                    # Example: (2mm)^3 = (2e-3 m)^3 = 8e-9 m^3.

A_SOURCE = 4.0e-6   # Effective Cross-sectional Area of the Neuronal Current Source (m^2)
                    # Example: (2mm)^2 = (2e-3 m)^2 = 4e-6 m^2.
                    # This assumes the source area is comparable to the US beam width.

# --- Lead Field Parameter ---
# This parameter represents the current density magnitude at the source location
# if 1 Ampere of current were injected into the recording electrode pair.
# It acts as a transfer coefficient for how well the electrodes "sense"
# the current at the focal spot. This is a highly simplified representation.
# Units: effectively m^-2 (A/m^2 at source per A injected).
J_L_NORM = 100.0    # Normalized Lead Field Current Density Magnitude (m^-2)
                    # This is a rough estimate for transcranial sensing of a deep source.
                    # Higher values mean better coupling. ECoG would have much higher values.

# --- Detection System Noise Parameters ---
T_KELVIN = 310.0    # Absolute Temperature (Kelvin, K)
                    # Example: 310 K is ~37 degrees Celsius (body temperature).

R_DET = 1000.0      # Effective Detection Resistance (Ohms)
                    # Includes electrode-skin interface resistance, tissue path,
                    # and input impedance of the amplifier if it contributes.
                    # Example: 1 kOhm for scalp EEG system.

DELTA_F_EFF = 250.0 # Effective Noise Bandwidth (Hertz, Hz)
                    # This is the bandwidth of the *final processed signal* after
                    # demodulation and filtering, not the raw RF bandwidth.
                    # Example: 3 kHz, allowing for detection of relatively fast
                    # neuronal current changes. A smaller value (e.g., 500 Hz)
                    # would reduce noise but limit temporal resolution.

# ==============================================================================
# --- 2. CALCULATIONS --------------------------------------------------------
# ==============================================================================

def calculate_effective_focal_length(v_focus, a_source):
    """
    Calculates an effective length scale related to the focal zone geometry.
    L_eff_focus = V_focus / A_source
    """
    if a_source <= 0:
        raise ValueError("A_source (Source Cross-sectional Area) must be positive.")
    return v_focus / a_source

def calculate_sensitivity(p0, k_ae, rho0, j_l_norm, l_eff_focus):
    """
    Calculates the system's sensitivity (S) to current.
    S = P0 * K_AE * rho0 * J_L_norm * L_eff_focus
    Units: Ohms
    """
    sensitivity = p0 * k_ae * rho0 * j_l_norm * l_eff_focus
    return sensitivity

def calculate_thermal_noise_voltage(k_b, t_kelvin, r_det, delta_f_eff):
    """
    Calculates the RMS thermal noise voltage (V_thermal).
    V_thermal = sqrt(4 * k_B * T * R_det * delta_f_eff)
    Units: Volts
    """
    if r_det < 0:
        raise ValueError("R_det (Effective Detection Resistance) cannot be negative.")
    if delta_f_eff <= 0:
        raise ValueError("DELTA_F_EFF (Effective Noise Bandwidth) must be positive.")
    thermal_noise_voltage = math.sqrt(4 * k_b * t_kelvin * r_det * delta_f_eff)
    return thermal_noise_voltage

def calculate_nec(v_thermal, sensitivity):
    """
    Calculates the Noise Equivalent Current (NEC).
    NEC = V_thermal / S
    Units: Amperes
    """
    if sensitivity == 0:
        # Avoid division by zero; indicates no sensitivity.
        return float('inf')
    nec = v_thermal / sensitivity
    return nec

# --- Perform the calculations ---
l_eff_focus = calculate_effective_focal_length(V_FOCUS, A_SOURCE)
sensitivity_S = calculate_sensitivity(P0, K_AE, RHO0, J_L_NORM, l_eff_focus)
v_thermal_rms = calculate_thermal_noise_voltage(K_B, T_KELVIN, R_DET, DELTA_F_EFF)
nec_A = calculate_nec(v_thermal_rms, sensitivity_S)

# ==============================================================================
# --- 3. DISPLAY RESULTS -----------------------------------------------------
# ==============================================================================
print("===================================================")
print("--- Acoustoelectric Imaging NEC Model Results ---")
print("===================================================")
print("\nInput Parameters:")
print(f"  Peak Ultrasound Pressure (P0):         {P0:.2e} Pa")
print(f"  AE Interaction Constant (K_AE):        {K_AE:.2e} Pa^-1")
print(f"  Tissue Resistivity (rho0):             {RHO0:.2f} Ohm*m")
print(f"  Focal Volume (V_focus):                {V_FOCUS:.2e} m^3")
print(f"  Source Area (A_source):                {A_SOURCE:.2e} m^2")
print(f"  Normalized Lead Field (J_L_norm):      {J_L_NORM:.2f} m^-2 (effective)")
print(f"  Temperature (T):                       {T_KELVIN:.2f} K")
print(f"  Detection Resistance (R_det):          {R_DET:.2f} Ohms")
print(f"  Effective Noise Bandwidth (delta_f_eff): {DELTA_F_EFF:.2f} Hz")
print("---------------------------------------------------")
print("Calculated Intermediate Values:")
print(f"  Effective Focal Length (L_eff_focus):  {l_eff_focus:.2e} m")
print(f"  System Sensitivity (S):                {sensitivity_S * 1e3:.4f} mOhms") # Display in mOhms
print(f"  RMS Thermal Noise Voltage (V_thermal): {v_thermal_rms * 1e6:.4f} microVolts") # Display in microVolts
print("---------------------------------------------------")
print("Estimated NEC:")
if nec_A == float('inf'):
    print("  Noise Equivalent Current (NEC):        Infinite (System has zero sensitivity)")
else:
    print(f"  Noise Equivalent Current (NEC):        {nec_A * 1e6:.2f} microAmperes (uA)") # Display in microAmperes
    print(f"                                         {nec_A * 1e3:.4f} milliAmperes (mA)")  # Display in milliAmperes
print("===================================================")

