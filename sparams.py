import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math # Import math for checking isnan

# --- Core Calculation Logic (adapted from Flask version) ---

def calculate_poles_and_params(omega_n, zeta):
    """Calculates pole locations and related parameters."""
    # Ensure inputs are valid numbers, default to safe values if not
    if not isinstance(omega_n, (int, float)) or omega_n <= 0:
        omega_n = 1.0 # Default safe value
    if not isinstance(zeta, (int, float)) or math.isnan(zeta):
         zeta = 1.0 # Default safe value


    poles = []
    sigma = zeta * omega_n
    omega_d = 0
    theta = 0
    damping_type = ""
    error_msg = None

    # Basic input validation handled by streamlit widgets mostly, but double check omega_n
    if omega_n <= 0:
         error_msg = "ωn (Natural Frequency) must be positive."
         # Use default values to prevent calculation errors
         omega_n = 1.0
         zeta = 1.0
         sigma = zeta * omega_n


    # Determine damping type and calculate poles/parameters
    try:
        if zeta < 0:
            damping_type = "Unstable"
            # Formula remains valid mathematically, sqrt of negative gives imaginary
            omega_d_unstable = omega_n * np.sqrt(1 - zeta**2 + 0j) # Add 0j for complex sqrt
            poles = [-sigma + 1j * omega_d_unstable.imag, -sigma - 1j * omega_d_unstable.imag]
            sigma = -abs(sigma) # Make sigma explicitly negative for unstable label
            omega_d = omega_d_unstable.imag
            # Theta calculation is less standard, angle from positive real axis
            if poles[0].imag >= 0 and poles[0].real != 0:
                 theta = np.degrees(np.arctan2(poles[0].imag, poles[0].real))
            elif poles[0].real != 0:
                 theta = np.degrees(np.arctan2(poles[1].imag, poles[1].real))
            else: # Purely imaginary unstable poles (not physically typical for this model)
                 theta = 90 if poles[0].imag > 0 else -90

        elif zeta == 0:
            damping_type = "Undamped (ζ = 0)"
            poles = [1j * omega_n, -1j * omega_n]
            sigma = 0
            omega_d = omega_n
            theta = 90 # Angle from negative real axis to upper pole

        elif 0 < zeta < 1:
            damping_type = f"Underdamped (0 < {zeta:.2f} < 1)"
            omega_d = omega_n * np.sqrt(1 - zeta**2)
            poles = [-sigma + 1j * omega_d, -sigma - 1j * omega_d]
            theta = np.degrees(np.arccos(zeta)) # Angle from negative real axis

        elif zeta == 1:
            damping_type = "Critically Damped (ζ = 1)"
            poles = [-omega_n, -omega_n] # Two real, equal poles
            sigma = omega_n
            omega_d = 0
            theta = 0 # Angle from negative real axis

        else: # zeta > 1
            damping_type = f"Overdamped (ζ = {zeta:.2f} > 1)"
            # Use np.sqrt with complex numbers enabled in case zeta^2-1 is tiny negative due to float precision
            sqrt_term = np.sqrt(zeta**2 - 1 + 0j)
            pole1 = -zeta * omega_n + omega_n * sqrt_term.real
            pole2 = -zeta * omega_n - omega_n * sqrt_term.real
            poles = [pole1, pole2] # Two distinct real poles
            sigma = zeta * omega_n # Still represents the ζωn term magnitude
            omega_d = 0
            theta = 0 # Angle from negative real axis
    except Exception as e:
        error_msg = f"Calculation Error: {e}. Please check inputs."
        # Reset to defaults on error
        omega_n, zeta = 1.0, 1.0
        poles, sigma, omega_d, theta, damping_type = calculate_poles_and_params(omega_n, zeta)


    return poles, sigma, omega_d, theta, damping_type, error_msg, omega_n, zeta

# --- Plotting Function ---

def generate_s_plane_plot(omega_n, zeta, poles, sigma, omega_d, theta):
    """Generates the S-plane plot using Matplotlib."""
    fig, ax = plt.subplots(figsize=(7, 6))

    # Plot poles
    pole_reals = [p.real for p in poles]
    pole_imags = [p.imag for p in poles]
    ax.plot(pole_reals, pole_imags, 'rx', markersize=10, markeredgewidth=2, label='Poles')

    # Determine plot limits dynamically
    max_real_abs = max(0.5, np.max(np.abs(pole_reals))) if pole_reals else 0.5
    max_imag_abs = max(0.5, np.max(np.abs(pole_imags))) if pole_imags else 0.5
    max_val = max(max_real_abs, max_imag_abs, omega_n if omega_n else 1.0) * 1.4 # Extend slightly beyond max magnitude

    # Adjust x-limits based on stability
    min_xlim = -max_val
    max_xlim = max_val * 0.2 # Default focus on LHP

    if any(p.real > 1e-9 for p in poles): # Check if any pole is significantly in RHP
        max_xlim = max_val
        if not any(p.real < -1e-9 for p in poles): # If ALL poles are in RHP or on axis
             min_xlim = -max_val * 0.2

    ax.set_xlim(min_xlim, max_xlim)
    ax.set_ylim(-max_val, max_val)


    # Draw axes
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('Re(s) (σ)')
    ax.set_ylabel('Im(s) (jω)')
    ax.set_title(f'S-Plane Plot (ωn={omega_n:.2f}, ζ={zeta:.2f})')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_aspect('equal', adjustable='box') # Crucial for correct angle visualization

    # --- Annotations ---
    # Get the primary pole for annotation (upper half plane or first real pole)
    primary_pole = None
    if poles:
        if any(p.imag > 1e-9 for p in poles):
             primary_pole = poles[np.argmax([p.imag for p in poles])]
        else:
             primary_pole = poles[np.argmax([p.real for p in poles])] # Rightmost pole if real

    # Default annotation positions
    text_offset_x = max_val * 0.03
    text_offset_y = max_val * 0.03

    if primary_pole is not None:
        # ωn - distance from origin (for stable/critically damped)
        if zeta >= 0 and zeta <= 1:
             ax.plot([0, primary_pole.real], [0, primary_pole.imag], 'b--', linewidth=1, alpha=0.7)
             # Position wn label along the line
             label_pos_x = primary_pole.real / 2 + text_offset_y * (-primary_pole.imag / omega_n if omega_n else 0)
             label_pos_y = primary_pole.imag / 2 + text_offset_x * (primary_pole.real / omega_n if omega_n else 0)
             ax.text(label_pos_x, label_pos_y, f'ωn={omega_n:.2f}', color='blue', ha='center', va='bottom')


        # σ - real part distance/value
        if abs(primary_pole.real) > 1e-9:
            ax.plot([primary_pole.real, primary_pole.real], [0, primary_pole.imag], 'g:', linewidth=1, alpha=0.7)
            label_sig = f'σ={abs(primary_pole.real):.2f}' if zeta >=0 else f'σ={primary_pole.real:.2f} (unstable)'
            ax.text(primary_pole.real, text_offset_y, label_sig , color='green', ha='center', va='bottom')

        # ωd - imaginary part
        if abs(primary_pole.imag) > 1e-9:
             ax.plot([0, primary_pole.real], [primary_pole.imag, primary_pole.imag], 'm:', linewidth=1, alpha=0.7)
             ax.text(text_offset_x if primary_pole.real>0 else -text_offset_x, primary_pole.imag, f'ωd={omega_d:.2f}', color='purple', va='center', ha='left' if primary_pole.real<=0 else 'right')
        elif abs(zeta) >= 1: # Overdamped/Critically Damped
             ax.text(primary_pole.real, -text_offset_y*2, f'ωd = 0', color='purple', ha='center', va='top')


        # θ - angle arc (Handle different cases carefully)
        if 0 < zeta < 1: # Underdamped
            # Angle from negative real axis
            angle_rad_start = np.pi
            angle_rad_end = np.arctan2(primary_pole.imag, primary_pole.real)
            arc_radius = omega_n * 0.3 # Adjust radius as needed
            arc = plt.matplotlib.patches.Arc((0,0), arc_radius*2, arc_radius*2, angle=0,
                                            theta1=np.degrees(angle_rad_end), theta2=180, # Arc from pole vector to neg real axis
                                            color='orange', linewidth=1.5)
            ax.add_patch(arc)
            # Label position requires calculating midpoint angle
            mid_angle_rad = (angle_rad_start + angle_rad_end) / 2
            label_radius = arc_radius * 1.15
            ax.text(label_radius * np.cos(mid_angle_rad), label_radius * np.sin(mid_angle_rad), f'θ={theta:.1f}°', color='orange', ha='center', va='center')
        elif zeta == 1: # Critically Damped
             ax.text(primary_pole.real * 0.8, text_offset_y*1.5, f'θ = 0°', color='orange', ha='center')
             ax.text(primary_pole.real, text_offset_y, '(x2)', color='red', ha='left', va='bottom') # Indicate double pole
        elif zeta == 0: # Undamped
             ax.text(-text_offset_x*1.5, primary_pole.imag * 0.8, f'θ = 90°', color='orange', ha='right')
        elif zeta > 1: # Overdamped
             ax.text(poles[0].real, -text_offset_y*2, f'θ = 0°', color='orange', ha='left', va='top')
             ax.text(poles[1].real, text_offset_y, f'θ = 0°', color='orange', ha='right', va='bottom')
             ax.text(poles[0].real, text_offset_y*1.5, f'Pole 1', color='red', ha='center', va='bottom', fontsize=8)
             ax.text(poles[1].real, text_offset_y*1.5, f'Pole 2', color='red', ha='center', va='bottom', fontsize=8)
        # Unstable theta annotation is less standard, often omitted or angle from positive real shown


    # Add legend if labels were added
    handles, labels = ax.get_legend_handles_labels()
    if handles:
      ax.legend(loc='upper right', fontsize='small')

    return fig # Return the figure object

# --- Streamlit App Layout ---

st.set_page_config(layout="wide") # Use wider layout

st.title("Interactive Second-Order System S-Plane Visualizer")
st.markdown("""
    Adjust the **Undamped Natural Frequency (ωn)** and **Damping Ratio (ζ)**
    using the sliders or input boxes below to see how the system's poles
    move in the S-Plane and how the key parameters (σ, ωd, θ) change.
""")

# --- Input Controls ---
col1, col2 = st.columns([1, 2]) # Input column smaller than plot column

with col1:
    st.subheader("System Parameters")
    # Use session state to remember last valid inputs
    if 'omega_n' not in st.session_state:
        st.session_state.omega_n = 10.0
    if 'zeta' not in st.session_state:
        st.session_state.zeta = 1.0 # Start critically damped

    # Sliders and Number Inputs for flexibility
    omega_n_input = st.number_input("ωn (Natural Frequency, rad/s)", min_value=0.01, value=st.session_state.omega_n, step=0.5, format="%.2f", key="omega_n_num")
    st.session_state.omega_n = omega_n_input # Update state

    zeta_input = st.number_input("ζ (Damping Ratio)", min_value=-2.0, max_value=5.0, value=st.session_state.zeta, step=0.05, format="%.2f", key="zeta_num")
    st.session_state.zeta = zeta_input # Update state

    st.markdown("---") # Separator
    st.subheader("Calculated Values")

    # --- Calculations ---
    poles, sigma, omega_d, theta, damping_type, error_msg, omega_n_used, zeta_used = calculate_poles_and_params(
        st.session_state.omega_n,
        st.session_state.zeta
    )

    if error_msg:
        st.error(error_msg)

    # Display calculated values
    st.metric(label="Damping Type", value=damping_type)

    st.write("**Poles (s):**")
    pole_strs = []
    for p in poles:
        if abs(p.imag) < 1e-9:
            pole_strs.append(f"{p.real:.3f}")
        else:
            pole_strs.append(f"{p.real:.3f} {'+' if p.imag >= 0 else '-'} j{abs(p.imag):.3f}")
    st.text("\n".join(pole_strs) if pole_strs else "N/A")

    # Display other parameters using columns for better layout
    sub_col1, sub_col2, sub_col3 = st.columns(3)
    with sub_col1:
        st.metric(label="σ (Damping Factor)", value=f"{abs(sigma):.3f}" if zeta_used>=0 else f"{sigma:.3f}!" )
        if zeta_used < 0: st.caption("(Negative indicates instability)")
    with sub_col2:
        st.metric(label="ωd (Damped Freq.)", value=f"{omega_d:.3f} rad/s")
    with sub_col3:
        st.metric(label="θ (Damping Angle)", value=f"{theta:.1f}°")

    st.caption("(θ relative to negative real axis for ζ ≥ 0)")


# --- Plot Display ---
with col2:
    st.subheader("S-Plane Plot")
    if not error_msg:
        fig = generate_s_plane_plot(omega_n_used, zeta_used, poles, sigma, omega_d, theta)
        st.pyplot(fig)
    else:
        st.warning("Plot cannot be generated due to input error.")


st.markdown("---")
st.markdown("Created by an AI Assistant based on Electrical Engineering principles.")