import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

st.set_page_config(layout="wide")

st.title("Interactive S-Plane Plot for Underdamped Second-Order System")
st.write("""
This app visualizes the pole locations of a stable, underdamped second-order system
in the S-plane (complex frequency plane). Adjust the parameters using the sliders
to see how they affect the pole positions and related geometric interpretations.
The system's characteristic equation is: `s² + 2ζωn s + ωn² = 0`.
""")

# --- User Inputs ---
st.sidebar.header("System Parameters")
omega_n_input = st.sidebar.slider(
    'Natural Frequency (ωn)',
    min_value=0.1,
    max_value=10.0,
    value=2.0,
    step=0.1,
    help="Undamped oscillation frequency (rad/s)"
)

zeta_input = st.sidebar.slider(
    'Damping Ratio (ζ)',
    min_value=0.01, # Strictly > 0 for damped frequency calculation
    max_value=0.99, # Strictly < 1 for underdamped
    value=0.3,
    step=0.01,
    help="Determines the level of damping (0 < ζ < 1 for underdamped)"
)

# --- Calculations ---
omega_n = omega_n_input
zeta = zeta_input

# Calculate derived parameters
sigma = zeta * omega_n          # Damping factor (magnitude of real part)
omega_d = omega_n * np.sqrt(1 - zeta**2) # Damped natural frequency
theta_rad = np.arccos(zeta)     # Angle in radians
theta_deg = np.degrees(theta_rad) # Angle in degrees

# Pole locations
pole1 = -sigma + 1j * omega_d
pole2 = -sigma - 1j * omega_d

# --- Plotting Function ---
def plot_s_plane(sigma_val, omega_d_val, omega_n_val, zeta_val, theta_deg_val, pole1_val, pole2_val):
    """Generates the S-plane matplotlib plot."""
    fig, ax = plt.subplots(figsize=(8, 8))

    # Determine plot limits dynamically with padding
    max_lim = omega_n_val * 1.4 # Ensure omega_n radius is visible
    ax.set_xlim(-max_lim, max_lim * 0.5) # More space on negative real axis
    ax.set_ylim(-max_lim, max_lim)
    ax.set_aspect('equal', adjustable='box') # Keep aspect ratio 1:1

    # Draw Axes
    ax.axhline(0, color='black', linewidth=0.75)
    ax.axvline(0, color='black', linewidth=0.75)

    # Add arrow heads to axes
    ax.plot(max_lim * 0.5, 0, ">k", markersize=5, clip_on=False)
    ax.plot(0, max_lim, "^k", markersize=5, clip_on=False)

    # Plot Poles
    ax.plot(pole1_val.real, pole1_val.imag, 'rx', markersize=10, markeredgewidth=2, label=f'Pole s1 ({pole1_val:.2f})')
    ax.plot(pole2_val.real, pole2_val.imag, 'bx', markersize=10, markeredgewidth=2, label=f'Pole s2 ({pole2_val:.2f})')

    # Plot Geometric Lines (Use different colors and styles)
    # omega_n (radius)
    ax.plot([0, pole1_val.real], [0, pole1_val.imag], 'm--', linewidth=1.5, label=f'ωn = {omega_n_val:.2f}')
    # sigma (damping factor - real part)
    ax.plot([pole1_val.real, pole1_val.real], [0, pole1_val.imag], 'g:', linewidth=1.5) # Vertical dashed from real axis
    ax.plot([0, pole1_val.real], [pole1_val.imag, pole1_val.imag], 'g:', linewidth=1.5) # Horizontal dashed from imag axis
    ax.plot([0, pole1_val.real], [0, 0], 'g-', linewidth=2.0, label=f'σ = {sigma_val:.2f}') # Line on real axis
    # omega_d (damped frequency - imag part)
    ax.plot([0, 0], [0, pole1_val.imag], 'c-', linewidth=2.0, label=f'ωd = {omega_d_val:.2f}') # Line on imag axis

    # Annotations
    ax.text(pole1_val.real / 2, pole1_val.imag * 1.05, f'ωn={omega_n_val:.2f}', color='m', ha='center', va='bottom')
    ax.text(pole1_val.real, pole1_val.imag / 2, f'ωd={omega_d_val:.2f}', color='c', ha='right', va='center', rotation=90)
    ax.text(pole1_val.real / 2, 0.05 * max_lim, f'σ={sigma_val:.2f}', color='g', ha='center', va='bottom')
    ax.text(pole1_val.real, pole1_val.imag + 0.05 * max_lim, 's1', color='r', ha='center', va='bottom', weight='bold')
    ax.text(pole2_val.real, pole2_val.imag - 0.1 * max_lim, 's2', color='b', ha='center', va='top', weight='bold')


    # Angle Theta
    arc_radius = omega_n_val * 0.3
    theta_arc = patches.Arc((0, 0), 2 * arc_radius, 2 * arc_radius, angle=0,
                            theta1=180 - theta_deg_val, theta2=180, color='purple', linewidth=1.5,
                            label=f'θ = {theta_deg_val:.1f}°')
    ax.add_patch(theta_arc)
    # Angle label position needs adjustment based on angle size
    angle_label_rad = np.radians(180 - theta_deg_val / 2)
    label_dist_factor = 0.35 # Adjust this factor as needed
    ax.text(arc_radius * label_dist_factor * np.cos(angle_label_rad),
            arc_radius * label_dist_factor * np.sin(angle_label_rad),
            f'θ={theta_deg_val:.1f}°', color='purple', ha='center', va='center', fontsize=9)

    # Labels and Title
    ax.set_xlabel('Real Axis (σ)')
    ax.set_ylabel('Imaginary Axis (jω)')
    ax.set_title('S-Plane Pole Locations')
    ax.grid(True, linestyle=':', alpha=0.6)
    # ax.legend(loc='upper right', fontsize='small') # Optional: Legend if labels aren't clear enough

    # Add zero point marker
    ax.plot(0,0, 'ko', markersize=4)

    return fig

# --- Display Results ---
col1, col2 = st.columns([2, 1]) # Plot takes more space

with col1:
    st.subheader("S-Plane Visualization")
    fig = plot_s_plane(sigma, omega_d, omega_n, zeta, theta_deg, pole1, pole2)
    st.pyplot(fig)

with col2:
    st.subheader("Calculated Values")
    st.metric(label="Damping Factor (σ)", value=f"{sigma:.4f}")
    st.metric(label="Damped Frequency (ωd)", value=f"{omega_d:.4f} rad/s")
    st.metric(label="Angle (θ)", value=f"{theta_deg:.2f} degrees")
    st.write("---")
    st.write("**Pole Locations:**")
    st.latex(f"s_1 = -\\sigma + j\\omega_d = {pole1.real:.4f} + {pole1.imag:.4f}j")
    st.latex(f"s_2 = -\\sigma - j\\omega_d = {pole2.real:.4f} + {pole2.imag:.4f}j")
    st.write("---")
    # Interpretation hints
    st.subheader("Interpretation")
    if zeta < 0.2:
        st.info("Low Damping (ζ < 0.2): Expect significant overshoot and ringing (oscillations) in the time-domain step response.")
    elif zeta < 0.7:
        st.info("Moderate Damping (0.2 ≤ ζ < 0.7): Some overshoot and ringing, but settles reasonably fast.")
    else:
        st.info("High Damping (0.7 ≤ ζ < 1.0): Little overshoot, approaches critical damping (fastest settling without oscillation). Poles are closer to the real axis.")

    st.write(f"Poles are located at a distance **ωn = {omega_n:.2f}** from the origin.")
    st.write(f"The real part **-σ = {-sigma:.2f}** determines the decay rate (settling time ≈ 4/σ).")
    st.write(f"The imaginary part **±ωd = ±{omega_d:.2f}** determines the oscillation frequency during transients.")


st.sidebar.write("---")
st.sidebar.write("Created by an AI Assistant")