import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import signal # Import signal processing library

st.set_page_config(layout="wide")

st.title("Interactive S-Plane & Time Domain Analysis of 2nd-Order System")
st.write("""
This app visualizes the S-plane pole locations and the corresponding time-domain step response
for a stable, underdamped second-order system. Adjust the parameters using the sliders
to observe the relationship between pole locations and transient behavior.
The system's transfer function is: `H(s) = ωn² / (s² + 2ζωn s + ωn²)`.
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

# Calculate derived parameters for S-plane
sigma = zeta * omega_n          # Damping factor (magnitude of real part)
omega_d = omega_n * np.sqrt(1 - zeta**2) # Damped natural frequency
theta_rad = np.arccos(zeta)     # Angle in radians
theta_deg = np.degrees(theta_rad) # Angle in degrees

# Pole locations
pole1 = -sigma + 1j * omega_d
pole2 = -sigma - 1j * omega_d

# --- Time Domain Calculations ---
# Define the transfer function H(s) = omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2)
num = [omega_n**2]
den = [1, 2 * zeta * omega_n, omega_n**2]
system = signal.TransferFunction(num, den)

# Calculate approximate settling time (2% criterion) for time vector generation
t_settle_approx = 4 / sigma if sigma > 1e-6 else 100 # Avoid division by zero, provide fallback

# Define time vector for simulation - go a bit beyond settling time
t_end = max(1.5 * t_settle_approx, 5 * (2 * np.pi / omega_n) if omega_n > 1e-6 else 10) # Ensure a few cycles or settling visible
t_end = min(t_end, 100) # Limit max simulation time
t = np.linspace(0, t_end, 500) # 500 points for smooth curve

# Calculate step response
t, y = signal.step(system, T=t)

# Calculate time-domain metrics
overshoot_pct = (np.max(y) / 1.0 - 1) * 100 if np.max(y) > 1 else 0.0

# --- Plotting Function: S-Plane ---
def plot_s_plane(sigma_val, omega_d_val, omega_n_val, zeta_val, theta_deg_val, pole1_val, pole2_val):
    """Generates the S-plane matplotlib plot."""
    fig, ax = plt.subplots(figsize=(7, 7)) # Adjusted size slightly

    max_lim = omega_n_val * 1.4
    ax.set_xlim(-max_lim, max_lim * 0.5)
    ax.set_ylim(-max_lim, max_lim)
    ax.set_aspect('equal', adjustable='box')

    ax.axhline(0, color='black', linewidth=0.75)
    ax.axvline(0, color='black', linewidth=0.75)
    ax.plot(max_lim * 0.5, 0, ">k", markersize=5, clip_on=False)
    ax.plot(0, max_lim, "^k", markersize=5, clip_on=False)

    ax.plot(pole1_val.real, pole1_val.imag, 'rx', markersize=10, markeredgewidth=2, label=f'Pole s1 ({pole1_val:.2f})')
    ax.plot(pole2_val.real, pole2_val.imag, 'bx', markersize=10, markeredgewidth=2, label=f'Pole s2 ({pole2_val:.2f})')

    ax.plot([0, pole1_val.real], [0, pole1_val.imag], 'm--', linewidth=1.5, label=f'ωn = {omega_n_val:.2f}')
    ax.plot([pole1_val.real, pole1_val.real], [0, pole1_val.imag], 'g:', linewidth=1.5)
    ax.plot([0, pole1_val.real], [pole1_val.imag, pole1_val.imag], 'g:', linewidth=1.5)
    ax.plot([0, pole1_val.real], [0, 0], 'g-', linewidth=2.0, label=f'σ = {sigma_val:.2f}')
    ax.plot([0, 0], [0, pole1_val.imag], 'c-', linewidth=2.0, label=f'ωd = {omega_d_val:.2f}')

    ax.text(pole1_val.real / 2, pole1_val.imag * 1.05, f'ωn={omega_n_val:.2f}', color='m', ha='center', va='bottom', fontsize=9)
    ax.text(pole1_val.real, pole1_val.imag / 2, f'ωd={omega_d_val:.2f}', color='c', ha='right', va='center', rotation=90, fontsize=9)
    ax.text(pole1_val.real / 2, 0.05 * max_lim, f'σ={sigma_val:.2f}', color='g', ha='center', va='bottom', fontsize=9)
    ax.text(pole1_val.real, pole1_val.imag + 0.05 * max_lim, 's1', color='r', ha='center', va='bottom', weight='bold', fontsize=9)
    ax.text(pole2_val.real, pole2_val.imag - 0.1 * max_lim, 's2', color='b', ha='center', va='top', weight='bold', fontsize=9)

    arc_radius = omega_n_val * 0.3
    theta_arc = patches.Arc((0, 0), 2 * arc_radius, 2 * arc_radius, angle=0,
                            theta1=180 - theta_deg_val, theta2=180, color='purple', linewidth=1.5,
                            label=f'θ = {theta_deg_val:.1f}°')
    ax.add_patch(theta_arc)
    angle_label_rad = np.radians(180 - theta_deg_val / 2)
    label_dist_factor = 0.35
    ax.text(arc_radius * label_dist_factor * np.cos(angle_label_rad),
            arc_radius * label_dist_factor * np.sin(angle_label_rad),
            f'θ={theta_deg_val:.1f}°', color='purple', ha='center', va='center', fontsize=9)

    ax.set_xlabel('Real Axis (σ)')
    ax.set_ylabel('Imaginary Axis (jω)')
    ax.set_title('S-Plane Pole Locations')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.plot(0,0, 'ko', markersize=4)

    return fig

# --- Plotting Function: Time Domain ---
def plot_time_response(t_vec, y_vec, t_settle_est, overshoot_val):
    """Generates the time-domain step response plot."""
    fig, ax = plt.subplots(figsize=(7, 4)) # Adjusted size slightly

    ax.plot(t_vec, y_vec, 'b-', linewidth=2, label='System Response y(t)')
    ax.axhline(1.0, color='grey', linestyle='--', linewidth=1.0, label='Final Value (1.0)')

    # Indicate Overshoot if significant
    if overshoot_val > 1.0: # Only show if overshoot > 1%
        max_y = np.max(y_vec)
        max_t = t_vec[np.argmax(y_vec)]
        ax.plot([max_t, max_t], [1.0, max_y], 'r--', linewidth=1.0, label=f'Overshoot: {overshoot_val:.1f}%')
        ax.plot(max_t, max_y, 'ro', markersize=5)

    # Indicate approx settling time
    ax.axvline(t_settle_est, color='orange', linestyle=':', linewidth=1.5, label=f'~Settling Time (Ts≈4/σ) = {t_settle_est:.2f}s')

    # Dynamic Y limits
    min_y = 0
    max_y_plot = max(1.1, np.max(y_vec) * 1.1) # Ensure final value is visible, add padding
    ax.set_ylim(min_y, max_y_plot)
    ax.set_xlim(0, t_vec[-1])

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude')
    ax.set_title('Step Response')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize='small', loc='lower right')

    return fig


# --- Display Results ---
col1, col2 = st.columns([2, 1]) # Plots take more space on left, details on right

with col1:
    st.subheader("S-Plane Visualization")
    fig_s_plane = plot_s_plane(sigma, omega_d, omega_n, zeta, theta_deg, pole1, pole2)
    st.pyplot(fig_s_plane)

    st.subheader("Time Domain Step Response")
    fig_time = plot_time_response(t, y, t_settle_approx, overshoot_pct)
    st.pyplot(fig_time)


with col2:
    st.subheader("S-Plane Parameters")
    st.metric(label="Damping Factor (σ)", value=f"{sigma:.4f}")
    st.metric(label="Damped Frequency (ωd)", value=f"{omega_d:.4f} rad/s")
    st.metric(label="Angle (θ)", value=f"{theta_deg:.2f} degrees")
    st.write("---")
    st.write("**Pole Locations:**")
    st.latex(f"s_1 = -\\sigma + j\\omega_d = {pole1.real:.4f} + {pole1.imag:.4f}j")
    st.latex(f"s_2 = -\\sigma - j\\omega_d = {pole2.real:.4f} + {pole2.imag:.4f}j")
    st.write("---")

    st.subheader("Time Domain Metrics")
    st.metric(label="Peak Overshoot", value=f"{overshoot_pct:.2f} %")
    st.metric(label="Approx. Settling Time (Ts ≈ 4/σ)", value=f"{t_settle_approx:.3f} s")
    st.write("*(Settling time is estimated for the response to stay within ±2% of the final value)*")
    st.write("---")

    # Interpretation hints
    st.subheader("Interpretation")
    st.markdown(f"""
    *   **Poles closer to the jω axis** (smaller `σ`, smaller `ζ`, `θ` closer to 90°) lead to **more oscillations (ringing)** and **higher overshoot** in the time response.
    *   **Poles further left** (larger `σ`) lead to **faster decay** of oscillations and a **shorter settling time**.
    *   The **imaginary part (`ωd`)** dictates the **frequency of the oscillations** seen in the time response.
    *   The **distance from the origin (`ωn`)** relates to the overall speed of the response.
    """)
    if zeta < 0.2:
        st.info("Low Damping (ζ < 0.2): Significant overshoot and ringing.")
    elif zeta < 0.707: # sqrt(2)/2 is often considered boundary for noticeable peaking
        st.info("Moderate Damping (0.2 ≤ ζ < 0.707): Noticeable overshoot and some ringing.")
    else:
        st.info("High Damping (0.707 ≤ ζ < 1.0): Little to no overshoot, response settles quickly without significant oscillation.")


st.sidebar.write("---")
st.sidebar.write("Created by an AI Assistant")
st.sidebar.write("[SciPy Signal Docs](https://docs.scipy.org/doc/scipy/reference/signal.html)")