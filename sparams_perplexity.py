# splane_visualizer.py
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc

def second_order_underdamped_system(zeta, wn):
    """
    Creates the mathematical representation of an under-damped second-order system.
    """
    if not (0 < zeta < 1):
        raise ValueError("Damping ratio must be between 0 and 1 for under-damped systems")
    
    sigma = zeta * wn
    omega_d = wn * np.sqrt(1 - zeta**2)
    theta = np.arccos(zeta)
    
    s1 = complex(-sigma, omega_d)
    s2 = complex(-sigma, -omega_d)
    
    return {
        "zeta": zeta,
        "wn": wn,
        "sigma": sigma,
        "omega_d": omega_d,
        "theta": theta,
        "poles": [s1, s2]
    }

def plot_splane(system):
    """Plot S-plane representation with key parameters."""
    wn = system['wn']
    zeta = system['zeta']
    sigma = system['sigma']
    omega_d = system['omega_d']
    theta = system['theta']
    pole = system['poles'][0]

    fig, ax = plt.subplots(figsize=(10, 8))
    limit = wn * 1.5
    ax.set_xlim(-limit, limit/2)
    ax.set_ylim(-limit, limit)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.axhline(0, color='k', lw=1.5)
    ax.axvline(0, color='k', lw=1.5)

    # Plot poles
    ax.plot(pole.real, pole.imag, 'rx', markersize=12)
    ax.plot(pole.real, -pole.imag, 'rx', markersize=12)

    # Natural frequency circle
    ax.add_patch(Circle((0, 0), wn, fill=False, color='blue', linestyle='--', lw=2))

    # Pole vector and annotations
    ax.plot([0, pole.real], [0, pole.imag], 'g-', lw=2)
    ax.plot([pole.real, 0], [pole.imag, pole.imag], 'r-', lw=2)
    ax.plot([0, 0], [0, omega_d], 'm-', lw=2)
    ax.add_patch(Arc((0, 0), wn, wn, theta1=0, theta2=np.degrees(theta), color='green', lw=2))

    # Annotations
    ax.annotate('ωn', (wn/2 * np.cos(np.pi/4), wn/2 * np.sin(np.pi/4)), 
               xytext=(wn/2 * np.cos(np.pi/4)+0.5, wn/2 * np.sin(np.pi/4)+0.5),
               fontsize=14, color='blue', weight='bold', arrowprops=dict(facecolor='blue', width=2))
    ax.annotate('ωd', (0, omega_d), xytext=(-1.5, omega_d), 
               fontsize=14, color='magenta', weight='bold', arrowprops=dict(facecolor='magenta', width=2))
    ax.annotate('σ', (pole.real, 0), xytext=(pole.real, -1), 
               fontsize=14, color='red', weight='bold', arrowprops=dict(facecolor='red', width=2))
    ax.annotate('θ', (wn/3 * np.cos(theta/2), wn/3 * np.sin(theta/2)), 
               xytext=(wn/3 * np.cos(theta/2)-0.5, wn/3 * np.sin(theta/2)+0.5),
               fontsize=14, color='green', weight='bold', arrowprops=dict(facecolor='green', width=2))

    # Equations
    plt.figtext(0.7, 0.9, r"$T(s) = \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2}$", fontsize=14)
    plt.figtext(0.7, 0.85, f"$s_{{1,2}} = -{sigma:.2f} \pm j{omega_d:.2f}$", fontsize=14)
    
    plt.tight_layout()
    return fig

def plot_time_response(system):
    """Plot step response with envelope."""
    wn = system['wn']
    zeta = system['zeta']
    sigma = system['sigma']
    omega_d = system['omega_d']

    fig, ax = plt.subplots(figsize=(10, 6))
    t = np.linspace(0, 5, 1000)
    y = 1 - np.exp(-sigma * t) * (np.cos(omega_d * t) + (zeta/np.sqrt(1-zeta**2)) * np.sin(omega_d * t))
    
    ax.plot(t, y, 'b-', lw=2)
    ax.plot(t, 1 + np.exp(-sigma * t), 'r--', alpha=0.6)
    ax.plot(t, 1 - np.exp(-sigma * t), 'r--', alpha=0.6)
    
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Amplitude', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    return fig

# Streamlit App Configuration
st.set_page_config(
    page_title="S-Plane Visualizer",
    page_icon="📊",
    layout="wide"
)

# Sidebar Information
with st.sidebar:
    st.header("Key Concepts")
    st.markdown("""
    **Transfer Function:**
    $T(s) = \\frac{\\omega_n^2}{s^2 + 2\\zeta\\omega_n s + \\omega_n^2}$
    
    **Parameters:**
    - ωn: Natural frequency
    - ζ: Damping ratio (0 < ζ < 1)
    - ωd: Damped frequency = ωn√(1-ζ²)
    - σ: Damping factor = ζωn
    - θ: Angle = arccos(ζ)
    """)

# Main App
st.title("Second-Order System S-Plane Visualizer")
st.markdown("Explore under-damped system behavior in the complex frequency domain")

# Controls
col1, col2 = st.columns(2)
with col1:
    wn = st.slider("Natural Frequency (ωn)", 1.0, 10.0, 5.0, 0.1)
with col2:
    zeta = st.slider("Damping Ratio (ζ)", 0.01, 0.99, 0.3, 0.01)

# Generate system
try:
    system = second_order_underdamped_system(zeta, wn)
    
    # Visualizations
    st.header("S-Plane Representation")
    st.pyplot(plot_splane(system))
    
    st.header("Time Domain Response")
    st.pyplot(plot_time_response(system))
    
    # Parameters
    st.header("System Parameters")
    st.markdown(f"""
    - Natural Frequency (ωn): `{system['wn']:.2f} rad/s`
    - Damping Ratio (ζ): `{system['zeta']:.3f}`
    - Damped Frequency (ωd): `{system['omega_d']:.2f} rad/s`
    - Damping Factor (σ): `{system['sigma']:.2f}`
    - Phase Angle (θ): `{np.degrees(system['theta']):.1f}°`
    """)
    
    # Signal Integrity Section
    st.header("Signal Integrity Implications")
    st.markdown("""
    - **Ringing**: Oscillations from stored energy in reactive components
    - **Overshoot**: Exceeding target voltage due to low damping
    - **Settling Time**: Governed by σ (decay rate)
    - **Practical Applications**: Transmission lines, RLC circuits, clock distribution networks
    """)

except ValueError as e:
    st.error(str(e))
    