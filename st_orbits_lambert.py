import numpy as np
import plotly.graph_objects as go
import streamlit as st
from orbit_elements_and_state_vector import compute_elements

from lamberts_problem import LambertAlgorithm

DICT_BODIES = {
    "Earth": {"mu": 398600.4418, "radius": 6371.0},
    "Mars": {"mu": 42828.3, "radius": 3389.5},
    "Sun": {"mu": 132712440041.939, "radius": 695700.0},
}


# Orbital equation in polar coordinates
def compute_orbit(a, e, omega_deg, num_points=500):
    omega = np.radians(omega_deg)
    nu = np.linspace(0, 2 * np.pi, num_points)
    r = a * (1 - e**2) / (1 + e * np.cos(nu))

    theta = omega + nu
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


# Compute satellite position at a given true anomaly
def compute_position_at_true_anomaly(a, e, omega_deg, nu_deg):
    omega = np.radians(omega_deg)
    nu = np.radians(nu_deg)
    r = a * (1 - e**2) / (1 + e * np.cos(nu))
    theta = omega + nu
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


# Streamlit app


planet_selected = st.sidebar.selectbox(
    "Which celestial body are we orbiting?",
    DICT_BODIES.keys(),
    index=0,
    # placeholder="Earth",
)
st.title(f"{planet_selected}-Centered Planar Orbit Visualizer: Lambert's Problem")

mu = DICT_BODIES[planet_selected]["mu"]
planet_radius = st.sidebar.number_input(
    "Planet Radius [km]", value=DICT_BODIES[planet_selected]["radius"]
)

# Delta time.
delta_time = st.sidebar.number_input(
    "Time [s]", min_value=1000, value=10000, step=1000, format="%d"
)


# ORBIT 1
st.sidebar.header("Orbit 1")
a1 = st.sidebar.number_input(
    "Semi-major axis a [km]", min_value=planet_radius + 100, value=planet_radius * 2
)
e1 = st.sidebar.slider("Eccentricity e", 0.0, 0.99, 0.3, step=0.01)
omega1 = st.sidebar.slider("Argument of periapsis ω [deg]", 0, 360, 0)
nu_current = st.sidebar.slider("True anomaly ν [deg]", 0, 360, 0)

# Compute orbit and object position
x1, y1 = compute_orbit(a1, e1, omega1)
x_obj1, y_obj1 = compute_position_at_true_anomaly(a1, e1, omega1, nu_current)

# ORBIT 2
st.sidebar.header("Orbit 2")
a2 = st.sidebar.number_input(
    "Semi-major axis a [km]",
    min_value=planet_radius + 100,
    value=planet_radius * 2.5,
    key="a_orbit2",
)
e2 = st.sidebar.slider("Eccentricity e", 0.0, 0.99, 0.3, step=0.01, key="e_orbit2")
omega2 = st.sidebar.slider("Argument of periapsis ω [deg]", 0, 360, 0, key="w_orbit2")
nu_desired = st.sidebar.slider("True anomaly ν [deg]", 0, 360, 20, key="nu_orbit2")

# Compute orbit and object position
x2, y2 = compute_orbit(a2, e2, omega2)
x_obj2, y_obj2 = compute_position_at_true_anomaly(a2, e2, omega2, nu_desired)


# Compute Lambert's problem
r1_vec = np.array([x_obj1, y_obj1, 0])
r2_vec = np.array([x_obj2, y_obj2, 0])

prograde = True

lambert_problem = LambertAlgorithm(r1_vec, r2_vec, delta_time, mu, prograde)
v1, v2 = lambert_problem.solve_it()
st.text(f"Transfer velocity in 1 is {v1}")
st.text(f"Transfer velocity in 2 is {v2}")

# compute the transfer orbit
_, e_tr, _, _, w_tr, _, _, a_tr = compute_elements(r1_vec, v1, mu)
x_tr, y_tr = compute_orbit(a_tr, e_tr, w_tr * 180 / np.pi)
print(w_tr)

# Create Plotly figure
fig = go.Figure()

# Plot orbit 1
fig.add_trace(
    go.Scatter(x=x1, y=y1, mode="lines", name="Orbit 1", line=dict(color="blue"))
)

# Plot current satellite position
fig.add_trace(
    go.Scatter(
        x=[x_obj1],
        y=[y_obj1],
        mode="markers+text",
        marker=dict(size=10, color="magenta"),
        text=["Current Position"],
        name="Current Position",
        textposition="top center",
    )
)

# Plot orbit 2
fig.add_trace(
    go.Scatter(x=x2, y=y2, mode="lines", name="Orbit 2", line=dict(color="green"))
)

# Plot NEW satellite position
fig.add_trace(
    go.Scatter(
        x=[x_obj2],
        y=[y_obj2],
        mode="markers+text",
        marker=dict(size=10, color="red"),
        text=["New Position"],
        name="New Position",
        textposition="top center",
    )
)

# Plot transfer orbit
fig.add_trace(
    go.Scatter(
        x=x_tr, y=y_tr, mode="lines", name="Orbit Transfer", line=dict(color="purple")
    )
)


# Plot Earth as a filled circle
theta_planet = np.linspace(0, 2 * np.pi, 200)
planet_x = planet_radius * np.cos(theta_planet)
planet_y = planet_radius * np.sin(theta_planet)
fig.add_trace(
    go.Scatter(
        x=planet_x,
        y=planet_y,
        fill="toself",
        name=planet_selected,
        fillcolor="rgba(24, 158, 149, 0.2)",
        line=dict(color="rgba(24, 158, 149, 1)"),
    )
)

# Plot settings
fig.update_layout(
    title="Orbital Trajectory",
    xaxis_title="x [km]",
    yaxis_title="y [km]",
    width=700,
    height=700,
    showlegend=True,
    xaxis=dict(scaleanchor="y", scaleratio=1),
    yaxis=dict(scaleanchor="x", scaleratio=1),
)

st.plotly_chart(fig)
