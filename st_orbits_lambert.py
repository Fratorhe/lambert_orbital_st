import numpy as np
import plotly.graph_objects as go
import streamlit as st

from lamberts_problem import LambertAlgorithm
from orbit_elements_and_state_vector import compute_elements

np.set_printoptions(precision=4)

DICT_BODIES = {
    "Earth": {"mu": 398600.4418, "radius": 6371.0, "z_guess": 0.0},
    "Mars": {"mu": 42828.3, "radius": 3389.5, "z_guess": 1},
    "Sun": {"mu": 132712440041.939, "radius": 695700.0, "z_guess": 20},
}

ORBIT_1_COLOR = "blue"
ORBIT_2_COLOR = "green"
ORBIT_TR_COLOR = "purple"
DELTA_V_COLOR = "pink"
ARROWSIZE = 1


def perifocal_to_eci_equatorial(v_pqw, omega_deg):
    omega = np.radians(omega_deg)
    R_z = np.array(
        [[np.cos(omega), -np.sin(omega), 0], [np.sin(omega), np.cos(omega), 0], [0, 0, 1]]
    )
    return R_z @ v_pqw


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


def compute_velocity_at_true_anomaly(a, e, nu_deg, mu, omega_deg):
    h = np.sqrt(mu * a * (1 - e**2))
    nu = np.radians(nu_deg)
    vx = -mu / h * np.sin(nu)
    vy = mu / h * (e + np.cos(nu))
    v = np.array([vx, vy, 0])
    v = perifocal_to_eci_equatorial(v, omega_deg)
    return v


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
st.text(f"r1 is {r1_vec} km")
st.text(f"r2 is {r2_vec} km")

prograde = True
zguess = DICT_BODIES[planet_selected]["z_guess"]
lambert_problem = LambertAlgorithm(r1_vec, r2_vec, delta_time, mu, prograde, zguess)
v1, v2 = lambert_problem.solve_it()
st.text(f"Vtr 1 is {v1} km/s")
st.text(f"Vtr 2 is {v2} km/s")

# compute the transfer orbit
_, e_tr, w_tr, _, a_tr = compute_elements(r1_vec, v1, mu)
x_tr, y_tr = compute_orbit(a_tr, e_tr, w_tr * 180 / np.pi)

# Create Plotly figure
fig = go.Figure()

# Plot orbit 1
fig.add_trace(
    go.Scatter(
        x=x1, y=y1, mode="lines", name="Orbit 1", line=dict(color="blue", dash="dashdot")
    )
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
    go.Scatter(
        x=x2, y=y2, mode="lines", name="Orbit 2", line=dict(color="green", dash="dash")
    )
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
        x=x_tr,
        y=y_tr,
        mode="lines",
        name="Orbit Transfer",
        line=dict(color="purple", dash="dash"),
    )
)


# Compute velocity vector components at current position
v_current = compute_velocity_at_true_anomaly(a1, e1, nu_current, mu, omega1)
st.text(f"	V at orbit 1 {v_current} km/s")
v_desired = compute_velocity_at_true_anomaly(a2, e2, nu_desired, mu, omega2)
st.text(f"	V at 0rbit 2 {v_desired} km/s")
v_current_mag = np.linalg.norm(v_current)
v_desired_mag = np.linalg.norm(v_desired)
v1_mag = np.linalg.norm(v1)  # vmag_transfer orbit at 1
v2_mag = np.linalg.norm(v2)  # vmag_transfer orbit at 2
velocity_scale = a1 / 2  # Scale factor to make the velocity vector visible on the plot


deltaV1 = v1 - v_current
deltaV2 = v_desired - v2
deltaV1_mag = np.linalg.norm(deltaV1)  # vmag_transfer orbit at 1
deltaV2_mag = np.linalg.norm(deltaV2)  # vmag_transfer orbit at 2

st.text(f"	ΔV1 {deltaV1} km/s")
st.text(f"	ΔV2 {deltaV2} km/s")
st.text(f"	ΔV TOTAL {(deltaV2_mag+deltaV1_mag):.3f} km/s")
st.text(f"	z={lambert_problem.z_solved:.3f} ")

# Arrow velocity vector orbit 1
fig.add_annotation(
    x=x_obj1 + velocity_scale * v_current[0] / v_current_mag,
    y=y_obj1 + velocity_scale * v_current[1] / v_current_mag,
    ax=x_obj1,
    ay=y_obj1,
    xref="x",
    yref="y",
    axref="x",
    ayref="y",
    showarrow=True,
    arrowhead=3,  # Arrowhead style (1=small, 3=large)
    arrowsize=ARROWSIZE,  # Size of the arrow
    arrowcolor=ORBIT_1_COLOR,
)

# Arrow velocity vector orbit 2
fig.add_annotation(
    x=x_obj2 + velocity_scale * v_desired[0] / v_desired_mag,
    y=y_obj2 + velocity_scale * v_desired[1] / v_desired_mag,
    ax=x_obj2,
    ay=y_obj2,
    xref="x",
    yref="y",
    axref="x",
    ayref="y",
    showarrow=True,
    arrowhead=3,  # Arrowhead style (1=small, 3=large)
    arrowsize=ARROWSIZE,  # Size of the arrow
    arrowcolor=ORBIT_2_COLOR,
)


# Arrow velocity vector transfer orbit 1
fig.add_annotation(
    x=x_obj1 + velocity_scale * v1[0] / v1_mag,
    y=y_obj1 + velocity_scale * v1[1] / v1_mag,
    ax=x_obj1,
    ay=y_obj1,
    xref="x",
    yref="y",
    axref="x",
    ayref="y",
    showarrow=True,
    arrowhead=3,  # Arrowhead style (1=small, 3=large)
    arrowsize=1,  # Size of the arrow
    arrowcolor=ORBIT_TR_COLOR,
)

# Arrow velocity vector transfer orbit 2
fig.add_annotation(
    x=x_obj2 + velocity_scale * v2[0] / v2_mag,
    y=y_obj2 + velocity_scale * v2[1] / v2_mag,
    ax=x_obj2,
    ay=y_obj2,
    xref="x",
    yref="y",
    axref="x",
    ayref="y",
    showarrow=True,
    arrowhead=3,  # Arrowhead style (1=small, 3=large)
    arrowsize=ARROWSIZE,  # Size of the arrow
    arrowcolor=ORBIT_TR_COLOR,
)

# Arrow velocity DeltaV1
fig.add_annotation(
    x=x_obj1 + velocity_scale * v1[0] / v1_mag,
    y=y_obj1 + velocity_scale * v1[1] / v1_mag,
    ax=x_obj1 + velocity_scale * v_current[0] / v_current_mag,
    ay=y_obj1 + velocity_scale * v_current[1] / v_current_mag,
    xref="x",
    yref="y",
    axref="x",
    ayref="y",
    showarrow=True,
    arrowhead=3,  # Arrowhead style (1=small, 3=large)
    arrowsize=ARROWSIZE,  # Size of the arrow
    arrowcolor=DELTA_V_COLOR,
)

# Arrow velocity deltaV2
fig.add_annotation(
    x=x_obj2 + velocity_scale * v_desired[0] / v_desired_mag,
    y=y_obj2 + velocity_scale * v_desired[1] / v_desired_mag,
    ax=x_obj2 + velocity_scale * v2[0] / v2_mag,
    ay=y_obj2 + velocity_scale * v2[1] / v2_mag,
    xref="x",
    yref="y",
    axref="x",
    ayref="y",
    showarrow=True,
    arrowhead=3,  # Arrowhead style (1=small, 3=large)
    arrowsize=ARROWSIZE,  # Size of the arrow
    arrowcolor=DELTA_V_COLOR,
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

fig2 = go.Figure()
# Plot F(z) function.
z_values = np.linspace(-10, 15, 100)
F_values = [lambert_problem.F_fun(z_value) for z_value in z_values]
fig2.add_trace(
    go.Scatter(
        x=z_values,
        y=F_values,
        mode="lines",
        name="Orbit Transfer",
        line=dict(color="purple", dash="dash"),
    )
)
st.plotly_chart(fig2)


print("done")
# Add credit sentence
st.markdown(
    """
---
This tool was created by Francisco Torres Herrador, Assistant Professor at New Mexico State University, 2025.
"""
)
