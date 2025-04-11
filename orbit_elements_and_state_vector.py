import numpy as np
from numpy import arctan, arctan2, cos, sin

to_rad = np.pi / 180


def print_type_orbit(eccentricity):
    if eccentricity < 1:
        print("Elliptical Orbit")
    else:
        print("Hyperbolic Orbit")


def rotx(angle):
    Rx = np.array(
        [
            [1, 0, 0],
            [0, cos(angle), sin(angle)],
            [0, -sin(angle), cos(angle)],
        ]
    )
    return Rx


def rotz(angle):
    Rz = np.array(
        [
            [cos(angle), sin(angle), 0],
            [-sin(angle), cos(angle), 0],
            [0, 0, 1],
        ]
    )
    return Rz


def compute_elements(r_vec, v_vec, mu):
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    v_r = np.dot(v_vec, r_vec) / r

    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    print(f"Angular momentum h is {h:.3e} km^2/s")

    e_vec = 1 / mu * ((v**2 - mu / r) * r_vec - r * v_r * v_vec)
    print(e_vec)
    e = np.linalg.norm(e_vec)

    w = arctan2(e_vec[1], e_vec[0])

    # true anomaly
    theta = np.arccos(np.dot(e_vec, r_vec) / e / r)
    if v_r < 0:
        print("Solved ambiguity for true anomaly")
        theta = 2 * np.pi - theta

    a = h**2 / mu / (1 - e**2)

    return h, e, w, theta, a


if __name__ == "__main__":
    # mu_sun = 1.32712440018e11  # km3/s2, for the sun

    # r_vec = np.array(
    #     [1.376978580204302e08, -6.238873772246556e07, 1.696893847996943e07]
    # )  # km
    # v_vec = np.array(
    #     [1.162927402282732e01, 2.724856969064504e01, 8.669287252356082e-01]
    # )  # km/s

    # r_vec = np.array([132.47e6, 67.533e6, 3.9022e6])  # km
    # v_vec = np.array([-14.010, 26.678, 0.10287])  # km/s

    # elements = compute_elements(r_vec, v_vec, mu_sun)
    # r_vec, v_vec, mu = compute_state_vector(*elements)

    mu_earth = 3.986e5  # km3/s2, for the sun

    r_vec = np.array([5657.83, 9799.64, 0])  # km
    v_vec = np.array([-7.28424433, 2.15804309, 0.0])  # km/s

    elements = compute_elements(r_vec, v_vec, mu_earth)
    # r_vec, v_vec, mu = compute_state_vector(*elements)

    print(elements)
