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

    i = np.arccos(h_vec[2] / h)
    print(f"Orbit inclination is {i:.3f} rad")
    print(f"Orbit inclination is {i/to_rad:.3f} degrees")

    if i > np.pi / 2:
        print(f"Since i greater than 90 degrees, this is a retrograde orbit.")
        print(f"Rotation of object is contrary to rotation of the body")
    else:
        print(f"Since i smaller than 90 degrees, this is a normal orbit.")
        print(f"Rotation of object is same to rotation of the body.")

    N_vec = np.cross([0, 0, 1], h_vec)
    print(N_vec)

    N = np.linalg.norm(N_vec)
    print(f"{N:.3e}")

    Omega = np.arccos(N_vec[0] / N)
    # Omega is returned from [0, pi].
    # However, we should check the quadrant based on N_Y
    if N_vec[1] < 0:
        print("y component of N is negative, therefore, we need to do 360-Omega")
        Omega = 2 * np.pi - Omega
    print(f"Omega is {Omega:.3f} rad")
    print(f"Omega is {Omega/to_rad:.3f} degrees")

    e_vec = 1 / mu * ((v**2 - mu / r) * r_vec - r * v_r * v_vec)
    print(e_vec)
    e = np.linalg.norm(e_vec)
    print(f"The eccentricity of the orbit is: {e:.3f}")
    print_type_orbit(e)

    # perigee argument
    w = np.arccos(np.dot(N_vec, e_vec) / N / e)
    if e_vec[2] < 0:
        print("Solved ambiguity for periapsis")
        w = 2 * np.pi - w
    if i == 0:
        w = arctan2(e_vec[1], e_vec[0])
    print(f"The periapsis argument is {w:.3f} rad")
    print(f"The periapsis argument is {w / to_rad:.3f} degrees")

    # true anomaly
    theta = np.arccos(np.dot(e_vec, r_vec) / e / r)
    if v_r < 0:
        print("Solved ambiguity for true anomaly")
        theta = 2 * np.pi - theta
    print(f"The true anomaly is {theta:.3f} rad")
    print(f"The true anomaly is {theta / to_rad:.3f} degrees")

    a = h**2 / mu / (1 - e**2)

    return h, e, i, Omega, w, theta, mu, a


def compute_state_vector(h, e, i, Omega, w, theta, mu, a):

    r = h**2 / mu / (1 + e * cos(theta))
    print(f"Distance is {r}km")

    # We first get the perifocal coordinates
    r_pf = np.array([r * cos(theta), r * sin(theta), 0])
    v_pf = np.array([-mu / h * sin(theta), mu / h * (e + cos(theta)), 0])

    print("The perifocal coordinates are: ")
    print(f"{r_pf} km")
    print(f"{v_pf} km/s")

    R_Omega = rotz(-Omega)
    R_incl = rotx(-i)
    R_omega = rotz(-w)

    r_vec = np.dot(np.dot(np.dot(R_Omega, R_incl), R_omega), r_pf)
    v_vec = np.dot(np.dot(np.dot(R_Omega, R_incl), R_omega), v_pf)

    print("The geocentric coordinates are: ")
    print(f"{r_vec} km")
    print(f"{v_vec} km/s")
    return r_vec, v_vec, mu


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
