## This file implements the methods required for Lambert's algorithm as a series of functions.


from dataclasses import dataclass

import numpy as np
from scipy.optimize import fsolve


def stumpff_S(z):
    if z == 0:
        return 1 / 6

    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (z ** (3 / 2))

    # if z < 0:
    return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / ((-z) ** (3 / 2))


def stumpff_C(z):
    if z == 0:
        return 1 / 2

    if z > 0:
        return (1 - np.cos(np.sqrt(z))) / z

    # if z < 0:
    return (np.cosh(np.sqrt(-z)) - 1) / (-z)


@dataclass
class LambertAlgorithm:
    r1_vec: list | np.ndarray
    r2_vec: list | np.ndarray
    delta_time: float
    mu: float
    prograde: bool = True
    z_guess: float = 0.0

    def __post_init__(self):
        if isinstance(self.r1_vec, list):
            self.r1_vec = np.array(self.r1_vec)
        if isinstance(self.r2_vec, list):
            self.r2_vec = np.array(self.r2_vec)
        self.r1 = np.linalg.norm(self.r1_vec)
        self.r2 = np.linalg.norm(self.r2_vec)
        self.z_solved = np.nan

    @property
    def delta_theta_degree(self):
        """
        Computes the angle theta between two vectors in degrees
        """
        return self.delta_theta * 180 / np.pi

    @property
    def delta_theta(self):
        """
        Computes the angle theta between two vectors provided if the orbit is prograde or retrograde
        """
        # compute the cross product of the two vectors.
        r1_cross_r2 = np.cross(self.r1_vec, self.r2_vec)
        z_cross = r1_cross_r2[2]  # take the z component
        r1_dot_r2 = np.dot(self.r1_vec, self.r2_vec)

        delta_theta = np.arccos(r1_dot_r2 / self.r1 / self.r2)
        # conditions for prograde and retrograde orbits.
        if self.prograde and z_cross < 0:
            delta_theta = 2 * np.pi - delta_theta
        if not self.prograde and z_cross > 0:
            delta_theta = 2 * np.pi - delta_theta
        return delta_theta

    @property
    def A(self):
        return np.sin(self.delta_theta) * np.sqrt(
            self.r1 * self.r2 / (1 - np.cos(self.delta_theta))
        )

    def y_fun(self, z):
        return self.r1 + self.r2 + self.A * (z * stumpff_S(z) - 1) / np.sqrt(stumpff_C(z))

    @property
    def y_zsolved(self):
        return self.y_fun(self.z_solved)

    def F_fun(self, z):
        return (
            (self.y_fun(z) / stumpff_C(z)) ** (3 / 2) * stumpff_S(z)
            + self.A * np.sqrt(self.y_fun(z))
            - np.sqrt(self.mu) * self.delta_time
        )

    def F_fun_prime(self, z):
        if z == 0:
            return np.sqrt(2) / 40 * self.y_fun(z) ** (3 / 2) + self.A / 8 * (
                np.sqrt(self.y_fun(z)) + self.A * np.sqrt(1 / 2 / self.y_fun(z))
            )

        first_summand = (
            (self.y_fun(z) / stumpff_C(z)) ** 3
            / 2
            * (
                (1 / 2 / z) * (stumpff_C(z) - 3 / 2 * stumpff_S(z) / stumpff_C(z))
                + 3 / 4 * stumpff_S(z) ** 2 / stumpff_C(z)
            )
        )
        second_summand = (
            self.A
            / 8
            * (
                3 * stumpff_S(z) / stumpff_C(z) * np.sqrt(self.y_fun(z))
                + self.A * np.sqrt(stumpff_C(z) / self.y_fun(z))
            )
        )
        return first_summand + second_summand

    def solve_z(self, z0):
        print("Iterative solver for z")
        zsol = fsolve(self.F_fun, z0, fprime=self.F_fun_prime)
        self.z_solved = zsol[0]  # save the solution for later use in other methods
        print("Iterative solver finished")
        return self.z_solved

    # LAGRANGE COEFFICIENTS
    @property
    def f_zsolved(self):
        return 1 - self.y_zsolved / self.r1

    @property
    def g_zsolved(self):
        return self.A * np.sqrt(self.y_zsolved / self.mu)

    @property
    def gdot_zsolved(self):
        return 1 - self.y_zsolved / self.r2

    @property
    def v1_vec(self):
        return 1 / self.g_zsolved * (self.r2_vec - self.f_zsolved * self.r1_vec)

    @property
    def v2_vec(self):
        return 1 / self.g_zsolved * (self.gdot_zsolved * self.r2_vec - self.r1_vec)

    def solve_it(self):
        z0 = self.z_guess
        self.solve_z(z0)
        return self.v1_vec, self.v2_vec


if __name__ == "__main__":
    # # EXERCISE 21.3.1 Morano
    # r1_vec = [4700, 9000, 2700]
    # r2_vec = [-24600, 3500, 6000]
    # l = LambertAlgorithm(
    #     r1_vec, r2_vec, delta_time=2 * 60 * 60, mu=3.986e5, prograde=True
    # )
    # print(f"Angle theta (in degrees): {l.delta_theta_degree}")
    # print(f"A is: {l.A}")

    # # # l.plot_F() # this can be used to give an estimate to start newton method
    # # zsol = l.solve_z(z0=0)
    # # print(f"Computed z is {zsol}")

    # v1, v2 = l.solve_it()
    # print(f"Velocity in 1 is {v1}")
    # print(f"Velocity in 2 is {v2}")

    # EXERCISE 21.4.1 Morano

    z = -2
    print(stumpff_S(z))
    print(stumpff_C(z))

    r1_vec = [5657.83, 9799.64, 0]
    r2_vec = [-18290.7, -2776.45, 0]
    delta_time = 4200
    l = LambertAlgorithm(r1_vec, r2_vec, delta_time=delta_time, mu=3.986e5, prograde=True)
    print(f"Angle theta (in degrees): {l.delta_theta_degree}")
    print(f"A is: {l.A}")

    print(f"{l.F_fun(2)=}")
    print(f"{l.F_fun_prime(2)=}")
    print(f"{l.y_fun(2)=}")

    print(f"{l.F_fun(0)=}")
    print(f"{l.F_fun_prime(0)=}")
    print(f"{l.y_fun(0)=}")

    v1, v2 = l.solve_it()
    print(f"Velocity in 1 is {v1}")
    print(f"Velocity in 2 is {v2}")
