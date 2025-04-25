function [v1, v2, z_solved] = lambert_solver(r1_vec, r2_vec, delta_t, mu, prograde, z_guess)
    if nargin < 6
        z_guess = 0;
    end
    if nargin < 5
        prograde = true;
    end

    r1 = norm(r1_vec);
    r2 = norm(r2_vec);
    z = z_guess;

    cross_r1_r2 = cross(r1_vec, r2_vec);
    z_cross = cross_r1_r2(3);
    theta = acos(dot(r1_vec, r2_vec) / (r1 * r2));
    if prograde && z_cross < 0
        theta = 2 * pi - theta;
    elseif ~prograde && z_cross > 0
        theta = 2 * pi - theta;
    end

    A = sin(theta) * sqrt(r1 * r2 / (1 - cos(theta)));

    y_fun = @(z) r1 + r2 + A * (z * stumpff_S(z) - 1) / sqrt(stumpff_C(z));
    F_fun = @(z) (y_fun(z) / stumpff_C(z))^(3/2) * stumpff_S(z) + A * sqrt(y_fun(z)) - sqrt(mu) * delta_t;

    F_prime = @(z) derivative_F(z, y_fun, A);
    options = optimset('Display','off');
    z_solved = fsolve(F_fun, z, options);

    y = y_fun(z_solved);
    f = 1 - y / r1;
    g = A * sqrt(y / mu);
    gdot = 1 - y / r2;

    v1 = (r2_vec - f * r1_vec) / g;
    v2 = (gdot * r2_vec - r1_vec) / g;
end

function dF = derivative_F(z, y_fun, A)
    y = y_fun(z);
    if z == 0
        dF = sqrt(2)/40 * y^(3/2) + A/8 * (sqrt(y) + A * sqrt(1 / (2 * y)));
    else
        C = stumpff_C(z);
        S = stumpff_S(z);
        term1 = (y / C)^3 / 2 * ((1 / (2 * z)) * (C - (3/2) * S / C) + (3/4) * S^2 / C);
        term2 = A / 8 * (3 * S / C * sqrt(y) + A * sqrt(C / y));
        dF = term1 + term2;
    end
end