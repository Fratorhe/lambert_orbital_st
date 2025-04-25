function [h, e, w, theta, a] = compute_elements(r_vec, v_vec, mu)
    r = norm(r_vec);
    v = norm(v_vec);
    v_r = dot(v_vec, r_vec) / r;

    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);
    fprintf('Angular momentum h = %.3e km^2/s\n', h);

    e_vec = (1/mu) * ((v^2 - mu/r) * r_vec - r * v_r * v_vec);
    e = norm(e_vec);
    disp('Eccentricity vector:');
    disp(e_vec);

    w = atan2(e_vec(2), e_vec(1));
    theta = acos(dot(e_vec, r_vec) / (e * r));
    if v_r < 0
        disp('Solved ambiguity for true anomaly');
        theta = 2 * pi - theta;
    end

    a = h^2 / mu / (1 - e^2);
end