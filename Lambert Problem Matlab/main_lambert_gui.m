function main_lambert_gui()
    bodies = struct( ...
        'Earth', struct('mu', 398600.4418, 'radius', 6371, 'z_guess', 0.0), ...
        'Mars', struct('mu', 42828.3, 'radius', 3389.5, 'z_guess', 1.0), ...
        'Sun', struct('mu', 132712440041.939, 'radius', 695700, 'z_guess', 20.0) ...
    );
    planets = fieldnames(bodies);
    f = figure('Position', [100, 100, 900, 600], 'Name', 'Lambert Orbital GUI', 'NumberTitle', 'off');

    uicontrol(f, 'Style', 'text', 'Position', [30, 560, 150, 20], 'String', 'Select Planet:');
    planetMenu = uicontrol(f, 'Style', 'popupmenu', 'String', planets, 'Position', [30, 540, 150, 25]);

    uicontrol(f, 'Style', 'text', 'Position', [30, 510, 150, 20], 'String', 'Time (s):');
    dtInput = uicontrol(f, 'Style', 'edit', 'String', '10000', 'Position', [30, 490, 150, 25]);

    uicontrol(f, 'Style', 'text', 'Position', [30, 450, 150, 20], 'String', 'Orbit 1: a, e, ω, ν');
    orb1Input = uicontrol(f, 'Style', 'edit', 'String', '13000, 0.3, 0, 0', 'Position', [30, 430, 150, 25]);

    uicontrol(f, 'Style', 'text', 'Position', [30, 400, 150, 20], 'String', 'Orbit 2: a, e, ω, ν');
    orb2Input = uicontrol(f, 'Style', 'edit', 'String', '16000, 0.3, 0, 20', 'Position', [30, 380, 150, 25]);

    uicontrol(f, 'Style', 'pushbutton', 'String', 'Solve Lambert', ...
        'Position', [30, 340, 150, 30], 'Callback', @(~, ~) solve_and_plot());

    ax = axes(f, 'Position', [0.35, 0.1, 0.6, 0.8]);
    axis equal; grid on;

    function solve_and_plot()
        cla(ax);
        planetName = planets{get(planetMenu, 'Value')};
        body = bodies.(planetName);
        mu = body.mu;
        planet_radius = body.radius;
        z_guess = body.z_guess;

        dt = str2double(get(dtInput, 'String'));
        orb1 = str2num(get(orb1Input, 'String')); %#ok<ST2NM>
        orb2 = str2num(get(orb2Input, 'String')); %#ok<ST2NM>

        [a1, e1, w1, nu1] = deal(orb1(1), orb1(2), deg2rad(orb1(3)), deg2rad(orb1(4)));
        [a2, e2, w2, nu2] = deal(orb2(1), orb2(2), deg2rad(orb2(3)), deg2rad(orb2(4)));

        [x1, y1] = orbit_xy(a1, e1, w1);
        [x2, y2] = orbit_xy(a2, e2, w2);
        [x_obj1, y_obj1] = true_anomaly_xy(a1, e1, w1, nu1);
        [x_obj2, y_obj2] = true_anomaly_xy(a2, e2, w2, nu2);

        r1 = [x_obj1, y_obj1, 0];
        r2 = [x_obj2, y_obj2, 0];

        [v1, v2, z_solved] = lambert_solver(r1, r2, dt, mu, true, z_guess);
        v1_actual = velocity_at_true_anomaly(a1, e1, nu1, mu, w1);
        v2_actual = velocity_at_true_anomaly(a2, e2, nu2, mu, w2);

        deltaV1 = norm(v1 - v1_actual);
        deltaV2 = norm(v2_actual - v2);
        deltaV_total = deltaV1 + deltaV2;

        fprintf('\n--- RESULTS ---\n');
        fprintf('v1 (transfer): [%f %f %f] km/s\n', v1);
        fprintf('v2 (transfer): [%f %f %f] km/s\n', v2);
        fprintf('ΔV1: %.3f km/s\n', deltaV1);
        fprintf('ΔV2: %.3f km/s\n', deltaV2);
        fprintf('ΔV TOTAL: %.3f km/s\n', deltaV_total);
        fprintf('z solved: %.3f\n', z_solved);

        plot(ax, x1, y1, '--b'); hold(ax, 'on');
        plot(ax, x2, y2, '--g');
        plot(ax, x_obj1, y_obj1, 'mo', 'MarkerFaceColor', 'm');
        plot(ax, x_obj2, y_obj2, 'ro', 'MarkerFaceColor', 'r');

        [~, ~, w_tr, ~, a_tr] = compute_elements(r1, v1, mu);
        [xtr, ytr] = orbit_xy(a_tr, norm((cross(r1, v1)))/mu, w_tr);
        plot(ax, xtr, ytr, '-m');

        theta = linspace(0, 2*pi, 100);
        fill(ax, planet_radius*cos(theta), planet_radius*sin(theta), [0.8 0.9 1], 'EdgeColor', 'k');
    end
end

function [x, y] = orbit_xy(a, e, omega)
    nu = linspace(0, 2*pi, 300);
    r = a * (1 - e^2) ./ (1 + e * cos(nu));
    theta = omega + nu;
    x = r .* cos(theta);
    y = r .* sin(theta);
end

function [x, y] = true_anomaly_xy(a, e, omega, nu)
    r = a * (1 - e^2) / (1 + e * cos(nu));
    theta = omega + nu;
    x = r * cos(theta);
    y = r * sin(theta);
end

function v = velocity_at_true_anomaly(a, e, nu, mu, omega)
    h = sqrt(mu * a * (1 - e^2));
    vr = -mu / h * sin(nu);
    vt = mu / h * (e + cos(nu));
    v_pqw = [vr, vt, 0];
    Rz = [cos(omega), -sin(omega), 0; sin(omega), cos(omega), 0; 0, 0, 1];
    v = (Rz * v_pqw')';
end