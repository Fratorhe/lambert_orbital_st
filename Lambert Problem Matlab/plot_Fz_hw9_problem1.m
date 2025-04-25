% AE 424 - HW9 - Problem 1 (b)
% Plot F(z) vs z for z in [-5, 5]

clc; clear;

% Parameters
a = 14300;
e = 0.3;
theta_R = deg2rad(60);
theta_Q = deg2rad(135);
mu = 398600.4418;
dt = 70 * 60; % 4200 s

% Positions R & Q
[xR, yR] = true_anomaly_xy(a, e, 0, theta_R);
[xQ, yQ] = true_anomaly_xy(a, e, 0, theta_Q);
r1 = [xR, yR, 0];
r2 = [xQ, yQ, 0];

theta = acos(dot(r1, r2)/(norm(r1)*norm(r2)));
A = sin(theta) * sqrt(norm(r1) * norm(r2)/(1 - cos(theta)));

% F(z)
z_vals = linspace(-5, 5, 1000);
F_vals = arrayfun(@(z) F_of_z(z, A, r1, r2, mu, dt), z_vals);

% Grph
figure;
plot(z_vals, F_vals, 'LineWidth', 1.5);
xlabel('z'); ylabel('F(z)');
title('F(z) vs z for HW9 Problem 1 (b)');
grid on;

% Función F(z)
function Fz = F_of_z(z, A, r1, r2, mu, dt)
    y = norm(r1) + norm(r2) + A * (z * stumpff_S(z) - 1) / sqrt(stumpff_C(z));
    Fz = (y / stumpff_C(z))^(3/2) * stumpff_S(z) + A * sqrt(y) - sqrt(mu) * dt;
end

% Función auxiliar de posición
function [x, y] = true_anomaly_xy(a, e, omega, nu)
    r = a * (1 - e^2) / (1 + e * cos(nu));
    theta = omega + nu;
    x = r * cos(theta);
    y = r * sin(theta);
end