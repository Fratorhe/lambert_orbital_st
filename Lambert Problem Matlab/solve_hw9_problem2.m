% AE 424 - HW9 - Problem 2
% Solves Lambert's problem for a 6-hour transfer between two positions in a geocentric orbit

clc; clear;

% Parámetros del problema
r1 = [7158.52, 2464.87, 0];      % km
r2 = [-28103.48, -31212.08, 0]; % km
dt = 6 * 3600;                  % 6 horas en segundos
mu = 398600.4418;               % Constante gravitacional Tierra (km^3/s^2)

% Resolver el problema de Lambert
[v1, v2, z_solved] = lambert_solver(r1, r2, dt, mu, true, 0);

% Mostrar resultados
fprintf('\n--- RESULTADOS PARA HW9 - PROBLEMA 2 ---\n');
fprintf('r1 = [%g, %g, %g] km\n', r1);
fprintf('r2 = [%g, %g, %g] km\n', r2);
fprintf('Tiempo de vuelo: %.0f segundos\n', dt);
fprintf('\nVelocidad inicial (v1_tr): [%f, %f, %f] km/s\n', v1);
fprintf('Velocidad final   (v2_tr): [%f, %f, %f] km/s\n', v2);
fprintf('z solucionado: %.4f\n', z_solved);