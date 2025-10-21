%% loadflow_decoupled_NR_history_Sbus.m
% Decoupled Newton–Raphson Load-Flow Solver
% Records iteration history and displays apparent power injection per bus

clc; clear; close all;

%% 1) SYSTEM DATA
% bus_data: [bus#, type, V_spec, δ_spec(deg), P_load, Q_load, P_gen, Q_gen]
bus_data = [
    % bus, type, Vsp, δsp, P_L, Q_L, P_G, Q_G
      1, 1, 1.05, 0, 0, 0, 100, 0; % slack
      2, 3, 1.00, 0, 50, 30, 0, 0; % load
      3, 2, 1.02, 0, 20, 10, 50, 0; % PV
      4, 3, 1.00, 0, 80, 50, 0, 0; % load
];
nb = size(bus_data,1);
slack = find(bus_data(:,2)==1); % type 1 is slack bus;
PVb = find(bus_data(:,2)==2); % type 2 is PV bus;
PQb = find(bus_data(:,2)==3); % type 3 is PQ bus;

% line_data: [fromBus, toBus, R+jX, Bc/2, tapRatio]
line_data = [
    % from to R+jX Bc/2 tap
      1, 2, 0.010+0.050j, 0.000j, 1;
      2, 3, 0.012+0.040j, 0.000j, 1;
      3, 4, 0.015+0.045j, 0.000j, 1;
      4, 1, 0.010+0.050j, 0.000j, 1;
];

%% 2) BUILD Ybus
Ybus = formYbus(line_data, nb);

%% 3) INITIAL GUESS
V = bus_data(:,3);
delta = bus_data(:,4)*pi/180;

%% 4) REDUCED JACOBIANS
NonSlack = setdiff(1:nb, slack);
Bp = -imag(Ybus(NonSlack, NonSlack));  
Bpp = -imag(Ybus(PQb, PQb));  
L1 = chol(Bp, 'lower');
L2 = chol(Bpp,'lower');

%% Preallocate history arrays to store V and delta each iteration
maxIter = 20;
V_hist = zeros(nb, maxIter+1);
delta_hist = zeros(nb, maxIter+1);
V_hist(:,1) = V;
delta_hist(:,1) = delta;

%% 5) ITERATIONS
tol = 1e-6; iter = 0;
while iter < maxIter
    % Compute P,Q injections
    [Pcalc, Qcalc] = calcPQ(Ybus, V, delta);
    Pspec = bus_data(:,7) - bus_data(:,5);
    Qspec = bus_data(:,8) - bus_data(:,6);

    % Mismatch
    dP = Pspec(NonSlack) - Pcalc(NonSlack);
    dQ = Qspec(PQb) - Qcalc(PQb);

    if max(abs([dP; dQ])) < tol, break; end

    % Solve updates
    y = L1 \ dP;     
    ddelta = L1' \ y;
    z = L2 \ dQ;     
    dV = L2' \ z;

    % Apply corrections
    delta(NonSlack) = delta(NonSlack) + ddelta;
    V(PQb) = V(PQb) + dV;

    iter = iter + 1;
    % Store history at iteration
    V_hist(:,iter+1) = V;
    delta_hist(:,iter+1) = delta;

    % Display results for this iteration
    fprintf('Iteration %d:\n', iter);
    fprintf('Bus | Voltage (p.u.) | Angle (deg)\n');
    for k = 1:nb
        fprintf('%3d | %8.4f | %8.4f\n', k, V(k), delta(k)*180/pi);
    end
    fprintf('-----------------------------------\n');
end

%% 6) FINAL SUMMARY
fprintf('\nConverged in %d iterations. Final max mismatch: %.2e p.u.\n', iter, max(abs([dP; dQ])));

%% 7) APPARENT POWER INJECTION AT BUSES
% Compute final P,Q injections and apparent power
[Pcalc, Qcalc] = calcPQ(Ybus, V, delta);
Sbus = Pcalc + 1j*Qcalc; % complex power in p.u.
Smag = abs(Sbus); % apparent power magnitude
Sang = angle(Sbus)*180/pi; % apparent power angle (deg)

fprintf('\nBus | P (p.u.) | Q (p.u.) | |S| (p.u.) | ∠S (deg)\n');
for i = 1:nb
    fprintf('%3d | %10.4f | %10.4f | %10.4f | %10.4f\n', i, Pcalc(i), Qcalc(i), Smag(i), Sang(i));
end

%% ————————————— Local Functions —————————————
function Ybus = formYbus(line_data, nb)
    nl = size(line_data,1);
    Y = zeros(nb);
    for k = 1:nl
        i = line_data(k,1); j = line_data(k,2);
        z = line_data(k,3); bsh = line_data(k,4);
        tap = line_data(k,5); y = 1/z;
        Y(i,j) = Y(i,j) - y/tap;
        Y(j,i) = Y(i,j);
        Y(i,i) = Y(i,i) + (y + bsh)/(tap^2);
        Y(j,j) = Y(j,j) + y + bsh;
    end
    Ybus = Y;
end

function [P, Q] = calcPQ(Ybus, V, delta)
    nb = length(V);
    P = zeros(nb,1); Q = zeros(nb,1);
    for i = 1:nb
        for j = 1:nb
            Ymag = abs(Ybus(i,j));
            thetaij = angle(Ybus(i,j));
            P(i) = P(i) + V(i)*V(j)*Ymag * cos(delta(i)-delta(j)-thetaij);
            Q(i) = Q(i) - V(i)*V(j)*Ymag * sin(delta(i)-delta(j)-thetaij);
        end
    end
end
