clc; clear; rng(42);

%% ── Robot parameters ──────────────────────────────────────────────────────
L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
radii  = repmat({[0.06, 0.057]}, 1, 7);
rho    = 2700;
m_act  = [2.2; 2.75; 1.54; 1.54; 0.715; 0.715; 0.715];
m_links = mass(L, radii, rho);

F_cases = {[0;  0; -50], '-Z 50N';
           [0; 50;   0], '+Y 50N'};

N       = 100000;
qd_val  = 1.0;   % rad/s  — fixed for all joints, all samples
qdd_val = 2.0;   % rad/s² — fixed for all joints, all samples

%% ── Joint limits ──────────────────────────────────────────────────────────
q_min = [-170; -80; -80; -170; -80; -75; -170] * pi/180;
q_max = [ 170;  80;  80;  170;  80;  75;  170] * pi/180;

%% ── Fixed qd and qdd vectors ──────────────────────────────────────────────
% All joints simultaneously at the specified values.
% Sign combinations: 2^7 = 128 possible sign patterns.
% Worst case is found by trying both +/- for each joint across samples.
% Simplest conservative approach: fix magnitude, randomise signs.
qd_pos  =  qd_val  * ones(7, 1);
qd_neg  = -qd_val  * ones(7, 1);
qdd_pos =  qdd_val * ones(7, 1);
qdd_neg = -qdd_val * ones(7, 1);

%% ── Preallocate ───────────────────────────────────────────────────────────
tau_peak_gravity  = zeros(7, 1);
tau_peak_payload  = zeros(7, 1);
tau_peak_coriolis = zeros(7, 1);
tau_peak_inertial = zeros(7, 1);
tau_peak_total    = zeros(7, 1);

q_star_total    = zeros(7, 7);
qd_star_total   = zeros(7, 7);
qdd_star_total  = zeros(7, 7);

%% ── Monte Carlo loop ──────────────────────────────────────────────────────
fprintf('Running %d Monte Carlo samples at qd = %.1f rad/s, qdd = %.1f rad/s²...\n', ...
    N, qd_val, qdd_val);

for s = 1:N

    % Random configuration only
    q = q_min + (q_max - q_min) .* rand(7, 1);

    % Random sign pattern for qd and qdd — covers all ± combinations
    sgn_qd  = sign(2*rand(7,1) - 1);   % each element +1 or -1
    sgn_qdd = sign(2*rand(7,1) - 1);

    qd  = qd_val  * sgn_qd;
    qdd = qdd_val * sgn_qdd;

    % ── Gravity only (qd=0, qdd=0) ──
    tau_g = dynamic_Torque(q, zeros(7,1), zeros(7,1), L, radii, rho, m_act);

    % ── Coriolis only ──
    tau_gc = dynamic_Torque(q, qd, zeros(7,1), L, radii, rho, m_act);
    tau_c  = tau_gc - tau_g;

    % ── Inertial only ──
    tau_gi = dynamic_Torque(q, zeros(7,1), qdd, L, radii, rho, m_act);
    tau_i  = tau_gi - tau_g;

    % ── Payload: worst case over both force directions ──
    J     = Geometric_Jacobian(q, L);
    tau_p = zeros(7, 1);
    for f = 1:size(F_cases, 1)
        tp    = J(1:3,:)' * F_cases{f,1};
        tau_p = sign(tp) .* max(abs(tau_p), abs(tp));
    end

    % ── Total ──
    tau_total = tau_g + tau_c + tau_i + tau_p;

    % ── Update per-joint peaks ──
    for j = 1:7
        if abs(tau_g(j)) > tau_peak_gravity(j)
            tau_peak_gravity(j) = abs(tau_g(j));
        end
        if abs(tau_p(j)) > tau_peak_payload(j)
            tau_peak_payload(j) = abs(tau_p(j));
        end
        if abs(tau_c(j)) > tau_peak_coriolis(j)
            tau_peak_coriolis(j) = abs(tau_c(j));
        end
        if abs(tau_i(j)) > tau_peak_inertial(j)
            tau_peak_inertial(j) = abs(tau_i(j));
        end
        if abs(tau_total(j)) > tau_peak_total(j)
            tau_peak_total(j)   = abs(tau_total(j));
            q_star_total(:, j)  = q;
            qd_star_total(:, j) = qd;
            qdd_star_total(:,j) = qdd;
        end
    end
end

fprintf('Done.\n\n');

%% ── DAF ───────────────────────────────────────────────────────────────────
tau_static_peak = tau_peak_gravity + tau_peak_payload;
DAF = tau_peak_total ./ (tau_static_peak + 1e-10);

%% ── Main results table ────────────────────────────────────────────────────
fprintf('%s\n', repmat('=',1,115));
fprintf('  Monte Carlo  |  N = %d  |  qd = %.1f rad/s (fixed)  |  qdd = %.1f rad/s² (fixed)  |  Payload = 50 N\n', ...
    N, qd_val, qdd_val);
fprintf('%s\n', repmat('=',1,115));
fprintf('  %-6s  %-12s  %-12s  %-12s  %-12s  %-12s  %-12s  %-8s  %s\n', ...
    'Joint','Peak Total','Gravity','Payload','Coriolis','Inertial','Static','DAF','Verdict');
fprintf('  %s\n', repmat('-',1,110));

for j = 1:7
    if tau_static_peak(j) < 0.01 * max(tau_static_peak)
        verdict = 'Purely dynamic';
        daf_str = 'N/A';
    else
        daf = DAF(j);
        if     daf > 1.5, verdict = 'Dynamic critical';
        elseif daf > 1.2, verdict = 'Dynamic significant';
        else,              verdict = 'Static sufficient';
        end
        daf_str = sprintf('%.3f', daf);
    end

    fprintf('  J%-5d  %-12.4f  %-12.4f  %-12.4f  %-12.4f  %-12.4f  %-12.4f  %-8s  %s\n', ...
        j, tau_peak_total(j), tau_peak_gravity(j), tau_peak_payload(j), ...
        tau_peak_coriolis(j), tau_peak_inertial(j), tau_static_peak(j), ...
        daf_str, verdict);
end
fprintf('%s\n\n', repmat('=',1,115));

%% ── Static vs Dynamic summary ─────────────────────────────────────────────
fprintf('%s\n', repmat('=',1,60));
fprintf('  Static vs Dynamic Peak Torques\n');
fprintf('%s\n', repmat('=',1,60));
fprintf('  %-6s  %-14s  %-14s  %-8s\n', ...
    'Joint','Static (Nm)','Dynamic (Nm)','DAF');
fprintf('  %s\n', repmat('-',1,48));
for j = 1:7
    fprintf('  J%-5d  %-14.4f  %-14.4f  %.3f\n', ...
        j, tau_static_peak(j), tau_peak_total(j), DAF(j));
end
fprintf('%s\n\n', repmat('=',1,60));

%% ── Critical configurations ───────────────────────────────────────────────
fprintf('%s\n', repmat('=',1,100));
fprintf('  Critical Configurations (q* at peak total torque per joint)\n');
fprintf('%s\n', repmat('=',1,100));
fprintf('  %-6s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Joint','q1','q2','q3','q4','q5','q6','q7');
fprintf('  %s\n', repmat('-',1,95));
for j = 1:7
    q_deg = q_star_total(:, j) * 180/pi;
    fprintf('  J%-5d  %-10.2f  %-10.2f  %-10.2f  %-10.2f  %-10.2f  %-10.2f  %-10.2f\n', ...
        j, q_deg(1), q_deg(2), q_deg(3), q_deg(4), q_deg(5), q_deg(6), q_deg(7));
end
fprintf('%s\n\n', repmat('=',1,100));

%% ── Bar chart ─────────────────────────────────────────────────────────────
figure('Color','w','Position',[100 100 900 500]);
joint_labels = {'J1 (Y)','J2 (P)','J3 (R)','J4 (P)','J5 (R)','J6 (Y)','J7 (P)'};
x = 1:7;

b1 = bar(x,        tau_static_peak, 0.35, 'FaceColor', [0.2 0.5 0.8]);
hold on;
b2 = bar(x + 0.38, tau_peak_total,  0.35, 'FaceColor', [0.8 0.2 0.2]);

xlabel('Joint'); ylabel('Peak Torque (Nm)');
title(sprintf('Peak Torque — Static vs Dynamic  |  qd = %.1f rad/s  |  qdd = %.1f rad/s²', ...
    qd_val, qdd_val));
set(gca, 'XTick', x + 0.19, 'XTickLabel', joint_labels);
legend([b1 b2], {'Static (gravity + payload)', 'Dynamic (total)'}, ...
    'Location', 'northeast');
grid on; box off;