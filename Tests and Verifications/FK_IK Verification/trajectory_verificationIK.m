clc; clear;

%% ── Robot parameters ─────────────────────────────────────────────────────
L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
radii = repmat({[0.06, 0.057]}, 1, 7);
rho = 2700;
m_links = mass(L,radii,rho);
m_act = [2.2000; 2.7500; 1.5400; 1.5400; 0.7150; 0.7150; 0.7150];
m_act = [2.2; 2.75; 1.54; 1.54; 0.715; 0.715; 0.715];

%% ── Joint limits ─────────────────────────────────────────────────────────
q_min = [-170; -80; -80; -170; -80; -75; -170] * pi/180;
q_max = [ 170;  80;  80;  170;  80;  75;  170] * pi/180;

%% ── Monte Carlo parameters ───────────────────────────────────────────────
N_samples = 200000;
F_pay     = [0;  0; -50];
F_lat     = [0; 50;   0];
qd_vec    = ones(7,1);
qdd_vec   = ones(7,1) * 2;

%% ── Sample joint space uniformly ─────────────────────────────────────────
Q = q_min + (q_max - q_min) .* rand(7, N_samples);

%% ── Evaluate torques ─────────────────────────────────────────────────────
% Three load cases:
%   LC1: vertical payload, static
%   LC2: vertical payload, dynamic (qd=1, qdd=2)
%   LC3: lateral payload,  static   (for J1)
%   LC4: lateral payload,  dynamic  (for J1)

tau_static_pay  = zeros(7, N_samples);   % LC1
tau_dyn_pay     = zeros(7, N_samples);   % LC2
tau_static_lat  = zeros(7, N_samples);   % LC3
tau_dyn_lat     = zeros(7, N_samples);   % LC4
tau_motion_only = zeros(7, N_samples);   % Coriolis + inertial (no payload)

fprintf('Running Monte Carlo: %d samples...\n', N_samples);
tic;
for k = 1:N_samples
    q = Q(:,k);

    tau_static_pay(:,k)  = dynamic_Torque(q, zeros(7,1), zeros(7,1), ...
                                           L, radii, rho, m_act, F_pay);
    tau_dyn_pay(:,k)     = dynamic_Torque(q, qd_vec,     qdd_vec,    ...
                                           L, radii, rho, m_act, F_pay);
    tau_static_lat(:,k)  = dynamic_Torque(q, zeros(7,1), zeros(7,1), ...
                                           L, radii, rho, m_act, F_lat);
    tau_dyn_lat(:,k)     = dynamic_Torque(q, qd_vec,     qdd_vec,    ...
                                           L, radii, rho, m_act, F_lat);

    % Pure motion increment (no payload)
    tg = dynamic_Torque(q, zeros(7,1), zeros(7,1), L, radii, rho, m_act, zeros(3,1));
    td = dynamic_Torque(q, qd_vec,     qdd_vec,    L, radii, rho, m_act, zeros(3,1));
    tau_motion_only(:,k) = td - tg;
end
fprintf('Done in %.1f s\n\n', toc);

%% ── Envelope: worst-case superposition ───────────────────────────────────
% For J1: use lateral load case (worst for yaw)
% For J2-J7: use vertical load case
tau_env_pay = abs(tau_static_pay) + abs(tau_motion_only);
tau_env_lat = abs(tau_static_lat) + abs(tau_motion_only);

% Combined envelope: take max of both load cases per joint
tau_envelope = max(tau_env_pay, tau_env_lat);

%% ── Find peak per joint ───────────────────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,115));
fprintf('  MONTE CARLO PEAK TORQUE ANALYSIS  |  %d samples  |  qd=1 rad/s  |  qdd=2 rad/s2\n', N_samples);
fprintf('%s\n', repmat('=',1,115));
fprintf('  %-4s  %-14s  %-14s  %-8s  %-60s\n', ...
    'Jnt','Static (Nm)','Envelope (Nm)','DAF','Critical q (deg)');
fprintf('  %s\n', repmat('-',1,110));

joint_labels = {'J1','J2','J3','J4','J5','J6','J7'};
q_critical   = zeros(7,7);   % row=joint, col=q values at peak

for j = 1:7
    % Peak static (vertical load)
    [tau_s_peak, idx_s] = max(abs(tau_static_pay(j,:)));

    % Check if lateral is worse for this joint
    [tau_s_lat_peak, idx_s_lat] = max(abs(tau_static_lat(j,:)));
    if tau_s_lat_peak > tau_s_peak
        tau_s_peak = tau_s_lat_peak;
        idx_s = idx_s_lat;
    end

    % Peak envelope
    [tau_e_peak, idx_e] = max(tau_envelope(j,:));
    q_crit = Q(:, idx_e) * 180/pi;
    q_critical(j,:) = q_crit';

    daf_j = tau_e_peak / (tau_s_peak + 1e-10);

    fprintf('  %-4s  %-14.4f  %-14.4f  %-8.3f  [%s]\n', ...
        joint_labels{j}, tau_s_peak, tau_e_peak, daf_j, ...
        sprintf('%7.1f', q_crit'));
end
fprintf('%s\n\n', repmat('=',1,115));

%% ── Per-joint detailed breakdown at critical config ──────────────────────
fprintf('\n%s\n', repmat('=',1,90));
fprintf('  TORQUE DECOMPOSITION AT CRITICAL CONFIGURATION PER JOINT\n');
fprintf('%s\n', repmat('=',1,90));

for j = 1:7
    q_c = q_critical(j,:)' * pi/180;

    tg_pay = dynamic_Torque(q_c, zeros(7,1), zeros(7,1), L, radii, rho, m_act, F_pay);
    tg_lat = dynamic_Torque(q_c, zeros(7,1), zeros(7,1), L, radii, rho, m_act, F_lat);
    tg_0   = dynamic_Torque(q_c, zeros(7,1), zeros(7,1), L, radii, rho, m_act, zeros(3,1));
    td_0   = dynamic_Torque(q_c, qd_vec,     qdd_vec,    L, radii, rho, m_act, zeros(3,1));

    tau_grav    = tg_0;                   % gravity only
    tau_pay_con = tg_pay - tg_0;          % payload contribution (vertical)
    tau_lat_con = tg_lat - tg_0;          % payload contribution (lateral)
    tau_mot     = td_0 - tg_0;            % Coriolis + inertial

    fprintf('\n  %s — Critical config: [%s] deg\n', ...
        joint_labels{j}, sprintf('%6.1f', q_critical(j,:)));
    fprintf('  %s\n', repmat('-',1,85));
    fprintf('  %-4s  %-12s  %-12s  %-12s  %-12s  %-12s\n', ...
        'Jnt','Gravity','Pay(vert)','Pay(lat)','Motion','Envelope');
    fprintf('  %s\n', repmat('-',1,75));
    for jj = 1:7
        env_j = abs(tg_pay(jj)) + abs(tau_mot(jj));
        fprintf('  %-4s  %-12.4f  %-12.4f  %-12.4f  %-12.4f  %-12.4f\n', ...
            joint_labels{jj}, tau_grav(jj), tau_pay_con(jj), ...
            tau_lat_con(jj), tau_mot(jj), env_j);
    end
end
fprintf('\n%s\n\n', repmat('=',1,90));

%% ── Plots ────────────────────────────────────────────────────────────────
figure('Color','w','Position',[50 50 1400 900]);

for j = 1:7
    subplot(4,2,j);

    % Distribution of envelope torque for this joint
    histogram(tau_envelope(j,:), 80, 'FaceColor', [0.2 0.4 0.8], ...
              'EdgeColor', 'none', 'Normalization', 'probability');
    hold on;

    % Mark peak
    [tau_e_peak, ~] = max(tau_envelope(j,:));
    xline(tau_e_peak, 'r-', 'LineWidth', 2, ...
          'Label', sprintf('Peak: %.1f Nm', tau_e_peak), ...
          'LabelVerticalAlignment','bottom');

    % Mark paper config value
    xlabel('Envelope Torque (Nm)');
    ylabel('Probability');
    title(sprintf('%s — Envelope Torque Distribution', joint_labels{j}));
    grid on; box off;
end
sgtitle(sprintf('Monte Carlo Envelope Torque Distributions  |  N = %d', N_samples), ...
    'FontSize', 12);

%% ── Heatmap: which joints most affect which joint torque ─────────────────
figure('Color','w','Position',[50 50 1200 500]);

% Correlation between joint angles and peak envelope torques
corr_mat = zeros(7,7);
for j = 1:7
    for jq = 1:7
        corr_mat(j,jq) = corr(Q(jq,:)', tau_envelope(j,:)');
    end
end

imagesc(abs(corr_mat));
colormap(hot); colorbar;
clim([0 1]);
xlabel('Joint Angle q_i');
ylabel('Joint Torque \tau_j');
title('|Correlation| between joint angles and envelope torques');
xticks(1:7); xticklabels({'q1','q2','q3','q4','q5','q6','q7'});
yticks(1:7); yticklabels({'τ1','τ2','τ3','τ4','τ5','τ6','τ7'});
set(gca,'FontSize',11);

%% ── Joint angle distributions at peak torque configs ─────────────────────
figure('Color','w','Position',[50 50 1200 600]);

% For each joint, find top 1% of envelope samples and plot q distribution
top_pct = ceil(0.01 * N_samples);
for j = 1:7
    subplot(2,4,j);
    [~, idx_sorted] = sort(tau_envelope(j,:), 'descend');
    idx_top = idx_sorted(1:top_pct);
    q_top   = Q(:, idx_top) * 180/pi;   % 7 × top_pct

    % Box plot of joint angles in top 1% samples
    boxplot(q_top', 'Labels', {'q1','q2','q3','q4','q5','q6','q7'});
    ylabel('Angle (deg)');
    title(sprintf('%s — q at top 1%% envelope', joint_labels{j}));
    grid on; box off;
end
sgtitle('Joint angle distributions at peak envelope torque configurations', ...
    'FontSize', 11);
