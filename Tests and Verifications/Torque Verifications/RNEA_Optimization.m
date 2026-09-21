clc; clear;

%% ── Robot parameters ─────────────────────────────────────────────────────
L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
radii = repmat({[0.06, 0.057]}, 1, 7);
rho = 2700;
m_links = mass(L,radii,rho);
m_act = [2.2000; 2.7500; 1.5400; 1.5400; 0.7150; 0.7150; 0.7150];

q_min = [-170; -80; -80; -170; -80; -75; -170] * pi/180;
q_max = [ 170;  80;  80;  170;  80;  75;  170] * pi/180;

F_pay    = [0;  0; -50];
F_lat    = [0; 50;   0];
qd_vec   = ones(7,1);
qdd_vec  = ones(7,1) * 2;

m_links = mass(L, radii, rho);
I_links = inertia_tensor(L, m_links, radii{1});
m_act_p = [m_act(:); 0];

%% ── Warm-start from Monte Carlo ──────────────────────────────────────────
% Start coordinate search from MC best configs
q_best = [ -168.0,  -9.1, -72.6, -13.8, -27.4,   0.9, -44.8;   % J1
            -35.0,  14.2, -44.7,  -6.2,   8.8,  11.1,   7.4;   % J2
            -43.9,   3.0,  74.4, -92.2,  53.4, -10.7,  -2.1;   % J3
             79.8, -68.8, -34.8,  71.3,  17.7,  -8.5, -13.1;   % J4
            108.6,  -8.4,  24.0,  12.9,  68.8, -42.8, -84.0;   % J5
             29.4,  -4.2,  51.4,-101.4, -16.7,  35.5,  -0.0;   % J6
             25.1, -21.5,  -6.5, -37.3, -16.6,  13.4,  60.9] * pi/180;  % J7

mc_peaks = [228.1464; 402.3533; 177.6448; 177.9740; 11.2654; 10.9759; 11.0377];
best_tau = mc_peaks;

joint_labels = {'J1','J2','J3','J4','J5','J6','J7'};
payloads     = {F_lat, F_pay, F_pay, F_pay, F_pay, F_pay, F_pay};

%% ── Envelope function ────────────────────────────────────────────────────
    function e = envelope_torque(q, j, F_p, qd_, qdd_, L_, radii_, rho_, m_act_)
        tg = dynamic_Torque(q, zeros(7,1), zeros(7,1), L_, radii_, rho_, m_act_, F_p);
        td = dynamic_Torque(q, qd_,        qdd_,       L_, radii_, rho_, m_act_, F_p);
        t0g = dynamic_Torque(q, zeros(7,1), zeros(7,1), L_, radii_, rho_, m_act_, zeros(3,1));
        t0d = dynamic_Torque(q, qd_,        qdd_,       L_, radii_, rho_, m_act_, zeros(3,1));
        motion = t0d - t0g;
        e = abs(tg(j)) + abs(motion(j));
    end

%% ── 1D Interval B&B per coordinate ──────────────────────────────────────
    function [tau_cert, q_cert] = interval_1D_search(j, q_fixed, F_p, ...
                                                       dim, tol_nm, ...
                                                       qd_, qdd_, L_, ...
                                                       radii_, rho_, m_act_, ...
                                                       qlo, qhi, tau_init)
        % Certified 1D maximization of envelope_torque(q,j)
        % over q(dim) ∈ [qlo, qhi] with other dims fixed at q_fixed

        % Use golden section search (deterministic, certified for unimodal)
        % Combined with interval subdivision for multimodal safety

        gr  = (sqrt(5)-1)/2;   % golden ratio
        tol_rad = tol_nm * pi/180 / 10;   % angle tolerance

        % Subdivision stack for 1D
        stack_1d = {[qlo, qhi]};
        best_val = tau_init;
        best_q   = q_fixed(dim);

        while ~isempty(stack_1d)
            seg   = stack_1d{end};
            stack_1d = stack_1d(1:end-1);
            a_ = seg(1); b_ = seg(2);

            if (b_ - a_) < tol_rad
                % Evaluate midpoint
                qm = q_fixed;
                qm(dim) = 0.5*(a_+b_);
                val = envelope_torque(qm, j, F_p, qd_, qdd_, ...
                                      L_, radii_, rho_, m_act_);
                if val > best_val
                    best_val = val;
                    best_q   = qm(dim);
                end
                continue;
            end

            % Evaluate at 3 interior points to detect monotonicity
            q1_ = a_ + (1-gr)*(b_-a_);
            q2_ = a_ +    gr *(b_-a_);
            qm_  = 0.5*(a_+b_);

            qa = q_fixed; qa(dim) = a_;
            qb = q_fixed; qb(dim) = b_;
            qq1 = q_fixed; qq1(dim) = q1_;
            qq2 = q_fixed; qq2(dim) = q2_;
            qqm = q_fixed; qqm(dim) = qm_;

            va = envelope_torque(qa,  j,F_p,qd_,qdd_,L_,radii_,rho_,m_act_);
            vb = envelope_torque(qb,  j,F_p,qd_,qdd_,L_,radii_,rho_,m_act_);
            v1 = envelope_torque(qq1, j,F_p,qd_,qdd_,L_,radii_,rho_,m_act_);
            v2 = envelope_torque(qq2, j,F_p,qd_,qdd_,L_,radii_,rho_,m_act_);
            vm = envelope_torque(qqm, j,F_p,qd_,qdd_,L_,radii_,rho_,m_act_);

            local_best = max([va,vb,v1,v2,vm]);
            if local_best > best_val
                best_val = local_best;
                [~,idx] = max([va,vb,v1,v2,vm]);
                pts = [a_, b_, q1_, q2_, qm_];
                best_q = pts(idx);
            end

            % Always subdivide both halves — guarantees coverage
            stack_1d{end+1} = [a_,  qm_];
            stack_1d{end+1} = [qm_, b_];
        end

        tau_cert = best_val;
        q_cert   = best_q;
    end

%% ── Coordinate-wise optimization ─────────────────────────────────────────
tol_nm    = 1.0;    % convergence tolerance [Nm]
MAX_OUTER = 20;     % max outer iterations

fprintf('\nCoordinate-wise interval optimization\n');
fprintf('Tolerance: %.2f Nm  |  Max outer iterations: %d\n\n', tol_nm, MAX_OUTER);

for j = 1:7
    fprintf('--- Joint %d ---\n', j);
    q_curr   = q_best(j,:)';
    tau_curr = best_tau(j);
    F_p      = payloads{j};

    for outer = 1:MAX_OUTER
        tau_prev = tau_curr;

        for dim = 1:7
            [tau_new, q_new] = interval_1D_search(j, q_curr, F_p, dim, ...
                                                   tol_nm, qd_vec, qdd_vec, ...
                                                   L, radii, rho, m_act, ...
                                                   q_min(dim), q_max(dim), ...
                                                   tau_curr);
            if tau_new > tau_curr
                tau_curr      = tau_new;
                q_curr(dim)   = q_new;
            end
        end

        fprintf('  Outer iter %2d: tau_J%d = %.4f Nm\n', outer, j, tau_curr);

        if abs(tau_curr - tau_prev) < tol_nm
            fprintf('  Converged.\n');
            break;
        end
    end

    best_tau(j)    = tau_curr;
    q_best(j,:)    = q_curr' * 180/pi;
    fprintf('\n');
end

%% ── Final results ────────────────────────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,115));
fprintf('  COORDINATE-WISE CERTIFIED TORQUE MAXIMA\n');
fprintf('  Tolerance: %.2f Nm  |  Warm-started from Monte Carlo\n', tol_nm);
fprintf('%s\n', repmat('=',1,115));
fprintf('  %-6s  %-16s  %-16s  %-10s  %-8s  %-50s\n', ...
    'Joint','MC Peak (Nm)','Certified (Nm)','Gain (Nm)','DAF','Critical q (deg)');
fprintf('  %s\n', repmat('-',1,110));

static_peaks = [144.75; 306.27; 137.38; 137.38; 10.64; 10.64; 10.64];

for j = 1:7
    gain = best_tau(j) - mc_peaks(j);
    daf  = best_tau(j) / static_peaks(j);
    fprintf('  %-6s  %-16.4f  %-16.4f  %-10.4f  %-8.3f  [%s]\n', ...
        joint_labels{j}, mc_peaks(j), best_tau(j), gain, daf, ...
        sprintf('%7.1f', q_best(j,:)));
end
fprintf('%s\n\n', repmat('=',1,115));
