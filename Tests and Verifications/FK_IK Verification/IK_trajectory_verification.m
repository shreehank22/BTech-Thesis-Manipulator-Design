%% IK_Trajectory_Test.m
clear; clc;

L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
radii = repmat({[0.06, 0.057]}, 1, 7);
rho = 2700;
m_links = mass(L,radii,rho);
m_act = [2.2000; 2.7500; 1.5400; 1.5400; 0.7150; 0.7150; 0.7150];
q_min = [-170,-80,-80,-170,-80,-75,-170]' * pi/180;
q_max = [ 170, 80, 80, 170, 80, 75, 170]' * pi/180;
psi_vals = linspace(0, 2*pi, 16);
N = 100;

function Q_all = collect_solutions(R_des, p_des, psi_vals, L)
    Q_all = [];
    warning('off','all');
    for k = 1:numel(psi_vals)
        Q_k = Inverse_Kinematics2(R_des, p_des, psi_vals(k), L);
        Q_all = [Q_all, Q_k];
    end
    warning('on','all');
end

function e = logSO3(R)
    tr = max(-1,min(1,(trace(R)-1)/2));
    th = acos(tr);
    if abs(th)<1e-7
        skew = 0.5*(R-R');
    elseif abs(th-pi)<1e-4
        A=(R+eye(3))/2; n=sqrt(max(0,diag(A)));
        if A(1,2)<0, n(2)=-n(2); end
        if A(1,3)<0, n(3)=-n(3); end
        e=pi*n; return
    else
        skew=(th/(2*sin(th)))*(R-R');
    end
    e=[skew(3,2);skew(1,3);skew(2,1)];
end

%% ── Define start and end configurations ──────────────────────────────────
q_start = [30, 45, -20,  60,  15, -30,  10]' * pi/180;
q_end   = [-30, 20,  45, -60, -15,  20, -10]' * pi/180;

[T_start,~,~] = Forward_Kinematics(q_start, L);
[T_end,~,~]   = Forward_Kinematics(q_end,   L);

p_start = T_start(1:3,4); R_start = T_start(1:3,1:3);
p_end   = T_end(1:3,4);   R_end   = T_end(1:3,1:3);

fprintf('Straight-line trajectory: %d waypoints\n', N);
fprintf('Start: [%.3f %.3f %.3f] m\n', p_start);
fprintf('End:   [%.3f %.3f %.3f] m\n\n', p_end);

%% ── Generate waypoints ────────────────────────────────────────────────────
t_vals = linspace(0,1,N);
waypoints_p = zeros(3,N);
waypoints_R = zeros(3,3,N);
e_rel = logSO3(R_start'*R_end); angle_rel = norm(e_rel);
for i = 1:N
    t = t_vals(i);
    waypoints_p(:,i) = (1-t)*p_start + t*p_end;
    if angle_rel < 1e-8
        waypoints_R(:,:,i) = R_start;
    else
        ax = e_rel/angle_rel;
        K  = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
        waypoints_R(:,:,i) = R_start*(eye(3)+sin(t*angle_rel)*K+(1-cos(t*angle_rel))*(K*K));
    end
end

%% ── Solve IK at each waypoint ─────────────────────────────────────────────
q_traj   = zeros(7, N);
solve_ok = false(1, N);
q_prev   = q_start;
w        = [1; 1; 1; 5; 10; 5; 10];

for i = 1:N
    p_i = waypoints_p(:,i);
    R_i = waypoints_R(:,:,i);
    Q_i = collect_solutions(R_i, p_i, psi_vals, L);
    if isempty(Q_i)
        continue;
    end
    best_score = inf; best_idx = 1;
    for s = 1:size(Q_i,2)
        dq        = Q_i(:,s) - q_prev;
        dq_w      = norm(w .* dq);
        violation = sum(max(0, q_min-Q_i(:,s)) + max(0, Q_i(:,s)-q_max));
        score     = dq_w + 1e3*violation;
        if score < best_score
            best_score = score; best_idx = s;
        end
    end
    q_traj(:,i) = Q_i(:,best_idx);
    dq = q_traj(:,i) - q_prev;
    q_traj(:,i) = q_traj(:,i) - 2*pi*round(dq/(2*pi));
    for j = 1:7
        if q_traj(j,i) < q_min(j) && q_traj(j,i)+2*pi <= q_max(j)
            q_traj(j,i) = q_traj(j,i) + 2*pi;
        elseif q_traj(j,i) > q_max(j) && q_traj(j,i)-2*pi >= q_min(j)
            q_traj(j,i) = q_traj(j,i) - 2*pi;
        end
    end
    solve_ok(i) = true;
    q_prev = q_traj(:,i);
end

%% ── TEST 1: IK solvability ────────────────────────────────────────────────
fprintf('=== TEST 1: IK solvability ===\n');
n_solved = sum(solve_ok);
fprintf('Solved: %d/%d waypoints\n', n_solved, N);
if n_solved==N
    fprintf('Status: PASS\n\n');
else
    fprintf('Failed at: %s\n', mat2str(find(~solve_ok)));
    fprintf('Status: FAIL\n\n');
end

%% ── TEST 2: FK round-trip ─────────────────────────────────────────────────
fprintf('=== TEST 2: FK round-trip error ===\n');
fprintf('%-10s  %-12s  %-12s  %s\n','Waypoint','pos_err','rot_err','Status');
fprintf('%s\n',repmat('-',1,48));
fk_pass=0; max_pe=0; max_re=0;
for i = 1:N
    if ~solve_ok(i), continue; end
    [Tc,~,~] = Forward_Kinematics(q_traj(:,i), L);
    pe = norm(Tc(1:3,4)-waypoints_p(:,i));
    re = norm(logSO3(Tc(1:3,1:3)*waypoints_R(:,:,i)'));
    max_pe=max(max_pe,pe); max_re=max(max_re,re);
    if pe<1e-4 && re<1e-4
        fk_pass=fk_pass+1;
    else
        fprintf('%-10d  %-12.2e  %-12.2e  FAIL\n',i,pe,re);
    end
end
fprintf('Max pos_err: %.2e m   Max rot_err: %.2e rad\n',max_pe,max_re);
fprintf('FK check: %d/%d passed\n',fk_pass,n_solved);
if fk_pass==n_solved, fprintf('Status: PASS\n\n');
else, fprintf('Status: FAIL\n\n'); end

%% ── TEST 3: Joint continuity ──────────────────────────────────────────────
fprintf('=== TEST 3: Joint continuity ===\n');
fprintf('%-10s  %-12s  %-12s  %s\n','Step','dq_norm','dq_max','Status');
fprintf('%s\n',repmat('-',1,48));
cont_pass=0; cont_total=0; max_jump=0;
dq_threshold=0.30;
for i = 2:N
    if ~solve_ok(i)||~solve_ok(i-1), continue; end
    dq=q_traj(:,i)-q_traj(:,i-1);
    dq_norm=norm(dq); dq_max=max(abs(dq));
    max_jump=max(max_jump,dq_norm);
    cont_total=cont_total+1;
    if dq_norm<dq_threshold
        cont_pass=cont_pass+1;
    else
        fprintf('%-10d  %-12.4f  %-12.4f  FAIL\n',i,dq_norm,dq_max);
    end
end
fprintf('Max step jump: %.4f rad\n',max_jump);
fprintf('Continuity: %d/%d steps passed (threshold %.2f rad)\n', ...
    cont_pass,cont_total,dq_threshold);
if cont_pass==cont_total, fprintf('Status: PASS\n\n');
else, fprintf('Status: FAIL\n\n'); end

%% ── DIAGNOSTIC: jumping steps ────────────────────────────────────────────
if cont_pass < cont_total
    fprintf('=== DIAGNOSTIC: jumping steps ===\n');
    for i = 2:N
        if ~solve_ok(i)||~solve_ok(i-1), continue; end
        dq=q_traj(:,i)-q_traj(:,i-1);
        if norm(dq)>=dq_threshold
            fprintf('Step %d  |  q_prev: %s deg\n',i, ...
                mat2str(round(rad2deg(q_traj(:,i-1))')));
            fprintf('         |  q_curr: %s deg\n', ...
                mat2str(round(rad2deg(q_traj(:,i))')));
            for j=1:7
                fprintf('  J%d: %+.4f rad (%+.2f deg)\n', ...
                    j,dq(j),rad2deg(dq(j)));
            end
        end
    end
    fprintf('\n');
end

%% ── TEST 4: Joint limit compliance ───────────────────────────────────────
fprintf('=== TEST 4: Joint limit compliance ===\n');
fprintf('%-8s  %-10s  %-10s  %-10s  %-10s  %s\n', ...
    'Joint','min(deg)','max(deg)','limit_min','limit_max','Status');
fprintf('%s\n',repmat('-',1,62));
jl_pass=0; Q_solved=q_traj(:,solve_ok);
for j = 1:7
    q_j=Q_solved(j,:); jmin=min(q_j); jmax=max(q_j);
    if jmin>=q_min(j) && jmax<=q_max(j)
        status='PASS'; jl_pass=jl_pass+1;
    else
        status='FAIL';
    end
    fprintf('J%-7d  %-10.2f  %-10.2f  %-10.2f  %-10.2f  %s\n',j, ...
        rad2deg(jmin),rad2deg(jmax),rad2deg(q_min(j)),rad2deg(q_max(j)),status);
end
fprintf('Joint limits: %d/7 joints fully compliant\n',jl_pass);
if jl_pass==7, fprintf('Status: PASS\n\n');
else, fprintf('Status: FAIL\n\n'); end

%% ── TEST 5: Singularity avoidance ────────────────────────────────────────
fprintf('=== TEST 5: Singularity avoidance (sigma_min) ===\n');
fprintf('%-10s  %-12s  %-12s  %s\n','Waypoint','sigma_min','sigma_max','Status');
fprintf('%s\n',repmat('-',1,48));
sigma_threshold=0.05;
sing_pass=0; sing_total=0;
sigma_min_all=zeros(1,N);
sigma_min_global=inf; sigma_min_wp=0;
for i = 1:N
    if ~solve_ok(i), continue; end
    J=Geometric_Jacobian(q_traj(:,i),L); sv=svd(J);
    sm=sv(end); sM=sv(1);
    sigma_min_all(i)=sm; sing_total=sing_total+1;
    if sm<sigma_min_global, sigma_min_global=sm; sigma_min_wp=i; end
    if sm>sigma_threshold
        sing_pass=sing_pass+1;
    else
        fprintf('%-10d  %-12.4f  %-12.4f  FAIL\n',i,sm,sM);
    end
end
fprintf('Global sigma_min: %.4f at waypoint %d\n',sigma_min_global,sigma_min_wp);
fprintf('Singularity avoidance: %d/%d waypoints above threshold %.2f\n', ...
    sing_pass,sing_total,sigma_threshold);
if sing_pass==sing_total, fprintf('Status: PASS\n\n');
else, fprintf('Status: FAIL\n\n'); end

%% ── SUMMARY ──────────────────────────────────────────────────────────────
fprintf('=== SUMMARY ===\n');
fprintf('IK solvability:       %d/%d\n',n_solved,N);
fprintf('FK round-trip:        %d/%d\n',fk_pass,n_solved);
fprintf('Joint continuity:     %d/%d\n',cont_pass,cont_total);
fprintf('Joint limit (joints): %d/7\n',jl_pass);
fprintf('Singularity avoid:    %d/%d\n',sing_pass,sing_total);

%% ── Plot ─────────────────────────────────────────────────────────────────
figure('Name','Trajectory','Position',[100 100 1200 800]);
subplot(2,2,1);
plot3(waypoints_p(1,:),waypoints_p(2,:),waypoints_p(3,:),'b-','LineWidth',2);
hold on;
plot3(p_start(1),p_start(2),p_start(3),'go','MarkerSize',10,'LineWidth',2);
plot3(p_end(1),p_end(2),p_end(3),'rs','MarkerSize',10,'LineWidth',2);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Cartesian Trajectory'); grid on;
legend('Path','Start','End');

subplot(2,2,2);
plot(t_vals,rad2deg(q_traj)','LineWidth',1.5);
xlabel('t'); ylabel('Joint angle (deg)');
title('Joint Trajectories');
legend('J1','J2','J3','J4','J5','J6','J7','Location','best');
grid on;

subplot(2,2,3);
dq_norms=zeros(1,N-1);
for i=2:N
    if solve_ok(i)&&solve_ok(i-1)
        dq_norms(i-1)=norm(q_traj(:,i)-q_traj(:,i-1));
    end
end
plot(2:N,dq_norms,'k-','LineWidth',1.5);
yline(dq_threshold,'r--','Threshold');
xlabel('Waypoint'); ylabel('||dq|| (rad)');
title('Joint Step Size'); grid on;

subplot(2,2,4);
plot(1:N,sigma_min_all,'b-','LineWidth',1.5);
yline(sigma_threshold,'r--','Threshold');
xlabel('Waypoint'); ylabel('\sigma_{min}');
title('Minimum Singular Value Along Path'); grid on;

sgtitle('NERO 7 — Straight-Line Trajectory Verification');
