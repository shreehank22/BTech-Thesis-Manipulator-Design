L = [0.1040, 0.1040, 0.9993, 0.9072, 0.6779, 0.1053, 0.1013];
N = 100000;
opts.elbow_exclusion = true;
opts.joint_limits = true;
Q = generate_joint_samples(N, L, opts);  

%% Global Conditioning Index analysis for Wrist configurations

% R-Y-P
LCI_vals = zeros(N, 1);
for i = 1:N
    q = Q(i,:)';
    J = Geometric_Jacobian(q, L);   
    J_wrist = J(4:6, 5:7);               
    sigma = svd(J_wrist);
    if max(sigma) < eps
        LCI_vals(i) = 0;
    else
        LCI_vals(i) = min(sigma) / max(sigma);
    end
end

fprintf('\n=== Wrist GCI — R-Y-P ===\n');
fprintf('GCI = %.4f\n', mean(LCI_vals));
fprintf('Std LCI = %.4f\n', std(LCI_vals));
fprintf('Min LCI = %.6f\n', min(LCI_vals));
fprintf('Max LCI = %.6f\n', max(LCI_vals));
fprintf('Near-singular = %.2f%%\n', 100*mean(LCI_vals < 0.01));

% R-P-Y
LCI_vals_rpy = zeros(N, 1);
for i = 1:N
    q = Q(i,:)';
    J = Geometric_Jacobian_RPY(q, L);   
    J_wrist = J(4:6, 5:7);               
    sigma = svd(J_wrist);
    if max(sigma) < eps
        LCI_vals_rpy(i) = 0;
    else
        LCI_vals_rpy(i) = min(sigma) / max(sigma);
    end
end

fprintf('\n=== Wrist GCI — R-P-Y ===\n');
fprintf('GCI = %.4f\n', mean(LCI_vals_rpy));
fprintf('Std LCI = %.4f\n', std(LCI_vals_rpy));
fprintf('Min LCI = %.6f\n', min(LCI_vals_rpy));
fprintf('Max LCI = %.6f\n', max(LCI_vals_rpy));
fprintf('Near-singular = %.2f%%\n', 100*mean(LCI_vals_rpy < 0.01));

MI_ryp = zeros(N,1);
for i = 1:N
    q = Q(i,:)';
    J = Geometric_Jacobian(q, L);           % R-Y-P topology
    J_wrist = J(4:6, 5:7);
    MI_ryp(i) = sqrt(det(J_wrist*J_wrist'));
end

MI_rpy = zeros(N,1);
for i = 1:N
    q = Q(i,:)';
    J = Geometric_Jacobian_RPY(q, L);       % R-P-Y topology
    J_wrist = J(4:6, 5:7);
    MI_rpy(i) = sqrt(det(J_wrist*J_wrist'));
end

%% Results
fprintf('\n=== Wrist MI Comparison ===\n');
fprintf('%-20s  %8s  %8s  %8s  %8s  %12s\n', ...
    'Topology', 'GMI', 'Std', 'Min', 'Max', 'Near-sing');
fprintf('%s\n', repmat('-',1,68));
fprintf('%-20s  %8.4f  %8.4f  %8.6f  %8.6f  %11.2f%%\n', ...
    'R-Y-P (current)', mean(MI_ryp), std(MI_ryp), min(MI_ryp), max(MI_ryp), 100*mean(MI_ryp<0.01));
fprintf('%-20s  %8.4f  %8.4f  %8.6f  %8.6f  %11.2f%%\n', ...
    'R-P-Y (alt)',     mean(MI_rpy), std(MI_rpy), min(MI_rpy), max(MI_rpy), 100*mean(MI_rpy<0.01));