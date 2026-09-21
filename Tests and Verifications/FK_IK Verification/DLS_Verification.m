%% DLS_Verification.m
clear; clc;

L = [0.200, 0.200, 1.000, 0.7424, 0.4635, 0.113, 0.106]
radii = repmat({[0.06, 0.057]}, 1, 7);
rho = 2700;
m_links = mass(L,radii,rho);

configs = {
    'General 1',   [30,  45, -20,  60,  15, -30,  10];
    'General 2',   [-45, 30,  90,  90, -15,  20,  -5];
    'J3 worst',    [0,    0,  90,  90,   0,   0,   0];
    'J5 worst',    [0,    0, -90, -90,   0,   0,  90];
    'Zero joints', [0,    0,   0,  60,   0,   0,   0];
    'Random 1',    [80,  -40,  60, 100, -30,  45, -20];
};

n_iter = 10;
rng(42);

fprintf('=== FK-IK Round-Trip Verification (10 iterations) ===\n\n');

all_means = zeros(size(configs,1),1);
all_max   = zeros(size(configs,1),1);

for c = 1:size(configs,1)
    fprintf('--- %s ---\n', configs{c,1});

    % Step 1: q_test -> FK -> T_home
    q_test = configs{c,2}' * pi/180;
    [T_home,~,~] = Forward_Kinematics(q_test, L);
    p_home = T_home(1:3,4);
    R_home = T_home(1:3,1:3);

    errors = zeros(n_iter,1);

    for i = 1:n_iter
        % Fresh random perturbation each iteration — this is the key fix
        q_init = q_test + (5*pi/180)*randn(7,1);

        % IK from perturbed init toward home pose
        [q_sol,~,~,~,~] = DLS(q_init, L, p_home, R_home, 1);

        % FK on IK solution
        [T_sol,~,~] = Forward_Kinematics(q_sol, L);
        p_sol = T_sol(1:3,4);
        R_sol = T_sol(1:3,1:3);

        % Position and orientation error
        pos_err = norm(p_sol - p_home);
        rot_err = norm(logSO3(R_sol * R_home'));

        errors(i) = pos_err;
        fprintf('Iteration %2d — pos_err: %.4e m   rot_err: %.4e rad\n', ...
            i, pos_err, rot_err);
    end

    all_means(c) = mean(errors);
    all_max(c)   = max(errors);

    fprintf('\n--- Round-Trip Statistics (%d iterations) ---\n', n_iter);
    fprintf('Mean Error : %.4e m\n', mean(errors));
    fprintf('Std Dev    : %.4e m\n', std(errors));
    fprintf('Min Error  : %.4e m\n', min(errors));
    fprintf('Max Error  : %.4e m\n\n', max(errors));
end

fprintf('=== Overall Summary ===\n');
fprintf('%-14s  %-12s  %-12s\n','Config','Mean err(m)','Max err(m)');
fprintf('%s\n', repmat('-',1,42));
for c = 1:size(configs,1)
    fprintf('%-14s  %-12.4e  %-12.4e\n', ...
        configs{c,1}, all_means(c), all_max(c));
end
