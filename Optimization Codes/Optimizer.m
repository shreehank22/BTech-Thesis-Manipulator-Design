clear; clc;
rng(42,'twister')

if isempty(gcp('nocreate'))
    parpool('threads');
end

joint_samples     = generate_joint_samples(5000);
workspace_samples = joint_samples;

% 6 variables: [l2, l3, l4, l5, l6, l7]
% l1 = l2 always by construction
lb = [0.10, 0.10, 0.10, 0.10, 0.10, 0.10];
ub = [1.00, 1.00, 1.00, 1.00, 1.00, 1.00];

% Expand 6-vector to 7-vector: l1 = l2 = x(1)
expand = @(x) [x(1), x(1), x(2), x(3), x(4), x(5), x(6)];

% Objective over 6 variables
objfun = @(x) Objective(expand(x), joint_samples, workspace_samples);

% Constraints over 6 variables — no ceq needed
confun = @(x) deal(constraints(expand(x)), []);

% Check feasibility of known good solution
L3 = [0.1085, 0.5115, 0.8982, 0.7683, 0.4549, 0.1524, 0.1062];
x3 = L3(2:7);   % drop l1 since l1=l2
c_check = constraints(expand(x3));
fprintf('L3 feasibility check — max(c) = %.6f (must <= 0)\n', max(c_check))

options_nsga = optimoptions('gamultiobj', ...
    'PopulationSize',     600, ...
    'MaxGenerations',     200, ...
    'CrossoverFraction',  0.8, ...
    'MutationFcn',        @mutationadaptfeasible, ...
    'SelectionFcn',       {@selectiontournament, 4}, ...
    'ConstraintTolerance',1e-6, ...
    'UseParallel',        true, ...
    'Display',            'iter', ...
    'PlotFcn',            {@gaplotpareto});

fprintf('\nRunning NSGA-II...\n')
tic
[x6_pareto, f_pareto] = gamultiobj( ...
    objfun, 6, [], [], [], [], lb, ub, confun, options_nsga);
t_nsga = toc;

% Expand Pareto solutions back to 7 variables
x_pareto = [x6_pareto(:,1), x6_pareto];
fprintf('\nNSGA-II finished in %.2f seconds\n', t_nsga)
fprintf('Pareto front size: %d solutions\n', size(x_pareto,1))