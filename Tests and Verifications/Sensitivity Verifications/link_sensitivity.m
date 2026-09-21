function sensitivity = link_sensitivity(x_best, joint_samples)

N = length(x_best);

% Reference metrics
[GMI_ref, GCI_ref, DWI_ref, SLI_ref] = compute_metrics(x_best, joint_samples);

% Tolerances
tol_GMI = 0.05;
tol_GCI = 0.05;
tol_DWI = 0.05;
tol_SLI = 0.10;

sensitivity = zeros(N,2);

for i = 1:N

    delta_vals = linspace(-0.1, 0.1, 80);
    valid = [];

    for d = delta_vals

        x_test = x_best;
        x_test(i) = x_test(i) + d;

        % Check constraints
        [c,~] = constraints(x_test);
        if any(c > 0)
            continue
        end

        % Compute metrics
        [GMI, GCI, DWI, SLI] = compute_metrics(x_test, joint_samples);

        % Check ALL conditions
        if (GMI >= (1 - tol_GMI)*GMI_ref) && ...
           (GCI >= (1 - tol_GCI)*GCI_ref) && ...
           (DWI >= (1 - tol_DWI)*DWI_ref) && ...
           (SLI <= (1 + tol_SLI)*SLI_ref)

            valid = [valid d];
        end

    end

    if ~isempty(valid)
        sensitivity(i,1) = min(valid);
        sensitivity(i,2) = max(valid);
    end

end
end