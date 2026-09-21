function results = monte_carlo_tolerance(x_best, joint_samples)

    N_mc     = 10000;   % number of Monte Carlo trials
    sigma    = 0.005;   % initial std dev for sampling (5 mm)
    tol_GMI  = 0.02;    % 2% degradation allowed
    tol_GCI  = 0.02;
    tol_DWI  = 0.02;
    tol_SLI  = 0.05;    % 5% increase allowed

    % Reference metrics at x_best
    [GMI_ref, GCI_ref, DWI_ref, SLI_ref] = compute_metrics(x_best, joint_samples);
    fprintf('\nReference metrics:\n');
    fprintf('  GMI=%.4f | GCI=%.4f | DWI=%.4f | SLI=%.4f\n\n', ...
        GMI_ref, GCI_ref, DWI_ref, SLI_ref);

    % ── Sample tolerances from truncated normal distribution ─────────────
    % Each trial: sample a 7-vector of perturbations ~ N(0, sigma^2)
    % Truncated to [-3sigma, +3sigma] to avoid extreme outliers
    n_links    = length(x_best);
    raw        = sigma * randn(N_mc, n_links);
    raw        = max(raw, -3*sigma);
    raw        = min(raw,  3*sigma);

    % ── Evaluate each trial ───────────────────────────────────────────────
    feasible_deltas = [];
    feasible_pct    = [];

    fprintf('Running %d Monte Carlo trials...\n', N_mc);
    for k = 1:N_mc
        delta   = raw(k, :);
        x_test  = x_best + delta;

        % 1. Hard constraint check
        [c, ~] = constraints(x_test);
        if any(c > 0)
            continue
        end

        % 2. Metric evaluation
        [GMI, GCI, DWI, SLI] = compute_metrics(x_test, joint_samples);

        pct_GMI = (GMI - GMI_ref) / GMI_ref * 100;
        pct_GCI = (GCI - GCI_ref) / GCI_ref * 100;
        pct_DWI = (DWI - DWI_ref) / DWI_ref * 100;
        pct_SLI = (SLI - SLI_ref) / SLI_ref * 100;

        % 3. Metric tolerance check
        if (pct_GMI >= -tol_GMI*100) && ...
           (pct_GCI >= -tol_GCI*100) && ...
           (pct_DWI >= -tol_DWI*100) && ...
           (pct_SLI <=  tol_SLI*100)

            feasible_deltas(end+1, :) = delta;
            feasible_pct(end+1, :)    = [pct_GMI, pct_GCI, pct_DWI, pct_SLI];
        end
    end

    n_feasible = size(feasible_deltas, 1);
    fprintf('Feasible trials: %d / %d  (%.1f%%)\n\n', ...
        n_feasible, N_mc, 100*n_feasible/N_mc);

    if n_feasible == 0
        fprintf('No feasible tolerance vectors found. Try reducing sigma.\n');
        results = [];
        return
    end

    % ── Derive recommended tolerance per link ────────────────────────────
    % For each link: find the 5th/95th percentile of feasible deltas
    % This gives a conservative tolerance interval
    tol_low  = prctile(feasible_deltas, 5,  1);   % 5th  percentile per link
    tol_high = prctile(feasible_deltas, 95, 1);   % 95th percentile per link
    tol_sym  = min(abs(tol_low), abs(tol_high));  % symmetric conservative bound

    % ── Find the single best tolerance vector ────────────────────────────
    % Best = feasible sample with smallest total perturbation norm
    % AND best metric preservation (weighted score)
    w = [1.0, 0.5, 0.5, 0.3];   % weights: GMI most important
    scores = zeros(n_feasible, 1);
    for k = 1:n_feasible
        metric_score = -sum(w .* abs(feasible_pct(k,:)));   % less change = higher score
        norm_score   = -norm(feasible_deltas(k,:));          % smaller delta = higher score
        scores(k)    = metric_score + norm_score;
    end
    [~, best_idx] = max(scores);
    best_delta    = feasible_deltas(best_idx, :);
    best_pct      = feasible_pct(best_idx, :);

    % ── Display results ───────────────────────────────────────────────────
    fprintf('Per-link tolerance bounds (5th-95th percentile of feasible samples):\n\n');
    fprintf('%-6s  %10s  %10s  %10s\n', 'Link', 'Lower(m)', 'Upper(m)', 'Symmetric(m)');
    fprintf('%s\n', repmat('-', 1, 44));
    for i = 1:n_links
        fprintf('l%-5d  %+10.4f  %+10.4f  %10.4f\n', i, ...
            tol_low(i), tol_high(i), tol_sym(i));
    end

    fprintf('\nBest tolerance vector (min perturbation + max metric preservation):\n\n');
    fprintf('%-6s  %10s  %8s  %8s  %8s  %8s\n', ...
        'Link', 'Delta(m)', 'dGMI', 'dGCI', 'dDWI', 'dSLI');
    fprintf('%s\n', repmat('-', 1, 58));
    for i = 1:n_links
        fprintf('l%-5d  %+10.4f  %7.2f%%  %7.2f%%  %7.2f%%  %7.2f%%\n', i, ...
            best_delta(i), best_pct(1), best_pct(2), best_pct(3), best_pct(4));
    end

    fprintf('\nMetric change at best tolerance vector:\n');
    fprintf('  dGMI=%+.2f%%  dGCI=%+.2f%%  dDWI=%+.2f%%  dSLI=%+.2f%%\n', ...
        best_pct(1), best_pct(2), best_pct(3), best_pct(4));

    % ── Pack results ─────────────────────────────────────────────────────
    results.feasible_deltas = feasible_deltas;
    results.feasible_pct    = feasible_pct;
    results.tol_low         = tol_low;
    results.tol_high        = tol_high;
    results.tol_sym         = tol_sym;
    results.best_delta      = best_delta;
    results.best_pct        = best_pct;
    results.acceptance_rate = n_feasible / N_mc;

    % ── Plot distribution of feasible deltas per link ────────────────────
    plot_tolerance_distributions(feasible_deltas, tol_low, tol_high, n_links);

end

% -------------------------------------------------------------------------
function plot_tolerance_distributions(feasible_deltas, tol_low, tol_high, n_links)

    figure('Name', 'Feasible Tolerance Distributions', ...
           'Position', [100 100 1100 500]);

    for i = 1:n_links
        subplot(2, 4, i);
        histogram(feasible_deltas(:,i)*1000, 30, ...
            'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        hold on;
        xline(tol_low(i)*1000,  'r--', 'LineWidth', 1.5);
        xline(tol_high(i)*1000, 'r--', 'LineWidth', 1.5);
        xline(0, 'k-', 'LineWidth', 1);
        xlabel('\Delta l (mm)');
        ylabel('Count');
        title(sprintf('l_%d  [%.1f, %.1f] mm', i, tol_low(i)*1000, tol_high(i)*1000));
        grid on; box off;
    end

    subplot(2, 4, 8);
    axis off;
    text(0.1, 0.6, 'Red dashed = 5th/95th', 'FontSize', 10);
    text(0.1, 0.4, 'percentile bounds',      'FontSize', 10);
    text(0.1, 0.2, 'Black = nominal (0)',    'FontSize', 10);

    sgtitle('Monte Carlo tolerance distributions — feasible samples only');

end

% -------------------------------------------------------------------------
function [GMI, GCI, DWI, SLI] = compute_metrics(x, joint_samples)

    GMI = compute_GMI(x, joint_samples);
    GMI = GMI / (GMI + 1);
    GCI = compute_GCI(joint_samples, x);
    DWI = compute_DWI(x, joint_samples);
    SLI = compute_SLI(x, joint_samples);

    if any(~isfinite([GMI GCI DWI SLI]))
        GMI = 0; GCI = 0; DWI = 0; SLI = inf;
    end

end