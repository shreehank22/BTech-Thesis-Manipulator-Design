function Qsamples = generate_joint_samples(N, link_lengths, opts)

    if nargin < 2 || isempty(link_lengths)
        link_lengths = [];
    end
    if nargin < 3
        opts = struct();
    end
    deg = pi/180;
    q_min_default = [-170, -80, -80, -170, -80, -75, -170] * pi/180;
    q_max_default = [ 170,  80,  80,  170,  80,  75,  170] * pi/180;
    q_min = double(get_option(opts, 'q_min', q_min_default));
    q_max = double(get_option(opts, 'q_max', q_max_default));
    verbose              = get_option(opts, 'verbose',                  false);
    use_joint_limits     = get_option(opts, 'joint_limits',             true);
    use_self_collision   = get_option(opts, 'self_collision',           false);
    use_workspace_bounds = get_option(opts, 'workspace_bounds',         false);
    use_singularity      = get_option(opts, 'singularity',              false);
    use_elbow_excl       = get_option(opts, 'elbow_exclusion',          true);
    sigma_thresh         = get_option(opts, 'singularity_threshold',    0.05);
    collision_clearance  = get_option(opts, 'self_collision_clearance', 0.06);
    max_attempt_factor   = get_option(opts, 'max_attempt_factor',       100);
    elbow_excl_rad       = get_option(opts, 'elbow_exclusion_deg',      15) * deg;
    if numel(q_min) ~= 7 || numel(q_max) ~= 7
        error('q_min and q_max must each have 7 elements.');
    end
    q_min = double(q_min(:).');   
    q_max = double(q_max(:).');   
    Qsamples     = zeros(N, 7);   
    accepted     = 0;
    attempts     = 0;
    max_attempts = max_attempt_factor * N;
    rej_limits      = 0;
    rej_elbow       = 0;
    rej_collision   = 0;
    rej_workspace   = 0;
    rej_singularity = 0;

    while accepted < N && attempts < max_attempts
        attempts = attempts + 1;
        q = q_min + (q_max - q_min) .* rand(1, 7);   
        if use_joint_limits && any(q < q_min | q > q_max)
            rej_limits = rej_limits + 1;
            continue;
        end
        if use_elbow_excl && abs(q(4)) < elbow_excl_rad
            rej_elbow = rej_elbow + 1;
            continue;
        end
        if ~isempty(link_lengths)
            q_col = double(q(:)); 
            [~, ee_pose, T_all] = Forward_Kinematics(q_col, link_lengths);
            if use_workspace_bounds && ee_pose{7}(3) < 0
                rej_workspace = rej_workspace + 1;
                continue;
            end
            if use_self_collision && has_self_collision(T_all, collision_clearance)
                rej_collision = rej_collision + 1;
                continue;
            end
            if use_singularity
                J         = Geometric_Jacobian(q_col, link_lengths);
                sigma_min = min(svd(J));
                if sigma_min < sigma_thresh
                    rej_singularity = rej_singularity + 1;
                    continue;
                end
            end
        end
        accepted = accepted + 1;
        Qsamples(accepted, :) = q;   
    end

    Qsamples = double(Qsamples(1:accepted, :)); 

    if accepted < N
        warning('generate_joint_samples:InsufficientSamples', ...
            'Only %d/%d samples generated after %d attempts.', ...
            accepted, N, attempts);
    end
    if verbose
        fprintf('generate_joint_samples: %d/%d accepted after %d attempts\n', ...
                accepted, N, attempts);
        fprintf('  Acceptance rate      : %.1f%%\n', 100*accepted/attempts);
        fprintf('  Rejected — limits    : %d\n',     rej_limits);
        fprintf('  Rejected — elbow     : %d  (|q4| < %.0f deg)\n', ...
                rej_elbow, elbow_excl_rad/deg);
        fprintf('  Rejected — workspace : %d\n',     rej_workspace);
        fprintf('  Rejected — collision : %d\n',     rej_collision);
        fprintf('  Rejected — SVD       : %d\n',     rej_singularity);
    end
end

function value = get_option(opts, name, default_value)
    if isfield(opts, name)
        value = opts.(name);
    else
        value = default_value;
    end
end

function tf = has_self_collision(T_all, clearance)
    p      = zeros(8, 3);
    p(1,:) = [0, 0, 0];
    for k  = 1:7
        p(k+1,:) = T_all{k}(1:3, 4).';
    end
    tf = false;
    for i = 1:6
        for j = i+2:8
            if norm(p(i,:) - p(j,:)) < clearance
                tf = true;
                return;
            end
        end
    end
end