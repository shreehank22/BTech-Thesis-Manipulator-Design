function [GMI, GCI, SLI] = compute_metrics(x, joint_samples)

GMI = compute_GMI(x, joint_samples);

GCI = compute_GCI(joint_samples, x);


SLI = compute_SLI(x,joint_samples);

% Safety
if any(~isfinite([GMI GCI SLI]))
    GMI = 0; GCI = 0; SLI = inf;
end

end