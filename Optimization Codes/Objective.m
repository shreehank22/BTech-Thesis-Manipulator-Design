function f = Objective(x, joint_samples, workspace_samples)
GMI = compute_GMI(x, joint_samples);
GCI = compute_GCI(joint_samples,x);
SLI = compute_SLI(x,workspace_samples);
if any(~isfinite([GCI GMI SLI]))
    f = [1e6,1e6,1e6];
return
end
f = [-GCI,-GMI,SLI];
end