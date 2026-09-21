function sli = compute_SLI(link_lengths, Qworkspace)
    L_total = sum(link_lengths);
    [V_workspace] = compute_workspace_volume(link_lengths, Qworkspace);
    if V_workspace <= 0
        sli = inf;
    else
        sli = L_total / (V_workspace^(1/3));
    end
end