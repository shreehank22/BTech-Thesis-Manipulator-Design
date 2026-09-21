function tau_g = gravity_Torque(q, m, l, m_act)
    gravity = [0; 0; -9.81];
    tau_g = zeros(length(q), 1);
    [link_com, ~, p_origins, z_axis] = linkCOMPositions(q, l, m);
    for i = 1:length(q)
        for j = i:length(q)
            r_ij = link_com(:,j) - p_origins{i};
            F = m(j) * gravity;
            tau_g(i) = tau_g(i) + dot(cross(r_ij, F), z_axis{i});
        end
        for j = i+1:length(q)
            r_ij     = p_origins{j} - p_origins{i}; 
            F_act    = m_act(j) * gravity;
            tau_g(i) = tau_g(i) + dot(cross(r_ij, F_act), z_axis{i});
        end
    end
end