function tau_total = payloadTorque_vec(L, q, radii, rho, F_vec,m_act)
    if iscell(q)
        q = cell2mat(q);
    end
    L = double(L(:));
    q  = double(q(:));
    F_vec = double(F_vec(:));
    N  = length(L);
    tau_p = zeros(N, 1);
    [~, ~, T_all] = Forward_Kinematics(q, L);
    ee_pos = T_all{N}(1:3, 4);
    for i = 1:N
        if i == 1
            p_i = zeros(3,1);
            z_i = [0;0;1];
        else
            p_i = T_all{i-1}(1:3,4);
            z_i = T_all{i-1}(1:3,3);
        end
        r        = ee_pos - p_i;
        tau_p(i) = dot(cross(r, F_vec), z_i);
    end
    m_links = mass(L, radii, rho);
    tau_g = gravity_Torque(q, m_links,L,m_act);
    tau_total = tau_p + tau_g;
end