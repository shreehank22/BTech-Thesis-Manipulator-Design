function tau = dynamic_Torque(q, qd, qdd, L, radii, rho, m_act, F_payload)
    if nargin < 8
        F_payload = zeros(3,1);
    end
    F_payload = F_payload(:);
    N  = length(q);
    g0 = [0; 0; -9.81];
    m_links = mass(L, radii, rho);
    I_links = inertia_tensor(L, m_links, radii{1});
    m_act_p = [m_act(:); 0];
    [~, ~, T_all] = Forward_Kinematics(q, L);
    R = cell(N+1, 1);
    p = cell(N+1, 1);
    z = cell(N+1, 1);
    r_com = cell(N, 1);
    R{1} = eye(3);
    p{1} = zeros(3,1);
    z{1} = [0;0;1];
    for i = 1:N
        R{i+1}   = T_all{i}(1:3,1:3);
        p{i+1}   = T_all{i}(1:3,4);
        z{i+1}   = T_all{i}(1:3,3);
        r_com{i} = 0.5*(p{i} + p{i+1});
    end
    omega = cell(N+1, 1);
    omegad = cell(N+1, 1);
    acc = cell(N+1, 1);
    acc_com = cell(N, 1);
    omega{1} = zeros(3,1);
    omegad{1} = zeros(3,1);
    acc{1} = g0;
    for i = 1:N
        zi = z{i};
        omega{i+1}  = omega{i}  + qd(i)  * zi;
        omegad{i+1} = omegad{i} + qdd(i) * zi ...
                    + cross(omega{i}, qd(i) * zi);
        r_i      = p{i+1} - p{i};
        acc{i+1} = acc{i} ...
                 + cross(omegad{i+1}, r_i) ...
                 + cross(omega{i+1}, cross(omega{i+1}, r_i));
        r_ci       = r_com{i} - p{i};
        acc_com{i} = acc{i} ...
                   + cross(omegad{i+1}, r_ci) ...
                   + cross(omega{i+1}, cross(omega{i+1}, r_ci));
    end
    F_out = cell(N+1, 1);
    M_out = cell(N+1, 1);
    F_out{N+1} = F_payload;        
    M_out{N+1} = zeros(3,1);      
    tau = zeros(N, 1);
    for i = N:-1:1
        mi = m_links(i);
        Ii = I_links(:,:,i);
        Ri = R{i+1};
        Ii_w = Ri * Ii * Ri';
        F_net = mi * acc_com{i};
        F_act = m_act_p(i+1) * acc{i+1};
        F_out{i} = F_net + F_act + F_out{i+1};
        alpha_i = omegad{i+1};
        omega_i = omega{i+1};
        tau_inertia = Ii_w * alpha_i + cross(omega_i, Ii_w * omega_i);
        r_ci   = r_com{i} - p{i};
        r_next = p{i+1}   - p{i};
        M_out{i} = tau_inertia ...
                 + M_out{i+1} ...
                 + cross(r_next, F_out{i+1} + F_act) ...
                 + cross(r_ci,   F_net);
        tau(i) = dot(M_out{i}, z{i});
    end
end