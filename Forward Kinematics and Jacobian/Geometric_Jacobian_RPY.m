function J = Geometric_Jacobian_RPY(theta, link_lengths)
J = zeros(6, 7); % Initialize the Jacobian matrix

[~, ~, T_all] = Forward_Kinematics_RPY(theta, link_lengths);
ee_pose = T_all{7}(1:3, 4);

% World frame
z_0 = [0; 0; 1];
p_0 = [0; 0; 0];

for i = 1:7
    if i == 1
        z_prev = z_0;
        p_prev = p_0;
    else
        z_prev = T_all{i-1}(1:3,3);
        p_prev = T_all{i-1}(1:3, 4);
    end
    r = ee_pose - p_prev;
    J(1:3, i) = cross(z_prev, r);
    J(4:6, i) = z_prev;
end
end