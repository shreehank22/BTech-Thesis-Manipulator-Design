function [link_com, overall_com, p_origins, z_axis, T_all] = COM_RST(q, L, masses, robot)
    N = length(L);
    body_names = robot.BodyNames;
    T_all = cell(1, N);
    for i = 1:N
        T_all{i} = getTransform(robot, q, body_names{i}, 'base');
    end
    p_origins    = cell(1, N+1);
    p_origins{1} = zeros(3,1); 
    for i = 1:N
        p_origins{i+1} = T_all{i}(1:3,4);
    end
    z_axis    = cell(1, N);
    z_axis{1} = [0;0;1]; 
    for i = 2:N
        z_axis{i} = T_all{i-1}(1:3,3);
    end
    link_com_cell = cell(1, N);
    for i = 1:N
        p_prox = p_origins{i};
        p_dist = p_origins{i+1};
        if norm(p_dist - p_prox) > 1e-6
            link_com_cell{i} = 0.5 * (p_prox + p_dist);
        else
            if i == 1
                R_prev = eye(3);
                T_local = T_all{1};
            else
                R_prev = T_all{i-1}(1:3,1:3);
                T_local = T_all{i-1} \ T_all{i};
            end
            local_disp = T_local(1:3,4);
            world_disp = R_prev * local_disp;
            link_com_cell{i} = p_prox + 0.5 * world_disp;
        end
    end
    link_com    = [link_com_cell{:}];
    overall_com = (link_com * masses(:)) / sum(masses);
end