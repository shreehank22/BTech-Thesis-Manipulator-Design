function [q_sol, converged,n_iter,pos_err,rot_err] = DLS(theta_curr, L, p_des, R_des,k)
    if nargin < 5, k = 1; end
    lambda_max = 0.1;
    epsilon  = 0.05;
    max_iter = 200;
    converged  = false;
    for iter = 1:max_iter
        n_iter=iter;
        [T_curr,~,~] = Forward_Kinematics(theta_curr, L);
        p_curr = T_curr(1:3,4);
        R_curr = T_curr(1:3,1:3);
        e_P = p_des - p_curr;
        e_R = logSO3(R_des * R_curr');
        e = [e_P; e_R];
        pos_err = norm(e_P);
        rot_err = norm(e_R);
        if norm(e) < 1e-7
            converged = true;
            break;
        end
        J = Geometric_Jacobian(theta_curr, L);
        [U,S,V] = svd(J);
        sigma = diag(S);
        lambda = lambda_max * exp(-(sigma/epsilon).^2);
        S_inv = zeros(7,6);
        for i = 1:6
            S_inv(i,i) = sigma(i) / (sigma(i)^2 + lambda(i)^2);
        end
        delta_theta = V * S_inv * U' * e;
        step_limit = 0.5;  
        if norm(k * delta_theta) > step_limit
            delta_theta = delta_theta * step_limit / norm(k * delta_theta);
        else
            delta_theta = k * delta_theta;
        end
        theta_curr = theta_curr + k*delta_theta;
    end
q_sol = theta_curr;
end