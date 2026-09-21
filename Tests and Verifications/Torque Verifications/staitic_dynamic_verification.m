clc; clear;

L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
radii  = repmat({[0.06, 0.057]}, 1, 7);
rho    = 2700;
m_act  = [2.2; 2.75; 1.54; 1.54; 0.715; 0.715; 0.715];
m_links = mass(L, radii, rho);

F_lat  = [0;  50;   0];   % 50 N lateral Y
F_down = [0;   0; -50];   % 50 N downward Z

%% --- Worst-case configurations ---
configs = {
    'J1 worst',  [0;  0;   0;   0;  0;  0;  0]*pi/180,  F_lat;
    'J2 worst',  [0;  0;   0;   0;  0;  0;  0]*pi/180,  F_down;
    'J4 worst',  [0;  0;   0;   0;  0;  0;  0]*pi/180,  F_down;
    'J7 worst',  [0;  0;   0;   0;  0;  0;  0]*pi/180,  F_down;
    'J3 worst',  [0;  0;  90;  90;  0;  0;  0]*pi/180,  F_down;
    'J5 worst',  [0;  0; -90; -90;  0;  0; 90]*pi/180,  F_down;
    'J6 worst',  [0;  0;   0;   0; 90;  0;  0]*pi/180,  F_down;
};

%% --- Static torque: gravity_Torque + Jacobian payload ---
fprintf('%s\n', repmat('=',1,100));
fprintf('  Static Torque Verification — gravity_Torque + Jacobian payload\n');
fprintf('%s\n', repmat('=',1,100));
fprintf('  %-12s  %-7s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Config','Force','J1','J2','J3','J4','J5','J6','J7');
fprintf('  %s\n', repmat('-',1,95));

for c = 1:size(configs,1)
    label = configs{c,1};
    q     = configs{c,2};
    F     = configs{c,3};

    tau_g     = gravity_Torque(q, m_links, L, m_act);
    J         = Geometric_Jacobian(q, L);
    tau_p     = J(1:3,:)' * F;
    tau_total = tau_g + tau_p;

    if F(2) ~= 0, flabel = '+Y 50N'; else, flabel = '-Z 50N'; end

    fprintf('  %-12s  %-7s  %-10.3f  %-10.3f  %-10.3f  %-10.3f  %-10.3f  %-10.3f  %-10.3f\n', ...
        label, flabel, tau_total(1), tau_total(2), tau_total(3), ...
        tau_total(4), tau_total(5), tau_total(6), tau_total(7));
end
fprintf('%s\n\n', repmat('=',1,100));

%% --- RNEA at zero qd and qdd ---
fprintf('%s\n', repmat('=',1,100));
fprintf('  RNEA Verification — dynamic_Torque(q, 0, 0) + Jacobian payload\n');
fprintf('%s\n', repmat('=',1,100));
fprintf('  %-12s  %-7s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Config','Force','J1','J2','J3','J4','J5','J6','J7');
fprintf('  %s\n', repmat('-',1,95));

for c = 1:size(configs,1)
    label = configs{c,1};
    q     = configs{c,2};
    F     = configs{c,3};

    tau_rnea  = dynamic_Torque(q, zeros(7,1), zeros(7,1), L, radii, rho, m_act);
    J         = Geometric_Jacobian(q, L);
    tau_p     = J(1:3,:)' * F;
    tau_total = tau_rnea + tau_p;

    if F(2) ~= 0, flabel = '+Y 50N'; else, flabel = '-Z 50N'; end

    fprintf('  %-12s  %-7s  %-10.3f  %-10.3f  %-10.3f  %-10.3f  %-10.3f  %-10.3f  %-10.3f\n', ...
        label, flabel, tau_total(1), tau_total(2), tau_total(3), ...
        tau_total(4), tau_total(5), tau_total(6), tau_total(7));
end
fprintf('%s\n\n', repmat('=',1,100));

%% --- Residual: gravity_Torque vs RNEA ---
fprintf('%s\n', repmat('=',1,100));
fprintf('  Residual — gravity_Torque vs RNEA (should be zero)\n');
fprintf('%s\n', repmat('=',1,100));
fprintf('  %-12s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-8s\n', ...
    'Config','J1','J2','J3','J4','J5','J6','J7','Status');
fprintf('  %s\n', repmat('-',1,95));

for c = 1:size(configs,1)
    label = configs{c,1};
    q     = configs{c,2};

    tau_g    = gravity_Torque(q, m_links, L, m_act);
    tau_rnea = dynamic_Torque(q, zeros(7,1), zeros(7,1), L, radii, rho, m_act);
    residual = tau_g - tau_rnea;

    if max(abs(residual)) < 1e-6
        status = 'PASS';
    else
        status = 'FAIL';
    end

    fprintf('  %-12s  %-10.2e  %-10.2e  %-10.2e  %-10.2e  %-10.2e  %-10.2e  %-10.2e  %-8s\n', ...
        label, residual(1), residual(2), residual(3), residual(4), ...
        residual(5), residual(6), residual(7), status);
end
fprintf('%s\n\n', repmat('=',1,100));