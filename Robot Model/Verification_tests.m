%% Loading robot parameters

q_reg = randomConfiguration(robot);
q = q_reg+dh(:,4); % RST is unable to consider offsets for 
L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
m_links = [0.3096,0.3096,2.9752,2.7010,2.0183,0.3135,0.3016];
m_act = [2.2000; 2.7500; 1.5400; 1.5400; 0.7150; 0.7150; 0.7150];
T = getTransform(robot,q,"end_effector","base");
Tf = Forward_Kinematics(q_reg,L);
diff = norm(T-Tf);
tolerance = 1e-6;
if diff < tolerance
    fprintf('[SUCCESS] Forward kinematics verified\n');
    fprintf('Error norm: %.3e\n', diff);
else
    fprintf('[WARNING] Forward kinematics mismatch ✖\n');
    fprintf('Error norm: %.3e (exceeds tolerance %.1e)\n', diff, tolerance);
end
%% Torque Verification
qd = ones(7,1);
qdd = ones(7,1)*2;

% Gravity Torque Verification
tau_rst = gravity_Torque_RST(q,m_links,L,m_act,robot);
tau_g = gravity_Torque(q_reg,m_links,L,m_act);
diff = tau_rst - tau_g;
disp(norm(diff));

% Dynamic Torque Verification
F = [0,0,-50];
tau_d_rst = dynamic_Torque_RST(q,qd,qdd,L,radii,rho,m_act,F,robot);
tau_d_g = dynamic_Torque(q_reg, qd, qdd, L, radii, rho, m_act, F);
diff_dynamic = tau_d_rst - tau_d_g;
disp(norm(diff_dynamic));

%% Jacobian verification
J = geometricJacobian(robot,q,"end_effector");
J_c = Geometric_Jacobian(q_reg,L);
J_rst = [J(4:6,:); J(1:3,:)];
disp(norm(J_rst - J_c)); 

%% Home pose verification
load("Arm.mat","robot","dh");
q_stow = [0; 80; 0; -165; 0; 0; 0] * pi/180;
L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
q = q_stow + dh(:,4);
T = getTransform(robot,q,"end_effector","base");
Tf = Forward_Kinematics(q_stow,L);
diff = norm(T-Tf);
tolerance = 1e-6;
if diff < tolerance
    fprintf('[SUCCESS] Forward kinematics verified\n');
    fprintf('Error norm: %.3e\n', diff);
else
    fprintf('[WARNING] Forward kinematics mismatch ✖\n');
    fprintf('Error norm: %.3e (exceeds tolerance %.1e)\n', diff, tolerance);
end

figure;
show(robot, q);
title('Minimal Stowable Pose');
view(0, 0);
