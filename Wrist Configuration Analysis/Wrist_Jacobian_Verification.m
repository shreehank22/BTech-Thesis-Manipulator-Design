load("Arm.mat","robot");
L = [0.1040, 0.1040, 0.9993, 0.9072, 0.6779, 0.1053, 0.1013];
q_init = zeros(7,1);
q = q_init + dh(:,4);
J = geometricJacobian(robot,q,"end_effector");


J_rst = [J(4:6,:); J(1:3,:)];
J_c = Geometric_Jacobian(q_init,L);

J_w_rst = J_rst(4:6, 5:7);
J_w_c = J_c(4:6, 5:7);
disp(norm(J_w_c-J_w_rst));