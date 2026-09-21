function FK_IK_tests(L,q_min,q_max)
   fprintf('\n===== FK/IK Validation Test =====\n\n');
   q_test = zeros(7,1);
   [T_fk,~,~] = Forward_Kinematics(q_test,L);
   fprintf('FK at q=zeros:\n');
   disp(T_fk);

   q_target = [pi/6, pi/4, -pi/3, pi/4, -pi/6, pi/3, -pi/4]';
   [T_home,~,~]=Forward_Kinematics(q_target,L);
   fprintf('FK at q_target:\n');
   disp(T_home);

   % Verifying IK
   q_ik = inverse_kinematics(q_test,T_home,L,q_min,q_max);
   [T_check] = Forward_Kinematics(q_ik,L);
   fprintf('IK solution FK result:\n');
   disp(T_check);
   fprintf('IK roundtrip pos error: %.3e m\n',   norm(T_check(1:3,4) - T_home(1:3,4)));
   fprintf('IK roundtrip rot error: %.3e rad\n', norm(logSO3(T_check(1:3,1:3)*T_home(1:3,1:3)')));
   fprintf('\nTarget EE position : [%.6f, %.6f, %.6f] m\n', T_home(1,4), T_home(2,4), T_home(3,4));
   fprintf('Target rotation matrix:\n');
   fprintf('  [%.6f  %.6f  %.6f]\n', T_home(1,1), T_home(1,2), T_home(1,3));
   fprintf('  [%.6f  %.6f  %.6f]\n', T_home(2,1), T_home(2,2), T_home(2,3));
   fprintf('  [%.6f  %.6f  %.6f]\n', T_home(3,1), T_home(3,2), T_home(3,3));
   n_iter = 10;
   fprintf('\n--- FK/IK Roundtrip Loop (%d iterations) ---\n', n_iter);fprintf('\n--- FK/IK Roundtrip Loop (10 iterations) ---\n');
   fprintf('%-6s %-14s %-14s\n', 'Iter', 'pos_error(m)', 'rot_error(rad)');
   errors = zeros(n_iter, 1);
   q_curr = q_test;    
for i = 1:n_iter
    q_curr = inverse_kinematics(q_curr, T_home, L, q_min, q_max);
    [T_curr, ~, ~] = Forward_Kinematics(q_curr, L);
    fprintf('\n--- Iteration %d ---\n', i);
    fprintf('EE position : [%.6f, %.6f, %.6f] m\n', T_curr(1,4), T_curr(2,4), T_curr(3,4));
    fprintf('Rotation matrix:\n');
    fprintf('  [%.6f  %.6f  %.6f]\n', T_curr(1,1), T_curr(1,2), T_curr(1,3));
    fprintf('  [%.6f  %.6f  %.6f]\n', T_curr(2,1), T_curr(2,2), T_curr(2,3));
    fprintf('  [%.6f  %.6f  %.6f]\n', T_curr(3,1), T_curr(3,2), T_curr(3,3));
    pos_err   = norm(T_curr(1:3,4) - T_home(1:3,4));
    rot_err   = norm(logSO3(T_curr(1:3,1:3) * T_home(1:3,1:3)'));
    errors(i) = pos_err;
    fprintf('%-6d %-14.3e %-14.3e\n', i, pos_err, rot_err);
end
fprintf('\n--- Roundtrip Statistics (%d iterations) ---\n', n_iter);
fprintf('Mean error : %.3e m\n', mean(errors));
fprintf('Std dev    : %.3e m\n', std(errors));
fprintf('Min error  : %.3e m\n', min(errors));
fprintf('Max error  : %.3e m\n', max(errors));