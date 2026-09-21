function run_trajectory()
%% ================================================================
%  Circular Trajectory Tracking — Proposed 7-DOF Arm
%  IK: Damped Least Squares (DLS)
%  Visualisation: MATLAB Animation + Dashboard
%% ================================================================

clear; clc; close all;
load('Arm.mat');

%% --- Robot Parameters ---
L = [0.2, 0.2, 1.0, 0.7424, 0.4635, 0.113, 0.106];

%% --- Circle Parameters ---
centre   = [1.0; 0.0; 0.5];
radius   = 0.3;
n_points = 100;
theta_c  = linspace(0, 2*pi, n_points);

%% --- Generate Circle in XZ Plane ---
p_traj = zeros(3, n_points);
for i = 1:n_points
    p_traj(:,i) = centre + ...
        radius * [cos(theta_c(i)); 0; sin(theta_c(i))];
end

%% --- Fixed End-Effector Orientation ---
R_des = eye(3);

%% --- Home Configuration ---
q_home = [0; pi/2; 0; 0; -pi/2; 0; 0];

%% --- Solve IK for Each Waypoint ---
fprintf('Solving IK for %d waypoints...\n', n_points);
Q_traj    = zeros(7, n_points);
pos_errs  = zeros(n_points, 1);
rot_errs  = zeros(n_points, 1);
converged = false(n_points, 1);
n_iters   = zeros(n_points, 1);

q_curr = q_home;
for i = 1:n_points
    [q_sol, conv, n_it, pe, re] = DLS(q_curr, L, p_traj(:,i), R_des);
    Q_traj(:,i)  = q_sol;
    pos_errs(i)  = pe;
    rot_errs(i)  = re;
    converged(i) = conv;
    n_iters(i)   = n_it;
    q_curr       = q_sol;
    if mod(i,10) == 0
        fprintf('  Waypoint %3d/%d — pos_err=%.2e  conv=%d\n', ...
            i, n_points, pe, conv);
    end
end

fprintf('Converged: %d/%d waypoints\n', sum(converged), n_points);

%% --- Compute FK for All Solutions ---
ee_actual = zeros(3, n_points);
for i = 1:n_points
    [T,~,~] = Forward_Kinematics(Q_traj(:,i), L);
    ee_actual(:,i) = T(1:3,4);
end

%% ================================================================
%  FIGURE 1 — Static Overview: Trajectory + Key Poses
%% ================================================================
key_idx = round(linspace(1, n_points, 8));

figure('Color','k','Position',[50 50 900 700]);
hold on;

plot3(p_traj(1,:), p_traj(2,:), p_traj(3,:), ...
    'w--', 'LineWidth', 2, 'DisplayName', 'Desired');
plot3(ee_actual(1,:), ee_actual(2,:), ee_actual(3,:), ...
    'c-', 'LineWidth', 1.5, 'DisplayName', 'Actual EE');

colors = jet(length(key_idx));
for k = 1:length(key_idx)
    plot_arm_links(Q_traj(:,key_idx(k)), L, colors(k,:));
end

scatter3(p_traj(1,1),   p_traj(2,1),   p_traj(3,1),   ...
    120, 'g', 'filled', 'DisplayName', 'Start');
scatter3(p_traj(1,end), p_traj(2,end), p_traj(3,end), ...
    120, 'r', 'filled', 'DisplayName', 'End');

xlabel('X (m)','Color','w','FontSize',11);
ylabel('Y (m)','Color','w','FontSize',11);
zlabel('Z (m)','Color','w','FontSize',11);
title('Circular Trajectory — Key Configurations', ...
    'Color','w','FontSize',13);
legend('Desired','Actual EE','Start','End', ...
    'TextColor','w','Color','k','Location','best');
ax = gca;
ax.Color = 'k'; ax.XColor = 'w';
ax.YColor = 'w'; ax.ZColor = 'w';
ax.GridColor = 'w'; grid on; axis equal; view(35,20);
saveas(gcf, 'Trajectory_Static.png');

%% ================================================================
%  FIGURE 2 — Dashboard: Errors + Joint Profiles
%% ================================================================
figure('Color','w','Position',[50 50 1200 500]);

subplot(1,3,1);
plot(1:n_points, pos_errs*1000, 'b-', 'LineWidth', 1.5);
xlabel('Waypoint','FontSize',11);
ylabel('Position Error (mm)','FontSize',11);
title('End-Effector Position Error','FontSize',12,'FontWeight','bold');
yline(1, 'r--', '1 mm', 'LineWidth', 1.2, 'FontSize', 9);
grid on; box off;

subplot(1,3,2);
plot(1:n_points, rot_errs, 'r-', 'LineWidth', 1.5);
xlabel('Waypoint','FontSize',11);
ylabel('Rotation Error (rad)','FontSize',11);
title('End-Effector Rotation Error','FontSize',12,'FontWeight','bold');
grid on; box off;

subplot(1,3,3);
colors_j = lines(7);
hold on;
labels = {'J1','J2','J3','J4','J5','J6','J7'};
for j = 1:7
    plot(1:n_points, rad2deg(Q_traj(j,:)), ...
        'Color', colors_j(j,:), 'LineWidth', 1.2, ...
        'DisplayName', labels{j});
end
xlabel('Waypoint','FontSize',11);
ylabel('Joint Angle (deg)','FontSize',11);
title('Joint Angle Profiles','FontSize',12,'FontWeight','bold');
legend('Location','best','FontSize',8);
grid on; box off;

sgtitle('DLS Inverse Kinematics — Circular Trajectory Tracking', ...
    'FontSize',13,'FontWeight','bold');
saveas(gcf, 'Trajectory_Dashboard.png');

%% ================================================================
%  FIGURE 3 — MATLAB Animation
%% ================================================================
figure('Color','k','Position',[100 100 900 700]);
hold on;

plot3(p_traj(1,:), p_traj(2,:), p_traj(3,:), ...
    'w--', 'LineWidth', 1.5);

h_links  = plot3(nan,nan,nan, 'c-',  'LineWidth', 3);
h_joints = plot3(nan,nan,nan, 'wo',  'MarkerSize', 8, ...
    'MarkerFaceColor','w');
h_ee     = plot3(nan,nan,nan, 'g.', 'MarkerSize', 18);
h_trail  = plot3(nan,nan,nan, 'y-', 'LineWidth', 1.5);

xlabel('X (m)','Color','w','FontSize',11);
ylabel('Y (m)','Color','w','FontSize',11);
zlabel('Z (m)','Color','w','FontSize',11);
title('Circular Trajectory Animation','Color','w','FontSize',13);
ax = gca;
ax.Color = 'k'; ax.XColor = 'w';
ax.YColor = 'w'; ax.ZColor = 'w';
ax.GridColor = 'w'; grid on;
view(35,20);

%% --- Compute full arm geometry bounds across all waypoints ---
fprintf('Computing arm geometry bounds...\n');
all_origins = zeros(3, 8*n_points);
for i = 1:n_points
    [~,~,T_all] = Forward_Kinematics(Q_traj(:,i), L);
    p_orig = zeros(3,8);
    p_orig(:,1) = [0;0;0];
    for k = 1:7
        p_orig(:,k+1) = T_all{k}(1:3,4);
    end
    all_origins(:, (i-1)*8+1 : i*8) = p_orig;
end

pad = 0.3;
x_all = [all_origins(1,:), p_traj(1,:)];
y_all = [all_origins(2,:), p_traj(2,:)];
z_all = [all_origins(3,:), p_traj(3,:)];

xlim([min(x_all)-pad, max(x_all)+pad]);
ylim([min(y_all)-pad, max(y_all)+pad]);
zlim([min(z_all)-pad, max(z_all)+pad]);

h_text = text(0.02, 0.95, '', 'Units','normalized', ...
    'Color','w', 'FontSize',10, 'FontWeight','bold');

trail_x = nan(n_points,1);
trail_y = nan(n_points,1);
trail_z = nan(n_points,1);

fprintf('\nRunning animation...\n');

for i = 1:n_points
    p_orig = all_origins(:, (i-1)*8+1 : i*8);

    set(h_links,  'XData', p_orig(1,:), ...
                  'YData', p_orig(2,:), ...
                  'ZData', p_orig(3,:));
    set(h_joints, 'XData', p_orig(1,:), ...
                  'YData', p_orig(2,:), ...
                  'ZData', p_orig(3,:));
    set(h_ee,     'XData', ee_actual(1,i), ...
                  'YData', ee_actual(2,i), ...
                  'ZData', ee_actual(3,i));

    trail_x(i) = ee_actual(1,i);
    trail_y(i) = ee_actual(2,i);
    trail_z(i) = ee_actual(3,i);
    set(h_trail,  'XData', trail_x(1:i), ...
                  'YData', trail_y(1:i), ...
                  'ZData', trail_z(1:i));

    set(h_text, 'String', sprintf(...
        'Waypoint : %d/%d\nPos error : %.3f mm\nConverged : %d', ...
        i, n_points, pos_errs(i)*1000, converged(i)));

    drawnow limitrate;
    pause(0.03);
end

fprintf('Animation complete.\n');
saveas(gcf, 'Trajectory_Animation_Final.png');

end  % <-- end of run_trajectory function

%% ================================================================
%  HELPER FUNCTION
%% ================================================================
function plot_arm_links(q, L, col)
    [~,~,T_all] = Forward_Kinematics(q, L);
    p = zeros(3,8);
    p(:,1) = [0;0;0];
    for k = 1:7
        p(:,k+1) = T_all{k}(1:3,4);
    end
    plot3(p(1,:), p(2,:), p(3,:), '-', ...
        'Color', col, 'LineWidth', 1.5);
    plot3(p(1,:), p(2,:), p(3,:), 'o', ...
        'Color', col, 'MarkerSize', 5, ...
        'MarkerFaceColor', col);
end