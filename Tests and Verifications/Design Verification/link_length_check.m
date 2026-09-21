L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
radii = repmat({[0.06, 0.057]}, 1, 7);
rho = 2700;
m_links = mass(L,radii,rho);
m_act = [2.2000; 2.7500; 1.5400; 1.5400; 0.7150; 0.7150; 0.7150];


% Mass of the arm
mass_arm = sum(m_links+m_act);
disp(mass_arm);

% Mass-weighted link index
weighted_sum = 0;
for i = 1:7
    weighted_sum = weighted_sum+i*m_links(i);
end
LI = sum((1:7)'.*m_links)/sum(m_links);
LI_norm = LI/7;
disp(LI_norm);

% Length norm
L_norm = zeros(7,1);
n=7;
for i=1:n
    L_norm(i,1) = L(i)/sum(L);
end
disp(L_norm);

% arm length fractions
wcr = (L(6) + L(7)) / (L(2) + L(3) + L(4));
arm_fraction = (L(1)+L(2)+L(3)+L(4)+L(5))/sum(L);
wrist_fraction = (L(6)+L(7))/sum(L);
disp(arm_fraction);disp(wrist_fraction);disp(wcr);

% Gravity Torque computation
q_home = [0; pi/2; 0; 0; -pi/2; 0; 0];
q_extend = zeros(7,1);
q_elbow = [0; 0; 0; pi/2; 0; 0; 0];

tau_home = gravity_Torque(q_home,m_links,L,m_act);
tau_extend = gravity_Torque(q_extend,m_links,L,m_act);
tau_elbow = gravity_Torque(q_elbow,m_links,L,m_act);

fprintf('\n%-10s  %-10s  %-10s  %-10s\n', 'Joint', 'Home', 'Extended', 'Elbow90');
fprintf('%s\n', repmat('-',1,44));
for i = 1:7
    fprintf('J%-9d  %-10.4f  %-10.4f  %-10.4f\n', i, tau_home(i), tau_extend(i), tau_elbow(i));
end

% COM
[~,com_home,~,~] = linkCOMPositions(q_home,L,m_links);
[~,com_extend,~,~] = linkCOMPositions(q_extend,L,m_links);
[~,com_elbow,~,~] = linkCOMPositions(q_elbow,L,m_links);

configs = {'Home', 'Extended', 'Elbow90'};
coms = {com_home, com_extend, com_elbow};

fprintf('\n%-10s  %-10s  %-10s  %-10s  %-15s  %-12s\n', ...
    'Config', 'COM_x(m)', 'COM_y(m)', 'COM_z(m)', 'r_radial(m)', 'r/L_total');
fprintf('%s\n', repmat('-',1,72));

for i = 1:3
    c = coms{i};
    r = norm(c(1:2));
    fprintf('%-10s  %-10.4f  %-10.4f  %-10.4f  %-15.4f  %-12.4f\n', ...
        configs{i}, c(1), c(2), c(3), r, r/sum(L));
end