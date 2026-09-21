% Robot Parameters
L = [0.2,0.2,1.0,0.7424,0.4635,0.113,0.106];
l1 = L(1);l2 = L(2);l3 = L(3);l4 = L(4);l5 = L(5);l6 = L(6);l7 = L(7);
rho   = 2700;
r_out = 0.06;
r_in  = 0.057;
radii   = repmat({[r_out, r_in]}, 1, 7);
m_links = mass(L, radii, rho);
m_act   = [2.20; 2.75; 1.54; 1.54; 0.715; 0.715; 0.715];

% Initializing robot model
robot = rigidBodyTree('DataFormat','column','MaxNumBodies',10);
robot.Gravity = [0; 0; -9.81];
dh = [0,+pi/2,l1,0;0,+pi/2,0,pi/2;0,-pi/2,l2+l3,0;0,+pi/2,0,0;0,+pi/2,l4+l5,pi/2;0,+pi/2,0,pi/2;l6+l7,0,0,0;];
joint_names = {'J1_Yaw','J2_Pitch','J3_Roll','J4_Pitch','J5_Yaw','J6_Pitch','J7_Roll'};
body_names  = {'link1','link2','link3','link4','link5','link6','link7'};
prev_body = 'base';
q_limits = [170, 80, 80, 170, 80, 75, 170];

% Assembling the body
for i = 1:7
    body = rigidBody(body_names{i});
    joint = rigidBodyJoint(joint_names{i}, 'revolute');
    setFixedTransform(joint, [dh(i,1), dh(i,2), dh(i,3), dh(i,4)], 'dh');
    joint.PositionLimits = [-q_limits(i), q_limits(i)] * pi/180;
    link_length_i = norm(dh(i,1)) + abs(dh(i,3));
    if link_length_i < 1e-6
        body.Mass = 1e-6;
        body.CenterOfMass = [0; 0; 0];
        Ixx = 0;
        Iyy = 0;
        Izz = 0;
    else
        body.Mass = m_links(i);
        com_z = dh(i,3)/2;
        body.CenterOfMass = [0; 0; com_z];
        Ixx = (1/12)*m_links(i)*(3*(r_out^2+r_in^2)/2 + L(i)^2);
        Iyy = Ixx;
        Izz = (1/2)*m_links(i)*(r_out^2+r_in^2)/2;
    end
    body.Inertia = [Ixx, Iyy, Izz, 0, 0, 0];
    body.Joint = joint;
    addBody(robot, body, prev_body);
    prev_body = body_names{i};
end

% End effector
ee = rigidBody('end_effector');
ee_jnt = rigidBodyJoint('ee_fixed', 'fixed');
setFixedTransform(ee_jnt, eye(4));
ee.Joint = ee_jnt;
addBody(robot, ee, 'link7');

% Saving
save('Arm.mat', 'robot', 'dh', 'L', 'm_links', 'm_act', 'joint_names', 'body_names');
fprintf('Model saved: Arm.mat\n');