function [T_matrix, ee_pose, T_all] = Forward_Kinematics(q, x)
q = double(q(:));
x = double(x(:));
q = q(:);

q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
q5 = q(5); q6 = q(6); q7 = q(7);

l1 = x(1); l2 = x(2); l3 = x(3);
l4 = x(4); l5 = x(5); l6 = x(6); l7 = x(7);

% Now Y-P-R-P
T1 = Transformation_Matrix(q1,l1,+pi/2,0);  % Y
T2 = Transformation_Matrix(q2+pi/2,0,+pi/2,0);  % P
T3 = Transformation_Matrix(q3,l2+l3,-pi/2,0);  % R
T4 = Transformation_Matrix(q4,0,+pi/2,0);  % P 
T5 = Transformation_Matrix(q5+pi/2,l4 + l5,+pi/2,0);  % R 
T6 = Transformation_Matrix(q6+pi/2,0,+pi/2,0);  % Y 
T7 = Transformation_Matrix(q7,0,0,l6 + l7);  % P 


T = {T1, T2, T3, T4, T5, T6, T7};
T_all = cell(7,1); %0T7

T_all{1} = T{1};
for i = 2:7
    T_all{i} = T_all{i-1} * T{i};
end

T_matrix = T_all{7};
ee_pose = cell(7,1);
for i=1:7
    ee_pose{i} = T_all{i}(1:3,4);
end

end
