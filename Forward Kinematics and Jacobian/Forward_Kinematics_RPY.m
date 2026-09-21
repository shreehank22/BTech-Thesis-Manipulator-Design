function [T_matrix, ee_pose, T_all] = Forward_Kinematics_RPY(q, x)
q = double(q(:));
x = double(x(:));

l1=x(1); l2=x(2); l3=x(3); l4=x(4); l5=x(5); l6=x(6); l7=x(7);

T{1} = Transformation_Matrix(q(1),l1,+pi/2, 0);   
T{2} = Transformation_Matrix(q(2)+pi/2,0,+pi/2,0);   
T{3} = Transformation_Matrix(q(3),l2+l3,-pi/2,0);   
T{4} = Transformation_Matrix(q(4),0,+pi/2, 0);   
T{5} = Transformation_Matrix(q(5)+pi/2,l4+l5,+pi/2,0);   
T{6} = Transformation_Matrix(q(6),0,-pi/2, 0);   
T{7} = Transformation_Matrix(q(7)+pi/2,0,0,l6+l7);

T_all = cell(7,1);
T_all{1} = T{1};
for i = 2:7
    T_all{i} = T_all{i-1} * T{i};
end

T_matrix = T_all{7};
ee_pose = cell(7,1);
for i = 1:7
    ee_pose{i} = T_all{i}(1:3,4);
end
end