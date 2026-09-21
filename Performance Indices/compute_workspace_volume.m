function [V, Pos] = compute_workspace_volume(link_lengths, joint_angles)
    N = size(joint_angles, 1);
    positions = zeros(N, 3);
    l1=link_lengths(1); l2=link_lengths(2); l3=link_lengths(3); l4=link_lengths(4);
    l5=link_lengths(5); l6=link_lengths(6); l7=link_lengths(7);

    for i = 1:N
        q = joint_angles(i,:);
        T1 = Transformation_Matrix(q(1),       l1,    +pi/2, 0);
        T2 = Transformation_Matrix(q(2)+pi/2,  0,     +pi/2, 0);
        T3 = Transformation_Matrix(q(3),       l2+l3, -pi/2, 0);
        T4 = Transformation_Matrix(q(4),       0,     +pi/2, 0);
        T5 = Transformation_Matrix(q(5)+pi/2,  l4+l5, +pi/2, 0);
        T6 = Transformation_Matrix(q(6)+pi/2,  0,     +pi/2, 0);
        T7 = Transformation_Matrix(q(7),       0,      0,    l6+l7);
        T = T1*T2*T3*T4*T5*T6*T7;
        positions(i,:) = T(1:3,4)';
    end

    try
        [~, V] = convhulln(positions);
    catch
        V = 0;
    end

    Pos = positions;
end