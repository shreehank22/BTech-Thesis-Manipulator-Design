function Inertia_Tensor = inertia_tensor(Body, Mass, Material)
    r = Material(2);
    Inertia_Tensor = zeros(3,3,length(Body));
    for i = 1:length(Body)
        Ixx=(1/12)*Mass(i)*(3*r^2+Body(i)^2);
        Iyy=(1/12)*Mass(i)*(3*r^2+Body(i)^2);
        Izz=(1/2)*Mass(i)*r^2;
        Inertia_Tensor(:,:,i)=diag([Ixx Iyy Izz]);
    end
end