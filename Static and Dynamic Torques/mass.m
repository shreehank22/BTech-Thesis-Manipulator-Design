function m = mass(L, radii, rho)
    N = length(L);
    m = zeros(N, 1);
    for i = 1:N
        r_out = radii{i}(1);
        r_in  = radii{i}(2);
        m(i)  = rho * pi * (r_out^2 - r_in^2) * L(i);
    end
end