function e = logSO3(R)
    trace_val = (trace(R) - 1) / 2;
    theta  = acos(max(-1, min(1, trace_val)));   
    if abs(theta) < 1e-7
        skew = 0.5 * (R - R');
    elseif abs(theta - pi) < 1e-4
        A = (R + eye(3)) / 2;
        n = sqrt(max(0, diag(A)));
        if A(1,2) < 0, n(2) = -n(2); end
        if A(1,3) < 0, n(3) = -n(3); end
        e = pi * n;
        return
    else
        skew = (theta / (2*sin(theta))) * (R - R');
    end
    e = [skew(3,2); skew(1,3); skew(2,1)];
end