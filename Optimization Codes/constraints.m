function [c, ceq] = constraints(x)
    l1=x(1); l2=x(2); l3=x(3); l4=x(4);
    l5=x(5); l6=x(6); l7=x(7);
    L = sum(x);
    c = zeros(4,1);
    c(1) = l4 - l3;       
    c(2) = l5 - l4;     
    c(3) = l6 - l5;     
    c(4) = l7 - l6;       
    c(5) = 2.5 - L;       
    c(6) = L - 3.0; 
    c(7) = 0.2 - l1; 
    ceq = l1 - l2;
end