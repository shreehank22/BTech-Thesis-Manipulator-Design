function mu = compute_MI(q, link_lengths)
J = Geometric_Jacobian(q,link_lengths);
mu = sqrt(det(J*J'));
end
