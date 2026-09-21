function lci_value = compute_LCI(theta, link_lengths)

J = Geometric_Jacobian(theta, link_lengths);

sigma = svd(J);

sigma_max = sigma(1);
sigma_min = sigma(end);

if sigma_min <= eps
    lci_value = 0;
else
    lci_value = sigma_min / sigma_max;
end

end