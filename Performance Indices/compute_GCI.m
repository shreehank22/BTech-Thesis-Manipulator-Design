function gci_value = compute_GCI(joint_angles, link_lengths)
N = size(joint_angles,1);
LCI = zeros(N,1);

for i = 1:N
    LCI(i) = compute_LCI(joint_angles(i,:), link_lengths);
end

gci_value = mean(LCI);
end
