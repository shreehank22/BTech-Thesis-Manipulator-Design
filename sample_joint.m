function q = sample_joint(N, q_min, q_max)
    q = q_min + (q_max - q_min) .* rand(N, 7);
end