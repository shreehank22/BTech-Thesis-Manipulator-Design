function GMI = compute_GMI(link_lengths, Qsamples)
    N = size(Qsamples, 1);   
    GMI = 0;
    for i = 1:N
        q = Qsamples(i,:);
        mu = compute_MI(q, link_lengths);
        GMI = GMI + mu;
    end
    GMI = GMI / N;
end