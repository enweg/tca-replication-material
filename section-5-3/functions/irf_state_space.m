function irfs = irf_state_space(A, B, C, D, horizon)
    n = size(C, 1);
    m = size(A, 1);
    ATilde = [A zeros(m, n); C zeros(n, n)];
    DTilde = [B; D];
    irfs = zeros(n, n, horizon+1);
    for h=0:horizon 
        tmp = ATilde^h*DTilde;
        irfs(:, :, h+1) = tmp((m+1):end, :);
    end
end
