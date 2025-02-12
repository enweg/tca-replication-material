function [shock_size, idx] = get_shock_size(M_, shock_name)
    shocks = dynare_cellarray_to_vec(M_.exo_names);
    idx = find(shocks == shock_name);
    shock_size = sqrt(M_.Sigma_e(idx, idx));
end
