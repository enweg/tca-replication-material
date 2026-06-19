function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = model_extensive_gap_no_investment.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 1
    T = [T; NaN(1 - size(T, 1), 1)];
end
end
