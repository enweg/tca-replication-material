%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes transmission effects for the setting in which the output-wage
% channel is defined in a narrow sense and the expectations channel in a
% broader sense. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear working environment and run the Smets & Wouters (2007) model 
% using Dynare.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
% the following lines needs to be adjusted to the system on which this 
% script is run.
addpath("/Applications/Dynare/5.5-arm64/matlab");
addpath("functions/")
cd Smets_Wouters_2007;
dynare Smets_Wouters_2007_45;
cd ..;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmission Channel Analsysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs_var = ["dy", "dc", "dinve", "piexp", "pinfobs", "robs", "dw"];
obs_var = get_obsvar_ids(M_, options_, oo_, obs_var);
[p, q] = determine_varma_order(M_, options_, oo_, 10, obs_var);
[A0, Phis, Psis, vars] = get_varma_coeffs(M_, options_, oo_, p, q, obs_var);
k = size(A0, 1);
ix_r = find(vars == "robs");
ix_w = find(vars == "dw");
ix_p = find(vars == "pinfobs");
ix_c = find(vars == "dc");
ix_inv = find(vars == "dinve");
ix_y = find(vars == "dy");
ix_e = find(vars == "piexp");
T = eye(k);
% order = [ix_r, ix_e, ix_w, ix_c, ix_inv, ix_y, setdiff(1:k, [ix_r, ix_e, ix_w, ix_c, ix_inv, ix_y, ix_p]), ix_p];
order = [ix_r, ix_w, ix_c, ix_inv, ix_y, ix_e, setdiff(1:k, [ix_r, ix_w, ix_c, ix_inv, ix_y, ix_e, ix_p]), ix_p];
T = T(order, :);


vars_original = vars;
horizons = [16, 32, 40, 52, 100];
for i=1:length(horizons)
    horizon = horizons(i);
    [B, Oomega, k, vars] = varma_to_static(vars_original, A0, Phis, Psis, horizon, T);
    ix_r = find(vars == "robs");
    ix_w = find(vars == "dw");
    ix_p = find(vars == "pinfobs");
    ix_c = find(vars == "dc");
    ix_inv = find(vars == "dinve");
    ix_y = find(vars == "dy");
    ix_e = find(vars == "piexp");

    shocks = dynare_cellarray_to_vec(M_.exo_names);
    ix_em = find(shocks == "em");

    shock_size = sqrt(M_.Sigma_e(ix_em, ix_em));
    irfs = irf_static_model(M_, B, Oomega, k) * shock_size;

    % decomposing the total effect into transmission effects
    effect_not_ywe = through_not_x(M_, B, Oomega, [ix_w, ix_c, ix_inv, ix_y, ix_e], k) * shock_size;
    effect_not_e = through_not_x(M_, B, Oomega, [ix_e], k) * shock_size;
    % effect_not_yw = through_not_x(M_, B, Oomega, [ix_w, ix_c, ix_inv, ix_y], k) * shock_size;
    effect_y_not_e = effect_not_e - effect_not_ywe;
    % effect_e_not_yw = effect_not_yw - effect_not_ywe;
    effect_e = irfs - effect_not_e;
    % effect_y_or_w = irfs - effect_not_yw;

    df = table(vec(irfs(ix_p, ix_em, :)), vec(effect_not_ywe(ix_p, ix_em, :)), vec(effect_y_not_e(ix_p, ix_em, :)), vec(effect_e(ix_p, ix_em, :)));
    df.Properties.VariableNames = {'total', 'interest_rate', 'output_wage', 'expectations'};
    writetable(df, sprintf("output/alternative-effects-horizons-%d.csv", horizon));
end

