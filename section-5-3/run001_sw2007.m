%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPLICATION INSTRUCTIONS: 
% 1. Adjust the path to dynare in line 12 to the dynare installation location 
%    on your computer.
% 2. Run the entire file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% ADJUST THE FOLLOWING LINE
addpath("/Applications/Dynare/5.5-arm64/matlab");

addpath("functions/")
cd SW2007;
dynare SW2007;
cd ..;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTHING NEEDS TO BE ADJUSTED IN THE LINES BELOW.
%
% Preparation for TCA. 
% Question: How important are wages in the transmission of monetary policy? 
% Do the transmission effects show patterns in line with the idea of 
% second round effects? 
%
% Transmission Matrix Decision: 
% - Order interest rates first because it is the naturally associated variable 
%   with a monetary policy shock. 
% - Order wages second. That way, we allow for the possibility that the other 
%   variables only partially adjust if effects cannot go through wages. 
%   If wages are mostly involved in second round effects, then ordering wages
%   before or after the other variables should not make a difference. 
%   ==> Will be investgated using `run003_sw2007.m` file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[As, Psis, p, q] = dynare_to_varma(M_, oo_, options_);
[shock_size, ix_em] = get_shock_size(M_, "em");

k = size(As{1}, 1);
vars = dynare_cellarray_to_vec(options_.varobs);
ix_r = find(vars == "robs");
ix_w = find(vars == "dw");
ix_p = find(vars == "pinfobs");
ix_c = find(vars == "dc");
ix_inv = find(vars == "dinve");
ix_y = find(vars == "dy");
ix_l = find(vars == "labobs");

T = eye(k);
order = [ix_r, ix_w, ix_l, ix_c, ix_inv, ix_y, ix_p];
T = T(order, :);

vars_original = vars;

horizon = 20;  
[B, Oomega] = varma_to_static(As, Psis, horizon, T);
vars = vars_original(order);
ix_r = find(vars == "robs");
ix_w = find(vars == "dw");
ix_p = find(vars == "pinfobs");
ix_c = find(vars == "dc");
ix_inv = find(vars == "dinve");
ix_y = find(vars == "dy");
ix_l = find(vars == "labobs");
irfs = irf_static_model(M_, B, Oomega, k) * shock_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTHING NEEDS TO BE ADJUSTED IN THE FOLLOWING LINES. 
%
% Definition of transmission channel and computation of transmission effects. 
% The wage channel will be defined as the effect going through wages in some 
% period. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

effect_not_w = through_not_x(M_, B, Oomega, [ix_w], k) * shock_size;
effect_w = irfs - effect_not_w; 

% saving the decomposition
df = table(vec(irfs(ix_p, ix_em, :)), vec(effect_not_w(ix_p, ix_em, :)), vec(effect_w(ix_p, ix_em, :)));
df.Properties.VariableNames = {'total', 'not_w', 'w'};
writetable(df, sprintf("output/effects-horizons-%d.csv", horizon));

% saving the total IRFs
irfs = squeeze(irfs(:, ix_em, :));
irfs = irfs';
irfs_table = array2table(irfs, 'VariableNames', vars);
writetable(irfs_table, sprintf("output/irfs-horizons-%d.csv", horizon));

