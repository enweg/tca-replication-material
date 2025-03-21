%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPLICATION INSTRUCTIONS: 
% 1. Adjust the path to dynare in line 12 to the dynare installation location 
%    on your computer.
% 2. Run the entire file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% ADJUST THE FOLLOWING LINE
addpath("/Applications/Dynare/6.3-arm64/matlab/");

addpath("functions/")
cd SW2007;
% The following line is the only part of the code that requires Dynare. 
% It finds the first-order approximation to the DSGE. 
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

% We work with VARMA representations, so we first transform the Dynare DSGE
% model to VARMA form.
[As, Psis, p, q] = dynare_to_varma(M_, oo_, options_);
% We are interested in the monetary policy shock, here denoted by em. 
[shock_size, ix_em] = get_shock_size(M_, "em");

% To create the transmission matrix, we first find the index that each
% variables has in the original ordering. 
k = size(As{1}, 1);
vars = dynare_cellarray_to_vec(options_.varobs);
ix_r = find(vars == "robs");     % interest rates
ix_w = find(vars == "dw");       % wages
ix_p = find(vars == "pinfobs");  % inflation
ix_c = find(vars == "dc");       % consumption
ix_inv = find(vars == "dinve");  % investments
ix_y = find(vars == "dy");       % output
ix_l = find(vars == "labobs");   % labour hours

% We then use these indices to create the transmission matrix. 
% The transmission matrix corresponds to the ordering 
% interest rates, wages, labour hours, consumptions, investments, output, inflation.
% As the paper states, [labour hours, consumption, investments, output] can be 
% re-ordered within the group. 
order = [ix_r, ix_w, ix_l, ix_c, ix_inv, ix_y, ix_p];
T = eye(k);
T = T(order, :);

vars_original = vars;

% Defining the maximum horizon for all IRFs and for TCA.
horizon = 20;  
% We transform the VARMA to the systems / static representation which we 
% then use to compute the transmission effects. 
[B, Oomega] = varma_to_static(As, Psis, horizon, T);
% The variables are in a different ordering than the original ordering. We thus
% need to find the indices for all variables again. 
vars = vars_original(order);
ix_r = find(vars == "robs");
ix_w = find(vars == "dw");
ix_p = find(vars == "pinfobs");
ix_c = find(vars == "dc");
ix_inv = find(vars == "dinve");
ix_y = find(vars == "dy");
ix_l = find(vars == "labobs");
% Computing the IRFs -- the total effects. We adjust for the desired shock size
% defined in the Dynare model.
irfs = irf_static_model(M_, B, Oomega, k) * shock_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTHING NEEDS TO BE ADJUSTED IN THE FOLLOWING LINES. 
%
% Definition of transmission channel and computation of transmission effects. 
% The wage channel will be defined as the effect going through wages in some 
% period. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We define the demand channel as the effect not going through wages in any 
% period. The effect along the demand channel can be found using `through_not_x`
% which takes the index of  the variables (ix_w) and computes the effect not 
% going through this variable in any period. We again adjust for the shock size.
effect_not_w = through_not_x(M_, B, Oomega, [ix_w], k) * shock_size;
% The wage channel is defined as the complement to the demand channel.
% According to the theory in the paper, this implies that the effect through 
% the wage channel is the difference between the total effect and the effect
% through the demand channel. 
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

