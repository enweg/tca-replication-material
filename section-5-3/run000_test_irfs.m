%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING WHETHER VARMA IRFS SAME AS DYNARE IRFS
% This file does not need to be run in the replication process, unless 
% you want to make sure that the IRFs that dynare provides are the same 
% as the IRFs of the VARMA representation, i.e. unless you want to test
% whether the found VARMA representation is correct. 
%
% If any of the tests fail, an error message will be returned. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
addpath("/Applications/Dynare/5.5-arm64/matlab");
addpath("functions/")
cd SW2007;
dynare SW2007;
cd ..;
clc;

tol = 1e-9;  % accepted tolerance throughout all tests

[As, Psis, p, q] = dynare_to_varma(M_, oo_, options_);
irfsVarma = irf_varma(As, Psis, 19);

[A, B, C, D] = get_ABCD(M_, oo_, options_);
irfsSS = irf_state_space(A, B, C, D, 19);
test = max(vec(abs(irfsSS - irfsVarma))) < tol;
if ~test
    error("irfsSS ~= irfsVarma")
end

n = size(C, 1);

irfsDSGEea = [dy_ea dc_ea dinve_ea pinfobs_ea dw_ea robs_ea labobs_ea];
irfsDSGEea = reshape(irfsDSGEea', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "ea");
test = max(vec(abs(irfsDSGEea - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for ea ~= VARMA IRFs")
end

irfsDSGEeb = [dy_eb dc_eb dinve_eb pinfobs_eb dw_eb robs_eb labobs_eb];
irfsDSGEeb = reshape(irfsDSGEeb', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "eb");
test = max(vec(abs(irfsDSGEeb - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for eb ~= VARMA IRFs")
end

irfsDSGEeg = [dy_eg dc_eg dinve_eg pinfobs_eg dw_eg robs_eg labobs_eg];
irfsDSGEeg = reshape(irfsDSGEeg', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "eg");
test = max(vec(abs(irfsDSGEeg - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for eg ~= VARMA IRFs")
end

irfsDSGEeqs = [dy_eqs dc_eqs dinve_eqs pinfobs_eqs dw_eqs robs_eqs labobs_eqs];
irfsDSGEeqs = reshape(irfsDSGEeqs', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "eqs");
test = max(vec(abs(irfsDSGEeqs - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for eqs ~= VARMA IRFs")
end

irfsDSGEem = [dy_em dc_em dinve_em pinfobs_em dw_em robs_em labobs_em];
irfsDSGEem = reshape(irfsDSGEem', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "em");
test = max(vec(abs(irfsDSGEem - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for em ~= VARMA IRFs")
end

irfsDSGEepinf = [dy_epinf dc_epinf dinve_epinf pinfobs_epinf dw_epinf robs_epinf labobs_epinf];
irfsDSGEepinf = reshape(irfsDSGEepinf', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "epinf");
test = max(vec(abs(irfsDSGEepinf - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for epinf ~= VARMA IRFs")
end

irfsDSGEew = [dy_ew dc_ew dinve_ew pinfobs_ew dw_ew robs_ew labobs_ew];
irfsDSGEew = reshape(irfsDSGEew', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "ew");
test = max(vec(abs(irfsDSGEew - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for ew ~= VARMA IRFs")
end

[B, Oomega] = varma_to_static(As, Psis, 19, eye(n));
irfsStatic = irf_static_model(M_, B, Oomega, n);
test = max(vec(abs(irfsStatic - irfsVarma))) < tol;
if ~test
    error("Static IRFs ~= VARMA IRFs")
end
