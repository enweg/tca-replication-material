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
addpath("/Applications/Dynare/6.3-arm64/matlab/");
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

irfsDSGEea = [oo_.irfs.dy_ea' oo_.irfs.dc_ea' oo_.irfs.dinve_ea' oo_.irfs.pinfobs_ea' oo_.irfs.dw_ea' oo_.irfs.robs_ea' oo_.irfs.labobs_ea'];
irfsDSGEea = reshape(irfsDSGEea', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "ea");
test = max(vec(abs(irfsDSGEea - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for ea ~= VARMA IRFs")
end

irfsDSGEeb = [oo_.irfs.dy_eb' oo_.irfs.dc_eb' oo_.irfs.dinve_eb' oo_.irfs.pinfobs_eb' oo_.irfs.dw_eb' oo_.irfs.robs_eb' oo_.irfs.labobs_eb'];
irfsDSGEeb = reshape(irfsDSGEeb', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "eb");
test = max(vec(abs(irfsDSGEeb - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for eb ~= VARMA IRFs")
end

irfsDSGEeg = [oo_.irfs.dy_eg' oo_.irfs.dc_eg' oo_.irfs.dinve_eg' oo_.irfs.pinfobs_eg' oo_.irfs.dw_eg' oo_.irfs.robs_eg' oo_.irfs.labobs_eg'];
irfsDSGEeg = reshape(irfsDSGEeg', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "eg");
test = max(vec(abs(irfsDSGEeg - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for eg ~= VARMA IRFs")
end

irfsDSGEeqs = [oo_.irfs.dy_eqs' oo_.irfs.dc_eqs' oo_.irfs.dinve_eqs' oo_.irfs.pinfobs_eqs' oo_.irfs.dw_eqs' oo_.irfs.robs_eqs' oo_.irfs.labobs_eqs'];
irfsDSGEeqs = reshape(irfsDSGEeqs', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "eqs");
test = max(vec(abs(irfsDSGEeqs - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for eqs ~= VARMA IRFs")
end

irfsDSGEem = [oo_.irfs.dy_em' oo_.irfs.dc_em' oo_.irfs.dinve_em' oo_.irfs.pinfobs_em' oo_.irfs.dw_em' oo_.irfs.robs_em' oo_.irfs.labobs_em'];
irfsDSGEem = reshape(irfsDSGEem', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "em");
test = max(vec(abs(irfsDSGEem - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for em ~= VARMA IRFs")
end

irfsDSGEepinf = [oo_.irfs.dy_epinf' oo_.irfs.dc_epinf' oo_.irfs.dinve_epinf' oo_.irfs.pinfobs_epinf' oo_.irfs.dw_epinf' oo_.irfs.robs_epinf' oo_.irfs.labobs_epinf'];
irfsDSGEepinf = reshape(irfsDSGEepinf', n, 1, []);
[shock_size, idx] = get_shock_size(M_, "epinf");
test = max(vec(abs(irfsDSGEepinf - irfsVarma(:, idx, :) * shock_size))) < tol;
if ~test
    error("DSGE IRFs for epinf ~= VARMA IRFs")
end

irfsDSGEew = [oo_.irfs.dy_ew' oo_.irfs.dc_ew' oo_.irfs.dinve_ew' oo_.irfs.pinfobs_ew' oo_.irfs.dw_ew' oo_.irfs.robs_ew' oo_.irfs.labobs_ew'];
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
