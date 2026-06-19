function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = model_extensive_gap_no_investment.static_resid_tt(T, y, x, params);
end
residual = zeros(33, 1);
    residual(1) = (y(8)) - (y(32)+y(7));
    residual(2) = (y(10)) - (y(32)+y(9));
    residual(3) = (y(31)) - (y(31)*params(14)+x(1));
    residual(4) = (y(32)) - (y(32)*params(15)+x(2));
    residual(5) = (y(1)) - (y(1)-T(1)*(y(15)-y(17))-y(31)*T(1)*(params(14)-1));
    residual(6) = (y(2)) - (y(2)-T(1)*y(16)-y(31)*T(1)*(params(14)-1));
    residual(7) = (params(4)*y(18)) - (y(21)-y(1)*params(3));
    residual(8) = (params(4)*y(19)) - (y(22)-params(3)*y(2));
    residual(9) = (y(24)) - (params(1)*params(5)*y(27)+y(24)*params(1)*(1-params(2))-(y(15)-y(17)));
    residual(10) = (y(25)) - (params(1)*params(5)*y(28)+params(1)*(1-params(2))*y(25)-y(16));
residual(11) = y(24);
residual(12) = y(25);
    residual(13) = (y(7)) - (params(2)*y(4)+y(7)*(1-params(2)));
    residual(14) = (y(9)) - (params(2)*y(5)+y(9)*(1-params(2)));
    residual(15) = (y(27)) - (y(30)+y(12)-y(7));
    residual(16) = (y(28)) - (y(13)-y(9));
    residual(17) = (y(21)) - (y(30)+y(12)-y(18));
    residual(18) = (y(22)) - (y(13)-y(19));
    residual(19) = (y(30)) - (y(27)*params(9)-y(32)*params(9)+y(21)*(1-params(9)));
    residual(20) = (0) - (y(28)*params(9)-y(32)*params(9)+y(22)*(1-params(9)));
    residual(21) = (y(17)) - (y(17)*params(1)+y(30)*(1-params(10))*(1-params(1)*params(10))/params(10));
    residual(22) = (y(15)) - (y(15)*params(17)+(1-params(17))*(y(17)*params(12)+params(13)*y(14))+y(33));
    residual(23) = (y(33)) - (y(33)*params(16)+x(3));
    residual(24) = (y(12)) - (y(1)*params(6)+y(4)*params(7));
    residual(25) = (y(13)) - (y(2)*params(6)+y(5)*params(7));
    residual(26) = (y(3)) - (y(1)-y(2));
    residual(27) = (y(6)) - (y(4)-y(5));
    residual(28) = (y(11)) - (y(7)-y(9));
    residual(29) = (y(14)) - (y(12)-y(13));
    residual(30) = (y(20)) - (y(18)-y(19));
    residual(31) = (y(23)) - (y(21)-y(22));
    residual(32) = (y(26)) - (y(24)-y(25));
    residual(33) = (y(29)) - (y(27)-y(28));

end
