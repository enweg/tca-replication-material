function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = model_extensive_gap_no_investment.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(33, 1);
    residual(1) = (y(16)) - (y(40)+y(3));
    residual(2) = (y(18)) - (y(40)+y(4));
    residual(3) = (y(39)) - (params(14)*y(6)+x(it_, 1));
    residual(4) = (y(40)) - (params(15)*y(7)+x(it_, 2));
    residual(5) = (y(9)) - (y(42)-T(1)*(y(23)-y(46))-y(39)*T(1)*(params(14)-1));
    residual(6) = (y(10)) - (y(43)-T(1)*y(24)-y(39)*T(1)*(params(14)-1));
    residual(7) = (params(4)*y(26)) - (y(29)-y(9)*params(3));
    residual(8) = (params(4)*y(27)) - (y(30)-params(3)*y(10));
    residual(9) = (y(32)) - (params(1)*params(5)*y(49)+params(1)*(1-params(2))*y(47)-(y(23)-y(46)));
    residual(10) = (y(33)) - (params(1)*params(5)*y(50)+params(1)*(1-params(2))*y(48)-y(24));
    residual(11) = (y(32)) - (params(8)*(y(12)-y(1)-params(1)*(y(44)-y(12))));
    residual(12) = (y(33)) - (params(8)*(y(13)-y(2)-params(1)*(y(45)-y(13))));
    residual(13) = (y(15)) - (params(2)*y(12)+(1-params(2))*y(3));
    residual(14) = (y(17)) - (params(2)*y(13)+(1-params(2))*y(4));
    residual(15) = (y(35)) - (y(38)+y(20)-y(3));
    residual(16) = (y(36)) - (y(21)-y(4));
    residual(17) = (y(29)) - (y(38)+y(20)-y(26));
    residual(18) = (y(30)) - (y(21)-y(27));
    residual(19) = (y(38)) - (y(35)*params(9)-y(40)*params(9)+y(29)*(1-params(9)));
    residual(20) = (0) - (y(36)*params(9)-y(40)*params(9)+y(30)*(1-params(9)));
    residual(21) = (y(25)) - (y(46)*params(1)+y(38)*(1-params(10))*(1-params(1)*params(10))/params(10));
    residual(22) = (y(23)) - (params(17)*y(5)+(1-params(17))*(y(25)*params(12)+params(13)*y(22))+y(41));
    residual(23) = (y(41)) - (params(16)*y(8)+x(it_, 3));
    residual(24) = (y(20)) - (y(9)*params(6)+y(12)*params(7));
    residual(25) = (y(21)) - (y(10)*params(6)+y(13)*params(7));
    residual(26) = (y(11)) - (y(9)-y(10));
    residual(27) = (y(14)) - (y(12)-y(13));
    residual(28) = (y(19)) - (y(3)-y(4));
    residual(29) = (y(22)) - (y(20)-y(21));
    residual(30) = (y(28)) - (y(26)-y(27));
    residual(31) = (y(31)) - (y(29)-y(30));
    residual(32) = (y(34)) - (y(32)-y(33));
    residual(33) = (y(37)) - (y(35)-y(36));

end
