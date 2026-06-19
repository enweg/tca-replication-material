function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = model_extensive_gap_no_investment.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(33, 53);
g1(1,3)=(-1);
g1(1,16)=1;
g1(1,40)=(-1);
g1(2,4)=(-1);
g1(2,18)=1;
g1(2,40)=(-1);
g1(3,6)=(-params(14));
g1(3,39)=1;
g1(3,51)=(-1);
g1(4,7)=(-params(15));
g1(4,40)=1;
g1(4,52)=(-1);
g1(5,9)=1;
g1(5,42)=(-1);
g1(5,23)=T(1);
g1(5,46)=(-T(1));
g1(5,39)=T(1)*(params(14)-1);
g1(6,10)=1;
g1(6,43)=(-1);
g1(6,24)=T(1);
g1(6,39)=T(1)*(params(14)-1);
g1(7,9)=params(3);
g1(7,26)=params(4);
g1(7,29)=(-1);
g1(8,10)=params(3);
g1(8,27)=params(4);
g1(8,30)=(-1);
g1(9,23)=1;
g1(9,46)=(-1);
g1(9,32)=1;
g1(9,47)=(-(params(1)*(1-params(2))));
g1(9,49)=(-(params(1)*params(5)));
g1(10,24)=1;
g1(10,33)=1;
g1(10,48)=(-(params(1)*(1-params(2))));
g1(10,50)=(-(params(1)*params(5)));
g1(11,1)=params(8);
g1(11,12)=(-(params(8)*(1+params(1))));
g1(11,44)=(-(params(8)*(-params(1))));
g1(11,32)=1;
g1(12,2)=params(8);
g1(12,13)=(-(params(8)*(1+params(1))));
g1(12,45)=(-(params(8)*(-params(1))));
g1(12,33)=1;
g1(13,12)=(-params(2));
g1(13,3)=(-(1-params(2)));
g1(13,15)=1;
g1(14,13)=(-params(2));
g1(14,4)=(-(1-params(2)));
g1(14,17)=1;
g1(15,3)=1;
g1(15,20)=(-1);
g1(15,35)=1;
g1(15,38)=(-1);
g1(16,4)=1;
g1(16,21)=(-1);
g1(16,36)=1;
g1(17,20)=(-1);
g1(17,26)=1;
g1(17,29)=1;
g1(17,38)=(-1);
g1(18,21)=(-1);
g1(18,27)=1;
g1(18,30)=1;
g1(19,29)=(-(1-params(9)));
g1(19,35)=(-params(9));
g1(19,38)=1;
g1(19,40)=params(9);
g1(20,30)=(-(1-params(9)));
g1(20,36)=(-params(9));
g1(20,40)=params(9);
g1(21,25)=1;
g1(21,46)=(-params(1));
g1(21,38)=(-((1-params(10))*(1-params(1)*params(10))/params(10)));
g1(22,22)=(-((1-params(17))*params(13)));
g1(22,5)=(-params(17));
g1(22,23)=1;
g1(22,25)=(-((1-params(17))*params(12)));
g1(22,41)=(-1);
g1(23,8)=(-params(16));
g1(23,41)=1;
g1(23,53)=(-1);
g1(24,9)=(-params(6));
g1(24,12)=(-params(7));
g1(24,20)=1;
g1(25,10)=(-params(6));
g1(25,13)=(-params(7));
g1(25,21)=1;
g1(26,9)=(-1);
g1(26,10)=1;
g1(26,11)=1;
g1(27,12)=(-1);
g1(27,13)=1;
g1(27,14)=1;
g1(28,3)=(-1);
g1(28,4)=1;
g1(28,19)=1;
g1(29,20)=(-1);
g1(29,21)=1;
g1(29,22)=1;
g1(30,26)=(-1);
g1(30,27)=1;
g1(30,28)=1;
g1(31,29)=(-1);
g1(31,30)=1;
g1(31,31)=1;
g1(32,32)=(-1);
g1(32,33)=1;
g1(32,34)=1;
g1(33,35)=(-1);
g1(33,36)=1;
g1(33,37)=1;

end
