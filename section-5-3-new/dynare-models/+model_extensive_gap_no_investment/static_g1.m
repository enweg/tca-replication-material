function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = model_extensive_gap_no_investment.static_g1_tt(T, y, x, params);
end
g1 = zeros(33, 33);
g1(1,7)=(-1);
g1(1,8)=1;
g1(1,32)=(-1);
g1(2,9)=(-1);
g1(2,10)=1;
g1(2,32)=(-1);
g1(3,31)=1-params(14);
g1(4,32)=1-params(15);
g1(5,15)=T(1);
g1(5,17)=(-T(1));
g1(5,31)=T(1)*(params(14)-1);
g1(6,16)=T(1);
g1(6,31)=T(1)*(params(14)-1);
g1(7,1)=params(3);
g1(7,18)=params(4);
g1(7,21)=(-1);
g1(8,2)=params(3);
g1(8,19)=params(4);
g1(8,22)=(-1);
g1(9,15)=1;
g1(9,17)=(-1);
g1(9,24)=1-params(1)*(1-params(2));
g1(9,27)=(-(params(1)*params(5)));
g1(10,16)=1;
g1(10,25)=1-params(1)*(1-params(2));
g1(10,28)=(-(params(1)*params(5)));
g1(11,24)=1;
g1(12,25)=1;
g1(13,4)=(-params(2));
g1(13,7)=1-(1-params(2));
g1(14,5)=(-params(2));
g1(14,9)=1-(1-params(2));
g1(15,7)=1;
g1(15,12)=(-1);
g1(15,27)=1;
g1(15,30)=(-1);
g1(16,9)=1;
g1(16,13)=(-1);
g1(16,28)=1;
g1(17,12)=(-1);
g1(17,18)=1;
g1(17,21)=1;
g1(17,30)=(-1);
g1(18,13)=(-1);
g1(18,19)=1;
g1(18,22)=1;
g1(19,21)=(-(1-params(9)));
g1(19,27)=(-params(9));
g1(19,30)=1;
g1(19,32)=params(9);
g1(20,22)=(-(1-params(9)));
g1(20,28)=(-params(9));
g1(20,32)=params(9);
g1(21,17)=1-params(1);
g1(21,30)=(-((1-params(10))*(1-params(1)*params(10))/params(10)));
g1(22,14)=(-((1-params(17))*params(13)));
g1(22,15)=1-params(17);
g1(22,17)=(-((1-params(17))*params(12)));
g1(22,33)=(-1);
g1(23,33)=1-params(16);
g1(24,1)=(-params(6));
g1(24,4)=(-params(7));
g1(24,12)=1;
g1(25,2)=(-params(6));
g1(25,5)=(-params(7));
g1(25,13)=1;
g1(26,1)=(-1);
g1(26,2)=1;
g1(26,3)=1;
g1(27,4)=(-1);
g1(27,5)=1;
g1(27,6)=1;
g1(28,7)=(-1);
g1(28,9)=1;
g1(28,11)=1;
g1(29,12)=(-1);
g1(29,13)=1;
g1(29,14)=1;
g1(30,18)=(-1);
g1(30,19)=1;
g1(30,20)=1;
g1(31,21)=(-1);
g1(31,22)=1;
g1(31,23)=1;
g1(32,24)=(-1);
g1(32,25)=1;
g1(32,26)=1;
g1(33,27)=(-1);
g1(33,28)=1;
g1(33,29)=1;

end
