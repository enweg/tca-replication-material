function [y, T, residual, g1] = static_7(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(11, 1);
  residual(1)=(params(4)*y(18))-(y(21)-y(1)*params(3));
  residual(2)=(y(24))-(params(1)*params(5)*y(27)+y(24)*params(1)*(1-params(2))-(y(15)-y(17)));
  residual(3)=(y(14))-(y(12)-y(13));
  residual(4)=(y(7))-(params(2)*y(4)+y(7)*(1-params(2)));
  residual(5)=(y(27))-(y(30)+y(12)-y(7));
  residual(6)=(y(21))-(y(30)+y(12)-y(18));
  residual(7)=(y(30))-(y(27)*params(9)-y(32)*params(9)+y(21)*(1-params(9)));
  residual(8)=(y(17))-(y(17)*params(1)+y(30)*(1-params(10))*(1-params(1)*params(10))/params(10));
  residual(9)=(y(15))-(y(15)*params(17)+(1-params(17))*(y(17)*params(12)+params(13)*y(14))+y(33));
  residual(10)=(y(1))-(y(1)-1/params(3)*(y(15)-y(17))-y(31)*1/params(3)*(params(14)-1));
  residual(11)=(y(12))-(y(1)*params(6)+y(4)*params(7));
if nargout > 3
    g1_v = NaN(31, 1);
g1_v(1)=params(4);
g1_v(2)=1;
g1_v(3)=(-(params(1)*params(5)));
g1_v(4)=1;
g1_v(5)=(-params(9));
g1_v(6)=1;
g1_v(7)=(-((1-params(17))*params(13)));
g1_v(8)=(-params(2));
g1_v(9)=(-params(7));
g1_v(10)=1-(1-params(2));
g1_v(11)=1;
g1_v(12)=(-1);
g1_v(13)=(-1);
g1_v(14)=(-1);
g1_v(15)=1;
g1_v(16)=(-1);
g1_v(17)=1;
g1_v(18)=(-(1-params(9)));
g1_v(19)=(-1);
g1_v(20)=(-1);
g1_v(21)=1;
g1_v(22)=(-((1-params(10))*(1-params(1)*params(10))/params(10)));
g1_v(23)=(-1);
g1_v(24)=1-params(1);
g1_v(25)=(-((1-params(17))*params(12)));
g1_v(26)=(-(1/params(3)));
g1_v(27)=1;
g1_v(28)=1-params(17);
g1_v(29)=1/params(3);
g1_v(30)=params(3);
g1_v(31)=(-params(6));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 11, 11);
end
end
