function [K]=calculate_K(x) %º∆À„«˙¬ 

d1=x(1);
d4=x(2);
y2=x(3);
syms u;


%% X
X=(-20)*(1-u)*u^3-5*u^4;
du_X=diff(X,u);
ddu_X=diff(du_X,u);

%% Y
Y=4*d1*(1-u)^3*u+6*y2*(1-u)^2*u^2+4*(10-d4)*(1-u)*u^3+10*u^4;
du_Y=diff(Y,u);
ddu_Y=diff(du_Y,u);

K_p=(du_X*ddu_Y-du_Y*ddu_X)/(du_X^2+du_Y^2)^(3/2);
K=matlabFunction(K_p);

end