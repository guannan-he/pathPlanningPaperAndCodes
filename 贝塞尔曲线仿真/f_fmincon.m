%This is the main function;
%Using SQP algorithm to optimize the parameter of the quintic polynomial;

function []=f_fmincon()
x0=[1 1 5];%初始值
A=[-1,0,0;%不等式约束
    0,-1,0;
    0,0,0];
b=[0;0;0];


options=optimoptions('fmincon','Display','iter','Algorithm','sqp');%SQP
x=fmincon(@bojecfun,x0,A,b,[],[],[],[],@nonlinear,options);

u=0:0.005:1;
d1=x(1);
d4=x(2);
y2=x(3);
X=(-20).*(1-u).*u.^3-5.*u.^4;
Y=4*d1*(1-u).^3.*u+6*y2*(1-u).^2.*u.^2+4*(10-d4)*(1-u).*u.^3+10*u.^4;
plot(X,Y)


end
