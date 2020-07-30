function [ceq,c]=nonlinear(x)%非线性约束
[K_max,K_min]=find_M(x);
c=[K_max-0.187;
    -K_min-0.187;
    ];
ceq=[];
end

