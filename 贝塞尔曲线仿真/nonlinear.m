function [ceq,c]=nonlinear(x)%������Լ��
[K_max,K_min]=find_M(x);
c=[K_max-0.187;
    -K_min-0.187;
    ];
ceq=[];
end

