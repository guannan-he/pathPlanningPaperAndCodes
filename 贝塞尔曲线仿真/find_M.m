function [K_max,K_min]=find_M(x)%�ҵ���������


[K]=calculate_K(x);
% 
% for i=0:0.005:1  
%     num=K(i);
% end
i=0:0.005:1;
num = K(i);

[K_max]=max(num);
[K_min]=min(num);

end