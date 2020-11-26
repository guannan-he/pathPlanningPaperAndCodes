function [J] = cost_func(x)
%x组成：[求解参数（3），初始点（4），末点（3）]，共10个元素
[min_k,max_k] = get_k(x);
J = max_k - min_k;
end