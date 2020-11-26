function [c_ieq_n,b_ieq_n] = ieq_constr(x)
b_ieq_n = [0.187;0.187];
[min_k,max_k] = get_k(x);
c_ieq_n = [max_k;-min_k];
end