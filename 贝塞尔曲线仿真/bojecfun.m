
function f=bojecfun(x) %% cost function

[K_max,K_min]=find_M(x);

f=K_max-K_min;

end






















