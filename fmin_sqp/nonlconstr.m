function [c,ceq] = nonlconstr(x)
[min_yd,max_yd] = yd(x);
[min_ydd,max_ydd] = ydd(x);
[min_phi_d,max_phi_d] = phi_d(x);
[min_spd,max_spd] = vel(x);
c = [max_yd - 3
     -min_yd - 3
     max_ydd - 1
     -min_ydd - 1
     max_phi_d - 0.15
     -min_phi_d - 0.15
     max_spd - 19.44
     -min_spd - 0];
 ceq = [];
end