function [J] = cost_func(x)
%x��ɣ�[��������3������ʼ�㣨4����ĩ�㣨3��]����10��Ԫ��
[min_k,max_k] = get_k(x);
J = max_k - min_k;
end