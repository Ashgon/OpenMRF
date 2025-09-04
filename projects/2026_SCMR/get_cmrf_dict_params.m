function [T1, T2, T1p] = get_cmrf_dict_params(T1_lim,T2_lim, T1p_lim, factor, T1_interest, T2_interest, T1p_interest, n_interest)

T1 = T1_lim(1);
while T1(end)<T1_lim(2)
    T1(end+1) = T1(end) * factor;
end
T1(end) = [];

T2 = T2_lim(1);
while T2(end)<T2_lim(2)
    T2(end+1) = T2(end) * factor;
end
T2(end) = [];

T1p = T1p_lim(1);
while T1p(end)<T1p_lim(2)
    T1p(end+1) = T1p(end) * factor;
end
T1p(end) = [];

[T1, T2, T1p] = ndgrid(T1, T2, T1p);
T1  = T1(:);
T2  = T2(:);
T1p = T1p(:);

temp_del = ((T1<T2) + (T1<T1p) + (T1p<T2)) > 0;
T1(temp_del)  = [];
T2(temp_del)  = [];
T1p(temp_del) = [];
clear temp_del;

temp_del =  (T1  >= T1_interest(1)  & T1  <= T1_interest(2)) ...
         &  (T2  >= T2_interest(1)  & T2  <= T2_interest(2)) ...
         &  (T1p >= T1p_interest(1) & T1p <= T1p_interest(2));
T1(temp_del)  = [];
T2(temp_del)  = [];
T1p(temp_del) = [];
clear temp_del;

T1_interest  = linspace(T1_interest(1),  T1_interest(2),  n_interest);
T2_interest  = linspace(T2_interest(1),  T2_interest(2),  n_interest);
T1p_interest = linspace(T1p_interest(1), T1p_interest(2), n_interest);

[T1_interest, T2_interest, T1p_interest] = ndgrid(T1_interest, T2_interest, T1p_interest);
T1_interest  = T1_interest(:);
T2_interest  = T2_interest(:);
T1p_interest = T1p_interest(:);

temp_del = ((T1_interest<=T2_interest) + (T1_interest<=T1p_interest) + (T1p_interest<=T2_interest)) > 0;
T1_interest(temp_del)  = [];
T2_interest(temp_del)  = [];
T1p_interest(temp_del) = [];
clear temp_del;

T1  = [T1; T1_interest];
T2  = [T2; T2_interest];
T1p = [T1p; T1p_interest];

figure();
subplot(1,3,1);
plot(T1, T2, '.');
xlabel('T1 [s]');
ylabel('T2 [s]');
axis square;
subplot(1,3,2);
plot(T1, T1p, '.');
xlabel('T1 [s]');
ylabel('T1p [s]');
axis square;
subplot(1,3,3);
plot(T1p, T2, '.');
xlabel('T1p [s]');
ylabel('T2 [s]');
axis square;

end