function [avg_cr_thinningRate_fit intercept_cr_thinningRate avg_cr_thinningRate_mean] = critical_dhdt(t, del, h_critical_start, h_critical_end);

for i = 1:length(del)-1
    if del(i) <= h_critical_start && del(i) >= h_critical_end
        del_rel(i) = del(i);
        t_rel(i) = t(i);
        dhmindt(i) = abs(del(i+1) - del(i))/(t(i+1) - t(i));
    end
end
del_rel(del_rel == 0) = [];
t_rel(t_rel == 0) = [];
dhmindt(dhmindt == 0) = [];

p = polyfit(t_rel, del_rel, 1);
avg_cr_thinningRate_fit = p(1);
avg_cr_thinningRate_mean = mean(dhmindt);
intercept_cr_thinningRate = p(2);

end