function [Press CapillaryPres DisjoiningPres] = calcPress(h, h_store, h_store_new, t_store, x, L_flat, deltaX, kappa);

for i = 2:length(h) - 1
    Pres(i) = -1./(x(i)*deltaX^2)*((x(i) + x(i+1))/2*h(i+1) - ((x(i) + x(i+1))/2 + (x(i) + x(i-1))/2)*h(i) + ...
                        (x(i) + x(i-1))/2*h(i-1)) + 1/(12*kappa)*1/h(i)^3;
end



[val_h_max idx_h_max] = max(h_store([1:round(length(L_flat)/deltaX)],end));
if idx_h_max == 1;
    Pres_h_max = -1./(x(idx_h_max+1)*deltaX^2)*((x(idx_h_max+1) + x(idx_h_max+2))/2*h(idx_h_max+2) - ((x(idx_h_max+1) + x(idx_h_max+2))/2 + (x(idx_h_max+1) + x(idx_h_max))/2)*h(idx_h_max+1) + ...
                    (x(idx_h_max+1) + x(idx_h_max))/2*h(idx_h_max));
else
    Pres_h_max = -1./(x(idx_h_max)*deltaX^2)*((x(idx_h_max) + x(idx_h_max+1))/2*h(idx_h_max+1) - ((x(idx_h_max) + x(idx_h_max+1))/2 + (x(idx_h_max) + x(idx_h_max-1))/2)*h(idx_h_max) + ...
                    (x(idx_h_max) + x(idx_h_max-1))/2*h(idx_h_max-1));
end

for i = 2:length(h) - 1
    for saver = 1:length(t_store)
        Pressure(i,saver) = -1./(x(i)*deltaX^2)*((x(i) + x(i+1))/2*h_store_new(i+1, saver) - ((x(i) + x(i+1))/2 + (x(i) + x(i-1))/2)*h_store_new(i, saver) +...
                        (x(i) + x(i-1))/2*h_store_new(i-1, saver)) + 1/(12*kappa)*1/h_store_new(i, saver)^3;
        CapillaryPres(i,saver) = -1./(x(i)*deltaX^2)*((x(i) + x(i+1))/2*h_store_new(i+1, saver) - ((x(i) + x(i+1))/2 + (x(i) + x(i-1))/2)*h_store_new(i, saver) +...
                        (x(i) + x(i-1))/2*h_store_new(i-1, saver));
        DisjoiningPres(i,saver) = 1/(12*kappa)*1/h_store(i, saver)^3;
    end
end

Press = [Pres' Pressure];

end

