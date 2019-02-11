function [locDimple_left locDimple_right x_centre] = keyLocationsInFilm(x, L_flat, deltaX);

for i = 1:length(x)-1
    dum_x(i) = x(i) + x(i+1);
    if (x(i)+L_flat) <= deltaX/2
        locDimple_left = i;
    elseif (x(i)- L_flat) <= deltaX/2
        locDimple_right = i;
    end
    if x(i) == 0
        x_centre = i;
    elseif abs(dum_x(i)) <= deltaX
        x_centre = i;
    end
end

end