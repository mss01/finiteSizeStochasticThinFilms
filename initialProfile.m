function [h x] = initialProfile(kappa,L_flat,L_curv,transitionLength,deltaX)

x1 = [-transitionLength-deltaX:-deltaX:-L_curv-deltaX];
x1 = fliplr(x1);
x2 = [-transitionLength:deltaX:transitionLength];
x3 = [transitionLength+deltaX:deltaX:L_flat];
x_left = [x1 x2 x3] - L_flat;
x_right = -fliplr(x_left);
% x_right = x_left + L_flat + L_curv + deltaX;
if x_left(end) == 0
    x_right(1) = [];
    x = [x_left x_right];
else
    x = [x_left 0 x_right];
end
h1 = 1+ transitionLength^2*kappa/3 + x1.^2*kappa;
h2 = 1 + transitionLength^2*kappa/6 - 0.5*transitionLength*x2*kappa + x2.^2*kappa/2 - x2.^3*kappa/(6*transitionLength); 
h3 = ones(size(x3));

h_left = [h1 h2 h3];
h_right = fliplr(h_left);
if x_left(end) == 0
    h_right(1) = [];
    h = [h_left h_right]';
else
    h = [h_left 1 h_right]';
end
end