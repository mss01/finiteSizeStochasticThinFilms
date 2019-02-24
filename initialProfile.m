function [h x] = initialProfile(kappa,L_flat,L_curv,transitionLength,deltaX, filmConfiguration)


switch filmConfiguration
    case 'flatFilms_PBC' 
        x = 0:deltaX:L_flat;
        h = ones(size(x)); 
    case 'semiInfiniteNonFlatFilms'
        x1 = (-L_curv - deltaX):deltaX:-deltaX;    % domain for the curved portion of the film along with the one ghost point
        x2 = 0:deltaX:(L_flat+2*deltaX);   % domain for the flat portion of the film along with the two ghost points
        x = [x1 x2];      % full domain
        h1_dom = 1 + kappa*x1.^2;          % height for the curved portion
        h2_dom = ones(max(size(x2)),1)';   % height of the flat potion
        h = [h1_dom h2_dom]';              % full domain 
    case 'finiteSizedNonFlatFilms'
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



end