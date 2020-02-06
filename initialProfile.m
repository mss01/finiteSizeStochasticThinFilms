function [h x] = initialProfile(kappa,L_flat,L_curv, R_f, Rc,transitionLength,deltaX, filmConfiguration, correctionLP_switch)


switch filmConfiguration
    case 'flatFilms_PBC' 
        x = -2*deltaX:deltaX:L_flat+2*deltaX; % add two ghost points each sides
        h = ones(size(x));
    case 'semiInfiniteNonFlatFilms'
%         x1 = (-L_curv - deltaX):deltaX:-deltaX;    % domain for the curved portion of the film along with the one ghost point
        x1 = -deltaX:-deltaX:(-L_curv - deltaX);    % domain for the curved portion of the film along with the one ghost point
        x1 = fliplr(x1);
        x2 = 0:deltaX:(L_flat+2*deltaX);   % domain for the flat portion of the film along with the two ghost points
        x = [x1 x2];      % full domain
        h1_dom = 1 + kappa*x1.^2; % + 10^-3.*(1 - 2.*rand(1,length(x1)));          % height for the curved portion
        h2_dom = ones(max(size(x2)),1)'; % + 10^-3.*(1 - 2.*rand(1,length(x2)));   % height of the flat potion
        h = [h1_dom h2_dom];              % full domain 
    case 'finiteSizedNonFlatFilms2D'
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
            h = [h_left h_right];
        else
            h = [h_left 1 h_right];
        end
    case 'axisSymmetricFilm'
%         cutOff_r = 3*deltaX;
%         transitionLength = 1.4;
%         FlatPortion = (cutOff_r - 2*deltaX):deltaX:L_flat - transitionLength;
%         actualTransition = (L_flat - transitionLength + deltaX):deltaX:L_flat + transitionLength;
%         if mod(L_flat, deltaX) ~= 0
%             CurvedPortion = (L_flat + transitionLength        ):deltaX:(L_flat + L_curv + deltaX);
%         else
%             CurvedPortion = (L_flat + transitionLength + deltaX):deltaX:(L_flat + L_curv + deltaX);
%         end
%         x = [FlatPortion actualTransition CurvedPortion];
%         hFlat = ones(size(FlatPortion));
%         dLeft =  L_flat - transitionLength;
%         dRight = L_flat + transitionLength;
%         B = -dLeft*(dLeft/2 + 1/(2*transitionLength)*(dLeft^2/3 - (L_flat^2 - transitionLength^2)/2));
%         A = B - dRight^3/(12*transitionLength);
%         D = -(dLeft^2/4 + (dLeft^3/9 - dRight*dLeft^2/4)./(2.*transitionLength) ...
%                 + B*log(dLeft));
%         C = -5*dRight^3/(72*transitionLength) + (B-A)*log(dRight) + D;
%         hTransition = 1 + actualTransition.^2./4 + (actualTransition.^3./9 - (dRight.*actualTransition.^2./4))./(2.*transitionLength) ...
%             + B.*log(actualTransition) + D;
%         switch correctionLP_switch
%             case 'on' % not modified for the transition length for now
%                 hCurv = 1 + (CurvedPortion.^2 - L_flat.^2)./4.*1./(1 - R_f.^2./Rc.^2) + L_flat^2/2.*log(L_flat./CurvedPortion).*1./(1 - R_f.^2./Rc.^2);
%             case 'off'
% %                 hCurv = 1 + (CurvedPortion.^2 - L_flat.^2)./4 + L_flat^2/2.*log(L_flat./CurvedPortion);
%                 hCurv = 1 + CurvedPortion.^2./4 + A.*log(CurvedPortion) + C;
%         end
%         h = [hFlat hTransition hCurv];


% 
%         ampl = 2.5*10^-2;
%         cutOff_r = 3*deltaX;
%         FlatPortion = cutOff_r - 2*deltaX:deltaX:L_flat;
%         if mod(L_flat, deltaX) ~= 0
%             CurvedPortion = (L_flat         ):deltaX:(L_flat + L_curv + deltaX);
%         else
%             CurvedPortion = (L_flat + deltaX):deltaX:(L_flat + L_curv + deltaX);
%         end
%         x = [FlatPortion CurvedPortion];
%         hFlat = ones(size(FlatPortion)) + ampl.*sin(6*pi.*FlatPortion/(L_flat));
%         switch correctionLP_switch
%             case 'on'
%                 hCurv = 1 + (CurvedPortion.^2 - L_flat.^2)./4.*1./(1 - R_f.^2./Rc.^2) + L_flat^2/2.*log(L_flat./CurvedPortion).*1./(1 - R_f.^2./Rc.^2);
%             case 'off'
%                 hCurv = 1 + (CurvedPortion.^2 - L_flat.^2)./4 + L_flat^2/2.*log(L_flat./CurvedPortion);
%         end
%         h = [hFlat hCurv];


        
        cutOff_r = 3*deltaX;
        FlatPortion = cutOff_r - 2*deltaX:deltaX:L_flat;
        if mod(L_flat, deltaX) ~= 0
            CurvedPortion = (L_flat         ):deltaX:(L_flat + L_curv + deltaX);
        else
            CurvedPortion = (L_flat + deltaX):deltaX:(L_flat + L_curv + deltaX);
        end
        x = [FlatPortion CurvedPortion];
        hFlat = ones(size(FlatPortion))*1;
        switch correctionLP_switch
            case 'on'
                hCurv = 1 + (CurvedPortion.^2 - L_flat.^2)./4.*1./(1 - R_f.^2./Rc.^2) + L_flat^2/2.*log(L_flat./CurvedPortion).*1./(1 - R_f.^2./Rc.^2);
            case 'off'
                hCurv = 1 + (CurvedPortion.^2 - L_flat.^2)./4 + L_flat^2/2.*log(L_flat./CurvedPortion);
        end
        h = [hFlat hCurv];
      
end

h = h';
x = x';

end