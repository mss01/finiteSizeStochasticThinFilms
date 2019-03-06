function [t_rupt x_rupt minH] = filmSolver(filmConfiguration, disjPress_switch, L_flat,transitionLength,L_curv,N,deltaX,deltaT,kappa,Tmp,gx,h_adjusted,A,p,endTime,seN, cutOff_thickness);


tic

%% equation
%                                                           ___
%  d H      d  /   3 d^3 H     1 dH   \      __    d  /    /           \
%  --- = - ---|   H  ----- +   -----  |  + _/2T   ---|    / H^3 N(X,T)  |
%  d T     d X \     d X^3     H dX   /           d X \ \/             / 

%% boundary conditions

% at X --> -L_curv, H = 1 + \kappa X^2;
% at X --> -L_curv, d^2H/dX^2 =  2 \kappa;

% at X --> L_flat, dH/dX = 0
% at X --> L_flat, d^3H/dX^3 = 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now here is the final discretization
%  _              
% |         j+1                        j+1                          j+1  
% | p2*h2* H        -p2*(3*h2 + h1)* H        ( 3*p2*(h1+h2) + 1 ) H     ... 
% |_         i-2                        i-1                          i
%                                          _
%                     j+1             j+1   |                
%   -p2*(3*h1 + h2)* H        p2*h1* H      |    =  LHS --> A matrix
%                     i+1             i+2  _|                       


% _                                                                               _          
%| 2*kappa*deltaX^2                                                                | 
%| kappa*x(2)^2+1                                                                  | 
%|                                                                                 | 
%|previous        disj pressure                   noise term                       |     
%| soln             vdW only                                                       |
%|  j                                                                              |
%| H       -     p1 (h2_r - h2_l)  +    p3*(sqrt(h1_r).*noi1 - sqrt(h1_l).*noi2)
%|  i                                                                              |        
%|                                                                                 |  
%|0                                                                                |  
%|_0                                                                              _| 



%%%%%----------- Recall the initial conditions ------------%%%%%

[h x] = initialProfile(kappa, L_flat,L_curv,transitionLength,deltaX, filmConfiguration);

if isequal(disjPress_switch, 'on') 
    p1 = deltaT/(deltaX^2);            % parameter for the explicit part (disj press)
elseif isequal(disjPress_switch, 'off') 
    p1 = 0;
end

%if disjPress_switch == 'on';
 %   p1 = deltaT/(deltaX^2);            % parameter for the explicit part (disj press)
%elseif disjPress_switch == 'off'
%    p1 = 0;
%end
p2 = deltaT/(deltaX^4);            % parameter for the implicit part (surf tension)
p3 = 1/(deltaX)*sqrt(2*deltaT*Tmp);     % parameter for the noise term, also explicit

%%  preallocation

h1_r = zeros(size(h));
h1_l = zeros(size(h));
h2_r = zeros(size(h));
h2_l = zeros(size(h));
R = zeros(size(x));
noi_r = zeros(size(x));
noi_l = zeros(size(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = linspace(1,h_adjusted,h_adjusted);     % vector used in vectorization
flag = 0;                          % counter for storing files
rng shuffle                        % change the seed every new realization

%%%% preallocating h_storage and t_storage. But it results into memory issues if you want to run the simulation until endTime of a large value
% h_store = zeros(length(h), floor(endTime/(seN*deltaT)));   
% t_store = zeros(1,floor(endTime/(seN*deltaT)));
t_range = deltaT:deltaT:endTime;
saver = 1;
for i = 1:length(t_range)      % time marching
    t = i*deltaT;

    R = randn(N+1,1);              % random numbers for each grid point (changes every time)
    gx_f = gx.*R;                  % use the gx matrix as per the model
    g = sum(gx_f);                 % Q-Wiener process definition - gives the noise term
    g = g';      % extending to the boundary points
    noi_r = (g(3:h_adjusted-2) + g(4:h_adjusted-1))./2;   % for the main domain
    noi_l = (g(2:h_adjusted-3) + g(3:h_adjusted-2))./2;   % for the main domain
%     noi1 = (g(2:h_size-2) + g(3:h_size-1))./2;
%     noi2 = (g(1:h_size-3) + g(2:h_size-2))./2;


    h1_r = 2*h(k(3:h_adjusted-2)).^2.*h(k(4:h_adjusted-1)).^2./(h(k(3:h_adjusted-2))+h(k(4:h_adjusted-1)));

    h1_l = 2*h(k(3:h_adjusted-2)).^2.*h(k(2:h_adjusted-3)).^2./(h(k(3:h_adjusted-2))+h(k(2:h_adjusted-3)));

    h2_r =  0.5.*(1./h(k(3:h_adjusted-2)) + 1./h(k(4:h_adjusted-1))).*(h(k(4:h_adjusted-1)) - h(k(3:h_adjusted-2)));

    h2_l =  0.5.*(1./h(k(3:h_adjusted-2)) + 1./h(k(2:h_adjusted-3))).*(h(k(3:h_adjusted-2)) - h(k(2:h_adjusted-3)));

    A(p) =[h1_l*p2; -(h1_r+3*h1_l)*p2; 3*(h1_r+h1_l)*p2+1; -(3*h1_r+h1_l)*p2; h1_r*p2];

    %% generate the b-vector in Ax = b
    switch filmConfiguration
        case 'flatFilms_PBC'
%             if Tmp == 0
%                 b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l); 0; 0;];
%             else
%                 b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi_r - sqrt(h1_l).*noi_l); 0; 0];
%             end
                b = [0; 0; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi_r - sqrt(h1_l).*noi_l); 0; 0];

        case 'semiInfiniteNonFlatFilms'
%             if Tmp == 0
%                 b = [2*kappa*deltaX^2; kappa*x(2)^2+1; h(3:h_adjusted-2)-p1*(h2_r - h2_l); 0; 0];
%             else
%                 b = [2*kappa*deltaX^2; kappa*x(2)^2+1; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi_r - sqrt(h1_l).*noi_l); 0; 0];
%             end
                b = [2*kappa*deltaX^2; kappa*x(2)^2+1; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi_r - sqrt(h1_l).*noi_l); 0; 0];
        case 'finiteSizedNonFlatFilms'
%             if Tmp == 0
%                 b = [2*kappa*deltaX^2; kappa*(x(2)+L_flat)^2+1; h(3:h_adjusted-2)-p1*(h2_r - h2_l); kappa*(x(end-1)-L_flat)^2+1 ; 2*kappa*deltaX^2];
%             else
%                 b = [2*kappa*deltaX^2; kappa*(x(2)+L_flat)^2+1; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi_r - sqrt(h1_l).*noi_l); kappa*(x(end-1)-L_flat)^2+1 ; 2*kappa*deltaX^2];
%             end
                b = [2*kappa*deltaX^2; kappa*(x(2)+L_flat)^2+1; h(3:h_adjusted-2)-p1*(h2_r - h2_l) + p3.*(sqrt(h1_r).*noi_r - sqrt(h1_l).*noi_l); kappa*(x(end-1)-L_flat)^2+1 ; 2*kappa*deltaX^2];
    end
    
    h = A\b;

%======================================================================
%               Save data after every few time steps
%======================================================================
    if i == length(t_range)
        error('Beware ---> endTime is not sufficient')
        break;
    end
    flag = flag + 1;
    if(flag==seN)
        h_store(:,saver) = h;
        t_store(saver) = t;
        flag=0;
%         dlmwrite(sprintf('Data_%1.15f.txt',t),h,'precision','%.16f','delimiter','\t')    % save it as a txt file to be used in matlab post processing script
        saver = saver + 1;
    end
    if min(h(:)) <= cutOff_thickness
%         dlmwrite(sprintf('Data_%1.15f.txt',t),h,'precision','%.16f','delimiter','\t')      % write the last file because it becomes important especially for high kappa values in the late regime
        h_store(:,saver) = h;
        t_store(saver) = t;
        break              % if the film height goes below a certain height stop the realization
    else
        continue
    end
    

        
%  some check to ensure if any realization has given singular matrix, Matlab gives this warning away anyway, so have commented it for now.       
%     if(rcond(full(A)) < 1e-12)
%         break
%     else
%         continue
%     end
end



% test1 = size(h_store)
% test2 = length(h)
% test3 = length(t_store)
t_store(t_store == 0) = [];
h_store(h_store == 0) = [];
h_store = reshape(h_store,length(h), length(t_store));
while min(h_store(:,end)) < 0
    h_store = h_store(:,[1:end-1]);
    t_store = t_store(1:end-1);
end
t_rupt = t;                     % rupture time obtained from this simulation
[minH{1} idx] = min(h_store);
x_min = x(idx);
x_rupt = x_min(end);

save('hData.mat','h_store','t_store','minH', 'x_min', 't_rupt', 'x_rupt', '-v7.3');

toc

end