clear all
close all
clc

tic

R_film = 100e-6;
Rc = 1.8e-3;
h0 = 300e-9; 
A_vw = 2e-21;
gamma = 0.034;
visc = 0.001;

kappa = pi*h0^3*gamma/A_vw/Rc;
l_scale = sqrt(2*pi*gamma/A_vw)*h0^2;
t_scale = 36*pi^2*gamma*visc*h0^5/A_vw^2;

R_film_sc = R_film/l_scale;
transLength = sqrt(1/(2*kappa));
R_PB = 15*transLength;

deltaR = 0.005;
ctimestep = 3.25;
endTime = 0.01;
seN = 20;
cutOff_thickness = 2e-5;
deltaT = deltaR^ctimestep;
cutOff_r = deltaR/2;
p2 = deltaT/(deltaR^2);
p1 = deltaT/(deltaR^4);
FlatPortion = cutOff_r - 2*deltaR:deltaR:R_film_sc;
CurvedPortion = R_film_sc+deltaR: deltaR : R_PB + deltaR;
gridR = [FlatPortion CurvedPortion]';
N = length(gridR)-1;
hFlat = ones(size(FlatPortion));
hCurv = 1 + kappa*(CurvedPortion - R_film_sc).^2;
h_init = [hFlat hCurv];
h = h_init';
plot(gridR, h_init)

h_size= N + 1; % extra 1+2 for the ghost point and extra 1 for the additional grid point
A=spalloc( h_size , h_size ,8 + h_size *5-20) ;% preallocate sparse matrix
%% generate q vector that is going to be used to fill up band in the sparse matrix
q1=linspace (3 , h_size-2,h_size-4) *(1 + h_size ) ;
q = [] ;
k = linspace(1,h_size,h_size);
for i =-3:1
    q=[q q1+ i*h_size] ;
end
%% boundary conditions
q_a=[1 h_size+2 h_size*(h_size-2) h_size*(h_size-2)+(h_size-1) h_size*(h_size-1) h_size^2];
A(q_a)=ones(1,6);
q_b=[3*h_size+2 4*h_size+1];
A(q_b)= ones (1,2) *-1;
A(h_size*(h_size-2)) = gridR(h_size - 2)/gridR(h_size - 1);
A(h_size*(h_size-1)) = -(gridR(h_size) + 2*gridR(h_size - 1) + gridR(h_size - 2))/gridR(h_size - 1);
A(h_size^2) = (gridR(h_size) + gridR(h_size - 1))/gridR(h_size - 1);
%% call the solver


k = linspace(1,h_size,h_size);     % vector used in vectorization
flag = 0;                          % counter for storing files
rng shuffle                        % change the seed every new realization

t_range = deltaT:deltaT:endTime;
saver = 1;
for i = 1:length(t_range)      % time marching
    t = i*deltaT;

    h1_r = 2*gridR(k(3:h_size-2)).*h(k(3:h_size-2)).^2.*gridR(k(4:h_size-1)).*h(k(4:h_size-1)).^2./(gridR(k(3:h_size-2)).*h(k(3:h_size-2))+ gridR(k(4:h_size-1)).*h(k(4:h_size-1)));

    h1_l = 2*gridR(k(3:h_size-2)).*h(k(3:h_size-2)).^2.*gridR(k(2:h_size-3)).*h(k(2:h_size-3)).^2./(gridR(k(3:h_size-2)).*h(k(3:h_size-2))+ gridR(k(2:h_size-3)).*h(k(2:h_size-3)));

    A(q) =[h1_l*p2; -(h1_r+3*h1_l)*p2; 3*(h1_r+h1_l)*p2+1; -(3*h1_r+h1_l)*p2; h1_r*p2];
    
    b = [0; 0; h(3:h_size-2) + p2./gridR(k(3:h_size-2)).*((h1_r./h(k(4:h_size-1)).^3) - (h1_r + h1_l)./h(k(3:h_size-2)).^3 + h1_l./h(k(2:h_size-3)).^3); 1+kappa*(gridR(end-1) - R_film_sc)^2; 8*kappa*deltaR^2];
    
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
    
end

t_store(t_store == 0) = [];
h_store(h_store == 0) = [];
h_store = reshape(h_store,length(h), length(t_store));
while min(h_store(:,end)) < 0
    h_store = h_store(:,[1:end-1]);
    t_store = t_store(1:end-1);
end
t_rupt = t;                     % rupture time obtained from this simulation
[minH{1} idx] = min(h_store);


save('hData.mat','h_store','t_store','minH', 't_rupt', '-v7.3');

toc


