function makeAnimation_det(filmConfiguration, animationSkip,kappa, L_flat, L_curv, transitionLength,deltaX, h_store, t_store)

[h x] = initialProfile(kappa,L_flat,L_curv,transitionLength,deltaX, filmConfiguration);
h_all = [h h_store];
t = [0 t_store];                    % corresponds to each time step, including the initial time step                    
formatSpec = 't = %0.3e,';          % to show time steps
t_string2 =  (sprintf(formatSpec,t));
t_new = strsplit(t_string2,',');

%% read the height profiles for every time steps
% for k = 1:1:b   % 1:s:(b-last)
%     str2 = strcat('.\', newFolder, '\', Out{k,1});
%     c(:,k)  = dlmread(str2);
% end
% c(:,all(c==0)) = [];
% c = [c(:,(1:end))];
q = max(size(t));

%% to specify animation conditions

v = VideoWriter('filmEvolution');
v.FrameRate = 5;  % Default 30
v.Quality = 100;    % Default 75
open(v)
hfig = figure;
i = 1;
j = 1;
while i <= q
    Y = h_all(:,i);
    M(:,:,:,j) = video_plot(L_flat, x,Y, t_new(i));
    i = i + animationSkip;
    j = j + 1;
end
% M = M(~cellfun(@isempty,{M.cdata}));
writeVideo(v,M)
close(v)


end