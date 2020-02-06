function makeAnimation_det(filmConfiguration, correctionLP_switch, animationSkip,kappa, L_flat, L_curv, R_f, Rc, transitionLength,deltaX, ...
                            h_store, t_store, h0_init, t_scale, l_scale, beginDrainageTime_right, endDrainageTime_right, x_dimple_loc_right, res_limit, h_drain_start, h_drain_end)

[h x] = initialProfile(kappa,L_flat,L_curv, R_f, Rc, transitionLength,deltaX, filmConfiguration, correctionLP_switch);
h_all = [h h_store];
t = [0 t_store]*t_scale;                    % corresponds to each time step, including the initial time step                    
% t = [0 t_store];                    % corresponds to each time step, including the initial time step                    
formatSpec = 't = %0.2f,';          % to show time steps
t_string2 =  (sprintf(formatSpec,t));
t_string3 = strsplit(t_string2,',');
t_new = strcat(t_string3, 's');

%% read the height profiles for every time steps
% for k = 1:1:b   % 1:s:(b-last)
%     str2 = strcat('.\', newFolder, '\', Out{k,1});
%     c(:,k)  = dlmread(str2);
% end
% c(:,all(c==0)) = [];
% c = [c(:,(1:end))];
q = max(size(t));

%% to specify animation conditions

v = VideoWriter('filmEvolution_dimensional');
v.FrameRate = 5;  % Default 30
v.Quality = 100;    % Default 75
open(v)
hfig = figure;
i = 1;
j = 1;
plot(x'*l_scale*10^6,h_all(:,1)*h0_init*10^9);
while i <= q
    Y = h_all(:,i);
    M(:,:,:,j) = video_plot(L_flat, x,Y, t_new(i), h0_init, t(i), t_scale, l_scale, beginDrainageTime_right, ...
                            endDrainageTime_right, x_dimple_loc_right(i), res_limit, deltaX, h_drain_start, h_drain_end);
    i = i + animationSkip;
    j = j + 1;
end
M(:,:,:,end) = video_plot(L_flat, x,h_all(:,end), t_new(end), h0_init, t(end), t_scale, l_scale, beginDrainageTime_right, ...
                            endDrainageTime_right, x_dimple_loc_right(end), res_limit, deltaX, h_drain_start, h_drain_end);
% M = M(~cellfun(@isempty,{M.cdata}));
writeVideo(v,M)
close(v)


end