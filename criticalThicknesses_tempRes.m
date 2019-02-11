function [h_cr h_cr_left h_cr_right] = criticalThicknesses_tempRes(deltaT, t, h_avg_j, del, cr_thickness, h_right_avg_j, h_left_avg_j, t_cr);

for i = 1:length(t)
    if (t(i) - (t(end) - t_cr)) <= deltaT
        t_crit = t(i);
        h_cr = h_avg_j(i);
        h_cr_right = h_right_avg_j(i);
        h_cr_left = h_left_avg_j(i);
    end
end

% [val_minHavg_j idx_minHavg_j] = min(h_avg_j);

% for i = 1:length(h_avg_j)
%         if del(i) <= cr_thickness && val_minHavg_j <= cr_thickness
%             h_cr_all(i) = h_avg_j(i);
%         elseif del(i) <= cr_thickness && min(h_avg_j) >= cr_thickness
%             h_cr_all_right(i) = h_right_avg_j(i); 
%             h_cr_all_left(i) = h_left_avg_j(i);
%         end
% end
% 
% if min(h_avg_j) < cr_thickness
%     h_cr_all(h_cr_all == 0) = [];
% elseif min(h_avg_j) >= cr_thickness
%     h_cr_all_right(h_cr_all_right == 0) = [];
%     h_cr_all_left(h_cr_all_left == 0) = [];
% end
% 
% if exist('h_cr_all','var') == 1;
%     h_cr = h_cr_all(1);
%     h_cr_right = 0;
%     h_cr_left = 0;
% elseif exist('h_cr_all_right','var') == 1 && exist('h_cr_all_left','var') == 1;
%     h_cr_right = h_cr_all_right(1);
%     h_cr_left = h_cr_all_left(1);
%     h_cr = 0;
% end

end