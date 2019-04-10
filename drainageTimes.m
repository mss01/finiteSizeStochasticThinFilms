function [drainageTime drainageTime_right drainageTime_right_rupt beginDrainageTime_right endDrainageTime_right] = drainageTimes(del,t,h_avg_j, h_right_avg_j, h_drain_start, h_drain_end, cr_thickness);
% function [drainageTime drainageTime_right drainageTime_left drainageTime_right_rupt drainageTime_left_rupt] = drainageTimes(del,t,h_avg_j, h_right_avg_j, h_left_avg_j, h_drain_start, h_drain_end, cr_thickness);

for i = 1:length(del)
    if h_avg_j(i) < h_drain_start && h_avg_j(i) > h_drain_end
        t_drain(i) = t(i);
    elseif (h_avg_j(i) < h_drain_start && min(h_avg_j) >= cr_thickness)
        t_drain(i) = t(i);
        t_drain_right = t((h_right_avg_j  < h_drain_start) & (h_right_avg_j > h_drain_end));
%         t_drain_left =  t((h_left_avg_j  < h_drain_start) & (h_left_avg_j > h_drain_end));
    elseif (h_avg_j(i) >= h_drain_start)
        t_drain(i) = 0;
        t_drain_right = t((h_right_avg_j  < h_drain_start) & (h_right_avg_j > h_drain_end));
%         t_drain_left =  t((h_left_avg_j  < h_drain_start) & (h_left_avg_j > h_drain_end));
    end
end

for i = 1:length(h_right_avg_j)
    if h_right_avg_j(i) < h_drain_start
        t_drain_right_rupt = t(h_right_avg_j < h_drain_start);
%         t_drain_left_rupt  = t(h_left_avg_j < h_drain_start);
    end
end

t_drain(t_drain == 0) = [];
if ~isempty(t_drain);
    drainageTime = t_drain(end) - t_drain(1); 
    drainageTime_right = t_drain_right(end) - t_drain_right(1);
%     drainageTime_left = t_drain_left(end) - t_drain_left(1);
    drainageTime_right_rupt = t_drain_right_rupt(end) - t_drain_right_rupt(1);
%     drainageTime_left_rupt = t_drain_left_rupt(end) - t_drain_left_rupt(1);
else
    drainageTime = 0;
    drainageTime_right = t_drain_right(end) - t_drain_right(1);
%     drainageTime_left = t_drain_left(end) - t_drain_left(1);
    drainageTime_right_rupt = t_drain_right_rupt(end) - t_drain_right_rupt(1);
%     drainageTime_left_rupt = t_drain_left_rupt(end) - t_drain_left_rupt(1);
end

beginDrainageTime_right = t_drain_right(1);
endDrainageTime_right = t_drain_right(end);

end