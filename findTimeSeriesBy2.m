function [idxTimeStamp] = findTimeSeriesBy2(t_store)

t_vector = t_store;
tes = length(t_vector);
tes2 = pow2(floor(log2(tes)));

N_power = log(tes2)/log(2);

for i = 1:N_power
    idxTimeStamp(i) = ceil(length(t_vector)/2);
    t_vector = t_store(1:idxTimeStamp(i));
end
idxTimeStamp = [length(t_store) idxTimeStamp 1];

end