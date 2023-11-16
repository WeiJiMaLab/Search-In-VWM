for ri = 1:30
%output = LL_EVPF(2,alldata(11).data.set_size,alldata(11).data.col_dist',alldata(11).data.response, 1000, [3.6 -3.8 2]);
output = LL_EVPF(2,alldata(11).data.set_size,alldata(11).data.col_dist',alldata(11).data.response, 5000, [3.6 -3.8 2]);
ll_ri(ri) = sum(output);
end

% 1000 samples -- 0.18
% 5000 samples --- 0.0799
%we dont know what would be for other parameter values, so in the code put
% 0.5