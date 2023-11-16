Jvecs = [ 3 3.5 4];% in log
alphas_vec = [ 0.8 1.3 1.8]; % not in log, it is a power - shall we even be log/exponentiating this?
taus_vec = [  1 2 3 ]; % in log

ki = 0;
for ji = 1: 3
    for ai = 1:3
        for ti = 1:3
            ki = ki + 1;
            params(ki,1) = Jvecs(ji);
            params(ki,2) = alphas_vec(ai);
            params(ki,3) = taus_vec(ti);
        end
    end
end

params