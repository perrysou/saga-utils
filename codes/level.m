function [tau_a] = level(z,meancc,lag0,lag,dt)
for ncorr = lag0+1:length(lag)
    if meancc(1,1,ncorr)<z
        tau_a = (z-meancc(1,1,ncorr))...
            /(meancc(1,1,ncorr-1)-meancc(1,1,ncorr));
        tau_a = dt*(lag(ncorr-1)-tau_a);
    end 
end
end

