% this function creates a kernel for low pass butterworth filter
function butterfilt = AJbutter(L,cutoff,order,dfreq)
%fnyquist = fsamp/2; %Hz nyquist freq

%create frequency array
% freq = dfreq:dfreq:L*dfreq;
freq = -L/2*dfreq+dfreq:dfreq:L/2*dfreq;%-fnyquist:dfreq:fnyquist;
% size(freq)

butterfilt = 1./ sqrt(1 + (freq/cutoff).^(2*order));

% plot(freq,butterfilt)
end