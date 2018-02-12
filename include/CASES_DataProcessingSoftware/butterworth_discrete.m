
% Equations used from:
% http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
% Cookbook formulae for audio EQ biquad filter coefficients
% by Robert Bristow-Johnson <rbj@audioimagination.com>

function outdata = butterworth_discrete(data, Fs, freqcut, order,...
    filter_type)

if mod(order,2) ~= 0 
    disp('This function can only do even-ordered filters')
end

if strcmp(filter_type,'lp') == 1,
    lp = 1;
elseif strcmp(filter_type,'hp') == 1,
    hp = 1;
else
    disp('Try writing lp/ hp for high pass or low pass butterworth filter')
end

    
n_prime = 20;

n_fudge = 0;
%n_fudge = round(Fs/(4.*freqcut))

%     Set up the output array
num = length(data);
outdata = zeros(1,num);

% For improved stability, the high-order filter will be constructed
%  as cascading 2nd order biquad filters.  We loop through order/2 and
%  apply each 2nd order biquad separately
for k=0:1:(order/2-1),
    if k == 0 
      temp = data;
    else
      temp = outdata;
    end

%      Pad the data to "prime" the IIR filter.  We're using n_prime
%      points at the beginning of the signal, plus n_fudge
%      additional points at the end.  The n_prime points are used to
%      get the filter going so it doesn't go unstable when it first
%      encounters data.  The n_fudge points at the end are
%      used to calibrate the time shift of the discrete filter
%      compared to a frequency-domain butterworth filter with the same
%      parameters, the output array will be shifted by these
%      positions, IE: for n_fudge=2
%      [input[0:100],input[100],input[100]] becomes
%      [trash,trash,output[0:100]]
%      Unfortunately an IIR filter has a frequency-dependent time
%      shift, so this can't be completely eliminated, but it can at
%      least be minimised in the pass-band
            
      tempdata = temp;  
      if n_prime > 0,
%           size(tempdata)
%           size(repmat(tempdata(1),1,n_prime))
          tempdata = [repmat(tempdata(1),1,n_prime), tempdata];
      end
      
      if n_fudge > 0,
          tempdata = [tempdata, repmat(tempdata(end-1),1,n_fudge)];
      end
         
%     Initialize outdata to the input array to prime the IIR filter
    outdata = tempdata;
    
%      Set up the filter coefficients
    coeff_a = zeros(1,3);
    coeff_b = zeros(1,3);
    
    w0 = 2*pi*(freqcut/Fs);
    
%     Since we're cascading multiple biquad filters, we must adjust
%     the Q for each of the biquads so the final result is a
%     butterworth (Q=.707).  Equation taken from
%     http://en.wikipedia.org/wiki/Butterworth_filter
    q = ((-2)*cos((2*(order/2.0-k)+order-1.0)/(2.0*order)*pi))^(-1);
    
    alpha = sin(w0)/(2*q);

%    Set up the coefficients for the filter
    if lp ==1, 
        coeff_b(1) = (1 - cos(w0))/2;
        coeff_b(2) = 1 - cos(w0);
        coeff_b(3) = coeff_b(1);
        coeff_a(1) = 1 + alpha;
        coeff_a(2) = -2*cos(w0);
        coeff_a(3) = 1 - alpha;
    elseif hp ==1,
        coeff_b(1) =  (1 + cos(w0))/2;
        coeff_b(2) = -(1 + cos(w0));
        coeff_b(3) =  coeff_b(1);
        coeff_a(1) =  1 + alpha;
        coeff_a(2) = -2*cos(w0);
        coeff_a(3) =  1 - alpha;
    end
    
%   Normalize the coefficients
    coeff_b = coeff_b/coeff_a(1);
    coeff_a = coeff_a/coeff_a(1);
    
%   Loop through the signal, applying the filter
    for i=3:1:num-1,
        temp1 = coeff_b(1)*tempdata(i) + coeff_b(2)*tempdata(i-1) +...
            coeff_b(3)*tempdata(i-2);
        temp2 = coeff_a(2)*outdata(i-1) + coeff_a(3)*outdata(i-2);
        
        outdata(i) = temp1 - temp2;
    end
    
%     size(outdata)
%  Remove the prime and fudge indices
    outdata = outdata((n_prime+n_fudge):end);
    
end %for k

end 

