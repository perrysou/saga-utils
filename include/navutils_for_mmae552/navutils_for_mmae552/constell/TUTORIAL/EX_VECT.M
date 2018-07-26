% Function to test LOS routine...                         
% Generate a constellation with four satellites
[constell_elems] = genconst([2 2 6378137+7e6 0.001 pi/4]);

% Convert to ephemeris format
constell_elems(:,8) = ones(size(constell_elems,1),1)*916;
constell_elems(:,9) = ones(size(constell_elems,1),1)*0; 
const_eph = kep2geph(constell_elems);

% Load the GPS almanac and convert to ephemeris format
alm = readyuma;               
%alm = alm(1:2,:);
gps_eph = alm2geph(alm);                            

skip_slow_way = 0;

% Propagate the orbits
t_start = [916 0];
t_stop = [916 100];
dt = 1;
[tc_gps,prnc,xc,vc]=propgeph(const_eph,t_start,t_stop,dt);
[tg_gps,prng,xg,vg]=propgeph(gps_eph,t_start,t_stop,dt);

tc = gpst2sec(tc_gps);
tg = gpst2sec(tg_gps);

% Find the unique times and satellite numbers
t_unique = unique(tg);
num_times = length(t_unique);                
const_sat_nums = unique(prnc);
gps_sat_nums = unique(prng);
num_const = length(const_sat_nums);
num_gps = length(gps_sat_nums);
                              
if skip_slow_way ~= 1
% Time the LOS computations
t0 = clock;
                              
% Initialize the outputs
los_vect = [];
los_t = [];
c_num = [];
g_num = [];

% Compute LOS by looping in time and SVs
for i = 1:num_times                
  for j = 1:num_gps                                        
    % Find the GPS index for this time and satellite number
    Ig = find(tg == t_unique(i) & prng == gps_sat_nums(j));        
    
    % Find all of the GPS satellites at this time
    Ic = find(tc == t_unique(i));
    
    % Just to make it really SLOW
    for k = 1:length(Ic)            
      if ~isempty(Ig)
        this_los = xg(Ig,:) - xc(Ic(k),:);
      
        % Add this LOS, time and sat numbers to the outputs
        los_vect = [los_vect; this_los];
        los_t = [los_t; tg(Ig)];
        c_num = [c_num; prnc(Ic(k))];
        g_num = [g_num; prng(Ig)];
      end
    end
  end
end                     

fprintf('Elapsed time for non-vectorized LOS was %g seconds.\n',etime(clock,t0));
end % if skip_slow_way ~= 1

t0 = clock;
% Compute LOS the new & vectorized way
[t_los_new, los_new, ind] = los(tc_gps,[prnc xc],tg_gps,[prng xg]);

fprintf('Elapsed time for vectorized LOS was %g seconds.\n',etime(clock,t0));

fprintf('There were %d GPS satellites and %d constellation satellites.\n',...
         num_gps,num_const);
fprintf('A total of %d time steps were evaluated.\n',length(unique(t_los_new(:,2))));
fprintf('For a total of %d LOS computations.\n',length(los_new));         

return

% Compute LOS from the output positions 
los_output_pos = xg(ind(:,2),:) - xc(ind(:,1),:);

% Compute LOS from the output indices
los_output_ind = xg(ind(:,2),:) - xc(ind(:,1),:);

if skip_slow_way ~= 1
% Plot the differences in the 3 LOS computations
figure
plot(los_new - los_vect);
title('Difference in vectorized and non-vectorized LOS computations')
ylabel('m')
xlabel('Should be 0')

figure
plot(los_new - los_output_pos);
title('Difference in vectorized and position output from LOS')
ylabel('m')
xlabel('Should be 0')

figure
plot(los_new - los_output_ind);
title('Difference in vectorized and index output from LOS')
ylabel('m')
xlabel('Should be 0')

end
    
