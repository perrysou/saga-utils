%% Fn_Plot_HighRate_CASESdata.m 
% This code is a part of CASES_Post_Processing_Software v1
% Copyright (C) 2012 - 2014, Kshitija Deshpande, Virginia Tech. 

% CASES_Post_Processing_Software v1 is free software: you can redistribute
% it and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% 
% CASES_Post_Processing_Software v1 is distributed in the hope that it will
% be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% ----------------------------------------------------------------------

% %% Fn_Plot_HighRate_CASESdata.m: Produces different types of high rate
% plots for each PRN.

% % Inputs: 
% 1. Needs to change "path" to enter the correct path of the CASES
% folder and folder name where all the CASES log files are stored.
% folder_path is the complete folder path including the folder name to be
% input to the read high rate data MATLAB function
% (Fn_ReadHighRate_CASESdata) and plotting high rate data function
% (Fn_Plot_HighRate_CASESdata). 
% 2. signal type:
% See Fn_ReadHighRate_CASESdata.m for details on signal type.
% Generally, it is 0 or 2 for AAL-PIP CASES. Thus, i=1 or 3. 
% 3. types of plots
% A = all & scintillating events, all segments
% of data -- processed
% B = S4 and Sigma_phi plots
% C = plot common scintillating times for PRNs
% 
% % Outputs: 
% PRN_files_<SignalType>/FilteredData_PRN<PRN>.mat with processed phase, 
% power data as a function of time for each PRN <PRN> and <SignalType>
% are saved.



function Fn_Plot_HighRate_CASESdata(folder_name,signal_type,set_plot)

%File separator for the selected operating system
sep = filesep;
tic

%Specify plot type: plot_type = 1 for eps, anything else for png
plot_type = 0;

set(0,'defaultlinelinewidth',1.8)
%set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

%Turn warnings off
warning off  all

%Specify signal type
%  signal_type = 0; 
switch signal_type
   case 0
      signal = 'L1CA';
   case 1
      signal = 'L2CM';
   case 2
      signal = 'L2CL';
   case 3
      signal = 'L2CLM';
   case 4
      signal = 'L5I';
   case 5
      signal = 'L5Q';
   case 6
      signal = 'L5IQ';
   case 7
      signal = 'L1CA-ALT1';
   case 8
      signal = 'CDMA-UHF-PILOT';
   case 9
      signal = 'CDMA-UHF-SYNC';
   otherwise
      error('Unknown signal.')
end

%specify types of plots 
% A = all & scintillating events, all segments of
%data -- processed 
% B = S4 and Sigma_phi plots
% C = plot common scintillating times for PRNs
%  set_plot = 'B';

%Specify path to save plots and other files
command = strcat('mkdir',{' '},folder_name,'plots_',signal);
system(cell2mat(command));
command = strcat('mkdir',{' '},folder_name,'plots_',signal,sep,...
    'HighRate_ScintAll',sep);
system(cell2mat(command));
command = strcat('mkdir',{' '},folder_name,'PRN_files_',signal);
system(cell2mat(command));
FiltDatadir = strcat(folder_name,'PRN_files_',signal,sep);
outdir = strcat(folder_name,'plots_',signal,sep,'HighRate_ScintAll',sep);

infilename = strcat(folder_name,'PRN_files_',signal,sep,...
    'Processed_HR_Data.mat'); 
load(infilename);

DATA = DATA;

%RGB color values
light_pink= [1, 0.71, 0.76];
hot_pink= [1,0.42,0.71];
light_skyblue= [0.53,0.81,0.98];
peacock_blue= [0.2, 0.63, 0.79];
forest_green= [0.133, 0.545, 0.133];

%Plotting intialization
fontsz = 12;
marksz = 4;
linestyle= '.';
fontsz_lg = 10; %font size of legends

%% Plots set A: Plot for each PRN with each segment zoomed in one all
%processed data
if strcmp(set_plot,'A') == 1,
%Plot scintillating events after Rx clock removal in all available PRNs
for kk = 1:1:32,
    sv = (DATA(:,8)==kk);
    svdata = DATA(sv,:);
    
    if(isempty(svdata)~=1),
        
                
        ort_iq = svdata(:,3); %ort seconds + fractional seconds
        % %Find and remove repeated entries
        %1st derivative to find the time slots (wherever there is a huge
        %difference, that's the BUG. There are repeated entries inside
        %iq.log)
        ort_iqd = diff(ort_iq);
        iqdn_ind=find(ort_iqd<0);
        for ii = 1:1:length(iqdn_ind),
            indii = find(ort_iq == ort_iq(iqdn_ind(ii)));
            ort_iq(indii(1):indii(end)-1)=0;
            %This will work well only if there are only 2 repetitions max.
            %ort_sc = ort_sc(ort_sc~=0) will remove the repeated values.
            %But better leave the zeros as they are and recreate the
            %original "scdata" array based on those zeros.
        end
        %Non-zero time stamps means without any repetition.
        nzero_ort = find(ort_iq~=0);
        svdata = svdata(nzero_ort,:);
        %Reread time for scint file(with no repeated entries).
        orts =svdata(:,3);
        
        %length of main data
        LM = length(svdata);
        
        %read all IQ data
        obstime_o = svdata(:,1); %original obstime
        ortw = svdata(:,2); % ORT weeks
%         orts = svdata(:,3); % ORT seconds
        pcdata = svdata(:,4); % processed carrier phase values
        piqphdata = svdata(:,5); % processed iq phase values
        piqpowdata = svdata(:,6); % processed iq power values
        refprn = svdata(:,7); %PRN of reference channel
        scintprn = svdata(:,8); % PRN for scintillating channel
        
        %gps to utc conversion
        [UTC_time, leap_sec, day_of_year]=gps2utc([ortw orts], 0);
        %UTC_time = year month day hour minutes seconds
        UTC_init = UTC_time(1,:);
        year = UTC_init(1);
        month = UTC_init(2);
        day = UTC_init(3);
        hour = UTC_init(4);
        
        %Observation time seconds after UTC_init hours:
        obstime = (UTC_time(:,4)-hour)*3600 + UTC_time(:,5)*60 ...
            + UTC_time(:,6);        
        
        % simply plot the high rate IQ data without the scintillating
        % events

        %process only if considerable amount of data present
        if length(obstime)>1000,
                    
        %Initialize the stacking variables for this PRN
        OBSTIME = zeros(1);
        PCDATA = zeros(1);
        IQPOW = zeros(1);
        IQPH = zeros(1);

        % % Zoom certain patches in time
        %find if there indeed are data chunks with such a time gap
        difference = 1;
        % Difference in two data 'patches' in seconds in the prn file
        diff_ind = find(diff(obstime) >= difference);

        fignum = 0;
        if(isempty(diff_ind)==0)
            for zmt = 1:1:length(diff_ind)+1, % finding segments in time
                %for all the high rate data

                if(zmt == 1)                 % first segment
                    zm_st = 1;
                else
                    zm_st = diff_ind(zmt-1)+1;
                end

                if(zmt == length(diff_ind)+1) % last segment
                    zm_end = LM;
                else
                    zm_end = diff_ind(zmt);
                end

                obstimez = obstime(zm_st:zm_end);
                fcdataz = pcdata(zm_st:zm_end);
                fiqpowz = piqpowdata(zm_st:zm_end);
                fiqphz = piqphdata(zm_st:zm_end);
                
                %process only if considerable amount of data present
                if length(obstimez)>1000,
                
                %length of each continuous segment in time
                LT = length(obstimez);
                % If there exists segments in phase data itself, handle
                % those first. Find discontinuities in phase data.
                rawphase = fiqphz;
                indphtemp = find( abs(diff(diff(rawphase))) > 5);

                % % Diagnosis
%                 figure
%                 plot(diff(diff(rawphase)))
                if isempty(indphtemp) ~= 1,
                    indph = indphtemp(2:2:end);
                    
%                     disp('Satellite PRN:')
%                     disp(kk)
%                     disp('IQ Phase discontinuities at time:')
%                     disp(obstime(indph))
%                     disp('After UT hours')
%                     disp(hour)
%                     disp('orts')
%                     disp(orts(indph))
                    
                    for zmph = 1:1:length(indph)+1,
                        %for all the high rate data
                        %  obstimez = obstime((diff_ind(zm))*(zm-1)+1:...
                        %   zm*diff_ind(zm));
                        if(zmph == 1)                 % first segment
                            zm_sts = 1;
                        else
                            zm_sts = indph(zmph-1)+1;
                        end

                        if(zmph == length(indph)+1) % last segment
                            zm_ende = LT;
                        else
                            zm_ende = indph(zmph);
                        end
                        fignum = fignum +1;

%                         disp(zmph)
%                         disp(obstime(zm_st))
%                         disp(obstime(zm_end))

                        %Phase differenced zoomed in variables
                        obstimezphd = obstimez(zm_sts:zm_ende);
                        fcdatazphd = fcdataz(zm_sts:zm_ende);
                        fiqpowzphd = fiqpowz(zm_sts:zm_ende);
                        fiqphzphd = fiqphz(zm_sts:zm_ende);
                        
                        %process only if considerable amount of data present
                        if length(obstimezphd)>1000,
                        
                        %Post processing goes here for the longer continuous
                        %lengths of data
                        % Polyfit subtraction
                        degree = 3;
                        rawphase = fiqphzphd;
                        %polynomial fitting and subtraction
                        poly_coef = polyfit(obstimezphd,...
                            rawphase,degree);
                        poly = polyval(poly_coef,obstimezphd);
                        poly_sub = rawphase - poly;
                        phaseph = detrend(poly_sub);

                        degree = 3;
                        rawphase = fcdatazphd;
                        %polynomial fitting and subtraction
                        poly_coef = polyfit(obstimezphd,...
                            rawphase,degree);
                        poly = polyval(poly_coef,obstimezphd);
                        poly_sub = rawphase - poly;
                        fcdatazphd = detrend(poly_sub);

                        % %AJ filtering of phase
                        % %Filter phase in frequency domain
                        tlength = obstimezphd(end) - obstimezphd(1);
                        dfreq = 1/tlength;
                        % disp(tlength)
                        % disp(dfreq)
                        Lzph = length(phaseph);
                        % NFFT = 2^nextpow2(L);
                        % fft of phase data
                        wdata_fft = fftshift(fft(ifftshift(phaseph))); 
                        cutoff = 0.1; %[Hz]
                        order = 6; %order of the butterworth filter
                        fsamp = 50; %[Hz]
                        %create a low pass butterworth filter Kernel
                        butterlow = AJbutter(Lzph,cutoff,order,dfreq);
                        butterhi = 1.0 - butterlow;
                        wdata_filt = ...
                            fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                        phase_fftfilt = real(wdata_filt);
                        fiqphzphd = phase_fftfilt;

                        %AJ way for processing of power
                        Tdata = obstimezphd(end) - obstimezphd(1);
                        if Tdata <=60,
                            fc = 0.2;
                        else
                            fc = 0.1; %Desired cut off freq in Hz
                        end
                        %sampling frequency
                        fsamp = round(1/(obstimezphd(4)-obstimezphd(3)));
                        norder = 6; %order of the butterworth filter
                        power = fiqpowzphd';
                        low_pass_power = ...
                            butterworth_discrete(power,fsamp,fc,...
                            norder,'lp');
                        low_pass_power = ...
                            low_pass_power(1:length(power));
                        power1 = power./low_pass_power;
                        fiqpowzphd = power1';
                        
                        edgeT = 10; %edge time in seconds to be thrown out
                        edsamp = edgeT * fsamp;
                        obstimezphd = obstimezphd(edsamp+1:end-edsamp);
                        fiqphzphd = fiqphzphd(edsamp+1:end-edsamp);
                        fcdatazphd = fcdatazphd(edsamp+1:end-edsamp);
                        fiqpowzphd = fiqpowzphd(edsamp+1:end-edsamp);

                        %plot these zoomed in segments
                        %Plot Carrier Phase
                        figure('Visible','Off')
                        set(gca,'FontSize',fontsz)
                        plot(obstimezphd,fcdatazphd*2*pi,'r',...
                            'MarkerSize',marksz,...
                            'linestyle',linestyle);
                        grid on
                        str = strcat('Detrended, clock error corrected carrier',...
                            {' (zoomed in) phase of the signal: '},signal,...
                            ', PRN:',num2str(kk),{', fig'},num2str(fignum));
                        title(str)
                        xstring = strcat({'Time [s] after '},num2str(hour),...
                            ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                            ,num2str(day),'/',num2str(year));
                        xlabel(xstring);
                        ylabel('Carrier Phase [rad]')
                        
                        str_name = strcat(signal,'_CorrCarrPh_PRN',num2str(kk),...
                            'zoom',num2str(fignum));
                        if(plot_type == 1)
                            plotfile = strcat(outdir,str_name,'.eps');
                            saveas(gcf,plotfile,'epsc2');
                        else
                            plotfile = strcat(outdir,str_name,'.png');
                            saveas(gca,plotfile);
                        end
                        clear gca
                        close

                        %     Plot IQ power and phase
                        figure('Visible','Off')
                        subplot(2,1,1)
                        set(gca,'FontSize',fontsz)
                        %To have both the powers overlapped
                        plot(obstimezphd,10*log10(fiqpowzphd),...
                            'b','MarkerSize',marksz);  %if power-mean(power)
                        grid on
                        str = strcat('CASES Detrended power (zoomed in)',...
                            {' for:'}, signal,', PRN:',num2str(kk),...
                            {', fig'},num2str(fignum));
                        title(str)
                        ylabel('Power [dB]')
                        axis([-inf inf -5 5]);
                        subplot(2,1,2)
                        set(gca,'FontSize',fontsz)
                        plot(obstimezphd,fiqphzphd,'g-','MarkerSize',marksz);
                        grid on
                        str = strcat('CASES Detrended & clock error',...
                            {' corrected phase (zoomed in) '},...
                            'for: ', signal,', PRN:',num2str(kk),...
                            {', fig'},num2str(fignum));
                        title(str)
                        xstring = strcat({'Time [s] after '},num2str(hour),...
                            ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                            ,num2str(day),'/',num2str(year));
                        xlabel(xstring);
                        ylabel('Phase [rad]')
                        if max(abs(fiqphzphd))>2,
                            axis([-inf inf -5 5]);
                        else
                            axis([-inf inf -2 2]);
                        end
                        str_name = strcat(signal,'_CorrIQ_PRN',num2str(kk),...
                            'zoom',num2str(fignum));
                        if(plot_type == 1)
                            plotfile = strcat(outdir,str_name,'.eps');
                            saveas(gcf,plotfile,'epsc2');
                        else
                            plotfile = strcat(outdir,str_name,'.png');
                            saveas(gca,plotfile);
                        end
                        clear gca
                        close

                        %save into the Mega array: Stitch into a bigger array
                        OBSTIME = [OBSTIME, obstimezphd'];
                        PCDATA = [PCDATA, fcdatazphd'];
                        IQPOW = [IQPOW, fiqpowzphd'];
                        IQPH = [IQPH, fiqphzphd'];

                        end% process only if considerable data present
                        
                    end %for zoom: zmph of all phase segments

                else %if no discontinutites in phase exist:

                    fignum = fignum + 1;

                    %Post processing goes here for the longer continuous
                    %lengths of data
                    % Polyfit subtraction
                    degree = 3;
                    rawphase = fiqphz;
                    %polynomial fitting and subtraction
                    poly_coef = polyfit(obstimez,...
                        rawphase,degree);
                    poly = polyval(poly_coef,obstimez);
                    poly_sub = rawphase - poly;
                    phase = detrend(poly_sub);

                    degree = 3;
                    rawphase = fcdataz;
                    %polynomial fitting and subtraction
                    poly_coef = polyfit(obstimez,...
                        rawphase,degree);
                    poly = polyval(poly_coef,obstimez);
                    poly_sub = rawphase - poly;
                    fcdataz = detrend(poly_sub);

                    % %AJ filtering of phase
                    % %Filter phase in frequency domain
                    tlength = obstimez(end) - obstimez(1);
                    dfreq = 1/tlength;
                    % disp(tlength)
                    % disp(dfreq)
                    Lz = length(phase);
                    % NFFT = 2^nextpow2(L);
                    wdata_fft = fftshift(fft(ifftshift(phase))); % fft of phase data
                    cutoff = 0.1; %[Hz]
                    order = 6; %order of the butterworth filter
                    fsamp = 50; %[Hz]
                    %create a low pass butterworth filter Kernel
                    butterlow = AJbutter(Lz,cutoff,order,dfreq);
                    butterhi = 1.0 - butterlow;
                    wdata_filt = fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                    phase_fftfilt = real(wdata_filt);
                    fiqphz = phase_fftfilt;
                    
                    %AJ way for processing of power
                    Tdata = obstimez(end) - obstimez(1);
                    if Tdata <=60,
                        fc = 0.2;
                    else
                        fc = 0.1; %Desired cut off freq in Hz
                    end
                    %sampling frequency
                    fsamp = round(1/(obstimez(4)-obstimez(3)));
                    norder = 6; %order of the butterworth filter
                    power = fiqpowz';
                    low_pass_power = ...
                        butterworth_discrete(power,fsamp,fc,...
                        norder,'lp');
                    low_pass_power = ...
                        low_pass_power(1:length(power));
                    power1 = power./low_pass_power;
                    fiqpowz = power1';
                    
                    edgeT = 10; %edge time in seconds to be thrown out
                    edsamp = edgeT * fsamp;
                    obstimez = obstimez(edsamp+1:end-edsamp);
                    fiqphz = fiqphz(edsamp+1:end-edsamp);
                    fcdataz = fcdataz(edsamp+1:end-edsamp);
                    fiqpowz = fiqpowz(edsamp+1:end-edsamp);

                    %plot these zoomed in segments
                    %Plot Carrier Phase
                    figure('Visible','Off')
                    set(gca,'FontSize',fontsz)
                    plot(obstimez,fcdataz*2*pi,'r','MarkerSize',marksz,...
                        'linestyle',linestyle);
                    grid on
                    str = strcat('Detrended, clock error corrected carrier',...
                        {' (zoomed in) phase of the signal: '},signal,...
                        ', PRN:',num2str(kk),{', fig'},num2str(fignum));
                    title(str)
                    xstring = strcat({'Time [s] after '},num2str(hour),...
                        ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                        ,num2str(day),'/',num2str(year));
                    xlabel(xstring);
                    ylabel('Carrier Phase [rad]')
                    %                         if (max(fcdataz)>10)
                    %                             axis([-inf inf -1 1]);
                    %                         end
                    %                         axis([-inf inf -2.5 2.5]);
                    str_name = strcat(signal,'_CorrCarrPh_PRN',num2str(kk),...
                        'zoom',num2str(fignum));
                    if(plot_type == 1)
                        plotfile = strcat(outdir,str_name,'.eps');
                        saveas(gcf,plotfile,'epsc2');
                    else
                        plotfile = strcat(outdir,str_name,'.png');
                        saveas(gca,plotfile);
                    end
                    clear gca
                    close

                    %     Plot IQ power and phase
                    figure('Visible','Off')
                    subplot(2,1,1)
                    set(gca,'FontSize',fontsz)
                    %To have both the powers overlapped
                    plot(obstimez,10*log10(fiqpowz),...
                        'b','MarkerSize',marksz);  %if power-mean(power)
                    grid on
                    str = strcat('CASES Detrended power (zoomed in)',...
                        {' for:'}, signal,', PRN:',num2str(kk),...
                        {', fig'},num2str(fignum));
                    title(str)
                    ylabel('Power [dB]')
                    axis([-inf inf -5 5]);
                    subplot(2,1,2)
                    set(gca,'FontSize',fontsz)
                    plot(obstimez,fiqphz,'g-','MarkerSize',marksz);
                    grid on
                    str = strcat('CASES Detrended & clock error',...
                        {' corrected phase (zoomed in) '},...
                        'for: ', signal,', PRN:',num2str(kk),...
                        {', fig'},num2str(fignum));
                    title(str)
                    xstring = strcat({'Time [s] after '},num2str(hour),...
                        ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                        ,num2str(day),'/',num2str(year));
                    xlabel(xstring);
                    ylabel('Phase [rad]')
                    if max(abs(fiqphz))>2,
                            axis([-inf inf -5 5]);
                        else
                            axis([-inf inf -2 2]);
                        end
                    str_name = strcat(signal,'_CorrIQ_PRN',num2str(kk),...
                        'zoom',num2str(fignum));
                    if(plot_type == 1)
                        plotfile = strcat(outdir,str_name,'.eps');
                        saveas(gcf,plotfile,'epsc2');
                    else
                        plotfile = strcat(outdir,str_name,'.png');
                        saveas(gca,plotfile);
                    end
                    clear gca
                    close

                    %save into the Mega array: Stitch into a bigger array
                    OBSTIME = [OBSTIME, obstimez'];
                    PCDATA = [PCDATA, fcdataz'];
                    IQPOW = [IQPOW, fiqpowz'];
                    IQPH = [IQPH, fiqphz'];

                    
                end % if discontinuities exist in the phase data
                end %process only if length(obstimez)>1000,
            end %for zoom: zmt of all segments

            % Save and Plot Stitched zoomed figures
            obstime = OBSTIME(2:end);
            pcdata = PCDATA(2:end);
            piqpowdata = IQPOW(2:end);
            piqphdata = IQPH(2:end);

            %save filtered data
            data_PRN = [obstime;piqpowdata;piqphdata];
            strprn = strcat('FilteredData_PRN',num2str(kk),'.mat');
            filename = strcat(FiltDatadir,strprn);
            save(filename,'data_PRN');
            strprn = strcat('FilteredData_PRN',num2str(kk),'.txt');
            filename = strcat(FiltDatadir,strprn);
%             csvwrite(filename,data_PRN');
            %Write with higher precision to capture time with 0.01s
            %resolution.
            dlmwrite(filename,data_PRN','Precision','%.10g',...
                'delimiter',',');
            
            %Plot Carrier Phase
            figure('Visible','Off')
            set(gca,'FontSize',fontsz)
            plot(obstime,pcdata*2*pi,'r','MarkerSize',marksz,'linestyle',...
                linestyle);
            grid on
            str = strcat({'Detrended, clock error corrected carrier '},...
                'phase of the signal: ',signal,', PRN:',num2str(kk));
            title(str)
            xstring = strcat({'Time [s] after '},num2str(hour),...
                ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                ,num2str(day),'/',num2str(year));
            xlabel(xstring);
            ylabel('Carrier Phase [rad]')
            %                 axis([-inf inf -2.5 2.5 ]);
            str_name = strcat(signal,'_CorrCarrPh_PRN',num2str(kk));
            if(plot_type == 1)
                plotfile = strcat(outdir,str_name,'.eps');
                saveas(gcf,plotfile,'epsc2');
            else
                plotfile = strcat(outdir,str_name,'.png');
                saveas(gca,plotfile);
            end
            clear gca
            close

            %Plot IQ power and phase
            figure('Visible','Off')
            subplot(2,1,1)
            set(gca,'FontSize',fontsz)
            plot(obstime,10*log10(piqpowdata),...
                'b-','MarkerSize',marksz);
            grid on
            str = strcat('Detrended IQ power and detrended, clock error',...
                {' corrected IQ phase for '},signal,', PRN:',num2str(kk));
            %str = strcat('Detrended IQ power and HPF phase of the signal: ',...
            % signal,', PRN:',num2str(kk));
            %str = strcat('Detrended, polyfit subtracted IQ power and ',...
            % ' phase of the signal: ',signal,', PRN:',num2str(kk));
            title(str)
            axis([-inf inf -5 5]);
            xstring = strcat({'Time [s] after '},num2str(hour),...
                ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                ,num2str(day),'/',num2str(year));
            xlabel(xstring);
            ylabel('Power [dB]')
            subplot(2,1,2)
            set(gca,'FontSize',fontsz)
            plot(obstime,piqphdata,'g-','MarkerSize',marksz);
            grid on
            xstring = strcat({'Time [s] after '},num2str(hour),...
                ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                ,num2str(day),'/',num2str(year));
            xlabel(xstring);
            ylabel('Phase [rad]')
            if max(abs(piqphdata))>2,
                axis([-inf inf -5 5]);
            else
                axis([-inf inf -2 2]);
            end
            str_name = strcat(signal,'_CorrIQ_PRN',num2str(kk));
            if(plot_type == 1)
                plotfile = strcat(outdir,str_name,'.eps');
                saveas(gcf,plotfile,'epsc2');
            else
                plotfile = strcat(outdir,str_name,'.png');
                saveas(gca,plotfile);
            end
            clear gca
            close

        else % no discontinuities in time
             

            % If no discontinuities exist, process here
            % Polyfit subtraction
            %disp('continuous segment of time data')

            %Initialize the stacking variables for this PRN
            OBSTIME = zeros(1);
            PCDATA = zeros(1);
            IQPOW = zeros(1);
            IQPH = zeros(1);

            %find discontinuities in phase data
            rawphase = piqphdata;
            indphtemp = find( abs(diff(diff(rawphase))) > 5);
            fignum = 0;
            if isempty(indphtemp) ~= 1,
                indph = indphtemp(2:2:end);
                
%                 disp('Satellite PRN:')
%                 disp(kk)
%                 disp('IQ Phase discontinuities at time:')
%                 disp(obstime(indph))
%                 disp('After UT hours')
%                 disp(hour)
%                 disp('orts')
%                 disp(orts(indph))
                
%                 figure
%                 plot(diff(diff(rawphase)))
                
                for zmph = 1:1:length(indph)+1,
                    %for all the high rate data

                    if(zmph == 1)                 % first segment
                        zm_st = 1;
                    else
                        zm_st = indph(zmph-1)+1;
                    end

                    if(zmph == length(indph)+1) % last segment
                        zm_end = LM;
                    else
                        zm_end = indph(zmph);
                    end

                    fignum = fignum +1;

%                     disp(zmph)
%                     disp(obstime(zm_st))
%                     disp(obstime(zm_end))

                    obstimez = obstime(zm_st:zm_end);
                    fcdataz = pcdata(zm_st:zm_end);
                    fiqpowz = piqpowdata(zm_st:zm_end);
                    fiqphz = piqphdata(zm_st:zm_end);

                    %process only if considerable amount of data present
                    if length(obstimez)>1000,
                    
                    %Post processing goes here for the longer continuous
                    %lengths of data
                    % Polyfit subtraction
                    degree = 3;
                    rawphase = fiqphz;
                    %polynomial fitting and subtraction
                    poly_coef = polyfit(obstimez,...
                        rawphase,degree);
                    poly = polyval(poly_coef,obstimez);
                    poly_sub = rawphase - poly;
                    phase = detrend(poly_sub);
                    
                    
                    degree = 3;
                    rawphase = fcdataz;
                    %polynomial fitting and subtraction
                    poly_coef = polyfit(obstimez,...
                        rawphase,degree);
                    poly = polyval(poly_coef,obstimez);
                    poly_sub = rawphase - poly;
                    fcdataz = detrend(poly_sub);

                    % %AJ filtering of phase
                    % %Filter phase in frequency domain
                    tlength = obstimez(end) - obstimez(1);
                    dfreq = 1/tlength;
                    % disp(tlength)
                    % disp(dfreq)
                    Lz = length(phase);
                    % NFFT = 2^nextpow2(L);
                    % fft of phase data
                    wdata_fft = fftshift(fft(ifftshift(phase))); 
                    cutoff = 0.1; %[Hz]
                    order = 6; %order of the butterworth filter
                    fsamp = 50; %[Hz]
                    %create a low pass butterworth filter Kernel
                    butterlow = AJbutter(Lz,cutoff,order,dfreq);
                    butterhi = 1.0 - butterlow;
                    wdata_filt = ...
                        fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                    phase_fftfilt = real(wdata_filt);
                    fiqphz = phase_fftfilt;
                    
                                    
                    %AJ way for processing of power
                    Tdata = obstimez(end) - obstimez(1);
                    if Tdata <=60,
                        fc = 0.2;
                    else
                        fc = 0.1; %Desired cut off freq in Hz
                    end
                    %sampling frequency
                    fsamp = round(1/(obstimez(4)-obstimez(3)));
                    norder = 6; %order of the butterworth filter
                    power = fiqpowz';
                    low_pass_power = ...
                        butterworth_discrete(power,fsamp,fc,...
                        norder,'lp');
                    low_pass_power = ...
                        low_pass_power(1:length(power));
                    power1 = power./low_pass_power;
                    fiqpowz = power1';
                    
                    edgeT = 10; %edge time in seconds to be thrown out
                    edsamp = edgeT * fsamp;
                    obstimez = obstimez(edsamp+1:end-edsamp);
                    fiqphz = fiqphz(edsamp+1:end-edsamp);
                    fcdataz = fcdataz(edsamp+1:end-edsamp);
                    fiqpowz = fiqpowz(edsamp+1:end-edsamp);

                    %plot these zoomed in segments
                    %Plot Carrier Phase
                    figure('Visible','Off')
                    set(gca,'FontSize',fontsz)
                    plot(obstimez,fcdataz*2*pi,'r','MarkerSize',marksz,...
                        'linestyle',linestyle);
                    grid on
                    str = strcat('Detrended, clock error corrected carrier',...
                        {' (zoomed in) phase of the signal: '},signal,...
                        ', PRN:',num2str(kk),{', fig'},num2str(fignum));
                    title(str)
                    xstring = strcat({'Time [s] after '},num2str(hour),...
                        ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                        ,num2str(day),'/',num2str(year));
                    xlabel(xstring);
                    ylabel('Carrier Phase [rad]')

                    str_name = strcat(signal,'_CorrCarrPh_PRN',num2str(kk),...
                        'zoom',num2str(fignum));
                    if(plot_type == 1)
                        plotfile = strcat(outdir,str_name,'.eps');
                        saveas(gcf,plotfile,'epsc2');
                    else
                        plotfile = strcat(outdir,str_name,'.png');
                        saveas(gca,plotfile);
                    end
                    clear gca
                    close

                    %     Plot IQ power and phase
                    figure('Visible','Off')
                    subplot(2,1,1)
                    set(gca,'FontSize',fontsz)
                    %To have both the powers overlapped
                    plot(obstimez,10*log10(fiqpowz),...
                        'b','MarkerSize',marksz);  %if power-mean(power)
                    grid on
                    str = strcat('CASES Detrended power (zoomed in)',...
                        {' for:'}, signal,', PRN:',num2str(kk),...
                        {', fig'},num2str(fignum));
                    title(str)
                    ylabel('Power [dB]')
                    axis([-inf inf -5 5]);
                    subplot(2,1,2)
                    set(gca,'FontSize',fontsz)
                    plot(obstimez,fiqphz,'g-','MarkerSize',marksz);
                    grid on
                    str = strcat('CASES Detrended & clock error',...
                        {' corrected phase (zoomed in) '},...
                        'for: ', signal,', PRN:',num2str(kk),...
                        {', fig'},num2str(fignum));
                    title(str)
                    xstring = strcat({'Time [s] after '},num2str(hour),...
                        ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                        ,num2str(day),'/',num2str(year));
                    xlabel(xstring);
                    ylabel('Phase [rad]')
                    if max(abs(fiqphz))>2,
                        axis([-inf inf -5 5]);
                    else
                        axis([-inf inf -2 2]);
                    end
                    str_name = strcat(signal,'_CorrIQ_PRN',num2str(kk),...
                        'zoom',num2str(fignum));
                    if(plot_type == 1)
                        plotfile = strcat(outdir,str_name,'.eps');
                        saveas(gcf,plotfile,'epsc2');
                    else
                        plotfile = strcat(outdir,str_name,'.png');
                        saveas(gca,plotfile);
                    end
                    clear gca
                    close

                    %save into the Mega array: Stitch into a bigger array
                    OBSTIME = [OBSTIME, obstimez'];
                    PCDATA = [PCDATA, fcdataz'];
                    IQPOW = [IQPOW, fiqpowz'];
                    IQPH = [IQPH, fiqphz'];

                    end %process only if length(obstimez)>1000,
                end %for zoom: zmph of all segments

                % Plot Stitch zoomed figures
                obstime = OBSTIME(2:end);
                pcdata = PCDATA(2:end);
                piqpowdata = IQPOW(2:end);
                piqphdata = IQPH(2:end);

                %save filtered data
                data_PRN = [obstime;piqpowdata;piqphdata];
                strprn = strcat('FilteredData_PRN',num2str(kk),'.mat');
                filename = strcat(FiltDatadir,strprn);
                save(filename,'data_PRN');
                strprn = strcat('FilteredData_PRN',num2str(kk),'.txt');
                filename = strcat(FiltDatadir,strprn);
%                 csvwrite(filename,data_PRN');
                dlmwrite(filename,data_PRN','Precision','%.10g',...
                    'delimiter',',');

                %Plot Carrier Phase
                figure('Visible','Off')
                set(gca,'FontSize',fontsz)
                plot(obstime,pcdata*2*pi,'r','MarkerSize',marksz,'linestyle',...
                    linestyle);
                grid on
                str = strcat({'Detrended, clock error corrected carrier '},...
                    'phase of the signal: ',signal,', PRN:',num2str(kk));
                title(str)
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ylabel('Carrier Phase [rad]')
                %                 axis([-inf inf -2.5 2.5 ]);
                str_name = strcat(signal,'_CorrCarrPh_PRN',num2str(kk));
                if(plot_type == 1)
                    plotfile = strcat(outdir,str_name,'.eps');
                    saveas(gcf,plotfile,'epsc2');
                else
                    plotfile = strcat(outdir,str_name,'.png');
                    saveas(gca,plotfile);
                end
                clear gca
                close

                %Plot IQ power and phase
                figure('Visible','Off')
                subplot(2,1,1)
                set(gca,'FontSize',fontsz)
                plot(obstime,10*log10(piqpowdata),...
                    'b-','MarkerSize',marksz);
                grid on
                str = strcat('Detrended IQ power and detrended, clock error',...
                    {' corrected IQ phase for '},signal,', PRN:',num2str(kk));
                title(str)
                axis([-inf inf -5 5]);
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ylabel('Power [dB]')
                subplot(2,1,2)
                set(gca,'FontSize',fontsz)
                plot(obstime,piqphdata,'g-','MarkerSize',marksz);
                grid on
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ylabel('Phase [rad]')
                if max(abs(piqphdata))>2,
                    axis([-inf inf -5 5]);
                else
                    axis([-inf inf -2 2]);
                end
                str_name = strcat(signal,'_CorrIQ_PRN',num2str(kk));
                if(plot_type == 1)
                    plotfile = strcat(outdir,str_name,'.eps');
                    saveas(gcf,plotfile,'epsc2');
                else
                    plotfile = strcat(outdir,str_name,'.png');
                    saveas(gca,plotfile);
                end
                clear gca
                close


            else
                % If there are no discontinuities in phase data

                degree = 3;
                rawphase = piqphdata;
                %polynomial fitting and subtraction
                poly_coef = polyfit(obstime,...
                    rawphase,degree);
                poly = polyval(poly_coef,obstime);
                poly_sub = rawphase - poly;
                phase = detrend(poly_sub);

                degree = 3;
                rawphase = pcdata;
                %polynomial fitting and subtraction
                poly_coef = polyfit(obstime,...
                    rawphase,degree);
                poly = polyval(poly_coef,obstime);
                poly_sub = rawphase - poly;
                pcdata = detrend(poly_sub);

                % %AJ filtering of phase
                % %Filter phase in frequency domain
                tlength = obstime(end) - obstime(1);
                dfreq = 1/tlength;
                % disp(tlength)
                % disp(dfreq)
                Lz = length(phase);
                % NFFT = 2^nextpow2(L);
                % fft of phase data
                wdata_fft = fftshift(fft(ifftshift(phase))); 
                cutoff = 0.1; %[Hz]
                order = 6; %order of the butterworth filter
                fsamp = 50; %[Hz]
                %create a low pass butterworth filter Kernel
                butterlow = AJbutter(Lz,cutoff,order,dfreq);
                butterhi = 1.0 - butterlow;
                wdata_filt = fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                phase_fftfilt = real(wdata_filt);
                piqphdata = phase_fftfilt;

                %AJ way for processing of power
                Tdata = obstime(end) - obstime(1);
                if Tdata <=60,
                    fc = 0.2;
                else
                    fc = 0.1; %Desired cut off freq in Hz
                end
                %sampling frequency
                fsamp = round(1/(obstime(4)-obstime(3)));
                norder = 6; %order of the butterworth filter
                power = piqpowdata';
                low_pass_power = ...
                    butterworth_discrete(power,fsamp,fc,...
                    norder,'lp');
                low_pass_power = ...
                    low_pass_power(1:length(power));
                power1 = power./low_pass_power;
                piqpowdata = power1';
                
                edgeT = 10; %edge time in seconds to be thrown out
                edsamp = edgeT * fsamp;
                obstime = obstime(edsamp+1:end-edsamp);
                piqphdata = piqphdata(edsamp+1:end-edsamp);
                pcdata = pcdata(edsamp+1:end-edsamp);
                piqpowdata = piqpowdata(edsamp+1:end-edsamp);

                %save filtered data:
                data_PRN = [obstime,piqpowdata,piqphdata]';
                strprn = strcat('FilteredData_PRN',num2str(kk),'.mat');
                filename = strcat(FiltDatadir,strprn);
                save(filename,'data_PRN');
                strprn = strcat('FilteredData_PRN',num2str(kk),'.txt');
                filename = strcat(FiltDatadir,strprn);
%                 csvwrite(filename,data_PRN');
                dlmwrite(filename,data_PRN','Precision','%.10g',...
                    'delimiter',',');
                
                %Plot Carrier Phase
                figure('Visible','Off')
                set(gca,'FontSize',fontsz)
                plot(obstime,pcdata*2*pi,'r','MarkerSize',marksz,'linestyle',...
                    linestyle);
                grid on
                str = strcat({'Detrended, clock error corrected carrier '},...
                    'phase of the signal: ',signal,', PRN:',num2str(kk));
                title(str)
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ylabel('Carrier Phase [rad]')
                %                 axis([-inf inf -2.5 2.5 ]);
                str_name = strcat(signal,'_CorrCarrPh_PRN',num2str(kk));
                if(plot_type == 1)
                    plotfile = strcat(outdir,str_name,'.eps');
                    saveas(gcf,plotfile,'epsc2');
                else
                    plotfile = strcat(outdir,str_name,'.png');
                    saveas(gca,plotfile);
                end
                clear gca
                close

                %Plot IQ power and phase
                figure('Visible','Off')
                subplot(2,1,1)
                set(gca,'FontSize',fontsz)
                plot(obstime,10*log10(piqpowdata),...
                    'b-','MarkerSize',marksz);
                grid on
                str = strcat('Detrended IQ power and detrended, clock error',...
                    {' corrected IQ phase for '},signal,', PRN:',num2str(kk));
                %str = strcat('Detrended IQ power and HPF phase of the signal: ',...
                % signal,', PRN:',num2str(kk));
                %str = strcat('Detrended, polyfit subtracted IQ power and ',...
                % ' phase of the signal: ',signal,', PRN:',num2str(kk));
                title(str)
                axis([-inf inf -5 5]);
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ylabel('Power [dB]')
                subplot(2,1,2)
                set(gca,'FontSize',fontsz)
                plot(obstime,piqphdata,'g-','MarkerSize',marksz);
                grid on
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ylabel('Phase [rad]')
                if max(abs(piqphdata))>2,
                    axis([-inf inf -5 5]);
                else
                    axis([-inf inf -2 2]);
                end
                str_name = strcat(signal,'_CorrIQ_PRN',num2str(kk));
                if(plot_type == 1)
                    plotfile = strcat(outdir,str_name,'.eps');
                    saveas(gcf,plotfile,'epsc2');
                else
                    plotfile = strcat(outdir,str_name,'.png');
                    saveas(gca,plotfile);
                end
                clear gca
                close

            end % if zoomed figures based on phase data
        end% if plotting zoomed figures based on time data
        end %process only if length(obstimez)>1000,
    end %isempty svdata
    %close all

end %for kk -- PRN 
% End Set A
end % set_plot == A

%% Plots set B: All Processed data: Plot S4, sigma_phi
if strcmp(set_plot,'B') == 1,
%Create necessary directories
command = strcat('mkdir',{' '},folder_name,'plots_',signal,sep,...
    'HighRate_S4_SPhi',sep);
system(cell2mat(command));
outdir = strcat(folder_name,'plots_',signal,sep,'HighRate_S4_SPhi',sep);

%% Gather data for all PRNs
ORTS = DATA(:,3);  %save all seconds to find out the length of data
L = length(DATA);

% Mega PRN arrays to carry information for all PRNs
OBSTIME_MPRN = zeros(L,32); % In UT
S4A_MPRN = zeros(L,32);
SPHA_MPRN = zeros(L,32);


for kk = 1:1:32,
    sv = (DATA(:,8)==kk);
    svdata = DATA(sv,:);
    
    if(isempty(svdata)~=1),
                
        ort_iq = svdata(:,3); %ort seconds + fractional seconds
        % %Find and remove repeated entries
        %1st derivative to find the time slots (wherever there is a huge
        %difference, that's the BUG. There are repeated entries inside
        %iq.log)
        ort_iqd = diff(ort_iq);
        iqdn_ind=find(ort_iqd<0);
        for ii = 1:1:length(iqdn_ind),
            indii = find(ort_iq == ort_iq(iqdn_ind(ii)));
            ort_iq(indii(1):indii(end)-1)=0;
            %This will work well only if there are only 2 repetitions max.
            %ort_sc = ort_sc(ort_sc~=0) will remove the repeated values.
            %But better leave the zeros as they are and recreate the
            %original "scdata" array based on those zeros.
        end
        %Non-zero time stamps means without any repetition.
        nzero_ort = find(ort_iq~=0);
        svdata = svdata(nzero_ort,:);
        %Reread time for scint file(with no repeated entries).
        orts =svdata(:,3);
        
        LM = length(svdata);
        
        %read all IQ data
        obstime_o = svdata(:,1); %original obstime
        ortw = svdata(:,2); % ORT weeks
        pcdata = svdata(:,4); % processed carrier phase values
        piqphdata = svdata(:,5); % processed iq phase values
        piqpowdata = svdata(:,6); % processed iq power values
        refprn = svdata(:,7); %PRN of reference channel
        scintprn = svdata(:,8); % PRN for scintillating channel
        
        %gps to utc conversion
        [UTC_time, leap_sec, day_of_year]=gps2utc([ortw orts], 0);
        %UTC_time = year month day hour minutes seconds
        UTC_init = UTC_time(1,:);
        year = UTC_init(1);
        month = UTC_init(2);
        day = UTC_init(3);
        hour = UTC_init(4);
        
        %Observation time seconds after UTC_init hours:
        obstime = (UTC_time(:,4)-hour)*3600 + UTC_time(:,5)*60 ...
            + UTC_time(:,6);
%         %absolute UTC observation time vector in seconds
%         obstime_abs = (UTC_time(:,4))*3600 + UTC_time(:,5)*60 ...
%             + UTC_time(:,6);

        %Observation time seconds (absolute) in UT:
    	abs_obstime = (UTC_time(:,4))*3600 + ...
            UTC_time(:,5)*60 + UTC_time(:,6);
        
        %Initialize the stacking variables for this PRN
        AOBSTIME = zeros(1);
        S4DATA = zeros(1);
        %sigma phi variables
        SPHDATA = zeros(1);
        
        % % Zoom certain patches in time
        %find if there indeed are data chunks with such a time gap
        difference = 5;
        % Difference in two data 'patches' in seconds in the prn file
        diff_ind = find(diff(obstime) > difference);

        fignum = 0;
        if(isempty(diff_ind)==0)
            for zmt = 1:1:length(diff_ind)+1, % finding segments in time
                %for all the high rate data
                if(zmt == 1)                 % first segment
                    zm_st = 1;
                else
                    zm_st = diff_ind(zmt-1)+1;
                end

                if(zmt == length(diff_ind)+1) % last segment
                    zm_end = LM;
                else
                    zm_end = diff_ind(zmt);
                end

                abs_obstimez = abs_obstime(zm_st:zm_end);
                obstimez = obstime(zm_st:zm_end);
                fiqpowz = piqpowdata(zm_st:zm_end);
                fiqphz = piqphdata(zm_st:zm_end);
                    
                %process only if considerable amount of data present
                if length(obstimez)>1000,
                
                %length of each continuous segment in time
                LT = length(obstimez);
                % If there exists segments in phase data itself, handle
                % those first
                %find discontinuities in phase data
                rawphase = fiqphz;
                indphtemp = find( abs(diff(diff(rawphase))) > 5);

                %                     disp(zmt)
                %                     figure
                %                     plot(diff(diff(rawphase)))
                if isempty(indphtemp) ~= 1,
                    indph = indphtemp(2:2:end);
                    for zmph = 1:1:length(indph)+1,
                        %for all the high rate data
                        if(zmph == 1)                 % first segment
                            zm_sts = 1;
                        else
                            zm_sts = indph(zmph-1)+1;
                        end

                        if(zmph == length(indph)+1) % last segment
                            zm_ende = LT;
                        else
                            zm_ende = indph(zmph);
                        end
                        fignum = fignum +1;

%                         disp(zmph)
%                         disp(obstime(zm_st))
%                         disp(obstime(zm_end))

                        %Phase differenced zoomed in variables
                        obstimezphd = obstimez(zm_sts:zm_ende);
                        abs_obstimezphd = abs_obstimez(zm_sts:zm_ende);
                        fiqpowzphd = fiqpowz(zm_sts:zm_ende);
                        fiqphzphd = fiqphz(zm_sts:zm_ende);

                        %process only if considerable amount of data present
                        if length(obstimezphd)>1000,
                        
                        %Post processing goes here for the longer continuous
                        %lengths of data
                        % Polyfit subtraction
                        degree = 3;
                        rawphase = fiqphzphd;
                        %polynomial fitting and subtraction
                        poly_coef = polyfit(obstimezphd,...
                            rawphase,degree);
                        poly = polyval(poly_coef,obstimezphd);
                        poly_sub = rawphase - poly;
                        phaseph = detrend(poly_sub);

                        % %AJ filtering of phase
                        % %Filter phase in frequency domain
                        tlength = obstimezphd(end) - obstimezphd(1);
                        dfreq = 1/tlength;
                        % disp(tlength)
                        % disp(dfreq)
                        Lzph = length(phaseph);
                        % NFFT = 2^nextpow2(L);
                        % fft of phase data
                        wdata_fft = fftshift(fft(ifftshift(phaseph))); 
                        cutoff = 0.1; %[Hz]
                        order = 6; %order of the butterworth filter
                        fsamp = 50; %[Hz]
                        %create a low pass butterworth filter Kernel
                        butterlow = AJbutter(Lzph,cutoff,order,dfreq);
                        butterhi = 1.0 - butterlow;
                        wdata_filt = ...
                            fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                        phase_fftfilt = real(wdata_filt);
                        fiqphzphd = phase_fftfilt;

                        %AJ way for processing of power
                        Tdata = obstimezphd(end) - obstimezphd(1);
                        if Tdata <=60,
                            fc = 0.2;
                        else
                            fc = 0.1; %Desired cut off freq in Hz
                        end
                        %sampling frequency
                        fsamp = round(1/(obstimezphd(4)-obstimezphd(3)));
                        norder = 6; %order of the butterworth filter
                        power = fiqpowzphd';
                        low_pass_power = ...
                            butterworth_discrete(power,fsamp,fc,...
                            norder,'lp');
                        low_pass_power = ...
                            low_pass_power(1:length(power));
                        power1 = power./low_pass_power;
                        fiqpowzphd = power1';
                        
                        edgeT = 10; %edge time in seconds to be thrown out
                        edsamp = edgeT * fsamp;
                        obstimezphd = obstimezphd(edsamp+1:end-edsamp);
                        fiqphzphd = fiqphzphd(edsamp+1:end-edsamp);
                        fiqpowzphd = fiqpowzphd(edsamp+1:end-edsamp);
                        abs_obstimezphd = abs_obstimezphd(edsamp+1:end-edsamp);
                        
                        %Get S4 and sigma phi data, decide the filter length,
                        %ususally 60 seconds
                        Tdata = obstimezphd(end) - obstimezphd(1);
                        if Tdata >=60,
                            Tfilt = 60; %Desired filter length for S4
                            %and sigma_phi in [s]

                            % Time length of data in [s]
                            DataLength = length(fiqphzphd);
                            %Number of samples in data
                            FL = floor(Tfilt/Tdata*DataLength);
                            [S4outphd, SPHoutphd] = ...
                                slidingS4SigmaPHI(fiqpowzphd,fiqphzphd,FL);

                            %save into the Mega array
                            AOBSTIME = [AOBSTIME, abs_obstimezphd'];
                            S4DATA = [S4DATA, S4outphd'];
                            SPHDATA = [SPHDATA, SPHoutphd'];

                            %plot S4 and sigma_phi these zoomed in segments
                            figure('Visible','Off')
                            subplot(2,1,1)
                            set(gca,'FontSize',fontsz)
                            %To have both the powers overlapped
                            plot(obstimezphd,S4outphd,...
                                'b-','MarkerSize',marksz);  %if power-mean(power)
                            grid on
                            str = strcat('Scintillation index and standard',...
                                {' deviation of phase (zoomed in) for '},...
                                'the signal: ', signal,', PRN:',num2str(kk),...
                                {', fig'},num2str(fignum));
                            %     str = strcat('Detrended IQ power and HPF phase of the
                            %       signal PRN:',num2str(kk));
                            %     str = strcat('Detrended, polyfit subtracted IQ power
                            %     and ',... ' phase of the signal: ',signal,', PRN:',
                            %       num2str(kk));
                            title(str)
                            xstring = strcat({'Time [s] after '},num2str(hour),...
                                ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                                ,num2str(day),'/',num2str(year));
                            xlabel(xstring);
                            ylabel('S_4')
                           
                            subplot(2,1,2)
                            set(gca,'FontSize',fontsz)
                            plot(obstimezphd,SPHoutphd,'g-','MarkerSize',marksz);
                            grid on
                            xstring = strcat({'Time [s] after '},num2str(hour),...
                                ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                                ,num2str(day),'/',num2str(year));
                            xlabel(xstring);
                            ylabel('\sigma_{\phi}')

                            str_name = strcat(signal,'_S4_Sphi_PRN',num2str(kk),...
                                'zoom',num2str(fignum));
                            if(plot_type == 1)
                                plotfile = strcat(outdir,str_name,'.eps');
                                saveas(gcf,plotfile,'epsc2');
                            else
                                plotfile = strcat(outdir,str_name,'.png');
                                saveas(gca,plotfile);
                            end
                            clear gca
                            close
                        end % only consider greater than 60 seconds long data

                        end% process only if considerable data present
                        
                    end %for zoom: zmph of all phase segments

                else %if no discontinutites in phase exist:

                    fignum = fignum + 1;

                    %Post processing goes here for the longer continuous
                    %lengths of data
                    % Polyfit subtraction
                    degree = 3;
                    rawphase = fiqphz;
                    %polynomial fitting and subtraction
                    poly_coef = polyfit(obstimez,...
                        rawphase,degree);
                    poly = polyval(poly_coef,obstimez);
                    poly_sub = rawphase - poly;
                    phase = detrend(poly_sub);

                    % %AJ filtering of phase
                    % %Filter phase in frequency domain
                    tlength = obstimez(end) - obstimez(1);
                    dfreq = 1/tlength;
                    % disp(tlength)
                    % disp(dfreq)
                    Lz = length(phase);
                    % NFFT = 2^nextpow2(L);
                    wdata_fft = fftshift(fft(ifftshift(phase))); % fft of phase data
                    cutoff = 0.1; %[Hz]
                    order = 6; %order of the butterworth filter
                    fsamp = 50; %[Hz]
                    %create a low pass butterworth filter Kernel
                    butterlow = AJbutter(Lz,cutoff,order,dfreq);
                    butterhi = 1.0 - butterlow;
                    wdata_filt = fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                    phase_fftfilt = real(wdata_filt);
                    fiqphz = phase_fftfilt;
                    
                    %AJ way for processing of power
                    Tdata = obstimez(end) - obstimez(1);
                    if Tdata <=60,
                        fc = 0.2;
                    else
                        fc = 0.1; %Desired cut off freq in Hz
                    end
                    %sampling frequency
                    fsamp = round(1/(obstimez(4)-obstimez(3)));
                    norder = 6; %order of the butterworth filter
                    power = fiqpowz';
                    low_pass_power = ...
                        butterworth_discrete(power,fsamp,fc,...
                        norder,'lp');
                    low_pass_power = ...
                        low_pass_power(1:length(power));
                    power1 = power./low_pass_power;
                    fiqpowz = power1';
                    
                    edgeT = 10; %edge time in seconds to be thrown out
                    edsamp = edgeT * fsamp;
                    obstimez = obstimez(edsamp+1:end-edsamp);
                    fiqphz = fiqphz(edsamp+1:end-edsamp);
                    fiqpowz = fiqpowz(edsamp+1:end-edsamp);
                    abs_obstimez = abs_obstimez(edsamp+1:end-edsamp);
                    
                    %Get S4 and sigma phi data, decide the filter length,
                    %ususally 60 seconds
                    Tdata = obstimez(end) - obstimez(1);
                    if Tdata >=60,
                        Tfilt = 60; %Desired filter length for S4 
                        %and sigma_phi in [s]
                        
                        % Time length of data in [s]
                        DataLength = length(fiqphz);
                        %Number of samples in data
                        FL = floor(Tfilt/Tdata*DataLength);
                        [S4out, SPHout] = slidingS4SigmaPHI(fiqpowz,fiqphz,FL);
                        
                        %save into the Mega array
                        AOBSTIME = [AOBSTIME, abs_obstimez'];
                        S4DATA = [S4DATA, S4out'];
                        SPHDATA = [SPHDATA, SPHout'];                        
                        
                        %plot S4 and sigma_phi these zoomed in segments
                        figure('Visible','Off') 
                        subplot(2,1,1)
                        set(gca,'FontSize',fontsz)
                        %To have both the powers overlapped
                        plot(obstimez,S4out,...
                            'b-','MarkerSize',marksz);  %if power-mean(power)
                        grid on
                        str = strcat('Scintillation index and standard',...
                            {' deviation of phase (zoomed in) for '},...
                            'the signal: ', signal,', PRN:',num2str(kk),...
                            {', fig'},num2str(fignum));
                        %     str = strcat('Detrended IQ power and HPF phase of the
                        %       signal PRN:',num2str(kk));
                        %     str = strcat('Detrended, polyfit subtracted IQ power
                        %     and ',... ' phase of the signal: ',signal,', PRN:',
                        %       num2str(kk));
                        title(str)
                        xstring = strcat({'Time [s] after '},num2str(hour),...
                            ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                            ,num2str(day),'/',num2str(year));
                        xlabel(xstring);
                        ylabel('S_4')

                        subplot(2,1,2)
                        set(gca,'FontSize',fontsz)
                        plot(obstimez,SPHout,'g-','MarkerSize',marksz);
                        grid on
                        xstring = strcat({'Time [s] after '},num2str(hour),...
                            ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                            ,num2str(day),'/',num2str(year));
                        xlabel(xstring);
                        ylabel('\sigma_{\phi}')
                        str_name = strcat(signal,'_S4_Sphi_PRN',num2str(kk),...
                            'zoom',num2str(fignum));
                        if(plot_type == 1)
                            plotfile = strcat(outdir,str_name,'.eps');
                            saveas(gcf,plotfile,'epsc2');
                        else
                            plotfile = strcat(outdir,str_name,'.png');
                            saveas(gca,plotfile);
                        end
                        clear gca
                        close
                    end % only consider greater than 60 seconds long data               
                end % if-else-ends discontinuities exist in the phase data
                end% process only if length(obstimecz)>1000,
            end %for zoom: zmt of all segments

        else % no discontinuities in time
             

            % If no discontinuities exist, process here
            % Polyfit subtraction
            %disp('continuous segment of time data')

            %find discontinuities in phase data
            rawphase = piqphdata;
            indphtemp = find( abs(diff(diff(rawphase))) > 100);
            fignum = 0;
            if isempty(indphtemp) ~= 1,
                indph = indphtemp(2:2:end);
                for zmph = 1:1:length(indph)+1,
                    %for all the high rate data
                    %  obstimez = obstime((diff_ind(zm))*(zm-1)+1:...
                    %   zm*diff_ind(zm));
                    if(zmph == 1)                 % first segment
                        zm_st = 1;
                    else
                        zm_st = indph(zmph-1)+1;
                    end

                    if(zmph == length(indph)+1) % last segment
                        zm_end = LM;
                    else
                        zm_end = indph(zmph);
                    end

                    fignum = fignum +1;

%                     disp(zmph)
%                     disp(obstime(zm_st))
%                     disp(obstime(zm_end))

                    abs_obstimez = abs_obstime(zm_st:zm_end);
                    obstimez = obstime(zm_st:zm_end);
                    fiqpowz = piqpowdata(zm_st:zm_end);
                    fiqphz = piqphdata(zm_st:zm_end);

                    %process only if considerable amount of data present
                    if length(obstimez)>1000,
                    
                    %Post processing goes here for the longer continuous
                    %lengths of data
                    % Polyfit subtraction
                    degree = 3;
                    rawphase = fiqphz;
                    %polynomial fitting and subtraction
                    poly_coef = polyfit(obstimez,...
                        rawphase,degree);
                    poly = polyval(poly_coef,obstimez);
                    poly_sub = rawphase - poly;
                    phase = detrend(poly_sub);
                    
                    % %AJ filtering of phase
                    % %Filter phase in frequency domain
                    tlength = obstimez(end) - obstimez(1);
                    dfreq = 1/tlength;
                    % disp(tlength)
                    % disp(dfreq)
                    Lz = length(phase);
                    % fft of phase data
                    wdata_fft = fftshift(fft(ifftshift(phase))); 
                    cutoff = 0.1; %[Hz]
                    order = 6; %order of the butterworth filter
                    fsamp = 50; %[Hz]
                    %create a low pass butterworth filter Kernel
                    butterlow = AJbutter(Lz,cutoff,order,dfreq);
                    butterhi = 1.0 - butterlow;
                    wdata_filt = fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                    phase_fftfilt = real(wdata_filt);
                    fiqphz = phase_fftfilt;
                    
                    %AJ way for processing of power
                    Tdata = obstimez(end) - obstimez(1);
                    if Tdata <=60,
                        fc = 0.2;
                    else
                        fc = 0.1; %Desired cut off freq in Hz
                    end
                    %sampling frequency
                    fsamp = round(1/(obstimez(4)-obstimez(3)));
                    norder = 6; %order of the butterworth filter
                    power = fiqpowz';
                    low_pass_power = ...
                        butterworth_discrete(power,fsamp,fc,...
                        norder,'lp');
                    low_pass_power = ...
                        low_pass_power(1:length(power));
                    power1 = power./low_pass_power;
                    fiqpowz = power1';
                    
                    edgeT = 10; %edge time in seconds to be thrown out
                    edsamp = edgeT * fsamp;
                    obstimez = obstimez(edsamp+1:end-edsamp);
                    fiqphz = fiqphz(edsamp+1:end-edsamp);
                    fiqpowz = fiqpowz(edsamp+1:end-edsamp);
                    abs_obstimez = abs_obstimez(edsamp+1:end-edsamp);

                     %Get S4 and sigma phi data, decide the filter length,
                    %ususally 60 seconds
                    Tdata = obstimez(end) - obstimez(1);
                    if Tdata >=60,
                        Tfilt = 60; %Desired filter length for S4 
                        %and sigma_phi in [s]
                        
                        % Time length of data in [s]
                        DataLength = length(fiqphz);
                        %Number of samples in data
                        FL = floor(Tfilt/Tdata*DataLength);
                        [S4out, SPHout] = slidingS4SigmaPHI(fiqpowz,fiqphz,FL);
                        
                        %save into the Mega array
                        AOBSTIME = [AOBSTIME, abs_obstimez'];
                        S4DATA = [S4DATA, S4out'];
                        SPHDATA = [SPHDATA, SPHout'];                        
                        
                        %plot S4 and sigma_phi these zoomed in segments
                        figure('Visible','Off') 
                        subplot(2,1,1)
                        set(gca,'FontSize',fontsz)
                        %To have both the powers overlapped
                        plot(obstimez,S4out,...
                            'b-','MarkerSize',marksz);  %if power-mean(power)
                        grid on
                        str = strcat('Scintillation index and standard',...
                            {' deviation of phase (zoomed in) for '},...
                            'the signal: ', signal,', PRN:',num2str(kk),...
                            {', fig'},num2str(fignum));
                        %     str = strcat('Detrended IQ power and HPF phase of the
                        %       signal PRN:',num2str(kk));
                        %     str = strcat('Detrended, polyfit subtracted IQ power
                        %     and ',... ' phase of the signal: ',signal,', PRN:',
                        %       num2str(kk));
                        title(str)
                        xstring = strcat({'Time [s] after '},num2str(hour),...
                            ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                            ,num2str(day),'/',num2str(year));
                        xlabel(xstring);
                        ylabel('S_4')

                        subplot(2,1,2)
                        set(gca,'FontSize',fontsz)
                        plot(obstimez,SPHout,'g-','MarkerSize',marksz);
                        grid on
                        xstring = strcat({'Time [s] after '},num2str(hour),...
                            ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                            ,num2str(day),'/',num2str(year));
                        xlabel(xstring);
                        ylabel('\sigma_{\phi}')

                        str_name = strcat(signal,'_S4_Sphi_PRN',num2str(kk),...
                            'zoom',num2str(fignum));
                        if(plot_type == 1)
                            plotfile = strcat(outdir,str_name,'.eps');
                            saveas(gcf,plotfile,'epsc2');
                        else
                            plotfile = strcat(outdir,str_name,'.png');
                            saveas(gca,plotfile);
                        end
                        clear gca
                        close
                    end % only consider greater than 60 seconds long data
                    end %process only if length(obstimecz)>1000,
                end %for zoom: zmph of all segments
            else
                
                %process only if considerable amount of data present
                if length(obstime)>1000,
                
                % If there are no discontinuities in phase data

                degree = 3;
                rawphase = piqphdata;
                %polynomial fitting and subtraction
                poly_coef = polyfit(obstime,...
                    rawphase,degree);
                poly = polyval(poly_coef,obstime);
                poly_sub = rawphase - poly;
                phase = detrend(poly_sub);

                % %AJ filtering of phase
                % %Filter phase in frequency domain
                tlength = obstime(end) - obstime(1);
                dfreq = 1/tlength;
                Lz = length(phase);
                % NFFT = 2^nextpow2(L);
                wdata_fft = fftshift(fft(ifftshift(phase))); % fft of phase data
                cutoff = 0.1; %[Hz]
                order = 6; %order of the butterworth filter
                fsamp = 50; %[Hz]
                %create a low pass butterworth filter Kernel
                butterlow = AJbutter(Lz,cutoff,order,dfreq);
                butterhi = 1.0 - butterlow;
                wdata_filt = fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                phase_fftfilt = real(wdata_filt);
                piqphdata = phase_fftfilt;
                
                %AJ way for processing of power
                Tdata = obstime(end) - obstime(1);
                if Tdata <=60,
                    fc = 0.2;
                else
                    fc = 0.1; %Desired cut off freq in Hz
                end
                %sampling frequency
                fsamp = round(1/(obstime(4)-obstime(3)));
                norder = 6; %order of the butterworth filter
                power = piqpowdata';
                low_pass_power = ...
                    butterworth_discrete(power,fsamp,fc,...
                    norder,'lp');
                low_pass_power = ...
                    low_pass_power(1:length(power));
                power1 = power./low_pass_power;
                piqpowdata = power1';    
                
                edgeT = 10; %edge time in seconds to be thrown out
                edsamp = edgeT * fsamp;
                obstime = obstime(edsamp+1:end-edsamp);
                piqphdata = piqphdata(edsamp+1:end-edsamp);
                piqpowdata = piqpowdata(edsamp+1:end-edsamp);
                abs_obstime = abs_obstime(edsamp+1:end-edsamp);

                % plot S4, sigma_phi only one segment existing in data
                %Plot Carrier Phase
                %Get S4 and sigma phi data, decide the filter length,
                %ususally 60 seconds
                Tdata = obstime(end) - obstime(1);
                if Tdata >=60,
                    Tfilt = 60; %filter for S4 and sigma phi is 60 sec
                    
                    % Time length of data in [s]
                    DataLength = length(piqphdata);
                    %Number of samples in data
                    FL = floor(Tfilt/Tdata*DataLength);
                    [S4out, SPHout] = slidingS4SigmaPHI(piqpowdata,piqphdata,FL);
                    
                    %save into the Mega array
                    AOBSTIME = [AOBSTIME, abs_obstime'];
                    S4DATA = [S4DATA, S4out'];
                    SPHDATA = [SPHDATA, SPHout'];
                                        
                    figure('Visible','Off') 
                    subplot(2,1,1)
                    set(gca,'FontSize',fontsz)
                    %To have both the powers overlapped
                    plot(obstime,S4out,...
                        'b-','MarkerSize',marksz);  %if power-mean(power)
                    grid on
                    str = strcat('Scintillation index and standard',...
                        {' deviation of phase for '},...
                        'the signal: ', signal,', PRN:',num2str(kk));
                    %     str = strcat('Detrended IQ power and HPF phase of the
                    %       signal PRN:',num2str(kk));
                    %     str = strcat('Detrended, polyfit subtracted IQ power
                    %     and ',... ' phase of the signal: ',signal,', PRN:',
                    %       num2str(kk));
                    title(str)
                    xstring = strcat({'Time [s] after '},num2str(hour),...
                        ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                        ,num2str(day),'/',num2str(year));
                    xlabel(xstring);
                    ylabel('S_4')
                    subplot(2,1,2)
                    set(gca,'FontSize',fontsz)
                    plot(obstime,SPHout,'g-','MarkerSize',marksz);
                    grid on
                    xstring = strcat({'Time [s] after '},num2str(hour),...
                        ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                        ,num2str(day),'/',num2str(year));
                    xlabel(xstring);
                    ylabel('\sigma_{\phi}')
                    str_name = strcat(signal,'_S4_Sphi_PRN',num2str(kk));
                    if(plot_type == 1)
                        plotfile = strcat(outdir,str_name,'.eps');
                        saveas(gcf,plotfile,'epsc2');
                    else
                        plotfile = strcat(outdir,str_name,'.png');
                        saveas(gca,plotfile);
                    end
                    clear gca
                    close
                    
                end % consider only longer than 60 s data
                end %process only if length(obstimecz)>1000,
            end % if zoomed figures based on phase data
        end% if plotting zoomed figures based on time data
          
        if (length(AOBSTIME)>10)
            %Save time, s4 and sigmaphi in the mega PRN arrays
            AOBSTIME = AOBSTIME(2:end);
            S4DATA = S4DATA(2:end);
            SPHDATA = SPHDATA(2:end);

            L_obs = length(AOBSTIME);
            OBSTIME_MPRN(1:L_obs,kk) = AOBSTIME;
            S4A_MPRN(1:L_obs,kk) = S4DATA;
            SPHA_MPRN(1:L_obs,kk) = SPHDATA;
            kk
        end
            
    end %isempty svdata
    close all
end %for kk -- PRN

% only if there is data for that signal go ahead
if(isempty(DATA)~=1),
    %The GPS week number will be the same in a day
    ORTW_ST = ortw(1);
    
    ORTS_min = min(ORTS); %Find the minimum value of orts: starting time for UT
    
    %gps to utc conversion
    [UTC_time, leap_sec, day_of_year]=gps2utc([ORTW_ST ORTS_min], 0);
    %UTC_time = year month day hour minutes seconds
    UTC_init = UTC_time(1,:);
    YEAR = UTC_init(1);
    MONTH = UTC_init(2);
    DAY = UTC_init(3);
    HOUR = UTC_init(4);
    
    %to completely randomize the colors of stack plot
%     RandStream.setDefaultStream ...
%      (RandStream('mt19937ar','seed',sum(100*clock)));
    

    %find the maximum scintillation that day for nomalization
    s4max = max(max(S4A_MPRN));
    sphimax = max(max(SPHA_MPRN));
    
    NScale = 2; %the values be scaled by this number

    %Displacement vector value on y axis
    S4_01 = 0.1/s4max;
    sigmaPhrad1 = 1/sphimax;
    
    % % Plotting Stacked high rate S4
    count = 1;
    figure
%     figure('Visible','Off') ;
%     set(gca,'FontSize',fontsz)
    for kk = 1:1:32,
        
        obstime = OBSTIME_MPRN(:,kk);
        s4adata = S4A_MPRN(:,kk);
        
        obstime = obstime(obstime~=0);
        s4adata = s4adata(s4adata~=0);
        
        if (length(obstime)>4),
            time = (obstime-HOUR*3600)/3600;
            S4 = s4adata/s4max*NScale+count;
            dt = time(4)-time(3);
            ind_diff= find(diff(time)>0.1);
            disc_ind = [1,ind_diff',length(time)];
            Ttime = 0;
            SS4 = 0;
            for idf = 1:1:length(disc_ind)-2,
                dTime = time(ind_diff(idf)):dt:time(ind_diff(idf)+1);
                timei = NaN*dTime';
                S4i = ones(size(timei));

                Ttime = [Ttime;time(disc_ind(idf)+1:disc_ind(idf+1));timei];
                SS4 = [SS4;S4(disc_ind(idf)+1:disc_ind(idf+1));S4i];
            end

            Ttime = [Ttime;time(disc_ind(end-1)+1:disc_ind(end))];
            SS4 = [SS4;S4(disc_ind(end-1)+1:disc_ind(end))];

            Ttime = Ttime(Ttime~=0);
            SS4 = SS4(SS4~=0);

            plot(Ttime,SS4,...
                'Color',[rand(1) rand(1) rand(1)]);
            PRN_label = (num2str(kk));
            text(time(1), S4(1), PRN_label, 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right')
            hold on
            count = count +1;
        end
    end
    
    %Change the fontsize of PRN numbers
    figureHandle = gca;
    set(findall(figureHandle,'type','text'),'fontSize',12)
    
    str = strcat({'High Rate Scintillation Index, S4, for '},num2str(count-1),...
    {' PRNs on '},signal, {' signal'});
    title(str)
    xstring = strcat({'Time [hours] after '},num2str(HOUR),...
        ':00 UT on: (mm/dd/yy)',num2str(MONTH),'/'...
        ,num2str(DAY),'/',num2str(YEAR));
    xlabel(xstring);
    set(gca,'ytick',[])
    ylabel('PRNs')
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',[-0.1 round(count/2) 0])
    %axis([-inf inf 0 count])
    axis tight
    
    %Draw an arrow annotation based on the x and y axis values from the plot
    
    axPos = get(gca,'Position'); %# gca gets the handle to the current axes

    % get the limits, i.e. min and max of the axes.
    xMinMax = xlim;
    yMinMax = ylim;

    %Compute x y values
    xPlot = [xMinMax(2)+0.4 xMinMax(2)+0.4];
    yPlot = [round(yMinMax(2)/2) round(yMinMax(2)/2)+S4_01*NScale];
    
    % calculate the annotation x and y from the plot x and y.
    xAn = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
    yAn = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);

    
    xAn
    yAn
    % If the Y-axis scale-arrow size goes beyond 1 (scaling is too high),
    % reduce it by 1 and replot the stacked plot for s4 figure
    if yAn(2)>1,
        close
        NScale = NScale -2;
        count = 1;
        figure
        %     figure('Visible','Off') ;
        %     set(gca,'FontSize',fontsz)
        for kk = 1:1:32,
            
            obstime = OBSTIME_MPRN(:,kk);
            s4adata = S4A_MPRN(:,kk);
            
            obstime = obstime(obstime~=0);
            s4adata = s4adata(s4adata~=0);
            
            if (length(obstime)>4),
                time = (obstime-HOUR*3600)/3600;
                S4 = s4adata/s4max*NScale+count;
                dt = time(4)-time(3);
                ind_diff= find(diff(time)>0.1);
                disc_ind = [1,ind_diff',length(time)];
                Ttime = 0;
                SS4 = 0;
                for idf = 1:1:length(disc_ind)-2,
                    dTime = time(ind_diff(idf)):dt:time(ind_diff(idf)+1);
                    timei = NaN*dTime';
                    S4i = ones(size(timei));
                    
                    Ttime = [Ttime;time(disc_ind(idf)+1:disc_ind(idf+1));timei];
                    SS4 = [SS4;S4(disc_ind(idf)+1:disc_ind(idf+1));S4i];
                end
                
                Ttime = [Ttime;time(disc_ind(end-1)+1:disc_ind(end))];
                SS4 = [SS4;S4(disc_ind(end-1)+1:disc_ind(end))];
                
                Ttime = Ttime(Ttime~=0);
                SS4 = SS4(SS4~=0);
                
                plot(Ttime,SS4,...
                    'Color',[rand(1) rand(1) rand(1)]);
                PRN_label = (num2str(kk));
                text(time(1), S4(1), PRN_label, 'VerticalAlignment','bottom', ...
                    'HorizontalAlignment','right')
                hold on
                count = count +1;
            end
        end
        
        %Change the fontsize of PRN numbers
        figureHandle = gca;
        set(findall(figureHandle,'type','text'),'fontSize',12)
        
        str = strcat({'High Rate Scintillation Index, S4, for '},num2str(count-1),...
            {' PRNs on '},signal, {' signal'});
        title(str)
        xstring = strcat({'Time [hours] after '},num2str(HOUR),...
            ':00 UT on: (mm/dd/yy)',num2str(MONTH),'/'...
            ,num2str(DAY),'/',num2str(YEAR));
        xlabel(xstring);
        set(gca,'ytick',[])
        ylabel('PRNs')
        ylabh = get(gca,'YLabel');
        set(ylabh,'Position',[-0.1 round(count/2) 0])
        %axis([-inf inf 0 count])
        axis tight
        
        %Draw an arrow annotation based on the x and y axis values from the plot
        
        axPos = get(gca,'Position'); %# gca gets the handle to the current axes
        
        % get the limits, i.e. min and max of the axes.
        xMinMax = xlim;
        yMinMax = ylim;
        
        %Compute x y values
        xPlot = [xMinMax(2)+0.4 xMinMax(2)+0.4];
        yPlot = [round(yMinMax(2)/2) round(yMinMax(2)/2)+S4_01*NScale];
        
        % calculate the annotation x and y from the plot x and y.
        xAn = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
        yAn = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
        
    end
    
    xAn
    yAn
    
%     txtar = annotation('textarrow',xAn,yAn,'String','S4=0.1','FontSize',12);
    
    str_name = strcat(signal,'_HighRateS4_StackPlot',num2str(count-1),'PRNs');
    if(plot_type == 1)
        %set(gcf, 'PaperPositionMode', 'auto');
        plotfile = strcat(outdir,str_name,'.eps');
        saveas(gcf,plotfile,'epsc2');
    else
        plotfile = strcat(outdir,str_name,'.png');
        saveas(gca,plotfile);
    end
    close
    
    %% Plotting Part of the code: sigma_phi
    count = 1;
    figure
%     figure('Visible','Off') ;
%     set(gca,'FontSize',fontsz)
    for kk = 1:1:32,
        
        obstime = OBSTIME_MPRN(:,kk);
        sphiadata = SPHA_MPRN(:,kk);
        
        obstime = obstime(obstime~=0);
        sphiadata = sphiadata(sphiadata~=0);
        
        %     sphiindex = find(sphiadata>10);
        %     sphiadata(sphiindex) =6;
        
        if (length(obstime)>4),
            time = (obstime-HOUR*3600)/3600;
            sigma_PH = sphiadata/sphimax*NScale+count;
            
            dt = time(4)-time(3);
            ind_diff= find(diff(time)>0.1);
            disc_ind = [1,ind_diff',length(time)];
            Ttime = 0;
            Ssigma_PH = 0;
            for idf = 1:1:length(disc_ind)-2,
                dTime = time(ind_diff(idf)):dt:time(ind_diff(idf)+1);
                timei = NaN*dTime';
                sigma_PHi = ones(size(timei));

                Ttime = [Ttime;time(disc_ind(idf)+1:disc_ind(idf+1));timei];
                Ssigma_PH = [Ssigma_PH;sigma_PH(disc_ind(idf)+1:disc_ind(idf+1));sigma_PHi];
            end

          Ttime = [Ttime;time(disc_ind(end-1)+1:disc_ind(end))];
          Ssigma_PH = [Ssigma_PH;sigma_PH(disc_ind(end-1)+1:disc_ind(end))];  
        
          Ttime = Ttime(Ttime~=0);
          Ssigma_PH = Ssigma_PH(Ssigma_PH~=0);
            
            plot(Ttime,Ssigma_PH,...
                'Color',[rand(1) rand(1) rand(1)]);
            PRN_label = (num2str(kk));
            text(time(1), sigma_PH(1), PRN_label, 'VerticalAlignment','bottom', ...
                'HorizontalAlignment','right')
            hold on
            count = count +1;
        end
    end
    
    %Change the fontsize of PRN numbers
    figureHandle = gca;
    set(findall(figureHandle,'type','text'),'fontSize',12)
    
    str = strcat({'High Rate Standard Deviation of Phase, \sigma_{\phi},'},...
    {' for '},num2str(count-1),...
    {' PRNs on '},signal,{' signal'});
    title(str)%,'FontSize',15)
    xstring = strcat({'Time [hours] after '},num2str(HOUR),...
        ':00 UT on: (mm/dd/yy)',num2str(MONTH),'/'...
        ,num2str(DAY),'/',num2str(YEAR));
    xlabel(xstring);
    set(gca,'ytick',[])
    ylabel('PRNs')
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',[-0.1 round(count/2) 0])
    %axis([-inf inf 0 count+2])
    axis tight    
%     axes('FontSize',15);
    
    %Draw an arrow annotation based on the x and y axis values from the plot
    xPlot = [xMinMax(2)+0.4 xMinMax(2)+0.4];
    yPlot = [round(yMinMax(2)/2)-1 round(yMinMax(2)/2)-1+sigmaPhrad1*NScale];

    axPos = get(gca,'Position'); %# gca gets the handle to the current axes

    % get the limits, i.e. min and max of the axes.
    xMinMax = xlim;
    yMinMax = ylim;

    % calculate the annotation x and y from the plot x and y.

    xAn = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
    yAn = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);

    % If the Y-axis scale-arrow size goes beyond 1 (scaling is too high),
    % reduce it by 1 and replot the stacked plot for sigma_phi figure
    if yAn(2)>1,
        close
        NScale = NScale -1;
        count = 1;
        figure
        %     figure('Visible','Off') ;
        %     set(gca,'FontSize',fontsz)
        for kk = 1:1:32,
            
            obstime = OBSTIME_MPRN(:,kk);
            sphiadata = SPHA_MPRN(:,kk);
            
            obstime = obstime(obstime~=0);
            sphiadata = sphiadata(sphiadata~=0);
            
            %     sphiindex = find(sphiadata>10);
            %     sphiadata(sphiindex) =6;
            
            if (length(obstime)>4),
                time = (obstime-HOUR*3600)/3600;
                sigma_PH = sphiadata/sphimax*NScale+count;
                
                dt = time(4)-time(3);
                ind_diff= find(diff(time)>0.1);
                disc_ind = [1,ind_diff',length(time)];
                Ttime = 0;
                Ssigma_PH = 0;
                for idf = 1:1:length(disc_ind)-2,
                    dTime = time(ind_diff(idf)):dt:time(ind_diff(idf)+1);
                    timei = NaN*dTime';
                    sigma_PHi = ones(size(timei));
                    
                    Ttime = [Ttime;time(disc_ind(idf)+1:disc_ind(idf+1));timei];
                    Ssigma_PH = [Ssigma_PH;sigma_PH(disc_ind(idf)+1:disc_ind(idf+1));sigma_PHi];
                end
                
                Ttime = [Ttime;time(disc_ind(end-1)+1:disc_ind(end))];
                Ssigma_PH = [Ssigma_PH;sigma_PH(disc_ind(end-1)+1:disc_ind(end))];
                
                Ttime = Ttime(Ttime~=0);
                Ssigma_PH = Ssigma_PH(Ssigma_PH~=0);
                
                plot(Ttime,Ssigma_PH,...
                    'Color',[rand(1) rand(1) rand(1)]);
                PRN_label = (num2str(kk));
                text(time(1), sigma_PH(1), PRN_label, 'VerticalAlignment','bottom', ...
                    'HorizontalAlignment','right')
                hold on
                count = count +1;
            end
        end
        
        %Change the fontsize of PRN numbers
        figureHandle = gca;
        set(findall(figureHandle,'type','text'),'fontSize',12)
        
        str = strcat({'High Rate Standard Deviation of Phase, \sigma_{\phi},'},...
            {' for '},num2str(count-1),...
            {' PRNs on '},signal,{' signal'});
        title(str)%,'FontSize',15)
        xstring = strcat({'Time [hours] after '},num2str(HOUR),...
            ':00 UT on: (mm/dd/yy)',num2str(MONTH),'/'...
            ,num2str(DAY),'/',num2str(YEAR));
        xlabel(xstring);
        set(gca,'ytick',[])
        ylabel('PRNs')
        ylabh = get(gca,'YLabel');
        set(ylabh,'Position',[-0.1 round(count/2) 0])
        %axis([-inf inf 0 count+2])
        axis tight
        %     axes('FontSize',15);
        
        %Draw an arrow annotation based on the x and y axis values from the plot
        xPlot = [xMinMax(2)+0.4 xMinMax(2)+0.4];
        yPlot = [round(yMinMax(2)/2)-1 round(yMinMax(2)/2)-1+sigmaPhrad1*NScale];
        
        axPos = get(gca,'Position'); %# gca gets the handle to the current axes
        
        % get the limits, i.e. min and max of the axes.
        xMinMax = xlim;
        yMinMax = ylim;
        
        % calculate the annotation x and y from the plot x and y.

        xAn = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
        yAn = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);

    end
    
%     txtar = annotation('textarrow',xAn,yAn,'String',...
%         '\sigma_{\phi}=1 rad','FontSize',12);

    str_name = strcat(signal,'_HighRateSPHI_StackPlot',num2str(count-1),'PRNs');
    if(plot_type == 1)
       % set(gcf, 'PaperPositionMode', 'auto');
        plotfile = strcat(outdir,str_name,'.eps');
        saveas(gcf,plotfile,'epsc2');
    else
        plotfile = strcat(outdir,str_name,'.png');
        saveas(gca,plotfile);
    end
    close
    
end % if DATA itself is not empty

%end of B set of plots
end % set_plot == B



%% Set Plot C: Plot common times
if strcmp(set_plot,'C') == 1,
% Create necessary directories
command = strcat('mkdir',{' '},folder_name,'plots_',signal,sep,...
    'HighRate_Common_times',sep);
system(cell2mat(command));
outdir = strcat(folder_name,'plots_',signal,sep,'HighRate_Common_times',sep);

%Create a common vector for all PRNs (Mega PRN)
orts_MPRN = zeros(floor(length(DATA)/2),32);
ortw_MPRN = zeros(floor(length(DATA)/2),32);
obstime_MPRN = zeros(floor(length(DATA)/2),32);
ph_MPRN = zeros(floor(length(DATA)/2),32);
pow_MPRN = zeros(floor(length(DATA)/2),32);
% length(DATA)

for kk = 1:1:32,
    sv = (DATA(:,8)==kk);
    svdata = DATA(sv,:);
    
    if(isempty(svdata)~=1),
                
        ort_iq = svdata(:,3); %ort seconds + fractional seconds
        % %Find and remove repeated entries
        %1st derivative to find the time slots (wherever there is a huge
        %difference, that's the BUG. There are repeated entries inside iq.log
        ort_iqd = diff(ort_iq);
        iqdn_ind=find(ort_iqd<0);
        for ii = 1:1:length(iqdn_ind),
            indii = find(ort_iq == ort_iq(iqdn_ind(ii)));
            ort_iq(indii(1):indii(end)-1)=0;
            %This will work well only if there are only 2 repetitions max. ort_sc =
            %ort_sc(ort_sc~=0) will remove the repeated values. But better leave
            %the zeros as they are and recreate the original "scdata" array based
            %on those zeros.
        end
        %Non-zero time stamps means without any repetition.
        nzero_ort = find(ort_iq~=0);
        svdata = svdata(nzero_ort,:);
        %Reread time for scint file(with no repeated entries).
        orts =svdata(:,3);
        
        L = length(svdata);
        
        %read all IQ data
        obstime_o = svdata(:,1); %original obstime
        ortw = svdata(:,2); % ORT weeks
%         orts = svdata(:,3); % ORT seconds
        pcdata = svdata(:,4); % processed carrier phase values
        piqphdata = svdata(:,5); % processed iq phase values
        piqpowdata = svdata(:,6); % processed iq power values
        refprn = svdata(:,7); %PRN of reference channel
        scintprn = svdata(:,8); % PRN for scintillating channel
        
        %gps to utc conversion
        [UTC_time, leap_sec, day_of_year]=gps2utc([ortw orts], 0);
        %UTC_time = year month day hour minutes seconds
        UTC_init = UTC_time(1,:);
        year = UTC_init(1);
        month = UTC_init(2);
        day = UTC_init(3);
        hour = UTC_init(4);
        
        %Observation time seconds after UTC_init hours:
        obstime = (UTC_time(:,4)-hour)*3600 + UTC_time(:,5)*60 ...
            + UTC_time(:,6);
        
        %absolute UTC observation time vector in seconds
        obstime_abs = (UTC_time(:,4))*3600 + UTC_time(:,5)*60 ...
            + UTC_time(:,6);
        
%         size(obstime_abs)
        ortw_MPRN(1:length(obstime_abs),kk) = ortw;
        orts_MPRN(1:length(obstime_abs),kk) = orts;
        obstime_MPRN(1:length(obstime_abs),kk) = obstime_abs;
        ph_MPRN(1:length(obstime_abs),kk) = piqphdata;
        pow_MPRN(1:length(obstime_abs),kk) = piqpowdata;
        
    end %isempty svdata
end %for kk -- PRN 

%Compare each PRN with another and note down common times
for kk = 1:1:32,
    %ortskk = orts_MPRN(:,kk);
    
    %having trouble with 0.001-0.009 difference in the times, thus intersec
    %t won't work. How about using something like: round(a*10^2)./10^2 to 
    %round the number to just close by number?
    
    obskk = round(obstime_MPRN(:,kk)*10^2)./10^2;
    for jj = kk+1:1:32,
        %ortsjj = orts_MPRN(:,jj);
        obsjj = round(obstime_MPRN(:,jj)*10^2)./10^2;
               
        [com_obstime, ijj, ikk] = intersect(obsjj,obskk);
        
        com_obstime=com_obstime(com_obstime~=0);
               
        phkk =  ph_MPRN(ikk,kk);
        phjj =  ph_MPRN(ijj,jj);
        
        phjj = phjj(phjj~=0);
        phkk = phkk(phkk~=0);
        
%         figure
%         plot(com_obstime,phjj)
%         hold on
%         pause
%         plot(com_obstime,phkk,'r')
        
        if(length(com_obstime)>10 && ~isempty(phjj)...
                && ~isempty(phkk))
            % If empty/ zero valued vectors, GPS2UTC fails giving
            % "undefined day_of_year variable" inside the gps2utc routine
            % itself.
            %gps to utc conversion
            ortsc = orts_MPRN(ikk,kk);
            ortwc = ortw_MPRN(ikk,kk);
            ortsc = ortsc(ortsc~=0);
            ortwc = ortwc(ortwc~=0);
            [UTC_time]=gps2utc([ortwc ortsc], 0);
            %UTC_time = year month day hour minutes seconds
            UTC_init = UTC_time(1,:);
            year = UTC_init(1);
            month = UTC_init(2);
            day = UTC_init(3);
            hour = UTC_init(4);
            
            %Observation time seconds after UTC_init hours:
            obstimec = (UTC_time(:,4)-hour)*3600 + UTC_time(:,5)*60 ...
                + UTC_time(:,6);
            
            L_obsc = length(obstimec);
            
            powkk =  pow_MPRN(ikk,kk);
            powjj =  pow_MPRN(ijj,jj);
            
            powkk =  powkk(powkk~=0);
            powjj =  powjj(powjj~=0);
            
            % Zoom certain patches
            %find if there indeed are data chunks with such a time gap
            difference = 5; % Difference in two data 'patches'
            %in seconds in the prn file
            diff_ind = find(diff(obstimec) > difference);
            
            if(isempty(diff_ind)==0)
                for zm = 1:1:length(diff_ind)+1,
                    %for all the high rate data
                    
                    if(zm == 1)                 % first segment
                        zm_st = 1;
                    else
                        zm_st = diff_ind(zm-1)+1;
                    end
                    
                    if(zm == length(diff_ind)+1) % last segment
                        zm_end = L_obsc;
                    else
                        zm_end = diff_ind(zm);
                    end
                    
                    obstimecz = obstimec(zm_st:zm_end);
                    
                    powkkz = powkk(zm_st:zm_end);
                    powjjz = powjj(zm_st:zm_end);
                    
                    phkkz = phkk(zm_st:zm_end);
                    phjjz = phjj(zm_st:zm_end);
                    
                    %process only if considerable amount of data present
                    if length(obstimecz)>1000,
                        %                     plot(powkkz)
                        %                     pause
                        %                     hold on
                        %                     plot(powjjz,'r')
                        %                     pause
                        
                        %Post processing goes here for the longer continuous
                        %lengths of data
                        % Polyfit subtraction
                        degree = 3;
                        rawphase = phkkz;
                        %polynomial fitting and subtraction
                        poly_coef = polyfit(obstimecz,...
                            rawphase,degree);
                        poly = polyval(poly_coef,obstimecz);
                        poly_sub = rawphase - poly;
                        phase = detrend(poly_sub);
                        
                        % %AJ filtering of phase
                        % %Filter phase in frequency domain
                        tlength = obstimecz(end) - obstimecz(1);
                        dfreq = 1/tlength;
                        % disp(tlength)
                        % disp(dfreq)
                        
                        Lz = length(phase);
                        % NFFT = 2^nextpow2(L);
                        wdata_fft = fftshift(fft(ifftshift(phase))); 
                        % fft of phase data
                        cutoff = 0.1; %[Hz]
                        order = 6; %order of the butterworth filter
                        %create a low pass butterworth filter Kernel
                        butterlow = AJbutter(Lz,cutoff,order,dfreq);
                        butterhi = 1.0 - butterlow;
                        wdata_filt = ...
                            fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                        phase_fftfilt = real(wdata_filt);
                        phkkz = phase_fftfilt;
                        
                        %AJ way for processing of power
                        Tdata = obstimecz(end) - obstimecz(1);
                        if Tdata <=60,
                            fc = 0.2;
                        else
                            fc = 0.1; %Desired cut off freq in Hz
                        end
                        %sampling frequency
                        fsamp = round(1/(obstimecz(4)-obstimecz(3)));
                        norder = 6; %order of the butterworth filter
                        power = powkkz';
                        low_pass_power = ...
                            butterworth_discrete(power,fsamp,fc,...
                            norder,'lp');
                        low_pass_power = ...
                            low_pass_power(1:length(power));
                        power1 = power./low_pass_power;
                        powkkz = power1';
                        
                        
                        %Post processing goes here for the longer continuous
                        %lengths of data
                        % Polyfit subtraction
                        degree = 3;
                        rawphase = phjjz;
                        %polynomial fitting and subtraction
                        poly_coef = polyfit(obstimecz,...
                            rawphase,degree);
                        poly = polyval(poly_coef,obstimecz);
                        poly_sub = rawphase - poly;
                        phase = detrend(poly_sub);
                        
                        % %AJ filtering of phase
                        % %Filter phase in frequency domain
                        tlength = obstimecz(end) - obstimecz(1);
                        dfreq = 1/tlength;
                        % disp(tlength)
                        % disp(dfreq)
                        
                        Lz = length(phase);
                        % NFFT = 2^nextpow2(L);
                        wdata_fft = fftshift(fft(ifftshift(phase))); 
                        % fft of phase data
                        cutoff = 0.1; %[Hz]
                        order = 6; %order of the butterworth filter
                        %create a low pass butterworth filter Kernel
                        butterlow = AJbutter(Lz,cutoff,order,dfreq);
                        butterhi = 1.0 - butterlow;
                        wdata_filt = ...
                            fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                        phase_fftfilt = real(wdata_filt);
                        phjjz = phase_fftfilt;
                        
                        
                        %AJ way for processing of power
                        Tdata = obstimecz(end) - obstimecz(1);
                        if Tdata <=60,
                            fc = 0.2;
                        else
                            fc = 0.1; %Desired cut off freq in Hz
                        end
                        %sampling frequency
                        fsamp = round(1/(obstimecz(4)-obstimecz(3)));
                        norder = 6; %order of the butterworth filter
                        power = powjjz';
                        low_pass_power = ...
                            butterworth_discrete(power,fsamp,fc,...
                            norder,'lp');
                        low_pass_power = ...
                            low_pass_power(1:length(power));
                        power1 = power./low_pass_power;
                        powjjz = power1';
                                
                        phjjztemp = phjjz(phjjz~=0);
                        phkkztemp = phkkz(phkkz~=0);
                        
                        edgeT = 10; %edge time in seconds to be thrown out
                        edsamp = edgeT * fsamp;
                        obstimecz = obstimecz(edsamp+1:end-edsamp);
                        phjjz = phjjz(edsamp+1:end-edsamp);
                        phkkz = phkkz(edsamp+1:end-edsamp);
                        powjjz = powjjz(edsamp+1:end-edsamp);
                        powkkz = powkkz(edsamp+1:end-edsamp);                        
                        
                        % plot only non-zero segments and non-negative time
                        % values!
                        if (min(obstimecz) > 0 && length(obstimecz)>10 && ...
                                ~isempty(phjjztemp) && ~isempty(phkkztemp))
                            %plot common power
                            figure('Visible','Off')
                            subplot(2,1,1)
                            set(gca,'FontSize',fontsz)
                            %To have both the powers overlapped
                            plot(obstimecz,10*log10(powkkz),'Color',...
                                peacock_blue,'MarkerSize',marksz);  
                            %if power-mean(power)
                            grid on
                            axis([-inf inf -5 5 ]);
                            str = strcat('Detrended, clock error corrected',...
                                {' (zoomed in) IQ power [dB] of '},...
                                'the signal: ', signal,', PRNs:',num2str(kk),...
                                '&',num2str(jj),{', fig'},num2str(zm));
                            title(str)
                            ystring = strcat('PRN',num2str(kk));
                            ylabel(ystring)
                            subplot(2,1,2)
                            set(gca,'FontSize',fontsz)
                            plot(obstimecz,10*log10(powjjz),'b',...
                                'MarkerSize',marksz);
                            grid on
                            axis([-inf inf -5 5 ]);
                            xstring = strcat({'Time [s] after '},...
                                num2str(hour),':00 UT on: (mm/dd/yy)',...
                                num2str(month),'/',num2str(day),...
                                '/',num2str(year));
                            xlabel(xstring);
                            ystring = strcat('PRN',num2str(jj));
                            ylabel(ystring)
                            str_name = strcat(signal,'_CorrIQPower_PRNs',...
                                num2str(kk),'AND',num2str(jj),'zoom',...
                                num2str(zm));
                            if(plot_type == 1)
                                plotfile = strcat(outdir,str_name,'.eps');
                                saveas(gcf,plotfile,'epsc2');
                            else
                                plotfile = strcat(outdir,str_name,'.png');
                                saveas(gca,plotfile);
                            end
                            clear gca
                            close
                            
                            figure('Visible','Off')
                            subplot(2,1,1)
                            set(gca,'FontSize',fontsz)
                            %To have both the phases overlapped
                            plot(obstimecz,phkkz,...
                                'g-','MarkerSize',marksz);
                            grid on
                            axis([-inf inf -2 2]);
                            str = strcat('Detrended, clock error corrected',...
                                {' (zoomed in) IQ phase [rad] of '},...
                                'the signal: ', signal,', PRNs:',num2str(kk),...
                                '&',num2str(jj),{', fig'},num2str(zm));
                            title(str)
                            ystring = strcat('PRN',num2str(kk));
                            ylabel(ystring)
                            subplot(2,1,2)
                            set(gca,'FontSize',fontsz)
                            plot(obstimecz,phjjz,'Color',...
                                forest_green,'MarkerSize',marksz);
                            grid on
                            axis([-inf inf -2 2]);
                            xstring = strcat({'Time [s] after '},...
                                num2str(hour),':00 UT on: (mm/dd/yy)',...
                                num2str(month),'/',num2str(day),...
                                '/',num2str(year));
                            xlabel(xstring);
                            ystring = strcat('PRN',num2str(jj));
                            ylabel(ystring)
                            str_name = strcat(signal,'_CorrIQPhase_PRNs',...
                                num2str(kk),'AND',num2str(jj),...
                                'zoom',num2str(zm));
                            if(plot_type == 1)
                                plotfile = strcat(outdir,str_name,'.eps');
                                saveas(gcf,plotfile,'epsc2');
                            else
                                plotfile = strcat(outdir,str_name,'.png');
                                saveas(gca,plotfile);
                            end
                            clear gca
                            close
                            
                        end % only non-zero zoomed in segments
                    end %process only considerable amounts of data
                end  %for all zoomed in parts
            else 
                %process only if considerable amount of data present
                if length(obstimec)>1000,
                %Process phase
                %Post processing goes here for the longer continuous
                %lengths of data
                % Polyfit subtraction
                degree = 3;
                rawphase = phkk;
                %polynomial fitting and subtraction
                poly_coef = polyfit(obstimec,...
                    rawphase,degree);
                poly = polyval(poly_coef,obstimec);
                poly_sub = rawphase - poly;
                phase = detrend(poly_sub);
                
                % %AJ filtering of phase
                % %Filter phase in frequency domain
                tlength = obstimec(end) - obstimec(1);
                dfreq = 1/tlength;
                % disp(tlength)
                % disp(dfreq)
                
                Lz = length(phase);
                % NFFT = 2^nextpow2(L);
                wdata_fft = fftshift(fft(ifftshift(phase))); % fft of phase data
                cutoff = 0.1; %[Hz]
                order = 6; %order of the butterworth filter
                %create a low pass butterworth filter Kernel
                butterlow = AJbutter(Lz,cutoff,order,dfreq);
                butterhi = 1.0 - butterlow;
                wdata_filt = fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                phase_fftfilt = real(wdata_filt);
                phkk = phase_fftfilt;
                
                %AJ way for processing of power
                Tdata = obstimec(end) - obstimec(1);
                if Tdata <=60,
                    fc = 0.2;
                else
                    fc = 0.1; %Desired cut off freq in Hz
                end
                %sampling frequency
                fsamp = round(1/(obstimec(4)-obstimec(3)));
                norder = 6; %order of the butterworth filter
                power = powkk';
                low_pass_power = ...
                    butterworth_discrete(power,fsamp,fc,...
                    norder,'lp');
                low_pass_power = ...
                    low_pass_power(1:length(power));
                power1 = power./low_pass_power;
                powkk = power1';
                
                %Post processing goes here for the longer continuous
                %lengths of data
                % Polyfit subtraction
                degree = 3;
                rawphase = phjj;
                %polynomial fitting and subtraction
                poly_coef = polyfit(obstimec,...
                    rawphase,degree);
                poly = polyval(poly_coef,obstimec);
                poly_sub = rawphase - poly;
                phase = detrend(poly_sub);
                
                % %AJ filtering of phase
                % %Filter phase in frequency domain
                tlength = obstimec(end) - obstimec(1);
                dfreq = 1/tlength;
                % disp(tlength)
                % disp(dfreq)
                
                Lz = length(phase);
                % NFFT = 2^nextpow2(L);
                wdata_fft = fftshift(fft(ifftshift(phase))); % fft of phase data
                cutoff = 0.1; %[Hz]
                order = 6; %order of the butterworth filter
                %create a low pass butterworth filter Kernel
                butterlow = AJbutter(Lz,cutoff,order,dfreq);
                butterhi = 1.0 - butterlow;
                wdata_filt = fftshift(ifft(ifftshift(wdata_fft.*butterhi')));
                phase_fftfilt = real(wdata_filt);
                phjj = phase_fftfilt;
                
                %AJ way for processing of power
                Tdata = obstimec(end) - obstimec(1);
                if Tdata <=60,
                    fc = 0.2;
                else
                    fc = 0.1; %Desired cut off freq in Hz
                end
                %sampling frequency
                fsamp = round(1/(obstimec(4)-obstimec(3)));
                norder = 6; %order of the butterworth filter
                power = powjj';
                low_pass_power = ...
                    butterworth_discrete(power,fsamp,fc,...
                    norder,'lp');
                low_pass_power = ...
                    low_pass_power(1:length(power));
                power1 = power./low_pass_power;
                powjj = power1';
                
                edgeT = 10; %edge time in seconds to be thrown out
                edsamp = edgeT * fsamp;
                Lobs = length(obstimec);
                obstimec = obstimec(edsamp+1:Lobs-edsamp);
                phjj = phjj(edsamp+1:Lobs-edsamp);
                phkk = phkk(edsamp+1:Lobs-edsamp);
                powjj = powjj(edsamp+1:Lobs-edsamp);
                powkk = powkk(edsamp+1:Lobs-edsamp);
                
                %plot common power
                figure('Visible','Off')
                subplot(2,1,1)
                set(gca,'FontSize',fontsz)
                %To have both the powers overlapped
                plot(obstimec,10*log10(powkk),'Color',...
                    peacock_blue,'MarkerSize',marksz);  %if power-mean(power)
                grid on
                axis([-inf inf -5 5 ]);
                str = strcat('Detrended, clock error corrected',...
                    {' IQ power [dB] of '},...
                    'the signal: ', signal,', PRNs:',num2str(kk),...
                    '&',num2str(jj));
                title(str)
                ystring = strcat('PRN',num2str(kk));
                ylabel(ystring)
                subplot(2,1,2)
                set(gca,'FontSize',fontsz)
                plot(obstimec,10*log10(powjj),'b',...
                    'MarkerSize',marksz);
                grid on
                axis([-inf inf -5 5 ]);
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ystring = strcat('PRN',num2str(jj));
                ylabel(ystring)
                str_name = strcat(signal,'_CorrIQPower_PRNs',...
                    num2str(kk),'AND',num2str(jj));
                if(plot_type == 1)
                    plotfile = strcat(outdir,str_name,'.eps');
                    saveas(gcf,plotfile,'epsc2');
                else
                    plotfile = strcat(outdir,str_name,'.png');
                    saveas(gca,plotfile);
                end
                clear gca
                close
                
                
                figure('Visible','Off') 
                subplot(2,1,1)
                set(gca,'FontSize',fontsz)
                %To have both the powers overlapped
                plot(obstimec,phkk,...
                    'g-','MarkerSize',marksz);  %if power-mean(power)
                grid on
                axis([-inf inf -2 2]);
                str = strcat('Detrended, clock error corrected',...
                    {' IQ phase [rad] of '},...
                    'the signal: ', signal,', PRNs:',num2str(kk),...
                    '&',num2str(jj));
                title(str)
                ystring = strcat('PRN',num2str(kk));
                ylabel(ystring)
                subplot(2,1,2)
                set(gca,'FontSize',fontsz)
                plot(obstimec,phjj,'Color',...
                    forest_green,'MarkerSize',marksz);
                grid on
                axis([-inf inf -2 2]);
                xstring = strcat({'Time [s] after '},num2str(hour),...
                    ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
                    ,num2str(day),'/',num2str(year));
                xlabel(xstring);
                ystring = strcat('PRN',num2str(jj));
                ylabel(ystring)
                str_name = strcat(signal,'_CorrIQPhase_PRNs',...
                    num2str(kk),'AND',num2str(jj));
                if(plot_type == 1)
                    plotfile = strcat(outdir,str_name,'.eps');
                    saveas(gcf,plotfile,'epsc2');
                else
                    plotfile = strcat(outdir,str_name,'.png');
                    saveas(gca,plotfile);
                end
                clear gca
                close
                end %process data only if of considerable length
            end % no segments to zoom on
        end %there exists common data
    end % for jj -- 2nd PRN
end %for kk -- 1st PRN
% End Set C
end % set_plot == C
toc
