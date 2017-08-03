function varargout = plotPFISRvs(ts,te,mflag,vflag)

dbstop if error;
[prop, op_path] = ver_chk;

if ~nargin
    ts = [2013 12 8 7 0 0];
    te = [2013 12 8 8 0 0];
    vflag = 'vmag_vang';
end
switch vflag
    case 'vmag_vang'
        %plot vgmag and vgang
        col = [20 21];
        errcol = [22 23];
    case 've_vn'
        %plot vge and vgn [14 15]
        col = [14 15];
        errcol = [17 18];
    case 'debug'
        col = [20 21];
        errcol = [22 23];
    otherwise
        megadata = {};
        lat = NaN;
        dtau = NaN;
        ebar = [];
        varargout = {megadata, lat, dtau, ebar};
        return;
end

load('PFISR_data1.txt');

PFISR_data = PFISR_data1;
data = [datenum(PFISR_data(:,1:6)) PFISR_data(:,7:end)];

data(:,11) = 90 - data(:,11) ;

data(data(:,11)<0,11) = data(data(:,11)<0,11) + 360;

% max(data,13)
data = data(data(:,1)<=datenum(te)+10/24/3600&data(:,1)>=datenum(ts)-10/24/3600,:);
if ~isempty(data)
    tdt = unique(data(:,1),'stable');
    dtau = unique(data(:,2),'rows');
    lat = unique(data(:,3),'rows');

    %rotation from geographic drift vg:[vge,vgn,vgz] to magnetic field drift vmf:[vipe,vipn,vi6]
    %actually a two-step rotation: 
    %geographic -> geomagnetic -> magnetic field line
    %declination angle delta, taken from Heinselman and Nicolls, 2008
    de = 19.0519*pi/180;
    %inclination (dip) angle I, taken from Heinselman and Nicolls, 2008
    in = 77.5175*pi/180;
    %rotation matrix
    R = [cos(de) -sin(de) 0;
        sin(in)*sin(de) cos(de)*sin(in) cos(in)
        -cos(in)*sin(de) -cos(in)*cos(de) sin(in)];

    %inverse transform to geographic velocities
    vg = R\data(:,4:6)';
    vge = vg(1,:);vgn = vg(2,:);
    
    %dvmmag and dvmang
    dvm = data(:,7:9);
    dvme = data(:,7);
    dvmn = data(:,8);
%     dvg = R\data(:,7:9)';
    for vi = 1:size(dvm,1)
        %assume the three velocity components are uncorrelated
        sigmavm = diag(dvm(vi,:).^2);
        sigmavg = inv(R)*sigmavm*inv(R)';
        dvge(vi,:) = sqrt(sigmavg(1,1));
        dvgn(vi,:) = sqrt(sigmavg(2,2));
        dvgu(vi,:) = sqrt(sigmavg(3,3));
        sigmavgen(vi,:) = sigmavg(1,2);
    end
        
    %compute the jacobian
    v = sqrt(vge.^2+vgn.^2);
    vgmag = v;
    th = 180/pi*atan2(vgn,vge);
    th(th<0) = th(th<0) + 360;
    vgang = th;
    dvde = vge./v;
    dvdn = vgn./v;
    dthde = -vgn./v.^2;
    dthdn = vge./v.^2;
    
%     Jv = [dvde dvdn];
%     Jth = [dthde dthdn];
    sigmavgmag = dvde'.^2.*dvge.^2 + ...
        dvdn'.^2.*dvgn.^2 + ...
        2*dvde'.*dvdn'.*sigmavgen;
    dvgmag = sqrt(sigmavgmag);
    sigmavgang = dthde'.^2.*dvge.^2 + ...
        dthdn'.^2.*dvgn.^2 + ...
        2*dthde'.*dthdn'.*sigmavgen;
    dvgang = sqrt(sigmavgang)*180/pi;
    
    dvgmagJ = dvgmag;
    dvgangJ = dvgang;

%     subplot(2,1,1);plot(data(:,3),dvgmag,data(:,3),dvgmagE,data(:,3),dvgmagJ)
%     subplot(2,1,2);plot(data(:,3),dvgang,data(:,3),dvgangE,data(:,3),dvgangJ)
%     legend('Easiest','Ensemble Average','Jacobian')
    vm = data(:,4:6)';
        
    data = [data vg' dvge dvgn dvgu vgmag' vgang' dvgmagJ dvgangJ];
    %geographic velocities vge vgn vgz dvge dvgn dvgu vgmag vgang dvgmag dvgang 
    %[14 15 16, 17 18 19, 20 21, 22 23]
    switch mflag
        case 'lat'
            ylabelstr = 'Latitude [\circ]';   
            ytick = 1:length(lat);
            yticklabel = lat;
%             dtau = dtau(1);        
%             titlestr = (['\Delta\tau = ',num2str(dtau),' s']);
        case 'dtau'
            ylabelstr = '\Delta\tau [s]';   
            ytick = 1:length(dtau);
            yticklabel = dtau;
%             lat = 66.25;
%             titlestr = (['Latitude = ',num2str(lat),'\circ']);
    end
    [megadata] = deal(cell(length(col),length(lat),length(dtau)));
    ebar = zeros(length(col),length(lat),length(dtau));
    for kk = 1:length(lat)
        for dt = 1:length(dtau)
            for subi = 1:length(col)
                latdata = data(data(:,3)==lat(kk)&data(:,2)==dtau(dt),[1 col(subi) errcol(subi)]);
                ebar(subi,kk,dt) = sqrt(meansqr(latdata(:,end)./latdata(:,end-1)))*100;
                latdata(isnan(latdata(:,2)),:) = [];            
                megadata{subi,kk,dt} = latdata;
            end
        end
    end
else   
    megadata = {};
    lat = NaN;
    dtau = NaN;
    ebar = [];
end
varargout = {megadata, lat, dtau, ebar};
end
