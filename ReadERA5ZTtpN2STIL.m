%%
ps = 100000;
R = 287;
Cp = 1004;
g = 9.81;
H = 7000;
FileList = dir(strcat('G:\Paper3Data\ERA5_T\','\*.nc'));
load('Lnsp');
TBh = -10000:100:10000;
load('hp','h');
midh = (h(2:end) + h(1:end-1))/2;
%%
Ztp = NaN(72,73,length(FileList),730);
Ttp = NaN(72,73,length(FileList),730);
N2 = NaN(137,73,length(FileList),730);
STIL = NaN(72,73,length(FileList),730);
N2TB = NaN(201,73,length(FileList),730);
%%
LnspSize = size(Lnsp);
for ListN = 1:length(FileList)
    filePath = strcat('G:\Paper3Data\ERA5_T\',FileList(ListN).name);
    ERA5T_full = ncread(filePath,'t');
    hyam = ncread(filePath,'hyam')';
    hybm = ncread(filePath,'hybm')';
    if size(ERA5T_full,4) == 732
        ERA5T_full(:,:,:,119:120) = [];
    end
    Lnsp_current = Lnsp(:,:,ListN,:);
    ERA5T_current = ERA5T_full;
    [~, ~, ~, nTime] = size(ERA5T_current);
    Ztp_temp = NaN(72,73,nTime);
    Ttp_temp = NaN(72,73,nTime);
    N2_temp = NaN(137,73,nTime);
    STIL_temp = NaN(72,73,nTime);
    N2TB_temp = NaN(201,73,nTime);
    parfor lat = 1:73
       abs_lat = abs(90 - (lat-1)*2.5);
       if abs_lat >= 60          
           h_min = 3000;       
           h_max = 15000;      
       elseif abs_lat >= 30    
           h_min = 5000;       
           h_max = 20000;    
       else                     
           h_min = 8000;     
           h_max = 22000;    
       end
        Ztp_lat = NaN(72,nTime);
        Ttp_lat = NaN(72,nTime);
        N2_lat = NaN(137,nTime);
        STIL_lat = NaN(72,nTime);
        N2TB_lat = NaN(201,nTime);
        for time = 1:nTime
            N2lon = NaN(72,137);
            N2TBlon = NaN(72,201);
            for lon = 1:72
                p = hyam + hybm .* exp(Lnsp_current(lon,lat,1,time));
                zlog = -H * log(p/ps);
                t_raw = squeeze(ERA5T_current(lon,lat,:,time));
                t = interp1(zlog, t_raw, h);
                tz = interp1(midh, diff(t)./diff(h), h);
                n2 = R/H * (tz + R/Cp * t/H);
                ztp = NaN;
                ttp = NaN;
                for height = 137:-1:1
                    if (h(height) < h_min) || (h(height) > h_max)
                        continue;
                    end
                    if height == 1
                        continue;
                    end
                    if (tz(height-1) > -2e-3) && (tz(height) < -2e-3)
                        ztp = interp1(tz(height-1:height), h(height-1:height), -2e-3);             
                        z_check = ztp+200:200:ztp+2000;
                        tz_check = interp1(h, tz, z_check, 'linear', 'extrap');
                        if all(cumsum(tz_check)./(1:length(tz_check)) >= -2e-3)
                            ttp = interp1(h, t, ztp);
                            break;
                        else
                            ztp = NaN;
                            ttp = NaN;
                        end
                    end
                end
                N2lon(lon,:) = n2;
                if ~isnan(ztp)
                    Ztp_lat(lon,time) = ztp;
                    Ttp_lat(lon,time) = ttp;
                    N2TBlon(lon,:) = interp1(h-ztp, n2, TBh);
                    STIL_lat(lon,time) = max(N2TBlon(lon,101:131));
                end
            end
            N2_lat(:,time) = mean(N2lon, 1, 'omitnan');
            if ~all(isnan(Ztp_lat(:,time)))
                N2TB_lat(:,time) = mean(N2TBlon, 1, 'omitnan');
            end
        end
        Ztp_temp(:,lat,:) = Ztp_lat;
        Ttp_temp(:,lat,:) = Ttp_lat;
        N2_temp(:,lat,:) = N2_lat;
        STIL_temp(:,lat,:) = STIL_lat;
        N2TB_temp(:,lat,:) = N2TB_lat;
    end
    Ztp(:,:,ListN,:) = Ztp_temp;
    Ttp(:,:,ListN,:) = Ttp_temp;
    N2(:,:,ListN,:) = N2_temp;
    STIL(:,:,ListN,:) = STIL_temp;
    N2TB(:,:,ListN,:) = N2TB_temp;
    fprintf('Completed ListN = %d/%d\n', ListN, length(FileList));
end
clearvars -except Ztp Ttp N2 N2TB STIL
Ztp = (Ztp(:,:,:,1:2:end-1) + Ztp(:,:,:,2:2:end))/2;
Ttp = (Ttp(:,:,:,1:2:end-1) + Ttp(:,:,:,2:2:end))/2;
N2 = (N2(:,:,:,1:2:end-1) + N2(:,:,:,2:2:end))/2;
N2TB = (N2TB(:,:,:,1:2:end-1) + N2TB(:,:,:,2:2:end))/2;
STIL = (STIL(:,:,:,1:2:end-1) + STIL(:,:,:,2:2:end))/2;
save('N2ZTtpSTILERA.mat')