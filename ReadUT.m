%%
ps = 100000;
R = 287;
Cp = 1004;
k = R/Cp;
g = 9.81;
H = 7000;
a = 6.371e6;
Omiga = 7.292e-5; 
load('Lnsp');
load('hp');
FileList = dir(strcat('G:\Paper3Data\ERA5_T\','\*.nc'));
Umean = NaN(137,73,47,730);
Tmean = NaN(137,73,47,730);
Thtamean = NaN(137,73,47,730);
for ListN = 1:length(FileList)
    T = ncread(strcat('G:\Paper3Data\ERA5_T\',FileList(ListN).name),'t');
    U = ncread(strcat('G:\Paper3Data\ERA5_U\',FileList(ListN).name),'u');
    hyam = ncread(strcat('G:\Paper3Data\ERA5_T\',FileList(ListN).name),'hyam')';
    hybm = ncread(strcat('G:\Paper3Data\ERA5_T\',FileList(ListN).name),'hybm')';
    if size(T,4) == 732
        T(:,:,:,119:120) = [];
        U(:,:,:,119:120) = [];
    end
    for lat = 1:73
        for time = 1:size(T,4)
            Ulon = NaN(72,137);
            Tlon = NaN(72,137);
            Thtalon = NaN(72,137);
            for lon = 1:72
                u = squeeze(U(lon,lat,:,time));
                t = squeeze(T(lon,lat,:,time));
                plev = squeeze(hyam + hybm.*exp(Lnsp(lon,lat,ListN,time)));
                zlev = -H*log(plev/ps);
                Ulon(lon,:) = interp1(zlev,u,h);
                Tlon(lon,:) = interp1(zlev,t,h);
                Thtalon(lon,:) = Tlon(lon,:).*(ps./p).^k;
            end
            Umean(:,lat,ListN,time) = mean(Ulon,1,'omitmissing');
            Tmean(:,lat,ListN,time) = mean(Tlon,1,'omitmissing');
            Thtamean(:,lat,ListN,time) = mean(Thtalon,1,'omitmissing');
        end
    end
end
clearvars -except Umean Tmean Thtamean
Umean = (Umean(:,:,:,1:2:end-1) + Umean(:,:,:,2:2:end))/2;
Tmean = (Tmean(:,:,:,1:2:end-1) + Tmean(:,:,:,2:2:end))/2;
Thtamean = (Thtamean(:,:,:,1:2:end-1) + Thtamean(:,:,:,2:2:end))/2;
save('UTmean.mat')