load('hp')
load('Lnsp')
ps = 100000;
R = 287;
Cp = 1004;
k = R/Cp;
g = 9.81;
H = 7000;
a = 6.371e6;
FileList = dir(strcat('G:\Paper3Data\ERA5_T\','\*.nc'));
VResidual = NaN(137,73,47,730);
WResidual = NaN(137,73,47,730);
for ListN = 1:length(FileList)
    T = ncread(strcat('G:\Paper3Data\ERA5_T\',FileList(ListN).name),'t');
    V = ncread(strcat('G:\Paper3Data\ERA5_V\',FileList(ListN).name),'v');
    W = ncread(strcat('G:\Paper3Data\ERA5_W\',FileList(ListN).name),'w');
    hyam = ncread(strcat('G:\Paper3Data\ERA5_T\',FileList(ListN).name),'hyam')';
    hybm = ncread(strcat('G:\Paper3Data\ERA5_T\',FileList(ListN).name),'hybm')';
    if size(T,4) == 732
        T(:,:,:,119:120) = [];
        V(:,:,:,119:120) = [];
        W(:,:,:,119:120) = [];
    end
    parfor time = 1:size(T, 4)
        vmean = NaN(137, 73);
        wmean = NaN(137, 73);
        thtamean = NaN(137, 73);
        roumean = NaN(137, 73);
        vthtapie = NaN(137, 73);
        thtaz = NaN(137, 73);
        cosphi = NaN(137, 73);
        rouvthtathtazz = NaN(137, 73);
        cosphivthtathtazphi = NaN(137, 73);
        tlon = NaN(72, 137);
        thtalon = NaN(72, 137);
        vlon = NaN(72, 137);
        wlon = NaN(72, 137);
        roulon = NaN(72, 137);
        for lat = 1:73
            tlon = NaN(72, 137);
            thtalon = NaN(72, 137);
            vlon = NaN(72, 137);
            wlon = NaN(72, 137);
            roulon = NaN(72, 137);
            for lon = 1:72
                plev = hyam + hybm .* exp(Lnsp(lon, lat, ListN, time));
                zlev = -H * log(plev / ps);
                tlon(lon, :) = T(lon, lat, :, time);
                vlon(lon, :) = V(lon, lat, :, time);
                wlon(lon, :) = W(lon, lat, :, time);
                thtalon(lon, :) = tlon(lon, :) .* (ps ./ plev).^k;
                roulon(lon, :) = plev ./ (R * tlon(lon, :));
                wlon(lon, :) = -wlon(lon, :) ./ (roulon(lon, :) * g);
                tlon(lon, :) = interp1(zlev, tlon(lon, :), h);
                vlon(lon, :) = interp1(zlev, vlon(lon, :), h);
                wlon(lon, :) = interp1(zlev, wlon(lon, :), h);
                thtalon(lon, :) = interp1(zlev, thtalon(lon, :), h);
                roulon(lon, :) = interp1(zlev, roulon(lon, :), h);
            end
            cosphi(:, lat) = cos(deg2rad(90 - (2.5) * (lat - 1)));
            vmean(:, lat) = mean(vlon, 1, 'omitmissing');
            wmean(:, lat) = mean(wlon, 1, 'omitmissing');
            thtamean(:, lat) = mean(thtalon, 1, 'omitmissing');
            roumean(:, lat) = mean(roulon, 1, 'omitmissing');
            vpie = vlon - vmean(:, lat)';
            thtapie = thtalon - thtamean(:, lat)';
            vthtapie(:, lat) = mean(vpie .* thtapie, 1, 'omitmissing');
            thtaz(:, lat) = interp1(Midh, diff(thtamean(:, lat)') ./ diff(h), h);
        end
        rouvthtathtaz = roumean .* vthtapie ./ thtaz;
        for lat = 1:73
            rouvthtathtazz(:, lat) = interp1(Midh, diff(rouvthtathtaz(:, lat)') ./ diff(h), h);
        end
        cosphivthtathtaz = cosphi .* vthtapie ./ thtaz;
        for height = 1:137
            cosphivthtathtazphi(height, :) = interp1(deg2rad(88.75:-2.5:-88.75), ...
                diff(cosphivthtathtaz(height, :)) ./ diff(deg2rad(90:-2.5:-90)), ...
                deg2rad(90:-2.5:-90));
        end
        VResidual_local = vmean - rouvthtathtazz ./ roumean;
        WResidual_local = wmean + cosphivthtathtazphi ./ (a * cosphi);
        VResidual(:, :, ListN, time) = VResidual_local;
        WResidual(:, :, ListN, time) = WResidual_local;
    end
end
clearvars -except VResidual WResidual
VResidual = (VResidual(:,:,:,1:2:end-1) + VResidual(:,:,:,2:2:end))/2;
WResidual = (WResidual(:,:,:,1:2:end-1) + WResidual(:,:,:,2:2:end))/2;
save('Residual.mat')