Str = ('G:\Paper3Data\ERA5_Lnsp\Lnsp.nc');
Lnsp(:,:,:) = ncread(Str,'lnsp');
lnsp = NaN(72,73,47,730);
Day = 0;
for year = 1979:2025
    day = 0;
    for month = 1:12
        if rem(year,4) ==0
            fday = 29*2;
        else
            fday = 28*2;
        end
        switch month
            case {1,3,5,7,8,10,12}
                lnsp(:,:,year-1978,(1:31*2)+day) = Lnsp(:,:,(1:31*2)+Day);
                Day = Day+31*2;
                day = day+31*2;
            case {4,6,9,11}
                lnsp(:,:,year-1978,(1:30*2)+day) = Lnsp(:,:,(1:30*2)+Day);
                Day = Day+30*2;
                day = day+30*2;
            case 2
                lnsp(:,:,year-1978,(1:28*2)+day) = Lnsp(:,:,(1:28*2)+Day);
                Day = Day+fday;
                day = day+28*2;
        end
    end
end
Lnsp = lnsp;
clearvars -except Lnsp
save('Lnsp.mat')