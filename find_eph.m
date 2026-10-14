function icol = find_eph(Eph,sat,time)

isat = find(Eph(30,:) ==  sat) ;
n=size(isat,2) ;

if  n==0
    icol=[ ];
    return
end

icol = isat(1) ;

% 注：weektow2time 已经为北斗把 +14 秒加进 Eph(32)，
%     所以 Eph(32) 对所有星座都是"GPS 连续秒"，
%     这里直接和 time（GPS 连续秒）比较即可，不再做北斗 -14 调整。
time_eph = Eph(32,icol);
dtmin = time_eph - time ;

for i= isat
    time_eph = Eph(32,i) ;
    dt = time_eph - time ;
    if abs(dt)<abs(dtmin)
        icol = i ;
        dtmin = dt;
    end
end

fit_interval = Eph(29,icol) ;
if fit_interval ~=0
    dtmax = fit_interval*3600/2 ;
else
    switch char(Eph(31,icol))
        case 'R' % GLONASS
            dtmax = 900;
        case 'J' % QZSS
            dtmax = 3600 ;
        otherwise
            dtmax = 7200;
    end
end

if fix(abs(dtmin)) > dtmax
    icol = [ ];
    return ;
end

%check satellite health
%the second and third conditions are temporary (QZSS and Galileo health flag is kept on for tests)
if (Eph(27,icol) ~= 0 && ~strcmp(char(Eph(31,icol)),'J') && ~strcmp(char(Eph(31,icol)),'E'))
    icol = [];
    return
end

end