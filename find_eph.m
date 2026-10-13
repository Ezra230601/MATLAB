function icol = find_eph(Eph,sat,time)

isat = find(Eph(30,:) ==  sat) ;
n=size(isat,2) ;

if  n==0
    icol=[ ];
    return
end

icol = isat(1) ;
delta =0 ;

%consider BeiDou time (BDT) for BeiDou satellites
if ( strcmp( char(Eph(31)) , 'C' ) )
    delta =14;
    time = time-delta ;
end

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

if fix(abs(dtmin)) -delta > dtmax
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