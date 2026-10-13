function  [time, datee,num_sat, sat, sat_types, tow] = RINEX_get_epoch(fid)


time = NaN ;
sat = [ ];
sat_types= [ ] ;

eof = 0;
tow = NaN ;


while (eof == 0)
    lin = fgets(fid) ;
    keywords ={'COMMENT', 'MARKER NAME', 'MARKER NUMBER', 'APPROX POSITION XYZ', 'ANTENNA: DELTA H/E/N'} ;
    answer = [ ];
    s = 1 ;
    
    while (s<=length(keywords) && isempty(answer))
        
        answer = strfind(lin, keywords{s});
        s = s+1 ;
        
    end
    
    while (~isempty(answer) && ~feof(fid))
        
        lin = fgetl(fid) ;
        answer = [ ] ;
        s = 1 ;
        while (s<=length(keywords) && isempty(answer))
            
            answer = strfind(lin, keywords{s});
            s = s+1 ;
            
        end
        
    end
    
    if (feof(fid)==1)
        return
    end
    
    if (~strcmp(lin(1),'>'))   % RINEX v2.XX
        
        if ( strcmp( lin(29),'0' ) || strcmp ( lin(29),'1' ) || strcmp( lin(29),'2' ) )
            
            data = textscan(lin(1:26), '%f%f%f%f%f%f') ;
            year = data{1};
            month = data{2} ;
            day = data{3};
            hour = data{4} ;
            minute = data{5} ;
            second = data{6} ;
            
            year = four_digit_year(year) ;
            [week, tow] = date2gps ([year, month, day, hour, minute, second]);
            [time] = weektow2time(week, tow, 'G');
            
            %number of visible satellites
            [num_sat] = sscanf(lin(30:32),'%d');
            lin = sat_data_sub(lin, 33, 68);
            
            %remove 'blank spaces' and unwanted characters at the end of the string
            lin = RemoveSpaces(lin);
            
            %read additional lines, depending on the number of satellites
            nlines = ceil(num_sat/12);
            for n = 1 : nlines - 1
                lin = [lin sat_data_sub(fgetl(fid), 33, 68)];
                lin = RemoveSpaces(lin);
            end
            
            pos = 1;
            sat = zeros(num_sat,1);
            sat_types = char(32*uint8(ones(num_sat,1))');
            for i = 1 : num_sat
                if (strcmp(lin(pos),' '))
                    type = 'G';
                else
                    type = lin(pos);
                end
                sat_types(i) = type;
                sat(i) = mod((lin(pos+1)-48)*10+(lin(pos+2)-48),160);
                pos = pos + 3;
            end
            
            eof = 1;
        end
    else %RINEX v3.xx
        if (strcmp(lin(32),'0') || strcmp(lin(32),'1') || strcmp(lin(32),'2'))
            
            %save time information
            data   = textscan(lin(2:29),'%f%f%f%f%f%f');
            year   = data{1};
            month  = data{2};
            day    = data{3};
            hour   = data{4};
            minute = data{5};
            second = data{6};
            
            %computation of the GPS time in weeks and seconds of week
            [week, tow] = date2gps([year, month, day, hour, minute, second]);
            [time] = weektow2time(week, tow, 'G');
            
            %number of visible satellites
            [num_sat] = sscanf(lin(33:35),'%d');
            
            eof = 1;
        end
    end
end

datee = [year month day hour minute second];