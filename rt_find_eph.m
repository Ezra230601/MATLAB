function   Eph_t = rt_find_eph(Eph_in,time,nsat) 

empty_col = zeros(33,1) ;
Eph_t = zeros(33,nsat) ;

for i = 1:nsat
    
    icol=find_eph(Eph_in,i,time) ;
    if (~isempty(icol))
        Eph_t(:,i) = Eph_in(:,icol) ;
    else
        Eph_t(:,i) = empty_col;
    end
    
end

end