function corr = relativistic_error_correction(time, Eph, XS, VS)
%   DESCRIPTION:
%   Computation of the relativistic correction term. From the
%   Interface Specification document revision E (IS-GPS-200E), page 86.

%INPUT:
%   time = GPS time
%   Eph = satellite ephemeris vector
%   XS = satellite position (X,Y,Z)
%   VS = satellite velocity (X,Y,Z)
global glp
if (Eph(4)~=0) %if not using SP3 ephemeris
    roota = Eph(4);
    ecc   = Eph(6);
    
    Ek = ecc_anomaly(time, Eph);
    corr = -4.442807633e-10 * ecc * roota * sin(Ek);
else
    corr = -2*dot(XS,VS)/(glp.vlight^2);
end