function [corr, distSR_corr] = relativistic_range_error_correction(posR, Satpos)

% SYNTAX:
%   [corr, distSR_corr] = relativistic_range_error_correction(XR, XS);
%
% INPUT:
%   XR = receiver position  (X,Y,Z)
%   XS = satellite position (X,Y,Z)
%
% OUTPUT:
%   corr = relativistic range error correction term
%   distSR_corr = corrected satellite-receiver distance
%
% DESCRIPTION:
%   Computation of the relativistic range error correction term (Shapiro
%   signal propagation delay).

global glp

distR = sqrt(sum(posR.^2 ,1));
distS = sqrt(sum(Satpos.^2 ,2));

posR_mat = posR(:,ones(size(Satpos,1),1))';
distSR = sqrt(sum((Satpos-posR_mat).^2 ,2));

% switch char(Eph(31))
%     case 'G'
%         %GM = goGNSS.GM_GPS;
        GM = 3.986005e14;
%     case 'R'
%         %GM = goGNSS.GM_GLO;
%         GM = 3.9860044e14;
%     case 'E'
%         %GM = goGNSS.GM_GAL;
%         GM = 3.986004418e14;
%     case 'C'
%         %GM = goGNSS.GM_BDS;
%         GM = 3.986004418e14;
%     case 'J'
%         %GM = goGNSS.GM_QZS;
%         GM = 3.986005e14;
%     otherwise
%         fprintf('Something went wrong in ecc_anomaly.m\nUnrecongized Satellite system!\n');
%         %GM = goGNSS.GM_GPS;
%         GM = 3.986005e14;
% end

corr = 2*GM/(glp.vlight^2)*log((distR + distS + distSR)./(distR + distS - distSR));

distSR_corr = distSR + corr;
