function glp_out = glp_ini

global glp   %SHL

global cut_off  n_sys


n_sys = 1 ;
cut_off = 10 ;  % elevation cutoff angle, 10 deg
glp.vlight= 299792458 ;

% CONSTELLATION REF -----------------------------------------------
% CRS parameters, according to each GNSS system CRS definition
% (ICD document in brackets):
%
% *_GPS --> WGS-84   (IS-GPS200E)
% *_GLO --> PZ-90    (GLONASS-ICD 5.1)
% *_GAL --> GTRF     (Galileo-ICD 1.1)
% *_BDS --> CGCS2000 (BeiDou-ICD 1.0)
% *_QZS --> GRS80   (IS-QZSS 1.5D)

gps.elA = 6378137;                                % GPS (WGS-84)      Ellipsoid semi-major axis [m]
glo.elAGLO = 6378136;                          % GLONASS (PZ-90)   Ellipsoid semi-major axis [m]
gal.elAGAL = 6378137;                           % Galileo (GTRF)    Ellipsoid semi-major axis [m]
bds.elABDS = 6378136;                          % BeiDou (CGCS2000) Ellipsoid semi-major axis [m]
qzs.elAQZS = 6378137;                          % QZSS (WGS-84)     Ellipsoid semi-major axis [m]

gps.elF = 1/298.257222101;                        % GPS (WGS-84)      Ellipsoid flattening
glo.elF = 1/298.257222101;                          % GLONASS (PZ-90)   Ellipsoid flattening
gal.elF = 1/298.257222101;                   % Galileo (GTRF)    Ellipsoid flattening
bds.elF = 1/298.257222101;                  % BeiDou (CGCS2000) Ellipsoid flattening
qzs.elF = 1/298.257222101;                  % QZSS (WGS-84)     Ellipsoid flattening

gps.elE = sqrt(1-(1-gps.elF)^2);             % GPS (WGS-84)      Eccentricity
glo.elE = sqrt(1-(1-glo.elF )^2);       % GLONASS (PZ-90)   Eccentricity
gal.elE = sqrt(1-(1-gal.elF)^2);         % Galileo (GTRF)    Eccentricity
bds.elE = sqrt(1-(1-bds.elF)^2);       % BeiDou (CGCS2000) Eccentricity
qzs.elE = sqrt(1-(1-qzs.elF)^2);         % QZSS (WGS-84)     Eccentricity

gps.GM = 3.986005e14;                        % GPS     Gravitational constant * (mass of Earth) [m^3/s^2]
glo.GM = 3.9860044e14;                   % GLONASS Gravitational constant * (mass of Earth) [m^3/s^2]
gal.GM = 3.986004418e14;                   % Galileo Gravitational constant * (mass of Earth) [m^3/s^2]
bds.GM = 3.986004418e14;                  % BeiDou  Gravitational constant * (mass of Earth) [m^3/s^2]
qzs.GM = 3.986005e14;                        % QZSS    Gravitational constant * (mass of Earth) [m^3/s^2]
% (NOTE: these values are not actually called from goGNSS.m
%        by ecc_anomaly.m for computation time reasons; if
%        it's needed to change them, please update also
%        the values in ecc_anomaly.m)

gps.OMEGAE_dot = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
glo.OMEGAE_dot = 7.292115e-5;                      % GLONASS Angular velocity of the Earth rotation [rad/s]
gal.OMEGAE_dot = 7.2921151467e-5;              % Galileo Angular velocity of the Earth rotation [rad/s]
bds.OMEGAE_dot = 7.292115e-5;                     % BeiDou  Angular velocity of the Earth rotation [rad/s]
qzs.OMEGAE_dot = 7.2921151467e-5;             % QZSS    Angular velocity of the Earth rotation [rad/s]

glo.J2= 1.0826257e-3;                        % GLONASS second zonal harmonic of the geopotential

glp.pi = 3.1415926535898;                   % pi value used for orbit computation
glp.circle_rad = 2*glp.pi ;              % 2 pi (NOTE: this value is not actually called from goGNSS.m
%             for computation time reasons; if it's needed to
%             change it, please update also ecc_anomaly.m and
%             satellite_orbits.m)



%                           CONSTELLATION SPECIFIC ------------------------------------------

gps.f1 = 1575.420;  % GPS [MHz]
gps.f2 = 1227.600;  %
gps.f5 = 1176.450;  %

gps.fr = [gps.f1 gps.f2 gps.f5]*1e6 ;         % GPS carriers frequencies [Hz]
gps.lambda = glp.vlight./gps.fr ;               % GPS carriers wavelengths [m] LAMBDAG
gps.alpha1  = gps.fr(1)^2/(gps.fr(1)^2 - gps.fr(2)^2); % GPS iono-free combination parameter
gps.alpha2 = gps.fr(2)^2/(gps.fr(1)^2 - gps.fr(2)^2); % GPS iono-free combination parameter
gps.alphaT = 77; % round(gps.f1/10.23/2);                   % GPS iono-free combination parameter
gps.alphaN = 60 ;% round(gps.f2/10.23/2);                   % GPS iono-free combination parameter


glo.f1_base  = 1602.000;  % GLONASS [MHz]
glo.f2_base  = 1246.000;  %
glo.f1_delta = 0.5625;    %
glo.f2_delta = 0.4375;    %
glo.fch = 6:-1:-7;

glo.fr_base = [glo.f1_base glo.f2_base];                             % GLONASS carriers base frequencies [Hz]
glo.fr_delta = [glo.f1_delta glo.f2_delta];                              % GLONASS carriers delta frequencies [Hz/n]
glo.fr1 = glo.fch' .* glo.fr_delta(1) + glo.fr_base(1);
glo.fr2 = glo.fch' .* glo.fr_delta(2) + glo.fr_base(2);
glo.fr = [glo.fr1 glo.fr2]*1e6;                                                 % GLONASS carriers frequencies [Hz]
glo.lambdaR = glp.vlight ./ glo.fr;                                         % GLONASS carriers wavelengths [m]
glo.alpha1 = glo.fr(:,1).^2./(glo.fr(:,1).^2 - glo.fr(:,2).^2);      % GLONASS iono-free combination parameter
glo.alpha2 = glo.fr(:,2).^2./(glo.fr(:,1).^2 - glo.fr(:,2).^2);      % GLONASS iono-free combination parameter
glo.alphaT = 9;                                                                     % GLONASS iono-free combination parameter
glo.alphaN = 7;                                                                     % GLONASS iono-free combination parameter


gal.f1  = gps.f1;     % Galileo [MHz]
gal.f5a = gps.f5;     %
gal.f5b = 1207.140;       %
gal.f5  = 1191.795;       %
gal.f6  = 1278.750;       %

gal.fr = [gal.f1 gal.f5a gal.f5b gal.f5 gal.f6]*1e6;            % Galileo carriers frequencies [Hz]
gal.lambda = glp.vlight  ./ gal.fr   ;                                  % Galileo carriers wavelengths [m]
gal.alpha1 = gal.fr(1)^2/(gal.fr(1)^2 - gal.fr(2)^2);          % Galileo iono-free combination parameter
gal.alpha2 = gal.fr(2)^2/(gal.fr(1)^2 - gal.fr(2)^2);          % Galileo iono-free combination parameter
gal.alphaT = 154;  % round(gal.f1/10.23);         % Galileo iono-free combination parameter
gal.alphaN = 115;  % round(gal.f5a/10.23);       % Galileo iono-free combination parameter


bds.f1  = 1589.740;       % BeiDou [MHz]
bds.f2  = 1561.098;       %
bds.f5b = gal.f5b;           %
bds.f6  = 1268.520;       %

bds.fr = [bds.f2 bds.f5b bds.f6 bds.f1]*1e6;                                 % BeiDou carriers frequencies [Hz]
bds.lambda = glp.vlight ./ bds.fr;                                       % BeiDou carriers wavelengths [m]
bds.alpha1 = bds.fr(1)^2/(bds.fr(1)^2 - bds.fr(2)^2);                       % BeiDou iono-free combination parameter
bds.alpha2 = bds.fr(2)^2/(bds.fr(1)^2 - bds.fr(2)^2);                 % BeiDou iono-free combination parameter
bds.alphaT = 763;% round(bds.f2/2.046);                             % BeiDou iono-free combination parameter
bds.alphaN = 590;% round(bds.f5b/2.046);                            % BeiDou iono-free combination parameter




qzss.f1 = gps.f1 ;      % QZSS [MHz]
qzss.f2 = gps.f2 ;     %
qzss.f5 = gps.f5 ;      %
qzss.f6 = gal.f6 ;      %
qzss.fr = [qzss.f1 qzss.f2 qzss.f5 qzss.f6]*1e6;
qzss.lambda = glp.vlight ./ qzss.fr;                      % QZSS carriers wavelengths [m]
qzss.alpha1 = qzss.fr(1)^2/(qzss.fr(1)^2 - qzss.fr(2)^2); % QZSS iono-free combination parameter
qzss.alpha2 = qzss.fr(2)^2/(qzss.fr(1)^2 - qzss.fr(2)^2); % QZSS iono-free combination parameter
qzss.alphaT = 77;% round(qzss.f1/10.23/2);                   % QZSS iono-free combination parameter
qzss.alphaN = 60;% round(qzss.f2/10.23/2);                   % QZSS iono-free combination parameter



sbas.f1 = gps.f1;      % SBAS [MHz]
sbas.f5 = gps.f5;      %
sbas.fr = [sbas.f1 sbas.f5]*1e6;                           % SBAS carriers frequencies [Hz]
sbas.lambda = glp.vlight ./ sbas.fr;                      % SBAS carriers wavelengths [m]
sbas.alpha1 = sbas.fr(1)^2/(sbas.fr(1)^2 - sbas.fr(2)^2); % SBAS iono-free combination parameter
sbas.alpha2 = sbas.fr(2)^2/(sbas.fr(1)^2 - sbas.fr(2)^2); % SBAS iono-free combination parameter
sbas.alphaT = 154;% round(gal.f1/10.23);                    % SBAS iono-free combination parameter
sbas.alphaN = 115;% round(gal.f5a/10.23);                   % SBAS iono-free combination parameter


glp.gps =gps ;
glp.qzs = qzs ;
glp.bds =bds ;
glp.gal = gal ;
glp.glo = glo ;
glp.qzss = qzss;
glp.sbas = sbas ;

% GNSS ID
glp.ID_GPS = 1; % Id of GPS constellation for goGPS internal use
glp.ID_GLONASS = 2; % Id of GLONASS constellation for goGPS internal use
glp.ID_GALILEO = 3; % Id of Galileo constellation for goGPS internal use
glp.ID_BEIDOU  = 4; % Id of BeiDou constellation for goGPS internal use
glp.ID_QZSS    = 5; % Id of QZSS constellation for goGPS internal use
glp.ID_SBAS    = 6; % Id of SBAS constellation for goGPS internal use



glp_out = glp ;


%% global parameters definition
% Fisher test table
%-------------------------------------------------------------------------------
global FTABLE;
FTABLE=finv(0.9995,1,1:200)';

global elea
elea = 10 ;
global flag_doppler_cs
global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_vel
global sigmaq_cod1 sigmaq_cod2 sigmaq_ph sigmaq0_N sigmaq_dtm
global min_nsat  snr_threshold cs_threshold  snr_a snr_0 snr_1 snr_A order o1 o2 o3
global h_antenna
% global tile_header tile_georef
% global master_ip master_port ntrip_user ntrip_pw ntrip_mountpoint
% global nmea_init
% global COMportR  COMportR ='NA'
global IAR_method P0 mu flag_auto_mu flag_default_P0

flag_doppler_cs = 0 ; %  use doppler for cycle slip
min_nsat = 2 ;
snr_threshold = 0 ; %this value could change due to the input
snr_a = 30;  snr_0 = 10;  snr_1 = 50; snr_A = 30;

order = 2 ; %
o1 = order;  o2 = order*2;  o3 = order*3;

sigmaq0 = 1 ;
sigmaq_vE = 0.25 ;   sigmaq_vN =0.25 ;   sigmaq_vU = 0.01 ; sigmaq_vel = 1 ;
sigmaq_cod1 = 9 ; sigmaq_cod2 = 0.16 ;
sigmaq_ph =   9.0000e-04 ; sigmaq0_N = 1000;  sigmaq_dtm=1e30 ;

IAR_method = 0 ;
P0= 0.001 ;
mu = 0.5 ;
flag_auto_mu = 1 ;  flag_default_P0 = 1;
h_antenna = 1;
cs_threshold= 3 ;

global clock_delay_thresh
global cond_num_threshold

clock_delay_thresh = 100;
cond_num_threshold = 1e6; %threshold on the condition number on the
% eigenvalues of the N matrix (least squares)

global Xhat_t_t
Xhat_t_t =[ ];

global azR elR distR
%azimuth, elevation and distance of satellites with respect to the ROVER
azR = []; elR = []; distR = [];

global azM elM distM
%azimuth, elevation and distance of satellites with respect to the MASTER
azM = []; elM = []; distM = [];

%estimation error covariance matrix at time t
global Cee
Cee = [];

global geoid
geoid.grid = 0;
geoid.cellsize = 0;
geoid.Xll = 0;
geoid.Yll = 0;
geoid.ncols = 0;
geoid.nrows = 0;



global max_code_residual flag_outlier flag_outlier_OLOO SPP_threshold
% Enable outlier detection (0/1)
flag_outlier = 0 ;
% Apply OLOO for outlier detection (0/1)
flag_outlier_OLOO = 0 ;
% Threshold on code LS estimation error [m]
pp_spp_thr = 4 ;
% Threshold on maximum residual of code obs [m]
max_code_residual = 30 ;
SPP_threshold = 4 ;



end



