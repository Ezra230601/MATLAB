%% ============================================================
%  主脚本：读广播星历 -> 画 GPS 卫星轨迹 -> 输出卫星位置
%  依赖：项目里已有的 Rinex_get_nav / find_eph / satellite_orbits /
%        glp_ini，以及我刚补上的 weektow2time / four_digit_year /
%        RemoveSpaces / check_t / ecc_anomaly
%  直接点运行即可
% ============================================================
clear; clc; close all;

%% 第 1 步：初始化全局常数（地球半径、光速、引力常数等）
glp_ini;

%% 第 2 步：构造"星座开关"结构体（本次只开 GPS，其他星座跳过不读）
constellations.GPS      = struct('numSat',32,'enabled',1,'indexes',1:32,'PRN',1:32);
constellations.GLONASS  = struct('numSat',24,'enabled',0,'indexes',1:24,'PRN',1:24);
constellations.Galileo   = struct('numSat',36,'enabled',0,'indexes',1:36,'PRN',1:36);
constellations.BeiDou    = struct('numSat',46,'enabled',0,'indexes',1:46,'PRN',1:46);
constellations.QZSS     = struct('numSat',4,'enabled',0,'indexes',1:4,'PRN',1:4);
constellations.SBAS     = struct('numSat',39,'enabled',0,'indexes',1:39,'PRN',1:39);

%% 第 3 步：读 N 文件（广播星历）
navfile = 'BRDC00IGS_R_20260010000_01D_MN(1).rnx';   % 你的 N 文件
[Eph, iono] = Rinex_get_nav(navfile, constellations);
fprintf('星历矩阵大小：%d 行 x %d 列（每列一组星历）\n', size(Eph,1), size(Eph,2));

%% 第 4 步：选卫星
prn = 30;                                   % 改成你想看的卫星号
idx = find(Eph(30,:) == prn);               % 这颗卫星有几组星历
if isempty(idx)
    error('文件里没有 GPS PRN %d，先运行 unique(Eph(30,:)'' ) 看看有哪些卫星', prn);
end

%% 第 5 步：在星历有效时间段内，每 5 分钟算一次卫星位置
t = min(Eph(32,idx)) : 300 : max(Eph(32,idx));
X = nan(size(t)); Y = nan(size(t)); Z = nan(size(t));

for k = 1:length(t)
    icol = find_eph(Eph, prn, t(k));        % 挑出最适合这个时刻的一组星历
    if isempty(icol)
        continue;                           % 这颗卫星此刻不可用（坏星/超出有效期）
    end
    satp = satellite_orbits(t(k), Eph(:,icol));   % 算地固系坐标（米）
    X(k) = satp(1);
    Y(k) = satp(2);
    Z(k) = satp(3);
end

valid = ~isnan(X);
fprintf('成功算出 %d 个时刻的位置\n', sum(valid));

%% 第 6 步：画两张图
figure('Color','w','Position',[80 80 1150 460]);

% 左：三维轨迹 + 地球（单位换成 km）
subplot(1,2,1); hold on; grid on; axis equal;
[Xs,Ys,Zs] = sphere(40);
surf(Xs*6371, Ys*6371, Zs*6371, 'EdgeColor','none', ...
     'FaceAlpha',0.25, 'FaceColor',[0.6 0.8 1]);
plot3(X(valid)/1000, Y(valid)/1000, Z(valid)/1000, 'r.-', ...
      'LineWidth',1.2, 'MarkerSize',8);
plot3(X(find(valid,1,'last'))/1000, Y(find(valid,1,'last'))/1000, ...
      Z(find(valid,1,'last'))/1000, 'kp','MarkerFaceColor','k','MarkerSize',14);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title(sprintf('GPS 卫星 PRN%d 三维轨迹', prn));
legend('地球','轨迹','最终位置','Location','best');

% 右：星下点轨迹（经纬度）
lat = zeros(1,length(t)); lon = zeros(1,length(t));
for k = find(valid)
    [lat(k), lon(k)] = ecef2geod(X(k), Y(k), Z(k));
end
subplot(1,2,2); hold on; grid on;
plot(lon(valid), lat(valid), 'b.-', 'LineWidth',1.2, 'MarkerSize',8);
kend = find(valid,1,'last');
plot(lon(kend), lat(kend), 'kp','MarkerFaceColor','k','MarkerSize',14);
xlabel('经度 (°)'); ylabel('纬度 (°)');
xlim([-180 180]); ylim([-90 90]);
title(sprintf('GPS 卫星 PRN%d 星下点轨迹', prn));
legend('星下点轨迹','最终位置','Location','best');

%% 第 7 步：打印最终位置
fprintf('\n===== 卫星 PRN%d 最终位置 =====\n', prn);
fprintf('ECEF 坐标：X = %.3f km, Y = %.3f km, Z = %.3f km\n', ...
        X(kend)/1000, Y(kend)/1000, Z(kend)/1000);
fprintf('经纬度  ：经度 %.4f°, 纬度 %.4f°\n', lon(kend), lat(kend));

%% ============================================================
%  局部函数：ECEF 直角坐标 -> 大地经纬度（WGS84）
% ============================================================
function [lat, lon, h] = ecef2geod(X, Y, Z)
a  = 6378137.0;
f  = 1/298.257223563;
e2 = f*(2-f);
lon = atan2(Y,X);
p   = sqrt(X^2 + Y^2);
lat = atan2(Z, p*(1-e2));
h   = 0;
for k = 1:10
    N   = a / sqrt(1 - e2*sin(lat)^2);
    h   = p/cos(lat) - N;
    lat = atan2(Z, p*(1 - e2*N/(N+h)));
end
lat = lat*180/pi;
lon = lon*180/pi;
end
