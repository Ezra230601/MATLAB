%% ============================================================
%  在惯性空间(ECI)里画卫星轨迹 —— 真正的"圆"
%  原理：把每个时刻的 ECEF 坐标绕 Z 轴反向旋转地球自转角，
%        去掉地球自转，就得到卫星在惯性空间的真实轨迹（圆）
% ============================================================
clear; clc; close all;
global glp;
glp_ini;

constellations.GPS      = struct('numSat',32,'enabled',1,'indexes',1:32,'PRN',1:32);
constellations.GLONASS  = struct('numSat',24,'enabled',0,'indexes',1:24,'PRN',1:24);
constellations.Galileo   = struct('numSat',36,'enabled',0,'indexes',1:36,'PRN',1:36);
constellations.BeiDou    = struct('numSat',46,'enabled',0,'indexes',1:46,'PRN',1:46);
constellations.QZSS     = struct('numSat',4,'enabled',0,'indexes',1:4,'PRN',1:4);
constellations.SBAS     = struct('numSat',39,'enabled',0,'indexes',1:39,'PRN',1:39);

navfile = 'BRDC00IGS_R_20260010000_01D_MN(1).rnx';
[Eph, ~] = Rinex_get_nav(navfile, constellations);

prn = 1;
idx = find(Eph(30,:) == prn);
t = min(Eph(32,idx)) : 300 : max(Eph(32,idx));
t0 = t(1);
omega_e = glp.gps.OMEGAE_dot;

X = nan(size(t)); Y = nan(size(t)); Z = nan(size(t));
Xe = X; Ye = Y; Ze = Z;   % ECEF 坐标
for k = 1:length(t)
    icol = find_eph(Eph, prn, t(k));
    if isempty(icol), continue; end
    satp = satellite_orbits(t(k), Eph(:,icol));
    Xe(k) = satp(1); Ye(k) = satp(2); Ze(k) = satp(3);

    % === 关键一步：把 ECEF 旋转回 ECI 惯性系 ===
    theta = omega_e * (t(k) - t0);          % 这段时间地球转了多少弧度
    X(k) = Xe(k)*cos(theta) - Ye(k)*sin(theta);
    Y(k) = Xe(k)*sin(theta) + Ye(k)*cos(theta);
    Z(k) = Ze(k);
end

valid = ~isnan(X);
fprintf('共画出 %d 个点\n', sum(valid));

%% 画图：左 ECEF（跟着地球转，不是圆），右 ECI（惯性空间，是圆）
figure('Color','w','Position',[80 80 1200 500]);

% 左：ECEF
subplot(1,2,1); hold on; grid on; axis equal;
[xs,ys,zs] = sphere(40);
surf(xs*6371, ys*6371, zs*6371, 'EdgeColor','none', ...
     'FaceAlpha',0.25, 'FaceColor',[0.6 0.8 1]);
plot3(Xe(valid)/1000, Ye(valid)/1000, Ze(valid)/1000, 'r.-', 'LineWidth',1.2);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('ECEF 地固系（跟着地球转，不是圆）');
legend('地球','卫星轨迹','Location','best');
view(3);

% 右：ECI 惯性系
subplot(1,2,2); hold on; grid on; axis equal;
surf(xs*6371, ys*6371, zs*6371, 'EdgeColor','none', ...
     'FaceAlpha',0.25, 'FaceColor',[0.6 0.8 1]);
plot3(X(valid)/1000, Y(valid)/1000, Z(valid)/1000, 'b.-', 'LineWidth',1.5);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('ECI 惯性系（去掉地球自转，完美的圆）');
legend('地球','卫星轨迹','Location','best');
view(3);

%% 打印验证：到地心距离应该基本不变
r = sqrt(X(valid).^2 + Y(valid).^2 + Z(valid).^2)/1000;
fprintf('\n到地心距离：平均 %.1f km，最大 %.1f km，最小 %.1f km（变化 %.1f km）\n', ...
    mean(r), max(r), min(r), max(r)-min(r));
