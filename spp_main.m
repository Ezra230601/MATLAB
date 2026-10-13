%% ============================================================
%  接收机单点定位（SPP）：用 O 文件伪距解算接收机坐标，画散点图
%  流程：读N文件 -> 读O文件GPS伪距 -> 逐历元最小二乘 -> 散点图
% ============================================================
clear; clc; close all;
fprintf('>>> 正在运行【修复版】SPP（钟差列=1，米制钟差） <<<\n');
global glp iono_model;
iono_model = 1;                       % 打开 Klobuchar 电离层模型

%% 第 1 步：初始化 + 读 N 文件
glp_ini;
constellations.GPS      = struct('numSat',32,'enabled',1,'indexes',1:32,'PRN',1:32);
constellations.GLONASS  = struct('numSat',24,'enabled',0,'indexes',1:24,'PRN',1:24);
constellations.Galileo   = struct('numSat',36,'enabled',0,'indexes',1:36,'PRN',1:36);
constellations.BeiDou    = struct('numSat',46,'enabled',0,'indexes',1:46,'PRN',1:46);
constellations.QZSS     = struct('numSat',4,'enabled',0,'indexes',1:4,'PRN',1:4);
constellations.SBAS     = struct('numSat',39,'enabled',0,'indexes',1:39,'PRN',1:39);

[Eph, iono] = Rinex_get_nav('BRDC00IGS_R_20260010000_01D_MN(1).rnx', constellations);

%% 第 2 步：读 O 文件（只取 GPS C1C 伪距）
[epoch_t, prn_cell, pr_cell] = read_obs_gps('abmf0010(1).26o');

%% 第 3 步：接收机近似坐标（来自 O 文件头 APPROX POSITION XYZ）
Xr0 = [2919785.7120; -5383745.0670; 1774604.6920];   % 已知近似坐标
Xr = Xr0;                             % 每历元都从这个初值开始迭代
dtR = 0;                             % 接收机钟差(米等效)
c  = glp.vlight;
omega_e = glp.gps.OMEGAE_dot;
cut_off = 10;                        % 截止高度角 10°

%% 第 4 步：逐历元最小二乘
step = 15;                           % 每 15 个历元算一次（15*30秒=7.5分钟）
lat_sol = []; lon_sol = []; h_sol = []; t_sol = []; n_used = [];

for e = 1:step:length(epoch_t)
    t_rx = epoch_t(e);
    Xr = Xr0; dtR = 0;              % 每个历元重新从初值开始
    for ls_iter = 1:5                % 最小二乘迭代次数
        H = []; b = [];
        for s = 1:length(prn_cell{e})
            prn = prn_cell{e}(s);
            P   = pr_cell{e}(s);
            icol = find_eph(Eph, prn, t_rx);
            if isempty(icol), continue; end

            % --- 迭代求信号发射时刻 t_tx 和卫星位置 ---
            tau = P/c;
            for inner = 1:3
                t_tx = t_rx - tau;
                dtS   = sat_clock_error_correction(t_tx, Eph(:,icol));
                dtrel = relativistic_error_correction(t_tx, Eph(:,icol), [], []);
                tgd   = Eph(28,icol);
                t_tx  = t_rx - P/c - dtS - dtrel + tgd;
                satp  = satellite_orbits(t_tx, Eph(:,icol));
                tau   = t_rx - t_tx;
            end
            dtS_total = dtS + dtrel - tgd;

            % --- 地球自转修正 ---
            satp = earth_rotation_correction(tau, satp, omega_e);

            % --- 几何距离、高度角 ---
            dr = satp - Xr;
            rho = norm(dr);
            [Az, El, ~] = topocent(Xr, satp');
            if El < cut_off, continue; end

            % --- 电离层延迟 ---
            [phi_r, lam_r, ~] = cart2geod(Xr(1), Xr(2), Xr(3));
            sod = mod(t_rx, 86400);
            I = iono_error_correction(phi_r*180/pi, lam_r*180/pi, Az, El, sod, iono);

            % --- 观测方程：b = P - rho - dtR + c*dtS_total - I ---
            % 注意：dtR 用"米"表示（即 c*钟差），所以设计矩阵最后一列是 1
            u = (Xr - satp)/rho;                       % 单位向量
            H(end+1, :) = [u(1), u(2), u(3), 1];
            b(end+1,1)  = P - rho - dtR + c*dtS_total - I;
        end

        if size(H,1) < 4, break; end
        dx = H\b;                       % 最小二乘解
        Xr  = Xr + dx(1:3);
        dtR = dtR + dx(4);             % dx(4) 就是米制钟差改正，直接加

        % ===== 诊断：只打印第一个历元第一次迭代的中间值 =====
        if e == 1 && ls_iter == 1
            fprintf('===== 诊断信息（首历元第1次迭代）=====\n');
            fprintf('首颗卫星 H 行 = [%.6f, %.6f, %.6f, %.1f]\n', H(1,:));
            fprintf('前4个 b 值(m) = [%.2f, %.2f, %.2f, %.2f]\n', b(1:4));
            fprintf('最小二乘 dx = [%.2f, %.2f, %.2f, %.2f]\n', dx);
            fprintf('迭代后 Xr = [%.2f, %.2f, %.2f], dtR=%.2f m\n', Xr, dtR);
            fprintf('用星数 = %d\n', size(H,1));
        end
    end

    if size(H,1) >= 4
        [phi_s, lam_s, h_s] = cart2geod(Xr(1), Xr(2), Xr(3));
        lat_sol(end+1) = phi_s*180/pi;
        lon_sol(end+1) = lam_s*180/pi;
        h_sol(end+1)   = h_s;
        t_sol(end+1)   = t_rx;
        n_used(end+1)  = size(H,1);
    end
end

fprintf('成功解算 %d 个历元\n', length(lat_sol));

%% 第 5 步：散点图
figure('Color','w','Position',[100 100 900 600]);
scatter(lon_sol, lat_sol, 25, t_sol - min(t_sol), 'filled'); hold on;
colorbar; colormap jet;
xlabel('经度 (°)'); ylabel('纬度 (°)');
title('ABMF 站接收机单点定位结果散点图（颜色=当天累计秒）');
axis equal; grid on;

% 画文件头给的已知近似位置作为参考点
[phi0, lam0, ~] = cart2geod(Xr0(1), Xr0(2), Xr0(3));
plot(lam0*180/pi, phi0*180/pi, 'kp', 'MarkerFaceColor','y', 'MarkerSize',15);
legend('解算位置','已知近似位置','Location','best');

%% 打印统计
fprintf('\n===== 定位结果 =====\n');
fprintf('平均经度：%.6f°，平均纬度：%.6f°，平均海拔：%.2f m\n', ...
        mean(lon_sol), mean(lat_sol), mean(h_sol));
fprintf('经纬度散布：%.3f\" / %.3f\"（1度=3600秒）\n', ...
        std(lon_sol)*3600, std(lat_sol)*3600);
fprintf('每个历元平均用星数：%.1f 颗\n', mean(n_used));
