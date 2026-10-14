%% ============================================================
%  导出老师要求的三个结果文件：
%    G*.txt      卫星位置（信号发射时刻的 ECEF 坐标）
%    G*Pr.txt    伪距观测值
%    Pos.txt     单点定位结果
%  直接运行本脚本即可，会在当前目录生成 txt 文件
% ============================================================
clear; clc; close all;
global glp iono_model;
iono_model = 1;
glp_ini;

%% 第 1 步：读 N 文件和 O 文件
constellations.GPS      = struct('numSat',32,'enabled',1,'indexes',1:32,'PRN',1:32);
constellations.GLONASS  = struct('numSat',24,'enabled',0,'indexes',1:24,'PRN',1:24);
constellations.Galileo   = struct('numSat',36,'enabled',0,'indexes',1:36,'PRN',1:36);
constellations.BeiDou    = struct('numSat',46,'enabled',0,'indexes',1:46,'PRN',1:46);
constellations.QZSS     = struct('numSat',4,'enabled',0,'indexes',1:4,'PRN',1:4);
constellations.SBAS     = struct('numSat',39,'enabled',0,'indexes',1:39,'PRN',1:39);

[Eph, iono] = Rinex_get_nav('BRDC00IGS_R_20260010000_01D_MN(1).rnx', constellations);
[epoch_t, prn_cell, pr_cell] = read_obs_gps('abmf0010(1).26o');

%% 第 2 步：初始化
Xr0 = [2919785.7120; -5383745.0670; 1774604.6920];
c  = glp.vlight;
omega_e = glp.gps.OMEGAE_dot;

% 收集所有出现过的 GPS 卫星号
all_prn = unique([prn_cell{:}]);

% 为每颗卫星打开 Gxx.txt 和 GxxPr.txt
sat_fid = containers.Map;      % Gxx.txt 文件句柄
pr_fid  = containers.Map;      % GxxPr.txt 文件句柄
for p = all_prn
    fname = sprintf('G%02d.txt', p);
    fid = fopen(fname, 'w');
    fprintf(fid, '# satpos  信号发射时刻的 ECEF 位置\n');
    fprintf(fid, '# sow_s          x_m            y_m            z_m\n');
    sat_fid(num2str(p)) = fid;

    fname2 = sprintf('G%02dPr.txt', p);
    fid2 = fopen(fname2, 'w');
    fprintf(fid2, '# prange  该频点原始伪距\n');
    fprintf(fid2, '# sow_s          P_m\n');
    pr_fid(num2str(p)) = fid2;
end

% 打开 Pos.txt
pos_fid = fopen('Pos.txt', 'w');
fprintf(pos_fid, '# spp  接收机天线相位中心的单点定位结果\n');
fprintf(pos_fid, '# sow_s          x_m            y_m            z_m\n');

%% 第 3 步：逐历元解算并写文件
% 注意：卫星位置 satp 只与发射时刻 t_tx 有关，而 t_tx = t_rx - P/c - 钟差项，
%       不依赖接收机坐标 Xr。所以每个历元对每颗卫星只需算一次 satp，
%       然后用预存的 satp 做最小二乘迭代。写入文件也只写一次。
for e = 1:length(epoch_t)
    t_rx = epoch_t(e);
    sow  = mod(t_rx, 604800);          % 周内秒

    % --- 预计算每颗卫星的位置和钟差（与接收机坐标无关） ---
    n_sat = length(prn_cell{e});
    pre = cell(1, n_sat);              % 每元素 = [prn, P, satp, dtS_total] 或 []
    for s = 1:n_sat
        prn = prn_cell{e}(s);
        P   = pr_cell{e}(s);
        icol = find_eph(Eph, prn, t_rx);
        if isempty(icol), continue; end

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
        satp = earth_rotation_correction(tau, satp, omega_e);
        pre{s} = struct('prn', prn, 'P', P, 'satp', satp, 'dtS_total', dtS_total);

        % 写 Gxx.txt：每颗卫星每历元只写一次（信号发射时刻的 ECEF 位置）
        fprintf(sat_fid(num2str(prn)), '%12.3f  %14.3f  %14.3f  %14.3f\n', ...
            sow, satp(1), satp(2), satp(3));
        % 写 GxxPr.txt：原始伪距
        fprintf(pr_fid(num2str(prn)), '%12.3f  %14.3f\n', sow, P);
    end

    % --- 最小二乘迭代（用预存的 satp，不再重算） ---
    Xr = Xr0; dtR = 0;
    for ls_iter = 1:5
        H = []; b = [];
        for s = 1:n_sat
            if isempty(pre{s}), continue; end
            satp      = pre{s}.satp;
            P         = pre{s}.P;
            dtS_total = pre{s}.dtS_total;

            dr = satp - Xr; rho = norm(dr);
            [Az, El, ~] = topocent(Xr, satp');
            if El < 10, continue; end
            [phi_r, lam_r, ~] = cart2geod(Xr(1), Xr(2), Xr(3));
            sod = mod(t_rx, 86400);
            I = iono_error_correction(phi_r*180/pi, lam_r*180/pi, Az, El, sod, iono);

            u = (Xr - satp)/rho;
            H(end+1,:) = [u(1) u(2) u(3) 1];
            b(end+1,1) = P - rho - dtR + c*dtS_total - I;
        end
        if size(H,1) < 4, break; end
        dx = H\b;
        Xr = Xr + dx(1:3); dtR = dtR + dx(4);
    end

    % 写 Pos.txt：这个历元的接收机位置
    if size(H,1) >= 4
        fprintf(pos_fid, '%12.3f  %14.3f  %14.3f  %14.3f\n', sow, Xr(1), Xr(2), Xr(3));
    end
end

%% 第 4 步：关闭所有文件
for p = all_prn
    fclose(sat_fid(num2str(p)));
    fclose(pr_fid(num2str(p)));
end
fclose(pos_fid);

fprintf('导出完成！已生成：\n');
fprintf('  - Gxx.txt（卫星位置，共 %d 颗卫星）\n', length(all_prn));
fprintf('  - GxxPr.txt（伪距观测值）\n');
fprintf('  - Pos.txt（单点定位结果，共 %d 个历元）\n', length(epoch_t));
