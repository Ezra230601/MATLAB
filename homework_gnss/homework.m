%% homework.m  ——  GNSS 作业：PRN=14 的 L1 观测 + ECEF 位置
% 任务 1：L1 载波相位 / 多普勒 / 信噪比 随时间变化（3 张图）
% 任务 2：X、Y、Z 坐标分量随时间变化（1 张图，3 个子图）
% 任务 3：X、Y、Z 三维轨迹（1 张图），标明历元顺序
% 输出 5 张图 + 1 个"时间+数据"表 (homework_table.txt)
% 时间系统：GPST（GPS 时）

clear; clc; close all;

%% ---------- 配置 ----------
prn       = 14;                                     % 学号末两位
obs_file  = 'abmf0010(1).26o';                       % RINEX 3 观测文件
nav_file  = 'BRDC00IGS_R_20260010000_01D_MN(1).rnx'; % RINEX 3 广播星历
out_table = 'homework_table.txt';                    % 输出数据表

glp = glp_ini;                                      % 加载常数（omega_e、c 等）
c       = glp.vlight;
omega_e = glp.gps.OMEGAE_dot;

%% ---------- 第 1 步：读广播星历 ----------
fprintf('读取广播星历 %s ...\n', nav_file);
constellations.GPS      = struct('numSat',32,'enabled',1,'indexes',1:32,'PRN',1:32);
constellations.GLONASS  = struct('numSat',24,'enabled',0,'indexes',1:24,'PRN',1:24);
constellations.Galileo   = struct('numSat',36,'enabled',0,'indexes',1:36,'PRN',1:36);
constellations.BeiDou    = struct('numSat',46,'enabled',0,'indexes',1:46,'PRN',1:46);
constellations.QZSS     = struct('numSat',4,'enabled',0,'indexes',1:4,'PRN',1:4);
constellations.SBAS     = struct('numSat',39,'enabled',0,'indexes',1:39,'PRN',1:39);
[Eph, iono] = Rinex_get_nav(nav_file, constellations);

%% ---------- 第 2 步：读 RINEX 3 观测文件，提取 PRN14 的 L1 观测 ----------
fprintf('读取观测文件 %s，提取 G%02d ...\n', obs_file, prn);
fid = fopen(obs_file);  if fid<0, error('打不开观测文件'); end

% --- 跳过文件头，找到 "END OF HEADER" ---
while ~feof(fid)
    lin = fgetl(fid);
    if contains(lin, 'END OF HEADER'), break; end
end

% --- 逐历元读取 ---
T      = [];     % GPST 连续秒（自 1980-01-06 起的秒数）
sod    = [];     % 日内秒 (0..86399.999)，用作画图横轴
tdstr  = {};     % "YYYY-MM-DD HH:MM:SS" 字符串
L1phas = [];     % L1 载波相位（cycles）
Dop    = [];     % L1 多普勒（Hz）
SNR    = [];     % L1 信噪比（dB-Hz）

buf_T = []; buf_sod = []; buf_td = {}; buf_L = []; buf_D = []; buf_S = [];

while ~feof(fid)
    lin = fgetl(fid);
    if ~ischar(lin) || isempty(lin), continue; end
    if lin(1) ~= '>', continue; end

    % 解析历元首行：> 2026 01 01 00 00  0.0000000  0 46
    data = sscanf(lin(2:end), '%f');
    if length(data) < 8, continue; end
    flag = data(7);
    if flag ~= 0, continue; end              % 跳过异常历元
    nsat = data(8);

    % 时间转 GPST 连续秒 + 日内秒 + 日期字符串
    [week, tow] = date2gps(data(1:6)');
    t_cont = week * 604800 + tow;
    sec_of_day = mod(t_cont, 86400);
    date_str = sprintf('%04d-%02d-%02d %02d:%02d:%05.2f', ...
        data(1), data(2), data(3), data(4), data(5), data(6));

    % --- 读取这一历元的所有卫星行，找 G{prn} ---
    found = false;
    for k = 1:nsat
        lin2 = fgetl(fid);
        if ~ischar(lin2) || length(lin2) < 5, continue; end
        if ~strcmp(lin2(1:1), 'G'), continue; end
        svprn = sscanf(lin2(2:3), '%d');
        if svprn ~= prn, continue; end

        % 按固定字符位提取 C1C / L1C / D1C / S1C
        C1C = parse_field(lin2,  5);    % 第 1 个 obs (C1C)
        L1C = parse_field(lin2, 21);    % 第 2 个 obs (L1C)
        D1C = parse_field(lin2, 37);    % 第 3 个 obs (D1C)
        S1C = parse_field(lin2, 53);    % 第 4 个 obs (S1C)

        % 存这一历元
        buf_T(end+1)   = t_cont;
        buf_sod(end+1) = sec_of_day;
        buf_td{end+1} = date_str;
        buf_L(end+1)  = L1C;
        buf_D(end+1)  = D1C;
        buf_S(end+1)  = S1C;
        found = true;
        break;     % 一个历元内 G{prn} 只出现一次
    end
    % 如果本历元没找到 G{prn}，"for k" 已读完 nsat 行；如果找到了，
    % 上面 break 跳出 for，但 fgetl 已读到 G{prn} 行，其余 nsat-k 行未读，
    % 需要继续 fgetl 把它们消耗掉，保证下次循环对齐到下一个 '>' 历元首行
    if found
        for j = k+1:nsat, fgetl(fid); end
    end
end
fclose(fid);

T      = buf_T; sod = buf_sod; tdstr = buf_td;
L1phas = buf_L; Dop = buf_D; SNR = buf_S;
fprintf('共读出 %d 个历元的 G%02d 观测。\n', length(T), prn);

if isempty(T), error('没找到 G%02d 的观测数据', prn); end

%% ---------- 第 3 步：计算卫星 ECEF 位置（X, Y, Z） ----------
fprintf('计算 G%02d 的 ECEF 位置 ...\n', prn);
X = nan(size(T)); Y = X; Z = X;
for k = 1:length(T)
    t_rx = T(k);
    icol = find_eph(Eph, prn, t_rx);
    if isempty(icol), continue; end

    % 用伪距估算信号传播时间（用 C1C 当伪距；缺失则用 0.07s 估值）
    if ~isnan(L1phas(k)) && ~isempty(L1phas(k))
        % 没存 C1C，用粗略估值 26000 km / c
        P = 2.6e7;
    else
        P = 2.6e7;
    end
    tau = P / c;
    for inner = 1:3
        t_tx  = t_rx - tau;
        dtS   = sat_clock_error_correction(t_tx, Eph(:,icol));
        dtrel = relativistic_error_correction(t_tx, Eph(:,icol), [], []);
        tgd   = Eph(28,icol);
        t_tx  = t_rx - P/c - dtS - dtrel + tgd;
        satp = satellite_orbits(t_tx, Eph(:,icol));
        tau  = t_rx - t_tx;
    end
    satp = earth_rotation_correction(tau, satp, omega_e);
    X(k) = satp(1); Y(k) = satp(2); Z(k) = satp(3);
end

%% ---------- 第 4 步：5 张图 ----------
% 在时间间隙处插入 NaN，让线段断开（卫星不可见时段不连线）
gap_idx = find(diff(T) > 60);   % 时间跳变 > 60s 的位置
insert_NaN = @(v) [v(1); repmat([NaN; v(2:end)], size(v,1), 1)];  % 占位
% 更稳的做法：逐个间隙插入 NaN
hhmm  = sod / 3600;
xx = X; yy = Y; zz = Z; ll = L1phas; dd = Dop; ss = SNR; tt = hhmm;
for g = gap_idx(end:-1:1)   % 从后往前插，避免索引乱
    tt = [tt(1:g), NaN, tt(g+1:end)];
    xx = [xx(1:g), NaN, xx(g+1:end)];
    yy = [yy(1:g), NaN, yy(g+1:end)];
    zz = [zz(1:g), NaN, zz(g+1:end)];
    ll = [ll(1:g), NaN, ll(g+1:end)];
    dd = [dd(1:g), NaN, dd(g+1:end)];
    ss = [ss(1:g), NaN, ss(g+1:end)];
end

% --- 图 1：L1 载波相位 ---
f1 = figure('Name','Fig 1: L1 carrier phase','Position',[100 100 900 400]);
plot(tt, ll, 'b-', 'LineWidth', 1.2); grid on;
xlabel('GPST 时间 (h, 2026-01-01)');
ylabel('L1 载波相位 \phi (cycles)');
title(sprintf('G%02d  L1 载波相位随时间变化', prn));

% --- 图 2：多普勒 ---
f2 = figure('Name','Fig 2: Doppler','Position',[120 100 900 400]);
plot(tt, dd, 'r-', 'LineWidth', 1.2); grid on;
xlabel('GPST 时间 (h, 2026-01-01)');
ylabel('L1 多普勒 f_d (Hz)');
title(sprintf('G%02d  L1 多普勒频移随时间变化', prn));

% --- 图 3：信噪比 ---
f3 = figure('Name','Fig 3: SNR','Position',[140 100 900 400]);
plot(tt, ss, 'g-', 'LineWidth', 1.2); grid on;
xlabel('GPST 时间 (h, 2026-01-01)');
ylabel('L1 信噪比 S/N_0 (dB-Hz)');
title(sprintf('G%02d  L1 信噪比随时间变化', prn));

% --- 图 4：X, Y, Z 分量 ---
f4 = figure('Name','Fig 4: X Y Z components','Position',[160 100 900 700]);
subplot(3,1,1); plot(tt, xx/1000, 'b-', 'LineWidth', 1.2); grid on;
ylabel('X (km)'); title(sprintf('G%02d  ECEF 坐标分量随时间变化', prn));
subplot(3,1,2); plot(tt, yy/1000, 'r-', 'LineWidth', 1.2); grid on;
ylabel('Y (km)');
subplot(3,1,3); plot(tt, zz/1000, 'k-', 'LineWidth', 1.2); grid on;
xlabel('GPST 时间 (h, 2026-01-01)');
ylabel('Z (km)');

% --- 图 5：3D 轨迹 + 历元顺序 ---
f5 = figure('Name','Fig 5: 3D trajectory','Position',[180 100 900 700]);
plot3(xx/1000, yy/1000, zz/1000, 'b-', 'LineWidth', 1.0); hold on; grid on;
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title(sprintf('G%02d  ECEF 三维轨迹（标号 = 历元序号）', prn));
% 每隔约 5% 的点画一个标号（用原始索引，不含 NaN 插入点）
N   = length(X);
step_anno = max(1, floor(N/20));
idx = 1:step_anno:N;
scatter3(X(idx)/1000, Y(idx)/1000, Z(idx)/1000, 30, 'ro', 'filled');
for ii = 1:length(idx)
    text(X(idx(ii))/1000, Y(idx(ii))/1000, Z(idx(ii))/1000, ...
        sprintf(' %d', idx(ii)), 'FontSize', 8, 'Color', 'k');
end
% 加一个箭头标方向（首→末）
quiver3(X(1)/1000, Y(1)/1000, Z(1)/1000, ...
        (X(end)-X(1))/1000*0.05, (Y(end)-Y(1))/1000*0.05, (Z(end)-Z(1))/1000*0.05, ...
        3, 'r', 'LineWidth', 1.5);

% 保存所有图
saveas(f1, 'fig1_L1phase.png');
saveas(f2, 'fig2_doppler.png');
saveas(f3, 'fig3_snr.png');
saveas(f4, 'fig4_xyz_components.png');
saveas(f5, 'fig5_3d_trajectory.png');
fprintf('5 张图已保存为 fig1~fig5 的 PNG。\n');

%% ---------- 第 5 步：写"时间 + 数据"表 ----------
fid = fopen(out_table, 'w');
fprintf(fid, '# G%02d  时间 + 数据表\n', prn);
fprintf(fid, '# 时间系统: GPST (GPS Time)\n');
fprintf(fid, '# 列含义:\n');
fprintf(fid, '#   1. GPST 连续秒 (自 1980-01-06 00:00:00 起的秒数)\n');
fprintf(fid, '#   2. GPST 日期时间 (YYYY-MM-DD HH:MM:SS.sss)\n');
fprintf(fid, '#   3. L1 载波相位 (cycles)\n');
fprintf(fid, '#   4. L1 多普勒 (Hz)\n');
fprintf(fid, '#   5. L1 信噪比 (dB-Hz)\n');
fprintf(fid, '#   6. X (m, ECEF)\n');
fprintf(fid, '#   7. Y (m, ECEF)\n');
fprintf(fid, '#   8. Z (m, ECEF)\n');
fprintf(fid, '# 每行对应一个观测历元 (历元间隔 30 s)\n');
fprintf(fid, '# ------------------------------------------------------------\n');
fclose(fid);

% 用循环逐行写（保证日期字符串不被破坏）
fid = fopen(out_table, 'a');
for k = 1:length(T)
    fprintf(fid, '%14.3f  %s  %14.4f  %12.4f  %8.3f  %14.3f  %14.3f  %14.3f\n', ...
        T(k), tdstr{k}, L1phas(k), Dop(k), SNR(k), X(k), Y(k), Z(k));
end
fclose(fid);
fprintf('数据表已保存到 %s （共 %d 行）。\n', out_table, length(T));

% ===== 局部函数：从 RINEX 3 行的固定字符位提取一个 obs 值 =====
function v = parse_field(lin, p0)
    % lin: 完整观测行
    % p0 : obs 值字段起点（14 个字符）
    if length(lin) < p0 + 13
        v = NaN; return;
    end
    s = lin(p0 : p0+13);
    s = strtrim(s);
    if isempty(s)
        v = NaN;
    else
        v = sscanf(s, '%f');
        if isempty(v), v = NaN; end
    end
end
