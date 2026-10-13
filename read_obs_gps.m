function [epoch_t, prn_cell, pr_cell] = read_obs_gps(filename)
% ============================================================
% 读 RINEX 3 观测文件(.o)，只提取 GPS 卫星的 C1C 伪距
%
% 输出：
%   epoch_t  : 每个历元的 GPS 连续时刻（秒）
%   prn_cell : 1xN cell，每个元素是该历元可见 GPS 卫星号数组
%   pr_cell  : 1xN cell，每个元素是对应的 C1C 伪距(米)
%
% 观测文件每行格式（RINEX 3）：
%   历元首行：> 年 月 日 时 分 秒 标志 卫星数
%   之后每颗卫星一行：G12  伪距值 ...（后面还有载波相位等，用不到）
% ============================================================

fid = fopen(filename, 'rt');
if fid == -1
    error('打不开文件：%s', filename);
end

% 跳过文件头
while true
    lin = fgetl(fid);
    if ~ischar(lin)
        error('文件格式不对');
    end
    if contains(lin, 'END OF HEADER')
        break;
    end
end

epoch_t = []; prn_cell = {}; pr_cell = {};

while ~feof(fid)
    lin = fgetl(fid);
    if ~ischar(lin), break; end
    if isempty(strtrim(lin)), continue; end
    if lin(1) ~= '>', continue;                 % 不是历元首行，跳过
    end

    % 解析历元首行：> 2026 01 01 00 00  0.0000000  0 46
    data = sscanf(lin(2:end), '%f');
    if length(data) < 8, continue; end
    flag = data(7);
    if flag ~= 0, continue; end                 % 异常历元跳过
    nsat = data(8);

    [week, tow] = date2gps(data(1:6)');         % 转 GPS 周、周内秒
    t = week*604800 + tow;                       % 连续秒

    % 逐颗卫星读一行
    prns = []; prs = [];
    for s = 1:nsat
        lin2 = fgetl(fid);
        if length(lin2) < 3, continue; end
        if lin2(1) ~= 'G', continue; end         % 只保留 GPS
        prn = str2double(lin2(2:3));
        v = sscanf(lin2(5:end), '%f', 1);        % 该行第一个数 = C1C 伪距
        if ~isempty(v) && v > 1e6                % 合理伪距 2万~3万公里
            prns(end+1) = prn; %#ok<AGROW>
            prs(end+1)  = v;   %#ok<AGROW>
        end
    end

    if length(prns) >= 5                          % 至少 5 颗才定位
        epoch_t(end+1) = t;
        prn_cell{end+1} = prns;
        pr_cell{end+1}  = prs;
    end
end
fclose(fid);
fprintf('共读取 %d 个有效历元（GPS）\n', length(epoch_t));
end
