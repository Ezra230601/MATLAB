function [Ek, n] = ecc_anomaly(t, Eph)
% ============================================================
% 解偏近点角 Ek（开普勒方程迭代）
% 这是 satellite_orbits.m 内部调用的小步骤：
%   给定一组星历 Eph 和时刻 t，先算平均角速度 n，
%   再解 M = E - e*sin(E) 求出 Ek。
% ============================================================

roota  = Eph(4);     % 半长轴的平方根
ecc    = Eph(6);     % 偏心率
deltan = Eph(5);     % 平均角速度修正
M0     = Eph(3);     % 平近点角初值
time_eph = Eph(32);  % 星历对应的连续时刻

mu = 3.986005e14;    % 地球引力常数 m^3/s^2

A  = roota^2;        % 半长轴
n0 = sqrt(mu / A^3); % 理想平均角速度
n  = n0 + deltan;    % 修正后的平均角速度

tk = check_t(t - time_eph);   % 距星历历元的时间差
Mk = M0 + n * tk;             % 当前平近点角

Ek = Mk;                      % 迭代解开普勒方程
for k = 1:10
    Ek = Mk + ecc * sin(Ek);
end

end
