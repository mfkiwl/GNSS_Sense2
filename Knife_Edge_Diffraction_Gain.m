% =========================================================================
% plot_knife_edge_diffraction.m
% 功能: 绘制标准的刀锋衍射增益图 (Knife-Edge Diffraction Gain)
% 理论基础: 菲涅尔积分 (Fresnel Integrals)
% =========================================================================

clear; clc; close all;

% 1. 定义菲涅尔衍射参数 v 的范围
% v < 0: 视距(LoS)畅通，但在第一菲涅尔区内 (有震荡)
% v = 0: 视距刚好被阻挡 (Grazing incidence, -6dB)
% v > 0: 视距被阻挡 (进入阴影区)
v = linspace(-5, 5, 1000); 

% 2. 计算菲涅尔积分 C(v) 和 S(v)
C = fresnelc(v);
S = fresnels(v);

% 3. 计算衍射增益 (Diffraction Gain/Loss)
% 公式推导: F(v) = (1+j)/2 * integral_v^inf exp(-j*pi/2*t^2) dt
% 模值: |F(v)| = sqrt( ((1/2 - C(v)).^2 + (1/2 - S(v)).^2) / 2 )
F_abs = sqrt( ((0.5 - C).^2 + (0.5 - S).^2) / 2 );

% 转换为分贝 (dB)
Gain_dB = 20 * log10(F_abs);

% 4. 绘图
figure('Color', 'w', 'Position', [300, 300, 800, 500]);
plot(v, Gain_dB, 'b-', 'LineWidth', 2); hold on;

% 标注关键点
grid on;
xlabel('Fresnel Diffraction Parameter, \nu', 'FontSize', 12);
ylabel('Diffraction Gain (dB)', 'FontSize', 12);
title('Knife-Edge Diffraction Gain', 'FontSize', 14, 'FontWeight', 'bold');

% 标注 v=0 (-6dB) 点
plot(0, -6.02, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(0.2, -6, 'v=0, Gain \approx -6 dB (LoS Grazing)', 'FontSize', 10);

% 标注区域
xline(0, 'k--', 'LineWidth', 1.5);
text(-3, 2, 'Interference Region (Line-of-Sight)', 'Color', [0 0.5 0], 'FontWeight', 'bold');
text(2, -20, 'Shadow Region (Obstructed)', 'Color', [0.5 0 0], 'FontWeight', 'bold');

% 标注第一菲涅尔区边界 (v ~ -1.414)
% 第一菲涅尔区边界对应路径差 lambda/2，即 v = -sqrt(2) ≈ -1.414
xline(-sqrt(2), 'm:', 'LineWidth', 1.5);
text(-sqrt(2), -15, ' 1st Fresnel Zone Boundary', 'Color', 'm', 'Rotation', 90);

ylim([-30, 6]);
xlim([-5, 5]);