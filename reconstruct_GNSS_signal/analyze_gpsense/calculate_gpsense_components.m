% ============== calculate_gpsense_components.m ==============
function [residual_doppler, pr_phase_cycles] = calculate_gpsense_components(current_epoch, sat_id, sat_states, receiver_pos, receiver_clk_err_sec)
% CALCULATE_GPSENSE_COMPONENTS - 为GPSense相位重建计算所需的组件。
%
% 输出:
%    residual_doppler - 残余多普勒 (测量多普勒 - 几何多普勒)，单位 Hz (周/秒)。
%    pr_phase_cycles  - 基于伪距估计的绝对相位，单位 周。

% --- 常数与数据提取 ---
c = 299792458.0;
f_L1 = 1575.42e6;
lambda_L1 = c / f_L1;

sat_state = sat_states.(sat_id);
P_raw = current_epoch.data.(sat_id).pseudorange.C1C;
measured_doppler = current_epoch.data.(sat_id).doppler.D1C;

if isnan(P_raw) || isnan(measured_doppler)
    residual_doppler = NaN;
    pr_phase_cycles = NaN;
    return;
end

% --- 计算组件 ---
% 1. 计算几何多普勒 f_D
sat_pos = sat_state.position;
sat_vel = sat_state.velocity;
v_rel = sat_vel; % 假设接收机静止
los_vec = (sat_pos - receiver_pos) / norm(sat_pos - receiver_pos);
fD_geometric = -(v_rel' * los_vec) / lambda_L1;

% 2. 计算残余多普勒 (dΦ/dt - f_D)
residual_doppler = measured_doppler - fD_geometric;

% 3. 计算基于伪距的相位项 (ρ - c*tb)/λ
pr_phase_cycles = (P_raw - c * receiver_clk_err_sec) / lambda_L1;

end