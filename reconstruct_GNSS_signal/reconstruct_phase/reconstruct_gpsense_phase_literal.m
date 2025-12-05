% ============== reconstruct_gpsense_phase_literal.m (适配GPS和北斗版) ==============
function recon_phase_rad = reconstruct_gpsense_phase_literal(dPhi_dt_calculated, P_raw, sat_state, receiver_pos, receiver_clk_err_sec, t_elapsed)
% RECONSTRUCT_GPSENSE_PHASE_LITERAL (覆盖修改版)
% 严格按GPSense论文公式(5)的字面意思构建相位，支持GPS和北斗。

% --- 常数 ---
c = 299792458.0;
f_GPS_L1 = 1575.42e6;
f_BDS_B1I = 1561.098e6;
lambda_GPS_L1 = c / f_GPS_L1;
lambda_BDS_B1I = c / f_BDS_B1I;

% 检查输入是否有效
if isnan(dPhi_dt_calculated) || isnan(P_raw)
    recon_phase_rad = NaN;
    return;
end

% --- 根据卫星系统选择正确的波长 ---
% 我们假设 sat_state 结构体中包含一个 'System' 字段
if strcmpi(sat_state.System, 'GPS')
    lambda = lambda_GPS_L1;
elseif strcmpi(sat_state.System, 'BeiDou')
    lambda = lambda_BDS_B1I;
else
    error('未知的卫星系统: %s', sat_state.System);
end

% === 几何多普勒频率 f_D (单位：Hz) ===
sat_pos = sat_state.position;
sat_vel = sat_state.velocity;
los_vec = (sat_pos - receiver_pos) / norm(sat_pos - receiver_pos);
fD_geometric = - (sat_vel' * los_vec) / lambda;

% === 伪距项 (ρ − c·tb) / λ，单位：cycles（周） ===
pseudorange_term_cycles = (P_raw - c * receiver_clk_err_sec) / lambda;

% === 构建总周数 (使用传入的dPhi_dt) ===
total_cycles = dPhi_dt_calculated - fD_geometric * t_elapsed + pseudorange_term_cycles;

% === 转为弧度 ===
recon_phase_rad = 2 * pi * total_cycles;
end