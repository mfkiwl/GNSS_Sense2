% ============== reconstruct_gps_signals.m (适配GPS和北斗版) ==============
function [recon_pseudorange, recon_phase] = reconstruct_gps_signals(method, current_epoch, sat_id, ...
                                                                    sat_states, receiver_pos, receiver_clk_err_sec)
% RECONSTRUCT_GPS_SIGNALS - 根据严谨的物理模型重建伪距和相位，支持GPS和北斗。

% --- 常数定义 ---
c = 299792458.0;
f_GPS_L1 = 1575.42e6;
f_BDS_B1I = 1561.098e6;
lambda_GPS_L1 = c / f_GPS_L1;
lambda_BDS_B1I = c / f_BDS_B1I;

% --- 根据卫星系统选择正确的观测码和波长 ---
sys_char = upper(sat_id(1));
if sys_char == 'G'
    pr_code = 'C1C';
    ph_code = 'L1C';
    lambda = lambda_GPS_L1;
elseif sys_char == 'C'
    % RINEX 3.x 中 B1I 信号的标准码
    pr_code = 'C2I'; 
    ph_code = 'L2I';
    lambda = lambda_BDS_B1I;
else
    error('不支持的卫星系统: %s', sat_id);
end

% --- 数据提取 ---
if ~(isfield(current_epoch.data, sat_id) && ...
     isfield(current_epoch.data.(sat_id), 'pseudorange') && isfield(current_epoch.data.(sat_id).pseudorange, pr_code) && ...
     isfield(current_epoch.data.(sat_id), 'carrier_phase') && isfield(current_epoch.data.(sat_id).carrier_phase, ph_code))
    recon_pseudorange = NaN; recon_phase = NaN;
    warning('卫星 %s 缺少 %s 或 %s 观测值，无法重建。', sat_id, pr_code, ph_code);
    return;
end

P_raw = current_epoch.data.(sat_id).pseudorange.(pr_code);
L_raw_cycles = current_epoch.data.(sat_id).carrier_phase.(ph_code);
sat_state = sat_states.(sat_id);
sat_clk_err = sat_state.clock_error;

if isnan(P_raw) || isnan(L_raw_cycles)
    recon_pseudorange = NaN; recon_phase = NaN;
    warning('卫星 %s 的 %s 或 %s 观测值为NaN。', sat_id, pr_code, ph_code);
    return;
end

% --- 1. 重建伪距 (通用) ---
recon_pseudorange = P_raw - c * (receiver_clk_err_sec - sat_clk_err);

% --- 2. 根据方法重建相位 ---
L_raw_meters = L_raw_cycles * lambda;
recon_phase_std_meters = L_raw_meters + c * receiver_clk_err_sec - c * sat_clk_err;

if strcmpi(method, 'standard')
    recon_phase = recon_phase_std_meters;
elseif strcmpi(method, 'gpsense')
    sat_pos = sat_state.position;
    geom_range = norm(sat_pos - receiver_pos);
    residual_phase_meters = recon_phase_std_meters - geom_range;
    recon_phase = (residual_phase_meters / lambda) * (2 * pi);
else
    error("未知的方法，请选择 'standard' 或 'gpsense'。");
end

end