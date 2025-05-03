# ========== 氯化锂溶解度参数 ==========                                                      添加

# ========== 辅助函数 ==========                                                               添加
function calc_sat_solubility(t_li::Float64)
    """计算氯化锂溶解度 (g/100g水)"""
    LiCl_a = 63.7
    LiCl_b = 0.356
    LiCl_c = 0.0012
    LiCl_d = 1.5e-6
    S_sat = LiCl_a + LiCl_b * t_li - LiCl_c * t_li^2 + LiCl_d * t_li^3
    return S_sat
end




#! 计算湿空气热力学参数
function calc_air_sp(t_C::Float64)
    """
    计算0~200℃范围内饱和水蒸气压力(Pa)
    """

    T = t_C + 273.15
    C8 = -5.8002206e3
    C9 = 1.3914993
    C10 = -4.8640239e-2
    C11 = 4.1764768e-5
    C12 = -1.4452093e-8
    C13 = 6.5459673

    ln_ps = (C8 / T) + C9 + (C10 * T) + (C11 * T^2) + (C12 * T^3) + (C13 * log(T))
    p_s = exp(ln_ps)
    return p_s
end

# t=30.0;RH=0.6695
function calc_air_humidity(t::Float64, RH::Float64)
    p_s = calc_air_sp(t)
    p_v = RH * p_s
    p_a = P_atm - p_v
    w_a = 0.622 * p_v / p_a
    return w_a
end

function calc_air_enthalpy(t::Float64, d::Float64)
    h_a = (1.006 * t + d * (2501 + 1.86 * t))
    return h_a
end







#! 海水热力学参数计算
function calc_w_p(t::Float64)
    """计算纯水的蒸汽压 (Pa) - Antoine方程"""
    pv_w = 610.78 * exp(17.2694 * t / (t + 238.3))
    return pv_w
end

function calc_p_ratio(S::Float64)
    """
    	计算 纯水的蒸汽压 / 海水的蒸汽压 的值
    	S (float): 盐浓度

    """
    ratio = 1 + 0.57357 * (S / (1000 - S))
    return ratio
end

function calc_sw_p(t::Float64, S::Float64)
    """计算海水蒸汽压 (Pa)"""
    pv_w = calc_w_p(t)
    ratio = calc_p_ratio(S)
    pv_sw = pv_w / ratio
    return pv_sw
end

function calc_w_enthalpy(t::Float64)
    """计算纯水的焓(kJ/kg)"""
    p = 610.78 * exp(17.2694 * t / (t + 238.3)) / 1e6

    h_w = 78 * log(100 * p) + 400 * √(p) - 200 * p^(1 / 3) - 15 * p + 200

    return h_w
end



#!  海水温度下,壁面饱和湿空气的湿度比和焓。
function calc_sw_humidity(t::Float64, S::Float64) #? 壁面饱和湿空气的含湿量
    """
    计算含湿量d
    pv_sw 计算海水饱和蒸汽压
    p_a干空气的分压力，即空气的压力减去水蒸气的压力
    """
    p_s = calc_sw_p(t, S)
    p_a = P_atm - p_s
    w_s_w = 0.622 * p_s / p_a
    return w_s_w
end

function calc_sw_enthalpy(t::Float64, d::Float64) #?壁面饱和湿空气的焓。
    h_a = (1.006 * t + d * (2501 + 1.86 * t))
    return h_a
end


function calc_sw_temp(h_a::Float64, d::Float64)
    """根据焓值和含湿量反求温度
    参数:
        h_a: 湿空气焓值 (kJ/kg)
        d: 含湿量 (kg/kg)
    返回:
        温度 t (°C)
    """
    t = (h_a - 2501 * d) / (1.006 + 1.86 * d)
    return t
end

#? 海水密度 rho_sw
function calc_water_density(t::Float64)
    """计算纯水密度(kg/m³)"""
    # 纯水密度系数
    a1 = 9.999e2
    a2 = 2.034e-2
    a3 = -6.162e-3
    a4 = 2.261e-5
    a5 = -4.657e-8

    # 计算纯水密度
    rho_w = (a1 + a2 * t + a3 * t^2 + a4 * t^3 + a5 * t^4)
    return rho_w
end

function cal_sw_density(t::Float64, S::Float64)
    """
    计算海水密度 ρ_sw (kg/m³) 公式(3)
    
    参数
    t (float) 温度，单位 °C
    S (float) 盐度，单位 g/kg
    
    返回
    float 海水密度值 (kg/m³)
    """
    # 海水密度修正系数
    b1 = 0.8020
    b2 = -2.001e-3
    b3 = 1.677e-5
    b4 = -3.060e-8
    b5 = -1.613e-11

    # 计算纯水密度
    rho_w = calc_w_density(t)

    # 计算海水密度修正项
    salt_correction = S * (b1 + b2 * t + b3 * t^2 + b4 * t^3 + b5 * S * t^2)

    # 计算海水密度x
    rho_sw = rho_w + salt_correction

    return rho_sw
end

#? 海水比热容
function calc_sw_cp(t::Float64, S::Float64)#? 已检查+1
    """K
    	计算海水比热容 cp_sw (KJ/kg·K) 的值
    	
    	参数
    	T (float) 温度，单位 K (范围 273-353K，即0-80°C)
    	S (float) 盐度，单位 g/kg (范围建议 0-40g/kg)

    """
    T = t + 273.15  # 转换为摄氏温度用于计算 - 正确定义变量t

    # 公式2
    A = 5.328 - 9.76E-2 * S + 4.04E-4 * S^2
    B = -6.913E-3 + 7.351E-4 * S - 3.15E-6 * S^2
    C = 9.6E-6 - 1.927E-6 * S + 8.23E-9 * S^2
    D = 2.5E-9 + 1.666E-9 * S - 7.125E-12 * S^2
    cp_sw = A + B * T + C * T^2 + D * T^3

    return cp_sw

end

#! 计算中间无量纲参数
function lewis_factor(t_sw, t_air, S)   #?没问题+1
    """计算刘易斯因子
    ω_s_w​：壁面处空气的饱和含湿量（由壁面温度决定）；
    ω_s_a​：主流空气的实际含湿量；
    
    """
    w_s_w = calc_sw_humidity(t_sw, S)
    w_s_air = calc_air_sp(t_air)

    term = ((w_s_w + 0.622) / (w_s_air + 0.622) - 1) / log((w_s_w + 0.622) / (w_s_air + 0.622))
    Le = 0.865^0.667 * term
    return Le
end
function mass_transfer_coefficient(t_sw, S)  #? 已检查+1
    """计算传质系数 (kg/m²·s)"""
    h_D0 = 0.01  # 淡水基准传质系数
    h_D = h_D0 * (1 - 0.05 * S / 120)
    return h_D
end
function merkel_number(t_sw, S) #? 已检查
    """计算默克尔数"""
    h_D = mass_transfer_coefficient(t_sw, S)
    Me = h_D * a * V / m_sw
    return Me
end
function mass_flow_ratio()  #? 已检查
    """计算质量流量比"""
    MR = m_sw / m_air
    return MR
end




#!======= 冷却塔微分方程组 =======

# 冷却塔微分方程组
function cooling_tower_equations(du, u, p, z)
	h_a, w, t_sw, S = u

	#MR = mass_flow_ratio()
	#Me = merkel_number(t_sw, S)
	MR = 0.51986   # 质量流量比                                                            修改 
	Me = 1.78896   # Merkel数
	Le = lewis_factor(t_sw, t_air_in, S)

	w_s_w = calc_sw_humidity(t_sw, S) #计算壁面处空气的饱和含湿量
	h_s_w = calc_sw_enthalpy(t_sw, w_s_w)    # 计算壁面处空气的饱和焓

	h_v = calc_air_enthalpy(t_sw, w)
	cp_sw = calc_sw_cp(t_sw, S)

	w_sw_change = w_s_w - w

	NTU = MR * Me

	du[1] = NTU * (Le * (h_s_w - h_a) + (1 - Le) * w_sw_change * h_v)
	du[2] = NTU * w_sw_change
	w_o=du[2]
	w_change = w_o - w

	du[3] = -(1 / (MR - w_change)) * ((1 / cp_sw) * du[1] - (t_sw - t_ref) * du[2])
	S_sat = calc_sat_solubility(du[3]) #对应浓度的食盐饱和浓度
	du[4] = (min(S, S_sat) / (MR - w_change)) * du[2]
end

function solve_cooling_tower(t_air_in, w_in, t_sw_in, S_in)
	"""求解冷却塔模型"""
	#? 初始化
	h_a_in = calc_air_enthalpy(t_air_in, w_in)
	u0 = [h_a_in, w_in, t_sw_in, S_in] # 初始值 

	tspan = (0.0, 15.0) # 时间范围
	prob = ODEProblem(cooling_tower_equations, u0, tspan)
	sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
	return sol
end

function plot_results(sol)
	# 创建4个子图
	p1 = plot(sol[1, :], sol[4, :], label = "Air Enthalpy (kJ/kg)", xlabel = "TimeStep", ylabel = "盐结晶量")
	p2 = plot(sol[2, :] .* 1000, sol[4, :], label = "Air Relative humidity (g/kg)", xlabel = "TimeStep", ylabel = "盐结晶量")
	p3 = plot(sol[3, :], sol[4, :], label = "Sea water temperature (°C)", xlabel = "TimeStep", ylabel = "盐结晶量")
	p4 = plot(sol[1, :], sol[4, :], label = "Salinity of seawater (g/kg)", xlabel = "TimeStep", ylabel = "盐结晶量")

	# 组合所有子图
	plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 800), legend = :bottomright)

	# 保存图像
	# savefig("cooling_tower_results.png")
	# println("\n结果图已保存为 cooling_tower_results.png")
end
function cooling_tower_equations(du, u, p, z)
	h_a, w, t_sw, S = u

	#MR = mass_flow_ratio()
	#Me = merkel_number(t_sw, S)
	MR = 0.51986   # 质量流量比                                                            修改 
	Me = 1.78896   # Merkel数
	Le = lewis_factor(t_sw, t_air_in, S)

	w_s_w = calc_sw_humidity(t_sw, S) #计算壁面处空气的饱和含湿量
	h_s_w = calc_sw_enthalpy(t_sw, w_s_w)    # 计算壁面处空气的饱和焓

	h_v = calc_air_enthalpy(t_sw, w)
	cp_sw = calc_sw_cp(t_sw, S)

	w_sw_change = w_s_w - w

	NTU = MR * Me

	du[1] = NTU * (Le * (h_s_w - h_a) + (1 - Le) * w_sw_change * h_v)
	du[2] = NTU * w_sw_change
	w_o=du[2]
	w_change = w_o - w

	du[3] = -(1 / (MR - w_change)) * ((1 / cp_sw) * du[1] - (t_sw - t_ref) * du[2])
	S_sat = calc_sat_solubility(du[3]) #对应浓度的食盐饱和浓度
	du[4] = (min(S, S_sat) / (MR - w_change)) * du[2]
end

function solve_cooling_tower(t_air_in, w_in, t_sw_in, S_in)
	"""求解冷却塔模型"""
	#? 初始化
	h_a_in = calc_air_enthalpy(t_air_in, w_in)
	u0 = [h_a_in, w_in, t_sw_in, S_in] # 初始值 

	tspan = (0.0, 15.0) # 时间范围
	prob = ODEProblem(cooling_tower_equations, u0, tspan)
	sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
	return sol
end

function plot_results(sol)
	# 创建4个子图
	p1 = plot(sol[1, :], sol[4, :], label = "Air Enthalpy (kJ/kg)", xlabel = "TimeStep", ylabel = "盐结晶量")
	p2 = plot(sol[2, :] .* 1000, sol[4, :], label = "Air Relative humidity (g/kg)", xlabel = "TimeStep", ylabel = "盐结晶量")
	p3 = plot(sol[3, :], sol[4, :], label = "Sea water temperature (°C)", xlabel = "TimeStep", ylabel = "盐结晶量")
	p4 = plot(sol[1, :], sol[4, :], label = "Salinity of seawater (g/kg)", xlabel = "TimeStep", ylabel = "盐结晶量")

	# 组合所有子图
	plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 800), legend = :bottomright)

	# 保存图像
	# savefig("cooling_tower_results.png")
	# println("\n结果图已保存为 cooling_tower_results.png")
end