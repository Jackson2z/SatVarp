using DifferentialEquations
using Plots
using LinearAlgebra

"""
	2025/5/3 14:11
	1. 完成了计算初始含湿量的函数
	2. 完成了计算冷却塔微分方程组的函数
	3. 完成了求解冷却塔微分方程组的函数

"""


# 物理常数
P_atm = 101325.0  # 大气压力 (Pa)
R = 8.314  # 通用气体常数 (J/mol·K)
M_w = 18.01528e-3  # 水的摩尔质量 (kg/mol)
M_a = 28.97e-3  # 空气的摩尔质量 (kg/mol)
M_s = 58.44e-3  # 盐(NaCl)的摩尔质量 (kg/mol)
h_fg = 2257e3  # 水的汽化潜热 (J/kg)
t_ref = 0.0  # 参考温度 (°C)

# 新疆典型气象参数 (简化版)
t_air_in = 35.0  # 入口空气温度,干球温度 (°C)
RH_in = 0.16  # 入口空气相对湿度

# 冷却塔运行参数
t_sw_in = 35.0  # 入口海水温度 (°C)
S_in = 74.4  # 入口海水盐度 (g/kg)
m_sw = 18  # 海水质量流量 (kg/s)
m_air = 34.625  # 空气质量流量 (kg/s).根据比热容计算
V = 100.0  # 冷却塔体积 (m³)
a = 145.0  # 单位体积传热传质面积 (m²/m³)

# 结晶参数
# S_sat = 74.75 
include("pkg.jl")

# 计算初始含湿量
w_in = calc_air_humidity(t_air_in, RH_in)



function main1(t_air_in, w_in, t_sw_in, S_in)

	sol = solve_cooling_tower(t_air_in, w_in, t_sw_in, S_in)# 求解
	# 1. 绘制结果变化曲线
	plot_results(sol)


	# 从解对象获取出口参数
	h_a_out, w_out, t_sw_out, S_out = sol[end]  # 出口空气焓,出口空气含湿量,出口海水温度,出口海水盐度

	t_air_out=calc_sw_temp(h_a_out, w_out) # 出口空气温度


	# 2. 典型工况计算
	println("\n===== 典型工况计算 =====")

	println("\n=== 冷却塔详细性能计算结果 ===")

	println("\n【空气参数】")
	println("- 入口温度: $(t_air_in) °C")
	println("- 出口温度: $(round(t_air_out, digits=2)) °C")
	println("- 入口湿度比: $(round(w_in*1e3, digits=2)) g/kg")
	println("- 出口湿度比: $(round(w_out*1e3, digits=2)) g/kg")

	println("\n【海水参数】")
	println("- 入口温度: $(t_sw_in) °C")
	println("- 出口温度: $(round(t_sw_out, digits=2)) °C")
	println("- 入口盐度: $(S_in) g/kg")
	println("- 出口盐度: $(round(S_out, digits=2)) g/kg")
	S_sat=calc_sat_solubility(t_sw_out)
	println("- 饱和盐度: $(round(S_sat, digits=2)) g/kg")

end

function main2(t_air_in, w_in, t_sw_in, S_in)

	sol = solve_cooling_tower(t_air_in, w_in, t_sw_in, S_in)# 求解
	# 1. 绘制结果变化曲线
	plot_results(sol)


	# 从解对象获取出口参数
	h_a_out, w_out, t_sw_out, S_out = sol[end]  # 出口空气焓,出口空气含湿量,出口海水温度,出口海水盐度
	t_air_out=calc_sw_temp(h_a_out, w_out) # 出口空气温度

return h_a_out, w_out, t_sw_out, S_out,t_air_out

end


#? 单次打印
main1(t_air_in, w_in, t_sw_in, S_in)

#? 循环运行
lnn=1

h_a_out=zeros(3,3)
w_out=zeros(3,3)
t_sw_out=zeros(3,3)
S_out=zeros(3,3)
t_air_out=zeros(3,3)

for (lnn, t_air_in) in enumerate([30.0, 35.0, 40.0])
    for (lnm, t_sw_in) in enumerate([30.0, 35.0, 40.0])
        h_a_out[lnn,lnm], w_out[lnn,lnm], t_sw_out[lnn,lnm], S_out[lnn,lnm],t_air_out[lnn,lnm]=main2(t_air_in, w_in, t_sw_in, S_in)
    end
end
