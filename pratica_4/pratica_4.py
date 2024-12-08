import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import math
from math import pi,sqrt,exp 
from uncertainties import ufloat
import backend_p4 as BE
#---------INCERTEZAS--------#
i_temperatura = 1 #°C
i_medidor_de_volume = 1*10**-6 #mL
i_comprimento = 0.03/1000 #mm -> m
#---------MEDIÇÕES do COBRE----------------#
T1 = np.mean([124.4,124.6]) #°C
T2 = np.mean([132.3,132.3]) #°C
T3 = np.mean([138.5,138.6]) #°C
T_amb = np.mean([25.2,25.9]) #°C
P_amb = np.mean([914.9,914.6]) #hPa

T1 = BE.convert_celsius_to_kelvin(T1,i_temperatura) #K
T2 = BE.convert_celsius_to_kelvin(T2,i_temperatura) #K
T3 = BE.convert_celsius_to_kelvin(T3,i_temperatura) #K
T_amb = BE.convert_celsius_to_kelvin(T_amb,i_temperatura) #K
P_amb = BE.convert_hPa_to_Pa(P_amb,0.15*P_amb)
d = ufloat(30/1000,i_comprimento) #mm -> m
x1 = ufloat(3/1000,i_comprimento) #mm -> m
x2 = ufloat(22/1000,i_comprimento) #mm -> m
x3 = ufloat(31/1000,i_comprimento) #mm -> m
#---------MEDIÇÕES CIRCUITO DE AGUA-----------#
T_eb = np.mean([96.8,97.2]) #°C
T_a1 = np.mean([24.5,24.4]) #°C 
T_a2 = np.mean([38.8,38.9]) #°C

V_a = np.mean([1,1]) #L
t_a = ['05:31.93', '05:40.42']
V_cond_1 = 50 #mL
t_cond_1 = '11:17.61' #min:s.ms
V_cond_2 = 20 #mL
t_cond_2 = '04:14.77' #min:s.ms

T_eb = BE.convert_celsius_to_kelvin(T_eb,i_temperatura)
T_a1 = BE.convert_celsius_to_kelvin(T_a1,i_temperatura)
T_a2 = BE.convert_celsius_to_kelvin(T_a2,i_temperatura)

V_a = BE.convert_liters_to_cubic_meters(V_a,i_medidor_de_volume) #m^3
t_a_arg_0 = BE.convert_time_to_seconds(t_a[0]) #s
t_a_arg_1 = BE.convert_time_to_seconds(t_a[1]) #s
t_a = np.mean([t_a_arg_0,t_a_arg_1]) #s
V_cond_1 = BE.convert_milliliters_to_cubic_meters(V_cond_1,i_medidor_de_volume) #m^3
V_cond_2 = BE.convert_milliliters_to_cubic_meters(V_cond_2,i_medidor_de_volume) #m^3
t_cond_1 = BE.convert_time_to_seconds(t_cond_1) #s
t_cond_2 = BE.convert_time_to_seconds(t_cond_2) #s

vazao_agua = V_a/t_a #m^3/s
vazao_cond_1 = V_cond_1/t_cond_1 #m^3/s
vazao_cond_2 = V_cond_2/t_cond_2 #m^3/s
vazao_cond = np.mean([vazao_cond_1,vazao_cond_2]) #m^3/s
#--------PROPRIEDADES FÍSICAS--------#
densidade_agua_fria = CP.PropsSI("D", "T", T_a1.nominal_value, "P", P_amb.nominal_value, 'Water') #kg/m^3
densidade_agua_quente = CP.PropsSI("D", "T", T_eb.nominal_value, "P", P_amb.nominal_value, 'Water') #kg/m^3
C_p_agua = CP.PropsSI('C', 'T', T_amb.nominal_value, 'P', P_amb.nominal_value, 'Water')
h_lv = PropsSI('H', 'T', T_eb.nominal_value, 'Q', 1, 'Water')  #J/kg
#--------CALCULOS PRELIMINARES-------#
Area = (pi*d**2)/4
vazao_massica_agua = densidade_agua_quente*vazao_agua
vazao_massica_cond = densidade_agua_fria*vazao_cond

calor_agua = vazao_massica_agua*C_p_agua*(T_a2-T_a1)
calor_cond = vazao_massica_cond*h_lv
#--------K EXP---------#
dTdx = (T3-T1)/(x3-x1)
K_exp = calor_agua/(Area*(dTdx + 32))
#--------RESULTADOS---------#
resultados=[['Densidade agua fria:',densidade_agua_fria,'kg/m^3'],
            ['Densidade agua quente:',densidade_agua_quente,'kg/m^3'],
            ['Cp da agua:',C_p_agua,'J/kg°C'],
            ['h da água quente:',h_lv,'J/kg'],
            ['Area:',Area,'m^2'],
            ['Vazão agua:',vazao_agua,'m^3/s'],
            ['Vazão cond:',vazao_cond,'m^3/s'],
            ['Vazão mássica agua:',vazao_massica_agua,'kg/s'],
            ['Vazão massica cond:',vazao_massica_cond,'kg/s'],
            ['Calor da água:',calor_agua,'W'],
            ['Calor de cond:',calor_cond,'W'],
            ['dT/dx:',dTdx,'K/m'],
            ['K_exp:',K_exp,'W/mK'],
            ['Temperatura de ebulição:',T_eb,'K']]

BE.print_results(resultados)