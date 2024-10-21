import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import math
from math import pi,sqrt,exp 
from uncertainties import ufloat

def show_results():
    print(f'Pressão Lab: {pressão_lab}')
    print(f'Temperatura Lab: {temperatura_lab}')
    print(f'Temperatura S: {Ts_media}')
    print(f'Tf: {Tf}')
    print(f'Rayleigh: {Ra}')
    print(f'Prandtl: {Pr}')
    print(f'h teórico laminar: {h_conv_teórico_laminar}')
    print(f'h radiação: {h_rad_teórico}')
    print(f'h total teórico: {h_conv_teórico_laminar + h_rad_teórico}')
    print(f'h experimental: {h_experimental}')

def calcular_densidade():
    densidade_mercurio = ufloat(13600,2720) #kg/m^3 -> PADRÃO
    gravidade = ufloat(9.78,0.05) #m/s^2 -> PADRÃO
    h_baro = ufloat(0.699,0.001) #m
    return densidade_mercurio*gravidade*h_baro

#----------------INCERTEZAS--------------#
i_paq = 0.03/1000 #m
#---------------DADOS DO LABORATÓRIO---------------#
Fluido = 'Air'
correção_kelvin = 273.15
temperatura_lab = ufloat(np.mean([27.6,27.8,27.6,27.7]) + correção_kelvin, 0.1) #K -> MEDIR TERMOMETRO
pressão_lab = ufloat(91500,0)
gravidade = ufloat(9.78,0.05) #m/s^2 -> PADRÃO

diametro_externo = ufloat(0.02575,0.005) #m -> MEDIR PAQ
comprimento_da_barra = ufloat(0.195,0.005) #m -> MEDIR PAQ
Ts1 = np.mean([72.7,72.8,72.6,72.6])+correção_kelvin #K -> MEDIR
Ts2 = np.mean([74,73.7,73.8,74])+correção_kelvin
Ts3 = np.mean([71.7,71.6,71.8,71.6])+correção_kelvin
V =  30 #Volts -> MEDIR
I =  0.293 #Amperes -> MEDIR

Ts1_ufloat = ufloat(Ts1, 1) #USANDO O CONJUNTO DE DADOS TS1, FAZER PARA CADA 
Ts2_ufloat = ufloat(Ts2, 1)
Ts3_ufloat = ufloat(Ts3, 1)
Ts_media = ufloat(np.mean([Ts1,Ts2,Ts3]), 1)

V_ufloat = ufloat(V, V*0.03)
I_ufloat = ufloat(I, I*0.03)

emissividade = 0.6
stephan_boltzman = 5.67*10**(-8)
#-------------PROPRIEDADES TÉRMICAS----------------#
Tf = (Ts_media + temperatura_lab)/2
betha=PropsSI("isobaric_expansion_coefficient","T",Tf.nominal_value,"P",pressão_lab.nominal_value,Fluido)      #coeficiente de expansao termica
Pr=PropsSI("Prandtl","T",Tf.nominal_value,"P",pressão_lab.nominal_value,Fluido)                                #Numero de Pr
mu=PropsSI("V","T",Tf.nominal_value,"P",pressão_lab.nominal_value,Fluido)                                      #vicosidade dinamica
rho=PropsSI("D","T",Tf.nominal_value,"P",pressão_lab.nominal_value,Fluido)                                     #densidade
ni=mu/rho                                                               #viscosidade cinematica

"Tipo de escoamento"
Gr=(gravidade*betha*(Ts_media-temperatura_lab)*comprimento_da_barra**3)/(ni**2)         #Numero de Gr
Ra=Gr*Pr                                    #Numero de Ra

#------------h PRÁTICO----------#
h_conv_teórico_laminar = 1.32*((Ts_media.nominal_value - temperatura_lab.nominal_value)/diametro_externo.nominal_value)**(1/4)
h_rad_teórico = emissividade*stephan_boltzman*(Ts_media.nominal_value**2 + temperatura_lab.nominal_value**2)*(Ts_media.nominal_value + temperatura_lab.nominal_value)
#-----------h EXPERIMENTAL--------------#
potência = V_ufloat*I_ufloat
h_experimental = potência/((math.pi*diametro_externo*comprimento_da_barra)*(Ts_media-temperatura_lab))
#-------------RESULTADOS------------#
show_results()