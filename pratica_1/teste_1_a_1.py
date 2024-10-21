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
    print(f'Prandtl: {Pr}')
    print(f'---------------PONTO S1------------------')
    print(f'Ts1: {Ts1_ufloat}')
    print(f'Rayleigh S1: {Ra1}')
    print(f'h teórico laminar S1: {h_conv_teórico_laminar1}')
    print(f'h radiação S1: {h_rad_teórico1}')
    print(f'h total teórico S1: {h_conv_teórico_laminar1 + h_rad_teórico1}')
    print(f'h experimental S1: {h_experimental1.nominal_value}\n')
    print(f'---------------PONTO S2------------------')
    print(f'Ts2: {Ts2_ufloat}')
    print(f'Rayleigh S2: {Ra2}')
    print(f'h teórico laminar S2: {h_conv_teórico_laminar2}')
    print(f'h radiação S2: {h_rad_teórico2}')
    print(f'h total teórico S2: {h_conv_teórico_laminar2 + h_rad_teórico2}')
    print(f'h experimental S2: {h_experimental2.nominal_value}\n')
    print(f'---------------PONTO S3------------------')
    print(f'Ts3: {Ts3_ufloat}')
    print(f'Rayleigh S3: {Ra3}')
    print(f'h teórico laminar S3: {h_conv_teórico_laminar3}')
    print(f'h radiação S3: {h_rad_teórico3}')
    print(f'h total teórico S3: {h_conv_teórico_laminar3 + h_rad_teórico3}')
    print(f'h experimental S3: {h_experimental3.nominal_value}\n')

def calcular_densidade():
    densidade_mercurio = ufloat(13600,2720) #kg/m^3 -> PADRÃO
    gravidade = ufloat(9.78,0.05) #m/s^2 -> PADRÃO
    h_baro = ufloat(0.699,0.001) #m
    return densidade_mercurio*gravidade*h_baro

def plot_curve(x,y,title,x_lable,y_lable):
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='o', linestyle='--', color='b', label=f'{y_lable} vs {x_lable}')
    # Adicionando títulos e rótulos
    plt.title(title)
    plt.xlabel(x_lable)
    plt.ylabel(y_lable)
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_series_of_data(x,y,title,x_lable,y_lable,y_lable_graf,colors):
    plt.figure(figsize=(10, 6))
    for a,b,c in zip(y,y_lable_graf,colors):
        plt.plot(x, a, marker='o', linestyle='-', color=f'{c}', label=f'{b} vs {x_lable}')
    # Adicionando títulos e rótulos
    plt.title(title)
    plt.xlabel(x_lable)
    plt.ylabel(y_lable)
    plt.legend()
    plt.grid(True)
    plt.show()
#----------------INCERTEZAS--------------#
i_paq = 0.03/1000 #m
#---------------DADOS DO LABORATÓRIO---------------#
Fluido = 'Air'
correção_kelvin = 273.15
temperatura_lab = ufloat(np.mean([27.6,27.8,27.6,27.7]) + correção_kelvin, 1) #K -> MEDIR TERMOMETRO
pressão_lab = ufloat(91500,0)
gravidade = ufloat(9.78,0.05) #m/s^2 -> PADRÃO

diametro_externo = ufloat(0.02,i_paq) #m -> MEDIR PAQ
comprimento_da_barra = ufloat(0.195,i_paq) #m -> MEDIR PAQ
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
Gr1=(gravidade*betha*(Ts1_ufloat-temperatura_lab)*comprimento_da_barra**3)/(ni**2)         #Numero de Gr
Ra1=Gr1*Pr                                    #Numero de Ra

Gr2=(gravidade*betha*(Ts2_ufloat-temperatura_lab)*comprimento_da_barra**3)/(ni**2)         #Numero de Gr
Ra2=Gr2*Pr                                    #Numero de Ra

Gr3=(gravidade*betha*(Ts3_ufloat-temperatura_lab)*comprimento_da_barra**3)/(ni**2)         #Numero de Gr
Ra3=Gr3*Pr                                    #Numero de Ra
#------------h PRÁTICO----------#
h_conv_teórico_laminar1 = 1.32*((Ts1_ufloat.nominal_value - temperatura_lab.nominal_value)/diametro_externo.nominal_value)**(1/4)
h_rad_teórico1 = emissividade*stephan_boltzman*(Ts1_ufloat.nominal_value**2 + temperatura_lab.nominal_value**2)*(Ts1_ufloat.nominal_value + temperatura_lab.nominal_value)

h_conv_teórico_laminar2 = 1.32*((Ts2_ufloat.nominal_value - temperatura_lab.nominal_value)/diametro_externo.nominal_value)**(1/4)
h_rad_teórico2 = emissividade*stephan_boltzman*(Ts2_ufloat.nominal_value**2 + temperatura_lab.nominal_value**2)*(Ts2_ufloat.nominal_value + temperatura_lab.nominal_value)

h_conv_teórico_laminar3 = 1.32*((Ts3_ufloat.nominal_value - temperatura_lab.nominal_value)/diametro_externo.nominal_value)**(1/4)
h_rad_teórico3 = emissividade*stephan_boltzman*(Ts3_ufloat.nominal_value**2 + temperatura_lab.nominal_value**2)*(Ts3_ufloat.nominal_value + temperatura_lab.nominal_value)
#-----------h EXPERIMENTAL--------------#
potência = V_ufloat*I_ufloat
h_experimental1 = potência/((math.pi*diametro_externo*comprimento_da_barra)*(Ts1_ufloat-temperatura_lab))
h_experimental2 = potência/((math.pi*diametro_externo*comprimento_da_barra)*(Ts2_ufloat-temperatura_lab))
h_experimental3 = potência/((math.pi*diametro_externo*comprimento_da_barra)*(Ts3_ufloat-temperatura_lab))
#-------------RESULTADOS------------#
show_results()
Ts_graf = [Ts1_ufloat.nominal_value,Ts2_ufloat.nominal_value,Ts3_ufloat.nominal_value]
pontos = [1,2,3]
plot_curve(pontos,Ts_graf,'Gráfico T em cada ponto', 'Ponto', 'Temperatura(K)')

h_conv_t_graf = [h_conv_teórico_laminar1,h_conv_teórico_laminar2,h_conv_teórico_laminar3]
h_rad_t_graf = [h_rad_teórico1,h_rad_teórico2,h_rad_teórico3]
h_teorico_graf = [a + b for a,b in zip(h_conv_t_graf,h_rad_t_graf)]
h_exp_graf = [h_experimental1.nominal_value,h_experimental2.nominal_value,h_experimental3.nominal_value]
y = [h_conv_t_graf,h_rad_t_graf,h_teorico_graf,h_exp_graf]
y_lable_graf = ['h_conv_teórico','h_rad_teórico','h_combinado_teórico','h_combinado_experimental']
colors = ['r','g','b','y']
plot_series_of_data(pontos,y,'h vs Pontos da barra','Ponto da barra','Coeficientes de transferência de calor(h[W/m^2*K])',y_lable_graf,colors)
