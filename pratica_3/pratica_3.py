import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import math
from math import pi,sqrt,exp 
from uncertainties import ufloat
import backend_p3 as BE
#-----------------MEDIÇÕES------------------#
comprimento_da_barra = 298/1000 #mm -> m
diametro_da_barra = 26.3/1000 #mm -> m
temperatura_ambiente = np.mean([23.4,24.5]) #°C
t_0 = np.mean([93.5,93.4,93.5]) #°C
t_45 = np.mean([78.8,78.8,79.2]) #°C
t_95 = np.mean([66.8,67.2,67.7]) #°C
t_145 = np.mean([59.1,59.5,59.6]) #°C
t_195 = np.mean([53.9,54,54.2]) #°C
t_245 = np.mean([50.5,50.6,50.7]) #°C
t_295 = np.mean([48.6,49.1,48.9]) #°C

#----------------CONSTANTES FORNECIDAS----------#
K_aço = 60 #W/m°C
emissividade_aleta = 0.6 

#---------------INCERTEZA NOS DADOS--------------#
i_temp = 1 #°C
i_comprimento = 0.1/1000 #mm -> m

temperatura_ambiente = BE.convert_celsius_to_kelvin(temperatura_ambiente,i_temp)
t_0 = BE.convert_celsius_to_kelvin(t_0,i_temp)
t_45 = BE.convert_celsius_to_kelvin(t_45,i_temp)
t_95 = BE.convert_celsius_to_kelvin(t_95,i_temp)
t_145 = BE.convert_celsius_to_kelvin(t_145,i_temp)
t_195 = BE.convert_celsius_to_kelvin(t_195,i_temp)
t_245 = BE.convert_celsius_to_kelvin(t_245,i_temp)
t_295 = BE.convert_celsius_to_kelvin(t_295,i_temp)
comprimento_da_barra = ufloat(comprimento_da_barra,i_comprimento)
diametro_da_barra = ufloat(diametro_da_barra,i_comprimento)

q_aleta, h = BE.iteration_for_q(12,comprimento_da_barra,diametro_da_barra,t_0,temperatura_ambiente,K_aço,show=False)

#----------EFETIVIDADE e EFICIÊNCIA-------------#
P = pi * diametro_da_barra  # Perímetro da aleta (m)
A_b = (pi * diametro_da_barra**2)/4 
A_aleta = P*comprimento_da_barra
efetividade = q_aleta/(h*A_b*(t_0 - temperatura_ambiente))
eficiencia = q_aleta/(h*A_aleta*(t_0 - temperatura_ambiente))

#-----------GRAFICOS--------------#
x_medições = [0,45,95,145,195,245,295]
t_medições_resultados = [t_0,t_45,t_95,t_145,t_195,t_245,t_295]
t_medições = [t_0.nominal_value,t_45.nominal_value,t_95.nominal_value,
                t_145.nominal_value,t_195.nominal_value,t_245.nominal_value,t_295.nominal_value]
x_array = []
i = 0

while i <= 298:
    x_array.append(i)
    i += 2

#Temperaturas usando a formula para os 150 pontos do x_array
T_x_resultados = [BE.temperature_distribution_on_bar(h,comprimento_da_barra,diametro_da_barra,t_0,
        temperatura_ambiente,K_aço,a/1000) for a in x_array]
T_x = [BE.temperature_distribution_on_bar(h,comprimento_da_barra,diametro_da_barra,t_0,
        temperatura_ambiente,K_aço,a/1000).nominal_value for a in x_array]

x_graf = [x_medições,x_array]
y_graf = [t_medições,T_x]
y_labels = ['Temperaturas obtidas experimentalmente','Temperaturas obtidas teoricamente(EQ 3.75)']
colors = ['red','blue']
BE.t_distribuitions(x_graf,y_graf,'Distribuição de temperatura ao longo da barra','X(mm)',
                    'Temperatura(K)',y_labels,colors)

#Temperaturas usando a formula para os 7 pontos de x_medições
T_x_comparation = [BE.temperature_distribution_on_bar(h,comprimento_da_barra,diametro_da_barra,t_0,
        temperatura_ambiente,K_aço,a/1000).nominal_value for a in x_medições]
#Temperatura 
t_diffs = [np.absolute(a-b) for a,b in zip(T_x_comparation,t_medições)]

BE.t_diffs(x_medições,t_diffs,'Diferenças entre as temperaturas obtidas experimentalmente e teoricamente','X(mm)',
        'ΔT(K)','green')

resultados=[['Pontos da barra:',x_medições,'mm'],
            ['Temperaturas obtidas experimentamente:',t_medições_resultados,'K'],
            ['Calor da aleta:',q_aleta,'W'],
            ['Coeficiente de transferência de calor:',h,''],
            ['Efetividade:',efetividade,''],
            ['Eficiência:',eficiencia,''],
            ['X para a função de destribuição:',x_array,'mm'],
            ['T(x) para os pontos X:',T_x_resultados,'K'],
            ['ΔT:',t_diffs,'K']]

BE.print_results(resultados)