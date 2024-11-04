import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import math
from math import pi,sqrt,exp 
from uncertainties import ufloat
import backend_p2 as BP

fluido = 'R12'

#AMBIENTE
temperatura_ambiente = 22.6 #°C
pressão_ambiente = 920.3 #hPa
#22/10 - 20:10 - Horário pro INMET

t_ambiente_kelvin = BP.convert_celsius_to_kelvin(temperatura_ambiente,1)
p_ambiente = ufloat(pressão_ambiente*100,0)

h_fluido_ambiente = BP.get_entalpy(p_ambiente,t_ambiente_kelvin,fluido) / 1000

#COMPRESSOR
pressao_compressor_entrada = 30.5 #psi
pressao_compressor_saida = 491 #psi
temperatura_compressor_entrada = 22 #°C
temperatura_compressor_saida = 85.7 #°C

p_entrada_compressor = BP.convert_psi_to_pascal(pressao_compressor_entrada,1) + p_ambiente
p_saida_compressor = BP.convert_psi_to_pascal(pressao_compressor_saida,1) + p_ambiente
t_entrada_compressor = BP.convert_celsius_to_kelvin(temperatura_compressor_entrada, 1)
t_saida_compressor = BP.convert_celsius_to_kelvin(temperatura_compressor_saida, 1)

h_entrada_compressor = BP.get_entalpy(p_entrada_compressor,t_entrada_compressor,fluido) /1000
h_saida_compressor = BP.get_entalpy(p_saida_compressor,t_saida_compressor,fluido) /1000

#EVAPORADOR
temperatura_de_entrada_evaporador = 2.1 # °C
temperatura_de_saida_evaporador = 8.5 # °C
pressao_entrada_evaporador = 32 #psi
pressao_saida_evaporador = 34 #psi

p_entrada_evaporador = BP.convert_psi_to_pascal(pressao_entrada_evaporador,1) + p_ambiente
p_saida_evaporador = BP.convert_psi_to_pascal(pressao_saida_evaporador,1) + p_ambiente
t_entrada_evaporador = BP.convert_celsius_to_kelvin(temperatura_de_entrada_evaporador, 1)
t_saida_evaporador = BP.convert_celsius_to_kelvin(temperatura_de_saida_evaporador, 1)

h_entrada_evaporador = BP.get_entalpy(p_entrada_evaporador,t_entrada_evaporador,fluido) /1000
h_saida_evaporador = BP.get_entalpy(p_saida_evaporador,t_saida_evaporador,fluido) /1000


#CONDENSADOR
temperatura_entrada_condensador = 71.1 #°C
pressao_entrada_condensador = 140 #psi
temperatura_saida_condensador = 42.7 #°C
pressao_saida_condensador = 140 #psi

p_entrada_condensador = BP.convert_psi_to_pascal(pressao_entrada_condensador,1) + p_ambiente
p_saida_condensador = BP.convert_psi_to_pascal(pressao_saida_condensador,1) + p_ambiente
t_entrada_condensador = BP.convert_celsius_to_kelvin(temperatura_entrada_condensador, 1)
t_saida_condensador = BP.convert_celsius_to_kelvin(temperatura_saida_condensador, 1)

h_entrada_condensador = BP.get_entalpy(p_entrada_condensador,t_entrada_condensador,fluido) /1000
h_saida_condensador = BP.get_entalpy(p_saida_condensador,t_saida_condensador,fluido) /1000

#VALVULA TERMOESTATICA (DE)
temperatura_entrada_valvula = 35.2 #°C
pressao_entrada_valvula = 130 #psi
temperatura_saida_valvula = 2.1 #°C
pressao_saida_valvula = 36 #psi
vazao = 1.3 #lb/min

p_entrada_valvula = BP.convert_psi_to_pascal(pressao_entrada_valvula,1) + p_ambiente
p_saida_valvula = BP.convert_psi_to_pascal(pressao_saida_valvula,1) + p_ambiente
t_entrada_valvula = BP.convert_celsius_to_kelvin(temperatura_entrada_valvula, 1)
t_saida_valvula = BP.convert_celsius_to_kelvin(temperatura_saida_valvula, 1)
vazao_valvula = BP.lb_min_to_kg_s(vazao)

h_entrada_valvula = BP.get_entalpy(p_entrada_valvula,t_entrada_valvula,fluido) /1000
h_saida_valvula = BP.get_entalpy(p_saida_valvula,t_saida_valvula,fluido) /1000

#TROCADOR DE CALOR
'''
A queda de pressão do trocador de calor é negligenciavel.
E por isso, a pressão na entrada deve ser igual a saida, sendo elas a pressão da saida do passo anterior
'''
temperatura_de_entrada_TC_liq = 41.6 #°C
temperatura_de_saida_TC_liq = 34.9 #°C
temperatura_de_entrada_TC_vap = 10.4 #°C
temperatura_de_saida_TC_vap = 16.7 #°C

t_entrada_TC_liq = BP.convert_celsius_to_kelvin(temperatura_de_entrada_TC_liq,1)
t_saida_TC_liq = BP.convert_celsius_to_kelvin(temperatura_de_saida_TC_liq,1)
t_entrada_TC_vap = BP.convert_celsius_to_kelvin(temperatura_de_entrada_TC_vap,1)
t_saida_TC_vap = BP.convert_celsius_to_kelvin(temperatura_de_saida_TC_vap,1)
p_entrada_TC_liq = p_saida_condensador 
p_saida_TC_liq = p_entrada_valvula
p_entrada_TC_vap = p_saida_evaporador
p_saida_TC_vap = p_entrada_compressor

h_entrada_TC_liq = BP.get_entalpy(p_entrada_TC_liq,t_entrada_TC_liq,fluido) /1000
h_saida_TC_liq = BP.get_entalpy(p_saida_TC_liq,t_saida_TC_liq,fluido) /1000
h_entrada_TC_vap = BP.get_entalpy(p_entrada_TC_vap,t_entrada_TC_vap,fluido) /1000
h_saida_TC_vap = BP.get_entalpy(p_saida_TC_vap,t_saida_TC_vap,fluido) /1000
#-----------RESULTADOS-------------#
resultados = [
        ['Entalpia do fluido ambiente:',h_fluido_ambiente,'kJ/kg'],
        ['Pressão de entrada no condensador:',p_entrada_condensador/1000,'kPa'],
        ['Pressão de saida no condensador:',p_saida_condensador/1000,'kPa'],
        ['Entalpia de entrada no compressor:',h_entrada_compressor,'kJ/kg'],
        ['Entalpia de saida no compressor:',h_saida_compressor,'kJ/kg'],
        ['Entalpia de entrada no evaporador:',h_entrada_evaporador,'kJ/kg'],
        ['Entalpia de saida no evaporador:',h_saida_evaporador,'kJ/kg'],
        ['Entalpia de entrada no condensador:',h_entrada_condensador,'kJ/kg'],
        ['Entalpia de saida no condensador:',h_saida_condensador,'kJ/kg']
        ]

BP.print_results(resultados)


#--------------GRAFICO----------------#
pressure_graf = [p_entrada_compressor.nominal_value,
                p_entrada_condensador.nominal_value, 
                p_saida_condensador.nominal_value,
                p_entrada_TC_liq.nominal_value,
                p_saida_TC_liq.nominal_value,
                p_entrada_valvula.nominal_value,
                p_saida_valvula.nominal_value,
                p_entrada_evaporador.nominal_value,
                p_saida_evaporador.nominal_value,
                p_entrada_TC_vap.nominal_value,
                p_saida_TC_vap.nominal_value,
                p_entrada_compressor.nominal_value] #eixo y
h_graf = [h_entrada_compressor,
        h_entrada_condensador,
        h_saida_condensador,
        h_entrada_TC_liq,
        h_saida_TC_liq,
        h_entrada_valvula,
        h_saida_valvula,
        h_entrada_evaporador,
        h_saida_evaporador,
        h_entrada_TC_vap,
        h_saida_TC_vap,
        h_entrada_compressor
        ] #eixo x

point_labels = ['Entrada compressor',
                'Entrada condensador',
                'Saida condensador',
                'Entrada TC liq',
                'Saida TC liq',
                'Entrada valvula',
                'Saida valvula',
                'Entrada evaporador',
                'Saida evaporador',
                ]
BP.plot_series_of_data(h_graf,pressure_graf,'Pxh','h','Pressão','b', point_labels)