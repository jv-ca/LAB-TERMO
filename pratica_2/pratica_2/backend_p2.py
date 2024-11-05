import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import math
from math import pi,sqrt,exp 
from uncertainties import ufloat

def convert_celsius_to_kelvin(temp, incerteza):
    '''
    Converte °C para K, retornando o ufloat
    '''
    temp = temp + 273.15
    ufloat_temp = ufloat(temp, incerteza)
    return ufloat_temp

def convert_psi_to_pascal(pressure, incerteza):
    '''
    Converte psi para Pascal, retornando ufloat
    '''
    pressure = pressure*6891.76
    ufloat_pressure = ufloat(pressure, incerteza)
    return ufloat_pressure

def get_entalpy(pressure, temperature, fluid = 'R12'):
    h = CP.PropsSI('H', 'T', temperature.nominal_value, 'P', pressure.nominal_value, fluid)
    return h

def print_results(results):
    for i in results:
        print(f'{i[0]} {i[1]} {i[2]}')

def lb_min_to_kg_s(mass_flow_lb_min):
    """
    Converte uma vazão mássica de lb/min para kg/s.
    
    Parameters:
    - mass_flow_lb_min: Vazão mássica em lb/min.

    Returns:
    - Vazão mássica em kg/s.
    """
    # Converte lb/min para kg/s (1 lb = 0.453592 kg)
    mass_flow_kg_s = mass_flow_lb_min * 0.453592 / 60
    return mass_flow_kg_s

def plot_series_of_data(x,y,title,x_lable,y_lable,color,point_labels):
    # Plotando a curva de sino no gráfico de Pressão vs Entalpia
    h_liq, h_vap, P_liq, P_vap = R12_curve_Pxh()
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker= 'o',linestyle='-', color=f'{color}', label=f'{y_lable} vs {x_lable}')
    plt.plot(h_liq, P_liq, linestyle='-', color=f'red', label=f'Curva de sino do R12')
    plt.plot(h_vap, P_vap, linestyle='-', color=f'red')
    # Adicionando uma anotação em cada ponto
    for i, (xi, yi) in enumerate(zip(x, y)):
        if i == len(x) - 1:  # Ignora o último ponto
            break
        # Alterna o valor de 'xytext' entre (5, 5) e (5, -15)
        offset = (0, 10) if i % 2 == 0 else (0, -8)
        if i == 11:
            offset = (0,-8)
        plt.annotate(str(i + 1),          # Rótulo numérico (1, 2, 3, ...)
                    (xi, yi),            # Coordenadas do ponto
                    textcoords="offset points", # Tipo de posicionamento
                    xytext=offset,       # Offset alternado
                    ha='center')         # Alinhamento horizontal
    # Adicionando títulos e rótulos
    plt.title(title)
    plt.xlabel(x_lable)
    plt.ylabel(y_lable)
    plt.legend()
    plt.grid(True)
    plt.show()

def R12_curve_Pxh():
    # Definir o fluido
    fluid = 'R12'

    # Intervalo de temperaturas para a curva de sino (Kelvin)
    T_min = CP.PropsSI('Tmin', fluid)    # Temperatura mínima possível para o R12
    T_crit = CP.PropsSI('Tcrit', fluid)  # Temperatura crítica do R12

    # Lista de temperaturas para amostrar a curva de sino
    temperatures = np.linspace(T_min+50, T_crit, 100)

    # Listas para armazenar entalpias e pressões para os estados de líquido e vapor saturados
    h_liq = []
    h_vap = []
    P_liq = []
    P_vap = []

    # Calculando entalpias e pressões para cada temperatura
    for T in temperatures:
        # Entalpia do líquido e do vapor saturado
        h_liq.append(CP.PropsSI('H', 'T', T, 'Q', 0, fluid) / 1000)  # Convertendo para kJ/kg
        h_vap.append(CP.PropsSI('H', 'T', T, 'Q', 1, fluid) / 1000)  # Convertendo para kJ/kg

        # Pressão do líquido e do vapor saturado (em MPa)
        P_liq.append(CP.PropsSI('P', 'T', T, 'Q', 0, fluid))  # Convertendo para MPa
        P_vap.append(CP.PropsSI('P', 'T', T, 'Q', 1, fluid))  # Convertendo para MPa

    # Convertendo listas para arrays do NumPy para facilitar o uso
    h_liq = np.array(h_liq).tolist()
    h_vap = np.array(h_vap).tolist()
    P_liq = np.array(P_liq).tolist()
    P_vap = np.array(P_vap).tolist()
    return h_liq,h_vap,P_liq,P_vap

