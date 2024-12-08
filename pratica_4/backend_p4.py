import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import math
from math import pi,sqrt,exp 
from uncertainties import ufloat

'''
Funções criadas para prática 2
'''

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
    incerteza = incerteza*6891.76
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

'''
Funções criadas para prática 3
'''

def iteration_for_q(h, L, d, Tb, Tamb, K, show = True):
    """
    Calcula a taxa de transferência de calor na aleta.
    
    Parâmetros:
    - h: Coeficiente de transferência de calor convectivo (W/m²K)
    - L: Comprimento da aleta (m)
    - d: Diâmetro da aleta (m)
    - Tb: Temperatura da base da aleta (K)
    - Tamb: Temperatura do ambiente (K)
    - K: Condutividade térmica do material da aleta (W/mK)
    
    Retorna:
    - q_aleta: Taxa de transferência de calor (W)
    """
    
    # Calcular o perímetro e a área da seção transversal da aleta
    P = pi * d  # Perímetro da aleta (m)
    A = (pi * d**2) / 4  # Área da seção transversal da aleta (m²)
    emissividade = 0.6
    stephan_boltzman = 5.67*10**(-8)
    erro = 10
    while erro >= 0.01:
        # Calcular m
        m = (h * P / (K * A))**(1/2)
        
        # Calcular q pela equação 3.77
        q_aleta = (h * P * K * A)**(1/2) * (Tb - Tamb) * \
                (np.sinh(m.nominal_value * L.nominal_value) + (h / (m * K)) * np.cosh(m.nominal_value * L.nominal_value)) / \
                (np.cosh(m.nominal_value * L.nominal_value) + (h / (m * K)) * np.sinh(m.nominal_value * L.nominal_value))
        
        T_media = Tamb + (q_aleta/(h*P*L))

        h_conv_teórico_laminar = 1.32*((T_media.nominal_value - Tamb.nominal_value)/d.nominal_value)**(0.25)
        h_rad_teórico = emissividade*stephan_boltzman*(T_media.nominal_value**2 + Tamb.nominal_value**2)*(T_media.nominal_value + Tamb.nominal_value)
        h_recalculado = h_conv_teórico_laminar + h_rad_teórico
        erro = np.absolute(h_recalculado - h)
        h = h_recalculado
        if show:
            # Exibindo resultados com precisão ajustada
            print(f"Erro: {erro:.4f}")
            print(f"h_recalculado: {h_recalculado:.4f}")
            print(f"q_aleta: {q_aleta:.4f}\n")  # Exibir mais casas decimais
    return q_aleta, h

def temperature_distribution_on_bar(h, L, d, Tb, Tamb, K, x):
    """
    Calcula a temperatura em um ponto x ao longo de uma aleta.
    
    Parâmetros:
    - h: Coeficiente de transferência de calor convectivo (W/m²K)
    - L: Comprimento da aleta (m)
    - d: Diâmetro da aleta (m)
    - Tb: Temperatura da base da aleta (K)
    - Tamb: Temperatura do ambiente (K)
    - K: Condutividade térmica do material da aleta (W/mK)
    - x: Distância ao longo da aleta (m)

    Retorna:
    - T(x): Temperatura no ponto x (K)
    """
    # Calcular o perímetro (P) e a área da seção transversal (A) da aleta
    P = pi * d  # Perímetro da aleta (m)
    A = (pi * d**2) / 4  # Área da seção transversal da aleta (m²)

    # Calcular m
    m = (h * P / (K * A))**0.5

    # Termos do numerador e denominador
    numerador = math.cosh(m.nominal_value * (L.nominal_value - x)) + (h / (m * K)) * math.sinh(m.nominal_value * (L.nominal_value - x))
    denominador = math.cosh(m.nominal_value * L.nominal_value) + (h / (m * K)) * math.sinh(m.nominal_value * L.nominal_value)

    # Cálculo da temperatura no ponto x
    T_x = Tamb + (Tb - Tamb) * (numerador / denominador)
    return T_x

def t_distribuitions(x,y,title,x_lable,y_lable,y_lable_graf,colors):
    plt.figure(figsize=(12, 8))
    for a,b,c,d in zip(y,y_lable_graf,colors,x):
        plt.plot(d, a, linestyle='-', color=f'{c}', label=f'{b} vs {x_lable}')
    # Adicionando títulos e rótulos
    plt.title(title)
    plt.xlabel(x_lable)
    plt.ylabel(y_lable)
    plt.legend()
    plt.grid(True)
    plt.show()

def t_diffs(x,y,title,x_lable,y_lable,colors):
    plt.figure(figsize=(12, 8))
    plt.plot(x, y, marker = 'o', linestyle='-', color=f'{colors}', label=f'{y_lable} vs {x_lable}')
    # Adicionando títulos e rótulos
    plt.title(title)
    plt.xlabel(x_lable)
    plt.ylabel(y_lable)
    plt.legend()
    plt.grid(True)
    plt.show()

'''
Funções criadas para prática 4
'''

def convert_hPa_to_Pa(pressure,incerteza):
    '''
    Rescales hPa pressures to Pa
    '''
    pressure = ufloat(pressure*100,incerteza)
    return pressure

def convert_time_to_seconds(time_string):
    """
    Converte um tempo no formato 'min:s.ms' em segundos.

    Args:
        time_string (str): Tempo no formato 'min:s.ms'.

    Returns:
        float: O tempo em segundos.
    """
    # Divide a string em minutos e segundos
    minutes, seconds = time_string.split(':')
    # Converte para float e calcula o total em segundos
    total_seconds = int(minutes) * 60 + float(seconds)
    return total_seconds

def convert_liters_to_cubic_meters(liters,incerteza):
    """
    Converte um volume de litros (L) para metros cúbicos (m^3).

    Args:
        liters (float): O volume em litros.

    Returns:
        float: O volume em metros cúbicos.
    """
    result = ufloat(liters/1000,incerteza)
    return result

def convert_milliliters_to_cubic_meters(milliliters,incerteza):
    """
    Converte um volume de mililitros (mL) para metros cúbicos (m^3).

    Args:
        milliliters (float): O volume em mililitros.

    Returns:
        float: O volume em metros cúbicos.
    """
    result = ufloat(milliliters * 10**-6,incerteza)
    return result