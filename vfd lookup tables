from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd

pumps = ['P01P02', 'P03', 'P04']
K_coeffs = np.array([7.67, 7.07, 7.67])
Sys_Len = np.array([10.93, 2.434, 2.434])

# Parameters
D_080 = 0.080
C_values = [120, 130, 140]

freq_range = np.append(np.arange(36, 58, 1), [58, 58.4, 59, 60])
static_head = np.arange(0.3, 2.45, 0.05)

pump_curve_freqs = np.array([36, 40, 45, 50, 55, 58.4, 60])
pump_curve_coeffs = np.array([
    [-0.175, 4.7],
    [-0.2, 5.8],
    [-0.2, 7.3],
    [-2/9, 9],
    [-1/4, 11],
    [-0.2595, 12.2],
    [-0.2857, 13]
])

interpolated_coeffs = interp1d(pump_curve_freqs, pump_curve_coeffs, axis=0, kind='linear')

flow_rates = []
tank_level = []

for Name, K, L_080 in zip(pumps, K_coeffs, Sys_Len):

    for freq in freq_range:
        for sh in static_head:
            if freq in pump_curve_freqs:
                coeff_idx = np.where(pump_curve_freqs == freq)[0][0]
                coeff = pump_curve_coeffs[coeff_idx]
            else:
                coeff = interpolated_coeffs(freq)  # Get the interpolated coefficients for the current frequency
            for C in C_values:
                def equations(p):
                    h , q = p
                    eq1 = coeff[0]*q*1000 + coeff[1] - h
                    eq2 = (10.67*(q**1.852)*L_080)/((C**1.852)*(D_080**4.8704)) + K * (q / (np.pi * (D_080 / 2)**2)) ** 2 / (2 * 9.81) - h + sh
                    return(eq1, eq2)

                h , q = fsolve(equations,(10 , 0.015))
                flow_rates.append({'Frequency (Hz)': freq, 'C': C, 'Flow Rate (L/s)': q*1000, 'Dynamic Head (m)': h, 'Tank Level (m)':sh})

    flow_rate_table = pd.DataFrame(flow_rates)

    flow_rate_table.to_csv('test.csv')
    flow_rate_dict = {}

    for C_Val in C_values:
        flow_rate_dict[C_Val] = flow_rate_table[flow_rate_table['C'] == C_Val]
        flow_rate_dict[C_Val] = np.round(flow_rate_dict[C_Val].pivot_table(index='Tank Level (m)', columns='Frequency (Hz)', values='Flow Rate (L/s)'), decimals=2)
        flow_rate_dict[C_Val].index = 2.8 - flow_rate_dict[C_Val].index
        flow_rate_dict[C_Val] = flow_rate_dict[C_Val][::-1].rename_axis(f'Tank Level/Hz {Name} C{C_Val}')
        flow_rate_dict[C_Val].to_csv(f'test_{Name}_system_C{C_Val}.csv')
        print(flow_rate_dict[C_Val])
