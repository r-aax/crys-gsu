from scipy.interpolate import make_interp_spline
from scipy.optimize import curve_fit

_KilocalorieDjoules = 4186.8;
_DjouleKilocalories = (1.0 / _KilocalorieDjoules);
_StefanBoltzmannConstant = 5.670374419e-8;

# Viscosity by temperature.
_water_dynamic_viscosity_data = {'x':
                                 [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
                                 'y':
                                 [1.787e-3, 1.519e-3, 1.307e-3, 1.002e-3, 0.798e-3,
                                  0.653e-3, 0.547e-3, 0.467e-3, 0.404e-3, 0.355e-3,
                                  0.315e-3, 0.282e-3]
                                 }

# Capacity by temperature.
_water_specific_heat_capacity_data = {'x':
                                      [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0,
                                       1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0, 11.0,
                                       12.0, 14.0, 15.0, 17.0, 18.0, 22.0, 23.0, 29.0],
                                      'y':
                                      [1.019, 1.018, 1.017, 1.016, 1.015, 1.014, 1.013, 1.012, 1.011, 1.01,
                                       1.009, 1.008, 1.008, 1.007, 1.006, 1.005, 1.004, 1.004, 1.003, 1.003,
                                       1.002, 1.002, 1.001, 1.001, 1.0, 1.0, 0.999, 0.999]}
_water_specific_heat_capacity_data['y'] = [v * _KilocalorieDjoules for v in _water_specific_heat_capacity_data['y']]

# Capacity by temperature.
_ice_specific_heat_capacity_data = {'x':
                                    [-28.0, -27.0, -26.0, -25.0, -24.0, -23.0, -22.0, -21.0, -20.0, -19.0,
                                     -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0, -10.0],
                                    'y':
                                    [0.452, 0.454, 0.455, 0.457, 0.459, 0.461, 0.463, 0.466, 0.467, 0.468,
                                     0.470, 0.472, 0.474, 0.476, 0.478, 0.480, 0.481, 0.483, 0.485]}
_ice_specific_heat_capacity_data['y'] = [v * _KilocalorieDjoules for v in _ice_specific_heat_capacity_data['y']]

# Capacity by temperature.
_air_specific_heat_capacity_data = {'x':
                                    [-50.0, -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0,  -5.0,
                                     0.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
                                     90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0,
                                     190.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0],
                                    'y':
                                    [1013.0, 1013.0, 1013.0, 1013.0, 1013.0, 1011.0, 1009.0, 1009.0, 1009.0, 1007.0,
                                     1005.0, 1005.0, 1005.0, 1005.0, 1005.0, 1005.0, 1005.0, 1005.0, 1009.0, 1009.0,
                                     1009.0, 1009.0, 1009.0, 1009.0, 1011.0, 1013.0, 1015.0, 1017.0, 1020.0, 1022.0,
                                     1024.0, 1026.0, 1037.0, 1047.0, 1058.0, 1068.0, 1081.0, 1093.0, 1104.0]}

# Heat by temperature.
_vaporization_specific_heat_data = {'x':
                                    [0.0, 10.0, 20.0, 30.0, 50.0, 70.0, 90.0, 100.0, 120.0, 150.0,
                                     180.0, 200.0, 220.0, 250.0, 300.0, 350.0, 370.0, 374.0, 374.15],
                                    'y':
                                    [2500, 2470, 2450, 2400, 2380, 2320, 2280, 2260, 2200, 2110,
                                     2010, 1940, 1860, 1700, 1400, 890, 440, 110, 0]}

# Heat by temperature.
_sublimation_specific_heat_data = {'x':
                                   [-39.0, -38.0, -37.0, -36.0, -35.0, -34.0, -33.0, -32.0, -31.0, -30.0,
                                    -29.0, -28.0, -27.0, -26.0, -25.0, -24.0, -23.0, -22.0, -21.0, -20.0,
                                    -19.0, -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0, -10.0,
                                    -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0],
                                   'y':
                                   [698.4, 697.8, 697.3, 696.8, 696.2, 695.6, 695.1, 694.5, 694.0, 693.4,
                                    692.9, 692.3, 691.8, 691.2, 690.6, 690.0, 689.5, 688.9, 688.6, 688.1,
                                    687.6, 687.0, 686.5, 685.9, 685.3, 684.7, 684.2, 683.6, 683.1, 682.5,
                                    681.9, 681.4, 680.8, 680.3, 679.7, 679.1, 678.5, 678.0, 677.6, 676.9]}

_sublimation_specific_heat_data['y'] = [v * _KilocalorieDjoules for v in _sublimation_specific_heat_data['y']]

class LinearLaw()
    def __init__(self, table):
        self.interpolate = make_interp_spline(table['x'], table['y'], 1)

    def __call__(self, x):
        return self.interpolate(x)

def _cubic(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

def _quadro(x, a, b, c, d, e):
    return a * x**4 + b * x**3 + c * x**2 + d * x + e

class FuncLaw():
    def __init__(self, table, func):
    """
    :param func: function with parameters to be optimized where 1st parameter is x and the rest are parameters to be optimized
    """
        self.func = func
        self.popt, self.pcov = curve_fit(self.func, table['x'], table['y'])

    def __call__(self, x):
        return self.func(x, *(self.popt))

class ConstLaw():
    def __init__(self, constant):
        self.y = constant

    def __call__(self, x):
        return self.y
