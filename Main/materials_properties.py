import numpy as np
from scipy import interpolate as scipy_interpolate
import inspect


def initiate_thermal_conductivity_interpolation():

    # "tungsten":
    # temperature_c =       [20, 100, 200, 300, 400, 500, 600, 700]
    temperature_k =        [293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15]
    thermal_conductivity = [172.8,  164.8,  155.5,  147.2,  139.8,  133.1,  127.2,  122.1]
    interpolated_tungsten_thermal_conductivity = scipy_interpolate.interp1d(temperature_k, thermal_conductivity) 

    # "lithium_lead":
    #temperature_c =       [20,     300,    350,    400,    450,    500,    550,    600,    650,    700]
    temperature_k =        [293.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15, 923.15, 973.15]
    thermal_conductivity = [7.69,   13.18,  14.16,  15.14,  16.12,  17.10,  18.08,  19.06,  20.04,  21.02]
    interpolated_lithium_lead_thermal_conductivity = scipy_interpolate.interp1d(temperature_k, thermal_conductivity) 

    # "eurofer":
    #temperature_c =       [20,     50,     100,    150,    200, 250, 300, 350, 400, 450, 500, 550,600]
    temperature_k =        [293.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15]
    thermal_conductivity = [27.63,  28.73,  29.87,  30.32,  30.28,  29.95,  29.51,  29.10,  28.84,  28.82,  29.08,  29.62,  30.38]
    interpolated_eurofer_thermal_conductivity = scipy_interpolate.interp1d(temperature_k, thermal_conductivity) 

    return {"tungsten": interpolated_tungsten_thermal_conductivity,
            "lithium_lead": interpolated_lithium_lead_thermal_conductivity,
            "eurofer": interpolated_eurofer_thermal_conductivity}



def initiate_specific_heat_interpolation():

    # "tungsten":
    # temperature_c = [20, 100, 200, 300, 400, 500, 600, 700]
    temperature_k =[293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15]
    specific_heat = [129, 131.6, 134.7, 137.8, 140.9, 133.1, 127.2, 122.1]
    interpolated_tungsten_specific_heat = scipy_interpolate.interp1d(temperature_k, specific_heat) # this object could be created once on inititation to speed up the code

    # "lithium_lead":
    # temperature_c = [20, 300, 350, 400, 450, 500, 550, 600, 650, 700]
    temperature_k = [293.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15, 923.15, 973.15]
    specific_heat = [192, 190, 189, 189, 188, 188, 187, 187, 187, 186]
    interpolated_lithium_lead_specific_heat = scipy_interpolate.interp1d(temperature_k, specific_heat) # this object could be created once on inititation to speed up the code

    # "eurofer":
    # temperature_c = [20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
    temperature_k =  [293.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15]
    specific_heat =[439, 462, 490, 509, 523, 534, 546, 562, 584, 616, 660, 721, 800]
    interpolated_eurofer_specific_heat = scipy_interpolate.interp1d(temperature_k, specific_heat) # this object could be created once on inititation to speed up the code

    return {"tungsten": interpolated_tungsten_specific_heat,
            "lithium_lead": interpolated_lithium_lead_specific_heat,
            "eurofer": interpolated_eurofer_specific_heat}


def initiate_density_interpolation():

    # "tungsten":
    # temperature_c = [20, 100, 200, 300, 400, 500, 600, 700]
    temperature_k =[293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15]
    density = [19298, 19279, 19254, 19229, 19205, 19178, 19152, 19125 ]
    interpolated_tungsten_density = scipy_interpolate.interp1d(temperature_k, density) 

    # "lithium_lead":
    # temperature_c = [20, 300, 350, 400, 450, 500, 550, 600, 650, 700]
    temperature_k = [293.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15, 923.15, 973.15]
    density = [10172, 9839, 9779, 9720, 9661, 9601, 9542, 9482, 9423, 9363]
    interpolated_lithium_lead_density = scipy_interpolate.interp1d(temperature_k, density)

    # "eurofer":
    # temperature_c = [20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
    temperature_k = [293.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15]
    density =       [7760,   7753,   7740,   7727,   7713,   7699,   7685,   7670,   7655,   7640,   7625,   7610, 7594]
    interpolated_eurofer_density = scipy_interpolate.interp1d(temperature_k, density) 

    return {"tungsten": interpolated_tungsten_density,
            "lithium_lead": interpolated_lithium_lead_density,
            "eurofer": interpolated_eurofer_density}


thermal_conductivities_dict = initiate_thermal_conductivity_interpolation()
specific_heat_dict = initiate_specific_heat_interpolation()
density_dict = initiate_density_interpolation()

def calculate_D(T, material_id):
    R = 8.314  # Perfect gas constant
    k_B = 8.6e-5
    if material_id == "concrete":  # Concrete
        return 2e-6  # 7.3e-7*np.exp(-6.3e3/T)
    elif material_id == "polymer":  # Polymer
        return 2.0e-7*np.exp(-29000.0/R/T)
    elif material_id == "steel":  # steel
        return 7.3e-7*np.exp(-6.3e3/T)
    elif material_id == "tungsten":
        return 4.1e-7*np.exp(-0.39/k_B/T)
    elif material_id == "eurofer":
        return 8.1e-8*np.exp(-14470/R/T)
    elif material_id == "lithium_lead":
        return 2.5e-7*np.exp(-27000/R/T)
    else:
        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))


def calculate_thermal_conductivity(T, material_id):
    if material_id == "concrete":
        return 0.5
    if material_id == "polymer":
        return 0.3
    if material_id == "steel":
        return 16
    if material_id not in thermal_conductivities_dict.keys():
        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))
    else:
        interpolated_object = thermal_conductivities_dict[material_id]
        return float(interpolated_object.__call__(T))


def calculate_specific_heat(T, material_id):
    if material_id == "concrete":
        return 880
    if material_id == "polymer":
        return 1000
    if material_id == "steel":
        return 500
    if material_id not in specific_heat_dict.keys():
        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))
    else:
        interpolated_object = specific_heat_dict[material_id]
        return float(interpolated_object.__call__(T))


def calculate_density(T, material_id):
    if material_id == "concrete":
        return 2400
    if material_id == "polymer":
        return 940
    if material_id == "steel":
        return 7700
    if material_id not in density_dict.keys():
        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))
    else:
        interpolated_object = density_dict[material_id]
        return float(interpolated_object.__call__(T))


def calculate_mu(T, material_id):
    if material_id == 'lithium_lead':
        return 0.187e-3*np.exp(1400/T)
