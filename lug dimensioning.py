import math
import random

# Material properties
materials = [
    {"name": "4130 Steel", "Fty": 460, "density": 7850, "curve": "curve_151"},  # MPa, kg/m3
    {"name": "8630 Steel", "Fty": 550, "density": 7850, "curve": "curve_151"},    #MPa, kg/m3
    {"name": "2014-T6", "Fty": 414, "density": 2800, "curve": "curve_158"},    #MPa, kg/m3
    {"name": "7075-T6", "Fty": 503, "density": 2810, "curve": "curve_158"}    #MPa, kg/m3

]

# Initial list of the designs
top_designs = []

# Functions for getting k values from the graphs
def curve_151(x):
    y = 0.122*x*6 - 0.5814*x5 + 0.8762*x4 - 0.4578*x3 - 0.1657*x*2 + 1.4777*x
    return y

def curve_158(x):
    if x < 0.44:
        y = 0.54/0.44 * x
    elif 0.44 <= x <= 1.2:
        y = 0.1252*math.log(x) + 0.66
    else:
        y = 0.68
    return y

def curve_121(x):
    if 1 <= x <= 2.9 :
        y = -0.08/1.9 * x + 1.04
    else:
        y = 1.1957 * (math.exp(-0.084 * x))
    return y


# Design variable ranges (in m)
D_min, D_max = 0.02, 0.1  # Pin diameter range
t_min, t_max = 0.005, 0.01  # Thickness range
w_min, w_max = 0.02, 0.2    #width

#Axial and transv forces
Pa = 445
Ptr = 1424.4

# Monte Carlo simulation parameters
num_samples = 1000000  # Number of random design configurations to evaluate

# Monte Carlo simulation
for _ in range(num_samples):
    # Randomly sample D, t, and material
    D = random.uniform(D_min, D_max)
    t = random.uniform(t_min, t_max)
    w = random.uniform(w_min,w_max)
    material = random.choice(materials)

    if w <= D or w-D < 0.005:
        continue  # Skip invalid configurations

    # Compute areas
    A1 = (w - 2*(D/2)*math.sin(45))/2
    A2 = (w-D)/2
    A3 = A2
    A4 = A1

    Abr = D * t
    Aav = 6/(3/A1 + 1/A2 + 1/A3 + 1/A4)
    Aav_Abr = Aav / Abr
    At = (w-D) * t
    w_D = w/D
    
    # Interpolate Ku and Ky for the current material (e.g., Curve 1)
    if material["curve"] == "curve_151":
        Ky = curve_151(Aav_Abr)
    elif material["curve"] == "curve_158":
        Ky = curve_158(Aav_Abr)
    
    Kt = curve_121(w_D)

    # Compute strengths (in N)
    Pty = Ky * Abr * material["Fty"] * 1000000 # Ultimate strength load
    Py = Kt * At * material["Fty"] * 1000000

    Ra = Pa/Py
    Rtr = Ptr/Pty

    MS = (1/((Ra*(1.6) + Rtr(1.6))*0.625))
    # Only continue with the design if MS is positive
    if MS<=0:
        continue

    # Calculate allowable transversal and axial loads
    load_allowable_tr = (Ptr/2) * MS
    load_allowable_a = (Pa/2) * MS

    # Check if design meets load requirements
    if Pty > load_allowable_tr and Py > load_allowable_a:
        # Compute weight
        volume = (w*2 + (w/2)2*math.pi - (D/2)*2*math.pi)*t  
        weight = volume * material["density"]  # kg

        # Define the design for the list
        current_design = {
            "D": D,
            "t": t,
            "w": w,
            "material": material["name"],
            "weight": weight,
            "Pty": Pty,
            "MS": MS
        }
            
        # Add the new design to the list
        top_designs.append(current_design)

# Sort the list by weight and keep the top 10
top_designs = sorted(top_designs, key=lambda x: x["weight"])[:10]

# Output the top 10
print("Top 10 Designs Found:")
for i, design in enumerate(top_designs, 1):
    print(f"Rank {i}: {design}")
