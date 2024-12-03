import numpy as np
from copy import deepcopy

# authors: Johannes Nilsson, __insert here__ of D04
# prupose: iteratively find best fastener setup for the lug
# version: 0.3
# date: 2024-11-28

#process:
# use a rough monte carlo method to get an approximate minimum.
# then use gradient descent to find the local minimum
class Fastener():
    def __init__(self, posx:float, posz:float, diameter:float) -> None:
        self.posx = posx
        self.posz = posz
        self.diameter = diameter

    def __repr__(self):
        return f"x position: {self.posx}, z position: {self.posz}, diameter: {self.diameter}"

class Setup():
    def __init__(self, n_pairs:int, thickness:float, fasteners:list[Fastener], mat_lug:int):
        self.n_pairs = n_pairs #number of pairs of fasteners
        self.thickness_lug = thickness #thickness of the lug
        self.fasteners = fasteners #position of the fasteners
        self.max_bearing_stress = None # to be calculated
        self.material_lug = mat_lug #material index of the attaching plate
        self.under_force_limit = True #assumed true unless changed

    def gen_monte_carlo():
        """generates a random setup"""
        #generate fasteners
        n_pairs = np.random.randint(2, max_pairs + 1)
        diameter = np.random.uniform(D_min,D_max)
        fasteners = []
        for i in range(n_pairs):
            fasteners.append(Fastener(np.random.uniform(x_min + diameter, x_max - diameter), np.random.uniform(z_min + diameter, z_max- diameter), diameter)) 
            #+- diameter to make sure it doesn't intersect with the edge
        new_setup = Setup(n_pairs,
                          np.random.uniform(*thickness_range),  # lug thickness (m)
                          fasteners,
                          mat_lug=np.random.randint(0,len(material_list))
                          )
        return new_setup


    def eval_mass(self):
        """evaluates the mass of the setup\n
        and the max force as well"""
        # This is all following WP4 instructions
        # Calculate areas
        for f in self.fasteners:
            f.area = np.pi * (f.diameter / 2) ** 2
        self.Area_fasteners_total = sum(f.area for f in self.fasteners)
        
       
        
        # Mass calculations
        lug_mass = W * H * self.thickness_lug * material_list[self.material_lug][1]  # (kg)
        total_fastener_mass = self.n_pairs* 2 * fastener_mass  # (kg)
        total_mass = lug_mass + total_fastener_mass
        return total_mass
    
    def force_check(self):
        """Checks if all the forces are under the limit"""
         # Find CG of the fasteners
        #x_cg = sum(f.area * f.posx for f in self.fasteners) / Area_fasteners_total
        x_cg = 0 #since we assume symmetry around z axis
        z_cg = sum(f.area * f.posz for f in self.fasteners) / self.Area_fasteners_total

        # Calculate distances to CG
        radius = lambda fastener : np.sqrt((fastener.posx-x_cg)**2 + (fastener.posz - z_cg)**2)
        
        # Forces and moments eq (4.2-4.4)
        F_cg = np.sqrt(F_x**2 + F_z**2)
        F_cgx = F_x
        F_cgz = F_z
        M_cgy = F_x * z_cg - F_z * x_cg


        denominator = sum(f.area * radius(f)**2 for f in self.fasteners)

        if denominator == 0:
            raise SyntaxError
   
        # Force calculations (bearing check continued)
        for f in self.fasteners:
            f.F_in_plane_x = F_cgx / self.n_pairs*2
            f.F_in_plane_z = F_cgz / self.n_pairs*2
            f.F_in_plane_My = (M_cgy * f.area * radius(f)) / denominator
            force_angle = np.arctan2(f.posz - z_cg,f.posx - x_cg) #angle of the moment force on the bearing
            f.F_total = np.sqrt((f.F_in_plane_x - np.cos(force_angle)*f.F_in_plane_My)** 2 + (f.F_in_plane_z + np.sin(force_angle)*f.F_in_plane_My) ** 2)/self.n_pairs*2 # total force magnitude on the fasteners
        
        # Bearing stress
        for f in self.fasteners:
            f.sigma_bearing = f.F_total / (f.diameter * self.thickness_lug)
        
        # Check if all fasteners allow the bearing stress
        max_fastener_bearing_stress = max(f.sigma_bearing for f in self.fasteners)
        self.max_bearing_stress = max_fastener_bearing_stress

        #is this above the limit
        if self.max_bearing_stress > material_list[self.material_lug][2] * safety_factor:
            self.under_force_limit = False

        #check pull through 

        F_pi = F_x/(self.n_pairs*2)
        denominator = sum(f.area * radius(f)**2 for f in self.fasteners)

        for f in self.fasteners: #(4.6-4.7)
            F_pMz = (M_x * f.area * radius(f)) / denominator
            f.F_pull_through = F_pi + F_pMz

        self.max_pull_through = max(f.F_pull_through for f in self.fasteners)

        #now find if the walls themselves resist the load
        #assuming the same diameter for the hole and the fastener itself
        for f in self.fasteners:

            crown_area = (np.pi/4)*((f.diameter*(1+fastener_head_ratio))**2-f.diameter**2)
            shaft_area = np.pi*f.diameter*self.thickness_lug
            stress = f.F_pull_through/crown_area
            strain = f.F_pull_through/shaft_area
            f.von_mises_stress = np.sqrt(stress**2 + 3*strain**2) #the tension yield stress.

        # Check if all fasteners allow the bearing stress
        self.max_wall_pull_through = max(f.von_mises_stress for f in self.fasteners)
        #is this above the limit
        if self.max_wall_pull_through > material_list[self.material_lug][2] * safety_factor:
            self.under_force_limit = False

        return

    def gradient_descent(self):
        """improves the design to minimize mass using gradient descent\n
        currently modifies thickness and the position, diameter of the fasteners"""

        #get the partial derivatives of each section and move the design in that direction
        total_delta_mass = 500
        step_size = 1e-2 # for actually moving the value
        delta_step = 1e-5 #for testing the derivative
        debug_count = 0 #counts number of iterations
        while (debug_count < 1000):

            total_old_mass = self.eval_mass()
            old_design = deepcopy(self)
            
            #step thickness:
            old_thickness = self.thickness_lug
            old_mass = self.eval_mass()

            #change setup:
            self.thickness_lug += delta_step
            new_mass = self.eval_mass()

            #move in direction of minimum mass
            self.thickness_lug = old_thickness + step_size*(old_mass-new_mass)

            #if we run into a limit stop it:
            if not (thickness_range[0] < self.thickness_lug < thickness_range[1]):
                self.thickness_lug = old_thickness

            #now do this with position of the fasteners:
            for f in self.fasteners:
                #in x
                old_posx = f.posx
                old_mass = self.eval_mass()

                #change x:
                f.posx += delta_step
                new_mass = self.eval_mass()

                f.posx = old_posx + step_size*(old_mass-new_mass)
                #if we run into a limit stop it:
                if not (x_min < f.posx < x_max):
                    f.posx = old_posx
                #in z
                old_posz = f.posz
                old_mass = self.eval_mass()

                #change z:
                f.posz += delta_step
                new_mass = self.eval_mass()

                f.posz = old_posz + step_size*(old_mass-new_mass)
                #if we run into a limit stop it:
                if not (z_min < f.posz < z_max):
                    f.posz = old_posz

                #diameter
                old_diameter = f.diameter
                old_mass = self.eval_mass()

                #change diameter:
                f.diameter += delta_step
                new_mass = self.eval_mass()

                f.diameter = old_diameter + step_size*(old_mass-new_mass)
                #if we run into a limit stop it:
                if not (D_min < f.diameter < D_max):
                    f.diameter = old_diameter
            #get total delta mass:
            total_delta_mass = total_old_mass - self.eval_mass()

            #check the forces:
            self.force_check()
            if self.under_force_limit == False:
                self = old_design

            #lower step size to slowly bring it to equilibrium
            step_size *= 0.5
            debug_count += 1
        #process complete
        print(f"Debug: {debug_count} iterations before completion")

        return #it edits in place so no need to return anything
        




    def __repr__(self) -> str:
        return f"lug mass: {self.eval_mass()}, thickness: {self.thickness_lug}, number of fasteners: {self.n_pairs*2}, max bearing stress: {self.max_bearing_stress}, material: {material_list[self.material_lug][0]}, max pull through: {self.max_pull_through}"



# Constants
F_x = 366.3/2
F_z = 1424.4/2
M_x = 57.2
fastener_mass = 0.04  # kg
density = 2810  # aluminum 7075-T6 (kg/m^3)
max_bearing_stress = 300e6  # 450 yield stress and S.F. = 1.5
max_pull_through_limit = 500
thickness_of_attachment_plate = 0.001
fastener_head_ratio = 0.5
safety_factor = 1.5
#name, density, max stress
material_list=[
    ["Aluminium 2014-T6", 2796, 414e6],
    ["Aluminium 7075-T6", 2823, 450e6],
    ["Steel 4130", 7850, 435e6],
    ["Steel 8630", 7850, 550e6],
    ["Magnesium AZ91C-T6", 1810, 145e6],
    ["Aluminium 356.0-T6", 2671, 152e6]
]
'''name,density,max stress'''
# Lug dimensions (back-up wall)
W = 0.5  # x-axis (m)
H = 0.3  # z-axis (m)

# Design Constraints
max_fasteners = 10
max_pairs = max_fasteners // 2
x_min, x_max = 0, W / 2    # m (the/2 for symmetry)
z_min, z_max = 0, H        # m
D_min, D_max = 0.001, 0.01  # m
thickness_range = [0.01, 0.05] # m

# Number of designs to evaluate (iterations)
iterations = 100000

# Optimization (Monte Carlo Method)
best_setup = None
best_mass = float('inf')

for _ in range(iterations):
    # Randomly generate a design
    tested_setup = Setup.gen_monte_carlo()
    mass = tested_setup.eval_mass()
    tested_setup.force_check()
    if tested_setup.max_bearing_stress > max_bearing_stress:
        continue #too high bearing stress

    # Check if this design is better than previous
    if mass < best_mass:
        best_mass = mass
        best_setup = tested_setup

# Print the best design

print("monte carlo best setup:")
print(best_setup)
for f in best_setup.fasteners:
    print(f)

best_setup.gradient_descent()
best_setup.force_check()
assert best_setup.under_force_limit == True
print("Gradient descent best setup:")
print(best_setup)
for f in best_setup.fasteners:
    print(f)
