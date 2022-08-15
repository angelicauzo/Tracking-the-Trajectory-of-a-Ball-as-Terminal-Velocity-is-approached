"""

author: Angelica Uzo
course: Chemical Engineering
school: University of Birmingham

"""

# This code will model the trajectory of a hard spherical steel ball with drag forces acting 
# travelling through different fluids provided the initial position, initial speed and angle of inclination 
# to the horizontal in degrees at an initial time, t0 in steps dt until final time, tf is reached
import numpy as np
import matplotlib.pyplot as plt
import seaborn; seaborn.set_style("whitegrid")

# Pre-defined parameters
# Initial position
rx0 = 0 #m
ry0 = 100 #m

# Initial speed and angle of inclination to the horizontal
vx = 50 #m s^-1
vy = 0 #m s^-1

# Gravitational acceleration
g = 9.81 #m s^-2

# Initial time
t0 = 0 #s

# Drag coefficient
cd_sphere = 0.47 #dimensionless 

# Density of fluids
rho_air = 1.2 #kg m^-3 
rho_water = 1000 #kg m^-3 
rho_glycerin = 1260 #kg m^-3 

# Nylon ball dimensions
diameter_nylon = 5e-3 #m
rho_nylon = 1140 #kg m^-3

# Steel ball dimensions
diameter_steel = 5e-3 #m
rho_steel = 7850 #kg m^-3

# User-defined parameters
# Time step                          
dt = 0.01 #s
# Final time
tf = 10000 #s 

# Current working parameters
rho_solid = rho_steel
rho_fluid = rho_air
diameter_solid = diameter_steel

# Spherical ball
radius_solid = diameter_solid / 2 #m
projected_area_solid = np.pi * radius_solid ** 2 #m^2
volume_solid = np.pi * 4 / 3 * radius_solid ** 3 #m^3

# Buoyancy
force_buoyancy = np.array([0, rho_fluid * volume_solid * g])

# Weight of ball
mass_ball = volume_solid * rho_solid
force_weight = np.array([0, mass_ball * -g])

# This fuction returns the drag force at the velocity
def drag_force(velocity):
    drag = cd_sphere * projected_area_solid * rho_fluid * velocity ** 2 / 2
    return drag

# t_current represents the time at the current position of the projectile
t_current = t0

# r_current is an array containing the current vertical and horizontal displacements of the projectile respectively 
r_current = np.array([rx0 , ry0]) 

# v_current is an array containing the current vertical and horizontal velocities of the projectile respectively
v_current = np.array([vx, vy])


# position, time and speed represent empty matrices into which the r_current, t_current and v_current
# values will be appended respectively
position = []
time = []
speed = []
terminal = []

# Euler's method
# This loop calculates r_current and v_current at t_current and appends it to the list 'position' and 'speed' 
# respectively until t_current is equal to tf after which the loop is terminated. 
while t_current <= tf:
    # Calculating acceleration
    # a[0] is the acceleration along the horizontal whereas a[1] is the acceleration along the vertical
    v0_magitude = np.sqrt(v_current[0]**2 + v_current[1]**2)
    force_drag = np.array([drag_force(v_current[0]) * (-v_current[0]/v0_magitude), drag_force(v_current[1]) * (-v_current[1]/v0_magitude)])
    resultant_force = force_weight + force_buoyancy + force_drag
    a = np.array([resultant_force[0] / mass_ball, resultant_force[1] / mass_ball])
    # r_current[0] represents the vertical displacement, r_current[1] represents the horizontal displacement
    # v_current[0] represents the vertical velocity, v_current[1] represents the horizontal velocity
    v_new = np.array([v_current[0] + dt * a[0] , v_current[1] + dt * a[1]])
    r_new = np.array([r_current[0] + v_current[0] * dt , r_current[1] + v_current[1] * dt])
    # Error margin of 0.001
    E = 0.001
    # B represents the magnitude of v_new
    B = np.sqrt(v_new[0]**2+v_new[1]**2)
    # C represents the magnitude of v_current
    C = np.sqrt(v_current[0]**2+v_current[1]**2)
    # Finding terminal velocity
    if B > C - E and B < C + E:
        terminal.append(B)
        v_terminal = np.array(terminal)
    # 'position.append(r_current)' modifies the list 'position' by adding r_current to the end of the list 
    # rx_and_ry represents an array of the entries within 'position'
    position.append(r_current)
    rx_and_ry = np.array(position)
    # 'speed.append(v_current)' modifies the list 'position' by adding r_current to the end of the list 
    # v represents an array of the entries within 'position'
    speed.append(v_current)
    v = np.array(speed)
    # r_new and v_new become the next timestep's r_current and v_current values
    v_current = v_new
    r_current = r_new
    # 'time.append(t_current)' modifies the list 'time' by adding t_current to the end of the list 
    # t represents an array of the entries within 'time'
    time.append(t_current)
    t = np.array(time)
    # This defines t_current at the new timestep and the loop repeats
    t_current = t_current + dt
    # To end the loop
    if r_current[1] < 0:
        break

terminal_velocity = max(v_terminal)
print("The Euler terminal velocity is", round(terminal_velocity, 4), "m/s")

# Analytical terminal velocity
force_drag = - force_weight + force_buoyancy
analytical_terminal_velocity = np.sqrt((2 * (force_drag) )/(cd_sphere * projected_area_solid * rho_fluid))
print("The analytical terminal velocity is", round(analytical_terminal_velocity[1],4), "m/s")

# Plotting the graphs
fig, axes = plt.subplots(1, 2)
fig.set_size_inches(18, 6)
# rx_and_ry[:,0] represents the horizontal displacement
# rx_and_ry[:,1] represents the vertical displacement
axes[0].scatter(rx_and_ry[:,0], rx_and_ry[:,1], c=t, label = "Euler Projectile Position")
axes[0].set_xlabel("Horizontal Displacement ($m$)")
axes[0].set_ylabel("Vertical Displacement ($m$)")
axes[0].set_title("Trajectory of Projectile")
# v[:,1] represents the vertical velocity
axes[1].plot(t, v[:,1], label = "vy(t) vs. t")
axes[1].set_xlabel("Time ($s$)")
axes[1].set_ylabel("Vertical Velocity ($m$)")
axes[1].set_title("Vertical Velocity vs. Time")




