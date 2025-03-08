import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import display, HTML
import rocketCEA
import math
dry_mass = 40
# Physics parameters
g = 9.81  # Gravity (m/s^2)
dt = 0.01  # Time step
velocity = 0  # Initial velocity
y_position = 0  # Initial position
x_position = 0
max_height = 100
# Exhaust trail
trail_x, trail_y = [], []
def anim(new_system):
    # Initialize figure
    fig, ax = plt.subplots(figsize=(5, 8))
    ax.set_xlim(-2, 2)  # Expand for fuel meters
    ax.set_ylim(0, 10)  # Rocket launch height
    ax.set_xlabel("Rocket")
    ax.set_ylabel("Altitude")
    trail, = ax.plot([], [], "orange", lw=2, alpha=0.6)

    # Rocket body
    rocket, = ax.plot([], [], marker="^", markersize=20, color="red", label="Rocket")
    marker, = ax.plot([0, 0], [10, 10], marker="X", markersize=20, color="black", label="max_height")



    # Fuel & Oxidizer meters (vertical bars)
    fuel_meter, = ax.plot([1.6, 1.6], [0, 10], "g", lw=8, label="Fuel")  # Green bar
    oxidizer_meter, = ax.plot([1.5, 1.5], [0, 10], "b", lw=8, label="Oxidizer")  # Blue bar

    plt.legend(loc = "lower left")




    #new_system = rocketCEA.on_button_click()
    
    title = ax.text(0.5,0.85, y_position, bbox={'facecolor':'w', 'alpha':0.5, 'pad':5}, transform=ax.transAxes, ha="center")
    # Fuel & Oxidizer levels

    wind_speed = 10
    def update(frame):
        global velocity, y_position, fuel, oxidizer, trail_x, trail_y, x_position, max_height
        fuel = new_system['Fuel_mass']
        oxidizer = new_system['Ox_Mass']
        thrust =new_system['Thrust']  # Thrust (m/s^2)

        if frame<len(thrust):
            thrust_acc = thrust[frame]/(dry_mass+oxidizer[frame]+fuel[frame])
        else:
            thrust_acc = 0
        # Newton's second law (Thrust - Gravity)
        acceleration = thrust_acc - g
        velocity += acceleration * dt
        y_position += velocity * dt
        x_position += wind_speed * dt
        # Decrease fuel & oxidizer as rocket moves
        if frame<len(fuel):
            fuel_new = fuel[frame]  # Simulating fuel burn
        else:
            fuel_new = 0

        if frame<len(oxidizer):
            oxidizer_new = oxidizer[frame]  # Oxidizer burns slightly faster
        else:
            oxidizer_new = 0
        #Dynamically adjust Y-axis limits
        ax.set_xlim(x_position - 10, x_position + 10)
        ax.set_ylim(max(0, y_position - 10), max(10, y_position + 2))


        # Update rocket position
        rocket.set_data(x_position, y_position)
        title.set_text("Current Position: {}. Max position: {}".format(round(y_position, 1), max_height))
        marker.set_data(x_position, min(max_height, max(10, y_position + 2)))
        # Update exhaust trail
        trail_x.append(x_position)
        trail_y.append(y_position - 0.5)  # Slightly below rocket
        
        trail.set_data(trail_x, trail_y)

        # Update fuel & oxidizer meters
        fuel_meter.set_data([max(10, x_position+ 10)-0.5, max(10, x_position+ 10)-0.5], [max(0, y_position - 10), fuel_new+max(0, y_position - 10)])
        oxidizer_meter.set_data([max(10, x_position+10)-0.1, max(10, x_position+10)-0.1], [max(0, y_position - 10), oxidizer_new+max(0, y_position - 10)])

        #if(y_position > max_height):
        #    ani.event_source.stop()

        return rocket, trail, fuel_meter, oxidizer_meter

    # Create looping animation
    ani = animation.FuncAnimation(fig, update, frames=400, interval=50, blit=False)
    ani.save("rocket_animation.gif", writer="pillow")
    print("Gif Updated!")
