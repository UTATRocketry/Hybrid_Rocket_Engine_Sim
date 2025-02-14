import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import display, HTML

# Initialize figure
fig, ax = plt.subplots(figsize=(5, 8))
ax.set_xlim(-2, 2)  # Expand for fuel meters
ax.set_ylim(0, 10)  # Rocket launch height
ax.set_xlabel("Rocket")
ax.set_ylabel("Altitude")

# Rocket body
rocket, = ax.plot([], [], marker="^", markersize=20, color="red", label="Rocket")

# Exhaust trail
trail_x, trail_y = [], []
trail, = ax.plot([], [], "orange", lw=2, alpha=0.6)

# Fuel & Oxidizer meters (vertical bars)
fuel_meter, = ax.plot([1.6, 1.6], [0, 10], "g", lw=8, label="Fuel")  # Green bar
oxidizer_meter, = ax.plot([1.5, 1.5], [0, 10], "b", lw=8, label="Oxidizer")  # Blue bar

# Physics parameters
g = 9.81  # Gravity (m/s^2)
thrust = 30  # Thrust (m/s^2)
dt = 0.1  # Time step
velocity = 0  # Initial velocity
position = 0  # Initial position

# Fuel & Oxidizer levels
fuel = 10  # Start full (max height 10)
oxidizer = 10

def update(frame):
    global velocity, position, fuel, oxidizer, trail_x, trail_y

    # Newton's second law (Thrust - Gravity)
    acceleration = thrust - g
    velocity += acceleration * dt
    position += velocity * dt

    # Decrease fuel & oxidizer as rocket moves
    fuel = max(0, fuel - 0.1)  # Simulating fuel burn
    oxidizer = max(0, oxidizer - 0.12)  # Oxidizer burns slightly faster

    #Dynamically adjust Y-axis limits
    ax.set_ylim(0, max(10, position + 2))

    # Reset when rocket reaches the top
    if position > 100:
        ani.event_source.stop()

    # Update rocket position
    rocket.set_data(0, position)

    # Update exhaust trail
    trail_x.append(0)
    trail_y.append(position - 0.5)  # Slightly below rocket
    if len(trail_y) > 15:  # Keep trail short
        trail_x.pop(0)
        trail_y.pop(0)
    
    trail.set_data(trail_x, trail_y)

    # Update fuel & oxidizer meters
    fuel_meter.set_data([1.6, 1.6], [0, fuel])
    oxidizer_meter.set_data([1.5, 1.5], [0, oxidizer])

    return rocket, trail, fuel_meter, oxidizer_meter

# Create looping animation
ani = animation.FuncAnimation(fig, update, frames=100, interval=50, blit=False)
ani.save("rocket_animation.gif", writer="pillow")