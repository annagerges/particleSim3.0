# Particle System

## Overview
A particle simulation engine that uses gravitational dynamics, spring forces, and collision detection using spatial partitioning for efficient computation. The system currently models **100–1000 particles** interacting on a 2D grid with energy conservation validation and adaptive timestep processing.

## Physics Model

### Forces and Acceleration

**Gravitational Force: for particles in freefall or flying upwards**:
```
a = -g = -9.8 m/s²
```

**Spring Force (when y ≤ spring height)**):
```
F_spring = -k(y - h)
a_total = -(k/m)(y - h) - g
acceleration is the force of the spring/mass minus gravity
```

where `k` is dynamically calculated based on particle count and assumed compression (0.2 m), scaled by factor of 4 for springiness:
```
k = (nP × 9.8 × m / 0.2) × 4
Multiplied by 4 for strength and springiness
```

### Numerical Integration

**RK4(Runge-Kutta 4th order) Approximation** with fixed timestep `dt = 0.1 s`:
```
		currentY = part[index].getY();
		currentVy = part[index].getVy();

		//Slope of postition and velosity respectively
		k1Y = part[index].getVy();
		k1Vy = checkAccel(currentY,s,part[index].getMass());

		//predicted y and velocity in 0.05 seconds
		testY2 = currentY + k1Y * (dt / 2);
		testVy2 = currentVy + k1Vy * (dt / 2);

		//Slope of postition and velosity respectively
		k2Y = testVy2;
		k2Vy= checkAccel(testY2, s, part[index].getMass());

		//predicted y and velocity in 0.05 seconds
		testY3 = testY2 + k2Y * (dt/2);
		testVy3 = testVy2 + k2Vy * (dt/2);

		//Slope of postition and velosity respectively
		k3Y = testVy3;
		k3Vy= checkAccel(testY3, s, part[index].getMass());

		//predicted y and velocity in 0.1 seconds
		testY4 = testY3 + k3Y * dt;
		testVy4 = testVy3 + k3Vy * dt;

		k4Y = testVy4;
		k4Vy = checkAccel(testY4, s, part[index].getMass());

		//weighted avg of y and vy to find final values (the midpoint values are twice as accurate because they provide better estimate of avg slope)
		finalY = currentY + (dt / 6) * (k1Y+(2*k2Y)+(2*k3Y)+k4Y);
		finalVy = currentVy + (dt / 6) * (k1Vy + (2 * k2Vy) + (2 * k3Vy) + k4Vy);
```

Updates apply at (`0.01 s` intervals) within each frame for stability.

## Program Architecture

### Spatial Partitioning

The simulation uses a dynamically sized grid(sqrt(num of particles/10)+2) of equal dimensions to reduce collision detection from O(n²) to O(n) for typical particle distributions:
- **Grid cell width**: 800 / cells per row
- **Collision checks**: Only particles in the same grid cell are tested
- **Grid update**: partial rebuild if any particles move cells via `clearAndFix()` to update the grid and manage resources effectively.

This optimization scales efficiently to 1000 particles without performance degradation.

### Collision Detection

#### Particle–Particle Collisions
- **Trigger**: Euclidean distance `d < 10` units (collision radius = 5 per particle)
- **Axis determination**: Compare `dx` vs. `dy` to resolve collision along correct axis
- **Response**:
  - **Y-axis collision**: Reverse `vy`, apply damping (−0.1), boost lower particle if `vy < 30`
  - **X-axis collision**: Reverse `vx` for both particles
  - **If X and Y axis are equal**: apply Y-axis collision conditions.

#### Wall Collisions
- **Bounds**: `[1, 799]` × `[1, 799]`
- **Response**: Reverse normal velocity component, clamp position to boundary
- **Energy boost**: Left wall adds `+0.1` to `vx` (if `vx < 30`) to maintain motion without overloading speed.

#### CSV Logging
Uses fstream C++ library to create a csv file and log the qualities of every particle (besides mechanical energy and mass (automatically 0.5kg to prevent division by 0 errors) for further analysis.

#### Python CSV Analysis
- **Libraries**: Uses matplotlib.pyplot and pandas to create particle trajectory graphs and rk4 energy conservation validation.
- **Error Handling**: errors='coerce' to safely convert invalid data points to NaN(not a number) without breaking the program and using isna() and continue to exclude any invalid points. Trajectory program prevents users from picking the same particle twice and picking an invalid amount of data points. It only graphs lines with valid data points.

```
for _, row in df.iterrows():
    #convert the particle id to a number and safely make it NaN if unsuccessful
    raw_pid=pd.to_numeric(row['Particle num'],errors='coerce')

    ##if the particle number is a number and is one of the target particles convert the number to an integer, the coordinates to number (with NaN error handling), and append it to the dictionary if the coordinates are numbers
    if not pd.isna(raw_pid) and int(raw_pid) in target_particles:
        p_id = int(raw_pid)
        px=pd.to_numeric(row['x'], errors='coerce')
        py=pd.to_numeric(row['y'], errors='coerce')

        if not (pd.isna(px) or pd.isna(py)):
            particle_paths[p_id]['x'].append(px)
            particle_paths[p_id]['y'].append(py)

```

- **Energy Conservation Validation Graph**: Graphs a line of ideal energy conservation (0% energy drift) compared to rk4 approximation energy drift % with time (Euler approximation line coming soon).
	

### Frame Timing

Uses `std::chrono::high_resolution_clock` to decouple frame rate from simulation timestep:
```
accumulator += frameDuration.count()
while (accumulator >= dt) {
    update()
    accumulator -= dt
}
```

Ensures consistent physics independent of frame rate or system load.

## Code Structure

| File | Purpose |
|------|---------|
| `particleSim2.0.cpp` | Initialization, user input validation, grid setup, main loop |
| `Particles.cpp` | Physics update, collision detection, grid management |
| `Particles.h` | Class definitions (`Particles`, `Spring`), function declarations |

### Key Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `updatePos()` | `void(vector<Particles>&, Spring&)` | Apply forces, update velocity and position |
| `wallCollis()` | `void(vector<Particles>&)` | Handle boundary collisions |
| `particleCollis()` | `void(vector<Particles*>&)` | Detect and resolve particle–particle collisions within grid cells with more than 1 particle |
| `clearAndFix()` | `void(vector<Particles>&, vector<vector<vector<Particles*>>>&, int width)` | Clears and puts any particle that moves cells in accordance to it's current position|
| `csvDump()` | `void csvDump(std::vector<Particles>&, std::fstream&);` | Logs particle state into csv file for futher analysis|


## Coming Soon!!
- **Rendering**: Visualize particles and grid using OpenGL or SDL2.
- **Damping Coefficient**: Make collision damping configurable; validate energy dissipation against expected mechanical loss.
- **Spring Stiffness Tuning**: Expose `k` scaling factor as command-line parameter.


---

**Date**: August 2026  
**License**: MIT
