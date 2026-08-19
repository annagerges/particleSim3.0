# Unified Multi-Physics Particle System

A C++ simulation modeling gravitational and spring forces on particles in 2D space, with spatial partitioning for efficient collision detection and RK4 numerical integration for physics accuracy. Python for data analysis.

---

## Overview

This project simulates a system of 100–1000 particles affected by gravity and spring forces, with particle-to-particle and wall collisions. The simulation uses:

- **RK4 integration** for <0.1% energy drift
- **Spatial partitioning grid** to reduce collision detection from O(n²) to O(n)
- **CSV export** for post-simulation analysis and validation
- **Python validation plots** to visualize energy conservation

---

## Physics Model

### Forces and Acceleration

**Gravitational Force (freefall)**:

```
a = -g = -9.8 m/s²
```

**Spring Force (when y ≤ spring height)**:

```
F_spring = -k(y - h)
a_total = -(k/m)(y - h) - g
acceleration is the force of the spring/mass minus gravity

The spring stiffness `k` is dynamically calculated to support all particles at equilibrium compression (0.2 m):
k = (nP × 9.8 × m / 0.2) × 4

Factor of 4 increases springiness and creates a responsive force.
```

### Numerical Integration

**RK4 (Runge-Kutta 4th Order)** with fixed timestep `dt = 0.1 s`. Computes four slope estimates (k1–k4) and combines them with weights (1:2:2:1) to approximate the solution:

```cpp
// RK4 update for a single particle's position and velocity
currentY = part[index].getY();
currentVy = part[index].getVy();

// First slope: current velocity and acceleration
k1Y = part[index].getVy();
k1Vy = checkAccel(currentY, s, part[index].getMass());

// Second slope: state at t + dt/2
testY2 = currentY + k1Y * (dt / 2);
testVy2 = currentVy + k1Vy * (dt / 2);
k2Y = testVy2;
k2Vy = checkAccel(testY2, s, part[index].getMass());

// Third slope: state at t + dt/2 (refined estimate)
testY3 = currentY + k2Y * (dt / 2);
testVy3 = currentVy + k2Vy * (dt / 2);
k3Y = testVy3;
k3Vy = checkAccel(testY3, s, part[index].getMass());

// Fourth slope: state at t + dt
testY4 = currentY + k3Y * dt;
testVy4 = currentVy + k3Vy * dt;
k4Y = testVy4;
k4Vy = checkAccel(testY4, s, part[index].getMass());

// Weighted average of slopes (midpoint estimates have weight 2)
finalY = currentY + (dt / 6) * (k1Y + 2*k2Y + 2*k3Y + k4Y);
finalVy = currentVy + (dt / 6) * (k1Vy + 2*k2Vy + 2*k3Vy + k4Vy);

```
Updates apply at (`0.01 s` intervals) within each frame for stability. Rk4 achieves <0.1% energy drift over long simulations compared to the approximately 2% drift with Euler integration.

## Program Architecture

### Spatial Partitioning

The simulation uses a dynamically sized to reduce collision detection from O(n²) to O(n) for typical particle distributions:
-**Grid Sizing**: Width and height are the maximum between 3 and  grid(sqrt(num of particles/10)+2)
- **Grid cell width**: 800 / cells per row
- **Collision checks**: Only particles in the same grid cell are tested
- **Grid update**: partial rebuild if any particles move cells via `clearAndFix()` to update the grid and manage resources effectively.

This optimization scales efficiently to 1000 particles without performance degradation.

### Collision Detection

#### Particle–Particle Collisions
- **Trigger**: Euclidean distance `d < 10` units (center to center. Radius of 5 per particle)
- **Axis determination**: Compare `dx` vs. `dy` to resolve collision along correct axis
	- **Response**:
	  - **Y-axis collision**: Reverse `vy`, apply damping (−0.1), boost lower particle if `vy < 30`
	  - **X-axis collision**: Reverse `vx` for both particles
	  - **If X and Y axis are equal**: apply Y-axis collision conditions.

#### Wall Collisions
- **Bounds**: `y values are from 6-max height and x values are from 0-max height (non-inclusive)`
- **Response**: Reverse normal velocity component, clamp position to boundary
- **Energy Injection**: Left wall adds `+0.1` to `vx` (if `vx < 30`) to counteract dampening and maintain motion. Prevents particles from stalling.

#### CSV Logging
Uses fstream C++ library to create a csv file and log the qualities of every particle (besides mechanical energy and mass (automatically 0.5kg to prevent division by 0 errors) for further analysis.

#### Python CSV Analysis
- **Method**: Reads CSV, computes kinetic, gravitational potential, and potential spring energy. Plots energy drift % over time.
- **Result**: RK4 method has <0.1% energy drift and demonstrates numerical stability

#### CSV Error Handling
Safely processes invalid data during analysis

```
for _, row in df.iterrows():
    # Convert particle ID; invalid entries become NaN
    raw\_pid = pd.to\_numeric(row['Particle num'], errors='coerce')
    
    # Only process valid particle IDs in target set
    if not pd.isna(raw\_pid) and int(raw\_pid) in target\_particles:
        p\_id = int(raw\_pid)
        px = pd.to\_numeric(row['x'], errors='coerce')
        py = pd.to\_numeric(row['y'], errors='coerce')
        
        # Skip invalid coordinates
        if not (pd.isna(px) or pd.isna(py)):
            particle\_paths[p\_id]['x'].append(px)
            particle\_paths[p\_id]['y'].append(py)
```
Prevents crashes and excludes invalid data.
	

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
| `clearAndFix()` | `void(vector<Particles>&, vector<vector<vector<Particles*>>>&, int width)` | Clears and puts any particle that moves cells in accordance to it's current position while maintaining O(1) time|
| `csvDump()` | `void csvDump(std::vector<Particles>&, std::fstream&);` | Logs particle state into csv file for futher analysis|

### Building and Running

#### Requirements
- **C++**: C++11 or later
- **Python**: 3.14+, with `pandas` and `matplotlib`


## Coming Soon!!
- **Damping Coefficient**: Make collision damping configurable; validate energy dissipation against expected mechanical loss.
- **Spring Stiffness Tuning**: Expose `k` scaling factor as command-line parameter.
- **And More!!**


---

**Date**: August 2026  
**License**: MIT
