# Relativistic Speed Simulator

A Three.js-based simulation that visualizes what a 3D object looks like when moving past an observer at relativistic speeds. The simulation includes:

- **Apparent position effects** (light-time delay)
- **Length contraction** (Lorentz contraction)
- **Doppler shift** (redshift/blueshift)
- **Geometric distortion** due to finite light speed

## Physics Background

### 1. Retarded Time
When an object moves at relativistic speeds, we don't see it where it actually is, but where it *was* when the light left it. For each point on the object, we solve:

```
|r(t') - r_observer| = c(t - t')
```

Where `t'` is the **retarded time** (when light left the point), and `t` is the observer's current time.

### 2. Apparent Position
Each vertex appears at its position at the retarded time:
```
r_apparent = r_0 + v * t'
```

Since different parts of the object have different retarded times, the object appears distorted.

### 3. Doppler Shift
The color of each vertex is shifted based on its velocity component along the line of sight:

```
f_observed = f_emitted / (γ(1 - β·cos(θ)))
```

Where:
- `γ = 1/√(1 - β²)` is the Lorentz factor
- `β = v/c` is the velocity as a fraction of light speed
- `θ` is the angle between velocity and line to observer

This causes:
- **Blueshift** (brighter, bluer) for parts moving toward you
- **Redshift** (dimmer, redder) for parts moving away

### 4. Length Contraction
Objects are contracted along their direction of motion by the Lorentz factor γ. This is automatically included in the apparent position calculations.

## Setup

### Quick Start

1. Open `index.html` in a modern web browser (Chrome, Firefox, Edge)
2. The simulation will start automatically with a default sphere model
3. Use the controls to adjust parameters

### Using Your Own 3D Model

#### From Polycam:

1. **Scan your subject** using the Polycam app on your phone
2. **Export the model**:
   - Open the scan in Polycam
   - Tap "Export"
   - Choose "GLB" format (recommended)
   - Quality: High
   - Download to your computer

3. **Add to project**:
   - Rename the file to `ilany.glb`
   - Drop it into the `model/` folder
   - Restart the simulation
   - The model will automatically load on startup!

4. **Alternative - Upload at runtime**:
   - Click "Upload Custom Model" button
   - Select any .glb or .gltf file
   - Great for testing different models

#### Model Requirements:
- Formats: `.glb` (preferred) or `.gltf`
- High-poly models (more vertices) show more detailed distortion
- Models are automatically centered and scaled to fit
- Falls back to procedural face if `model/ilany.glb` is not found

## Controls

### Velocity Slider
- Range: 0 to 0.999c (where c is the speed of light)
- Higher velocities show more dramatic effects
- The Lorentz factor (γ) is displayed below the slider

### Closest Pass Distance
- Range: 0.5 to 20 meters
- Distance of closest approach to the observer
- Closer passes show more extreme distortion

### Time Scale
- Range: 0.1x to 5x
- Controls animation speed
- Useful for studying effects in slow motion

### Reset Animation
- Restarts the animation from the beginning

## What to Look For

### At Low Speeds (< 0.3c)
- Minimal visual distortion
- Slight color shifts as object approaches/recedes

### At Moderate Speeds (0.3c - 0.8c)
- Noticeable compression along direction of motion
- Clear blueshift when approaching, redshift when receding
- Asymmetric appearance due to light-time effects

### At High Speeds (> 0.9c)
- Extreme distortion and warping
- Object appears "smeared" due to different retarded times
- Very pronounced color shifts
- The "trailing" side of the object is still visible from the front due to light-time delay

### Specific Effects to Notice:

1. **Headlight Effect**: The object appears brighter and bluer when approaching
2. **Searchlight Effect**: As the object passes, you briefly see parts of it that should be on the "back"
3. **Terrell Rotation**: At very high speeds, the object appears rotated rather than simply contracted
4. **Doppler Gradient**: Different parts of the object show different colors simultaneously

## Technical Details

### Coordinate System
- Observer is at the origin (0, 0, 0)
- Object moves along the x-axis
- Closest approach is at distance `d` along the y-axis

### Simulation vs Reality
- Speed of light is scaled down for visualization (C_VISUAL = 10 m/s in simulation)
- All relativistic effects are correctly scaled
- The actual speed of light is ~3×10⁸ m/s

### Performance
- Uses per-vertex calculations (CPU-based)
- Higher poly models may reduce framerate
- Optimized for models with 10k-100k vertices

## Troubleshooting

### Model doesn't appear
- Check browser console for errors
- Ensure model is in .glb or .gltf format
- Try the default sphere first to verify setup

### Simulation is slow
- Try a lower-poly model
- Reduce time scale
- Close other browser tabs

### Colors look wrong
- This is correct! Relativistic Doppler shift changes colors dramatically
- Try different velocities to see the range of effects

## Educational Use

This simulator is perfect for:
- Physics classes studying special relativity
- Visualizing length contraction and time dilation effects
- Understanding the difference between "seeing" and "measuring" in relativity
- Demonstrating that relativistic effects aren't just mathematical—they're visual!

## Further Reading

- Einstein's Special Relativity (1905)
- "The Terrell Effect" - why objects don't appear contracted the way you'd expect
- Relativistic ray tracing and rendering

## Credits

Built with Three.js and modern web standards.

Physics implementation based on:
- Lorentz transformations
- Retarded time calculations
- Relativistic Doppler formula
