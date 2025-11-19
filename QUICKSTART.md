# Quick Start Guide

## Running the Simulator

### Option 1: Simple Double-Click (Recommended)
1. Simply open `index.html` in your web browser
2. The simulation will start immediately with a default model

### Option 2: Local Server (Better for development)
Run a local server to avoid CORS issues with model loading:

```bash
# Using Python 3 (built into most systems)
python3 -m http.server 8000

# Or using npm
npm run dev
```

Then open: `http://localhost:8000`

## Adding Your Teacher's 3D Model from Polycam

### Step 1: Capture the Scan
1. Open **Polycam** on your phone
2. Create a new scan (Photo mode recommended for faces)
3. Take 20-50 photos around your teacher's head from different angles
4. Let Polycam process the 3D model

### Step 2: Export from Polycam
1. Open your completed scan
2. Tap the **Share** icon
3. Select **Export**
4. Choose format: **GLB** (recommended - single file)
5. Quality: Choose **High** for more detail
6. Download to your computer

### Step 3: Add to Project
1. Rename the exported file to `ilany.glb`
2. Drop it into the `model/` folder in this project
3. Run the simulator (open `index.html`)
4. Watch your teacher's face get relativistically distorted!

**Note:** If no model is found, a procedural face will be used as fallback. You can also use the "Upload Custom Model" button to try different models without replacing the file.

## First-Time Experience

### Try These Settings:

**Subtle Effects (Start Here):**
- Velocity: 0.5c
- Distance: 5m
- Watch the face compress and shift colors

**Moderate Effects:**
- Velocity: 0.85c
- Distance: 3m
- Notice the asymmetric distortion

**Extreme Effects:**
- Velocity: 0.99c
- Distance: 1m
- See extreme warping and the "Terrell rotation" effect

## What You'll See

### As the face approaches:
- ✨ **Blueshift**: Face gets bluer and brighter
- 🔷 **Compression**: Face appears squished along motion direction
- 👁️ **Weird angles**: You can see parts of the back of the head

### As the face passes by:
- 🌈 **Color gradient**: Front is blue, back is red
- 🌀 **Maximum distortion**: Most warped at closest point

### As the face recedes:
- 🔴 **Redshift**: Face gets redder and dimmer
- 📏 **Stretching**: Face appears elongated
- 👻 **Ghost effect**: Trailing edge is still visible

## Tips for Best Results

### For Cool Screenshots:
1. Set velocity to 0.95c
2. Set distance to 2m
3. Pause at interesting moments (adjust time scale to 0.1x)
4. Use your browser's screenshot tool

### For Understanding the Physics:
1. Start at 0.3c and slowly increase
2. Pay attention to how the distortion changes
3. Notice that it's NOT just simple squishing!

### Model Quality:
- **More vertices = More detail** in the distortion
- Polycam's high-quality exports work best
- Lower quality is faster but less impressive

## Troubleshooting

**Model appears inside-out:**
- This can happen with some exports
- Try re-exporting from Polycam

**Model is too small/large:**
- The simulator auto-scales, but if it looks wrong
- Try re-centering in Polycam before export

**Performance is slow:**
- Try reducing model quality in Polycam
- Close other tabs/applications
- Use time scale slider to slow down animation

**Colors look weird:**
- This is CORRECT! Relativistic Doppler shift is dramatic
- Real objects at these speeds would look like this

## Making It Your Own

### Modify the trajectory:
Edit `main.js`, line ~195:
```javascript
const x0 = -50; // Starting position
const y0 = this.closestDistance; // Height of trajectory
```

### Change background:
Edit `main.js`, line ~51:
```javascript
this.renderer.setClearColor(0x000510); // Dark blue space
```

### Adjust visual speed of light:
Edit `main.js`, line ~8:
```javascript
const C_VISUAL = 10; // Higher = faster animation
```

## Next Steps

1. Try different objects (not just faces!)
2. Experiment with different camera angles (drag to rotate view)
3. Screenshot at different velocities and compare
4. Show your physics teacher 😄

Enjoy warping reality!
