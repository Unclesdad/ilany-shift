import * as THREE from 'three';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// Physical constants
const C = 299792458; // Speed of light in m/s
const C_VISUAL = 10; // Visual speed of light for simulation (m/s in simulation units)

class RelativisticSimulator {
    constructor() {
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true });

        // Simulation parameters
        this.velocity = 0.9; // Fraction of c
        this.closestDistance = 5; // Meters
        this.timeScale = 1.0;
        this.currentTime = -30; // Start time (negative so object approaches from distance)
        this.animationRunning = true;

        // Model
        this.originalGeometry = null;
        this.relativisticMesh = null;
        this.defaultModelLoaded = false;

        this.init();
        this.setupControls();
        this.loadDefaultModel();
        this.animate();
    }

    init() {
        // Setup renderer
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setClearColor(0x000510);
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);

        // Setup camera
        this.camera.position.set(0, 0, 0);
        this.camera.lookAt(0, 0, -1);

        // Add orbit controls for better viewing
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        // Add lighting
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        this.scene.add(ambientLight);

        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(5, 5, 5);
        this.scene.add(directionalLight);

        // Add reference grid
        const gridHelper = new THREE.GridHelper(50, 50, 0x444444, 0x222222);
        gridHelper.position.y = -3;
        this.scene.add(gridHelper);

        // Add axes helper
        const axesHelper = new THREE.AxesHelper(5);
        this.scene.add(axesHelper);

        // Handle window resize
        window.addEventListener('resize', () => this.onWindowResize());

        document.getElementById('loading').style.display = 'none';
    }

    setupControls() {
        // Velocity control
        const velocitySlider = document.getElementById('velocity');
        const velocityValue = document.getElementById('velocity-value');
        const gammaValue = document.getElementById('gamma-value');

        velocitySlider.addEventListener('input', (e) => {
            this.velocity = parseFloat(e.target.value);
            const gamma = this.getLorentzFactor(this.velocity);
            velocityValue.textContent = `${this.velocity.toFixed(3)}c`;
            gammaValue.textContent = `γ = ${gamma.toFixed(3)}`;
        });

        // Distance control
        const distanceSlider = document.getElementById('distance');
        const distanceValue = document.getElementById('distance-value');

        distanceSlider.addEventListener('input', (e) => {
            this.closestDistance = parseFloat(e.target.value);
            distanceValue.textContent = `${this.closestDistance.toFixed(1)} m`;
        });

        // Time scale control
        const timeScaleSlider = document.getElementById('time-scale');
        const timeScaleValue = document.getElementById('time-scale-value');

        timeScaleSlider.addEventListener('input', (e) => {
            this.timeScale = parseFloat(e.target.value);
            timeScaleValue.textContent = `${this.timeScale.toFixed(1)}x`;
        });

        // Reset button
        document.getElementById('reset').addEventListener('click', () => {
            this.currentTime = -30;
            this.animationRunning = true;
        });

        // Upload button
        document.getElementById('upload').addEventListener('click', () => {
            document.getElementById('file-input').click();
        });

        document.getElementById('file-input').addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (file) {
                this.loadModelFromFile(file);
            }
        });

        // Initialize displays
        velocitySlider.dispatchEvent(new Event('input'));
        distanceSlider.dispatchEvent(new Event('input'));
        timeScaleSlider.dispatchEvent(new Event('input'));
    }

    getLorentzFactor(beta) {
        return 1.0 / Math.sqrt(1.0 - beta * beta);
    }

    // Calculate the retarded time for a point
    // Given observer time t, find t' such that the light travel time matches
    calculateRetardedTime(t, vertexPos, observerPos) {
        const beta = this.velocity;

        // Actual position of object: x(t') = x0 + v*t'
        // For object moving along x-axis, starting far left
        // Distance from observer: |r(t') - r_obs| = c * (t - t')

        // Object trajectory: moves along x-axis at height y = closestDistance
        // x(t') = -50 + v*C_VISUAL*t' (starts at x = -50)
        // y(t') = closestDistance
        // z(t') = 0

        // Observer at origin (0, 0, 0)

        // We need to solve: |r(t') - r_obs| = C_VISUAL * (t - t')
        // This is a quadratic equation in t'

        const v = beta * C_VISUAL;
        const x0 = -50; // Starting position

        // Vertex position in object frame (relative to object center)
        const vx = vertexPos.x;
        const vy = vertexPos.y;
        const vz = vertexPos.z;

        // Full 3D calculation
        // x(t') = x0 + v*t' + vx
        // y(t') = closestDistance + vy
        // z(t') = vz

        const y0 = this.closestDistance;

        // Distance squared: (x0 + v*t' + vx)^2 + (y0 + vy)^2 + vz^2 = C_VISUAL^2 * (t - t')^2

        // Expand: (x0 + vx)^2 + 2(x0 + vx)*v*t' + v^2*t'^2 + (y0 + vy)^2 + vz^2 = C_VISUAL^2 * (t^2 - 2*t*t' + t'^2)

        // Rearrange: (v^2 - C_VISUAL^2)*t'^2 + (2*v*(x0 + vx) + 2*C_VISUAL^2*t)*t' + ((x0 + vx)^2 + (y0 + vy)^2 + vz^2 - C_VISUAL^2*t^2) = 0

        const a = v * v - C_VISUAL * C_VISUAL;
        const b = 2 * v * (x0 + vx) + 2 * C_VISUAL * C_VISUAL * t;
        const c = (x0 + vx) * (x0 + vx) + (y0 + vy) * (y0 + vy) + vz * vz - C_VISUAL * C_VISUAL * t * t;

        // Solve quadratic equation
        const discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return t; // No solution, use current time
        }

        const t1 = (-b - Math.sqrt(discriminant)) / (2 * a);
        const t2 = (-b + Math.sqrt(discriminant)) / (2 * a);

        // Choose the retarded time (past time, t' < t)
        return Math.min(t1, t2);
    }

    // Calculate apparent position of vertex at observer time t
    getApparentPosition(t, vertexPosLocal) {
        const tRetarded = this.calculateRetardedTime(t, vertexPosLocal, new THREE.Vector3(0, 0, 0));

        // Position at retarded time
        const v = this.velocity * C_VISUAL;
        const x0 = -50;

        const x = x0 + v * tRetarded + vertexPosLocal.x;
        const y = this.closestDistance + vertexPosLocal.y;
        const z = vertexPosLocal.z;

        return new THREE.Vector3(x, y, z);
    }

    // Calculate Doppler factor for a vertex
    getDopplerFactor(t, vertexPosLocal) {
        const tRetarded = this.calculateRetardedTime(t, vertexPosLocal, new THREE.Vector3(0, 0, 0));

        const v = this.velocity * C_VISUAL;
        const x0 = -50;

        const x = x0 + v * tRetarded + vertexPosLocal.x;
        const y = this.closestDistance + vertexPosLocal.y;
        const z = vertexPosLocal.z;

        // Vector from observer to vertex
        const r = new THREE.Vector3(x, y, z);
        const rMag = r.length();

        if (rMag < 0.001) return 1.0;

        // Velocity vector (along x-axis)
        const vel = new THREE.Vector3(v, 0, 0);

        // Component of velocity along line of sight (radial velocity)
        const vRadial = vel.dot(r) / rMag;

        // Doppler factor: f_observed / f_emitted = sqrt((1 - beta) / (1 + beta)) for recession
        // More generally: f_obs / f_emit = 1 / (gamma * (1 + beta * cos(theta)))
        // where theta is angle between velocity and line to observer

        const beta = this.velocity;
        const gamma = this.getLorentzFactor(beta);

        // cos(theta) = -vRadial / (v * 1) (negative because we want angle from velocity to observer)
        const cosTheta = -vRadial / (v * 1.0 + 0.0001);

        // Relativistic Doppler factor
        const dopplerFactor = 1.0 / (gamma * (1.0 - beta * cosTheta));

        return dopplerFactor;
    }

    // Apply Doppler shift to color
    applyDopplerShift(originalColor, dopplerFactor) {
        // dopplerFactor > 1: blueshift (approaching)
        // dopplerFactor < 1: redshift (receding)

        // Convert color to HSL and shift hue
        const color = new THREE.Color(originalColor);
        const hsl = { h: 0, s: 0, l: 0 };
        color.getHSL(hsl);

        // Map Doppler factor to hue shift
        // Factor of 2 = full blueshift, factor of 0.5 = full redshift
        const logFactor = Math.log2(dopplerFactor);
        const hueShift = logFactor * 0.15; // Adjust sensitivity

        // Shift hue (0.66 = blue, 0 = red in HSL)
        let newHue = hsl.h + hueShift;
        newHue = Math.max(0, Math.min(1, newHue));

        // Also adjust brightness based on Doppler factor
        const newL = hsl.l * Math.pow(dopplerFactor, 0.5);

        color.setHSL(newHue, hsl.s, Math.max(0.1, Math.min(0.9, newL)));

        return color;
    }

    updateRelativisticMesh() {
        if (!this.relativisticMesh || !this.originalGeometry) return;

        const t = this.currentTime;
        const positions = this.relativisticMesh.geometry.attributes.position.array;
        const colors = this.relativisticMesh.geometry.attributes.color.array;
        const originalPositions = this.originalGeometry.attributes.position.array;
        const originalColors = this.originalGeometry.attributes.color?.array;

        // Update each vertex
        for (let i = 0; i < originalPositions.length; i += 3) {
            const originalPos = new THREE.Vector3(
                originalPositions[i],
                originalPositions[i + 1],
                originalPositions[i + 2]
            );

            // Get apparent position
            const apparentPos = this.getApparentPosition(t, originalPos);

            positions[i] = apparentPos.x;
            positions[i + 1] = apparentPos.y;
            positions[i + 2] = apparentPos.z;

            // Calculate Doppler shift
            const dopplerFactor = this.getDopplerFactor(t, originalPos);

            // Get base color from original geometry or use gray default
            let baseColor;
            if (originalColors) {
                baseColor = new THREE.Color(
                    originalColors[i],
                    originalColors[i + 1],
                    originalColors[i + 2]
                );
            } else {
                baseColor = new THREE.Color(0xcccccc);
            }

            const shiftedColor = this.applyDopplerShift(baseColor, dopplerFactor);

            colors[i] = shiftedColor.r;
            colors[i + 1] = shiftedColor.g;
            colors[i + 2] = shiftedColor.b;
        }

        this.relativisticMesh.geometry.attributes.position.needsUpdate = true;
        this.relativisticMesh.geometry.attributes.color.needsUpdate = true;
        this.relativisticMesh.geometry.computeVertexNormals();
    }

    loadDefaultModel() {
        // Create a detailed procedural face model
        const geometry = this.createProceduralFace();
        this.setupRelativisticMesh(geometry);
        this.defaultModelLoaded = true;
    }

    createProceduralFace() {
        // Start with a high-resolution sphere as base
        const baseGeometry = new THREE.SphereGeometry(1.5, 128, 128);
        const positions = baseGeometry.attributes.position.array;
        const colors = new Float32Array(positions.length);

        // Define face features in local coordinates
        // Front of face is in +z direction

        for (let i = 0; i < positions.length; i += 3) {
            let x = positions[i];
            let y = positions[i + 1];
            let z = positions[i + 2];

            // Normalize to get direction
            const len = Math.sqrt(x*x + y*y + z*z);
            const nx = x / len;
            const ny = y / len;
            const nz = z / len;

            // Make head shape more elongated (ellipsoid)
            y *= 1.2; // Taller head
            x *= 0.9; // Narrower head
            z *= 1.0;

            // Flatten the back of the head
            if (nz < -0.3) {
                z *= 0.7;
            }

            // Add eyes (two spherical bulges)
            const eyeLeftCenter = new THREE.Vector3(-0.5, 0.4, 1.2);
            const eyeRightCenter = new THREE.Vector3(0.5, 0.4, 1.2);
            const eyeRadius = 0.35;

            const distToLeftEye = Math.sqrt(
                Math.pow(x - eyeLeftCenter.x, 2) +
                Math.pow(y - eyeLeftCenter.y, 2) +
                Math.pow(z - eyeLeftCenter.z, 2)
            );

            const distToRightEye = Math.sqrt(
                Math.pow(x - eyeRightCenter.x, 2) +
                Math.pow(y - eyeRightCenter.y, 2) +
                Math.pow(z - eyeRightCenter.z, 2)
            );

            // Eye bulges
            if (distToLeftEye < eyeRadius && nz > 0.5) {
                const bulge = (1 - distToLeftEye / eyeRadius) * 0.3;
                z += bulge;
                // Pupil (darker color)
                if (distToLeftEye < 0.15) {
                    colors[i] = 0.1;
                    colors[i + 1] = 0.1;
                    colors[i + 2] = 0.1;
                } else {
                    colors[i] = 0.9;
                    colors[i + 1] = 0.9;
                    colors[i + 2] = 0.9;
                }
            } else if (distToRightEye < eyeRadius && nz > 0.5) {
                const bulge = (1 - distToRightEye / eyeRadius) * 0.3;
                z += bulge;
                // Pupil (darker color)
                if (distToRightEye < 0.15) {
                    colors[i] = 0.1;
                    colors[i + 1] = 0.1;
                    colors[i + 2] = 0.1;
                } else {
                    colors[i] = 0.9;
                    colors[i + 1] = 0.9;
                    colors[i + 2] = 0.9;
                }
            } else {
                // Skin tone base
                colors[i] = 0.9;
                colors[i + 1] = 0.75;
                colors[i + 2] = 0.65;
            }

            // Add nose (triangular protrusion in center)
            if (Math.abs(x) < 0.3 && y > -0.2 && y < 0.4 && nz > 0.6) {
                const noseFactor = (0.3 - Math.abs(x)) / 0.3;
                const heightFactor = 1 - Math.abs(y - 0.1) / 0.5;
                const noseBulge = noseFactor * heightFactor * 0.4;
                z += noseBulge;

                // Nostrils (darker)
                if (y < 0.0 && Math.abs(x) > 0.1 && Math.abs(x) < 0.25) {
                    colors[i] *= 0.7;
                    colors[i + 1] *= 0.7;
                    colors[i + 2] *= 0.7;
                }
            }

            // Add mouth (indentation)
            if (Math.abs(x) < 0.5 && y > -0.8 && y < -0.3 && nz > 0.5) {
                const mouthWidth = Math.pow(1 - Math.abs(x) / 0.5, 2);
                const mouthHeight = 1 - Math.abs(y + 0.55) / 0.25;
                const mouthDepth = mouthWidth * mouthHeight * 0.15;
                z -= mouthDepth;

                // Lips (slightly reddish)
                if (Math.abs(y + 0.55) < 0.1) {
                    colors[i] = 0.8;
                    colors[i + 1] = 0.4;
                    colors[i + 2] = 0.4;
                }
            }

            // Add ears
            const earLeftCenter = new THREE.Vector3(-1.2, 0.0, 0.0);
            const earRightCenter = new THREE.Vector3(1.2, 0.0, 0.0);

            if (Math.abs(y) < 0.6) {
                // Left ear
                if (x < -0.9 && nz > -0.5 && nz < 0.5) {
                    const earBulge = (1.2 - Math.abs(x)) * 0.3;
                    x -= earBulge;
                    z += earBulge * 0.2 * Math.sin(y * 2);
                }
                // Right ear
                if (x > 0.9 && nz > -0.5 && nz < 0.5) {
                    const earBulge = (1.2 - Math.abs(x)) * 0.3;
                    x += earBulge;
                    z += earBulge * 0.2 * Math.sin(y * 2);
                }
            }

            // Add some eyebrows (color only)
            if (y > 0.6 && y < 0.8 && nz > 0.7) {
                if ((x > 0.2 && x < 0.7) || (x < -0.2 && x > -0.7)) {
                    colors[i] *= 0.5;
                    colors[i + 1] *= 0.4;
                    colors[i + 2] *= 0.3;
                }
            }

            // Add some hair texture on top
            if (y > 0.8 && nz > -0.5) {
                const hairNoise = Math.sin(x * 15) * Math.cos(z * 15) * 0.05;
                y += Math.abs(hairNoise);
                // Hair color
                colors[i] = 0.3;
                colors[i + 1] = 0.2;
                colors[i + 2] = 0.15;
            }

            // Add facial structure (cheekbones)
            if (Math.abs(x) > 0.5 && Math.abs(x) < 1.0 && y > -0.2 && y < 0.3 && nz > 0.3) {
                const cheekFactor = (Math.abs(x) - 0.5) / 0.5;
                z += cheekFactor * 0.15;
            }

            // Chin
            if (Math.abs(x) < 0.4 && y < -0.7 && nz > 0.4) {
                const chinBulge = (1 - Math.abs(x) / 0.4) * 0.2;
                z += chinBulge;
            }

            // Update positions
            positions[i] = x;
            positions[i + 1] = y;
            positions[i + 2] = z;
        }

        // Set the colors attribute
        baseGeometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        baseGeometry.computeVertexNormals();

        return baseGeometry;
    }

    loadModelFromFile(file) {
        const loader = new GLTFLoader();
        const url = URL.createObjectURL(file);

        document.getElementById('loading').style.display = 'block';

        loader.load(
            url,
            (gltf) => {
                // Extract geometry from loaded model
                let geometry = null;

                gltf.scene.traverse((child) => {
                    if (child.isMesh && !geometry) {
                        geometry = child.geometry.clone();
                    }
                });

                if (geometry) {
                    // Center and scale the geometry
                    geometry.center();
                    geometry.computeBoundingSphere();
                    const scale = 2.0 / geometry.boundingSphere.radius;
                    geometry.scale(scale, scale, scale);

                    this.setupRelativisticMesh(geometry);
                    document.getElementById('loading').style.display = 'none';
                } else {
                    alert('No mesh found in the model!');
                    document.getElementById('loading').style.display = 'none';
                }

                URL.revokeObjectURL(url);
            },
            (progress) => {
                console.log('Loading:', (progress.loaded / progress.total * 100).toFixed(2) + '%');
            },
            (error) => {
                console.error('Error loading model:', error);
                alert('Error loading model: ' + error.message);
                document.getElementById('loading').style.display = 'none';
            }
        );
    }

    setupRelativisticMesh(geometry) {
        // Remove old mesh
        if (this.relativisticMesh) {
            this.scene.remove(this.relativisticMesh);
            this.relativisticMesh.geometry.dispose();
            this.relativisticMesh.material.dispose();
        }

        // Store original geometry
        this.originalGeometry = geometry.clone();

        // Create geometry with vertex colors
        const newGeometry = geometry.clone();
        const colors = new Float32Array(newGeometry.attributes.position.count * 3);
        for (let i = 0; i < colors.length; i += 3) {
            colors[i] = 0.8;
            colors[i + 1] = 0.8;
            colors[i + 2] = 0.8;
        }
        newGeometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

        // Create material that uses vertex colors
        const material = new THREE.MeshPhongMaterial({
            vertexColors: true,
            side: THREE.DoubleSide,
            flatShading: false
        });

        this.relativisticMesh = new THREE.Mesh(newGeometry, material);
        this.scene.add(this.relativisticMesh);

        // Add wireframe for better visualization
        const wireframe = new THREE.WireframeGeometry(newGeometry);
        const line = new THREE.LineSegments(wireframe);
        line.material.opacity = 0.1;
        line.material.transparent = true;
        line.material.color = new THREE.Color(0x00ffff);
        this.relativisticMesh.add(line);
    }

    updateInfo() {
        document.getElementById('time-display').textContent = this.currentTime.toFixed(2);

        const v = this.velocity * C_VISUAL;
        const x = -50 + v * this.currentTime;
        document.getElementById('position-display').textContent = x.toFixed(2);
        document.getElementById('rel-velocity').textContent = this.velocity.toFixed(3);
    }

    onWindowResize() {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }

    animate() {
        requestAnimationFrame(() => this.animate());

        if (this.animationRunning) {
            this.currentTime += 0.016 * this.timeScale; // Approximately 60 FPS

            // Reset if object has passed far away
            if (this.currentTime > 30) {
                this.currentTime = -30;
            }

            this.updateRelativisticMesh();
            this.updateInfo();
        }

        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
}

// Start the simulation
const simulator = new RelativisticSimulator();
