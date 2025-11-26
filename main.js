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
        this.staticRotation = Math.PI / 2; // Static rotation angle in radians (90 degrees)
        this.angularVelocity = 0; // rad/s (positive = CW, negative = CCW)
        this.spinningEnabled = false; // Whether spinning mode is enabled
        this.timeScale = 1.0;
        this.currentTime = -10; // Start time (negative so object approaches from distance)
        this.animationRunning = true;
        this.manualTimeControl = false; // Track if user is manually controlling time

        // Model
        this.originalGeometry = null;
        this.relativisticMesh = null;
        this.actualPositionMesh = null; // Transparent mesh showing actual position
        this.defaultModelLoaded = false;

        // Display options
        this.showActualPosition = false;

        this.init();
        this.setupControls();
        this.loadDefaultModel();
        this.animate();
    }

    init() {
        // Setup renderer
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setClearColor(0xffffff); // White background
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);

        // Setup camera
        // Observer effectively at origin, rotated 180 degrees to look backward initially
        // Position camera very slightly offset so OrbitControls can rotate
        this.camera.position.set(0, 0, -0.01);
        this.camera.lookAt(0, 0, 0);

        // Add orbit controls for rotation (camera orbits very close to origin)
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;
        this.controls.target.set(0, 0, 0); // Target at origin
        this.controls.enableZoom = false; // No zooming
        this.controls.enablePan = false; // No panning
        this.controls.enableRotate = true; // Allow rotation
        this.controls.minDistance = 0.01; // Lock distance at 1cm (effectively at origin)
        this.controls.maxDistance = 0.01;

        // Add lighting - reduced intensity to not wash out vertex colors
        const ambientLight = new THREE.AmbientLight(0xffffff, 1.0);
        this.scene.add(ambientLight);

        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.3);
        directionalLight.position.set(5, 5, 5);
        this.scene.add(directionalLight);

        // Add reference grid (horizontal)
        const gridHelper = new THREE.GridHelper(50, 50, 0xcccccc, 0xeeeeee);
        gridHelper.position.y = -2;
        this.scene.add(gridHelper);

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

        // Spinning mode checkbox
        const enableSpinningCheckbox = document.getElementById('enable-spinning');
        const rotationGroup = document.getElementById('rotation-group');
        const angularVelocityGroup = document.getElementById('angular-velocity-group');

        enableSpinningCheckbox.addEventListener('change', (e) => {
            this.spinningEnabled = e.target.checked;
            if (this.spinningEnabled) {
                rotationGroup.style.display = 'none';
                angularVelocityGroup.style.display = 'block';
            } else {
                rotationGroup.style.display = 'block';
                angularVelocityGroup.style.display = 'none';
            }
        });

        // Static rotation control
        const rotationSlider = document.getElementById('rotation');
        const rotationValue = document.getElementById('rotation-value');

        rotationSlider.addEventListener('input', (e) => {
            const degrees = parseFloat(e.target.value);
            this.staticRotation = degrees * Math.PI / 180; // Convert to radians
            rotationValue.textContent = `${degrees.toFixed(0)}°`;
        });

        // Angular velocity control
        const angularVelocitySlider = document.getElementById('angular-velocity');
        const angularVelocityValue = document.getElementById('angular-velocity-value');

        angularVelocitySlider.addEventListener('input', (e) => {
            this.angularVelocity = parseFloat(e.target.value);
            angularVelocityValue.textContent = `${this.angularVelocity.toFixed(1)} rad/s`;
        });

        // Time slider control
        const timeSlider = document.getElementById('time-slider');
        const timeSliderValue = document.getElementById('time-slider-value');

        timeSlider.addEventListener('input', (e) => {
            this.currentTime = parseFloat(e.target.value);
            timeSliderValue.textContent = `${this.currentTime.toFixed(1)} s`;
            this.manualTimeControl = true;
        });

        timeSlider.addEventListener('change', (e) => {
            // When user releases the slider, resume animation
            this.manualTimeControl = false;
        });

        // Time scale control
        const timeScaleSlider = document.getElementById('time-scale');
        const timeScaleValue = document.getElementById('time-scale-value');

        timeScaleSlider.addEventListener('input', (e) => {
            this.timeScale = parseFloat(e.target.value);
            timeScaleValue.textContent = `${this.timeScale.toFixed(1)}x`;
        });

        // Play/Pause button
        const playPauseButton = document.getElementById('play-pause');
        playPauseButton.addEventListener('click', () => {
            this.animationRunning = !this.animationRunning;
            playPauseButton.textContent = this.animationRunning ? 'Pause Animation' : 'Resume Animation';
        });

        // Reset button
        document.getElementById('reset').addEventListener('click', () => {
            this.currentTime = -10;
            this.animationRunning = true;
            playPauseButton.textContent = 'Pause Animation';
            timeSlider.value = -10;
            timeSliderValue.textContent = '-10.0 s';
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

        // Checkbox controls
        document.getElementById('show-actual').addEventListener('change', (e) => {
            this.showActualPosition = e.target.checked;
            if (this.actualPositionMesh) {
                this.actualPositionMesh.visible = this.showActualPosition;
            }
        });

        // Initialize displays
        velocitySlider.dispatchEvent(new Event('input'));
        distanceSlider.dispatchEvent(new Event('input'));
        rotationSlider.dispatchEvent(new Event('input'));
        angularVelocitySlider.dispatchEvent(new Event('input'));
        timeSlider.dispatchEvent(new Event('input'));
        timeScaleSlider.dispatchEvent(new Event('input'));
    }

    getLorentzFactor(beta) {
        return 1.0 / Math.sqrt(1.0 - beta * beta);
    }

    // Rotate a vertex around the Y-axis
    rotateVertexY(vertex, angle) {
        const cos = Math.cos(angle);
        const sin = Math.sin(angle);

        return new THREE.Vector3(
            vertex.x * cos + vertex.z * sin,
            vertex.y,
            -vertex.x * sin + vertex.z * cos
        );
    }

    // Calculate the delayed time for a point with rotation
    // Given observer time t, find t' such that the light travel time matches
    // vertexPosRotated should be the vertex position AFTER rotation has been applied
    calculateDelayedTime(t, vertexPosRotated, observerPos) {
        const beta = this.velocity;
        const v = beta * C_VISUAL;
        const x0 = 0; // Object center passes through x=0 at t=0
        const z0 = this.closestDistance;

        // If spinning is disabled OR angular velocity is zero, use static rotation
        // (Using quadratic solution for accuracy and consistency)
        if (!this.spinningEnabled || Math.abs(this.angularVelocity) < 0.001) {
            // Solve for delayed time with static rotation (quadratic solution)
            // The vertex position passed in is already rotated
            const vx = vertexPosRotated.x;
            const vy = vertexPosRotated.y;
            const vz = vertexPosRotated.z;

            const a = v * v - C_VISUAL * C_VISUAL;
            const b = 2 * v * (x0 + vx) + 2 * C_VISUAL * C_VISUAL * t;
            const c = (x0 + vx) * (x0 + vx) + vy * vy + (z0 + vz) * (z0 + vz) - C_VISUAL * C_VISUAL * t * t;

            const discriminant = b * b - 4 * a * c;
            if (discriminant < 0) {
                // No real solution - return approximate delayed time
                const approxDist = Math.sqrt((x0 + vx) * (x0 + vx) + vy * vy + (z0 + vz) * (z0 + vz));
                return t - approxDist / C_VISUAL;
            }

            const sqrtDisc = Math.sqrt(discriminant);
            const t1 = (-b - sqrtDisc) / (2 * a);
            const t2 = (-b + sqrtDisc) / (2 * a);

            // Choose the solution that represents light from the past (t' < t)
            const tDelayed = Math.min(t1, t2);

            // Clamp to reasonable range
            const maxDelay = 20;
            return Math.max(tDelayed, t - maxDelay);
        }

        // Spinning mode with non-zero angular velocity: iterative solution
        // For spinning objects, we need to find the delayed time iteratively
        // because the rotation angle depends on the delayed time itself
        const omega = this.angularVelocity;

        // Note: vertexPosRotated here is the LOCAL unrotated position
        // We'll rotate it based on the delayed time we're solving for
        const vertexPosLocal = vertexPosRotated;

        // Better initial guess: estimate based on object center position at time t
        const centerX = x0 + v * t;
        const centerDist = Math.sqrt(centerX * centerX + z0 * z0);
        let tDelayed = t - centerDist / C_VISUAL;

        // Clamp to reasonable range to prevent extreme values
        const maxDelay = 20; // Don't look more than 20 seconds into the past
        tDelayed = Math.max(tDelayed, t - maxDelay);

        // Iteratively refine tDelayed with damping for stability
        let prevTDelayed = tDelayed;
        for (let iter = 0; iter < 30; iter++) {
            // Calculate rotated position at this delayed time
            const angle = omega * tDelayed;
            const rotatedPos = this.rotateVertexY(vertexPosLocal, angle);

            // World position at delayed time
            const worldX = x0 + v * tDelayed + rotatedPos.x;
            const worldY = rotatedPos.y;
            const worldZ = z0 + rotatedPos.z;

            // Distance from observer to this position
            const distance = Math.sqrt(worldX * worldX + worldY * worldY + worldZ * worldZ);

            // Update delayed time based on light travel time
            let newTDelayed = t - distance / C_VISUAL;

            // Clamp to reasonable range
            newTDelayed = Math.max(newTDelayed, t - maxDelay);

            // Apply damping to prevent oscillation (blend 80% new, 20% old)
            newTDelayed = 0.8 * newTDelayed + 0.2 * tDelayed;

            // Check for convergence
            if (Math.abs(newTDelayed - tDelayed) < 0.0001) {
                return newTDelayed;
            }

            prevTDelayed = tDelayed;
            tDelayed = newTDelayed;
        }

        return tDelayed;
    }

    // Calculate apparent position of vertex at observer time t
    getApparentPosition(t, vertexPosLocal) {
        // Determine rotation angle
        let angle;
        let rotatedPos;

        if (!this.spinningEnabled || Math.abs(this.angularVelocity) < 0.001) {
            // Static rotation: apply rotation first, then solve for delayed time
            angle = this.spinningEnabled ? 0 : this.staticRotation;
            rotatedPos = this.rotateVertexY(vertexPosLocal, angle);

            // Calculate delayed time using the rotated position
            const tDelayed = this.calculateDelayedTime(t, rotatedPos, new THREE.Vector3(0, 0, 0));

            // Position at delayed time (rotation already applied)
            const v = this.velocity * C_VISUAL;
            const x0 = 0; // Object passes through x=0 at t=0

            const x = x0 + v * tDelayed + rotatedPos.x;
            const y = rotatedPos.y;
            const z = this.closestDistance + rotatedPos.z;

            return new THREE.Vector3(x, y, z);
        } else {
            // Spinning: calculateDelayedTime handles rotation internally
            const tDelayed = this.calculateDelayedTime(t, vertexPosLocal, new THREE.Vector3(0, 0, 0));

            // Calculate rotation at the delayed time
            angle = this.angularVelocity * tDelayed;
            rotatedPos = this.rotateVertexY(vertexPosLocal, angle);

            // Position at delayed time
            const v = this.velocity * C_VISUAL;
            const x0 = 0; // Object passes through x=0 at t=0

            const x = x0 + v * tDelayed + rotatedPos.x;
            const y = rotatedPos.y;
            const z = this.closestDistance + rotatedPos.z;

            return new THREE.Vector3(x, y, z);
        }
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

            // Get apparent position (rotation is handled inside based on delayed time)
            const apparentPos = this.getApparentPosition(t, originalPos);

            positions[i] = apparentPos.x;
            positions[i + 1] = apparentPos.y;
            positions[i + 2] = apparentPos.z;

            // Apply colors
            if (originalColors) {
                colors[i] = originalColors[i];
                colors[i + 1] = originalColors[i + 1];
                colors[i + 2] = originalColors[i + 2];
            } else {
                // Default gray color
                colors[i] = 0.8;
                colors[i + 1] = 0.8;
                colors[i + 2] = 0.8;
            }
        }

        this.relativisticMesh.geometry.attributes.position.needsUpdate = true;
        this.relativisticMesh.geometry.attributes.color.needsUpdate = true;
        this.relativisticMesh.geometry.computeVertexNormals();
    }

    addMaterialColorsToGeometry(geometry, material) {
        // Extract color from material and apply as vertex colors
        const vertexCount = geometry.attributes.position.count;
        const colors = new Float32Array(vertexCount * 3);

        if (Array.isArray(material)) {
            material = material[0];
        }

        // Check if material has a texture map
        if (material && material.map && material.map.image) {
            const texture = material.map;

            // Make sure texture image is fully loaded
            // Check if texture has valid dimensions (works for both HTMLImageElement and other image types)
            if ((!texture.image.width || texture.image.width === 0) &&
                (!texture.image.naturalWidth || texture.image.naturalWidth === 0)) {
                console.warn('Texture not fully loaded, using base color', {
                    width: texture.image.width,
                    height: texture.image.height,
                    naturalWidth: texture.image.naturalWidth,
                    complete: texture.image.complete
                });
                this.applyBaseColor(colors, vertexCount, material);
                geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
                return;
            }

            console.log('Extracting colors from texture map');
            const canvas = document.createElement('canvas');
            const ctx = canvas.getContext('2d', { willReadFrequently: true });

            canvas.width = texture.image.width;
            canvas.height = texture.image.height;

            try {
                ctx.drawImage(texture.image, 0, 0);
                const imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
                const uvs = geometry.attributes.uv;

                if (uvs) {
                    // Sample texture at UV coordinates for each vertex
                    for (let i = 0; i < vertexCount; i++) {
                        let u = uvs.getX(i);
                        let v = uvs.getY(i); // Don't flip for GLTF

                        // Clamp UV coordinates to valid range
                        u = Math.max(0, Math.min(1, u));
                        v = Math.max(0, Math.min(1, v));

                        const x = Math.floor(u * (canvas.width - 1));
                        const y = Math.floor(v * (canvas.height - 1));

                        const pixelIndex = (y * canvas.width + x) * 4;

                        colors[i * 3] = imageData.data[pixelIndex] / 255;
                        colors[i * 3 + 1] = imageData.data[pixelIndex + 1] / 255;
                        colors[i * 3 + 2] = imageData.data[pixelIndex + 2] / 255;
                    }
                    console.log('Applied texture colors to', vertexCount, 'vertices');
                    console.log('Sample colors:', {
                        first: [colors[0], colors[1], colors[2]],
                        mid: [colors[Math.floor(vertexCount/2)*3], colors[Math.floor(vertexCount/2)*3+1], colors[Math.floor(vertexCount/2)*3+2]],
                        last: [colors[vertexCount*3-3], colors[vertexCount*3-2], colors[vertexCount*3-1]]
                    });
                } else {
                    console.warn('No UV coordinates found, using base color');
                    this.applyBaseColor(colors, vertexCount, material);
                }
            } catch (error) {
                console.error('Error sampling texture:', error);
                this.applyBaseColor(colors, vertexCount, material);
            }
        } else {
            // No texture, use material base color
            this.applyBaseColor(colors, vertexCount, material);
        }

        geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    }

    applyBaseColor(colors, vertexCount, material) {
        let baseColor = new THREE.Color(0.8, 0.8, 0.8);

        if (material && material.color) {
            baseColor = material.color;
        }

        for (let i = 0; i < vertexCount; i++) {
            colors[i * 3] = baseColor.r;
            colors[i * 3 + 1] = baseColor.g;
            colors[i * 3 + 2] = baseColor.b;
        }

        console.log('Applied base color:', baseColor);
    }

    loadDefaultModel() {
        // Try to load the model from model/ilany.glb first
        const loader = new GLTFLoader();

        // Add cache-busting parameter to force reload if file changes
        const cacheBuster = '?t=' + new Date().getTime();

        loader.load(
            'model/ilany.glb' + cacheBuster,
            (gltf) => {
                // Successfully loaded ilany.glb
                console.log('Loaded model/ilany.glb');
                let geometry = null;
                let material = null;

                gltf.scene.traverse((child) => {
                    if (child.isMesh && !geometry) {
                        geometry = child.geometry.clone();
                        material = child.material;
                    }
                });

                if (geometry) {
                    // Handle texture loading
                    if (!geometry.attributes.color && material && material.map) {
                        const texture = material.map;

                        // Wait for texture image to load if it's not ready
                        if (texture.image && !texture.image.complete) {
                            console.log('Waiting for texture to load...', {
                                complete: texture.image.complete,
                                width: texture.image.width,
                                height: texture.image.height
                            });

                            const processGeometry = () => {
                                console.log('Processing geometry with texture');
                                this.addMaterialColorsToGeometry(geometry, material);

                                // Center and scale the geometry
                                geometry.center();
                                geometry.computeBoundingSphere();
                                const scale = 2.0 / geometry.boundingSphere.radius;
                                geometry.scale(scale, scale, scale);

                                this.setupRelativisticMesh(geometry);
                                this.defaultModelLoaded = true;
                            };

                            texture.image.onload = processGeometry;

                            // Check multiple conditions after delays
                            setTimeout(() => {
                                console.log('Timeout check:', {
                                    complete: texture.image.complete,
                                    width: texture.image.width,
                                    height: texture.image.height,
                                    defaultModelLoaded: this.defaultModelLoaded
                                });

                                if (!this.defaultModelLoaded) {
                                    // Process if texture has dimensions (image data loaded) OR if complete flag is set
                                    if ((texture.image.width > 0 && texture.image.height > 0) || texture.image.complete) {
                                        console.log('Texture has data, processing now');
                                        processGeometry();
                                    } else {
                                        // Try again after another delay
                                        setTimeout(() => {
                                            if (!this.defaultModelLoaded) {
                                                console.log('Second timeout check, forcing processing');
                                                processGeometry();
                                            }
                                        }, 500);
                                    }
                                }
                            }, 100);
                            return;
                        }
                    }

                    // Texture already loaded or no texture
                    if (!geometry.attributes.color && material) {
                        this.addMaterialColorsToGeometry(geometry, material);
                    }

                    // Center and scale the geometry
                    geometry.center();
                    geometry.computeBoundingSphere();
                    const scale = 2.0 / geometry.boundingSphere.radius;
                    geometry.scale(scale, scale, scale);

                    this.setupRelativisticMesh(geometry);
                    this.defaultModelLoaded = true;
                } else {
                    console.warn('No mesh found in ilany.glb, using procedural face');
                    this.loadProceduralFace();
                }
            },
            (progress) => {
                console.log('Loading ilany.glb:', (progress.loaded / progress.total * 100).toFixed(2) + '%');
            },
            (error) => {
                // Failed to load - use procedural face as fallback
                console.log('model/ilany.glb not found, using procedural face');
                this.loadProceduralFace();
            }
        );
    }

    loadProceduralFace() {
        // Create a detailed procedural face model as fallback
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
                let material = null;

                gltf.scene.traverse((child) => {
                    if (child.isMesh && !geometry) {
                        geometry = child.geometry.clone();
                        material = child.material;
                    }
                });

                if (geometry) {
                    // Wait for texture to load before processing
                    if (!geometry.attributes.color && material && material.map) {
                        const texture = material.map;

                        // Wait for texture image to load
                        if (texture.image && !texture.image.complete) {
                            console.log('Waiting for texture to load...', {
                                complete: texture.image.complete,
                                width: texture.image.width,
                                height: texture.image.height
                            });

                            let processed = false;
                            const processGeometry = () => {
                                if (processed) return;
                                processed = true;
                                console.log('Processing geometry with texture');
                                this.addMaterialColorsToGeometry(geometry, material);

                                // Center and scale the geometry
                                geometry.center();
                                geometry.computeBoundingSphere();
                                const scale = 2.0 / geometry.boundingSphere.radius;
                                geometry.scale(scale, scale, scale);

                                this.setupRelativisticMesh(geometry);
                                document.getElementById('loading').style.display = 'none';
                            };

                            texture.image.onload = processGeometry;

                            // Check multiple conditions after delays
                            setTimeout(() => {
                                console.log('Timeout check:', {
                                    complete: texture.image.complete,
                                    width: texture.image.width,
                                    height: texture.image.height,
                                    processed: processed
                                });

                                if (!processed) {
                                    // Process if texture has dimensions (image data loaded) OR if complete flag is set
                                    if ((texture.image.width > 0 && texture.image.height > 0) || texture.image.complete) {
                                        console.log('Texture has data, processing now');
                                        processGeometry();
                                    } else {
                                        // Try again after another delay
                                        setTimeout(() => {
                                            if (!processed) {
                                                console.log('Second timeout check, forcing processing');
                                                processGeometry();
                                            }
                                        }, 500);
                                    }
                                }
                            }, 100);
                            return;
                        }
                    }

                    // Texture already loaded or no texture
                    if (!geometry.attributes.color && material) {
                        this.addMaterialColorsToGeometry(geometry, material);
                    }

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

        // Only create grey colors if the geometry doesn't already have colors
        if (!newGeometry.attributes.color) {
            const colors = new Float32Array(newGeometry.attributes.position.count * 3);
            for (let i = 0; i < colors.length; i += 3) {
                colors[i] = 0.8;
                colors[i + 1] = 0.8;
                colors[i + 2] = 0.8;
            }
            newGeometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        }

        // Create material that uses vertex colors
        // MeshLambertMaterial handles vertex colors better than Phong
        const material = new THREE.MeshLambertMaterial({
            vertexColors: true,
            side: THREE.DoubleSide
        });

        this.relativisticMesh = new THREE.Mesh(newGeometry, material);
        this.scene.add(this.relativisticMesh);

        // Create actual position mesh (transparent, monochrome)
        if (this.actualPositionMesh) {
            this.scene.remove(this.actualPositionMesh);
            this.actualPositionMesh.geometry.dispose();
            this.actualPositionMesh.material.dispose();
        }

        const actualGeometry = geometry.clone();
        const actualMaterial = new THREE.MeshLambertMaterial({
            color: 0x888888, // Gray color
            transparent: true,
            opacity: 0.3,
            side: THREE.DoubleSide,
            wireframe: false
        });

        this.actualPositionMesh = new THREE.Mesh(actualGeometry, actualMaterial);
        this.actualPositionMesh.visible = this.showActualPosition;
        this.scene.add(this.actualPositionMesh);

        // Update mesh position immediately to avoid ghost at origin
        this.updateRelativisticMesh();
        this.updateActualPositionMesh();
    }

    updateActualPositionMesh() {
        if (!this.actualPositionMesh || !this.originalGeometry) return;

        const t = this.currentTime;
        const positions = this.actualPositionMesh.geometry.attributes.position.array;
        const originalPositions = this.originalGeometry.attributes.position.array;

        const beta = this.velocity;
        const gamma = this.getLorentzFactor(beta);
        const v = beta * C_VISUAL;

        // Object's center position at time t (no delayed time, just actual position)
        const x0 = 0; // Object passes through x=0 at t=0
        const centerX = x0 + v * t;
        const centerY = 0;
        const centerZ = this.closestDistance;

        // Rotation angle at current time
        let angle;
        if (this.spinningEnabled) {
            angle = this.angularVelocity * t;
        } else {
            angle = this.staticRotation;
        }

        // Update each vertex
        for (let i = 0; i < originalPositions.length; i += 3) {
            const originalPos = new THREE.Vector3(
                originalPositions[i],
                originalPositions[i + 1],
                originalPositions[i + 2]
            );

            // Apply rotation first
            const rotatedPos = this.rotateVertexY(originalPos, angle);

            // Apply length contraction in the direction of motion (x-direction)
            // Length contraction: L = L0 / gamma
            const contractedX = rotatedPos.x / gamma;
            const contractedY = rotatedPos.y; // No contraction perpendicular to motion
            const contractedZ = rotatedPos.z;

            // Position relative to object center
            positions[i] = centerX + contractedX;
            positions[i + 1] = centerY + contractedY;
            positions[i + 2] = centerZ + contractedZ;
        }

        this.actualPositionMesh.geometry.attributes.position.needsUpdate = true;
        this.actualPositionMesh.geometry.computeVertexNormals();
    }

    updateInfo() {
        document.getElementById('time-display').textContent = this.currentTime.toFixed(2);

        const v = this.velocity * C_VISUAL;
        const x = 0 + v * this.currentTime; // x=0 at t=0
        document.getElementById('position-display').textContent = x.toFixed(2);
        document.getElementById('rel-velocity').textContent = this.velocity.toFixed(3);

        // Update time slider if not being manually controlled
        if (!this.manualTimeControl) {
            const timeSlider = document.getElementById('time-slider');
            const timeSliderValue = document.getElementById('time-slider-value');
            timeSlider.value = this.currentTime;
            timeSliderValue.textContent = `${this.currentTime.toFixed(1)} s`;
        }
    }

    onWindowResize() {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }

    animate() {
        requestAnimationFrame(() => this.animate());

        if (this.animationRunning && !this.manualTimeControl) {
            this.currentTime += 0.016 * this.timeScale; // Approximately 60 FPS

            // Reset if object has passed far away
            if (this.currentTime > 10) {
                this.currentTime = -10;
            }
        }

        // Always update mesh and info, even when paused or manually controlled
        this.updateRelativisticMesh();
        this.updateActualPositionMesh();
        this.updateInfo();

        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
}

// Start the simulation
const simulator = new RelativisticSimulator();
