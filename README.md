# Optimization Framework for Kinematic Design of Redundant Manipulators

MATLAB implementation of the kinematic design optimization framework developed
for my B.Tech thesis (NIT Tiruchirappalli). Covers forward/inverse kinematics,
Jacobian-based performance indices, dynamics, and optimization routines for
kinematic design of redundant manipulator architectures, along with
sensitivity and torque verification studies.

---

## 1. Clone the Repository

```bash
git clone git@github.com:shreehank22/BTech-Thesis-Manipulator-Design.git
cd BTech-Thesis-Manipulator-Design
```

---

## 2. Requirements

- MATLAB (tested on R2023a or later)
- Optimization Toolbox

No external submodules — everything runs from the MATLAB scripts/functions in
this repo.

---

## 3. Setup

Run the setup script first to add all subfolders to the MATLAB path and load
the robot model:

```matlab
setup
```

This loads `Arm.mat` (link/DH parameters) and adds all subfolders below to the
MATLAB path.

---

## 4. Workflow

### Step 1 — Sample the Joint Space

```matlab
generate_joint_samples   % batch sampling for evaluation/optimization
q = sample_joint();       % single random joint configuration
```

### Step 2 — Forward Kinematics & Jacobian
*(`Forward Kinematics and Jacobian/`)*

- `Transformation_Matrix.m` — homogeneous transforms from DH parameters
- `Forward_Kinematics.m` / `Forward_Kinematics_RPY.m` — end-effector pose (RPY variant included)
- `Geometric_Jacobian.m` / `Geometric_Jacobian_RPY.m` — manipulator Jacobian

### Step 3 — Inverse Kinematics
*(`Inverse Kinematics/`)*

- `DLS.m` — Damped Least Squares IK solver

### Step 4 — Performance Indices
*(`Performance Indices/`)*

- `compute_MI.m` — Manipulability Index
- `compute_GMI.m` — Global Manipulability Index
- `compute_LCI.m` — Local Conditioning Index
- `compute_GCI.m` — Global Conditioning Index
- `compute_SLI.m` — Service/Link-length-related Index
- `compute_workspace_volume.m` — reachable workspace volume
- `compute_metrics.m` — aggregates the above for a given design

### Step 5 — Optimization
*(`Optimization Codes/`)*

- `Objective.m` — objective function (built on the performance indices above)
- `constraints.m` — design/kinematic constraints
- `Optimizer.m` — optimization driver
- `helper.m` — shared utilities

### Step 6 — Robot Model & Dynamics
*(`Robot Model/`, `Static and Dynamic Torques/`)*

- `robotmodel2.m` — manipulator model construction
- `COM_RST.m`, `dynamic_Torque_RST.m`, `gravity_Torque_RST.m`, `payloadTorque_RST.m` — dynamics for the optimized design (RST formulation)
- `run_trajectory.m` — trajectory execution/simulation
- `Static and Dynamic Torques/` — supporting torque, mass/inertia, and COM utilities (`dynamic_Torque.m`, `gravity_Torque.m`, `inertia_tensor.m`, `linkCOMPositions.m`, `mass.m`, `payloadTorque_vec.m`, `logSO3.m`)

### Step 7 — Wrist Analysis
*(`Wrist Configuration Analysis/`)*

- `wrist_analysis.m` — wrist configuration kinematics
- `Wrist_Jacobian_Verification.m` — wrist Jacobian checks

### Step 8 — Verification
*(`Tests and Verifications/`, `tests.mlx`)*

- **Design Verification** — `link_length_check.m`
- **FK_IK Verification** — `FK_IK_tests.m`, `DLS_Verification.m`, `IK_trajectory_verification.m`, `trajectory_verificationIK.m`
- **Sensitivity Verifications** — `sensitivity_analysis.m`, `link_sensitivity.m`
- **Torque Verifications** — `staitic_dynamic_verification.m`, `dynamic_Torque_verification.m`, `RNEA_Optimization.m`
- `tests.mlx` — Live Script consolidating key verification tests with inline plots/outputs

---

## 5. Project Structure

```
BTech-Thesis-Manipulator-Design/
├── setup.m                              — adds subfolders to path, loads Arm.mat
├── generate_joint_samples.m             — batch joint-space sampling
├── sample_joint.m                       — single random joint configuration sampler
├── tests.mlx                            — Live Script with consolidated verification tests
├── Arm.mat                              — manipulator model (link parameters / DH table)
├── Forward Kinematics and Jacobian/     — FK, Jacobian, transformation matrices
├── Inverse Kinematics/                  — DLS-based IK solver
├── Performance Indices/                 — manipulability/conditioning/workspace metrics
├── Optimization Codes/                  — objective, constraints, optimizer
├── Robot Model/                         — model construction, dynamics (RST), trajectory execution
├── Static and Dynamic Torques/          — torque, inertia, mass, COM utilities
├── Wrist Configuration Analysis/        — wrist kinematics and Jacobian checks
└── Tests and Verifications/
    ├── Design Verification/             — link length checks
    ├── FK_IK Verification/              — FK/IK/trajectory verification
    ├── Sensitivity Verifications/       — link sensitivity analysis
    └── Torque Verifications/            — static/dynamic torque and RNEA verification
```

---

## 6. Notes

This repo contains the code accompanying my B.Tech thesis and is kept
independent of any ongoing/organizational research repositories.
