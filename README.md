# ⚡ Load Flow Analysis using Decoupled Newton–Raphson Method 

## 📖 Project Overview
This project develops a **generalized load flow analysis program** for power systems using the **Decoupled Newton–Raphson (DNR) method** in MATLAB.  
The program efficiently computes the steady-state operating point of multi-bus power systems by calculating bus voltages, phase angles, and power injections.  
It also records iteration history and displays the apparent power injection at each bus for detailed analysis.

---

## 🧠 Key Features
- ✅ **Generalized Implementation** — Works for any n-bus power system.
- ⚙️ **Automated Y-Bus Formation** from line data.
- ⚡ **Fast Convergence** using the Decoupled Newton–Raphson algorithm.
- 📈 **Iteration History Recording** for voltage magnitudes and angles.
- 🔍 **Computation of Apparent Power** (|S| and ∠S) at each bus.
- 🖥️ **Formatted Command Window Output** for iteration-wise tracking.

---

## 🧩 System Model
### Bus Types
| Type | Description |
|------|--------------|
| 1 | Slack (Reference) Bus |
| 2 | PV (Generator) Bus |
| 3 | PQ (Load) Bus |

🎯 Learning Outcomes
Understand the Decoupled Newton–Raphson load flow algorithm.
Gain hands-on experience with Y-bus matrix formation and power equations.
Learn MATLAB-based power system modeling and simulation.
Visualize convergence behavior across iterations.
