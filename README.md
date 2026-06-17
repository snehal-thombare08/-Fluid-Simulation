# 💧 Fluid Simulation

A real-time **Navier-Stokes fluid simulator** built with **C++ and SFML 3.0**, featuring interactive painting with realistic fluid dynamics.

![Fluid Simulation](https://raw.githubusercontent.com/snehal-thombare08/-Fluid-Simulation/main/Screenshot%202026-06-17%20164302.png)

## ✨ Features

- **Navier-Stokes physics** — diffusion, advection, pressure projection
- **4 Visual Modes** — Cyan, Lava, White, Rainbow
- **Interactive painting** — drag mouse to add fluid with velocity
- **Velocity field overlay** — see the fluid flow vectors
- **Erase mode** — RMB to clear fluid areas
- Real-time pixel rendering at 1200×800

## 🎮 Controls

| Key / Button | Action |
|---|---|
| LMB (drag) | Paint fluid |
| RMB (drag) | Erase fluid |
| `1` | Cyan mode |
| `2` | Lava mode 🌋 |
| `3` | White mode |
| `4` | Rainbow mode 🌈 |
| `V` | Toggle velocity field |
| `Space` | Clear all fluid |
| `Esc` | Quit |

## 🚀 How to Run

1. Download **FluidSimulation-v1.0-Windows.zip**
2. Extract all files to a folder
3. Run **FluidSimulation.exe**

> No installation required. All DLLs included.

## 🛠️ Build from Source

**Requirements:**
- C++17 compiler (MinGW-w64)
- SFML 3.0 (via vcpkg)
- CMake 3.16+

```bash
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -G "MinGW Makefiles"
mingw32-make
```

## 🔬 Physics

Implements **Jos Stam's stable fluid solver**:
- **Diffusion** — fluid spreads over time
- **Advection** — fluid moves with velocity field
- **Pressure projection** — ensures incompressibility (∇·v = 0)
- **Gauss-Seidel iterations** for linear solve

## 📁 Project Structure

```
Fluid Simulation/
├── src/
│   └── main.cpp
├── CMakeLists.txt
└── .gitignore
```

## 🧰 Tech Stack

- **Language:** C++17
- **Graphics:** SFML 3.0
- **Build:** CMake + MinGW
- **Physics:** Navier-Stokes (Jos Stam method)

---
*Part of a C++ graphics & simulation portfolio*
