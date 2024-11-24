# Firefly Algorithm Parallelization Project

## Overview

This project implements the Firefly Algorithm (FFA) in both **sequential** and **parallelized** versions using C++. The algorithm simulates the behavior of fireflies to solve optimization problems based on their attractiveness and light intensity.

The parallelized version leverages **OpenMP** to enhance performance on multi-core processors, allowing the evaluation of larger populations and more generations within reduced computation times.

---

## Project Files

1. **`luciernagasSec.cpp`**  
   Contains the sequential implementation of the Firefly Algorithm.  
   Key features:
   - Simulates firefly behavior based on a defined objective function.
   - Evaluates light intensity and movement iteratively for each generation.
   - Outputs firefly positions, intensities, and execution time.

2. **`luciernagasPar.cpp`**  
   Contains the parallelized implementation using OpenMP.  
   Key enhancements:
   - Parallelized initialization, movement, and intensity evaluation.
   - Optimized firefly sorting and movement using OpenMP directives.
   - Enables processing of larger populations and higher iterations efficiently.

---

## How It Works

### Algorithm
1. **Objective Function:**  
   The fitness of fireflies is determined using predefined objective functions. Three test functions (`f_0`, `f_1`, `f_2`) are included in the parallel version for flexibility.
   
2. **Initialization:**  
   Fireflies are randomly distributed within a defined search range.

3. **Movement:**  
   Fireflies move towards brighter ones based on attractiveness, distance, and randomness. The position updates are governed by:
   - Attraction strength (\(\beta\)) inversely proportional to the square of the distance.
   - A randomness parameter (\(\alpha\)) to explore new areas.

4. **Convergence:**  
   Over multiple generations, fireflies converge to regions with higher light intensity, corresponding to optimal solutions of the objective function.

### Parallelization Techniques
- **`#pragma omp parallel for`:** Distributes tasks such as initialization and intensity evaluation across threads.  
- **`#pragma omp parallel for collapse(2)`:** Enhances nested loops for tasks like grid evaluation.  
- **Thread-Safe Operations:** Ensures race conditions are avoided during sorting and updates.

---

## Compilation & Execution

### Prerequisites
- C++ compiler with OpenMP support (e.g., GCC or Clang).  
- Compatible with Windows, Linux, and macOS systems.

### Compilation Commands
1. Sequential Version:
   ```bash
   g++ luciernagasSec.cpp -o fireflySeq -std=c++11
   ```
2. Parallelized Version:
   ```bash
   g++ luciernagasPar.cpp -o fireflyPar -fopenmp -std=c++11
   ```

### Execution Commands
1. Sequential:
   ```bash
   ./fireflySeq
   ```
2. Parallel:
   ```bash
   ./fireflyPar
   ```

---

## Input Parameters

Modify the following parameters in the source code to experiment with different configurations:

- **Population Size (`n`)**: Number of fireflies.  
  - Sequential: `n = 25`  
  - Parallel: `n = 6000`  
- **Generations (`MaxGeneration`)**: Number of iterations.  
- **Search Range (`range`)**: Specifies the x and y limits for firefly movement.  
  ```cpp
  double range[] = { -5, 5, -5, 5 };
  ```
- **Randomness (`alpha`)** and **Attraction Decay (`gamma`)**: Control firefly movement dynamics.

---

## Output

1. **Firefly Positions:** Final positions \((x, y)\) of all fireflies.  
2. **Intensity:** Corresponding intensity values based on the objective function.  
3. **Execution Time:** Total time taken for the computation.

Example Output:
```plaintext
Firefly 0 (1.23, 2.34),  intensity 15.62
Firefly 1 (-3.12, 4.56), intensity 13.14
...
Execution Time = 0.412 seconds
```

---

## Performance Comparison

| Metric                | Sequential (`luciernagasSec.cpp`) | Parallel (`luciernagasPar.cpp`) |
|-----------------------|-----------------------------------|---------------------------------|
| Population Size       | Small (<1000)                   | Large (>1000)                  |
| Execution Time        | Longer                          | Significantly shorter          |
| Scalability           | Limited                         | High (with multi-core CPUs)    |

---

## Customization

- **Add Objective Functions:**  
  Define new fitness functions by modifying the `f_*` methods. Example:
  ```cpp
  double f_custom(double x, double y) {
      return sin(x) + cos(y);
  }
  ```

- **Change Optimization Behavior:**  
  Tweak parameters such as `alpha` and `gamma` to influence exploration and exploitation.

---
