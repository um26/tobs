# tobs

A repository for **Topology Optimization of Binary Structures (TOBS)**.

---

## Table of Contents

1. [Introduction to Topology Optimization](#introduction-to-topology-optimization)
   - [What Is Topology Optimization?](#what-is-topology-optimization)
   - [Key Concepts](#key-concepts)
   - [Common Methods](#common-methods)
   - [Applications](#applications)
2. [TOBS — Topology Optimization of Binary Structures](#tobs--topology-optimization-of-binary-structures)
   - [Overview](#overview)
   - [How TOBS Works](#how-tobs-works)
   - [Advantages](#advantages)
   - [Limitations](#limitations)
3. [References](#references)

---

## Introduction to Topology Optimization

### What Is Topology Optimization?

**Topology optimization** is a mathematical method that optimizes material layout within a given design space, for a given set of loads, boundary conditions, and constraints, with the goal of maximizing the performance of the structure. Unlike size or shape optimization—which work on predefined geometries—topology optimization can freely change the fundamental form (topology) of a design, potentially discovering entirely new structural concepts.

The technique has its roots in the pioneering work of Bendsøe and Kikuchi (1988), and has since become a standard tool in engineering design, particularly for lightweight and high-performance structures.

### Key Concepts

| Concept | Description |
|---|---|
| **Design Domain** | The volume of material from which the optimal structure is carved out |
| **Design Variables** | Per-element quantities (typically density or presence/absence of material) that the optimizer adjusts |
| **Objective Function** | The quantity to be minimized or maximized (e.g., structural compliance, weight, natural frequency) |
| **Constraints** | Restrictions on the design, such as a maximum allowable volume fraction of material |
| **Sensitivity Analysis** | Computation of gradients of the objective and constraints with respect to design variables |
| **Filter / Regularization** | Techniques to avoid numerical issues like checkerboard patterns and mesh-dependency |

### Common Methods

- **SIMP (Solid Isotropic Material with Penalization)**: The most widely used method. Design variables are continuous densities between 0 and 1, penalized to drive solutions toward 0–1 (void–solid) distributions.
- **BESO (Bi-directional Evolutionary Structural Optimization)**: An evolutionary approach that iteratively adds and removes elements based on sensitivity numbers.
- **Level-Set Methods**: The structural boundary is represented implicitly as the zero-level of a scalar field, which is evolved to minimize the objective.
- **TOBS (Topology Optimization of Binary Structures)**: Uses strictly binary (0/1) design variables and solves the problem via sequential integer linear programming (see below).
- **Phase-Field Methods**: Model the solid–void interface as a diffuse region using a phase-field variable, evolved through a Cahn–Hilliard-type equation.

### Applications

Topology optimization is used across a wide range of industries and disciplines:

- **Aerospace & Automotive**: Lightweight structural components (brackets, ribs, bulkheads) that achieve the best stiffness-to-weight ratio.
- **Civil Engineering**: Bridge girders, building frames, and foundation designs.
- **Biomedical Engineering**: Patient-specific implants and bone scaffolds that mimic natural trabecular architecture.
- **Additive Manufacturing (3D Printing)**: Topology-optimized geometries—often organic and lattice-like—are increasingly manufactured directly without traditional tooling constraints.
- **Acoustic and Thermal Design**: Optimizing material layouts for noise reduction, heat dissipation, and thermal management.
- **Microelectromechanical Systems (MEMS)**: Designing compliant mechanisms and actuators at the micro-scale.

---

## TOBS — Topology Optimization of Binary Structures

### Overview

**TOBS** (Topology Optimization of Binary Structures) is a topology optimization method that treats the design variables as **strictly binary integers** (0 = void, 1 = solid). This stands in contrast to density-based methods like SIMP, where design variables are relaxed to a continuous range [0, 1] and intermediate densities are penalized but not eliminated.

The TOBS method was introduced by **Sivapuram and Picelli (2018)** as a way to produce clean, black-and-white designs without the need for post-processing thresholding or filtering.

### How TOBS Works

TOBS solves the binary topology optimization problem through a **sequential integer linear programming (SILP)** approach:

1. **Initialization**: The design domain is discretized into finite elements. Each element is assigned an initial binary state (solid or void).

2. **Finite Element Analysis (FEA)**: The structural response (displacements, stresses, etc.) is computed for the current design.

3. **Sensitivity Analysis**: Gradients of the objective function and constraints are evaluated with respect to each element's binary design variable.

4. **Integer Linear Program (ILP) Subproblem**: The nonlinear optimization problem is linearized around the current design point. The resulting ILP is solved to find the best update (which elements to flip from 0→1 or 1→0).

   The linearized subproblem for a compliance minimization problem can be written as:

   ```
   minimize    ∑ᵢ (∂c/∂xᵢ) · Δxᵢ
   subject to  ∑ᵢ (∂V/∂xᵢ) · Δxᵢ ≤ ΔV_max
               xᵢ + Δxᵢ ∈ {0, 1}   ∀ i
   ```

   where `c` is compliance, `V` is the volume, and `Δxᵢ` is the proposed change to element `i`.

5. **Update and Repeat**: The design is updated and steps 2–4 are repeated until convergence.

A **move limit** controls how many elements can be flipped per iteration, providing stability similar to trust-region methods in continuous optimization.

### Advantages

- **Crisp boundaries**: Produces genuinely binary designs without gray-scale intermediate regions, eliminating the need for post-processing.
- **No penalization required**: Avoids the artificial SIMP penalization exponent, which can distort sensitivities for intermediate densities.
- **Direct physical interpretation**: Every element is either fully solid or fully void at every iteration, so intermediate results are physically meaningful.
- **Constraint handling**: Multiple constraints (volume, stress, frequency, etc.) can be incorporated naturally into the ILP subproblem.
- **Compatible with standard FEA**: Can be implemented on top of any existing finite element solver.

### Limitations

- **Combinatorial complexity**: The underlying binary problem is NP-hard in general; practical scalability depends on the move-limit strategy and the size of the ILP subproblem.
- **Local optima**: Like most topology optimization methods, TOBS can converge to local rather than global optima, and the final result may depend on initialization.
- **ILP solver dependency**: Requires an efficient integer linear programming solver (e.g., CPLEX, Gurobi, or open-source alternatives) to handle the subproblem at each iteration.
- **Slower convergence**: Compared to gradient-based continuous methods, TOBS may require more iterations to converge on fine meshes.

---

## References

1. Bendsøe, M. P., & Kikuchi, N. (1988). *Generating optimal topologies in structural design using a homogenization method*. Computer Methods in Applied Mechanics and Engineering, 71(2), 197–224.

2. Bendsøe, M. P., & Sigmund, O. (2003). *Topology Optimization: Theory, Methods and Applications*. Springer, Berlin.

3. Sivapuram, R., & Picelli, R. (2018). *Topology optimization of binary structures using integer linear programming*. Finite Elements in Analysis and Design, 139, 49–61.

4. Querin, O. M., Steven, G. P., & Xie, Y. M. (1998). *Evolutionary structural optimisation (ESO) using a bidirectional algorithm*. Engineering Computations, 15(8), 1031–1048.

5. Sigmund, O., & Maute, K. (2013). *Topology optimization approaches: A comparative review*. Structural and Multidisciplinary Optimization, 48(6), 1031–1055.

6. van Dijk, N. P., Maute, K., Langelaar, M., & van Keulen, F. (2013). *Level-set methods for structural topology optimization: a review*. Structural and Multidisciplinary Optimization, 48(3), 437–472.