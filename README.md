# Approximate Lie Symmetries & Differential Geometry Tools

This repository contains a computational and numeric framework designed to analyze systems of ordinary differential equations by finding approximate Lie symmetries (continuous symmetries) of the corresponding solution set of a differential equation.

---

## 🗂️ Repository Structure & File Descriptions

* **`1dquanticNLS.mpl`**  
  Implements the computational routine for analyzing the 1-dimensional quintic Nonlinear Schrödinger (NLS) equation framework.
  
* **`1dquanticNLSreal.mpl`**  
  Handles the real-valued formulation and decomposition of the 1D quintic NLS system for symmetry reduction and analysis.

* **`HybridInvolutiveFormLHPDE.mpl`**  
  A specialized script focused on bringing linear homogeneous partial differential equations (LHPDEs) or related systems into a hybrid involutive form, facilitating exact and approximate analysis.

* **`NumericDiffGeometryTools.mpl`**  
  Core computational tools for numeric differential geometry calculations, mapping out geometric structures on solution manifolds.

* **`NumericDiffGeometryToolsComplex.mpl`**  
  An extension of the differential geometry toolkit tailored for complex-valued manifolds and systems (crucial for wave equations like the NLS).

---

## 🛠️ Technical Scope & Methodology
* **Lie Symmetry Analysis:** Determining infinitesimal generators of continuous transformation groups that leave the differential equation's solution set invariant under perturbation or exact parameters.
* **Involutive Systems:** Using differential algebra to analyze involutive forms of differential systems.
* **Symbolic-Numeric Computation:** Written in **Maple (`.mpl`)**, bridging rigorous algebraic derivation with computational geometry.

---
*Authored by [Sahaib Rakhra](https://sahaibrakhra.neocities.org) with Dr. Greg Reid and Dr. Siyuan Deng*
