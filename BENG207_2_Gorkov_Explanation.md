<pre>
Title: Gor'kov Acoustic Radiation Force — Explanation for BENG207 Students
Author: AJM
Context: Supplement to BENG207_2 Rigid Cavity Resonator Model
</pre>

---

## Why Beads and Cells Are Different

This is the single most important practical point for BENG207 students to understand:

| Property         | Polystyrene bead | iPSC organoid | Why it matters                              |
| ---------------- | ---------------- | ------------- | ------------------------------------------- |
| $\rho_p$ (kg/m³) | 1050             | 1050          | Same — density contrast is identical        |
| $c_p$ (m/s)      | 2350             | 1550          | **Very different** — beads are much stiffer |
| $f_1$ (monopole) | 0.62             | 0.128         | **5× smaller for cells**                    |
| $f_2$ (dipole)   | 0.032            | 0.032         | Same                                        |
| **Φ**            | **0.22**         | **0.059**     | **Cells have ~3× lower contrast**           |

The entire difference comes from sound speed — polystyrene is a hard plastic ($c_p$ = 2350 m/s), while cells are basically structured water with some protein ($c_p$ ≈ 1550 m/s, only ~5% above water's 1483 m/s).

**Practical consequence:** If a paper demonstrates acoustic trapping with polystyrene beads at some acoustic pressure $p_0$, you cannot assume the same $p_0$ will work for cells. Since force scales as $\Phi \cdot p_0^2$, to get the same trapping force with cells (Φ = 0.059) as with beads (Φ = 0.22), you need:

$$
p_{0,\text{cells}} = p_{0,\text{beads}} \times \sqrt{\frac{0.22}{0.059}} \approx 1.9 \times p_{0,\text{beads}}
$$

Or in terms of acoustic power (which scales as $p_0^2$):

$$
P_{\text{cells}} \approx 3.7 \times P_{\text{beads}}
$$

This is not a minor correction — you need almost **4× more power** to trap cells as effectively as beads.

## Order-of-Magnitude Force Estimate

Let us calculate the actual force on a 200 µm diameter organoid in our device.

**Given:**

- Φ = 0.059
- a = 100 µm = 10⁻⁴ m
- ρ_f = 1000 kg/m³
- c_f = 1483 m/s
- Assume p₀ = 100 kPa (a reasonable acoustic pressure for a SAW device)

**Gor'kov prefactor:**

$$
\frac{\Phi \cdot \pi \cdot a^3}{3 \cdot \rho_f \cdot c_f^2} = \frac{0.059 \times \pi \times (10^{-4})^3}{3 \times 1000 \times 1483^2} \approx 2.8 \times 10^{-23} \;\text{m}^3/\text{Pa}
$$

**Pressure gradient** (for a perfect standing wave, the maximum gradient of ⟨p²⟩ is $k \cdot p_0^2$ where $k = 2\pi f/c_f$):

$$
k = \frac{2\pi \times 741500}{1483} \approx 3142 \;\text{m}^{-1}
$$

$$
\left|\frac{\partial \langle p^2 \rangle}{\partial x}\right|_{\max} = k \cdot p_0^2 = 3142 \times (10^5)^2 = 3.1 \times 10^{13} \;\text{Pa}^2/\text{m}
$$

**Maximum force:**

$$
F_{\max} = 2.8 \times 10^{-23} \times 3.1 \times 10^{13} \approx 0.9 \;\text{pN}
$$

So the force on a 200 µm organoid at p₀ = 100 kPa is on the order of **~1 piconewton**. For context:

| Force                                   | Magnitude  |
| --------------------------------------- | ---------- |
| Acoustic radiation force (our device)   | ~1 pN      |
| Optical trap on a cell                  | ~1–100 pN  |
| Single molecular motor (kinesin)        | ~5 pN      |
| Cell adhesion (single integrin bond)    | ~50–100 pN |
| Gravity minus buoyancy on a 200 µm cell | ~0.02 pN   |

The acoustic force is well above gravity/buoyancy (so cells will be positioned even against sedimentation) and in the same range as optical trapping. However, it is below typical cell adhesion forces, which means cells can be trapped and positioned before they adhere to a substrate, but once they adhere (after minutes to hours in culture), the acoustic force alone will not detach them. This is actually desirable for our application — we want the acoustic field to position organoids initially, and then have them adhere and grow in place after the field is turned off.

## Summary for the Model Code

In the BENG207_2 C code (`rigid_cavity_model.c`), the Gor'kov addition does the following:

1. **New parameters** defined at the top: `rho_p`, `c_p`, `a_organoid`
2. **Derived in main()**: `f1_monopole`, `f2_dipole`, `Phi_contrast` — computed from the equations above
3. **New function `gorkov_force()`**: takes the existing `radiation_force()` output (which is −∂⟨p²⟩/∂x) and multiplies by the Gor'kov prefactor Φ·π·a³/(3·ρ_f·c_f²)
4. **Fig 4 CSV** now includes two extra columns: `F_gorkov_glass_N` and `F_gorkov_silicon_N` — force in Newtons at p₀ = 1 Pa normalization. To get actual force, multiply these values by your actual p₀².

All additions are marked with `/* BEGIN ADDITION */` and `/* END ADDITION */` comments in the code.

## References

- Gor'kov, L. P. (1962). On the forces acting on a small particle in an acoustical field in an ideal fluid. *Soviet Physics Doklady*, 6, 773–775.
- Bruus, H. (2012). Acoustofluidics 7: The acoustic radiation force on small particles. *Lab on a Chip*, 12, 1014–1021. ← **Best reference; recommended reading for BENG207 students.**
- Settnes, M., & Bruus, H. (2012). Forces acting on a small particle in an acoustical field in a viscous fluid. *Physical Review E*, 85, 016327. ← Extension to viscous fluids.
- Cohen, S., et al. (2020). Large-scale acoustic-driven neuronal patterning and directed outgrowth. *Scientific Reports*, 10, 4932.
- Barnkob, R., Augustsson, P., Laurell, T., & Bruus, H. (2010). Measuring the local pressure amplitude in microchannel acoustophoresis. *Lab on a Chip*, 10, 563–570. ← Experimental measurement of acoustic contrast for cells.
