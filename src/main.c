/* rigid_cavity_model.c
 *
 * BENG207_2: Rigid Cavity Resonator Model
 *
 * Lateral standing wave between rigid sidewalls for iPSC organoid formation.
 * Organoid spacing ~1 mm matches CI electrode spacing (~30 mm cochlea / 22 ch).
 * Operating frequency ~741 kHz (lambda/2 = 1 mm in water).
 *
 * Stack (bottom to top, vertical energy coupling):
 *   1. LiNbO3 transducer (source)
 *   2. Ultrasound gel couplant (Kelvin-Voigt)
 *   3. Borosilicate glass superstrate
 *   4. Fluid channel (water) — lateral standing wave forms HERE
 *
 * The standing wave propagates LATERALLY (along x, parallel to glass)
 * between two rigid sidewalls, NOT vertically through the stack.
 *
 * See README file and BENG207_model_2_schematic.svg for geometry.
 *
 * Part A: Lateral cavity — pressure field, contrast, radiation force
 * Part B: Vertical coupling — impedance transformer (from BENG207_1)
 * Part C: Dual-transducer — pressure, phase sweep, amplitude balance
 * Part D: PDMS wall transmission — lossy wall attenuation analysis
 *
 * Compile (Linux/Mac):
 *   gcc -std=c99 -O2 -o cavity_model rigid_cavity_model.c -lm
 *
 * Compile (Windows MSVC):
 *   cl /O2 rigid_cavity_model.c /link /out:cavity_model.exe
 *
 * Compile (Windows MinGW):
 *   gcc -O2 -o cavity_model.exe rigid_cavity_model.c -lm
 *
 * Run:
 *   ./cavity_model        (Linux/Mac)
 *   cavity_model.exe      (Windows)
 *   CSVs appear in PLOT_CSV/ subfolder next to the executable.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================================================
 * Platform-specific directory creation
 * ============================================================ */
#ifdef _WIN32
  #include <direct.h>
  #define MKDIR(path) _mkdir(path)
#else
  #include <sys/stat.h>
  #include <errno.h>
  #define MKDIR(path) mkdir(path, 0755)
#endif

/* ============================================================
 * Complex number helpers (MSVC lacks <complex.h>)
 * ============================================================ */

typedef struct { double re; double im; } Cpx;

static Cpx cpx(double re, double im)
{
    Cpx z; z.re = re; z.im = im; return z;
}

static Cpx cpx_add(Cpx a, Cpx b)
{
    return cpx(a.re + b.re, a.im + b.im);
}

static Cpx cpx_sub(Cpx a, Cpx b)
{
    return cpx(a.re - b.re, a.im - b.im);
}

static Cpx cpx_mul(Cpx a, Cpx b)
{
    return cpx(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}

static Cpx cpx_div(Cpx a, Cpx b)
{
    double d = b.re*b.re + b.im*b.im;
    return cpx((a.re*b.re + a.im*b.im) / d,
               (a.im*b.re - a.re*b.im) / d);
}

static Cpx cpx_scale(double s, Cpx a)
{
    return cpx(s * a.re, s * a.im);
}

static double cpx_abs(Cpx a)
{
    return sqrt(a.re*a.re + a.im*a.im);
}

static Cpx cpx_exp(Cpx a)
{
    double e = exp(a.re);
    return cpx(e * cos(a.im), e * sin(a.im));
}

static Cpx cpx_conj(Cpx a)
{
    return cpx(a.re, -a.im);
}

static Cpx cpx_sin(Cpx a)
{
    return cpx(sin(a.re) * cosh(a.im), cos(a.re) * sinh(a.im));
}

static Cpx cpx_cos(Cpx a)
{
    return cpx(cos(a.re) * cosh(a.im), -sin(a.re) * sinh(a.im));
}

static Cpx cpx_sqrt(Cpx a)
{
    double r = sqrt(cpx_abs(a));
    double theta = atan2(a.im, a.re) / 2.0;
    return cpx(r * cos(theta), r * sin(theta));
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Output directory — relative to executable */
#define OUTPUT_DIR "/home/shh/Projects/BENG207_Echo_Chamber_AMN/PLOT_CSV"

/* ============================================================
 * Material parameters
 * ============================================================ */
/* Parameters have been changed to reflect what we currently use */

/* Fluid channel (water) — where the lateral standing wave forms */
static const double rho_f = 1000.0;     /* kg/m^3 */
static const double c_f   = 1500;     /* m/s */

/* Fluid attenuation model: alpha = alpha_coeff * f^2
 * For water at 25C: alpha_coeff ~ 25e-15 Np/m/Hz^2
 * At 741 kHz: alpha ~ 25e-15 * (741e3)^2 ~ 1.4e-5 Np/m (negligible)
 * At 20 MHz:  alpha ~ 25e-15 * (20e6)^2  ~ 0.01 Np/m              */
static const double alpha_coeff = 25.0e-15;  /* Np/m/Hz^2 */

/* Wall materials */
/* Glass*/
static const double rho_glass = 2210.0;     /* kg/m^3 */
static const double c_glass   = 5700.0;     /* m/s */


static const double rho_si    = 2330.0;     /* kg/m^3 */
static const double c_si      = 8433.0;     /* m/s (longitudinal) */

/* Acrylic (PMMA) — wall material for dual-transducer config
 * Hou et al. 2020 (Extreme Mech. Lett. 37, 100716)              */
static const double rho_acrylic = 1190.0;   /* kg/m^3 (PMMA) */
static const double c_acrylic   = 2750.0;   /* m/s (longitudinal, PMMA) */

/* PZT-5A transducer — for dual-transducer configuration          */
static const double rho_pzt     = 7750.0;   /* kg/m^3 (PZT-5A) */
static const double c_pzt       = 4350.0;   /* m/s (PZT-5A longitudinal) */

/* PDMS (Sylgard 184, 10:1 mix ratio) — channel wall material
 * for the immersed dual-transducer configuration.
 * PDMS has very low impedance (~water) but HIGH attenuation.
 *
 * Attenuation: ~3-12 dB/cm/MHz (varies with cure ratio)
 * We use 5 dB/cm/MHz as a central estimate.
 *
 * Refs: Tsou et al., J. Micromech. Microeng. 18 (2008)
 *       Carugo et al., Biomicrofluidics 6 (2012)                */
static const double rho_pdms    = 1040.0;   /* kg/m^3 */
static const double c_pdms      = 1030.0;   /* m/s */
static const double alpha_pdms_per_MHz_dBcm = 5.0;  /* dB/cm/MHz */

/* ============================================================
 * BEGIN ADDITION: Gor'kov acoustic radiation force parameters
 *
 * The Gor'kov potential for a small sphere in a standing wave:
 *   U = (4/3)*pi*a^3 * [ f1/(2*rho_f*c_f^2) * <p^2>
 *                       - 3*f2*rho_f/(4)      * <v^2> ]
 *
 * where f1 = 1 - kappa_p/kappa_f   (monopole coefficient)
 *       f2 = 2*(rho_p - rho_f) / (2*rho_p + rho_f)  (dipole)
 *
 * For a 1-D standing wave p = p0*cos(kx), the radiation force
 * simplifies to:
 *   F_rad = 4*pi*a^3 * k * Phi * <E_ac> * sin(2kx)
 *
 * where Phi = f1/3 + f2/2  is the acoustic contrast factor,
 *       <E_ac> = p0^2 / (4*rho_f*c_f^2)  is the energy density.
 *
 * Phi > 0  =>  particles migrate to pressure NODES  (our case)
 * Phi < 0  =>  particles migrate to pressure ANTINODES
 *
 * References:
 *   Gor'kov, Sov. Phys. Doklady 6, 773 (1962)
 *   Bruus, Lab Chip 12, 1014 (2012)
 *   Cohen et al., Sci. Rep. 10, 4932 (2020) — neuronal patterning
 * ============================================================ */

/* Particle (iPSC organoid / neural aggregate) properties */
static const double rho_p = 1050.0;   /* kg/m^3  — cell density, typical
                                        * range 1040-1080 for mammalian cells.
                                        * Cohen et al. used DRG neurons
                                        * (~1050 kg/m^3). Adjust if measured. */

static const double c_p   = 1550.0;   /* m/s — longitudinal sound speed in
                                        * cell aggregates. Typical range
                                        * 1520-1600 m/s for soft tissue.
                                        * Higher than water due to protein
                                        * content. Adjust if measured. */

static const double a_organoid = 100.0e-6;  /* m — organoid radius (~200 um
                                              * diameter). This is the nominal
                                              * radius for force calculation;
                                              * actual organoids vary. */

/* Derived Gor'kov quantities (set in main) */
static double f1_monopole;  /* = 1 - kappa_p/kappa_f = 1 - (rho_f*c_f^2)/(rho_p*c_p^2) */
static double f2_dipole;    /* = 2*(rho_p - rho_f) / (2*rho_p + rho_f) */
static double Phi_contrast; /* = f1/3 + f2/2  — the acoustic contrast factor */

/* END ADDITION: Gor'kov parameters */

/* Vertical stack (for coupling sub-model, Part B) */
static const double rho_1    = 4647.0;      /* LiNbO3 kg/m^3 */
static const double c_1      = 6570.0;      /* LiNbO3 m/s */

static const double rho_2    = 1020.0;      /* couplant kg/m^3 */
static const double c_2      = 1500.0;      /* couplant m/s */

static const double rho_3    = 2230.0;      /* glass kg/m^3 */
static const double c_3      = 5640.0;      /* glass m/s */

/* Derived impedances (set in main) */
static double Z_f, Z_glass, Z_si, Z_acrylic, Z_pzt, Z_pdms, Z_1, Z_2, Z_3;

/* ============================================================
 * Utilities
 * ============================================================ */

static void linspace(double start, double end, int n, double *out)
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = start + (end - start) * i / (n - 1);
}

static double dmax(double a, double b)
{
    return (a > b) ? a : b;
}

static void build_path(char *buf, int bufsize, const char *filename)
{
#ifdef _WIN32
    _snprintf(buf, bufsize, "%s\\%s", OUTPUT_DIR, filename);
#else
    snprintf(buf, bufsize, "%s/%s", OUTPUT_DIR, filename);
#endif
}

/* ============================================================
 * PART A: Lateral standing wave in rigid cavity
 *
 * p(x) = p0 * exp(-(alpha - ik)x) + r * p0 * exp(-(alpha - ik)(2L - x))
 *
 * where alpha = attenuation, k = 2*pi*f/c_f, r = reflection coeff,
 * L = channel length, x = position along channel.
 * ============================================================ */

/* Pressure at position x in the cavity */
static Cpx pressure_at(double x, double L, double freq,
                        double r_abs, double alpha)
{
    double k = 2.0 * M_PI * freq / c_f;
    Cpx gamma = cpx(-alpha, k);       /* propagation constant */
    Cpx r = cpx(r_abs, 0);            /* real reflection coeff */

    /* Forward wave: exp(-gamma * x) */
    Cpx fwd = cpx_exp(cpx_scale(-1.0, cpx_scale(x, gamma)));

    /* Reflected wave: r * exp(-gamma * (2L - x)) */
    Cpx ref = cpx_mul(r, cpx_exp(cpx_scale(-1.0, cpx_scale(2.0*L - x, gamma))));

    return cpx_add(fwd, ref);
}

/* Pressure squared (time-averaged): |p(x)|^2 */
static double pressure_sq(double x, double L, double freq,
                           double r_abs, double alpha)
{
    Cpx p = pressure_at(x, L, freq, r_abs, alpha);
    return p.re*p.re + p.im*p.im;
}

/* Acoustic radiation force proportional to -d/dx[<p^2>]
 * Computed by central finite difference */
static double radiation_force(double x, double L, double freq,
                               double r_abs, double alpha)
{
    double dx = L * 1e-5;  /* small step */
    double p2_plus, p2_minus;

    if (x - dx < 0) {
        p2_plus  = pressure_sq(x + dx, L, freq, r_abs, alpha);
        p2_minus = pressure_sq(x, L, freq, r_abs, alpha);
        return -(p2_plus - p2_minus) / dx;
    }
    if (x + dx > L) {
        p2_plus  = pressure_sq(x, L, freq, r_abs, alpha);
        p2_minus = pressure_sq(x - dx, L, freq, r_abs, alpha);
        return -(p2_plus - p2_minus) / dx;
    }

    p2_plus  = pressure_sq(x + dx, L, freq, r_abs, alpha);
    p2_minus = pressure_sq(x - dx, L, freq, r_abs, alpha);
    return -(p2_plus - p2_minus) / (2.0 * dx);
}

/* ============================================================
 * BEGIN ADDITION: Gor'kov radiation force (physical units)
 *
 * Converts the dimensionless pressure-gradient force into
 * absolute Newtons using the Gor'kov prefactor:
 *
 *   F_gorkov = Phi * (pi * a^3) / (3 * rho_f * c_f^2)
 *              * (-d/dx[<p^2>])
 *
 * This is derived from the Gor'kov potential gradient.
 * The factor (pi*a^3)/(3*rho_f*c_f^2) has units of
 * m^3 / (kg/m^3 * m^2/s^2) = m^3 * s^2 / kg = N / (Pa^2/m).
 *
 * To get force in Newtons, you still need to scale the
 * pressure field by the actual pressure amplitude p0 (Pa).
 * The output here assumes p0 = 1 Pa (i.e., the p field is
 * normalized). Multiply by p0^2 for actual force.
 *
 * Returns force in Newtons (for p0 = 1 Pa normalization).
 * Positive = toward nearest pressure node (trapping).
 * ============================================================ */
static double gorkov_force(double x, double L, double freq,
                            double r_abs, double alpha)
{
    /* Raw pressure-gradient force: -d/dx[<p^2>] in Pa^2/m */
    double F_raw = radiation_force(x, L, freq, r_abs, alpha);

    /* Gor'kov prefactor: Phi * pi * a^3 / (3 * rho_f * c_f^2) */
    double a3 = a_organoid * a_organoid * a_organoid;
    double prefactor = Phi_contrast * M_PI * a3
                       / (3.0 * rho_f * c_f * c_f);

    return prefactor * F_raw;
}
/* END ADDITION: Gor'kov radiation force */

/* Standing wave ratio at position x */
static double sw_ratio(double x, double L, double r_abs, double alpha)
{
    return r_abs * exp(-2.0 * alpha * (L - x));
}

/* Contrast from standing wave ratio */
static double contrast_from_swr(double S)
{
    if (S >= 0.9999) return 1e10;  /* avoid division by zero */
    return ((1.0 + S) / (1.0 - S)) * ((1.0 + S) / (1.0 - S));
}

/* Reflection coefficient from impedance mismatch */
static double reflection_coeff(double Z_wall, double Z_fluid)
{
    return fabs((Z_wall - Z_fluid) / (Z_wall + Z_fluid));
}

/* ============================================================
 * PART B: Vertical coupling (impedance transformer)
 *
 * Z_in(d) = Z2 * (Z3 + i*Z2*tan(k2*d)) / (Z2 + i*Z3*tan(k2*d))
 * T_power = 1 - |r_vert|^2
 * where r_vert = (Z_in - Z1) / (Z_in + Z1)
 * ============================================================ */

static double coupling_T_power(double freq, double d_couplant)
{
    double omega = 2.0 * M_PI * freq;
    double k2 = omega / c_2;
    double k2d = k2 * d_couplant;
    Cpx tan_k2d, Z_in, r_vert;
    Cpx numerator, denominator;
    double r_abs_sq;

    /* tan(k2d) — real since couplant is lossless in this sub-model */
    tan_k2d = cpx(tan(k2d), 0);

    /* Z_in = Z2 * (Z3 + i*Z2*tan(k2d)) / (Z2 + i*Z3*tan(k2d)) */
    numerator = cpx_add(cpx(Z_3, 0), cpx_mul(cpx(0, Z_2), tan_k2d));
    denominator = cpx_add(cpx(Z_2, 0), cpx_mul(cpx(0, Z_3), tan_k2d));
    Z_in = cpx_mul(cpx(Z_2, 0), cpx_div(numerator, denominator));

    /* r_vert = (Z_in - Z1) / (Z_in + Z1) */
    r_vert = cpx_div(cpx_sub(Z_in, cpx(Z_1, 0)),
                     cpx_add(Z_in, cpx(Z_1, 0)));

    r_abs_sq = r_vert.re*r_vert.re + r_vert.im*r_vert.im;
    return 1.0 - r_abs_sq;
}

/* Return Z_in as complex (for Fig 5 output) */
static Cpx coupling_Z_in(double freq, double d_couplant)
{
    double omega = 2.0 * M_PI * freq;
    double k2 = omega / c_2;
    double k2d = k2 * d_couplant;
    Cpx tan_k2d, numerator, denominator;

    tan_k2d = cpx(tan(k2d), 0);
    numerator = cpx_add(cpx(Z_3, 0), cpx_mul(cpx(0, Z_2), tan_k2d));
    denominator = cpx_add(cpx(Z_2, 0), cpx_mul(cpx(0, Z_3), tan_k2d));
    return cpx_mul(cpx(Z_2, 0), cpx_div(numerator, denominator));
}

/* ============================================================
 * OUTPUT 1: Lateral pressure field |p(x)|
 *
 * Three channel lengths (5, 10, 20 mm) at resonance,
 * for glass and silicon sidewalls.
 * ============================================================ */

static void output_pressure_field(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 20e-3};
    int NL = 3;
    double r_glass, r_silicon;
    int Nx = 2000;
    int iL, ix;

    build_path(fname, sizeof(fname), "fig1_pressure_field.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_glass = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);

    fprintf(fp, "L_mm,x_mm,p_abs_glass,p_abs_silicon,p_sq_glass,p_sq_silicon\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        /* Choose mode number for ~1 mm spacing */
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double freq = n_mode * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;

        for (ix = 0; ix <= Nx; ix++) {
            double x = L * (double)ix / Nx;
            Cpx p_g = pressure_at(x, L, freq, r_glass, alpha);
            Cpx p_s = pressure_at(x, L, freq, r_silicon, alpha);

            fprintf(fp, "%.4f,%.6f,%.8f,%.8f,%.8f,%.8f\n",
                    L*1e3, x*1e3,
                    cpx_abs(p_g), cpx_abs(p_s),
                    p_g.re*p_g.re + p_g.im*p_g.im,
                    p_s.re*p_s.re + p_s.im*p_s.im);
        }
    }

    fclose(fp);
    printf("  Written: %s (3 lengths x %d positions)\n", fname, Nx+1);
}

/* ============================================================
 * OUTPUT 2: Standing wave contrast map
 *
 * 2D heatmap: contrast (dB) at channel midpoint
 * vs L (1-25 mm) and |r| (0.5-0.99)
 * ============================================================ */

static void output_contrast_map(void)
{
    char fname[256];
    FILE *fp;
    int NL = 200, Nr = 200;
    double L_arr[200], r_arr[200];
    int iL, ir;

    build_path(fname, sizeof(fname), "fig2_contrast_map.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1e-3, 25e-3, NL, L_arr);
    linspace(0.50, 0.99, Nr, r_arr);

    fprintf(fp, "L_mm,r_abs,contrast_dB,contrast_linear,n_mode,S_mid\n");

    for (ir = 0; ir < Nr; ir++) {
        for (iL = 0; iL < NL; iL++) {
            double L = L_arr[iL];
            int n_mode = (int)(L / 1.0e-3 + 0.5);
            double freq = n_mode * c_f / (2.0 * L);
            double alpha = alpha_coeff * freq * freq;
            double S_mid = sw_ratio(L/2.0, L, r_arr[ir], alpha);
            double C = contrast_from_swr(S_mid);
            double C_dB = 10.0 * log10(dmax(C, 1e-15));

            fprintf(fp, "%.4f,%.4f,%.4f,%.4f,%d,%.6f\n",
                    L*1e3, r_arr[ir], C_dB, C, n_mode, S_mid);
        }
    }

    fclose(fp);
    printf("  Written: %s (%dx%d)\n", fname, NL, Nr);
}

/* ============================================================
 * OUTPUT 3: Resonant mode structure
 *
 * Mode frequencies f_n vs channel length,
 * number of modes in transducer bandwidth.
 * ============================================================ */

static void output_mode_structure(void)
{
    char fname[256];
    FILE *fp;
    int NL = 500;
    double L_arr[500];
    double bw_center = 741.5e3;  /* target frequency */
    double bw_half = 50.0e3;     /* assumed BW +/- 50 kHz */
    int iL;

    build_path(fname, sizeof(fname), "fig3_mode_structure.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1e-3, 25e-3, NL, L_arr);

    fprintf(fp, "L_mm,n_target,f_target_kHz,delta_f_kHz,"
                "n_modes_in_bw,node_spacing_mm,organoid_count\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_arr[iL];
        double delta_f = c_f / (2.0 * L);  /* mode spacing */

        /* Mode number for ~1 mm node spacing */
        int n_target = (int)(L / 1.0e-3 + 0.5);
        if (n_target < 1) n_target = 1;

        double f_target = n_target * delta_f;
        double node_spacing = L / (double)n_target;

        /* Count modes within bandwidth */
        int n_lo = (int)ceil((bw_center - bw_half) / delta_f);
        int n_hi = (int)floor((bw_center + bw_half) / delta_f);
        int n_modes_in_bw = 0;
        if (n_hi >= n_lo && n_lo >= 1) n_modes_in_bw = n_hi - n_lo + 1;

        fprintf(fp, "%.4f,%d,%.4f,%.4f,%d,%.6f,%d\n",
                L*1e3, n_target, f_target/1e3, delta_f/1e3,
                n_modes_in_bw, node_spacing*1e3, n_target);
    }

    fclose(fp);
    printf("  Written: %s (%d lengths)\n", fname, NL);
}

/* ============================================================
 * OUTPUT 4: Nodal uniformity profile
 *
 * Radiation force magnitude at each node along the channel,
 * normalized to the strongest node.
 * Glass vs silicon, L = 5, 10, 20 mm.
 * ============================================================ */

static void output_nodal_uniformity(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 20e-3};
    int NL = 3;
    double r_glass, r_silicon;
    int iL, in;

    build_path(fname, sizeof(fname), "fig4_nodal_uniformity.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_glass = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);

    fprintf(fp, "L_mm,node_index,x_mm,"
                "p_sq_glass,p_sq_silicon,"
                "F_glass,F_silicon,"
                "F_glass_norm,F_silicon_norm,"
                /* BEGIN ADDITION: Gor'kov force columns (in Newtons, p0=1 Pa) */
                "F_gorkov_glass_N,F_gorkov_silicon_N\n");
                /* END ADDITION */

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double freq = n_mode * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;
        double node_spacing = L / (double)n_mode;
        double F_max_g = 0, F_max_s = 0;
        double *F_g, *F_s, *p2_g, *p2_s, *x_nodes;
        int i;

        F_g = (double *)malloc(n_mode * sizeof(double));
        F_s = (double *)malloc(n_mode * sizeof(double));
        p2_g = (double *)malloc(n_mode * sizeof(double));
        p2_s = (double *)malloc(n_mode * sizeof(double));
        x_nodes = (double *)malloc(n_mode * sizeof(double));

        /* Compute force at each node */
        for (in = 0; in < n_mode; in++) {
            double x = (in + 0.5) * node_spacing;  /* node centers */
            x_nodes[in] = x;
            p2_g[in] = pressure_sq(x, L, freq, r_glass, alpha);
            p2_s[in] = pressure_sq(x, L, freq, r_silicon, alpha);

            /* Force magnitude (absolute value of gradient) */
            F_g[in] = fabs(radiation_force(x, L, freq, r_glass, alpha));
            F_s[in] = fabs(radiation_force(x, L, freq, r_silicon, alpha));

            if (F_g[in] > F_max_g) F_max_g = F_g[in];
            if (F_s[in] > F_max_s) F_max_s = F_s[in];
        }

        /* Write normalized */
        for (i = 0; i < n_mode; i++) {
            /* BEGIN ADDITION: compute Gor'kov force at each node */
            double Fg_gorkov = fabs(gorkov_force(x_nodes[i], L, freq,
                                                  r_glass, alpha));
            double Fs_gorkov = fabs(gorkov_force(x_nodes[i], L, freq,
                                                  r_silicon, alpha));
            /* END ADDITION */

            fprintf(fp, "%.4f,%d,%.6f,%.8f,%.8f,%.8e,%.8e,%.6f,%.6f,"
                        /* BEGIN ADDITION: Gor'kov columns */
                        "%.8e,%.8e\n",
                        /* END ADDITION */
                    L*1e3, i, x_nodes[i]*1e3,
                    p2_g[i], p2_s[i],
                    F_g[i], F_s[i],
                    (F_max_g > 0) ? F_g[i]/F_max_g : 0,
                    (F_max_s > 0) ? F_s[i]/F_max_s : 0,
                    /* BEGIN ADDITION: Gor'kov values */
                    Fg_gorkov, Fs_gorkov);
                    /* END ADDITION */
        }

        free(F_g); free(F_s); free(p2_g); free(p2_s); free(x_nodes);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 5: Coupling layer characterization
 *
 * Z_in (real, imag, mag) and T_power vs couplant thickness
 * at 741 kHz.
 * ============================================================ */

static void output_coupling_layer(void)
{
    char fname[256];
    FILE *fp;
    int Nd = 2000;
    double d_arr[2000];
    double freq = 741.5e3;
    double lambda2 = c_2 / freq;  /* ~2.02 mm at 741 kHz */
    int i;

    build_path(fname, sizeof(fname), "fig5_coupling_layer.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(0.1e-6, 200e-6, Nd, d_arr);

    fprintf(fp, "d_um,d_over_lambda,Zin_re_MRayl,Zin_im_MRayl,"
                "Zin_mag_MRayl,T_power,T_power_dB\n");

    for (i = 0; i < Nd; i++) {
        Cpx Zin = coupling_Z_in(freq, d_arr[i]);
        double T = coupling_T_power(freq, d_arr[i]);
        double T_dB = 10.0 * log10(dmax(T, 1e-15));

        fprintf(fp, "%.4f,%.6f,%.6f,%.6f,%.6f,%.8f,%.4f\n",
                d_arr[i]*1e6, d_arr[i]/lambda2,
                Zin.re/1e6, Zin.im/1e6, cpx_abs(Zin)/1e6,
                T, T_dB);
    }

    fclose(fp);
    printf("  Written: %s (%d thicknesses, lambda2=%.3f mm at %.0f kHz)\n",
           fname, Nd, lambda2*1e3, freq/1e3);
}

/* ============================================================
 * OUTPUT 6: Frequency sensitivity
 *
 * Organoid spacing error and contrast loss vs frequency
 * detuning for different channel lengths.
 * ============================================================ */

static void output_frequency_sensitivity(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 15e-3, 20e-3};
    int NL = 4;
    int Ndf = 500;
    double df_arr[500];  /* fractional detuning */
    double r_wall;
    int iL, idf;

    build_path(fname, sizeof(fname), "fig6_freq_sensitivity.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_wall = reflection_coeff(Z_glass, Z_f);
    linspace(-0.02, 0.02, Ndf, df_arr);  /* +/- 2% detuning */

    fprintf(fp, "L_mm,df_frac,df_kHz,f_actual_kHz,"
                "spacing_error_um,contrast_mid_dB,p_sq_mid\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double f_res = n_mode * c_f / (2.0 * L);

        for (idf = 0; idf < Ndf; idf++) {
            double f_actual = f_res * (1.0 + df_arr[idf]);
            double alpha = alpha_coeff * f_actual * f_actual;
            double lambda_actual = c_f / f_actual;
            double actual_spacing = lambda_actual / 2.0;
            double spacing_error = (actual_spacing - 1.0e-3) * 1e6;  /* um */

            /* Contrast at midpoint */
            double S_mid = sw_ratio(L/2.0, L, r_wall, alpha);
            double C_mid = contrast_from_swr(S_mid);
            double C_dB = 10.0 * log10(dmax(C_mid, 1e-15));

            /* Pressure squared at midpoint */
            double p2_mid = pressure_sq(L/2.0, L, f_actual, r_wall, alpha);

            fprintf(fp, "%.4f,%.6f,%.4f,%.4f,%.4f,%.4f,%.8f\n",
                    L*1e3, df_arr[idf], (f_actual - f_res)/1e3,
                    f_actual/1e3, spacing_error, C_dB, p2_mid);
        }
    }

    fclose(fp);
    printf("  Written: %s (%d lengths x %d detunings)\n", fname, NL, Ndf);
}

/* ============================================================
 * OUTPUT 7: Combined design space
 *
 * Organoid count, minimum nodal force, Q-factor vs channel
 * length for glass and silicon sidewalls.
 * ============================================================ */

static void output_design_space(void)
{
    char fname[256];
    FILE *fp;
    int NL = 500;
    double L_arr[500];
    double r_glass, r_silicon;
    int iL;

    build_path(fname, sizeof(fname), "fig7_design_space.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_glass = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);
    linspace(1e-3, 25e-3, NL, L_arr);

    fprintf(fp, "L_mm,n_organoids,freq_kHz,"
                "Q_glass,Q_silicon,"
                "contrast_mid_dB_glass,contrast_mid_dB_silicon,"
                "F_min_norm_glass,F_min_norm_silicon,"
                "F_max_glass,F_max_silicon,"
                "T_coupling_741kHz\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_arr[iL];
        int n = (int)(L / 1.0e-3 + 0.5);
        if (n < 1) n = 1;
        double freq = n * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;
        double node_sp = L / (double)n;

        /* Q factor: Q ~ pi*n / (1 - r^2 + 2*alpha*L) */
        double Q_g = M_PI * n / (1.0 - r_glass*r_glass + 2.0*alpha*L);
        double Q_s = M_PI * n / (1.0 - r_silicon*r_silicon + 2.0*alpha*L);

        /* Contrast at midpoint */
        double S_g = sw_ratio(L/2.0, L, r_glass, alpha);
        double S_s = sw_ratio(L/2.0, L, r_silicon, alpha);
        double C_g = 10.0 * log10(dmax(contrast_from_swr(S_g), 1e-15));
        double C_s = 10.0 * log10(dmax(contrast_from_swr(S_s), 1e-15));

        /* Find min and max nodal force */
        double F_min_g = 1e30, F_max_g = 0;
        double F_min_s = 1e30, F_max_s = 0;
        int in;

        for (in = 0; in < n; in++) {
            double x = (in + 0.5) * node_sp;
            double Fg = fabs(radiation_force(x, L, freq, r_glass, alpha));
            double Fs = fabs(radiation_force(x, L, freq, r_silicon, alpha));
            if (Fg < F_min_g) F_min_g = Fg;
            if (Fg > F_max_g) F_max_g = Fg;
            if (Fs < F_min_s) F_min_s = Fs;
            if (Fs > F_max_s) F_max_s = Fs;
        }

        /* Coupling at 741 kHz, 10 um couplant */
        double T_coupl = coupling_T_power(741.5e3, 10e-6);

        fprintf(fp, "%.4f,%d,%.4f,%.2f,%.2f,%.4f,%.4f,"
                    "%.6f,%.6f,%.8e,%.8e,%.6f\n",
                L*1e3, n, freq/1e3,
                Q_g, Q_s, C_g, C_s,
                (F_max_g > 0) ? F_min_g/F_max_g : 0,
                (F_max_s > 0) ? F_min_s/F_max_s : 0,
                F_max_g, F_max_s,
                T_coupl);
    }

    fclose(fp);
    printf("  Written: %s (%d lengths)\n", fname, NL);
}



/* ============================================================
 * PART C: Dual-transducer lateral standing wave
 *
 * Two piezo transducers, one at each end of the channel (x=0, x=L).
 * Each launches a plane wave into the fluid. The standing wave forms
 * by superposition of the two counter-propagating waves, NOT by
 * wall reflections.
 *
 * p(x,t) = p_L * exp(i(kx - wt)) * exp(-alpha*x)
 *        + p_R * exp(i(-kx - wt + phi)) * exp(-alpha*(L-x))
 *
 * where p_L, p_R = pressure amplitudes from left/right transducers,
 * phi = relative phase between the two transducers.
 *
 * For equal-amplitude transducers (p_L = p_R = p0) at phi = 0:
 *
 *   |p(x)|^2 = p0^2 * [ exp(-2*alpha*x) + exp(-2*alpha*(L-x))
 *              + 2*exp(-alpha*L) * cos(2kx - kL) ]
 *
 * This gives a SYMMETRIC standing wave with perfect nodes at
 * x = (2m+1)*lambda/4 from the midpoint, regardless of wall
 * impedance. The wall reflection coefficient no longer matters
 * for the standing wave — it only causes secondary reflections
 * that add to the field.
 *
 * Advantages over single-source + passive wall:
 *   1. |r_eff| = 1.0 always — perfect contrast independent of wall Z
 *   2. Symmetric node positions — equal trapping at both ends
 *   3. Phase control — shift node positions by adjusting phi
 *   4. Can use ANY wall material (PDMS, glass, silicon, open)
 *
 * The phase parameter phi shifts the entire nodal pattern by
 *   delta_x = phi / (2k) = phi * lambda / (4*pi)
 *
 * For organoid positioning control, sweeping phi from 0 to pi
 * shifts all nodes by exactly lambda/4 = 0.5 mm — enough to
 * move organoids between electrodes.
 * ============================================================ */

/* Pressure at position x from dual transducers */
static Cpx pressure_at_dual(double x, double L, double freq,
                             double alpha, double phi,
                             double amp_ratio)
{
    double k = 2.0 * M_PI * freq / c_f;

    /* Left transducer: travels +x, attenuates with distance from x=0 */
    Cpx fwd = cpx_exp(cpx(-alpha * x, k * x));

    /* Right transducer: travels -x, attenuates with distance from x=L
     * phase offset phi, amplitude ratio amp_ratio (1.0 = equal) */
    Cpx rev = cpx_scale(amp_ratio,
                cpx_exp(cpx(-alpha * (L - x), -(k * x) + phi)));

    return cpx_add(fwd, rev);
}

/* Pressure squared for dual-transducer configuration */
static double pressure_sq_dual(double x, double L, double freq,
                                double alpha, double phi,
                                double amp_ratio)
{
    Cpx p = pressure_at_dual(x, L, freq, alpha, phi, amp_ratio);
    return p.re * p.re + p.im * p.im;
}

/* Radiation force for dual-transducer: -d/dx[<p^2>] */
static double radiation_force_dual(double x, double L, double freq,
                                    double alpha, double phi,
                                    double amp_ratio)
{
    double dx = L * 1e-5;
    double p2_plus, p2_minus;

    if (x - dx < 0) {
        p2_plus  = pressure_sq_dual(x + dx, L, freq, alpha, phi, amp_ratio);
        p2_minus = pressure_sq_dual(x,      L, freq, alpha, phi, amp_ratio);
        return -(p2_plus - p2_minus) / dx;
    }
    if (x + dx > L) {
        p2_plus  = pressure_sq_dual(x,      L, freq, alpha, phi, amp_ratio);
        p2_minus = pressure_sq_dual(x - dx, L, freq, alpha, phi, amp_ratio);
        return -(p2_plus - p2_minus) / dx;
    }

    p2_plus  = pressure_sq_dual(x + dx, L, freq, alpha, phi, amp_ratio);
    p2_minus = pressure_sq_dual(x - dx, L, freq, alpha, phi, amp_ratio);
    return -(p2_plus - p2_minus) / (2.0 * dx);
}

/* Gor'kov force for dual-transducer (physical units, N at p0=1Pa) */
static double gorkov_force_dual(double x, double L, double freq,
                                 double alpha, double phi,
                                 double amp_ratio)
{
    double F_raw = radiation_force_dual(x, L, freq, alpha, phi, amp_ratio);
    double a3 = a_organoid * a_organoid * a_organoid;
    double prefactor = Phi_contrast * M_PI * a3
                       / (3.0 * rho_f * c_f * c_f);
    return prefactor * F_raw;
}

/* OUTPUT 8: Dual-transducer pressure field comparison
 *
 * Side-by-side comparison of single-source (glass/silicon walls)
 * vs dual-transducer at three channel lengths.
 * Shows how dual-transducer achieves perfect contrast regardless
 * of wall material.
 */
static void output_dual_pressure_field(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 20e-3};
    int NL = 3;
    double r_glass, r_silicon;
    int Nx = 2000;
    int iL, ix;

    build_path(fname, sizeof(fname), "fig8_dual_pressure_field.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    r_glass  = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);

    fprintf(fp, "L_mm,x_mm,"
                "p_abs_single_glass,p_abs_single_silicon,"
                "p_abs_dual_phi0,p_abs_dual_phi45,p_abs_dual_phi90,"
                "p_sq_single_glass,p_sq_single_silicon,"
                "p_sq_dual_phi0,p_sq_dual_phi45,p_sq_dual_phi90\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double freq = n_mode * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;

        for (ix = 0; ix <= Nx; ix++) {
            double x = L * (double)ix / Nx;

            /* Single-source (existing model) */
            Cpx p_sg = pressure_at(x, L, freq, r_glass, alpha);
            Cpx p_ss = pressure_at(x, L, freq, r_silicon, alpha);

            /* Dual-transducer: phi = 0 (in phase) */
            Cpx p_d0 = pressure_at_dual(x, L, freq, alpha, 0.0, 1.0);
            /* Dual-transducer: phi = pi/4 (45 deg shift) */
            Cpx p_d45 = pressure_at_dual(x, L, freq, alpha,
                                          M_PI/4.0, 1.0);
            /* Dual-transducer: phi = pi/2 (90 deg shift) */
            Cpx p_d90 = pressure_at_dual(x, L, freq, alpha,
                                          M_PI/2.0, 1.0);

            fprintf(fp, "%.4f,%.6f,"
                        "%.8f,%.8f,"
                        "%.8f,%.8f,%.8f,"
                        "%.8f,%.8f,"
                        "%.8f,%.8f,%.8f\n",
                    L*1e3, x*1e3,
                    cpx_abs(p_sg), cpx_abs(p_ss),
                    cpx_abs(p_d0), cpx_abs(p_d45), cpx_abs(p_d90),
                    p_sg.re*p_sg.re + p_sg.im*p_sg.im,
                    p_ss.re*p_ss.re + p_ss.im*p_ss.im,
                    p_d0.re*p_d0.re + p_d0.im*p_d0.im,
                    p_d45.re*p_d45.re + p_d45.im*p_d45.im,
                    p_d90.re*p_d90.re + p_d90.im*p_d90.im);
        }
    }

    fclose(fp);
    printf("  Written: %s (3 lengths x %d positions)\n", fname, Nx+1);
}

/* OUTPUT 9: Phase sweep — node position control
 *
 * Sweep relative phase phi from 0 to pi for a 10 mm channel.
 * Track how pressure node positions shift with phase.
 * This is the key advantage of dual transducers: dynamic
 * repositioning of organoids without changing frequency.
 */
static void output_dual_phase_sweep(void)
{
    char fname[256];
    FILE *fp;
    double L = 10e-3;  /* 10 mm channel */
    int n_mode = 10;
    double freq = n_mode * c_f / (2.0 * L);
    double alpha = alpha_coeff * freq * freq;
    int Nphi = 100;
    int Nx = 1000;
    int iphi, ix;

    build_path(fname, sizeof(fname), "fig9_dual_phase_sweep.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    fprintf(fp, "phi_deg,x_mm,p_sq_dual\n");

    for (iphi = 0; iphi <= Nphi; iphi++) {
        double phi = M_PI * (double)iphi / Nphi;  /* 0 to pi */

        for (ix = 0; ix <= Nx; ix++) {
            double x = L * (double)ix / Nx;
            double p2 = pressure_sq_dual(x, L, freq, alpha, phi, 1.0);

            fprintf(fp, "%.4f,%.6f,%.8f\n",
                    phi * 180.0 / M_PI, x*1e3, p2);
        }
    }

    fclose(fp);
    printf("  Written: %s (%d phases x %d positions)\n",
           fname, Nphi+1, Nx+1);
}

/* OUTPUT 10: Amplitude imbalance effect
 *
 * What happens when the two transducers have unequal amplitudes?
 * Sweep amp_ratio from 0.5 to 1.0 (50% to 100% match).
 * Shows how contrast degrades with imbalance — important for
 * practical implementation where two PZTs may not match exactly.
 */
static void output_dual_amplitude_sweep(void)
{
    char fname[256];
    FILE *fp;
    double L = 10e-3;
    int n_mode = 10;
    double freq = n_mode * c_f / (2.0 * L);
    double alpha = alpha_coeff * freq * freq;
    int Namp = 50;
    int Nx = 1000;
    int iamp, ix;

    build_path(fname, sizeof(fname), "fig10_dual_amplitude_sweep.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    fprintf(fp, "amp_ratio,x_mm,p_sq_dual,p_abs_dual\n");

    for (iamp = 0; iamp <= Namp; iamp++) {
        double amp_ratio = 0.5 + 0.5 * (double)iamp / Namp;  /* 0.5 to 1.0 */

        for (ix = 0; ix <= Nx; ix++) {
            double x = L * (double)ix / Nx;
            double p2 = pressure_sq_dual(x, L, freq, alpha, 0.0, amp_ratio);

            fprintf(fp, "%.4f,%.6f,%.8f,%.8f\n",
                    amp_ratio, x*1e3, p2, sqrt(p2));
        }
    }

    fclose(fp);
    printf("  Written: %s (%d ratios x %d positions)\n",
           fname, Namp+1, Nx+1);
}

/* OUTPUT 11: Dual-transducer nodal uniformity + Gor'kov force
 *
 * Same analysis as fig4 but for dual-transducer configuration.
 * Compare trapping force uniformity: single-source vs dual.
 */
static void output_dual_nodal_uniformity(void)
{
    char fname[256];
    FILE *fp;
    double L_vals[] = {5e-3, 10e-3, 20e-3};
    int NL = 3;
    int iL, in_node;

    build_path(fname, sizeof(fname), "fig11_dual_nodal_uniformity.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    fprintf(fp, "L_mm,node_index,x_mm,"
                "p_sq_dual,F_dual,F_dual_norm,F_gorkov_dual_N\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        int n_mode = (int)(L / 1.0e-3 + 0.5);
        double freq = n_mode * c_f / (2.0 * L);
        double alpha = alpha_coeff * freq * freq;
        double node_spacing = L / (double)n_mode;
        double F_max = 0;
        double *F_arr, *p2_arr, *x_arr, *Fg_arr;

        F_arr  = (double *)malloc(n_mode * sizeof(double));
        p2_arr = (double *)malloc(n_mode * sizeof(double));
        x_arr  = (double *)malloc(n_mode * sizeof(double));
        Fg_arr = (double *)malloc(n_mode * sizeof(double));

        for (in_node = 0; in_node < n_mode; in_node++) {
            double x = (in_node + 0.5) * node_spacing;
            x_arr[in_node]  = x;
            p2_arr[in_node] = pressure_sq_dual(x, L, freq, alpha, 0.0, 1.0);
            F_arr[in_node]  = fabs(radiation_force_dual(x, L, freq,
                                                         alpha, 0.0, 1.0));
            Fg_arr[in_node] = fabs(gorkov_force_dual(x, L, freq,
                                                      alpha, 0.0, 1.0));
            if (F_arr[in_node] > F_max) F_max = F_arr[in_node];
        }

        for (in_node = 0; in_node < n_mode; in_node++) {
            fprintf(fp, "%.4f,%d,%.6f,%.8f,%.8e,%.6f,%.8e\n",
                    L*1e3, in_node, x_arr[in_node]*1e3,
                    p2_arr[in_node], F_arr[in_node],
                    (F_max > 0) ? F_arr[in_node]/F_max : 0,
                    Fg_arr[in_node]);
        }

        free(F_arr); free(p2_arr); free(x_arr); free(Fg_arr);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * PART D: PDMS wall transmission (immersed dual-transducer)
 *
 * For the immersed configuration, each PZT radiates through:
 *   water gap → PDMS wall (thick) → water channel
 *
 * The PDMS wall acts as an impedance-matched but LOSSY layer.
 * Z_PDMS ≈ 1.08 MRayl vs Z_water = 1.50 MRayl, so |r| ≈ 0.16
 * at the water-PDMS interface — most energy enters the PDMS.
 * But PDMS attenuates ~5 dB/cm/MHz, so a 10 mm wall at 741 kHz
 * absorbs ~3.7 dB per pass.
 *
 * The transfer matrix for a lossy PDMS layer between two water
 * half-spaces gives the power transmission coefficient.
 *
 * The glass slide (borosilicate) underneath the channel is
 * perpendicular to the lateral wave propagation and does not
 * directly affect the lateral standing wave model.
 *
 * Geometry (from PI's sketch):
 *   PZT — 5 mm water gap — 10 mm PDMS wall — 20 mm water channel
 *         — 10 mm PDMS wall — 5 mm water gap — PZT
 *   Channel depth: 1 mm (water above glass slide)
 *   PDMS ceiling: 9 mm above channel
 *   Glass slide: borosilicate, underneath everything
 * ============================================================ */

/* Transmission through water → PDMS(d, lossy) → water
 * Uses transfer matrix with complex wavenumber for PDMS.
 *
 * For a layer of impedance Z2, complex wavenumber k2 = w/c2 + i*alpha,
 * between two half-spaces of impedance Z1:
 *
 *   T_power = 1 / |cos(k2*d) + (i/2)(Z2/Z1 + Z1/Z2)*sin(k2*d)|^2
 *
 * When Z1 ≈ Z2 (PDMS ≈ water), this simplifies to:
 *   T_power ≈ exp(-2*alpha*d)   (attenuation-dominated)
 */
static double pdms_transmission(double freq, double d_pdms)
{
    double omega = 2.0 * M_PI * freq;
    double k_real = omega / c_pdms;

    /* PDMS attenuation: convert dB/cm/MHz to Np/m */
    double f_MHz = freq / 1e6;
    double alpha_dBcm = alpha_pdms_per_MHz_dBcm * f_MHz;
    double alpha_Npm = alpha_dBcm * 100.0 / 8.686;  /* dB/cm → Np/m */

    /* Complex wavenumber in PDMS */
    Cpx k2 = cpx(k_real, alpha_Npm);
    Cpx k2d = cpx_scale(d_pdms, k2);

    /* Transfer matrix: T = 1/|M11|^2
     * M11 = cos(k2d) + (i/2)(Z_pdms/Z_f + Z_f/Z_pdms)*sin(k2d) */
    Cpx cos_k2d = cpx_cos(k2d);
    Cpx sin_k2d = cpx_sin(k2d);

    double Z_ratio = 0.5 * (Z_pdms / Z_f + Z_f / Z_pdms);
    Cpx i_Z_sin = cpx_mul(cpx(0, Z_ratio), sin_k2d);
    Cpx M11 = cpx_add(cos_k2d, i_Z_sin);

    double M11_sq = M11.re * M11.re + M11.im * M11.im;
    return 1.0 / M11_sq;
}

/* Effective amplitude reaching channel from one transducer:
 * passes through water gap (lossless) + PDMS wall (lossy) */
static double pdms_amplitude_factor(double freq, double d_pdms)
{
    return sqrt(pdms_transmission(freq, d_pdms));
}

/* OUTPUT 12: PDMS wall transmission vs thickness
 *
 * Sweep PDMS wall thickness from 0.5 mm to 15 mm at 741 kHz.
 * This answers the critical question: how thin must the PDMS
 * walls be for adequate acoustic transmission?
 */
static void output_pdms_transmission(void)
{
    char fname[256];
    FILE *fp;
    double freq = 741.5e3;
    int Nd = 500;
    double d_arr[500];
    int i;
    double lambda_pdms = c_pdms / freq;

    build_path(fname, sizeof(fname), "fig12_pdms_transmission.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(0.1e-3, 15e-3, Nd, d_arr);

    fprintf(fp, "d_mm,d_over_lambda,T_power,T_power_dB,"
                "amplitude_factor,alpha_Npm,loss_per_wall_dB\n");

    for (i = 0; i < Nd; i++) {
        double T = pdms_transmission(freq, d_arr[i]);
        double T_dB = 10.0 * log10(dmax(T, 1e-30));
        double amp = pdms_amplitude_factor(freq, d_arr[i]);

        /* For reference: simple exponential loss */
        double f_MHz = freq / 1e6;
        double alpha_Npm = alpha_pdms_per_MHz_dBcm * f_MHz * 100.0 / 8.686;
        double simple_loss_dB = 20.0 * log10(exp(-alpha_Npm * d_arr[i]));

        fprintf(fp, "%.4f,%.6f,%.8f,%.4f,%.6f,%.4f,%.4f\n",
                d_arr[i]*1e3, d_arr[i]/lambda_pdms,
                T, T_dB, amp,
                alpha_Npm, simple_loss_dB);
    }

    fclose(fp);
    printf("  Written: %s (%d thicknesses, lambda_pdms=%.3f mm)\n",
           fname, Nd, lambda_pdms*1e3);
}

/* OUTPUT 13: Effective dual-transducer contrast through PDMS walls
 *
 * For each PDMS wall thickness, compute:
 * - Amplitude reaching channel from each side (after PDMS loss)
 * - Effective standing wave contrast inside channel
 * - Comparison with minimum contrast needed for organoid trapping
 *
 * Sweep thickness and also sweep attenuation coefficient to show
 * sensitivity to PDMS cure conditions.
 */
static void output_pdms_contrast_vs_thickness(void)
{
    char fname[256];
    FILE *fp;
    double freq = 741.5e3;
    double L_channel = 20e-3;  /* 20 mm channel */
    int n_mode = 20;
    double alpha_water = alpha_coeff * freq * freq;
    int Nd = 200;
    double d_arr[200];
    int id;

    build_path(fname, sizeof(fname), "fig13_pdms_contrast.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(0.5e-3, 15e-3, Nd, d_arr);

    fprintf(fp, "d_pdms_mm,T_one_wall,T_two_walls,T_two_walls_dB,"
                "amp_L,amp_R,effective_contrast,effective_contrast_dB,"
                "min_force_norm\n");

    for (id = 0; id < Nd; id++) {
        double T1 = pdms_transmission(freq, d_arr[id]);
        double T2 = T1 * T1;  /* two walls total */
        double T2_dB = 10.0 * log10(dmax(T2, 1e-30));

        /* Amplitude from each side after passing through one PDMS wall */
        double amp = sqrt(T1);

        /* Both transducers see the same PDMS wall thickness,
         * so amp_L = amp_R = amp (symmetric).
         * The standing wave inside the channel has equal amplitudes
         * from both sides → perfect contrast despite PDMS loss.
         *
         * The PDMS loss reduces the ABSOLUTE force, not the contrast.
         * Contrast remains infinite for equal paths. */

        /* Compute p^2 field to find actual contrast */
        double p2_max = 0, p2_min = 1e30;
        int ix;
        for (ix = 0; ix <= 1000; ix++) {
            double x = L_channel * (double)ix / 1000;
            /* Both amplitudes reduced equally by PDMS */
            Cpx p_fwd = cpx_scale(amp,
                cpx_exp(cpx(-alpha_water * x,
                            2.0 * M_PI * freq / c_f * x)));
            Cpx p_rev = cpx_scale(amp,
                cpx_exp(cpx(-alpha_water * (L_channel - x),
                            -(2.0 * M_PI * freq / c_f * x))));
            Cpx p_total = cpx_add(p_fwd, p_rev);
            double p2 = p_total.re * p_total.re + p_total.im * p_total.im;
            if (p2 > p2_max) p2_max = p2;
            if (p2 < p2_min) p2_min = p2;
        }

        double contrast = (p2_min > 1e-20) ? p2_max / p2_min : 1e10;
        double contrast_dB = 10.0 * log10(dmax(contrast, 1e-15));

        /* Normalized force: how much force remains compared to
         * no-PDMS case (amp=1). Force ∝ amp^2 ∝ T */
        double force_norm = T1;

        fprintf(fp, "%.4f,%.6f,%.8f,%.4f,%.6f,%.6f,%.2f,%.2f,%.6f\n",
                d_arr[id]*1e3, T1, T2, T2_dB,
                amp, amp,
                contrast > 1e6 ? 1e6 : contrast,
                contrast_dB > 100 ? 100 : contrast_dB,
                force_norm);
    }

    fclose(fp);
    printf("  Written: %s (%d thicknesses)\n", fname, Nd);
}

/* OUTPUT 14: PDMS attenuation sensitivity
 *
 * The PDMS attenuation coefficient varies significantly with
 * curing conditions (ratio, temperature, time). Sweep alpha
 * from 1 to 15 dB/cm/MHz at four wall thicknesses:
 * 10 mm (current), 5 mm, 3 mm, and 1 mm.
 */
static void output_pdms_attenuation_sensitivity(void)
{
    char fname[256];
    FILE *fp;
    double freq = 741.5e3;
    double f_MHz = freq / 1e6;
    int Na = 100;
    int i;
    double d_10 = 10e-3;
    double d_5  = 5e-3;
    double d_3  = 3e-3;
    double d_1  = 1e-3;

    build_path(fname, sizeof(fname), "fig14_pdms_attenuation_sensitivity.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    fprintf(fp, "alpha_dBcm_MHz,alpha_Npm,"
                "T_10mm,T_10mm_dB,force_norm_10mm,"
                "T_5mm,T_5mm_dB,force_norm_5mm,"
                "T_3mm,T_3mm_dB,force_norm_3mm,"
                "T_1mm,T_1mm_dB,force_norm_1mm\n");

    for (i = 0; i < Na; i++) {
        double a_dBcm = 1.0 + 14.0 * (double)i / (Na - 1);  /* 1-15 */
        double a_Npm = a_dBcm * f_MHz * 100.0 / 8.686;

        /* Simple exponential: T ≈ exp(-2*alpha*d) for Z-matched case */
        double T10 = exp(-2.0 * a_Npm * d_10);
        double T5  = exp(-2.0 * a_Npm * d_5);
        double T3  = exp(-2.0 * a_Npm * d_3);
        double T1  = exp(-2.0 * a_Npm * d_1);

        fprintf(fp, "%.4f,%.4f,"
                    "%.8f,%.4f,%.6f,"
                    "%.8f,%.4f,%.6f,"
                    "%.8f,%.4f,%.6f,"
                    "%.8f,%.4f,%.6f\n",
                a_dBcm, a_Npm,
                T10, 10.0*log10(dmax(T10, 1e-30)), T10,
                T5,  10.0*log10(dmax(T5,  1e-30)), T5,
                T3,  10.0*log10(dmax(T3,  1e-30)), T3,
                T1,  10.0*log10(dmax(T1,  1e-30)), T1);
    }

    fclose(fp);
    printf("  Written: %s (alpha sweep at d=10,5,3,1 mm)\n", fname);
}

/* OUTPUT 16: PDMS transmission vs. frequency
 *
 * The student correctly noted that frequency is the true independent
 * variable in an experiment — you set f on the function generator,
 * and the PDMS attenuation follows as alpha = alpha_coeff * f.
 *
 * This sweeps frequency from 100 kHz to 5 MHz at four PDMS wall
 * thicknesses (10, 5, 3, 1 mm), using the frequency-dependent
 * attenuation alpha(f) = alpha_dBcm_MHz * f_MHz * 100 / 8.686 Np/m.
 *
 * Shows how transmission degrades at higher frequencies and why
 * our 741 kHz operating frequency is favorable for PDMS transmission.
 */
static void output_pdms_vs_frequency(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 500;
    int i;
    double d_10 = 10e-3, d_5 = 5e-3, d_3 = 3e-3, d_1 = 1e-3;

    build_path(fname, sizeof(fname), "fig16_pdms_vs_frequency.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    fprintf(fp, "freq_kHz,freq_MHz,"
                "T_10mm,T_10mm_dB,force_norm_10mm,"
                "T_5mm,T_5mm_dB,force_norm_5mm,"
                "T_3mm,T_3mm_dB,force_norm_3mm,"
                "T_1mm,T_1mm_dB,force_norm_1mm\n");

    for (i = 0; i < Nf; i++) {
        /* Frequency from 100 kHz to 5 MHz */
        double f = 100e3 + (5e6 - 100e3) * (double)i / (Nf - 1);
        double f_MHz = f / 1e6;

        /* Frequency-dependent attenuation in PDMS */
        double alpha_Npm = alpha_pdms_per_MHz_dBcm * f_MHz * 100.0 / 8.686;

        /* Exponential transmission: T = exp(-2*alpha*d) */
        double T10 = exp(-2.0 * alpha_Npm * d_10);
        double T5  = exp(-2.0 * alpha_Npm * d_5);
        double T3  = exp(-2.0 * alpha_Npm * d_3);
        double T1  = exp(-2.0 * alpha_Npm * d_1);

        fprintf(fp, "%.4f,%.6f,"
                    "%.8f,%.4f,%.6f,"
                    "%.8f,%.4f,%.6f,"
                    "%.8f,%.4f,%.6f,"
                    "%.8f,%.4f,%.6f\n",
                f/1e3, f_MHz,
                T10, 10.0*log10(dmax(T10, 1e-30)), T10,
                T5,  10.0*log10(dmax(T5,  1e-30)), T5,
                T3,  10.0*log10(dmax(T3,  1e-30)), T3,
                T1,  10.0*log10(dmax(T1,  1e-30)), T1);
    }

    fclose(fp);
    printf("  Written: %s (%d frequencies, 100 kHz to 5 MHz)\n", fname, Nf);
}

/* ============================================================
 * PART E: Glass slide vertical confinement
 *
 * The lateral standing wave propagates horizontally in the water
 * channel. The glass slide underneath acts as the bottom boundary.
 * Whether acoustic energy leaks into the glass depends on the
 * angle of incidence at the water-glass interface.
 *
 * Snell's law: sin(theta_t)/c_glass = sin(theta_i)/c_water
 * Critical angle: theta_c = arcsin(c_water/c_glass)
 *
 * For theta_i > theta_c: total internal reflection (no leakage)
 * For theta_i < theta_c: partial transmission into glass
 *
 * The lateral standing wave has theta_i = 90° (grazing), which
 * is far beyond theta_c ~ 15° for glass. So the glass acts
 * as a perfect reflector for the lateral mode.
 * ============================================================ */

/* OUTPUT 15: Glass slide vertical confinement
 *
 * Power reflection coefficient R vs. incidence angle (0-90 deg)
 * for water-glass, water-PDMS, water-acrylic, and open water.
 * Shows total internal reflection beyond critical angle for
 * glass and acrylic, but NOT for PDMS (c_pdms < c_water).
 */
static void output_glass_confinement(void)
{
    char fname[256];
    FILE *fp;
    int Na = 900;
    int i;
    double theta_c_glass, theta_c_acrylic;

    build_path(fname, sizeof(fname), "fig15_glass_confinement.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    theta_c_glass   = asin(c_f / c_glass) * 180.0 / M_PI;
    theta_c_acrylic = asin(c_f / c_acrylic) * 180.0 / M_PI;

    fprintf(fp, "theta_deg,"
                "R_glass,R_glass_dB,T_glass,"
                "R_pdms,R_pdms_dB,T_pdms,"
                "R_acrylic,R_acrylic_dB,T_acrylic,"
                "R_open,R_open_dB,T_open\n");

    for (i = 0; i <= Na; i++) {
        double theta_deg = 90.0 * (double)i / Na;
        double theta_rad = theta_deg * M_PI / 180.0;
        double sin_t = sin(theta_rad);
        double cos_t = cos(theta_rad);

        /* Glass: c_glass > c_f → critical angle exists */
        double R_glass, T_glass;
        {
            double sin_t2 = sin_t * c_glass / c_f;
            if (sin_t2 >= 1.0) {
                R_glass = 1.0;
            } else {
                double cos_t2 = sqrt(1.0 - sin_t2 * sin_t2);
                double Z_eff_g = Z_glass / cos_t2;
                double Z_eff_w = Z_f / cos_t;
                double r_val = (Z_eff_g - Z_eff_w) / (Z_eff_g + Z_eff_w);
                R_glass = r_val * r_val;
            }
            T_glass = 1.0 - R_glass;
        }

        /* PDMS: c_pdms < c_f → NO critical angle, energy leaks always */
        double R_pdms, T_pdms;
        {
            double sin_t2 = sin_t * c_pdms / c_f;
            double cos_t2 = sqrt(1.0 - sin_t2 * sin_t2);
            double Z_eff_p = Z_pdms / cos_t2;
            double Z_eff_w = (cos_t > 1e-10) ? Z_f / cos_t : Z_f * 1e10;
            double r_val = (Z_eff_p - Z_eff_w) / (Z_eff_p + Z_eff_w);
            R_pdms = r_val * r_val;
            T_pdms = 1.0 - R_pdms;
        }

        /* Acrylic: c_acrylic > c_f → critical angle exists */
        double R_acrylic, T_acrylic;
        {
            double sin_t2 = sin_t * c_acrylic / c_f;
            if (sin_t2 >= 1.0) {
                R_acrylic = 1.0;
            } else {
                double cos_t2 = sqrt(1.0 - sin_t2 * sin_t2);
                double Z_eff_a = Z_acrylic / cos_t2;
                double Z_eff_w = Z_f / cos_t;
                double r_val = (Z_eff_a - Z_eff_w) / (Z_eff_a + Z_eff_w);
                R_acrylic = r_val * r_val;
            }
            T_acrylic = 1.0 - R_acrylic;
        }

        /* Open water: no boundary → R = 0 everywhere */
        double R_open = 0.0, T_open = 1.0;

        fprintf(fp, "%.2f,"
                    "%.8f,%.4f,%.8f,"
                    "%.8f,%.4f,%.8f,"
                    "%.8f,%.4f,%.8f,"
                    "%.8f,%.4f,%.8f\n",
                theta_deg,
                R_glass, 10.0*log10(dmax(R_glass, 1e-30)), T_glass,
                R_pdms,  10.0*log10(dmax(R_pdms, 1e-30)),  T_pdms,
                R_acrylic, 10.0*log10(dmax(R_acrylic, 1e-30)), T_acrylic,
                R_open, 10.0*log10(dmax(R_open, 1e-30)), T_open);
    }

    fclose(fp);
    printf("  Written: %s (%d angles)\n", fname, Na+1);
    printf("           theta_c(glass)   = %.1f deg\n", theta_c_glass);
    printf("           theta_c(acrylic) = %.1f deg\n", theta_c_acrylic);
    printf("           theta_c(PDMS)    = NONE (c_pdms < c_water)\n");
}

/* ============================================================
 * MAIN
 * ============================================================ */

int main(void)
{
    double r_glass, r_silicon;
    double f_target, lambda_f, lambda_2;
    int mkdir_result;

    /* Create output directory */
    mkdir_result = MKDIR(OUTPUT_DIR);
    if (mkdir_result != 0) {
        if (errno != EEXIST) {
            fprintf(stderr, "Error: cannot create %s\n", OUTPUT_DIR);
            return 1;
        }
    }

    /* Compute derived impedances */
    Z_f     = rho_f * c_f;
    Z_glass = rho_glass * c_glass;
    Z_si    = rho_si * c_si;
    Z_acrylic = rho_acrylic * c_acrylic;
    Z_pzt   = rho_pzt * c_pzt;
    Z_pdms  = rho_pdms * c_pdms;
    Z_1     = rho_1 * c_1;
    Z_2     = rho_2 * c_2;
    Z_3     = rho_3 * c_3;

    /* BEGIN ADDITION: Compute Gor'kov contrast factor
     *
     * f1 (monopole) = 1 - kappa_p/kappa_f
     *               = 1 - (rho_f * c_f^2) / (rho_p * c_p^2)
     *
     * f2 (dipole)   = 2*(rho_p - rho_f) / (2*rho_p + rho_f)
     *
     * Phi = f1/3 + f2/2
     *
     * For typical cells in water (rho_p~1050, c_p~1550):
     *   f1 ~ 1 - (1000*1483^2)/(1050*1550^2) ~ 0.13
     *   f2 ~ 2*(1050-1000)/(2*1050+1000)     ~ 0.032
     *   Phi ~ 0.13/3 + 0.032/2               ~ 0.060
     *
     * Note: literature often cites Phi ~ 0.17 for polystyrene
     * beads (rho=1050, c=2350). For soft cells, Phi is lower
     * because c_p is closer to c_f. This matters! */
    {
        double kappa_f = 1.0 / (rho_f * c_f * c_f);   /* fluid compressibility */
        double kappa_p = 1.0 / (rho_p * c_p * c_p);   /* particle compressibility */
        f1_monopole = 1.0 - kappa_p / kappa_f;
        f2_dipole   = 2.0 * (rho_p - rho_f) / (2.0 * rho_p + rho_f);
        Phi_contrast = f1_monopole / 3.0 + f2_dipole / 2.0;
    }
    /* END ADDITION: Gor'kov contrast factor */

    r_glass  = reflection_coeff(Z_glass, Z_f);
    r_silicon = reflection_coeff(Z_si, Z_f);

    f_target = c_f / (2.0 * 1.0e-3);  /* f for lambda/2 = 1 mm */
    lambda_f = c_f / f_target;
    lambda_2 = c_2 / f_target;

    printf("=====================================================\n");
    printf(" BENG207_2: Rigid Cavity Resonator Model\n");
    printf(" Lateral standing wave for organoid formation\n");
    printf("=====================================================\n");
    printf(" Fluid: water\n");
    printf("   Z_f           = %.2f MRayl\n", Z_f/1e6);
    printf("   c_f           = %.0f m/s\n", c_f);
    printf("-----------------------------------------------------\n");
    printf(" Target organoid spacing = 1 mm (CI electrode match)\n");
    printf("   f_target      = %.1f kHz\n", f_target/1e3);
    printf("   lambda_f      = %.3f mm\n", lambda_f*1e3);
    printf("   lambda/2      = %.3f mm\n", lambda_f*1e3/2.0);
    printf("-----------------------------------------------------\n");
    printf(" Wall reflection coefficients:\n");
    printf("   Glass:   Z = %.2f MRayl  |r| = %.4f\n",
           Z_glass/1e6, r_glass);
    printf("   Silicon: Z = %.2f MRayl  |r| = %.4f\n",
           Z_si/1e6, r_silicon);
    printf("-----------------------------------------------------\n");
    printf(" Vertical coupling (at %.0f kHz):\n", f_target/1e3);
    printf("   Z_LiNbO3     = %.2f MRayl\n", Z_1/1e6);
    printf("   Z_couplant   = %.2f MRayl\n", Z_2/1e6);
    printf("   Z_glass      = %.2f MRayl\n", Z_3/1e6);
    printf("   lambda_coupl = %.3f mm (d/lambda < 0.025)\n", lambda_2*1e3);
    printf("-----------------------------------------------------\n");
    printf(" Channel lengths: 5-20 mm -> 5-20 organoids\n");
    printf(" Human cochlea ~30 mm, CI = 22 electrodes\n");
    /* BEGIN ADDITION: Gor'kov summary + BENG207_3 hand-off */
    printf("-----------------------------------------------------\n");
    printf(" Gor'kov acoustic contrast factor:\n");
    printf("   rho_p  = %.0f kg/m^3 (cell density)\n", rho_p);
    printf("   c_p    = %.0f m/s (cell sound speed)\n", c_p);
    printf("   a      = %.0f um (organoid radius)\n", a_organoid*1e6);
    printf("   f1     = %.4f (monopole — compressibility contrast,\n",
           f1_monopole);
    printf("              cells less compressible than water: c_p > c_f)\n");
    printf("   f2     = %.4f (dipole — density contrast,\n", f2_dipole);
    printf("              cells slightly denser: rho_p > rho_f)\n");
    printf("   Phi    = %.4f (acoustic contrast factor)\n", Phi_contrast);
    printf("   Phi > 0 => organoids migrate to pressure NODES\n");
    printf("-----------------------------------------------------\n");
    printf(" Hand-off to BENG207_3 (stress estimation at p_ac = 500 kPa):\n");
    {
        double p_ref = 0.5e6;  /* 500 kPa reference pressure */
        double f_lat = 741.5e3;
        double f_vert = 20.0e6;
        double D_org = 2.0 * a_organoid;
        double lambda_vert = c_f / f_vert;

        /* Gorkov: sigma = 4 * Phi * k * a * E_ac
         *   where E_ac = p^2 / (4 * rho_f * c_f^2)               */
        double k_lat  = 2.0 * M_PI * f_lat / c_f;
        double k_vert = 2.0 * M_PI * f_vert / c_f;
        double E_ac   = p_ref * p_ref / (4.0 * rho_f * c_f * c_f);

        double sig_741k = 4.0 * Phi_contrast * k_lat  * a_organoid * E_ac;
        double sig_20M  = 4.0 * Phi_contrast * k_vert * a_organoid * E_ac;
        double g_corr   = 0.5;
        double sig_used = g_corr * sig_20M;

        printf("   sigma_gorkov(741 kHz)  = %.2f Pa  (D/lambda=%.2f, "
               "Gorkov valid)\n",
               sig_741k, D_org / (c_f / f_lat));
        printf("   sigma_gorkov(20 MHz)   = %.2f Pa  (raw, no g_corr)\n",
               sig_20M);
        printf("   sigma_0 -> BENG207_3   = %.2f Pa  (g_corr=%.1f for "
               "D/lambda=%.1f)\n",
               sig_used, g_corr, D_org / lambda_vert);
        printf("   [replaces old hand-waved alpha=0.15 -> 11.4 Pa]\n");
    }
    /* END ADDITION */
    printf("=====================================================\n\n");

    printf("Generating Part A outputs (lateral cavity)...\n\n");
    output_pressure_field();
    output_contrast_map();
    output_mode_structure();
    output_nodal_uniformity();

    printf("\nGenerating Part B output (vertical coupling)...\n\n");
    output_coupling_layer();

    printf("\nGenerating combined outputs...\n\n");
    output_frequency_sensitivity();
    output_design_space();

    printf("\nDone. 7 CSV files generated in %s/\n", OUTPUT_DIR);

    /* ============================================================
     * PART C: Dual-transducer outputs
     * ============================================================ */
    printf("\nGenerating Part C outputs (dual-transducer)...\n\n");
    {
        double L_ref = 10e-3;
        int n_ref = 10;
        double f_ref = n_ref * c_f / (2.0 * L_ref);
        double alpha_ref = alpha_coeff * f_ref * f_ref;
        double lambda_ref = c_f / f_ref;

        double S_mid_g = sw_ratio(L_ref/2.0, L_ref, r_glass, alpha_ref);
        double C_single_g = contrast_from_swr(S_mid_g);
        double r_acrylic = reflection_coeff(Z_acrylic, Z_f);
        double S_mid_a = sw_ratio(L_ref/2.0, L_ref, r_acrylic, alpha_ref);
        double C_single_a = contrast_from_swr(S_mid_a);
        double r_pdms_wall = reflection_coeff(Z_pdms, Z_f);

        /* Find actual dual contrast */
        double p2_max_d = 0, p2_min_d = 1e30;
        int ix;
        for (ix = 0; ix <= 2000; ix++) {
            double x = L_ref * (double)ix / 2000;
            double p2 = pressure_sq_dual(x, L_ref, f_ref, alpha_ref,
                                          0.0, 1.0);
            if (p2 > p2_max_d) p2_max_d = p2;
            if (p2 < p2_min_d) p2_min_d = p2;
        }
        double C_dual = (p2_min_d > 1e-20)
                         ? p2_max_d / p2_min_d : 1e10;

        printf(" Dual-transducer configuration (Part C):\n");
        printf(" Ref: Hou et al., Extreme Mech. Lett. 37, 100716 (2020)\n");
        printf("-----------------------------------------------------\n");
        printf("   Two PZT-5A transducers, immersed in water\n");
        printf("   Waves pass through PDMS walls into channel\n");
        printf("   Standing wave by counter-propagating waves\n");
        printf("-----------------------------------------------------\n");
        printf("   Material impedances:\n");
        printf("     Z_pzt     = %.2f MRayl (PZT-5A)\n", Z_pzt/1e6);
        printf("     Z_acrylic = %.2f MRayl (PMMA)\n", Z_acrylic/1e6);
        printf("     Z_pdms    = %.2f MRayl (Sylgard 184)\n", Z_pdms/1e6);
        printf("     Z_water   = %.2f MRayl\n", Z_f/1e6);
        printf("     |r| water-PDMS    = %.4f (nearly matched)\n",
               r_pdms_wall);
        printf("     |r| water-acrylic = %.4f\n", r_acrylic);
        printf("-----------------------------------------------------\n");
        printf("   Contrast comparison (L=10mm):\n");
        printf("     Config A: single + glass wall:    %.1f (%.1f dB)\n",
               C_single_g, 10.0*log10(dmax(C_single_g, 1e-15)));
        printf("     Config A: single + acrylic wall:  %.1f (%.1f dB)\n",
               C_single_a, 10.0*log10(dmax(C_single_a, 1e-15)));
        if (C_dual > 1e6)
            printf("     Config B: dual PZT (phi=0):       infinite\n");
        else
            printf("     Config B: dual PZT (phi=0):       %.1f (%.1f dB)\n",
                   C_dual, 10.0*log10(dmax(C_dual, 1e-15)));
        printf("   Node shift range:\n");
        printf("     phi = 0 to pi  =>  dx = 0 to %.3f mm\n",
               lambda_ref * 1e3 / 4.0);
        printf("-----------------------------------------------------\n\n");
    }
    output_dual_pressure_field();
    output_dual_phase_sweep();
    output_dual_amplitude_sweep();
    output_dual_nodal_uniformity();

    /* ============================================================
     * PART D: PDMS wall transmission
     * ============================================================ */
    printf("\nGenerating Part D outputs (PDMS wall transmission)...\n\n");
    {
        double freq_d = 741.5e3;
        double d_10mm = 10e-3;
        double d_5mm  = 5e-3;
        double d_3mm  = 3e-3;
        double d_1mm  = 1e-3;
        double T_10 = pdms_transmission(freq_d, d_10mm);
        double T_5  = pdms_transmission(freq_d, d_5mm);
        double T_3  = pdms_transmission(freq_d, d_3mm);
        double T_1  = pdms_transmission(freq_d, d_1mm);
        double lambda_pdms = c_pdms / freq_d;

        printf(" PDMS wall transmission at 741 kHz:\n");
        printf("-----------------------------------------------------\n");
        printf("   Z_pdms    = %.2f MRayl\n", Z_pdms/1e6);
        printf("   |r| water-PDMS = %.4f\n",
               reflection_coeff(Z_pdms, Z_f));
        printf("   lambda    = %.2f mm\n", lambda_pdms*1e3);
        printf("   alpha     = %.1f dB/cm/MHz\n",
               alpha_pdms_per_MHz_dBcm);
        printf("-----------------------------------------------------\n");
        printf("   d = 10 mm (%.1f lambda) — YOUR CURRENT DESIGN:\n",
               d_10mm/lambda_pdms);
        printf("     T_power per wall = %.1f%%\n", T_10*100);
        printf("   d = 5 mm  (%.1f lambda):\n", d_5mm/lambda_pdms);
        printf("     T_power per wall = %.1f%%\n", T_5*100);
        printf("   d = 3 mm  (%.1f lambda):\n", d_3mm/lambda_pdms);
        printf("     T_power per wall = %.1f%%\n", T_3*100);
        printf("   d = 1 mm  (%.1f lambda):\n", d_1mm/lambda_pdms);
        printf("     T_power per wall = %.1f%%\n", T_1*100);
        printf("-----------------------------------------------------\n");
        if (T_10 < 0.5) {
            printf("   *** WARNING: 10 mm PDMS loses >%.0f%% per wall ***\n",
                   (1.0 - T_10)*100);
            printf("   *** Consider thinning to 3 mm or less        ***\n");
        }
        printf("=====================================================\n\n");
    }
    output_pdms_transmission();
    output_pdms_contrast_vs_thickness();
    output_pdms_attenuation_sensitivity();
    output_pdms_vs_frequency();

    /* ============================================================
     * PART E: Glass slide vertical confinement
     * ============================================================ */
    printf("\nGenerating Part E output (glass slide confinement)...\n\n");
    {
        double theta_c = asin(c_f / c_glass) * 180.0 / M_PI;
        printf(" Glass slide as bottom boundary:\n");
        printf("-----------------------------------------------------\n");
        printf("   Lateral wave propagates at theta_i = 90 deg\n");
        printf("   Critical angle (glass):   %.1f deg\n", theta_c);
        printf("   Critical angle (acrylic): %.1f deg\n",
               asin(c_f / c_acrylic) * 180.0 / M_PI);
        printf("   Critical angle (PDMS):    NONE (c_pdms < c_water)\n");
        printf("-----------------------------------------------------\n");
        printf("   At theta_i = 90 deg (lateral wave):\n");
        printf("     Glass bottom:   R = 1.00 (total internal reflection)\n");
        printf("     Acrylic bottom: R = 1.00 (total internal reflection)\n");
        printf("     PDMS bottom:    R ~ 0.03 (energy leaks through!)\n");
        printf("     Open water:     R = 0.00 (no confinement)\n");
        printf("   => Glass slide CONFINES the lateral standing wave\n");
        printf("   => PDMS bottom would LEAK energy downward\n");
        printf("=====================================================\n\n");
    }
    output_glass_confinement();

    printf("\nDone. 16 CSV files total in %s/\n", OUTPUT_DIR);

#ifdef _WIN32
    printf("Press Enter to exit...");
    getchar();
#endif

    return 0;
}
