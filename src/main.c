/*
 * echo_chamber_model.c
 *
 * Standing-Wave "Echo Chamber" Model + Coupling Layer Transmission
 * for acoustic tweezer channel design.
 *
 * Physics from the analytical framework:
 *   A) Standing wave in channel: p(x) = p0*exp(...) + r2*p0*exp(...)
 *   B) Coupling layer: Zin(d) impedance transformer model
 *
 * Uses same material parameters as the transfer-matrix model.
 *
 * Compile:  gcc -std=c99 -O2 -o echo_chamber echo_chamber_model.c -lm
 * Run:      ./echo_chamber
 * Output:   CSV files for plotting
 *
 * License:  do not care
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <errno.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Alias for imaginary unit to avoid conflict with variable names */
#define IMAG_UNIT (_Complex_I)

/* Output directory — change YOUR_USERNAME to your actual username */
#define OUTPUT_DIR "/home/shh/Projects/BENG207_2/PLOT_CSV/"

/* ============================================================
 * Material parameters (same stack as transfer-matrix model) Anna: Talk to Alexi about these.
 * ============================================================ */

/* Transducer: LiNbO3 */
static const double rho_1    = 4647.0;     /* kg/m^3 */
static const double c_1_long = 6570.0;     /* m/s */

/* Couplant: US gel */
static const double rho_2   = 1020.0;      /* kg/m^3 */
static const double c_2     = 1500.0;      /* m/s (bulk sound speed in gel) */

/* Glass superstrate: borosilicate */
static const double rho_3 = 2230.0;        /* kg/m^3 */
static const double c_3   = 5640.0;        /* m/s */

/* Fluid: water */
static const double rho_f = 1000.0;        /* kg/m^3 */
static const double c_f   = 1483.0;        /* m/s */

/* Frequency */
static const double f0 = 20.0e6;           /* 20 MHz */

/* Derived (set in main) */
static double Z_1, Z_2, Z_3, Z_f;
static double omega0, k_f, lambda_f, lambda_2;

/* ============================================================
 * Utility
 * ============================================================ */

static void linspace(double a, double b, int n, double *out)
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = a + (b - a) * i / (n - 1);
}

/* (logspace removed — unused in this model) */

/* ============================================================
 * A) STANDING WAVE MODEL
 *
 * p(x) = p0 * exp(-(alpha - ik)*x)
 *       + r2 * p0 * exp(-(alpha - ik)*(2L - x))
 *
 * S(x) = |r2| * exp(-2*alpha*(L - x))
 * contrast = ((1+S)/(1-S))^2
 * ============================================================ */

/* Compute |p(x)|^2 normalized to p0^2 */
static double pressure_field(double x, double L, double alpha,
                              double k, double complex r2)
{
    double complex gamma = alpha - IMAG_UNIT * k;
    double complex p_fwd = cexp(-gamma * x);
    double complex p_ref = r2 * cexp(-gamma * (2.0*L - x));
    double complex p_tot = p_fwd + p_ref;
    return cabs(p_tot);
}

/* Standing wave ratio S(x) */
static double sw_ratio(double x, double L, double alpha, double r_mag)
{
    return r_mag * exp(-2.0 * alpha * (L - x));
}

/* Standing wave contrast at position x */
static double sw_contrast(double S)
{
    if (S >= 1.0) return 1e6;  /* clamp */
    return ((1.0 + S) / (1.0 - S)) * ((1.0 + S) / (1.0 - S));
}

/* ============================================================
 * OUTPUT 1: Pressure field p(x) for various L and |r|
 * ============================================================ */

static void output_pressure_field(void)
{
    int iL, ir, ix;
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "echo_fig1_pressure_field.csv");
    FILE *fp = fopen(fname, "w");

    /* Attenuation in water at 20 MHz: ~0.1 Np/m (conservative) */
    double alpha = 0.1;  /* Np/m — adjustable */

    /* Channel lengths */
    double L_vals[] = {2.5e-3, 5.0e-3, 10.0e-3};  /* 2.5, 5, 10 mm */
    int NL = 3;

    /* Reflection coefficients */
    double r_vals[] = {0.95, 0.7, 0.3};  /* rigid, moderate, leaky */
    int Nr = 3;

    int Nx = 500;

    fprintf(fp, "x_mm,L_mm,r_mag,p_norm\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        double x_arr[500];
        linspace(0, L, Nx, x_arr);

        for (ir = 0; ir < Nr; ir++) {
            double complex r2 = r_vals[ir];  /* real, positive (rigid) */

            for (ix = 0; ix < Nx; ix++) {
                double p = pressure_field(x_arr[ix], L, alpha,
                                           k_f, r2);
                fprintf(fp, "%.6f,%.2f,%.2f,%.8f\n",
                        x_arr[ix]*1e3, L*1e3, r_vals[ir], p);
            }
        }
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 2: Standing wave contrast vs L and |r|
 * ============================================================ */

static void output_contrast_map(void)
{
    int iL, ir;
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "echo_fig2_constant.csv");
    FILE *fp = fopen(fname, "w");


    double alpha = 0.1;  /* Np/m */

    int NL = 100, Nr = 100;
    double L_arr[100], r_arr[100];
    linspace(0.5e-3, 15e-3, NL, L_arr);
    linspace(0.05, 0.99, Nr, r_arr);

    fprintf(fp, "L_mm,r_mag,S_mid,contrast_mid,contrast_dB\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_arr[iL];
        double x_mid = L / 2.0;  /* evaluate at channel midpoint */

        for (ir = 0; ir < Nr; ir++) {
            double S = sw_ratio(x_mid, L, alpha, r_arr[ir]);
            double C = sw_contrast(S);
            double C_dB = 10.0 * log10(C);

            fprintf(fp, "%.4f,%.4f,%.6f,%.4f,%.2f\n",
                    L*1e3, r_arr[ir], S, C, C_dB);
        }
    }

    fclose(fp);
    printf("  Written: %s (%dx%d)\n", fname, NL, Nr);
}

/* ============================================================
 * OUTPUT 3: Resonance spacing vs channel length
 * ============================================================ */

static void output_resonance_spacing(void)
{
    int i;
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "echo_fig3_resonance_spacing.csv");
    FILE *fp = fopen(fname, "w");


    /* Assume transducer bandwidth ~ 2 MHz around 20 MHz */
    double BW = 2.0e6;
    int NL = 200;
    double L_arr[200];
    linspace(0.5e-3, 20e-3, NL, L_arr);

    fprintf(fp, "L_mm,delta_f_kHz,n_modes_in_BW\n");

    for (i = 0; i < NL; i++) {
        double df = c_f / (2.0 * L_arr[i]);
        double n_modes = BW / df;
        fprintf(fp, "%.4f,%.4f,%.2f\n",
                L_arr[i]*1e3, df/1e3, n_modes);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 4: Reflected/forward ratio S(x) along channel
 *           for various alpha and L
 * ============================================================ */

static void output_sw_ratio_profile(void)
{
    int iL, ia, ix;
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR,"echo_fig4_sw_ratio_profile.csv");
    FILE *fp = fopen(fname, "w");

    double r_mag = 0.8;  /* decent reflection */
    double alpha_vals[] = {0.01, 0.1, 0.5, 2.0};  /* Np/m */
    int Na = 4;

    double L_vals[] = {2.5e-3, 5.0e-3, 10.0e-3};
    int NL = 3;
    int Nx = 300;

    fprintf(fp, "x_mm,L_mm,alpha_Npm,S,contrast\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_vals[iL];
        double x_arr[300];
        linspace(0, L, Nx, x_arr);

        for (ia = 0; ia < Na; ia++) {
            for (ix = 0; ix < Nx; ix++) {
                double S = sw_ratio(x_arr[ix], L,
                                     alpha_vals[ia], r_mag);
                double C = sw_contrast(S);
                fprintf(fp, "%.6f,%.2f,%.3f,%.6f,%.4f\n",
                        x_arr[ix]*1e3, L*1e3,
                        alpha_vals[ia], S, C);
            }
        }
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * B) COUPLING LAYER MODEL
 *
 * Zin(d) = Z2 * (Z3 + i*Z2*tan(k2*d)) / (Z2 + i*Z3*tan(k2*d))
 *
 * Transmission coefficient through 3-layer stack:
 * T = 2*Z3 / (Zin + Z1) * ... (simplified for power)
 *
 * We compute:
 *   - |Zin(d)| vs d/lambda
 *   - Power transmission coefficient vs d/lambda
 * ============================================================ */

/* Input impedance looking into couplant backed by glass */
static double complex Z_input(double d, double freq __attribute__((unused)),
                               double Z2, double Z3_load,
                               double k2)
{
    double complex tan_kd = ctan(k2 * d);  /* real k2 for lossless */
    double complex num = Z3_load + IMAG_UNIT * Z2 * tan_kd;
    double complex den = Z2 + IMAG_UNIT * Z3_load * tan_kd;
    return Z2 * num / den;
}

/* Power transmission through transducer → couplant → glass */
static double power_transmission(double d, double freq)
{
    double k2 = 2.0 * M_PI * freq / c_2;
    double complex Zin = Z_input(d, freq, Z_2, Z_3, k2);

    /* Reflection at transducer-couplant interface */
    double complex r = (Zin - Z_1) / (Zin + Z_1);
    double R = cabs(r) * cabs(r);
    return 1.0 - R;  /* transmitted power fraction */
}

static void output_coupling_layer(void)
{
    int i;
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "echo_fig5_coupling_layer.csv");
    FILE *fp = fopen(fname, "w");

    int Nd = 1000;
    double d_norm[1000];  /* d/lambda_2 */
    double k2 = 2.0 * M_PI * f0 / c_2;
    linspace(0.001, 2.0, Nd, d_norm);

    fprintf(fp, "d_over_lambda,d_um,Zin_real_MRayl,Zin_imag_MRayl,"
                "Zin_mag_MRayl,T_power,T_dB\n");

    for (i = 0; i < Nd; i++) {
        double d = d_norm[i] * lambda_2;
        double complex Zin = Z_input(d, f0, Z_2, Z_3, k2);
        double T = power_transmission(d, f0);
        double T_dB = 10.0 * log10(fmax(T, 1e-15));

        fprintf(fp, "%.6f,%.4f,%.6f,%.6f,%.6f,%.8f,%.4f\n",
                d_norm[i], d*1e6,
                creal(Zin)/1e6, cimag(Zin)/1e6, cabs(Zin)/1e6,
                T, T_dB);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 6: Coupling sensitivity — T vs d for multiple freqs
 * ============================================================ */

static void output_coupling_vs_freq(void)
{
    int i, j;
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "echo_fig6_coupling_vs_freq.csv");
    FILE *fp = fopen(fname, "w");

    double freqs[] = {5e6, 10e6, 20e6, 50e6};
    int Nf = 4;
    int Nd = 500;

    /* Thickness in µm (absolute, not normalized) */
    double d_um[500];
    linspace(0.1, 100.0, Nd, d_um);

    fprintf(fp, "d_um");
    for (j = 0; j < Nf; j++)
        fprintf(fp, ",T_%.0fMHz,d_over_lam_%.0fMHz",
                freqs[j]/1e6, freqs[j]/1e6);
    fprintf(fp, "\n");

    for (i = 0; i < Nd; i++) {
        double d = d_um[i] * 1e-6;  /* convert to meters */
        fprintf(fp, "%.4f", d_um[i]);

        for (j = 0; j < Nf; j++) {
            double lam = c_2 / freqs[j];
            double d_norm = d / lam;
            double T = power_transmission(d, freqs[j]);
            fprintf(fp, ",%.8f,%.6f", T, d_norm);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 7: Combined figure — design space
 *           Contrast vs L for rigid/leaky ends +
 *           coupling T vs d at 20 MHz
 * ============================================================ */

static void output_design_space(void)
{
    int j;
    char fname[256];
    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "echo_fig7_design_space.csv");
    FILE *fp = fopen(fname, "w");

    /* Part A: contrast vs L for different end conditions */
    double alpha = 0.1;
    int NL = 200;
    double L_arr[200];

    /* r values representing: rigid Si/glass, moderate seal, open/leaky */
    double r_design[] = {0.95, 0.80, 0.50, 0.20};
    const char *r_labels[] = {"rigid", "good_seal", "moderate", "leaky"};
    int Nr = 4;
    int iL, ir;

    linspace(0.5e-3, 15e-3, NL, L_arr);

    fprintf(fp, "L_mm");
    for (j = 0; j < Nr; j++)
        fprintf(fp, ",contrast_dB_%s_r%.2f", r_labels[j], r_design[j]);
    fprintf(fp, ",delta_f_kHz\n");

    for (iL = 0; iL < NL; iL++) {
        double L = L_arr[iL];
        double x_mid = L / 2.0;
        double df;
        fprintf(fp, "%.4f", L*1e3);

        for (ir = 0; ir < Nr; ir++) {
            double S = sw_ratio(x_mid, L, alpha, r_design[ir]);
            double C = sw_contrast(S);
            double C_dB = 10.0 * log10(C);
            fprintf(fp, ",%.4f", C_dB);
        }

        df = c_f / (2.0 * L);
        fprintf(fp, ",%.4f\n", df/1e3);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * MAIN
 * ============================================================ */
int main(void)
{

    /* Derived quantities */
    Z_1 = rho_1 * c_1_long;
    Z_2 = rho_2 * c_2;
    Z_3 = rho_3 * c_3;
    Z_f = rho_f * c_f;

    omega0   = 2.0 * M_PI * f0;
    k_f      = omega0 / c_f;
    lambda_f = c_f / f0;
    lambda_2 = c_2 / f0;

    printf("=====================================================\n");
    printf(" Echo Chamber + Coupling Layer Model\n");
    printf("=====================================================\n");
    printf(" f0           = %.1f MHz\n", f0/1e6);
    printf(" Z_transducer = %.2f MRayl (LiNbO3)\n", Z_1/1e6);
    printf(" Z_couplant   = %.2f MRayl (US gel)\n", Z_2/1e6);
    printf(" Z_glass      = %.2f MRayl (borosilicate)\n", Z_3/1e6);
    printf(" Z_fluid      = %.2f MRayl (water)\n", Z_f/1e6);
    printf(" lambda_fluid = %.1f um\n", lambda_f*1e6);
    printf(" lambda_coupl = %.1f um\n", lambda_2*1e6);
    printf(" lambda/4 gel = %.1f um\n", lambda_2*1e6/4);
    printf(" lambda/10 gel= %.1f um (thin-layer limit)\n", lambda_2*1e6/10);
    printf("=====================================================\n\n");

    printf("Generating output files...\n\n");

    output_pressure_field();
    output_contrast_map();
    output_resonance_spacing();
    output_sw_ratio_profile();
    output_coupling_layer();
    output_coupling_vs_freq();
    output_design_space();

    printf("\nDone. %d CSV files generated.\n", 7);

    return 0;
}
