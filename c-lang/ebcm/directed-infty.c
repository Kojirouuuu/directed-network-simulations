#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static void build_poisson_pmf(double avg, int kmax, double *P) {
    P[0] = exp(-avg);
    for (int k = 1; k <= kmax; k++) {
        P[k] = P[k - 1] * (avg / (double)k);
    }

    double s = 0.0;
    for (int k = 0; k <= kmax; k++)
        s += P[k];
    for (int k = 0; k <= kmax; k++)
        P[k] /= s;
}

static void build_powerlaw_pmf(double gamma, int kmin, int kmax, double *P) {
    for (int k = 0; k <= kmax; k++) {
        P[k] = 0.0;
    }

    double Z = 0.0;
    for (int k = kmin; k <= kmax; k++) {
        Z += pow((double)k, -gamma);
    }
    for (int k = kmin; k <= kmax; k++) {
        P[k] = pow((double)k, -gamma) / Z;
    }
}

static void build_pmf(double avg, double gamma, int kmin, int kmax, const char *type, double *P) {
    if (strcmp(type, "Poi") == 0) {
        build_poisson_pmf(avg, kmax, P);
    } else if (strcmp(type, "Pow") == 0) {
        if (kmin < 0) {
            fprintf(stderr, "Error: kmin must be greater than 0 for power law "
                            "distribution\n");
            exit(1);
        }
        build_powerlaw_pmf(gamma, kmin, kmax, P);
    } else {
        fprintf(stderr, "Error: Invalid type: %s\n", type);
        exit(1);
    }
}

double calculate_mean_powerlaw(int min, int max, double gamma) {
    double Z = 0.0;
    double num = 0.0;
    for (int k = min; k <= max; k++) {
        double w = pow((double)k, -gamma);
        Z += w;
        num += (double)k * w;
    }
    return num / Z;
}

// 各次数分布（入/出/無向）の設定
typedef struct {
    double mean;
    int min;
    int max;
    double gamma;
    const char *type;
} DegreeConfig;

// EBCM 全体の次数分布設定
typedef struct {
    DegreeConfig ki;
} EBCMConfig;

// ダイナミクスパラメータ（lambda, mu, rho0, T）
typedef struct {
    int T;
    double rho0;
    double lambda_d;
    double mu;
} DynamicsConfig;

typedef struct {
    int ki_min, ki_max;
    double *Pi;
    double mean_ki;
    const char *type_i;
    double gamma_i;
} DegreeDist;

/* 分布タイプに応じた実効的な最大次数（ループ範囲に使う） */
static int get_effective_degree_max(const char *type, double mean, int max) {
    if (strcmp(type, "Poi") == 0) {
        return (int)(mean + 3.0 * sqrt(mean));
    }
    return max;
}

// ポアソン分布か冪分布
double *create_degree_dist(double mean, int min, int max, double gamma, const char *type) {
    double *P = calloc((size_t)max + 1, sizeof(double));
    if (!P) {
        fprintf(stderr, "Error: calloc P failed\n");
        exit(1);
    }
    build_pmf(mean, gamma, min, max, type, P);
    return P;
}

DegreeDist *build_all_degree_dist(const EBCMConfig *cfg) {
    const DegreeConfig *ki = &cfg->ki;

    DegreeDist *D = calloc(1, sizeof *D);
    if (!D) {
        fprintf(stderr, "Error: calloc DegreeDist failed\n");
        exit(1);
    }

    D->type_i = ki->type;
    D->gamma_i = ki->gamma;

    D->ki_min = (strcmp(D->type_i, "Pow") == 0) ? ki->min : 0;
    D->ki_max = get_effective_degree_max(D->type_i, ki->mean, ki->max);

    D->mean_ki = (strcmp(D->type_i, "Pow") == 0)
                     ? calculate_mean_powerlaw(D->ki_min, D->ki_max, D->gamma_i)
                     : ki->mean;

    // 配列生成（Poisson のとき min は無視される設計）
    D->Pi = create_degree_dist(D->mean_ki, D->ki_min, D->ki_max, D->gamma_i, D->type_i);

    return D;
}

static double *build_binom(int kmax) {
    double *Binom = calloc((size_t)(kmax + 1) * (kmax + 1), sizeof(double));
    if (!Binom) {
        fprintf(stderr, "Error: calloc Binom failed\n");
        exit(1);
    }
    for (int k = 0; k <= kmax; k++) {
        Binom[k * (kmax + 1) + 0] = 1.0;
        for (int m = 1; m <= kmax; m++) {
            Binom[k * (kmax + 1) + m] =
                Binom[k * (kmax + 1) + m - 1] * (double)(k - m + 1) / (double)m;
        }
    }
    return Binom;
}

static double Theta_d(const double *Binom, int ki, int kmax, int T, double theta_d) {
    if (theta_d == 1.0) {
        return 1.0;
    }
    if (T == 1) {
        return pow(theta_d, (double)ki);
    }
    if (ki == 0) {
        return 1.0;
    }
    int md_max = (ki < T - 1) ? ki : T - 1;
    double sum = 0.0;
    double q = 1.0 - theta_d;

    for (int md = 0; md <= md_max; md++) {
        sum += Binom[ki * (kmax + 1) + md] * pow(theta_d, (double)(ki - md)) * pow(q, (double)md);
    }
    return sum;
}

static double Theta_d_prime(const double *Binom, int ki, int kmax, int T, double theta_d) {
    if (ki == 0)
        return 0.0;

    if (T == 1) {
        // d/dθ θ^ki = ki θ^(ki-1)
        return (ki == 0) ? 0.0 : (double)ki * pow(theta_d, (double)(ki - 1));
    }

    int md_max = (ki < T - 1) ? ki : T - 1;
    double q = 1.0 - theta_d;
    double sum = 0.0;

    for (int md = 0; md <= md_max; md++) {
        int a = ki - md; // power of theta
        int b = md;      // power of q

        double term = 0.0;

        // a * theta^(a-1) * q^b
        if (a >= 1) {
            term += (double)a * pow(theta_d, (double)(a - 1)) * pow(q, (double)b);
        }
        // - b * theta^a * q^(b-1)
        if (b >= 1) {
            term -= (double)b * pow(theta_d, (double)a) * pow(q, (double)(b - 1));
        }

        sum += Binom[ki * (kmax + 1) + md] * term;
    }
    return sum;
}

static double Theta_d_prime_prime(const double *Binom, int ki, int kmax, int T, double theta_d) {
    if (T == 1) {
        fprintf(stderr, "Theta_d_prime_prime not needed for T=1\n");
        exit(1);
    }
    if (ki == 0)
        return 0.0;

    int md_max = (ki < T - 1) ? ki : T - 1;
    double q = 1.0 - theta_d;
    double sum = 0.0;

    for (int md = 0; md <= md_max; md++) {
        int a = ki - md; // power of theta
        int b = md;      // power of q
        double term = 0.0;

        if (a >= 2) {
            term += (double)a * (a - 1) * pow(theta_d, (double)(a - 2)) * pow(q, (double)b);
        }
        if (a >= 1 && b >= 1) {
            term += -2.0 * (double)a * (double)b * pow(theta_d, (double)(a - 1)) *
                    pow(q, (double)(b - 1));
        }
        if (b >= 2) {
            term += (double)b * (b - 1) * pow(theta_d, (double)a) * pow(q, (double)(b - 2));
        }

        sum += Binom[ki * (kmax + 1) + md] * term;
    }

    return sum;
}

static double xiS_directed(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                           double theta_d) {
    double s = 0.0;

    for (int k = D->ki_min; k <= D->ki_max; k++) {
        double pk = D->Pi[k];
        if (pk == 0.0)
            continue;
        double xiS_d = Theta_d(Binom, k, D->ki_max, p->T, theta_d);
        s += pk * xiS_d;
    }
    return (1.0 - p->rho0) * s;
}

static double xiS_directed_prime(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                                 double theta_d) {
    double s = 0.0;

    for (int k = D->ki_min; k <= D->ki_max; k++) {
        double pk = D->Pi[k];
        if (pk == 0.0)
            continue;
        double xiS_d_prime = Theta_d_prime(Binom, k, D->ki_max, p->T, theta_d);
        s += pk * xiS_d_prime;
    }
    return (1.0 - p->rho0) * s;
}

static double xiS_directed_prime_prime(const DegreeDist *D, const DynamicsConfig *p,
                                       const double *Binom, double theta_d) {
    double s = 0.0;

    for (int k = D->ki_min; k <= D->ki_max; k++) {
        double pk = D->Pi[k];
        if (pk == 0.0)
            continue;
        double xiS_d_prime_prime = Theta_d_prime_prime(Binom, k, D->ki_max, p->T, theta_d);
        s += pk * xiS_d_prime_prime;
    }
    return (1.0 - p->rho0) * s;
}

/* 感受性ノードの割合 S(t) = (1-ρ0) Σ P_i（因数分解時は xiS_directed と同じ） */
static double compute_S(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                        double theta_d) {
    return xiS_directed(D, p, Binom, theta_d);
}

static double rhs_Phi(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                      double theta_d) {
    double s = 0.0;

    for (int k = D->ki_min; k <= D->ki_max; k++) {
        double pk = D->Pi[k];
        if (pk == 0.0)
            continue;
        if (theta_d == 1.0) {
            if (p->T == 1) {
                return 1.0;
            }
            return 0.0;
        }
        if (p->T == 1) {
            s += pk * pow(theta_d, (double)k);
            continue;
        }
        if (k == 0) {
            continue;
        }
        double q = 1.0 - theta_d;
        double Phi = Binom[k * (D->ki_max + 1) + p->T - 1] *
                     pow(theta_d, (double)(k) - (double)(p->T - 1)) * pow(q, (double)(p->T - 1));
        s += pk * Phi;
    }
    return (1.0 - p->rho0) * s;
}

static double g_d(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                  double theta_d) {
    double xiS_d = xiS_directed(D, p, Binom, theta_d);
    return -p->lambda_d * (theta_d - xiS_d) + p->mu * (1.0 - theta_d);
}

static double g_d_prime(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                        double theta_d) {
    double xiS_d_prime = xiS_directed_prime(D, p, Binom, theta_d);
    return -p->lambda_d * (1 - xiS_d_prime) - p->mu;
}

static double g_d_prime_prime(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                              double theta_d) {
    double xiS_d_prime_prime = xiS_directed_prime_prime(D, p, Binom, theta_d);
    return p->lambda_d * xiS_d_prime_prime;
}

typedef double (*Func)(const DegreeDist *, const DynamicsConfig *, const double *, double);

static int bisect_interval(Func f, const DegreeDist *D, const DynamicsConfig *p,
                           const double *Binom, double a, double b, double tol, int itmax,
                           double *root) {
    double fa = f(D, p, Binom, a);
    double fb = f(D, p, Binom, b);

    if (fabs(fa) < tol) {
        *root = a;
        return 0;
    }
    if (fabs(fb) < tol) {
        *root = b;
        return 0;
    }
    if (fa * fb > 0.0)
        return -1;

    for (int it = 0; it < itmax; it++) {
        double m = 0.5 * (a + b);
        double fm = f(D, p, Binom, m);

        if (fabs(fm) < tol || fabs(b - a) < tol) {
            *root = m;
            return 0;
        }
        if (fa * fm <= 0.0) {
            b = m;
            fb = fm;
        } else {
            a = m;
            fa = fm;
        }
    }
    return -1;
}

static void append_root(double *roots, int *n, int max_roots, double r, double merge_eps) {
    // 近接根の重複を潰す（簡易）
    for (int i = 0; i < *n; i++) {
        if (fabs(roots[i] - r) < merge_eps)
            return;
    }
    if (*n < max_roots)
        roots[(*n)++] = r;
}

static int find_roots(Func f, const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                      double tmin, double tmax, double step, double tol, int itmax, double *roots,
                      int max_roots) {
    int n = 0;

    double a = tmin;
    double fa = f(D, p, Binom, a);

    for (double b = a + step; b <= tmax + 1e-12; b += step) {
        double fb = f(D, p, Binom, b);

        // 符号反転 → 二分法
        if (fa * fb < 0.0) {
            double r;
            if (bisect_interval(f, D, p, Binom, a, b, tol, itmax, &r) == 0) {
                append_root(roots, &n, max_roots, r, 0);
            }
        } else {
            // サンプル点がほぼ0なら拾う（接線根の取りこぼし補助）
            if (fabs(fa) < tol)
                append_root(roots, &n, max_roots, a, 0);
            if (fabs(fb) < tol)
                append_root(roots, &n, max_roots, b, 0);
        }

        a = b;
        fa = fb;
    }
    if (fabs(f(D, p, Binom, tmax)) < tol) {
        append_root(roots, &n, max_roots, tmax, 0.0);
    }
    return n;
}

#define MAX_GD_ROOTS 32

int main(void) {
    EBCMConfig cfg = {
        // .ki = {.mean = 12.0, .min = 5, .max = 500000, .gamma = 2.5, .type = "Pow"},
        .ki = {.mean = 12.0, .min = 5, .max = 707, .gamma = 2.5, .type = "Pow"},
    };
    double mu = 1.0;

    const int T_list[] = {3};
    const int T_count = (int)(sizeof(T_list) / sizeof(T_list[0]));

    const double lambda_d_min = 0.0;
    const double lambda_d_max = 4.0;
    const double lambda_d_step = 0.004;

    const double rho0_min = 0.0;
    const double rho0_max = 0.4;
    const double rho0_step = 0.0004;

    const double theta_search_step = 0.005; /* g_d=0 の根探索の刻み */

    int total_tasks = 0;
    int rho0_tasks = 0;
    for (double rho0 = rho0_min; rho0 <= rho0_max + 1e-6; rho0 += rho0_step) {
        rho0_tasks++;
    }
    int lambda_d_tasks = 0;
    for (double lambda_d = lambda_d_min; lambda_d <= lambda_d_max + 1e-6;
         lambda_d += lambda_d_step) {
        lambda_d_tasks++;
    }
    total_tasks = T_count * rho0_tasks * lambda_d_tasks;

    DegreeDist *D = build_all_degree_dist(&cfg);
    double *Binom = build_binom(D->ki_max);

    const int progress_width = 100;
    char dirbuf[256];
    char pathbuf[256];
    mkdir("out", 0755);
    mkdir("out/ebcm", 0755);
    mkdir("out/ebcm/directed-infty", 0755);
    if (strcmp(cfg.ki.type, "Pow") == 0) {
        snprintf(dirbuf, sizeof dirbuf, "out/ebcm/directed-infty/Pow/gamma=%.2f/kmin=%d/kmax=%d",
                 cfg.ki.gamma, D->ki_min, D->ki_max);
        mkdir("out/ebcm/directed-infty/Pow", 0755);
        char subdir[256];
        snprintf(subdir, sizeof subdir, "out/ebcm/directed-infty/Pow/gamma=%.2f", cfg.ki.gamma);
        mkdir(subdir, 0755);
        snprintf(subdir, sizeof subdir, "out/ebcm/directed-infty/Pow/gamma=%.2f/kmin=%d",
                 cfg.ki.gamma, D->ki_min);
        mkdir(subdir, 0755);
        mkdir(dirbuf, 0755);
    } else if (strcmp(cfg.ki.type, "Poi") == 0) {
        snprintf(dirbuf, sizeof dirbuf, "out/ebcm/directed-infty/Poi/kuave=%.2f", D->mean_ki);
        mkdir("out/ebcm/directed-infty/Poi", 0755);
        mkdir(dirbuf, 0755);
    } else {
        snprintf(dirbuf, sizeof dirbuf, "out/ebcm/directed-infty/%s", cfg.ki.type);
        mkdir(dirbuf, 0755);
    }

    snprintf(pathbuf, sizeof pathbuf, "%s/gd_zero.csv", dirbuf);
    FILE *fp_gd = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gdp_zero.csv", dirbuf);
    FILE *fp_gdp = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gdpp_zero.csv", dirbuf);
    FILE *fp_gdpp = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gd_gdp_zero.csv", dirbuf);
    FILE *fp_gd_gdp_zero = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gd_gdp_gdpp_zero.csv", dirbuf);
    FILE *fp_gd_gdp_gdpp_zero = fopen(pathbuf, "w");
    if (fp_gd == NULL || fp_gdp == NULL || fp_gdpp == NULL || fp_gd_gdp_zero == NULL ||
        fp_gd_gdp_gdpp_zero == NULL) {
        perror("fopen");
        return 1;
    }
    fprintf(fp_gd, "T,rho0,lambda_d,mu,theta_d,R,Phi\n");
    fprintf(fp_gdp, "T,rho0,lambda_d,mu,theta_d\n");
    fprintf(fp_gdpp, "T,rho0,lambda_d,mu,theta_d\n");
    fprintf(fp_gd_gdp_zero, "T,rho0,lambda_d,mu,theta_d\n");
    fprintf(fp_gd_gdp_gdpp_zero, "T,rho0,lambda_d,mu\n");
    /* 実行開始時刻・進捗 */
    time_t t_start = time(NULL);
    struct tm *tm_start = localtime(&t_start);
    char buf_start[64], buf_end[64];
    strftime(buf_start, sizeof buf_start, "%Y-%m-%d %H:%M:%S", tm_start);
#ifdef _OPENMP
    fprintf(stderr, "start: %s | threads: %d | tasks: %d\n", buf_start, omp_get_max_threads(),
            total_tasks);
#else
    fprintf(stderr, "start: %s | tasks: %d\n", buf_start, total_tasks);
#endif
    int current_task = 0;
    fprintf(stderr, "\r[");
    for (int i = 0; i < progress_width; i++)
        fputc('-', stderr);
    fprintf(stderr, "] 0/%d (0.0%%)", total_tasks);
    fflush(stderr);

#pragma omp parallel for schedule(dynamic)
    for (int task = 0; task < total_tasks; task++) {
        /* T, rho, lambda の順: task = ii*(rho*lambda) + jj*lambda + kk */
        int ii = task / (rho0_tasks * lambda_d_tasks); /* T index */
        int jj = (task / lambda_d_tasks) % rho0_tasks; /* rho index */
        int kk = task % lambda_d_tasks;                /* lambda index */
        int T = T_list[ii];
        double rho0 = rho0_min + jj * rho0_step;
        double lambda_d = lambda_d_min + kk * lambda_d_step;
        DynamicsConfig dynamics = {
            .T = T,
            .rho0 = rho0,
            .lambda_d = lambda_d,
            .mu = mu,
        };

        /* g_d(theta)=0 を満たす theta in [0,1] の根を列挙し、最大を採用 */
        double roots[MAX_GD_ROOTS];
        double roots_prime[MAX_GD_ROOTS];
        double roots_prime_prime[MAX_GD_ROOTS];

        int n_roots = 0;
        int n_roots_prime = 0;
        int n_roots_prime_prime = 0;

        double valid_theta_d = 0.0;

        double tmin, tmax, step;
        if (T == 1) {
            tmin = 0.0;
            tmax = 1.0 - 1e-10;
            step = 0.5;
        } else {
            tmin = 0.0;
            tmax = 1.0;
            step = 0.005;
        }

        if (T == 1 && rho0 == 0.0) {
            double m = g_d_prime(D, &dynamics, Binom, 1.0);
            if (fabs(m) < 1e-10) {
                n_roots = 1;
                roots[0] = 0.0;
                n_roots_prime = 1;
                roots_prime[0] = 1.0;
            } else if (m < 0.0) {
                n_roots = 1;
                roots[0] = 1.0;
            } else {
                n_roots = find_roots(g_d, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000, roots,
                                     MAX_GD_ROOTS);
                while (n_roots == 0 && tmin > 0.0) {
                    tmin *= 0.9;
                    n_roots = find_roots(g_d, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000,
                                         roots, MAX_GD_ROOTS);
                }
                for (int r = 0; r < n_roots; r++) {
                    if (roots[r] < 1.0 && roots[r] >= valid_theta_d) {
                        valid_theta_d = roots[r];
                    }
                }
            }
        } else {
            n_roots = find_roots(g_d, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000, roots,
                                 MAX_GD_ROOTS);
            while (n_roots == 0 && tmin > 0.0) {
                tmin *= 0.9;
                n_roots = find_roots(g_d, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000, roots,
                                     MAX_GD_ROOTS);
            }
            for (int r = 0; r < n_roots; r++) {
                if (roots[r] <= 1.0 && roots[r] >= valid_theta_d) {
                    valid_theta_d = roots[r];
                }
            }
        }

        if (valid_theta_d == 0.0 && T == 1) {
            valid_theta_d = 1.0;
        }
        /* 平衡時は A=0 なので R_inf = 1 - S_inf */
        double S_inf = compute_S(D, &dynamics, Binom, valid_theta_d);
        double R_inf = 1.0 - S_inf;
        double Phi_inf = rhs_Phi(D, &dynamics, Binom, valid_theta_d);

        n_roots_prime = find_roots(g_d_prime, D, &dynamics, Binom, 0.0, 1.0, theta_search_step,
                                   1e-6, 100, roots_prime, MAX_GD_ROOTS);

        if (T > 1) {
            n_roots_prime_prime =
                find_roots(g_d_prime_prime, D, &dynamics, Binom, 0.0, 1.0, theta_search_step, 1e-6,
                           100, roots_prime_prime, MAX_GD_ROOTS);
        }

        bool prime_found = false;
        double delta_min = 0.02; // InPow200*200: 0.030
        for (int r = 0; r < n_roots_prime; r++) {
            double root = roots_prime[r];
            double delta = fabs(root - valid_theta_d);
            double g_d_check = g_d(D, &dynamics, Binom, root);
            if (root <= 1.0 && root >= 0.0 && delta < delta_min && fabs(g_d_check) < 1e-4) {
                prime_found = true;
                delta_min = delta;
            }
        }

        bool prime_prime_found = false;
        double delta_min_prime_prime = 0.08; // OutPow200*200: 0.04
        if (prime_found && T > 1) {
            for (int r = 0; r < n_roots_prime_prime; r++) {
                double root = roots_prime_prime[r];
                double delta = fabs(root - valid_theta_d);
                double g_d_check = g_d(D, &dynamics, Binom, root);
                double g_d_prime_check = g_d_prime(D, &dynamics, Binom, root);
                if (root <= 1.0 && root >= 0.0 && delta < delta_min_prime_prime &&
                    fabs(g_d_prime_check) < 1e-2 && fabs(g_d_check) < 1e-2) {
                    prime_prime_found = true;
                    delta_min_prime_prime = delta;
                }
            }
        }

#pragma omp critical
        {
            fprintf(fp_gd, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_d, mu,
                    valid_theta_d, R_inf, Phi_inf);
            for (int r = 0; r < n_roots_prime; r++) {
                double root = roots_prime[r];
                if (root <= 1.0 && root >= 0.0) {
                    fprintf(fp_gdp, "%d,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_d, mu, root);
                }
            }

            for (int r = 0; r < n_roots_prime_prime; r++) {
                if (roots_prime_prime[r] <= 1.0 && roots_prime_prime[r] >= 0.0) {
                    fprintf(fp_gdpp, "%d,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_d, mu,
                            roots_prime_prime[r]);
                }
            }

            if (prime_found) {
                fprintf(fp_gd_gdp_zero, "%d,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_d, mu,
                        valid_theta_d);
            }
            if (prime_prime_found) {
                fprintf(fp_gd_gdp_gdpp_zero, "%d,%.6f,%.6f,%.6f\n", T, rho0, lambda_d, mu);
            }
            current_task++;
            if (current_task % 10 == 0 || current_task == total_tasks) {
                fprintf(stderr, "\r[");
                int filled = progress_width * current_task / total_tasks;
                for (int p = 0; p < progress_width; p++)
                    fputc(p < filled ? '#' : '-', stderr);
                fprintf(stderr, "] %d/%d (%.1f%%)   ", current_task, total_tasks,
                        100.0 * current_task / total_tasks);
            }
        }
    }

    time_t t_end = time(NULL);
    struct tm *tm_end = localtime(&t_end);
    strftime(buf_end, sizeof buf_end, "%Y-%m-%d %H:%M:%S", tm_end);
    double elapsed = difftime(t_end, t_start);
    fprintf(stderr, "\nend: %s | total: %.1f seconds (%.1f minutes)\n", buf_end, elapsed,
            elapsed / 60.0);

    fclose(fp_gd);
    fclose(fp_gdp);
    fclose(fp_gdpp);
    fclose(fp_gd_gdp_zero);
    fclose(fp_gd_gdp_gdpp_zero);
    free(Binom);
    free(D);
    return 0;
}