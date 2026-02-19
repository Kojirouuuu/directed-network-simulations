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
    DegreeConfig ku;
} EBCMConfig;

// ダイナミクスパラメータ（lambda, mu, rho0, T）
typedef struct {
    int T;
    double rho0;
    double lambda_u;
    double mu;
} DynamicsConfig;

typedef struct {
    int ku_min, ku_max;
    double *Pu;
    double mean_ku;
    const char *type_u;
    double gamma_u;
} DegreeDist;

/* 分布タイプに応じた実効的な最大次数（ループ範囲に使う） */
static int get_effective_degree_max(const char *type, double mean, int max) {
    if (strcmp(type, "Pow") == 0) {
        return (int)(pow((double)max, 0.5));
    }
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
    const DegreeConfig *ku = &cfg->ku;

    DegreeDist *D = calloc(1, sizeof *D);
    if (!D) {
        fprintf(stderr, "Error: calloc DegreeDist failed\n");
        exit(1);
    }

    D->type_u = ku->type;
    D->gamma_u = ku->gamma;

    D->ku_min = (strcmp(D->type_u, "Pow") == 0) ? ku->min : 0;
    D->ku_max = get_effective_degree_max(D->type_u, ku->mean, ku->max);

    D->mean_ku = (strcmp(D->type_u, "Pow") == 0)
                     ? calculate_mean_powerlaw(D->ku_min, D->ku_max, D->gamma_u)
                     : ku->mean;

    // 配列生成（Poisson のとき min は無視される設計）
    D->Pu = create_degree_dist(D->mean_ku, D->ku_min, D->ku_max, D->gamma_u, D->type_u);

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

static double Theta_u(const double *Binom, int ku_minus_1, int kmax, int T, double theta_u) {
    if (theta_u == 1.0) {
        return 1.0;
    }
    if (T == 1) {
        return pow(theta_u, (double)ku_minus_1);
    }
    if (ku_minus_1 == 0) {
        return 1.0;
    }
    int md_max = (ku_minus_1 < T - 1) ? ku_minus_1 : T - 1;
    double sum = 0.0;
    double q = 1.0 - theta_u;

    for (int md = 0; md <= md_max; md++) {
        sum += Binom[ku_minus_1 * (kmax + 1) + md] * pow(theta_u, (double)(ku_minus_1 - md)) *
               pow(q, (double)md);
    }
    return sum;
}

static double Theta_u_prime(const double *Binom, int ku_minus_1, int kmax, int T, double theta_u) {
    if (ku_minus_1 == 0)
        return 0.0;

    if (T == 1) {
        // d/dθ θ^(k-1) = (k-1) θ^(k-2)
        return (double)(ku_minus_1)*pow(theta_u, (double)(ku_minus_1 - 1));
    }

    int md_max = (ku_minus_1 < T - 1) ? ku_minus_1 : T - 1;
    double q = 1.0 - theta_u;
    double sum = 0.0;

    for (int md = 0; md <= md_max; md++) {
        int a = ku_minus_1 - md; // power of theta
        int b = md;              // power of q

        double term = 0.0;

        // a * theta^(a-1) * q^b
        if (a >= 1) {
            term += (double)a * pow(theta_u, (double)(a - 1)) * pow(q, (double)b);
        }
        // - b * theta^a * q^(b-1)
        if (b >= 1) {
            term -= (double)b * pow(theta_u, (double)a) * pow(q, (double)(b - 1));
        }

        sum += Binom[ku_minus_1 * (kmax + 1) + md] * term;
    }
    return sum;
}

static double Theta_u_prime_prime(const double *Binom, int ku_minus_1, int kmax, int T,
                                  double theta_u) {
    if (T == 1) {
        fprintf(stderr, "Theta_u_prime_prime not needed for T=1\n");
        exit(1);
    }
    if (ku_minus_1 == 0)
        return 0.0;

    int md_max = (ku_minus_1 < T - 1) ? ku_minus_1 : T - 1;
    double q = 1.0 - theta_u;
    double sum = 0.0;

    for (int md = 0; md <= md_max; md++) {
        int a = ku_minus_1 - md; // power of theta
        int b = md;              // power of q
        double term = 0.0;

        if (a >= 2) {
            term += (double)a * (a - 1) * pow(theta_u, (double)(a - 2)) * pow(q, (double)b);
        }
        if (a >= 1 && b >= 1) {
            term += -2.0 * (double)a * (double)b * pow(theta_u, (double)(a - 1)) *
                    pow(q, (double)(b - 1));
        }
        if (b >= 2) {
            term += (double)b * (b - 1) * pow(theta_u, (double)a) * pow(q, (double)(b - 2));
        }

        sum += Binom[ku_minus_1 * (kmax + 1) + md] * term;
    }

    return sum;
}

static double xiS_undirected(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                             double theta_u) {
    double s = 0.0;

    for (int k = D->ku_min; k <= D->ku_max; k++) {
        if (k == 0) {
            continue;
        }
        double pk = D->Pu[k];
        if (pk == 0.0)
            continue;
        double xiS_u = Theta_u(Binom, k - 1, D->ku_max, p->T, theta_u);
        s += (double)k * pk * xiS_u;
    }
    return (1.0 - p->rho0) * s / (double)D->mean_ku;
}

static double xiS_undirected_prime(const DegreeDist *D, const DynamicsConfig *p,
                                   const double *Binom, double theta_u) {
    double s = 0.0;

    for (int k = D->ku_min; k <= D->ku_max; k++) {
        if (k == 0) {
            continue;
        }
        double pk = D->Pu[k];
        if (pk == 0.0)
            continue;
        double xiS_u_prime = Theta_u_prime(Binom, k - 1, D->ku_max, p->T, theta_u);
        s += (double)k * pk * xiS_u_prime;
    }
    return (1.0 - p->rho0) * s / (double)D->mean_ku;
}

static double xiS_undirected_prime_prime(const DegreeDist *D, const DynamicsConfig *p,
                                         const double *Binom, double theta_u) {
    double s = 0.0;

    for (int k = D->ku_min; k <= D->ku_max; k++) {
        if (k == 0) {
            continue;
        }
        double pk = D->Pu[k];
        if (pk == 0.0)
            continue;
        double xiS_u_prime_prime = Theta_u_prime_prime(Binom, k - 1, D->ku_max, p->T, theta_u);
        s += (double)k * pk * xiS_u_prime_prime;
    }
    return (1.0 - p->rho0) * s / (double)D->mean_ku;
}

/* 感受性ノードの割合 S(t) = (1-ρ0) Σ k*P(k)/⟨k⟩ * Theta_u */
static double compute_S(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                        double theta_u) {
    return xiS_undirected(D, p, Binom, theta_u);
}

static double rhs_Phi(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                      double theta_u) {
    int m = p->T - 1; // subcritical: m = T-1
    if (m < 0)
        return 0.0;

    double q = 1.0 - theta_u;
    double s = 0.0;
    int stride = D->ku_max + 1;

    for (int k = D->ku_min; k <= D->ku_max; k++) {
        double pk = D->Pu[k];
        if (pk == 0.0)
            continue;

        if (k < m)
            continue; // C(k,m)=0 の領域を明示的に除外

        double Ckm = Binom[k * stride + m]; // binom(k,m)
        s += pk * Ckm * pow(theta_u, (double)(k - m)) * pow(q, (double)m);
    }

    return (1.0 - p->rho0) * s;
}

static double g_u(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                  double theta_u) {
    double xiS_u = xiS_undirected(D, p, Binom, theta_u);
    return -p->lambda_u * (theta_u - xiS_u) + p->mu * (1.0 - theta_u);
}

static double g_u_prime(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                        double theta_u) {
    double xiS_u_prime = xiS_undirected_prime(D, p, Binom, theta_u);
    return -p->lambda_u * (1 - xiS_u_prime) - p->mu;
}

static double g_u_prime_prime(const DegreeDist *D, const DynamicsConfig *p, const double *Binom,
                              double theta_u) {
    double xiS_u_prime_prime = xiS_undirected_prime_prime(D, p, Binom, theta_u);
    return p->lambda_u * xiS_u_prime_prime;
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
    int N = 100000;
    EBCMConfig cfg = {
        .ku = {.mean = 6.0, .min = 5, .max = N, .gamma = 2.5, .type = "Pow"},
    };
    double mu = 1.0;

    const int T_list[] = {1, 2, 3, 4};
    const int T_count = (int)(sizeof(T_list) / sizeof(T_list[0]));

    const double lambda_u_min = 0.0;
    const double lambda_u_max = 0.2;
    const double lambda_u_step = 0.001;

    const double rho0_min = 0.0;
    const double rho0_max = 0.4;
    const double rho0_step = 0.002;

    const double theta_search_step = 0.005; /* g_u=0 の根探索の刻み */

    int total_tasks = 0;
    int rho0_tasks = 0;
    for (double rho0 = rho0_min; rho0 <= rho0_max + 1e-6; rho0 += rho0_step) {
        rho0_tasks++;
    }
    int lambda_u_tasks = 0;
    for (double lambda_u = lambda_u_min; lambda_u <= lambda_u_max + 1e-6;
         lambda_u += lambda_u_step) {
        lambda_u_tasks++;
    }
    total_tasks = T_count * rho0_tasks * lambda_u_tasks;

    DegreeDist *D = build_all_degree_dist(&cfg);
    double *Binom = build_binom(D->ku_max);

    const int progress_width = 100;
    char dirbuf[256];
    char pathbuf[256];
    snprintf(dirbuf, sizeof dirbuf, "out/ebcm/undirected-infty/%s", cfg.ku.type);
    mkdir("out/ebcm/undirected-infty", 0755); /* 親が無いと子が作れない */
    mkdir(dirbuf, 0755);

    snprintf(pathbuf, sizeof pathbuf, "%s/gu_zero.csv", dirbuf);
    FILE *fp_gu = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gup_zero.csv", dirbuf);
    FILE *fp_gup = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gupp_zero.csv", dirbuf);
    FILE *fp_gupp = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gu_gup_zero.csv", dirbuf);
    FILE *fp_gu_gup_zero = fopen(pathbuf, "w");
    snprintf(pathbuf, sizeof pathbuf, "%s/gu_gup_gupp_zero.csv", dirbuf);
    FILE *fp_gu_gup_gupp_zero = fopen(pathbuf, "w");
    if (fp_gu == NULL || fp_gup == NULL || fp_gupp == NULL || fp_gu_gup_zero == NULL ||
        fp_gu_gup_gupp_zero == NULL) {
        perror("fopen");
        return 1;
    }
    fprintf(fp_gu, "T,rho0,lambda_u,mu,theta_u,R,Phi\n");
    fprintf(fp_gup, "T,rho0,lambda_u,mu,theta_u\n");
    fprintf(fp_gupp, "T,rho0,lambda_u,mu,theta_u\n");
    fprintf(fp_gu_gup_zero, "T,rho0,lambda_u,mu,theta_u\n");
    fprintf(fp_gu_gup_gupp_zero, "T,rho0,lambda_u,mu\n");
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
        int ii = task / (rho0_tasks * lambda_u_tasks); /* T index */
        int jj = (task / lambda_u_tasks) % rho0_tasks; /* rho index */
        int kk = task % lambda_u_tasks;                /* lambda index */
        int T = T_list[ii];
        double rho0 = rho0_min + jj * rho0_step;
        double lambda_u = lambda_u_min + kk * lambda_u_step;
        DynamicsConfig dynamics = {
            .T = T,
            .rho0 = rho0,
            .lambda_u = lambda_u,
            .mu = mu,
        };

        /* g_u(theta)=0 を満たす theta in [0,1] の根を列挙し、最大を採用 */
        double roots[MAX_GD_ROOTS];
        double roots_prime[MAX_GD_ROOTS];
        double roots_prime_prime[MAX_GD_ROOTS];

        int n_roots = 0;
        int n_roots_prime = 0;
        int n_roots_prime_prime = 0;

        double valid_theta_u = 0.0;

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
            double m = g_u_prime(D, &dynamics, Binom, 1.0);
            if (fabs(m) < 1e-10) {
                n_roots = 1;
                roots[0] = 0.0;
                n_roots_prime = 1;
                roots_prime[0] = 1.0;
            } else if (m < 0.0) {
                n_roots = 1;
                roots[0] = 1.0;
            } else {
                n_roots = find_roots(g_u, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000, roots,
                                     MAX_GD_ROOTS);
                while (n_roots == 0 && tmin > 0.0) {
                    tmin *= 0.9;
                    n_roots = find_roots(g_u, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000,
                                         roots, MAX_GD_ROOTS);
                }
                for (int r = 0; r < n_roots; r++) {
                    if (roots[r] < 1.0 && roots[r] >= valid_theta_u) {
                        valid_theta_u = roots[r];
                    }
                }
            }
        } else {
            n_roots = find_roots(g_u, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000, roots,
                                 MAX_GD_ROOTS);
            while (n_roots == 0 && tmin > 0.0) {
                tmin *= 0.9;
                n_roots = find_roots(g_u, D, &dynamics, Binom, tmin, tmax, step, 1e-12, 1000, roots,
                                     MAX_GD_ROOTS);
            }
            for (int r = 0; r < n_roots; r++) {
                if (roots[r] <= 1.0 && roots[r] >= valid_theta_u) {
                    valid_theta_u = roots[r];
                }
            }
        }

        if (valid_theta_u == 0.0 && T == 1) {
            valid_theta_u = 1.0;
        }
        /* 平衡時は A=0 なので R_inf = 1 - S_inf */
        double S_inf = compute_S(D, &dynamics, Binom, valid_theta_u);
        double R_inf = 1.0 - S_inf;
        double Phi_inf = rhs_Phi(D, &dynamics, Binom, valid_theta_u);

        n_roots_prime = find_roots(g_u_prime, D, &dynamics, Binom, 0.0, 1.0, theta_search_step,
                                   1e-6, 100, roots_prime, MAX_GD_ROOTS);

        if (T > 1) {
            n_roots_prime_prime =
                find_roots(g_u_prime_prime, D, &dynamics, Binom, 0.0, 1.0, theta_search_step, 1e-6,
                           100, roots_prime_prime, MAX_GD_ROOTS);
        }

        bool prime_found = false;
        double delta_min = 0.02; // InPow200*200: 0.030
        for (int r = 0; r < n_roots_prime; r++) {
            double root = roots_prime[r];
            double delta = fabs(root - valid_theta_u);
            double g_u_check = g_u(D, &dynamics, Binom, root);
            if (root <= 1.0 && root >= 0.0 && delta < delta_min && fabs(g_u_check) < 1e-3) {
                prime_found = true;
                delta_min = delta;
            }
        }

        bool prime_prime_found = false;
        double delta_min_prime_prime = 0.08; // OutPow200*200: 0.04
        if (prime_found && T > 1) {
            for (int r = 0; r < n_roots_prime_prime; r++) {
                double root = roots_prime_prime[r];
                double delta = fabs(root - valid_theta_u);
                double g_u_check = g_u(D, &dynamics, Binom, root);
                double g_u_prime_check = g_u_prime(D, &dynamics, Binom, root);
                if (root <= 1.0 && root >= 0.0 && delta < delta_min_prime_prime &&
                    fabs(g_u_prime_check) < 1e-2 * 4 && fabs(g_u_check) < 1e-2 * 4) {
                    prime_prime_found = true;
                    delta_min_prime_prime = delta;
                }
            }
        }

#pragma omp critical
        {
            fprintf(fp_gu, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_u, mu,
                    valid_theta_u, R_inf, Phi_inf);
            for (int r = 0; r < n_roots_prime; r++) {
                double root = roots_prime[r];
                if (root <= 1.0 && root >= 0.0) {
                    fprintf(fp_gup, "%d,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_u, mu, root);
                }
            }

            for (int r = 0; r < n_roots_prime_prime; r++) {
                if (roots_prime_prime[r] <= 1.0 && roots_prime_prime[r] >= 0.0) {
                    fprintf(fp_gupp, "%d,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_u, mu,
                            roots_prime_prime[r]);
                }
            }

            if (prime_found) {
                fprintf(fp_gu_gup_zero, "%d,%.6f,%.6f,%.6f,%.6f\n", T, rho0, lambda_u, mu,
                        valid_theta_u);
            }
            if (prime_prime_found) {
                fprintf(fp_gu_gup_gupp_zero, "%d,%.6f,%.6f,%.6f\n", T, rho0, lambda_u, mu);
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

    fclose(fp_gu);
    fclose(fp_gup);
    fclose(fp_gupp);
    fclose(fp_gu_gup_zero);
    fclose(fp_gu_gup_gupp_zero);
    free(Binom);
    free(D);
    return 0;
}
