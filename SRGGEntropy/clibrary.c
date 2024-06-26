#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// RNG
double unif(double low, double high) {
    return low + (high - low) * rand() / (RAND_MAX + 1.0);
}

// Connection Functions

// 1
double hard(double r, double r0) {
    return r < r0 ? 1 : 0;
}

// 2
double waxman(double r, double r0) {
    return exp(-r / r0);
}

// 3
double rayleigh(double r, double r0, double eta) {
    return exp(-pow(r / r0, eta));
}

// 4
double triangular(double r, double r0) {
    return 1 - r / r0;
}

// 5
double exclusion(double r, double r0) {
    return r > r0 ? 1 : 0;
}

// Utils
int binary_to_int(int *arr, int n) {
    int res = 0;
    for (int i = 0; i < n; i++) {
        res += arr[i] * pow(2, n - i - 1);
    }
    return res;
}

double torus_distance_1d(double x, double y) {
    double dx = fabs(x - y);
    if (dx > 0.5) {
        dx = 1 - dx;
    }
}

double torus_distance_2d(double x[], double y[]) {
    double dx = fabs(x[0] - y[0]);
    double dy = fabs(x[1] - y[1]);

    if (dx > 0.5) {
        dx = 1 - dx;
    }
    if (dy > 0.5) {
        dy = 1 - dy;
    }

    return pow(pow(dx, 2) + pow(dy, 2), 0.5);
}

// Exact Geometric Functions

double f_square(double r) {
    double x = r;
    double u = sqrt(x*x-1);
    double v = sqrt(x*x-1);
    double gamma;

    if (x <= 1) {
        gamma = M_PI-4*x+x*x;
    } else if (x <= 1) {
        gamma = 2*asin(1/x)-1-2*(x-u);
    } else if (x <= sqrt(2)) {
        gamma = 2*asin((1-u*v)/(x*x))+2*u+2*v-1-(1+x*x);
    }
    else {
        return 0;
    }

    return 2*r*gamma;
}

double f_line(double r) {
    return 2-2*r;
}

double f_torus_1D(double r) {
    return 2;
}

double f_torus_2D(double r) {
    if (r <= 0.5) {
        return 2*M_PI*r;
    } else if (r <= sqrt(2)/2) {
        return 2*r*(M_PI - 4*acos(1.0/(2*r)));
    } else {
        return 0;
    }
}

double f_disc(double r) {
    return 4*r/M_PI*acos(r/2)-2*r*r/M_PI*sqrt(1-r*r/4);
}

// Entropy Functions per Geometry
double torus_1d_entropy_hard(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_1d(nodes[j], nodes[k]);
                graph[count] = hard(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double torus_1d_entropy_rayleigh(const int n, double r0, double eta, int lim) {
        const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_1d(nodes[j], nodes[k]);
                if (unif(0,1) < rayleigh(dist, r0, eta)) {
                    graph[count] = 1;
                } else {
                    graph[count] = 0;
                }
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double torus_1d_entropy_triangular(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_1d(nodes[j], nodes[k]);
                if (unif(0,1) < triangular(dist, r0)) {
                    graph[count] = 1;
                } else {
                    graph[count] = 0;
                }
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double torus_1d_entropy_exclusion(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_1d(nodes[j], nodes[k]);
                graph[count] = exclusion(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}


double line_entropy_hard(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = fabs(nodes[j] - nodes[k]);
                graph[count] = hard(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;

}

double line_entropy_rayleigh(const int n, double r0, double eta, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = fabs(nodes[j] - nodes[k]);
                if (unif(0,1) < rayleigh(dist, r0, eta)) {
                    graph[count] = 1;
                } else {
                    graph[count] = 0;
                }
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double line_entropy_triangular(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = fabs(nodes[j] - nodes[k]);
                if (unif(0,1) < triangular(dist, r0)) {
                    graph[count] = 1;
                } else {
                    graph[count] = 0;
                }
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double line_entropy_exclusion(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = fabs(nodes[j] - nodes[k]);
                graph[count] = exclusion(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}


double torus_2d_entropy_hard(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_2d(nodes[j], nodes[k]);
                graph[count] = hard(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;

}

double torus_2d_entropy_rayleigh(const int n, double r0, double eta, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_2d(nodes[j], nodes[k]);
                if (unif(0,1) < rayleigh(dist, r0, eta)) {
                    graph[count] = 1;
                } else {
                    graph[count] = 0;
                }
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;

}

double torus_2d_entropy_triangular(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_2d(nodes[j], nodes[k]);
                if (unif(0,1) < triangular(dist, r0)) {
                    graph[count] = 1;
                } else {
                    graph[count] = 0;
                }
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;

}

double torus_2d_entropy_exclusion(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = torus_distance_2d(nodes[j], nodes[k]);
                graph[count] = exclusion(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;

}


double square_entropy_hard(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = hard(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double square_entropy_rayleigh(const int n, double r0, double eta, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = unif(0,1) < rayleigh(dist, r0, eta) ? 1 : 0;
                count++;
            }
        }  
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double square_entropy_triangular(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = unif(0,1) < triangular(dist, r0) ? 1 : 0;
                count++;
            }
        }  
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double square_entropy_exclusion(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            nodes[j_][0] = unif(0, 1);
            nodes[j_][1] = unif(0, 1);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = exclusion(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}


double disc_entropy_hard(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            double n1r = unif(0,1);
            double n1theta = unif(0,2*M_PI);

            nodes[j_][0] = n1r*cos(n1theta);
            nodes[j_][1] = n1r*sin(n1theta);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = hard(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double disc_entropy_rayleigh(const int n, double r0, double eta, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            double n1r = unif(0,1);
            double n1theta = unif(0,2*M_PI);

            nodes[j_][0] = n1r*cos(n1theta);
            nodes[j_][1] = n1r*sin(n1theta);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = unif(0,1) < rayleigh(dist, r0, eta) ? 1 : 0;
                count++;
            }
        }  
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double disc_entropy_triangular(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            double n1r = unif(0,1);
            double n1theta = unif(0,2*M_PI);

            nodes[j_][0] = n1r*cos(n1theta);
            nodes[j_][1] = n1r*sin(n1theta);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = unif(0,1) < triangular(dist, r0) ? 1 : 0;
                count++;
            }
        }  
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}

double disc_entropy_exculsion(const int n, double r0, int lim) {
    const int n_graphs = pow(2, 0.5*n*(n-1));
    int graph_map[n_graphs];

    for (int i=0; i<n_graphs; i++) {
        graph_map[i] = 0;
    }

    for (int i = 0; i < lim; i++) {
        int graph[n+1];
        double nodes[n][2];

        for (int j_ = 0; j_ < n; j_++) {
            double n1r = unif(0,1);
            double n1theta = unif(0,2*M_PI);

            nodes[j_][0] = n1r*cos(n1theta);
            nodes[j_][1] = n1r*sin(n1theta);
        }

        int count = 0;
        for (int j = 0; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                double dist = pow(pow(nodes[j][0] - nodes[k][0], 2) + pow(nodes[j][1] - nodes[k][1], 2), 0.5);
                graph[count] = exclusion(dist, r0);
                count++;
            }
        }
        graph_map[binary_to_int(graph, n)] += 1;
    }
        
    // calculate entropy
    double entropy = 0;

    for (int i = 0; i < n_graphs; i++) {
        if (graph_map[i] > 0) {
            double p = (double)graph_map[i] / lim;
            entropy += -p * log(p);
        }
    }

    return entropy;
}


// Quadrature for exact f(r)

double quad_p_bar_square_rayleigh(double r0, double eta, int lim) {
    // integrating f_square(r)*rayleigh(r, r0, eta) from 0 to sqrt(2) by Trapezoidal Rule
    double sum = 0;
    double h = sqrt(2) / lim;
    for (int i=0; i<lim; i++) {
        double x = i * h;
        sum += 0.5 * h * (f_square(x)*rayleigh(x, r0, eta) + f_square(x + h)*rayleigh(x + h, r0, eta));
    }
    return sum;
}

double p_bar_rayleigh(char geometry[], double r0, double eta, int lim) {
    if (strcmp(geometry, "square") == 0) {
        double sum = 0;
        double h = sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_square(x)*rayleigh(x, r0, eta) + f_square(x + h)*rayleigh(x + h, r0, eta));
        return sum; 
        }
    } else if (strcmp(geometry, "1d_torus") == 0) {
        double sum = 0;
        double h = 0.5 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_1D(x)*rayleigh(x, r0, eta) + f_torus_1D(x + h)*rayleigh(x + h, r0, eta));
        }
        return sum;
    } else if (strcmp(geometry, "2d_torus") == 0)  {
        double sum = 0;
        double h = 0.5*sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_2D(x)*rayleigh(x, r0, eta) + f_torus_2D(x + h)*rayleigh(x + h, r0, eta));
        }
        return sum;
    } else if (strcmp(geometry, "line") == 0) {
        double sum = 0;
        double h = 1.0 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_line(x)*rayleigh(x, r0, eta) + f_line(x + h)*rayleigh(x + h, r0, eta));
        }
        return sum;
    } else if (strcmp(geometry, "disc") == 0) {
        double sum = 0;
        double h = 2.0 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_disc(x)*rayleigh(x, r0, eta) + f_disc(x + h)*rayleigh(x + h, r0, eta));
        }
        return sum;
    }
    return 0;
}

double p_bar_hard(char geometry[], double r0, int lim) {
    if (strcmp(geometry, "square") == 0) {
        double sum = 0;
        double h = sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_square(x)*hard(x, r0) + f_square(x + h)*hard(x + h, r0));
        return sum; 
        }
    } else if (strcmp(geometry, "1d_torus") == 0) {
        double sum = 0;
        double h = 0.5 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_1D(x)*hard(x, r0) + f_torus_1D(x + h)*hard(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "2d_torus") == 0)  {
        double sum = 0;
        double h = 0.5*sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_2D(x)*hard(x, r0) + f_torus_2D(x + h)*hard(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "line") == 0) {
        double sum = 0;
        double h = 1 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_line(x)*hard(x, r0) + f_line(x + h)*hard(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "disc") == 0) {
        double sum = 0;
        double h = 2 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_disc(x)*hard(x, r0) + f_disc(x + h)*hard(x + h, r0));
        }
        return sum;
    }
    return 0;
}

double p_bar_triangular(char geometry[], double r0, int lim) {
    if (strcmp(geometry, "square") == 0) {
        double sum = 0;
        double h = sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_square(x)*triangular(x, r0) + f_square(x + h)*triangular(x + h, r0));
        return sum; 
        }
    } else if (strcmp(geometry, "1d_torus") == 0) {
        double sum = 0;
        double h = 0.5 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_1D(x)*triangular(x, r0) + f_torus_1D(x + h)*triangular(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "2d_torus") == 0)  {
        double sum = 0;
        double h = 0.5*sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_2D(x)*triangular(x, r0) + f_torus_2D(x + h)*triangular(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "line") == 0) {
        double sum = 0;
        double h = 1 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_line(x)*triangular(x, r0) + f_line(x + h)*triangular(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "disc") == 0) {
        double sum = 0; 
        double h = 2 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_disc(x)*triangular(x, r0) + f_disc(x + h)*triangular(x + h, r0));
        }
        return sum;
    }
    return 0;
}

double p_bar_exclusion(char geometry[], double r0, int lim) {
    if (strcmp(geometry, "square") == 0) {
        double sum = 0;
        double h = sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_square(x)*exclusion(x, r0) + f_square(x + h)*exclusion(x + h, r0));
        return sum; 
        }
    } else if (strcmp(geometry, "1d_torus") == 0) {
        double sum = 0;
        double h = 0.5 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_1D(x)*exclusion(x, r0) + f_torus_1D(x + h)*exclusion(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "2d_torus") == 0)  {
        double sum = 0;
        double h = 0.5*sqrt(2) / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_torus_2D(x)*exclusion(x, r0) + f_torus_2D(x + h)*exclusion(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "line") == 0) {
        double sum = 0;
        double h = 1 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_line(x)*exclusion(x, r0) + f_line(x + h)*exclusion(x + h, r0));
        }
        return sum;
    } else if (strcmp(geometry, "disc") == 0) {
        double sum = 0;
        double h = 2 / lim;
        for (int i=0; i<lim; i++) {
            double x = i * h;
            sum += 0.5 * h * (f_disc(x)*exclusion(x, r0) + f_disc(x + h)*exclusion(x + h, r0));
        }
        return sum;
    }
}
                