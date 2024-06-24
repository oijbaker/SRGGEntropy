#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

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
double rayleigh(double r, double r0, double eta) {
    return 1 - exp(-pow(r / r0, eta));
}

// 3
double triangular(double r, double r0) {
    return r < r0 ? r / r0 : 1 - (r - r0) / (1 - r0);
}

// Utils
int binary_to_int(int *arr, int n) {
    int res = 0;
    for (int i = 0; i < n; i++) {
        res += arr[i] * pow(2, n - i - 1);
    }
    return res;
}

// Entropy Functions per Geometry
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