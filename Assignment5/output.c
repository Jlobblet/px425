#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "output.h"

void DecompResults_add_times(DecompResults* decomp_results, double start, double search, double merge, double span) {
    decomp_results->cluster_search_time = search - start;
    decomp_results->cluster_merge_time = merge - search;
    decomp_results->cluster_span_time = span - merge;
}

void create_RunResults(RunResults* run_results) {
    run_results->decomp_results = calloc(run_results->n_cell_sizes, sizeof(DecompResults));
    if (run_results->decomp_results == NULL) {
        fprintf(stderr, "Failed to allocate memory for decomposition results array.");
        abort();
    }
}

void destroy_RunResults(RunResults* run_results) {
    if (run_results->decomp_results) {
        free(run_results->decomp_results);
    }
}

void create_Results(Results* results) {
    results->run_results = calloc(results->n_runs, sizeof(RunResults));
    if (results->run_results == NULL) {
        fprintf(stderr, "Failed to allocate memory for run results array.");
        abort();
    }
}

void destroy_Results(Results* results) {
    if (results->run_results) {
        for (int i = 0; i < results->n_runs; i++) {
            destroy_RunResults(&results->run_results[i]);
        }
        free(results->run_results);
    }
}

void print_Results(const Results* results) {
    for (int i = 0; i < results->n_runs; i++) {
        RunResults* run_results = &results->run_results[i];
        print_RunResults(run_results);
    }
}

void print_RunResults(const RunResults* run_results) {
    bool single_size = run_results->n_cell_sizes == 1;
    printf("S = %6.2f P = %8.6f ", run_results->S, run_results->P);
    if (!single_size) { printf("\n"); }
    for (int j = 0; j < run_results->n_cell_sizes; j++) {
        DecompResults decomp_results = run_results->decomp_results[j];
        if (!single_size) { printf("ncells = %3d ", decomp_results.ncells); }
        printf("(%8.4f %8.4f %8.4f sec) ",
               decomp_results.cluster_search_time,
               decomp_results.cluster_merge_time,
               decomp_results.cluster_span_time);
        printf(": %6d clusters, ", decomp_results.n_clusters);
        if (decomp_results.spanning_cluster) {
            printf("%7d spans", decomp_results.spanning_cluster);
        } else {
            printf("none spanning");
        }
        printf("\n");
    }
}
