#ifndef PX425_OUTPUT_H
#define PX425_OUTPUT_H

typedef struct DecompResults {
    int ncells;
    double cluster_search_time;
    double cluster_merge_time;
    double cluster_span_time;
    int n_clusters;
    int spanning_cluster;
} DecompResults;

void DecompResults_add_times(DecompResults* decomp_results, double start, double search, double merge, double span);

typedef struct RunResults {
    double S;
    double P;
    int n_cell_sizes;
    DecompResults* decomp_results;
} RunResults;

void create_RunResults(RunResults* run_results);
void destroy_RunResults(RunResults* run_results);

typedef struct Results {
    int nruns;
    RunResults* run_results;
} Results;

void create_Results(Results* results);
void destroy_Results(Results* results);

void print_Results(const Results* results);

#endif //PX425_OUTPUT_H
