#include "graph.h"
#include "wtime.h"

#include "scc_common.h"
#include <unistd.h>

int
main(int args, char** argv)
{
  printf("args = %d\n", args);
  if (args != 11) {
    std::cout << "Usage: ./scc_cpu <fw_beg_file> <fw_csr_file> <bw_beg_file> <bw_csr_file> <thread_count> <alpha> <beta> <gamma> <theta> <run_times>\n";
    exit(-1);
  }

  const char* fw_beg_file = argv[1];
  const char* fw_csr_file = argv[2];
  const char* bw_beg_file = argv[3];
  const char* bw_csr_file = argv[4];
  const index_t thread_count = atoi(argv[5]);
  const int alpha = atof(argv[6]);
  const int beta = atof(argv[7]);
  const int gamma = atof(argv[8]);
  const double theta = atof(argv[9]);
  const index_t run_times = atoi(argv[10]);
  printf("Thread = %lu, alpha = %d, beta = %d, gamma = %d, theta = %g, run_times = %lu\n", thread_count, alpha, beta, gamma, theta, run_times);

  auto* avg_time = new double[15];

  graph* g = graph_load(fw_beg_file,
                        fw_csr_file,
                        bw_beg_file,
                        bw_csr_file,
                        avg_time);
  index_t i = 0;

  while (i++ < run_times) {
    printf("\nRuntime: %lu\n", i);
    scc_detection(g,
                  alpha,
                  beta,
                  gamma,
                  theta,
                  thread_count,
                  avg_time);
  }

  print_time_result(run_times,
                    avg_time);

  return 0;
}
