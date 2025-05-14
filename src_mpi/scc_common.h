#ifndef SCC_COMMON_H
#define SCC_COMMON_H
#include "graph.h"
#include "util.h"
#include "wtime.h"

#include <memory>

inline auto
prepare_assignment(const graph* g)
{
  const auto entries = static_cast<std::size_t>(g->vert_count + 1);
  return std::make_unique<vertex_t[]>(entries);
}

graph*
graph_load(
  const char* fw_beg_file,
  const char* fw_csr_file,
  const char* bw_beg_file,
  const char* bw_csr_file,
  double* avg_time);

void
scc_detection(
  const graph* g,
  int alpha,
  int beta,
  int gamma,
  double theta,
  index_t thread_count,
  double* avg_time,
  int world_rank,
  int world_size,
  int run_time,
  vertex_t* assignment = nullptr);

void
get_scc_result(
  index_t* scc_id,
  index_t vert_count);

void
print_time_result(
  index_t run_times,
  double* avg_time);

inline index_t
pivot_selection(
  index_t* scc_id,
  long_t* fw_beg_pos,
  long_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  index_t* max_pivot_list,
  index_t* max_degree_list,
  index_t tid,
  index_t thread_count)
{
  index_t max_pivot_thread = 0;
  index_t max_degree_thread = 0;
  for (vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id) {
    if (scc_id[vert_id] == 0) {
      long_t out_degree = fw_beg_pos[vert_id + 1] - fw_beg_pos[vert_id];
      long_t in_degree = bw_beg_pos[vert_id + 1] - bw_beg_pos[vert_id];
      long_t degree_mul = out_degree * in_degree;

      if (degree_mul > max_degree_thread) {
        max_degree_thread = degree_mul;
        max_pivot_thread = vert_id;
      }
    }
  }
  index_t max_pivot = max_pivot_thread;
  index_t max_degree = max_degree_thread;

  {
    if (tid == 0)
      printf("max_pivot, %lu, max_degree, %lu\n", max_pivot, max_degree);
  }
  return max_pivot;
}

inline index_t
pivot_selection_from_fq(
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  index_t* max_pivot_list,
  index_t* max_degree_list,
  index_t tid,
  index_t thread_count,
  index_t* small_queue

)
{
  index_t max_pivot_thread = 0;
  index_t max_degree_thread = 0;
  for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id) {
    vertex_t vert_id = small_queue[fq_vert_id];

    index_t degree_mul = (fw_beg_pos[vert_id + 1] - fw_beg_pos[vert_id]) * (bw_beg_pos[vert_id + 1] - bw_beg_pos[vert_id]);

    if (degree_mul > max_degree_thread) {
      max_degree_thread = degree_mul;
      max_pivot_thread = vert_id;
    }
  }
  max_pivot_list[tid] = max_pivot_thread;
  max_degree_list[tid] = max_degree_thread;

#pragma omp barrier

  index_t max_pivot = 0;
  index_t max_degree = 0;

  for (index_t i = 0; i < thread_count; ++i) {
    if (max_degree_list[i] > max_degree) {
      max_degree = max_degree_list[i];
      max_pivot = max_pivot_list[i];
    }
  }
  if (DEBUG) {
    if (tid == 0)
      printf("max_pivot, %lu, max_degree, %lu\n", max_pivot, max_degree);
  }
  return max_pivot;
}

#endif
