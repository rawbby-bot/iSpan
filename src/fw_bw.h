#ifndef FW_BW_H
#define FW_BW_H

#include "graph.h"
#include "util.h"
#include "wtime.h"
#include <iostream>
#include <set>
#include <trim_1_gfq.h>

inline void
fw_bfs(
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  vertex_t* fw_sa,
  index_t* front_comm,
  index_t* work_comm,
  vertex_t root,
  index_t tid,
  index_t thread_count,
  double alpha,
  double beta,
  vertex_t edge_count,
  vertex_t vert_count)
{
  depth_t level = 0;
  fw_sa[root] = 0;
  bool is_top_down = true;
  bool is_top_down_queue = false;
  index_t queue_size = vert_count / thread_count;
  while (true) {
    double ltm = wtime();
    index_t front_count = 0;
    index_t my_work_next = 0;
    index_t my_work_curr = 0;

    if (is_top_down) {
      for (vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++) {
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] == level) {
          index_t my_beg = fw_beg_pos[vert_id];
          index_t my_end = fw_beg_pos[vert_id + 1];

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = fw_csr[my_beg];
            if (scc_id[nebr] == 0 && fw_sa[nebr] == -1) {
              fw_sa[nebr] = level + 1;
              my_work_next += fw_beg_pos[nebr + 1] - fw_beg_pos[nebr];

              front_count++;
            }
          }
        }
      }
      work_comm[tid] = my_work_next;
    } else if (!is_top_down_queue) {
      for (vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++) {
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] == -1) {
          index_t my_beg = bw_beg_pos[vert_id];
          index_t my_end = bw_beg_pos[vert_id + 1];
          my_work_curr += my_end - my_beg;

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = bw_csr[my_beg];
            if (scc_id[vert_id] == 0 && fw_sa[nebr] != -1) {
              fw_sa[vert_id] = level + 1;
              front_count++;
              break;
            }
          }
        }
      }
      work_comm[tid] = my_work_curr;

    } else {
      index_t* q = new index_t[queue_size];
      index_t head = 0;
      index_t tail = 0;

      for (vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++) {
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] == level) {
          q[tail++] = vert_id;
        }
      }
      while (head != tail) {
        vertex_t temp_v = q[head++];
        if (head == queue_size)
          head = 0;
        index_t my_beg = fw_beg_pos[temp_v];
        index_t my_end = fw_beg_pos[temp_v + 1];

        for (; my_beg < my_end; ++my_beg) {
          index_t w = fw_csr[my_beg];

          if (scc_id[w] == 0 && fw_sa[w] == -1) {
            q[tail++] = w;
            if (tail == queue_size)
              tail = 0;
            fw_sa[w] = level + 1;
          }
        }
      }
      delete[] q;
      break;
    }

    front_comm[tid] = front_count;

#pragma omp barrier
    front_count = 0;
    my_work_next = 0;

    for (index_t i = 0; i < thread_count; ++i) {
      front_count += front_comm[i];
      my_work_next += work_comm[i];
    }

    if (front_count == 0)
      break;

    if (is_top_down && my_work_next > (alpha * edge_count)) {
      is_top_down = false;
    }
    if (!is_top_down && my_work_next < (beta * edge_count)) {
      is_top_down = true;
    }

#pragma omp barrier

    level++;
  }
}

inline void
bw_bfs(
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  vertex_t* fw_sa,
  vertex_t* bw_sa,
  index_t* front_comm,
  index_t* work_comm,
  vertex_t root,
  index_t tid,
  index_t thread_count,
  double alpha,
  double beta,
  vertex_t edge_count,
  vertex_t vert_count)
{
  bw_sa[root] = 0;
  bool is_top_down = true;
  bool is_top_down_queue = false;
  index_t level = 0;
  scc_id[root] = -9;
  index_t queue_size = vert_count / thread_count;
  while (true) {
    double ltm = wtime();
    index_t front_count = 0;
    index_t my_work_next = 0;
    index_t my_work_curr = 0;

    if (is_top_down) {
      for (vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++) {
        if (scc_id[vert_id] == -9 && bw_sa[vert_id] == level) {
          index_t my_beg = bw_beg_pos[vert_id];
          index_t my_end = bw_beg_pos[vert_id + 1];

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = bw_csr[my_beg];
            if (scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1) {
              bw_sa[nebr] = level + 1;
              my_work_next += bw_beg_pos[nebr + 1] - bw_beg_pos[nebr];
              front_count++;
              scc_id[nebr] = -9;
            }
          }
        }
      }
      work_comm[tid] = my_work_next;
    } else if (!is_top_down_queue) {
      for (vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++) {
        if (scc_id[vert_id] == 0 && bw_sa[vert_id] == -1 && fw_sa[vert_id] != -1) {
          index_t my_beg = fw_beg_pos[vert_id];
          index_t my_end = fw_beg_pos[vert_id + 1];
          my_work_curr += my_end - my_beg;

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = fw_csr[my_beg];
            if (scc_id[nebr] == -9 && bw_sa[nebr] != -1) {
              bw_sa[vert_id] = level + 1;
              front_count++;
              scc_id[vert_id] = -9;
              break;
            }
          }
        }
      }
      work_comm[tid] = my_work_curr;
    } else {

      index_t* q = new index_t[queue_size];
      index_t head = 0;
      index_t tail = 0;

      for (vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++) {
        if (scc_id[vert_id] == -9 && bw_sa[vert_id] == level) {

          q[tail++] = vert_id;
        }
      }

      while (head != tail) {

        vertex_t temp_v = q[head++];
        if (head == queue_size)
          head = 0;
        index_t my_beg = bw_beg_pos[temp_v];
        index_t my_end = bw_beg_pos[temp_v + 1];

        for (; my_beg < my_end; ++my_beg) {
          index_t w = bw_csr[my_beg];

          if (scc_id[w] == 0 && bw_sa[w] == -1 && fw_sa[w] != -1) {

            q[tail++] = w;
            if (tail == queue_size)
              tail = 0;
            scc_id[w] = -9;
            bw_sa[w] = level + 1;
          }
        }
      }

      delete[] q;
    }

    front_comm[tid] = front_count;

#pragma omp barrier
    front_count = 0;
    my_work_next = 0;

    for (index_t i = 0; i < thread_count; ++i) {
      front_count += front_comm[i];
      my_work_next += work_comm[i];
    }

    if (front_count == 0) {
      break;
    }

    if (is_top_down && my_work_next > (alpha * edge_count)) {
      is_top_down = false;
    }
    if (!is_top_down && my_work_next < (beta * edge_count)) {
      is_top_down = true;
    }

#pragma omp barrier

    level++;
  }
}

inline void
fw_bfs_fq_queue(
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  vertex_t* fw_sa,
  index_t* vertex_cur,
  index_t* vertex_front,
  vertex_t root,
  index_t tid,
  index_t thread_count,
  const int alpha,
  const int beta,
  const int gamma,
  vertex_t* frontier_queue,
  vertex_t fq_size,
  const double avg_degree,
  vertex_t vertex_visited,
  vertex_t* temp_queue,
  index_t* prefix_sum,
  vertex_t upper_bound,
  vertex_t* thread_queue)
{
  depth_t level = 0;
  fw_sa[root] = 0;
  temp_queue[0] = root;
  vertex_t queue_size = 1;

  bool is_top_down = true;

  bool is_top_down_queue = false;

#pragma omp barrier
  while (true) {
    double ltm = wtime();
    vertex_t vertex_frontier = 0;

#pragma omp barrier
    if (is_top_down) {
      vertex_t step = queue_size / thread_count;
      vertex_t queue_beg = tid * step;
      vertex_t queue_end = (tid == thread_count - 1 ? queue_size : queue_beg + step);
      for (vertex_t q_vert_id = queue_beg; q_vert_id < queue_end; q_vert_id++) {
        vertex_t vert_id = temp_queue[q_vert_id];

        if (scc_id[vert_id] == 0 && fw_sa[vert_id] != -1) {
          index_t my_beg = fw_beg_pos[vert_id];
          index_t my_end = fw_beg_pos[vert_id + 1];
          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = fw_csr[my_beg];
            if (scc_id[nebr] == 0 && fw_sa[nebr] == -1) {
              fw_sa[nebr] = level + 1;
              thread_queue[vertex_frontier] = nebr;
              vertex_frontier++;
            }
          }
        }
      }
    } else if (!is_top_down_queue) {
      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] == -1) {
          index_t my_beg = bw_beg_pos[vert_id];
          index_t my_end = bw_beg_pos[vert_id + 1];

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = bw_csr[my_beg];
            if (scc_id[nebr] == 0 && fw_sa[nebr] != -1)

            {
              fw_sa[vert_id] = level + 1;
              vertex_frontier++;

              break;
            }
          }
        }
      }
    } else {
      vertex_t end_queue = upper_bound;
      index_t head = 0;
      index_t tail = 0;

      vertex_t step = queue_size / thread_count;
      vertex_t queue_beg = tid * step;
      vertex_t queue_end = (tid == thread_count - 1 ? queue_size : queue_beg + step);

      for (vertex_t q_vert_id = queue_beg; q_vert_id < queue_end; q_vert_id++) {
        thread_queue[tail] = temp_queue[q_vert_id];
        tail++;
      }
      while (head != tail) {
        vertex_t temp_v = thread_queue[head++];

        if (head == end_queue)
          head = 0;
        index_t my_beg = fw_beg_pos[temp_v];
        index_t my_end = fw_beg_pos[temp_v + 1];

        for (; my_beg < my_end; ++my_beg) {
          index_t w = fw_csr[my_beg];

          if (scc_id[w] == 0 && fw_sa[w] == -1) {
            thread_queue[tail++] = w;
            if (tail == end_queue)
              tail = 0;
            fw_sa[w] = level + 1;
          }
        }
      }
    }

    vertex_front[tid] = vertex_frontier;
#pragma omp barrier
    vertex_frontier = 0;

    for (index_t i = 0; i < thread_count; ++i) {
      vertex_frontier += vertex_front[i];
    }
    vertex_visited += vertex_frontier;

    if (VERBOSE) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
      if (tid == 0 && level < 50)
        std::cout << "Level-" << (int)level << " "

                  << vertex_frontier << " "
                  << fq_size << " "
                  << (double)(fq_size) / vertex_frontier << " "
                  << (wtime() - ltm) * 1000 << "ms "
                  << vertex_visited << " "
                  << edge_frontier << " "
                  << edge_remaider << " "
                  << edge_remaider / edge_frontier << "\n";
    }

    if (vertex_frontier == 0)
      break;

    if (is_top_down) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;

      if (!is_top_down_queue && (edge_remainder / alpha) < edge_frontier) {
        is_top_down = false;
        if (VERBOSE) {
          if (tid == 0) {

            std::cout << "--->Switch to bottom up\n";
          }
        }
      }

    } else if ((!is_top_down && !is_top_down_queue && (fq_size * 1.0 / beta) > vertex_frontier) || (!is_top_down && !is_top_down_queue && level > gamma))

    {

      vertex_frontier = 0;
      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] == level + 1) {
          thread_queue[vertex_frontier] = vert_id;
          vertex_frontier++;
        }
      }
      vertex_front[tid] = vertex_frontier;

      is_top_down = false;
      is_top_down_queue = true;

      if (VERBOSE) {
        if (tid == 0)
          std::cout << "--->Switch to top down queue\n";
      }
    }

#pragma omp barrier

    if (is_top_down || is_top_down_queue) {
      get_queue(thread_queue,
                vertex_front,
                prefix_sum,
                tid,
                temp_queue);
      queue_size = prefix_sum[thread_count - 1] + vertex_front[thread_count - 1];
    }
#pragma omp barrier

    level++;
  }
}

inline void
fw_bfs_fq(
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  vertex_t* fw_sa,
  index_t* vertex_cur,
  index_t* vertex_front,
  vertex_t root,
  index_t tid,
  index_t thread_count,
  double alpha,
  double beta,
  vertex_t* frontier_queue,
  vertex_t fq_size,
  const double avg_degree,
  vertex_t vertex_visited)
{
  depth_t level = 0;
  fw_sa[root] = 0;
  vertex_t root_out_degree = fw_beg_pos[root + 1] - fw_beg_pos[root];
  bool is_top_down = true;
  bool is_top_down_async = false;
  if (VERBOSE) {
    if (tid == 0) {
      printf("out_degree, %lu, limit, %.3lf\n", root_out_degree, alpha * beta * fq_size);
    }
  }
  if (root_out_degree < alpha * beta * fq_size) {
    is_top_down_async = true;
  }
  bool is_top_down_queue = false;
  index_t queue_size = fq_size / thread_count;
#pragma omp barrier

  while (true) {

    double ltm = wtime();

    vertex_t vertex_frontier = 0;

    if (is_top_down) {
      if (is_top_down_async) {
        for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
          vertex_t vert_id = frontier_queue[fq_vert_id];

          if (scc_id[vert_id] == 0 && (fw_sa[vert_id] == level || fw_sa[vert_id] == level + 1)) {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id + 1];

            for (; my_beg < my_end; my_beg++) {
              vertex_t nebr = fw_csr[my_beg];
              if (scc_id[nebr] == 0 && fw_sa[nebr] == -1) {
                fw_sa[nebr] = level + 1;

                vertex_frontier++;
              }
            }
          }
        }
      } else {
        for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
          vertex_t vert_id = frontier_queue[fq_vert_id];

          if (scc_id[vert_id] == 0 && fw_sa[vert_id] == level)

          {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id + 1];

            for (; my_beg < my_end; my_beg++) {
              vertex_t nebr = fw_csr[my_beg];
              if (scc_id[nebr] == 0 && fw_sa[nebr] == -1) {
                fw_sa[nebr] = level + 1;

                vertex_frontier++;
              }
            }
          }
        }
      }

    } else if (!is_top_down_queue) {
      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] == -1) {
          index_t my_beg = bw_beg_pos[vert_id];
          index_t my_end = bw_beg_pos[vert_id + 1];

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = bw_csr[my_beg];
            if (scc_id[vert_id] == 0 && fw_sa[nebr] != -1) {
              fw_sa[vert_id] = level + 1;
              vertex_frontier++;

              break;
            }
          }
        }
      }
    } else {
      index_t* q = new index_t[queue_size];
      index_t head = 0;
      index_t tail = 0;

      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] == level) {
          q[tail++] = vert_id;
        }
      }
      while (head != tail) {
        vertex_t temp_v = q[head++];

        if (head == queue_size)
          head = 0;
        index_t my_beg = fw_beg_pos[temp_v];
        index_t my_end = fw_beg_pos[temp_v + 1];

        for (; my_beg < my_end; ++my_beg) {
          index_t w = fw_csr[my_beg];

          if (scc_id[w] == 0 && fw_sa[w] == -1) {
            q[tail++] = w;
            if (tail == queue_size)
              tail = 0;
            fw_sa[w] = level + 1;
          }
        }
      }
      delete[] q;
    }

    vertex_front[tid] = vertex_frontier;

#pragma omp barrier
    vertex_frontier = 0;

    for (index_t i = 0; i < thread_count; ++i) {
      vertex_frontier += vertex_front[i];
    }
    vertex_visited += vertex_frontier;

    if (VERBOSE) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
      if (tid == 0)
        std::cout << "Level-" << (int)level << " "

                  << vertex_frontier << " "
                  << fq_size << " "
                  << (double)(fq_size) / vertex_frontier << " "
                  << (wtime() - ltm) * 1000 << "ms "
                  << vertex_visited << " "
                  << edge_frontier << " "
                  << edge_remaider << " "
                  << edge_remaider / edge_frontier << "\n";
    }
    if (vertex_frontier == 0)
      break;

#pragma omp barrier

    if (is_top_down) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
      if ((edge_remainder / alpha) < edge_frontier) {
        is_top_down = false;
        if (VERBOSE) {
          if (tid == 0)
            std::cout << "--->Switch to bottom up\n";
        }
      }

    } else if (!is_top_down && !is_top_down_queue && (fq_size * 1.0 / beta) > vertex_frontier) {
      is_top_down_queue = true;
      if (VERBOSE) {
        if (tid == 0)
          std::cout << "--->Switch to top down queue\n";
      }
    }
#pragma omp barrier
    level++;
  }
}

inline void
bw_bfs_fq_queue(
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  vertex_t* fw_sa,
  vertex_t* bw_sa,
  index_t* vertex_cur,
  index_t* vertex_front,
  vertex_t root,
  index_t tid,
  index_t thread_count,
  const int alpha,
  const int beta,
  const int gamma,
  vertex_t* frontier_queue,
  vertex_t fq_size,
  const double avg_degree,
  vertex_t vertex_visited,
  vertex_t* temp_queue,
  index_t* prefix_sum,
  vertex_t upper_bound,
  vertex_t* thread_queue)
{
  depth_t level = 0;
  bw_sa[root] = 0;
  scc_id[root] = -9;
  temp_queue[0] = root;
  vertex_t queue_size = 1;

  if (VERBOSE) {
    if (tid == 0) {
      printf("upperbound, %lu\n", upper_bound);
    }
  }

  bool is_top_down = true;

  bool is_top_down_queue = false;

#pragma omp barrier
  while (true) {
    double ltm = wtime();
    vertex_t vertex_frontier = 0;

#pragma omp barrier
    if (is_top_down) {
      vertex_t step = queue_size / thread_count;
      vertex_t queue_beg = tid * step;
      vertex_t queue_end = (tid == thread_count - 1 ? queue_size : queue_beg + step);

      for (vertex_t q_vert_id = queue_beg; q_vert_id < queue_end; q_vert_id++) {
        vertex_t vert_id = temp_queue[q_vert_id];

        if (scc_id[vert_id] == -9) {
          index_t my_beg = bw_beg_pos[vert_id];
          index_t my_end = bw_beg_pos[vert_id + 1];

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = bw_csr[my_beg];
            if (scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1) {
              bw_sa[nebr] = level + 1;
              thread_queue[vertex_frontier] = nebr;
              vertex_frontier++;
              scc_id[nebr] = -9;
            }
          }
        }
      }
    } else if (!is_top_down_queue) {
      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if (scc_id[vert_id] == 0 && fw_sa[vert_id] != -1) {
          index_t my_beg = fw_beg_pos[vert_id];
          index_t my_end = fw_beg_pos[vert_id + 1];

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = fw_csr[my_beg];
            if (scc_id[nebr] == -9)

            {
              bw_sa[vert_id] = level + 1;
              vertex_frontier++;
              scc_id[vert_id] = -9;

              break;
            }
          }
        }
      }
    } else {
      vertex_t end_queue = upper_bound;
      index_t head = 0;
      index_t tail = 0;

      vertex_t step = queue_size / thread_count;
      vertex_t queue_beg = tid * step;
      vertex_t queue_end = (tid == thread_count - 1 ? queue_size : queue_beg + step);

      for (vertex_t q_vert_id = queue_beg; q_vert_id < queue_end; q_vert_id++) {
        thread_queue[tail] = temp_queue[q_vert_id];
        tail++;
      }
      while (head != tail) {
        vertex_t temp_v = thread_queue[head++];

        if (head == end_queue)
          head = 0;
        index_t my_beg = bw_beg_pos[temp_v];
        index_t my_end = bw_beg_pos[temp_v + 1];

        for (; my_beg < my_end; ++my_beg) {
          index_t w = bw_csr[my_beg];

          if (scc_id[w] == 0 && bw_sa[w] == -1 && fw_sa[w] != -1) {
            thread_queue[tail++] = w;
            if (tail == end_queue)
              tail = 0;
            scc_id[w] = -9;
            bw_sa[w] = level + 1;
          }
        }
      }
    }

    vertex_front[tid] = vertex_frontier;
#pragma omp barrier
    vertex_frontier = 0;

    for (index_t i = 0; i < thread_count; ++i) {
      vertex_frontier += vertex_front[i];
    }
    vertex_visited += vertex_frontier;

    if (VERBOSE) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
      if (tid == 0 && level < 50)
        std::cout << "Level-" << (int)level << " "

                  << vertex_frontier << " "
                  << fq_size << " "
                  << (double)(fq_size) / vertex_frontier << " "
                  << (wtime() - ltm) * 1000 << "ms "
                  << vertex_visited << " "
                  << edge_frontier << " "
                  << edge_remaider << " "
                  << edge_remaider / edge_frontier << "\n";
    }

    if (vertex_frontier == 0)
      break;

    if (is_top_down) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
      if (!is_top_down_queue && (edge_remainder / alpha) < edge_frontier) {
        is_top_down = false;
        if (VERBOSE) {
          if (tid == 0)
            std::cout << "--->Switch to bottom up\n";
        }
      }

    } else if ((!is_top_down && !is_top_down_queue && (fq_size * 1.0 / beta) > vertex_frontier) || (!is_top_down && !is_top_down_queue && level > gamma))

    {
      vertex_frontier = 0;
      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if (scc_id[vert_id] == -9 && bw_sa[vert_id] == level + 1) {
          thread_queue[vertex_frontier] = vert_id;
          vertex_frontier++;
        }
      }
      vertex_front[tid] = vertex_frontier;

      is_top_down = false;
      is_top_down_queue = true;
      if (VERBOSE) {
        if (tid == 0)
          std::cout << "--->Switch to top down queue\n";
      }
    }

#pragma omp barrier

    if (is_top_down || is_top_down_queue) {
      get_queue(thread_queue,
                vertex_front,
                prefix_sum,
                tid,
                temp_queue);
      queue_size = prefix_sum[thread_count - 1] + vertex_front[thread_count - 1];
    }

#pragma omp barrier

    level++;
  }
}

inline void
bw_bfs_fq(
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  vertex_t* fw_sa,
  vertex_t* bw_sa,
  index_t* vertex_cur,
  index_t* vertex_front,
  vertex_t root,
  index_t tid,
  index_t thread_count,
  double alpha,
  double beta,
  vertex_t* frontier_queue,
  vertex_t fq_size,
  const double avg_degree,
  vertex_t vertex_visited)
{
  bw_sa[root] = 0;
  vertex_t root_in_degree = bw_beg_pos[root + 1] - bw_beg_pos[root];
  bool is_top_down = true;
  bool is_top_down_queue = false;
  bool is_top_down_async = false;
  if (DEBUG) {
    if (tid == 0) {
      printf("in_degree, %lu, limit, %.3lf\n", root_in_degree, alpha * beta * fq_size);
    }
  }
  if (root_in_degree < alpha * beta * fq_size) {
    is_top_down_async = true;
  }
  index_t level = 0;
  scc_id[root] = -9;
  index_t queue_size = fq_size / thread_count;
  while (true) {
    double ltm = wtime();

    vertex_t vertex_frontier = 0;

    if (is_top_down) {
      if (is_top_down_async) {
        for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
          vertex_t vert_id = frontier_queue[fq_vert_id];
          if (scc_id[vert_id] == -9 && (bw_sa[vert_id] == level || bw_sa[vert_id] == level + 1)) {
            index_t my_beg = bw_beg_pos[vert_id];
            index_t my_end = bw_beg_pos[vert_id + 1];

            for (; my_beg < my_end; my_beg++) {
              vertex_t nebr = bw_csr[my_beg];
              if (scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1) {
                bw_sa[nebr] = level + 1;

                vertex_frontier++;
                scc_id[nebr] = -9;
              }
            }
          }
        }
      } else {
        for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
          vertex_t vert_id = frontier_queue[fq_vert_id];
          if (scc_id[vert_id] == -9 && bw_sa[vert_id] == level) {
            index_t my_beg = bw_beg_pos[vert_id];
            index_t my_end = bw_beg_pos[vert_id + 1];

            for (; my_beg < my_end; my_beg++) {
              vertex_t nebr = bw_csr[my_beg];
              if (scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1) {
                bw_sa[nebr] = level + 1;

                vertex_frontier++;

                scc_id[nebr] = -9;
              }
            }
          }
        }
      }

    } else if (!is_top_down_queue) {
      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];

        if (scc_id[vert_id] == 0) {
          index_t my_beg = fw_beg_pos[vert_id];
          index_t my_end = fw_beg_pos[vert_id + 1];

          for (; my_beg < my_end; my_beg++) {
            vertex_t nebr = fw_csr[my_beg];

            if (scc_id[nebr] == -9) {
              bw_sa[vert_id] = level + 1;

              vertex_frontier++;
              scc_id[vert_id] = -9;
              break;
            }
          }
        }
      }

    } else {

      index_t* q = new index_t[queue_size];
      index_t head = 0;
      index_t tail = 0;

      for (vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; fq_vert_id++) {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if (scc_id[vert_id] == -9 && bw_sa[vert_id] == level) {
          q[tail++] = vert_id;
        }
      }

      while (head != tail) {
        vertex_t temp_v = q[head++];

        if (head == queue_size)
          head = 0;
        index_t my_beg = bw_beg_pos[temp_v];
        index_t my_end = bw_beg_pos[temp_v + 1];

        for (; my_beg < my_end; ++my_beg) {
          index_t w = bw_csr[my_beg];

          if (scc_id[w] == 0 && bw_sa[w] == -1 && fw_sa[w] != -1) {
            q[tail++] = w;
            if (tail == queue_size)
              tail = 0;
            scc_id[w] = -9;
            bw_sa[w] = level + 1;
          }
        }
      }
      delete[] q;
      if (!DEBUG)
        break;
    }

    vertex_front[tid] = vertex_frontier;

#pragma omp barrier
    vertex_frontier = 0;

    for (index_t i = 0; i < thread_count; ++i) {
      vertex_frontier += vertex_front[i];
    }
    vertex_visited += vertex_frontier;

    if (VERBOSE) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
      if (tid == 0)
        std::cout << "Level-" << (int)level << " "

                  << vertex_frontier << " "
                  << fq_size << " "
                  << (double)(fq_size) / vertex_frontier << " "
                  << (wtime() - ltm) * 1000 << "ms "
                  << vertex_visited << " "
                  << edge_frontier << " "
                  << edge_remaider << " "
                  << edge_remaider / edge_frontier << "\n";
    }
    if (vertex_frontier == 0)
      break;

#pragma omp barrier

    if (is_top_down) {
      double edge_frontier = (double)vertex_frontier * avg_degree;
      double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
      if ((edge_remainder / alpha) < edge_frontier) {
        is_top_down = false;
        if (VERBOSE) {
          if (tid == 0)
            std::cout << "--->Switch to bottom up\n";
        }
      }

    } else if (!is_top_down && !is_top_down_queue && (fq_size * 1.0 / beta) > vertex_frontier) {
      is_top_down_queue = true;
      if (VERBOSE) {
        if (tid == 0)
          std::cout << "--->Switch to top down queue\n";
      }
    }
#pragma omp barrier
    level++;
  }
  if (tid == 0)
    printf("bw_level, %lu\n", level);
}

inline void
init_fw_sa(
  index_t vert_beg,
  index_t vert_end,
  vertex_t* fw_sa,
  vertex_t* frontier_queue,
  vertex_t fq_size,
  vertex_t* wcc_fq,
  index_t tid,
  vertex_t* color,
  vertex_t& wcc_fq_size)
{
  for (index_t i = vert_beg; i < vert_end; ++i) {
    fw_sa[i] = -1;
  }
  if (tid == 0) {
    std::set<int> s_fq;
    for (vertex_t i = 0; i < fq_size; ++i) {

      s_fq.insert(color[frontier_queue[i]]);
    }
    wcc_fq_size = s_fq.size();
    std::set<int>::iterator it;
    int i = 0;
    for (it = s_fq.begin(); it != s_fq.end(); ++it, ++i) {

      wcc_fq[i] = *it;
    }
  }
}

inline void
mice_fw_bw(
  color_t* wcc_color,
  index_t* scc_id,
  index_t* fw_beg_pos,
  index_t* bw_beg_pos,
  vertex_t* fw_csr,
  vertex_t* bw_csr,
  vertex_t* fw_sa,
  index_t tid,
  index_t thread_count,
  vertex_t* frontier_queue,
  vertex_t fq_size,
  vertex_t* wcc_fq,
  vertex_t wcc_fq_size)
{
  index_t step = wcc_fq_size / thread_count;
  index_t wcc_beg = tid * step;
  index_t wcc_end = (tid == thread_count - 1 ? wcc_fq_size : wcc_beg + step);

  index_t* q = new index_t[fq_size];
  index_t head = 0;
  index_t tail = 0;

  for (index_t fq_i = 0; fq_i < fq_size; ++fq_i) {
    vertex_t v = frontier_queue[fq_i];
    if (scc_id[v] == 0) {
      vertex_t cur_wcc = wcc_color[v];
      bool in_wcc = false;
      for (index_t i = wcc_beg; i < wcc_end; ++i) {
        if (wcc_fq[i] == cur_wcc) {
          in_wcc = true;
          break;
        }
      }
      if (in_wcc) {

        fw_sa[v] = v;
        q[tail++] = v;
        if (tail == fq_size)
          tail = 0;
        while (head != tail) {
          vertex_t temp_v = q[head++];

          if (head == fq_size)
            head = 0;
          index_t my_beg = fw_beg_pos[temp_v];
          index_t my_end = fw_beg_pos[temp_v + 1];

          for (; my_beg < my_end; ++my_beg) {
            index_t w = fw_csr[my_beg];

            if (scc_id[w] == 0 && fw_sa[w] != v) {
              q[tail++] = w;
              if (tail == fq_size)
                tail = 0;
              fw_sa[w] = v;
            }
          }
        }

        scc_id[v] = v;
        q[tail++] = v;
        if (tail == fq_size)
          tail = 0;

        while (head != tail) {
          vertex_t temp_v = q[head++];

          if (head == fq_size)
            head = 0;
          index_t my_beg = bw_beg_pos[temp_v];
          index_t my_end = bw_beg_pos[temp_v + 1];

          for (; my_beg < my_end; ++my_beg) {
            index_t w = bw_csr[my_beg];

            if (scc_id[w] == 0 && fw_sa[w] == v) {
              q[tail++] = w;
              if (tail == fq_size)
                tail = 0;
              scc_id[w] = v;
            }
          }
        }
      }
    }
  }
  delete[] q;
}
#endif
