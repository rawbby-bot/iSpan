#include "../src_mpi/graph.h"
#include "../src_mpi/scc_common.h"
#include "./test_util.h"

#include <filesystem>
#include <vector>

int
main()
{
  // forward graph
  constexpr auto fn_fw_adjacent = "toy2_fw_adjacent.bin";
  constexpr auto fn_fw_begin = "toy2_fw_begin.bin";
  constexpr auto fn_fw_head = "toy2_fw_head.bin";
  constexpr auto fn_out_degree = "toy2_out_degree.bin";
  ASSERT(std::filesystem::exists(fn_fw_adjacent));
  ASSERT(std::filesystem::exists(fn_fw_begin));
  ASSERT(std::filesystem::exists(fn_fw_head));
  ASSERT(std::filesystem::exists(fn_out_degree));

  // backward graph
  constexpr auto fn_bw_adjacent = "toy2_bw_adjacent.bin";
  constexpr auto fn_bw_begin = "toy2_bw_begin.bin";
  constexpr auto fn_bw_head = "toy2_bw_head.bin";
  constexpr auto fn_in_degree = "toy2_in_degree.bin";
  ASSERT(std::filesystem::exists(fn_bw_adjacent));
  ASSERT(std::filesystem::exists(fn_bw_begin));
  ASSERT(std::filesystem::exists(fn_bw_head));
  ASSERT(std::filesystem::exists(fn_in_degree));

  MPI_Init(nullptr, nullptr);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto g = graph(fn_fw_begin, fn_fw_adjacent, fn_bw_begin, fn_bw_adjacent);
  ASSERT(g.vert_count == 7);
  ASSERT(g.edge_count == 11);

  ASSERT(g.fw_beg_pos[0] == 0);
  ASSERT(g.fw_beg_pos[1] == g.fw_beg_pos[0] + 1);
  ASSERT(g.fw_beg_pos[2] == g.fw_beg_pos[1] + 3);
  ASSERT(g.fw_beg_pos[3] == g.fw_beg_pos[2] + 1);
  ASSERT(g.fw_beg_pos[4] == g.fw_beg_pos[3] + 2);
  ASSERT(g.fw_beg_pos[5] == g.fw_beg_pos[4] + 2);
  ASSERT(g.fw_beg_pos[6] == g.fw_beg_pos[5] + 1);
  ASSERT(g.fw_beg_pos[7] == g.fw_beg_pos[6] + 1);

  ASSERT(g.fw_csr[0x0] == 2); // 0 2
  ASSERT(g.fw_csr[0x1] == 0); // 1 0
  ASSERT(g.fw_csr[0x2] == 2); // 1 2
  ASSERT(g.fw_csr[0x3] == 3); // 1 3
  ASSERT(g.fw_csr[0x4] == 0); // 2 0
  ASSERT(g.fw_csr[0x5] == 2); // 3 2
  ASSERT(g.fw_csr[0x6] == 4); // 3 4
  ASSERT(g.fw_csr[0x7] == 1); // 4 1
  ASSERT(g.fw_csr[0x8] == 3); // 4 3
  ASSERT(g.fw_csr[0x9] == 6); // 5 6
  ASSERT(g.fw_csr[0xa] == 5); // 6 5

  ASSERT(g.bw_beg_pos[0] == 0);
  ASSERT(g.bw_beg_pos[1] == g.bw_beg_pos[0] + 2);
  ASSERT(g.bw_beg_pos[2] == g.bw_beg_pos[1] + 1);
  ASSERT(g.bw_beg_pos[3] == g.bw_beg_pos[2] + 3);
  ASSERT(g.bw_beg_pos[4] == g.bw_beg_pos[3] + 2);
  ASSERT(g.bw_beg_pos[5] == g.bw_beg_pos[4] + 1);
  ASSERT(g.bw_beg_pos[6] == g.bw_beg_pos[5] + 1);
  ASSERT(g.bw_beg_pos[7] == g.bw_beg_pos[6] + 1);

  ASSERT(g.bw_csr[0x0] == 1); // 1 0
  ASSERT(g.bw_csr[0x1] == 2); // 2 0
  ASSERT(g.bw_csr[0x2] == 4); // 4 1
  ASSERT(g.bw_csr[0x3] == 0); // 0 2
  ASSERT(g.bw_csr[0x4] == 1); // 1 2
  ASSERT(g.bw_csr[0x5] == 3); // 3 2
  ASSERT(g.bw_csr[0x6] == 1); // 1 3
  ASSERT(g.bw_csr[0x7] == 4); // 4 3
  ASSERT(g.bw_csr[0x8] == 3); // 3 4
  ASSERT(g.bw_csr[0x9] == 6); // 6 5
  ASSERT(g.bw_csr[0xa] == 5); // 5 6

  std::vector<double> avg_time(15);

  if (!world_rank) {

    const auto assignment = prepare_assignment(&g);
    scc_detection(&g, 30, 200, 10, 0.01, 1, avg_time, world_rank, world_size, 1, assignment.get());

    ASSERT(assignment[0] == 0);
    ASSERT(assignment[1] == 1);
    ASSERT(assignment[2] == 0);
    ASSERT(assignment[3] == 1);
    ASSERT(assignment[4] == 1);
    ASSERT(assignment[5] == 5);
    ASSERT(assignment[6] == 5);

  } else {

    scc_detection(&g, 30, 200, 10, 0.01, 1, avg_time, world_rank, world_size, 1, nullptr);
  }

  MPI_Finalize();
}
