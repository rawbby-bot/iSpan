#include "../src/graph.h"
#include "./test_util.h"

#include <filesystem>

int
main()
{
  // forward graph
  constexpr auto fn_fw_adjacent = "fw_adjacent.bin";
  constexpr auto fn_fw_begin = "fw_begin.bin";
  constexpr auto fn_fw_head = "fw_head.bin";
  constexpr auto fn_out_degree = "out_degree.bin";
  ASSERT(std::filesystem::exists(fn_fw_adjacent));
  ASSERT(std::filesystem::exists(fn_fw_begin));
  ASSERT(std::filesystem::exists(fn_fw_head));
  ASSERT(std::filesystem::exists(fn_out_degree));

  // backward graph
  constexpr auto fn_bw_adjacent = "bw_adjacent.bin";
  constexpr auto fn_bw_begin = "bw_begin.bin";
  constexpr auto fn_bw_head = "bw_head.bin";
  constexpr auto fn_in_degree = "in_degree.bin";
  ASSERT(std::filesystem::exists(fn_bw_adjacent));
  ASSERT(std::filesystem::exists(fn_bw_begin));
  ASSERT(std::filesystem::exists(fn_bw_head));
  ASSERT(std::filesystem::exists(fn_in_degree));

  const auto g = graph(fn_fw_begin, fn_fw_adjacent, fn_bw_begin, fn_bw_adjacent);
  ASSERT(g.vert_count == 7);
  ASSERT(g.edge_count == 16);

  ASSERT(g.fw_beg_pos[0] == 0);
  ASSERT(g.fw_beg_pos[1] == g.fw_beg_pos[0] + 2);
  ASSERT(g.fw_beg_pos[2] == g.fw_beg_pos[1] + 4);
  ASSERT(g.fw_beg_pos[3] == g.fw_beg_pos[2] + 3);
  ASSERT(g.fw_beg_pos[4] == g.fw_beg_pos[3] + 3);
  ASSERT(g.fw_beg_pos[5] == g.fw_beg_pos[4] + 2);
  ASSERT(g.fw_beg_pos[6] == g.fw_beg_pos[5] + 1);
  ASSERT(g.fw_beg_pos[7] == g.fw_beg_pos[6] + 1);

  ASSERT(g.fw_csr[0x0] == 1); // 0 1
  ASSERT(g.fw_csr[0x1] == 2); // 0 2
  ASSERT(g.fw_csr[0x2] == 0); // 1 0
  ASSERT(g.fw_csr[0x3] == 2); // 1 2
  ASSERT(g.fw_csr[0x4] == 3); // 1 3
  ASSERT(g.fw_csr[0x5] == 4); // 1 4
  ASSERT(g.fw_csr[0x6] == 0); // 2 0
  ASSERT(g.fw_csr[0x7] == 1); // 2 1
  ASSERT(g.fw_csr[0x8] == 3); // 2 3
  ASSERT(g.fw_csr[0x9] == 1); // 3 1
  ASSERT(g.fw_csr[0xa] == 2); // 3 2
  ASSERT(g.fw_csr[0xb] == 4); // 3 4
  ASSERT(g.fw_csr[0xc] == 1); // 4 1
  ASSERT(g.fw_csr[0xd] == 3); // 4 3
  ASSERT(g.fw_csr[0xe] == 6); // 5 6
  ASSERT(g.fw_csr[0xf] == 5); // 6 5

  ASSERT(g.bw_beg_pos[0] == 0);
  ASSERT(g.bw_beg_pos[1] == g.bw_beg_pos[0] + 2);
  ASSERT(g.bw_beg_pos[2] == g.bw_beg_pos[1] + 4);
  ASSERT(g.bw_beg_pos[3] == g.bw_beg_pos[2] + 3);
  ASSERT(g.bw_beg_pos[4] == g.bw_beg_pos[3] + 3);
  ASSERT(g.bw_beg_pos[5] == g.bw_beg_pos[4] + 2);
  ASSERT(g.bw_beg_pos[6] == g.bw_beg_pos[5] + 1);
  ASSERT(g.bw_beg_pos[7] == g.bw_beg_pos[6] + 1);

  ASSERT(g.bw_csr[0x0] == 1); // 0 1
  ASSERT(g.bw_csr[0x1] == 2); // 0 2
  ASSERT(g.bw_csr[0x2] == 0); // 1 0
  ASSERT(g.bw_csr[0x3] == 2); // 1 2
  ASSERT(g.bw_csr[0x4] == 3); // 1 3
  ASSERT(g.bw_csr[0x5] == 4); // 1 4
  ASSERT(g.bw_csr[0x6] == 0); // 2 0
  ASSERT(g.bw_csr[0x7] == 1); // 2 1
  ASSERT(g.bw_csr[0x8] == 3); // 2 3
  ASSERT(g.bw_csr[0x9] == 1); // 3 1
  ASSERT(g.bw_csr[0xa] == 2); // 3 2
  ASSERT(g.bw_csr[0xb] == 4); // 3 4
  ASSERT(g.bw_csr[0xc] == 1); // 4 1
  ASSERT(g.bw_csr[0xd] == 3); // 4 3
  ASSERT(g.bw_csr[0xe] == 6); // 5 6
  ASSERT(g.bw_csr[0xf] == 5); // 6 5
}
