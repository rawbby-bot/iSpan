#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <omp.h>
#include <algorithm>
#include <mpi.h>
#include "wtime.h"
#include "graph.h"

#include "scc_common.h"

#include "trim_1_gfq.h"
#include "trim_2_3.h"
#include "color_propagation.h"
#include "fw_bw.h"
#include "openmp_wcc.hpp"
#define INF -1

void scc_detection(
        const graph *g,
        const int alpha,
        const int beta,
        const int gamma,
        const double theta,
        const index_t thread_count,
        double *avg_time,
        int world_rank,
        int world_size,
        int run_time
        )
{
    const index_t vert_count = g->vert_count;
    const long_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;
    if(DEBUG)
        printf("vert_count = %d, edge_count = %ld, avg_degree = %.3lf\n", vert_count, edge_count, avg_degree);

    long_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    long_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    if(VERBOSE)
    {
        for(int i=fw_beg_pos[vert_count]; i<fw_beg_pos[vert_count+1]; ++i)
            printf("%d\n", fw_csr[i]);
    }

    index_t *max_pivot_list = new index_t[thread_count];
    index_t *max_degree_list = new index_t[thread_count];

    index_t *thread_bin = new index_t[thread_count];

    index_t *front_comm = (index_t *)calloc(thread_count, sizeof(index_t));

    long_t *work_comm = (long_t *)calloc(thread_count, sizeof(long_t));

    index_t *size_bw_beg = (index_t *)calloc(thread_count + 1, sizeof(index_t));
    index_t *size_bw_csr = (index_t *)calloc(thread_count + 1, sizeof(index_t));

    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);

    vertex_t wcc_fq_size = 0;

    const index_t pid = world_rank;
    const index_t tid = pid;
    const index_t p_count = world_size;
    index_t s = vert_count / 32;
    if(vert_count % 32 != 0)
        s += 1;
    index_t t = s / p_count;
    if(s % p_count != 0)
        t += 1;
    index_t step = t * 32;
    index_t virtual_count = t * p_count * 32;

    index_t vert_beg = pid * step;
    index_t vert_end = (pid == p_count - 1 ? vert_count : vert_beg + step);
    unsigned int *sa_compress = (unsigned int *)calloc(s, sizeof(unsigned int));

    index_t *small_queue = new index_t[virtual_count + 1];
    index_t *wcc_fq= new index_t[virtual_count + 1];
    vertex_t *vert_map = (vertex_t *)calloc(vert_count + 1, sizeof(vertex_t));
    vertex_t *sub_fw_beg = (vertex_t *)calloc(vert_count + 1, sizeof(vertex_t));
    vertex_t *sub_fw_csr = (vertex_t *)calloc(edge_count + 1, sizeof(vertex_t));
    vertex_t *sub_bw_beg = (vertex_t *)calloc(vert_count + 1, sizeof(vertex_t));
    vertex_t *sub_bw_csr = (vertex_t *)calloc(edge_count + 1, sizeof(vertex_t));

    vertex_t *buf_fw_beg = (vertex_t *)calloc(vert_count / 10 + 1, sizeof(vertex_t));
    vertex_t *buf_fw_csr = (vertex_t *)calloc(vert_count + 1, sizeof(vertex_t));
    vertex_t *buf_bw_beg = (vertex_t *)calloc(vert_count / 10 + 1, sizeof(vertex_t));
    vertex_t *buf_bw_csr = (vertex_t *)calloc(vert_count + 1, sizeof(vertex_t));

    depth_t *fw_sa;
	depth_t *bw_sa;

	if(posix_memalign((void **)&fw_sa,getpagesize(),
		sizeof(depth_t)*(virtual_count + 1)))
		perror("posix_memalign");

    if(posix_memalign((void **)&bw_sa,getpagesize(),
		sizeof(depth_t)*(virtual_count + 1)))
		perror("posix_memalign");
    vertex_t *fq_comm = (vertex_t *)calloc(virtual_count + 1, sizeof(vertex_t));
    vertex_t *scc_id = new vertex_t[virtual_count + 1];
    vertex_t *scc_id_mice = (vertex_t *)calloc(virtual_count + 1, sizeof(vertex_t));
    vertex_t *fw_sa_temp = (vertex_t *)calloc(virtual_count + 1, sizeof(vertex_t));
    vertex_t *color = (vertex_t *)calloc(virtual_count + 1, sizeof(vertex_t));

    for(long_t i=0; i<virtual_count + 1; ++i)
    {

        fw_sa[i] = -1;

        bw_sa[i] = -1;
        scc_id[i] = 0;
        scc_id_mice[i] = -1;

    }

    vertex_t vertex_fw = 0;
    vertex_t vertex_bw = 0;
    index_t size_3_1 = 0;
    index_t size_3_2 = 0;
    bool changed = false;
    double start_time = wtime();
    double end_time;

    {

        double time_size_1_first;
        double time_size_1;
        double time_fw;
        double time_bw;
        double time_size_2;
        double time_size_3;
        double time_gfq;
        double time_color_1;
        double time_color_2;
        double time_color;
        double pivot_time;
        double time_color_init;
        double time_wcc;
        double time_mice_fw_bw;
        const vertex_t upper_bound = vert_count / thread_count * 5;
        vertex_t *thread_queue = new vertex_t[upper_bound];

        double time_comm;

        MPI_Barrier(MPI_COMM_WORLD);
        double time = wtime();

        index_t trim_times = 1;

        trim_1_first(scc_id,
                fw_beg_pos,
                bw_beg_pos,
                vert_beg,
                vert_end);
        std::cout<<tid<<",Computing size_1_first cost,"<<(wtime() - time) * 1000 <<" ms\n";

        double temp_time = wtime();
        MPI_Allgather(MPI_IN_PLACE,
                0,
                MPI_INT,
                scc_id,
                step,
                MPI_INT,
                MPI_COMM_WORLD);

        double time_comm_trim_1 = wtime() - temp_time;
        std::cout<<tid<<",trim-1 comm time,"<<time_comm_trim_1 * 1000 <<",ms\n";

        if(pid == 0)
        {
            time_size_1_first = wtime() - time;
        }

        time = wtime();
        vertex_t root = pivot_selection(scc_id,
                        fw_beg_pos,
                        bw_beg_pos,
                        0,
                        vert_count,
                        fw_csr,
                        bw_csr,
                        max_pivot_list,
                        max_degree_list,
                        pid,
                        thread_count);

        pivot_time = wtime() - time;

        time = wtime();
        fw_bfs(scc_id,
                fw_beg_pos,
                bw_beg_pos,

                vert_beg,
                vert_end,
                fw_csr,
                bw_csr,
                fw_sa,
                front_comm,
                work_comm,
                root,
                tid,
                thread_count,
                alpha,
                beta,
                edge_count,
                vert_count,
                world_size,
                world_rank,
                step,
                fq_comm,
                sa_compress,
                virtual_count);

        time_fw = wtime() - time;
        for(vertex_t i = 0; i < virtual_count / 32; ++i)
        {
            sa_compress[i] = 0;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        time = wtime();
        bw_bfs(scc_id,
                fw_beg_pos,
                bw_beg_pos,
                vert_beg,
                vert_end,
                fw_csr,
                bw_csr,
                fw_sa,
                bw_sa,
                front_comm,
                work_comm,
                root,
                tid,
                thread_count,
                alpha,
                beta,
                edge_count,
                vert_count,
                world_size,
                world_rank,
                step,
                fq_comm,
                sa_compress
                );

        time_bw = wtime() - time;

        printf("bw_bfs done\n");

            time = wtime();

            trim_1_normal(scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr);

            if(tid == 0)
            {
                time_size_1 += wtime() - time;
            }

            time = wtime();
            trim_1_normal(scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr);

            time_size_1 += wtime() - time;

            time = wtime();

        MPI_Allreduce(MPI_IN_PLACE,
                scc_id,
                vert_count,
                MPI_INT,
                MPI_MAX,
                MPI_COMM_WORLD);

            temp_time = wtime();

            gfq_origin(vert_count,
                        scc_id,
                        small_queue,
                        vert_beg,
                        vert_end,
                        fw_beg_pos,
                        fw_csr,
                        bw_beg_pos,
                        bw_csr,
                        sub_fw_beg,
                        sub_fw_csr,
                        sub_bw_beg,
                        sub_bw_csr,
                        front_comm,
                        work_comm,
                        world_rank,
                        world_size,
                        vert_map);

            vertex_t sub_v_count = front_comm[world_rank];
            vertex_t sub_e_count = work_comm[world_rank];

            if(sub_v_count > 0)
            {
                step = sub_v_count / p_count;
                if(sub_v_count % p_count != 0)
                    step += 1;
                vert_beg = pid * step;
                vert_end = (pid == p_count - 1 ? sub_v_count : vert_beg + step);
                for(index_t i = 0; i < sub_v_count; ++i)
                {
                    color[i] = i;
                }

                time = wtime();
                coloring_wcc(
                        color,
                        sub_fw_beg,
                        sub_fw_csr,
                        sub_bw_beg,
                        sub_bw_csr,
                        step,
                        world_size,
                        0,
                        sub_v_count,
                        sub_v_count);

                if(tid == 0)
                {
                    time_wcc += wtime() - time;
                }
                time = wtime();

                process_wcc(0,
                        sub_v_count,
                        wcc_fq,
                        color,
                        wcc_fq_size);

                if(tid == 0)
                {
                    printf("color time (ms), %lf, wcc_fq, %d, time (ms), %lf\n", time_wcc * 1000, wcc_fq_size, 1000 * (wtime() - time));
                }

                mice_fw_bw(color,
                        scc_id_mice,
                        sub_fw_beg,
                        sub_bw_beg,
                        sub_fw_csr,
                        sub_bw_csr,
                        fw_sa_temp,
                        world_rank,
                        world_size,
                        small_queue,
                        sub_v_count,
                        wcc_fq,
                        wcc_fq_size
                        );

                temp_time = wtime();
                MPI_Allreduce(MPI_IN_PLACE,
                        scc_id_mice,
                        vert_count,
                        MPI_LONG,
                        MPI_MAX,
                        MPI_COMM_WORLD);
                time_comm = wtime() - temp_time;

                printf("%d,final comm time,%.3lf\n", tid, time_comm * 1000);

                for(int i = 0; i < sub_v_count; ++i)
                {
                    vertex_t actual_v = small_queue[i];
                    scc_id[actual_v] = small_queue[scc_id_mice[i]];
                }

            }
            if(tid == 0)
            {
                time_mice_fw_bw = wtime() - time;

            }

        if(tid == 0 && run_time != 1)
        {
            avg_time[0] += time_size_1_first + time_size_1 + time_size_2 + time_size_3;

            avg_time[1] += time_fw + time_bw;

            avg_time[2] += time_wcc + time_mice_fw_bw;
            avg_time[4] += time_size_1_first + time_size_1;

            avg_time[5] += time_size_2;
            avg_time[6] += pivot_time;
            avg_time[7] += time_fw;
            avg_time[8] += time_bw;
            avg_time[9] += time_wcc;
            avg_time[10] += time_mice_fw_bw;

            avg_time[13] += time_size_3;
            avg_time[14] += time_gfq;

        }
        if(OUTPUT_TIME)
        {
            if(tid == 0)
            {
                printf("\ntime size_1_first, %.3lf\ntime size_1, %.3lf\ntime pivot, %.3lf\nlargest fw, %.3lf\nlargest bw, %.3lf\nlargest fw/bw, %.3lf\ntrim size_2, %.3lf\ntrim size_3, %.3lf\nwcc time, %.3lf\nmice fw-bw time, %.3lf\nmice scc time, %.3lf\ntotal time, %.3lf\n", time_size_1_first * 1000, time_size_1 * 1000, pivot_time * 1000, time_fw * 1000, time_bw * 1000, (pivot_time + time_fw + time_bw) * 1000, time_size_2 * 1000, time_size_3 * 1000, time_wcc * 1000, time_mice_fw_bw * 1000, (time_wcc + time_mice_fw_bw) * 1000, (time_size_1_first + time_size_1 + pivot_time + time_fw + time_bw + time_size_2 + time_size_3 + time_wcc + time_mice_fw_bw) * 1000);
            }
        }
        #pragma omp barrier
    }
    end_time = wtime() - start_time;
    avg_time[3] += end_time;
    if(DEBUG)
        printf("total time, %.3lf\n", end_time * 1000);

    get_scc_result(scc_id,
            vert_count);

}

