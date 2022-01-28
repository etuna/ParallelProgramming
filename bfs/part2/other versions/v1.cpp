int mpi_frontier(graph_t *graph, int start_vertex, int *result)
{
    int num_vertices = graph->num_vertices;
    fill_n(result, num_vertices, MAX_DIST);

    int depth = 0;
    result[start_vertex] = depth;

    int *frontier_in = new int[num_vertices];
    int *frontier_out = new int[num_vertices];
    frontier_in[0] = start_vertex;
    int front_in_size = 1;
    int front_out_size = 0;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int *front_out_sizes = new int[size];
    auto start_time = Time::now();
    while (front_in_size != 0)
    {
        front_out_size = 0;

        for (int v = 0; v < front_in_size; v += size)
        {
            int mV = v + rank;
            int vertex = frontier_in[mV];

            for (int n = graph->v_adj_begin[vertex]; n < graph->v_adj_begin[vertex] + graph->v_adj_length[vertex]; n++)
            {
                int neighbor = graph->v_adj_list[n];

                if (result[neighbor] > depth + 1)
                {
                    result[neighbor] = depth + 1;
                    frontier_out[front_out_size] = neighbor;
                    front_out_size++;
                }
            }
        }

        MPI_Allgather(&front_out_size, 1, MPI_INT, front_out_sizes, 1, MPI_INT, MPI_COMM_WORLD);

        int max = 0;
        int front_out_temp = 0;
        for (int i = 0; i < size; i++)
        {
            if (max > front_out_sizes[i])
            {
                max = front_out_sizes[i];
                front_out_temp += front_out_sizes[i];
            }
        }

        int *update = new int[max * 2];
        fill_n(update, max * 2, -1);
        for (int i = 0; i < front_out_size; i++)
        {
            update[i] = frontier_out[i];
            update[i + max] = result[frontier_out[i]];
        }

        int *front_out_update = new int[size * max * 2];
        MPI_Allgather(update, max * 2, MPI_INT, front_out_update, max * 2, MPI_INT, MPI_COMM_WORLD);
        int front_out_ind = 0;
        for (int i = 0; i < size; i++)
        {
            for (int m = 0; m < max; m++)
            {
                if (front_out_update[m + 2 * max * i] != -1)
                {
                    frontier_out[front_out_ind] = front_out_update[m + 2 * max * i];
                    if (result[front_out_update[m + 2 * max * i]] > front_out_update[m + 2 * max * i + max])
                    {
                        result[front_out_update[m + 2 * max * i]] = front_out_update[m + 2 * max * i + max];
                    }
                    front_out_ind++;
                }
                else
                {
                    continue;
                }
            }
        }
        front_out_size = front_out_temp;

        /* MPI_Allreduce(MPI_IN_PLACE, &front_out_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, result, num_vertices, MPI_INT, MPI_MIN, MPI_COMM_WORLD); */

        front_in_size = front_out_size;
        int *temp = frontier_in;
        frontier_in = frontier_out;
        frontier_out = temp;
        depth++;
    }

    return std::chrono::duration_cast<us>(Time::now() - start_time).count();
}
