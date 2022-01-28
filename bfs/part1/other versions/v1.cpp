#include <mpi.h>
#include <list>
int mpi_vertex_dist2(graph_t *graph, int start_vertex, int *result)
{
    int num_vertices = graph->num_vertices;
    fill_n(result, num_vertices, MAX_DIST);

    auto start_time = Time::now();

    int depth = 0;
    result[start_vertex] = depth;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int complete = 0;
    int keep_going = true;
    int *tmpResult = new int[graph->num_vertices];
    while (keep_going)
    {
        keep_going = false;

        for (int vertex = 0; vertex < num_vertices;)
        {
            int mVertex = vertex + rank;
            printf("Rank:%d handling %d vertex\n", rank, mVertex);
            if (result[mVertex] == depth)
            {
                for (int n = graph->v_adj_begin[mVertex];
                     n < graph->v_adj_begin[mVertex] + graph->v_adj_length[mVertex];
                     n++)
                {
                    int neighbor = graph->v_adj_list[n];

                    if (result[neighbor] > depth + 1)
                    {
                        printf("hey rank %d!which vertex??? %d\n", rank, mVertex);
                        result[neighbor] = depth + 1;
                        keep_going = true;
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            vertex += size;
            //MPI_Barrier(MPI_COMM_WORLD);
        }

        depth++;
    }

    //print_result(graph, result, depth);
    return std::chrono::duration_cast<us>(Time::now() - start_time).count();
}

struct dist_data
{
    int *result;
    int startVertex;
    int keepGoing;
} dist_data;

int findOwner(int **vertexOwnerList, int vertex, int size, int numVertices)
{

    int i;
    for (i = 0; i < size; i++)
    {
        int *vertices = vertexOwnerList[i];
        int m;
        for (m = 0; m < numVertices; m++)
        {
            if (vertices[m] == vertex)
            {
                return i;
            }
        }
    }
    return -1;
}

int mpi_vertex_dist(graph_t *graph, int start_vertex, int *result)
{
    int num_vertices = graph->num_vertices;
    fill_n(result, num_vertices, MAX_DIST);
    int depth = 0;
    int i;
    result[start_vertex] = depth;
    int keep_going = true;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int M;
    list<int> update_list;
    list<int> prev_update_list;

    auto start_time = Time::now();
    while (keep_going)
    {

        keep_going = false;

        update_list.clear();
        for (int vertex = 0; vertex < num_vertices; vertex += size)
        {
            int mVertex = vertex + rank;

            if (result[mVertex] == depth)
            {
                for (int n = graph->v_adj_begin[mVertex];
                     n < graph->v_adj_begin[mVertex] + graph->v_adj_length[mVertex]; n++)
                {
                    int neighbor = graph->v_adj_list[n];

                    if (result[neighbor] > depth + 1)
                    {
                        result[neighbor] = depth + 1;
                        keep_going = true;
                        update_list.push_back(neighbor);
                    }
                }
            }
        }
        if (prev_update_list.size() > 0)
        {
            keep_going = true;
        }
        MPI_Allreduce(MPI_IN_PLACE, &keep_going, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        int p;
        int N = update_list.size();
        int K = prev_update_list.size();

        printf("rank:%d, size:%d\n", rank, K);

        int *update = new int[(K + N) * 2];
        for (p = 0; p < K; p++)
        {
            update[p] = prev_update_list.front();
            if (rank == 2)
            {
                printf("prev upd list->%d\n", prev_update_list.front());
            }
            prev_update_list.pop_front();
            update[p + K + N] = result[update[p]];
        }
        for (p = K; p < K + N; p++)
        {
            update[p] = update_list.front();
            update_list.pop_front();
            update[p + K + N] = result[update[p]];
        }
        int *upd_buf;
        int total_send = N + K;
        for (int rn = 0; rn < size; rn++)
        {

            if (rank % 2 == 0)
            {
                MPI_Send(&total_send, 1, MPI_INT, (rank + 1) % size, 1, MPI_COMM_WORLD);
                MPI_Send(update, total_send * 2, MPI_INT, (rank + 1) % size, 1, MPI_COMM_WORLD);
                MPI_Recv(&M, 1, MPI_INT, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                upd_buf = new int[M * 2];
                MPI_Recv(upd_buf, M * 2, MPI_INT, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Recv(&M, 1, MPI_INT, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                upd_buf = new int[M * 2];
                MPI_Recv(upd_buf, M * 2, MPI_INT, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&total_send, 1, MPI_INT, (rank + 1) % size, 1, MPI_COMM_WORLD);
                MPI_Send(update, total_send * 2, MPI_INT, (rank + 1) % size, 1, MPI_COMM_WORLD);
            }
            int r;
            for (r = 0; r < M; r++)
            {
                printf("updating.. rank:%d, upd_el:%d, upd_val%d\n ", rank, upd_buf[r], upd_buf[r + M]);
                if (result[upd_buf[r]] > upd_buf[r + M])

                {
                    result[upd_buf[r]] = upd_buf[r + M];
                    prev_update_list.push_back(upd_buf[r]);
                    if (rank == 3)
                    {
                        printf("upd buf->%d\n", upd_buf[r]);
                    }
                }
            }
        }

        /*         MPI_Allreduce(MPI_IN_PLACE, &keep_going, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
 */
        /*MPI_Allreduce(MPI_IN_PLACE, result, num_vertices, MPI_INT, MPI_MIN, MPI_COMM_WORLD); */
        depth++;
    }

    int kkk;
    if (rank == 1)
    {
        for (kkk = 0; kkk < num_vertices; kkk++)
        {
            printf("result of %d : %d\n", kkk, result[kkk]);
        }
    }
    //print_result(graph, result, depth);
    return std::chrono::duration_cast<us>(Time::now() - start_time).count();
}