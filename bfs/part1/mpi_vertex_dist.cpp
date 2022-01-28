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
    int *emptyArr = new int[2];
    fill_n(emptyArr, 2, MAX_DIST);
    emptyArr[0] = 0;
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
    int *upd_size = new int[size * 2];
/*     auto t_time1 = Time::now(), t_time2 = Time::now(), t_time3 = Time::now(), t_time4 = Time::now();
 */    int *displs = new int[size];
    fill_n(displs, size, 0);
    int *substatus = new int[2];
    auto start_time = Time::now();
    while (keep_going)
    {
        keep_going = false;
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
/*         if (rank == 0)
        {
            t_time1 = Time::now();
        }
 */
        int num_upd = update_list.size();
        substatus[0] = keep_going;
        substatus[1] = num_upd;

        //MPI_Allreduce(MPI_IN_PLACE, &keep_going, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allgather(substatus, 2, MPI_INT, upd_size, 2, MPI_INT, MPI_COMM_WORLD);
/*         if (rank == 0)
        {
            printf("Duration1:%d\n", std::chrono::duration_cast<us>(Time::now() - t_time1).count());
        }

        if (rank == 0)
        {
            t_time2 = Time::now();
        } */
        int max = 0;
        for (int ss = 0; ss < size; ss++)
        {
            if (upd_size[ss * 2])
            {
                keep_going = true;
            } 
            if (upd_size[ss * 2 + 1] > max)
            {
                max = upd_size[ss * 2 + 1];
            }
        }

        int *local_upd = new int[max * 2];
        fill_n(local_upd, max * 2, -1);
        for (int nu = 0; nu < num_upd; nu++)
        {
            local_upd[nu] = update_list.front();
            update_list.pop_front();
            local_upd[nu + max] = result[local_upd[nu]];
        }

        int *upd_list_sync = new int[max * 2 * size];

/*         if (rank == 0)
        {
            printf("Duration2:%d\n", std::chrono::duration_cast<us>(Time::now() - t_time2).count());
        }

        if (rank == 0)
        {
            t_time3 = Time::now();
        } */
        MPI_Allgather(local_upd, max * 2, MPI_INT, upd_list_sync, max * 2, MPI_INT, MPI_COMM_WORLD);

        for (int ns = 0; ns < size; ns++)
        {
            if(ns == rank){
                continue;
            }
            for (int ins = ns * max * 2; ins < ns * max * 2 + max; ins++)
            {


                    if (upd_list_sync[ins] != -1)
                    {
                        if (result[upd_list_sync[ins]] > upd_list_sync[ins + max])
                        {
                            result[upd_list_sync[ins]] = upd_list_sync[ins + max];
                        }
                    }
                
            }
        }
/*         if (rank == 0)
        {
            printf("Duration3:%d\n", std::chrono::duration_cast<us>(Time::now() - t_time3).count());
        } */
        /*         MPI_Allreduce(MPI_IN_PLACE, &keep_going, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
 */
        /*MPI_Allreduce(MPI_IN_PLACE, result, num_vertices, MPI_INT, MPI_MIN, MPI_COMM_WORLD); */

        depth++;
    }

    //print_result(graph, result, depth);
    return std::chrono::duration_cast<us>(Time::now() - start_time).count();
}