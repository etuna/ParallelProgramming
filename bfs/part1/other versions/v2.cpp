#include <mpi.h>
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
                        //updatedVals[neighbor] = depth + 1;
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &keep_going, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, result, num_vertices, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        depth++;
    }

    //print_result(graph, result, depth);
    return std::chrono::duration_cast<us>(Time::now() - start_time).count();
}