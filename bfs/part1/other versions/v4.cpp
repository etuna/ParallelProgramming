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

    auto start_time = Time::now();

    int depth = 0;
    int i;
    /*     for (i = 0; i < num_vertices; i++)
    {
        result[i] = -999999;
    } */
    result[start_vertex] = depth;
    int keep_going = true;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int *keep_goings = new int[size];

    int *main_result = new int[num_vertices * size];
    int *main_upd = new int[num_vertices * size];
    int *main_keep_going = new int[size];

    int main_res_size = num_vertices * size * 2 + size;
    int *main_updvals = new int[num_vertices * size];

    // list<int> main_upd;
    int p;
    /*     for (p = 0; p < size * num_vertices; p++)
    {
        main_result[p] = 99999;
    } */
    int **results = new int *[size];
    int nind;
    for (i = 0; i < size; i++)
    {
        results[size] = new int[num_vertices];
    }
    int *updatedVals = new int[num_vertices];
    int *tmp_result = new int[num_vertices * 2 + 1];
    int *tmp_upd = new int[num_vertices];
    int *upd_ref = new int[num_vertices];
    while (keep_going)
    {
        keep_going = false;
        fill_n(main_upd, num_vertices, -1);
        fill_n(upd_ref, num_vertices, -1);
        fill_n(updatedVals, num_vertices, -1);
        int ind = 0;
        int kk;
        for (int vertex = 0; vertex < num_vertices;)
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
                        if(upd_ref[neighbor]==-1){
                        updatedVals[ind] = neighbor;
                        ind++;
                        upd_ref[neighbor]=1;
                        }
                    }
                }
            }
            vertex += size;
        }
        int c;


        MPI_Allgather(&keep_going, 1, MPI_INT, main_keep_going, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(result, num_vertices , MPI_INT, main_result, num_vertices, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(updatedVals, num_vertices, MPI_INT, main_upd, num_vertices, MPI_INT, MPI_COMM_WORLD);

        int km;
        for (c=0; c<size; c++){
            if(main_keep_going[c]){
                keep_going = true;
            }
            for(km=0; km<num_vertices; km++){
                if(main_upd[c*num_vertices+km]==-1){
                    continue;
                }else 
                {
                    int val = main_upd[c*num_vertices+km];
                    if(result[val]>main_result[c*num_vertices+val]){
                        result[val] = main_result[c*num_vertices+val];
                    }
                }
            }
        }
        depth++;
    }
    int h;

    //print_result(graph, result, depth);
    return std::chrono::duration_cast<us>(Time::now() - start_time).count();
}