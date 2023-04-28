#include<iostream>
#include<random>
#include <chrono>
#include <cstdlib>
using namespace std;

/*Determines the number of edges in the multigraph, which are directed into the vertex labelled 'vertex_id'*/
int incoming_edges(int **multigraph, int vertex_id)
{
	int incoming_edges = 0;
	for(int w = 0; w < num_vertices; w++)
	{
		incoming_edges += multigraph[w][vertex_id];
	}

	return incoming_edges;
}

/*Determines the number of edges in the multigraph, which are directed away from the vertex labelled 'vertex_id'*/
int outgoing_edges(int **multigraph, int vertex_id)
{
	int outgoing_edges = 0;
	for(int w = 0; w < num_vertices; w++)
	{
		outgoing_edges += multigraph[vertex_id][w];
	}

	return outgoing_edges;
}

/*Determines the degree of the vertex labelled 'vertex_id' in the underlying static graph*/
int degree(int **underlying_static_graph, int vertex_id)
{
	int degree = 0;
	for(int w = 0; w < num_vertices; w++)
	{
		degree += underlying_static_graph[w][vertex_id]
	}

	return degree;
}

/*Determines the number of 3-paths of the form 1->2->3->4 in the input multigraph*/
int count_three_paths(int **multigraph, int num_vertices)
{
	int three_path_count = 0;

	for(int u = 0; u < num_vertices, u++)
	{
		for(int v = 0; v < num_vertices; v++)
		{
			three_path_count += (incoming_edges(multigraph,u))*(multigraph[u][v])*(outgoing_edges(multigraph,v));
		}
	}

	return three_path_count;
}
/*Determines the number of 3-paths of the form 1->2->3->4 in the input multigraph G w.r.t to a degree based ordering of the underlying static graph.*/
/*In particular, imagine sorting the vertices of G w.r.t their degree in the underlying static undirected graph. For each edge (u->v) in G, if w is an in-neighbour of u that follows v, x is an outneighbour of of v that follows u in the ordering, and (x->w) forms a directed edge, only then we increment the number of 3-paths. Note that every chordal 4-cycle will contain such a 3-path.*/
int count_centered_three_paths(int **multigraph, int **underlying_static_graph, int num_vertices)
{
	int centered_three_path_count = 0;

	for(u = 0; u < num_vertices; u++)
	{
		for(v = 0; v < num_vertices; v++)
		{
			if(multigraph[u][v] > 0)
			{
				centered_three_path_count += count_centered_three_paths_per_edge(u,v,multigraph,underlying_static_graph,num_vertices);
			}
		}
	}

	return centered_three_path_count;
}

int count_centered_three_paths_per_edge(int tail_id, int head_id, int **multigraph, int **underlying_static_graph, int num_vertices)
{

	centered_three_path_pereedge_count = 0;
	for(w = 0; w < num_vertices; w++)
	{
		for(x = 0; x < num_vertices; x++)
		{
			if(degree(underlying_static_graph, w) > degree(tail_id) and degree(underlying_static_graph, x) > degree(head_id))
			{
				centered_three_path_pereedge_count += (multigraph[w][head_id])*(multigraph[head_id][tail_id])*(multigraph[tail_id][x]);
			}
		}
	}
	return centered_three_path_pereedge_count;
}
/*Determines the number of 3-paths 1->2->3->4 in the input multigraph w.r.t to an ordering on the vertices, resulting from the number of edges incident on every vertex.*/
/*For each edge (u->v) in G, if w is an in-neighbour of u that follows v, x is an outneighbour of of v that follows u in the ordering, and (x->w) forms a directed edge, only then we increment the number of 3-paths. Note that every chordal 4-cycle will contain such a 3-path.*/


int main(int argc, char *argv[]) /*argv[0]: number of vertices, argv[1]: edge multiplicity paramter for binomial distribution*/
{
	int edge_multiplicity, direction;
	int num_vertices = atoi(argv[1]);
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	int **multigraph; /*The entry in the u-th row and v-th column of this matrix will store the number of outgoing edges from u to v in the multigraph.*/
	int **underlying_static_graph; /*The entry in the u-th row and v-th column of this matrix  will be 1 if anf only if the input multigraph contains an edge between u and v, irrespective of its direction.*/

	default_random_engine generator(seed);
	binomial_distribution<int> multiplicity_distribution(atoi(argv[2]), 0.5);
	bernoulli_distribution direction_distribution(0.5);
	

	multigraph = new int*[num_vertices];
	underlying_static_graph = new int*[num_vertices];

	for(int i = 0; i < num_vertices; i++)
	{	
		multigraph[i] = new int[num_vertices];
		underlying_static_graph[i] = new int[num_vertices];
	}


	/*Initialise the input multigraph and its underlying static undirected graph*/

	for(int u = 0; u < num_vertices; u++)
	{
		for(int v = 0; v < num_vertices; v++)
		{
			multigraph[u][v] = 0;
			underlying_static_graph[u][v] = 0;

			if(u != v)
			{
				edge_multiplicity = multiplicity_distribution(generator);
				if(edge_multiplicity > 0)
				{
					underlying_static_graph[u][v] = 1;
					underlying_static_graph[v][u] = 1;
				}
				for(int i = 0; i < edge_multiplicity; i++)
				{
					direction = direction_distribution(generator);
					if(direction == 1) /*Represents outgoing edges*/
					{
						multigraph[u][v] += 1;
					}
					else
					{
						multigraph[v][u] += 1;
					}
				}
			}

		}
	}

	return 0;
}