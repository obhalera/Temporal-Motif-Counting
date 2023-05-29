#include<iostream>
#include<utility>
#include<list>
#include<iterator>
#include<random>
#include <chrono>
#include <cstdlib>
#include<map>
#include<algorithm>
using namespace std;

const double BINOMIAL_PARAMETER = 0.7;
const double UNIFORM_LOWER_BOUND = 1.0;
const double UNIFORM_UPPER_BOUND= 50.0;
const int MAX_MULTIPLICITY = 100;

struct CSR_edge{
	int second_end;
	double timestamp;
};

struct temporal_edge{
	int first_end, second_end;
	double timestamp;

	bool operator<(const temporal_edge &edge) const{									/*Need this function to facilitate the insertion of keys in map, since C++ map always stores the keys and their corresponding values in sorted order of the keys.*/
		return first_end < edge.first_end || (first_end == edge.first_end) 
			   && (second_end < edge.second_end) || (first_end == edge.first_end) 
			   && (second_end == edge.second_end) && (timestamp <= edge.timestamp);
	}
};

struct static_edge{
	int first_end, second_end;
};

int *Initialise(int arr_length, int *array)
{
	for(int i = 0; i < arr_length; i++)
		array[i] = 0;

	return array;
}

/*Output: 1.A randomly generated list of temporal edges. Each entry of this list is a triple (u,v,w), where:
		  		- 'u' is the tail of the edge (u,v)
		  		- 'v' is the head of the edge (u,v)
		  		- 'w' is the timestamp associated with this edge
		  2. The number of temporal edges inserted between a given pair of vertices is decided based on a Bernoulli distribution.
		  3. The timestamp assigned to every temporal edge is drawn uniformly from the interval [UNIFORM_LOWER_BOUND, UNIFORM_UPPER_BOUND].*/

list<temporal_edge> generateRandomGraph(int num_vertices)
{
	int edge_multiplicity,edge_count = 0;
	double timestamp;
	list<temporal_edge> list_of_temporal_edges;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	binomial_distribution<int> multiplicity_distribution(MAX_MULTIPLICITY, BINOMIAL_PARAMETER);
	uniform_real_distribution<double> time_distribution(UNIFORM_LOWER_BOUND, UNIFORM_UPPER_BOUND);	
	
	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = 0; j < num_vertices; j++)
		{
			if(i == j)
				continue;
			else
			{
				edge_multiplicity = multiplicity_distribution(generator);
				for(int k = 0; k < edge_multiplicity; k++)
				{
					temporal_edge edge;
					timestamp = time_distribution(generator);
					edge.first_end = i;
					edge.second_end = j;
					edge.timestamp = timestamp;
					list_of_temporal_edges.push_back(edge);
				}
				edge_count += edge_multiplicity;
			}
		}
	}
	cout<<"The number of edges in the randomly generated temporal network is: "<<edge_count<<"\n";

	return list_of_temporal_edges;
}

/*Ouput: CSR representation of the temporal network T*/

CSR_edge* temporalConvertToCSR(int num_vertices, list<temporal_edge> list_of_temporal_edges,  int *temporal_starting_indices)
{
	int j;
	list<temporal_edge> :: iterator itr;
	CSR_edge *temporal_graph_CSR = new CSR_edge[list_of_temporal_edges.size()];

	for(itr = list_of_temporal_edges.begin(); itr != list_of_temporal_edges.end(); ++itr)
		itr->timestamp = 0.0;

	for(itr = list_of_temporal_edges.begin(); itr != list_of_temporal_edges.end(); ++itr)
	{
		j = temporal_starting_indices[itr->first_end];
		while(temporal_graph_CSR[j].timestamp != 0.0)   								/*This loop finds the next non-empty position in the subarray graph_CSR[start: start + num_temporal_edges_per_vertex[*itr.first_end]]  */
			j++;
		temporal_graph_CSR[j].second_end = itr->second_end;
		temporal_graph_CSR[j].timestamp = itr->timestamp;
	}

	return temporal_graph_CSR;

}

/*Output: CSR representation of the underlying static graph G.*/

int* staticConvertToCSR(int num_vertices, int num_static_edges, list<temporal_edge> list_of_temporal_edges, int* static_starting_indices)
{
	int start_index_first_end, start_index_second_end, i ,j, static_edge_count = 0;
	int *static_graph_CSR = new int[2*num_static_edges];
	int *visited = new int[num_vertices];
	visited = Initialise(num_vertices, visited);
	static_graph_CSR = Initialise(num_vertices, static_graph_CSR);

	for(auto itr = list_of_temporal_edges.begin(); itr != list_of_temporal_edges.end(); ++itr)
	{
		if(visited[itr->first_end] == 1 and visited[itr->second_end] == 1)
			continue;
		else
		{
			visited[itr->first_end] = 1;
			visited[itr->second_end] = 1;
			i = static_starting_indices[itr->first_end];
			j = static_starting_indices[itr->second_end];
			while(1)
			{
				if(static_graph_CSR[i] != 0)
					i += 1;
				if(static_graph_CSR[j] != 0)
					j += 1;
				if(static_graph_CSR[i] == 0 and static_graph_CSR[j] == 0)
				{
					static_graph_CSR[i] = itr->second_end;
					static_graph_CSR[j] = itr->first_end;
					break;
				}
			}
			static_edge_count += 1;
		}
	}

	cout<<"The number of edges in the underlying static, undirected graph is: "<<static_edge_count<<"\n";
	return static_graph_CSR;
}

/*Output: The lower degree endpoint of e = {first_end, second_end} in the underlying static graph G.*/

int static_edge_degree(int first_end, int second_end,  int *num_static_edges_per_vertex)
{
	return num_static_edges_per_vertex[first_end] <= num_static_edges_per_vertex[second_end] ? first_end : second_end;
}

/*Output: A list of temporal edges directed from first_end to second_end in T.*/

list<CSR_edge> temporal_edges_list(int first_end, int second_end,  CSR_edge *temporal_graph_CSR, int* temporal_starting_indices, int *num_temporal_edges_per_vertex)
{
	int start;
	list<CSR_edge> specific_temporal_edges;
	start = temporal_starting_indices[first_end];
	for(int  i = start; i <= start + num_temporal_edges_per_vertex[first_end]; i++)
	{
		if(temporal_graph_CSR[i].second_end == second_end)
			specific_temporal_edges.push_back(temporal_graph_CSR[i]);
	}

	return specific_temporal_edges;
}

/*Output: A list of temporal edges directed from first_end to second_end with timestamp within the interval [timestamp - delta, timestamp + delta].*/

list<CSR_edge> temporal_edges_within_interval(int first_end, int second_end, CSR_edge *temporal_graph_CSR, int* temporal_starting_indices, int *num_temporal_edges_per_vertex, double timestamp, double delta)
{
	CSR_edge edge;
	list<CSR_edge> :: iterator itr;
	list<CSR_edge> candidate_temporal_edges = temporal_edges_list(first_end,second_end,temporal_graph_CSR,temporal_starting_indices,num_temporal_edges_per_vertex);
	list<CSR_edge> temporal_edges_within_delta;

	for(itr = candidate_temporal_edges.begin(); itr != candidate_temporal_edges.end(); ++itr)
	{
		if(itr->timestamp >= timestamp - delta and itr->timestamp <= timestamp + delta)
			{
				edge.second_end = second_end;
				edge.timestamp = itr->timestamp;
			}
			temporal_edges_within_delta.push_back(edge);
	}

	return temporal_edges_within_delta;
}

/*Output: 'true' if second_vertex is a neighbour of the first_vertex.*/

bool is_neighbour(int first_vertex, int second_vertex, int *static_graph_CSR, int *static_starting_indices, int *num_static_edges_per_vertex)
{
	int start;
	bool output = false;

	start = static_starting_indices[first_vertex];
	for(int i = start; i < start + num_static_edges_per_vertex[first_vertex]; i++)
	{
		if(static_graph_CSR[i] == second_vertex)
			output = true;
	}
	return output;
}

/*Output: The labels of the common neighbours of 'first_vertex' and 'second_vertex'.
		  The method looks for neighbours of 'first_vertex' in the adjacency list of 'second_vertex'.
		  Note that we would like the 'first_vertex' to be the one, that determines the degree of the edge {first_vertex, second_vertex} in the underlying static graph. */

list<int> common_neighbours(int first_vertex, int second_vertex, int *static_starting_indices, int *num_static_edges_per_vertex, int *static_graph_CSR)
{
	bool found;
	auto begin_second_neighbours = static_graph_CSR + static_starting_indices[second_vertex];
	auto end_second_neighbours = begin_second_neighbours + num_static_edges_per_vertex[second_vertex];
	list<int> common_neighbours;

	for(int i = static_starting_indices[first_vertex]; i < static_starting_indices[first_vertex] + num_static_edges_per_vertex[first_vertex]; i++)
	{
		found = (find(begin_second_neighbours, end_second_neighbours, static_graph_CSR[i]) != end_second_neighbours);
		if(found)
			common_neighbours.push_back(static_graph_CSR[i]);
	}

	return common_neighbours;
}

/*Output: A map which, for each temporal edge, stores the number of (undirected) delta-temporal triangles incident on that edge.*/

map<temporal_edge, int> num_triangles_within_delta(int num_vertices, int *static_starting_indices, int *num_static_edges_per_vertex, int *temporal_starting_indices, int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, int *static_graph_CSR, double delta)
{
	int first_end, second_end, delta_triangle_count = 0, delta_edge_count_one = 0, delta_edge_count_two = 0;
	map<temporal_edge, int> temporal_triangle_index;
	list<int> shared_neighbours;
	list<CSR_edge> temporal_edge_list;
	
	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = static_starting_indices[i]; j < static_starting_indices[i] + num_static_edges_per_vertex[i]; j++)
		{

			first_end = static_edge_degree(i,j,num_static_edges_per_vertex);   											/*Finds the vertex with smaller degree*/
			second_end = (i == first_end ? j : i);
			shared_neighbours = common_neighbours(first_end, second_end, static_starting_indices, num_static_edges_per_vertex, static_graph_CSR); 		/*Returns the list of common neighbours*/
			temporal_edge_list = temporal_edges_list(i, static_graph_CSR[j], temporal_graph_CSR, temporal_starting_indices, num_temporal_edges_per_vertex);				/*Enumerate the temporal edges on {i,j}*/

			for(auto itr = temporal_edge_list.begin(); itr != temporal_edge_list.end(); ++itr)					/*Iterate over temporal edges on (i,static_graph_CSR[j])*/
			{
				temporal_edge edge;
				edge.first_end = i;
				edge.second_end = j;
				edge.timestamp = itr->timestamp;

				for(auto ngbr_itr = shared_neighbours.begin(); ngbr_itr != shared_neighbours.end(); ++ngbr_itr) 			/*Iterate over the common neighbours of (i,static_graph_CSR[j])*/
				{
					delta_edge_count_one = ((temporal_edges_within_interval(i, *ngbr_itr,temporal_graph_CSR,temporal_starting_indices, num_temporal_edges_per_vertex, itr->timestamp,delta)).size() 
										+  (temporal_edges_within_interval( *ngbr_itr, i,temporal_graph_CSR,temporal_starting_indices, num_temporal_edges_per_vertex, itr->timestamp,delta)).size());		/*Outputs the number of temporal edges on (i, *ngbr_itr), which are within delta-interval of the timestamp itr->timestamp.*/

					delta_edge_count_two = ((temporal_edges_within_interval(j, *ngbr_itr,temporal_graph_CSR,temporal_starting_indices, num_temporal_edges_per_vertex, itr->timestamp,delta)).size() 
										+  (temporal_edges_within_interval( *ngbr_itr, j,temporal_graph_CSR,temporal_starting_indices, num_temporal_edges_per_vertex, itr->timestamp,delta)).size());

					delta_triangle_count = delta_edge_count_one * delta_edge_count_two;
				}
				temporal_triangle_index[edge] = delta_triangle_count;
			}
		}
	}

	return temporal_triangle_index;
}

int main(int argc, char *argv[])
{
	int num_vertices = atoi(argv[1]); 
	int num_static_edges = 0;													
	int *num_temporal_edges_per_vertex = new int[num_vertices];
	int *temporal_starting_indices = new int[num_vertices];
	int *num_static_edges_per_vertex = new int[num_vertices];
	int *static_starting_indices = new int[num_vertices];
	int *visited = new int[num_vertices];
	list<temporal_edge> list_of_temporal_edges;																					
	map<temporal_edge, int> temporal_triangle_index; 

	list_of_temporal_edges = generateRandomGraph(num_vertices);
																		

	num_temporal_edges_per_vertex = Initialise(num_vertices, num_temporal_edges_per_vertex);				/*For every vertex v, num_temporal_edges_per_vertex[v] stores the number of temporal edges directed away from v in T.*/
	temporal_starting_indices = Initialise(num_vertices, temporal_starting_indices);						/*For every vertex v, temporal_starting_indices[v] stores the starting point of v's adjacency list (which consists of x for every (v,x)) in CSR representation of the temporal network T.*/
	num_static_edges_per_vertex = Initialise(num_vertices, num_static_edges_per_vertex);					/*For every vertex v, num_static_edges_per_vertex[v] stores the number of edges incident on v in the underlying static undirected graph G.*/
	static_starting_indices = Initialise(num_vertices, static_starting_indices);							/*For every vertex v, static_starting_indices[v] stores the starting point of v's adjacency list in the in CSR representation of G.*/
	visited = Initialise(num_vertices, visited);

	for(auto itr = list_of_temporal_edges.begin(); itr != list_of_temporal_edges.end(); ++itr) 		
	{
		num_temporal_edges_per_vertex[itr->first_end] += 1;

		if(visited[itr->first_end] == 1 and visited[itr->second_end] == 1)
			continue;
		else
		{
			visited[itr->first_end] = 1;
			visited[itr->second_end] = 1;
			num_static_edges_per_vertex[itr->first_end] += 1;
			num_static_edges_per_vertex[itr->second_end] += 1;
			num_static_edges += 1;
		}

	}

	for(int i = 1; i < num_vertices; i++)
	{
		temporal_starting_indices[i] = temporal_starting_indices[i - 1] + num_temporal_edges_per_vertex[i - 1]; 	
		static_starting_indices[i] = static_starting_indices[i-1] + num_static_edges_per_vertex[i-1];
	}
	

	return 0;
}