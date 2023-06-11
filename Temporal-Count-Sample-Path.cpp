/*--------------------------------------
Author: Omkar Bhalerao
Description: The code implements and analyses 3-path sampling on randomly generated temporal multigraphs
----------------------------------------*/

#include<iostream>
#include<string>
#include<utility>
#include<array>
#include<list>
#include<vector>
#include<set>
#include<iterator>
#include<random>
#include <chrono>
#include <cstdlib>
#include<map>
#include<algorithm>

using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

const double BINOMIAL_PARAMETER = 0.8;
const double UNIFORM_LOWER_BOUND = 1.0;
const double UNIFORM_UPPER_BOUND= 50.0;
const int MAX_MULTIPLICITY = 100;
const long int MAX_NUM_SAMPLES = 1000;

struct CSR_edge{
	long int second_end;
	double timestamp;
};

struct temporal_edge{
	long int first_end, second_end;
	double timestamp;

	bool operator<(const temporal_edge &edge) const{									/*Need this function to facilitate the insertion of keys in map, since C++ map always stores the keys and their corresponding values in sorted order of the keys.*/
		return (first_end < edge.first_end) || ((first_end == edge.first_end) 
			   && (second_end < edge.second_end)) || ((first_end == edge.first_end) 
			   && (second_end == edge.second_end) && (timestamp <= edge.timestamp));
	}
};

struct static_edge{
	long int first_end, second_end;
};

struct three_path{
	temporal_edge left_edge, middle_edge, right_edge;
};

long int *Initialise(long int arr_length, long int *array)
{
	for(long int i = 0; i < arr_length; i++)
		array[i] = 0;

	return array;
}

map<long int, set<long int>> Initialise_set(long int num_vertices, map<long int, set<long int>> input_map)
{

	for(long int i = 0; i < num_vertices; i++)
	{
		set<long int> visited_vertices;
		input_map[i] = visited_vertices;
	}

	return input_map;
}

/*Output: 1.A randomly generated list of temporal edges. Each entry of this list is a triple (u,v,w), where:
		  		- 'u' is the tail of the edge (u,v)
		  		- 'v' is the head of the edge (u,v)
		  		- 'w' is the timestamp associated with this edge
		  2. The number of temporal edges inserted between a given pair of vertices is decided based on a Bernoulli distribution.
		  3. The timestamp assigned to every temporal edge is drawn uniformly from the interval [UNIFORM_LOWER_BOUND, UNIFORM_UPPER_BOUND].*/

list<temporal_edge> generateRandomGraph(long int num_vertices)
{
	int edge_multiplicity;
	long int edge_count = 0;
	double timestamp;
	list<temporal_edge> list_of_temporal_edges;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	binomial_distribution<int> multiplicity_distribution(MAX_MULTIPLICITY, BINOMIAL_PARAMETER);
	uniform_real_distribution<double> time_distribution(UNIFORM_LOWER_BOUND, UNIFORM_UPPER_BOUND);	
	
	for(long int i = 0; i < num_vertices; i++)
	{
		for(long int j = i; j < num_vertices; j++)
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
	cout<<"The number of temporal edges in the randomly generated temporal network is: "<<edge_count<<"\n";

	return list_of_temporal_edges;
}

/*Ouput: CSR representation of the temporal network T
		 Note that for every vertex v. you store the end points the edges directed away from v in the CSR representation.*/

CSR_edge* temporalConvertToCSR(long int num_vertices, list<temporal_edge> list_of_temporal_edges,  long int *temporal_starting_indices)
{
	long int j; 
	list<temporal_edge> :: iterator itr;
	CSR_edge *temporal_graph_CSR = new CSR_edge[list_of_temporal_edges.size()];

	for(long int i = 0; i < list_of_temporal_edges.size(); i++)
		temporal_graph_CSR[i].timestamp = 0.0;

	for(itr = list_of_temporal_edges.begin(); itr != list_of_temporal_edges.end(); ++itr)
	{
		j = temporal_starting_indices[itr->first_end];
		while(temporal_graph_CSR[j].timestamp != 0.0)   								/*This loop finds the next non-empty position in the subarray graph_CSR[start: start + num_temporal_edges_per_vertex[*itr.first_end]]  */
			j++;
		temporal_graph_CSR[j].second_end = itr->second_end;
		temporal_graph_CSR[j].timestamp = itr->timestamp;
		//cout<<"Added "<< temporal_graph_CSR[j].second_end <<" in the adjacency list of "<<itr->first_end<<" at index "<<j<<" with timestamp: "<<temporal_graph_CSR[j].timestamp<<" while original timestamp is: " <<itr->timestamp<<"\n";
	}

	return temporal_graph_CSR;

}

/*Output: CSR representation of the underlying static graph G.*/

long int* staticConvertToCSR(long int num_vertices, long int num_static_edges, list<temporal_edge> list_of_temporal_edges, long int* static_starting_indices)
{
	long int start_index_first_end, start_index_second_end, i ,j, static_edge_count = 0;
	long int *static_graph_CSR = new long int[2*num_static_edges];
	map<long int, set<long int>> visited;
	visited = Initialise_set(num_vertices, visited);
	static_graph_CSR = Initialise(num_vertices, static_graph_CSR);

	for(auto itr = list_of_temporal_edges.begin(); itr != list_of_temporal_edges.end(); ++itr)
	{
		if(((visited[itr->first_end]).find(itr->second_end) != (visited[itr->first_end]).end()) and 
		   ((visited[itr->second_end]).find(itr->first_end) != (visited[itr->second_end]).end()))
			continue;
		else
		{
			(visited[itr->first_end]).insert(itr->second_end);
			(visited[itr->second_end]).insert(itr->first_end);
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
					//cout<<"Added "<<static_graph_CSR[j]<<" at index "<<j<<" in "<<itr->second_end<<"'s adjacency list and added "<<static_graph_CSR[i]<<" at index "<<i<<" in "<<itr->first_end<<"'s adjacency list \n";
					break;
				}
			}
			static_edge_count += 1;
		}
	}

	cout<<"The number of edges in the underlying static, undirected graph is: "<<static_edge_count<<"\n";
	return static_graph_CSR;
}

/*Output: A list of temporal edges directed from first_end to second_end in T.*/

list<CSR_edge> temporal_edges_list(long int first_end, long int second_end,  CSR_edge *temporal_graph_CSR, long int* temporal_starting_indices, long int *num_temporal_edges_per_vertex)
{
	long int start;
	list<CSR_edge> specific_temporal_edges;
	start = temporal_starting_indices[first_end];
	for(long int  i = start; i <= start + num_temporal_edges_per_vertex[first_end]; i++)
	{
		if(temporal_graph_CSR[i].second_end == second_end)
			specific_temporal_edges.push_back(temporal_graph_CSR[i]);
	}

	return specific_temporal_edges;
}

/*Output: A list of temporal edges directed from first_end to second_end with timestamp within the interval [left_timestamp, right_timestamp].*/

list<CSR_edge> temporal_edges_within_interval(long int first_end, long int second_end, CSR_edge *temporal_graph_CSR, long int* temporal_starting_indices, long int *num_temporal_edges_per_vertex, double start_timestamp, double end_timestamp)
{
	CSR_edge edge;
	list<CSR_edge> :: iterator itr;
	list<CSR_edge> candidate_temporal_edges = temporal_edges_list(first_end,second_end,temporal_graph_CSR,temporal_starting_indices,num_temporal_edges_per_vertex);
	list<CSR_edge> temporal_edges_within_delta;

	for(itr = candidate_temporal_edges.begin(); itr != candidate_temporal_edges.end(); ++itr)
	{
		if(itr->timestamp >= start_timestamp and itr->timestamp <= end_timestamp)
			{
				edge.second_end = second_end;
				edge.timestamp = itr->timestamp;
			}
			temporal_edges_within_delta.push_back(edge);
	}

	return temporal_edges_within_delta;
}

/*Output: For a vertex u and timestamp t, list of all the temporal edges with timestamp at least t - delta and at most t, directed into u, except for the ones coming into u from destination_vertex.*/

list<temporal_edge> left_candidate_edges(long int candidate_vertex, long int destination_vertex, double candidate_timestamp, long int num_vertices, long int *static_starting_indices, long int *num_static_edges_per_vertex, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, long int*static_graph_CSR, double delta)
{
	long int neighbour;
	list<temporal_edge> left_edges;

	long int start = static_starting_indices[candidate_vertex];
	for(long int i = start; i < start + num_static_edges_per_vertex[candidate_vertex]; i++)
	{
		neighbour = static_graph_CSR[i];
		if(neighbour != destination_vertex){
			for(long int j = temporal_starting_indices[neighbour]; j < temporal_starting_indices[neighbour] + num_temporal_edges_per_vertex[neighbour]; j++)
			{
				if(temporal_graph_CSR[j].second_end == candidate_vertex and ((candidate_timestamp - delta < temporal_graph_CSR[j].timestamp or candidate_timestamp - delta == temporal_graph_CSR[j].timestamp) and (candidate_timestamp  > temporal_graph_CSR[j].timestamp or candidate_timestamp  == temporal_graph_CSR[j].timestamp)))
				{
					temporal_edge edge;
					edge.first_end = neighbour;
					edge.second_end = candidate_vertex;
					edge.timestamp = temporal_graph_CSR[j].timestamp;
					left_edges.push_back(edge);
				}
			}
		}
	}
		
	return left_edges;
}

/*Output: For a vertex u and timestamp t, list of all the temporal edges with timestamp at least t and at most t + delta, directed away from u, except for the ones going into 'destination_vertex'.*/

list<temporal_edge> right_candidate_edges(long int candidate_vertex, long int destination_vertex, double candidate_timestamp, long int num_vertices, long int *static_starting_indices, long int *num_static_edges_per_vertex, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, long int*static_graph_CSR, double delta)
{
	list<temporal_edge> right_edges;

	long int start = temporal_starting_indices[candidate_vertex];

	for(long int j = start; j < start + num_temporal_edges_per_vertex[candidate_vertex]; j++)
	{
		if(temporal_graph_CSR[j].second_end != destination_vertex and  (candidate_timestamp + delta > temporal_graph_CSR[j].timestamp or candidate_timestamp + delta == temporal_graph_CSR[j].timestamp) and (candidate_timestamp  < temporal_graph_CSR[j].timestamp or candidate_timestamp  == temporal_graph_CSR[j].timestamp))
		{
			temporal_edge edge;
			edge.first_end = candidate_vertex;
			edge.second_end = temporal_graph_CSR[j].second_end;
			edge.timestamp = temporal_graph_CSR[j].timestamp;
			right_edges.push_back(edge);
		}
	}

	return right_edges;
}

/*Output: A map, that for each temporal edge (e2), stores the number of 3-paths e1->e2->e3 centered at the edge e2, which satisfy t(e1) <= t(e2) <= t(e3).*/

pair<map<temporal_edge, long int>, long int> setting_weight_distribution(long int num_vertices, long int *static_starting_indices, long int *num_static_edges_per_vertex, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, long int *static_graph_CSR, double delta)
{
	long int total_weight = 0, current_weight = 0;
	pair<map<temporal_edge, long int>, long int> distribution_pair;
	map<temporal_edge, long int> weight_assignment_per_temporal_edge;
	list<temporal_edge> left_candidates;
	list<temporal_edge> right_candidates;

	for(long int i = 0; i < num_vertices; i++)
	{
		for(long int j = temporal_starting_indices[i]; j < temporal_starting_indices[i] + num_temporal_edges_per_vertex[i]; j++)
		{
			temporal_edge edge;
			edge.first_end = i;
			edge.second_end = temporal_graph_CSR[j].second_end;
			edge.timestamp = temporal_graph_CSR[j].timestamp;
			left_candidates = left_candidate_edges(i, temporal_graph_CSR[j].second_end, temporal_graph_CSR[j].timestamp,  num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);
			right_candidates = right_candidate_edges(temporal_graph_CSR[j].second_end, i, temporal_graph_CSR[j].timestamp,num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);
			current_weight = (left_candidates.size()) * (right_candidates.size());
			weight_assignment_per_temporal_edge[edge] = current_weight;
			total_weight += current_weight;
		}
	}

	distribution_pair.first = weight_assignment_per_temporal_edge;
	distribution_pair.second = total_weight;

	return distribution_pair;
}

/**/

long int  num_induced_chordal_four_cycles(three_path path, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, double delta)
{
	long int count = 0;
	list<CSR_edge> diagonal_edges;
	list<CSR_edge> completing_edges;
	list<CSR_edge>:: iterator diagonal_edges_itr;
	list<CSR_edge>:: iterator completing_edges_itr;

	diagonal_edges = temporal_edges_within_interval((path.right_edge).second_end, (path.middle_edge).first_end, temporal_graph_CSR, temporal_starting_indices, num_temporal_edges_per_vertex, (path.right_edge).timestamp, (path.left_edge).timestamp + delta);
	completing_edges = temporal_edges_within_interval((path.right_edge).second_end, (path.left_edge).first_end, temporal_graph_CSR, temporal_starting_indices, num_temporal_edges_per_vertex, (path.right_edge).timestamp, (path.left_edge).timestamp + delta);

	for(diagonal_edges_itr = diagonal_edges.begin(); diagonal_edges_itr != diagonal_edges.end(); ++diagonal_edges_itr)
	{
		for(completing_edges_itr = completing_edges.begin(); completing_edges_itr != completing_edges.end(); ++completing_edges_itr)
		{
			if(completing_edges_itr->timestamp >= diagonal_edges_itr -> timestamp)
				count += 1;
		}
	}

	return count;
}

long int num_four_cycles(three_path path, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, double delta)
{
	long int count = 0;
	list<CSR_edge> completing_edges;
	list<CSR_edge>::iterator completing_edges_itr;

	completing_edges = temporal_edges_within_interval((path.right_edge).second_end, (path.left_edge).first_end, temporal_graph_CSR, temporal_starting_indices, num_temporal_edges_per_vertex, (path.right_edge).timestamp, (path.left_edge).timestamp + delta);

	count = completing_edges.size();
	return count;

}
long int subgraph_count_estimate(long int num_vertices, long int num_static_edges, long int *static_starting_indices, long int *num_static_edges_per_vertex, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, long int *static_graph_CSR, double delta)
{
	long int total_weight = 0, estimate = 0, chordal_count = 0, count = 0, hits = 0;

	temporal_edge left_edge, right_edge;
	list<temporal_edge> left_edges, right_edges;
	map<temporal_edge, long int> assigned_weights;
	pair<map<temporal_edge, long int>, long int> distribution_pair;

	map<temporal_edge, long int>:: iterator map_itr;
	list<temporal_edge>:: iterator left_list_itr, right_list_itr;

	vector<long int> weight_distribution;
	vector<temporal_edge> temporal_edges;


	three_path path;
	std::random_device rd;
    std::mt19937 gen(rd());

    auto setup_start_time = high_resolution_clock::now();

	distribution_pair = setting_weight_distribution(num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);

	auto setup_end_time = high_resolution_clock::now();
	duration<double,milli> setup_weights = setup_end_time - setup_start_time;

	assigned_weights= distribution_pair.first;
	total_weight = distribution_pair.second;

	for(map_itr = assigned_weights.begin(); map_itr != assigned_weights.end(); ++map_itr)
	{
		temporal_edges.push_back(map_itr->first);
		weight_distribution.push_back(map_itr->second);
	}
	discrete_distribution<long int> weighted_edge_sampler(weight_distribution.begin(), weight_distribution.end());

	auto start_time_to_sample_and_count = high_resolution_clock::now();

	cout<<"SAMPLING 3-PATHS AND COUNTING CHORDAL 4-CYCLES ............ \n";

	for(long int i = 0; i < MAX_NUM_SAMPLES; i++)
	{

		path.middle_edge = temporal_edges[weighted_edge_sampler(gen)];
		left_edges = left_candidate_edges((path.middle_edge).first_end, (path.middle_edge).second_end, (path.middle_edge).timestamp,num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);
		right_edges = right_candidate_edges((path.middle_edge).second_end, (path.middle_edge).first_end, (path.middle_edge).timestamp,num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);
		

		if(left_edges.size() == 0 or right_edges.size() == 0)
		{
			chordal_count = 0;
			continue;
		}
		else
		{
			count += 1;
			uniform_int_distribution<long int> left_edge_distribution(0, left_edges.size()-1);
			uniform_int_distribution<long int> right_edge_distribution(0, right_edges.size()-1);

			left_list_itr = left_edges.begin();
			advance(left_list_itr, left_edge_distribution(gen));
			left_edge = *left_list_itr;

			right_list_itr = right_edges.begin();
			advance(right_list_itr, right_edge_distribution(gen));
			right_edge = *right_list_itr;

			path.left_edge = left_edge;
			path.right_edge = right_edge;

			chordal_count = num_induced_chordal_four_cycles(path, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, delta);
		}		

		if(chordal_count != 0)
			hits += 1;


		estimate += chordal_count;
	}
	auto end_time_to_sample_and_count = high_resolution_clock::now();
	duration<double, milli> time_to_sample_and_count = end_time_to_sample_and_count - start_time_to_sample_and_count;

	cout<<"Time to set up weights: "<<setup_weights.count()<<"\n";
	cout<<"Number of samples: "<<count<<"\n";
	cout<<"Number of hits: "<<hits<<"\n";
	cout<<"Time to sample and count across "<<MAX_NUM_SAMPLES<<" chordal 4-cycles is: "<<time_to_sample_and_count.count()<<"\n";
	
	return estimate;
}

vector<three_path> sample_3_paths(long int num_vertices, long int num_static_edges, long int *static_starting_indices, long int *num_static_edges_per_vertex, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, long int *static_graph_CSR, double delta)
{
	long int total_weight = 0, count = 0;

	temporal_edge left_edge, right_edge;
	list<temporal_edge> left_edges, right_edges;
	map<temporal_edge, long int> assigned_weights;
	pair<map<temporal_edge, long int>, long int> distribution_pair;

	map<temporal_edge, long int>:: iterator map_itr;
	list<temporal_edge>:: iterator left_list_itr, right_list_itr;

	vector<long int> weight_distribution;
	vector<temporal_edge> temporal_edges;
	vector<three_path> sampled_three_paths;


	three_path path;
	std::random_device rd;
    std::mt19937 gen(rd());

    auto setup_start_time = high_resolution_clock::now();
	distribution_pair = setting_weight_distribution(num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);
	auto setup_end_time = high_resolution_clock::now();
	duration<double,milli> setup_weights = setup_end_time - setup_start_time;

	assigned_weights= distribution_pair.first;
	total_weight = distribution_pair.second;

	for(map_itr = assigned_weights.begin(); map_itr != assigned_weights.end(); ++map_itr)
	{
		temporal_edges.push_back(map_itr->first);
		weight_distribution.push_back(map_itr->second);
	}

	discrete_distribution<long int> weighted_edge_sampler(weight_distribution.begin(), weight_distribution.end());

	cout<<"SAMPLING 3-PATHS ......... \n";

	for(long int i = 0; i < MAX_NUM_SAMPLES; i++)
	{
		path.middle_edge = temporal_edges[weighted_edge_sampler(gen)];
		left_edges = left_candidate_edges((path.middle_edge).first_end, (path.middle_edge).second_end, (path.middle_edge).timestamp,num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);
		right_edges = right_candidate_edges((path.middle_edge).second_end, (path.middle_edge).first_end, (path.middle_edge).timestamp,num_vertices, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_graph_CSR, static_graph_CSR, delta);
		

		if(left_edges.size() == 0 or right_edges.size() == 0)
			continue;

		else
		{
			count += 1;
			uniform_int_distribution<long int> left_edge_distribution(0, left_edges.size()-1);
			uniform_int_distribution<long int> right_edge_distribution(0, right_edges.size()-1);

			left_list_itr = left_edges.begin();
			advance(left_list_itr, left_edge_distribution(gen));
			left_edge = *left_list_itr;

			right_list_itr = right_edges.begin();
			advance(right_list_itr, right_edge_distribution(gen));
			right_edge = *right_list_itr;

			path.left_edge = left_edge;
			path.right_edge = right_edge;

			sampled_three_paths.push_back(path);

		}		

	}
	cout<<"Number of samples: "<<count<<"\n";
	cout<<"Time to set up weights: "<<setup_weights.count()<<"\n";
	
	return sampled_three_paths;
}

long int cumulative_subgraph_estimate(vector<three_path> sampled_three_paths, long int *temporal_starting_indices, long int *num_temporal_edges_per_vertex, CSR_edge *temporal_graph_CSR, double delta)
{
	long int count = 0, hits = 0;
	list<CSR_edge> diagonal_edges;
	list<CSR_edge> completing_edges;
	list<CSR_edge>:: iterator diagonal_edges_itr;
	list<CSR_edge>:: iterator completing_edges_itr;

	cout<<"COUNTING CHORDAL 4-CYCLES ............... \n";

	for(long int i = 0; i < sampled_three_paths.size(); i++)
	{
		diagonal_edges = temporal_edges_within_interval((sampled_three_paths[i].right_edge).second_end, (sampled_three_paths[i].middle_edge).first_end, temporal_graph_CSR, temporal_starting_indices, num_temporal_edges_per_vertex, (sampled_three_paths[i].right_edge).timestamp, (sampled_three_paths[i].left_edge).timestamp + delta);
		completing_edges = temporal_edges_within_interval((sampled_three_paths[i].right_edge).second_end, (sampled_three_paths[i].left_edge).first_end, temporal_graph_CSR, temporal_starting_indices, num_temporal_edges_per_vertex, (sampled_three_paths[i].right_edge).timestamp, (sampled_three_paths[i].left_edge).timestamp + delta);

		if(diagonal_edges.size() == 0 or completing_edges.size() == 0)
			continue;
		else
		{
			hits += 1;

			for(diagonal_edges_itr = diagonal_edges.begin(); diagonal_edges_itr != diagonal_edges.end(); ++diagonal_edges_itr)
			{
				for(completing_edges_itr = completing_edges.begin(); completing_edges_itr != completing_edges.end(); ++completing_edges_itr)
				{
					if(completing_edges_itr->timestamp >= diagonal_edges_itr -> timestamp)
						count += 1;
				}
			}
		
		}
	}

	cout<<"Number of hits: "<<hits<<"\n"; /*Number of times the sampled 3-path managed to lie on a chordal 4-cycle.*/
	return count;
}

int main(int argc, char *argv[])
{
	long int num_vertices = atoi(argv[1]);
	double delta = stod(argv[2]);
	long int num_static_edges = 0, estimate = 0, alternate_estimate = 0;	

	map<long int, set<long int>> visited;
	list<temporal_edge> list_of_temporal_edges;	
	CSR_edge* temporal_network_CSR;		
	long int *static_graph_CSR;

	long int *num_temporal_edges_per_vertex = new long int[num_vertices];
	long int *temporal_starting_indices = new long int[num_vertices];
	long int *num_static_edges_per_vertex = new long int[num_vertices];
	long int *static_starting_indices = new long int[num_vertices];

	vector<three_path> sampled_paths;
																		
	/*map<temporal_edge, int> temporal_triangle_index;*/

	cout<<"GENERATING RANDOM GRAPH ...............\n";
	list_of_temporal_edges = generateRandomGraph(num_vertices);

	temporal_network_CSR = new CSR_edge[list_of_temporal_edges.size()];														

	num_temporal_edges_per_vertex = Initialise(num_vertices, num_temporal_edges_per_vertex);				/*For every vertex v, num_temporal_edges_per_vertex[v] stores the number of temporal edges directed away from v in T.*/
	temporal_starting_indices = Initialise(num_vertices, temporal_starting_indices);						/*For every vertex v, temporal_starting_indices[v] stores the starting point of v's adjacency list (which consists of x for every (v,x)) in CSR representation of the temporal network T.*/
	num_static_edges_per_vertex = Initialise(num_vertices, num_static_edges_per_vertex);					/*For every vertex v, num_static_edges_per_vertex[v] stores the number of edges incident on v in the underlying static undirected graph G.*/
	static_starting_indices = Initialise(num_vertices, static_starting_indices);							/*For every vertex v, static_starting_indices[v] stores the starting point of v's adjacency list in the in CSR representation of G.*/
	visited = Initialise_set(num_vertices, visited);

	auto start_preprocessing_time = high_resolution_clock::now();

	for(auto itr = list_of_temporal_edges.begin(); itr != list_of_temporal_edges.end(); ++itr) 		
	{
		num_temporal_edges_per_vertex[itr->first_end] += 1;

		if(((visited[itr->first_end]).find(itr->second_end) != (visited[itr->first_end]).end()) and 
		   ((visited[itr->second_end]).find(itr->first_end) != (visited[itr->second_end]).end()))
			continue;
		else
		{
			(visited[itr->first_end]).insert(itr->second_end);
			(visited[itr->second_end]).insert(itr->first_end);
			num_static_edges_per_vertex[itr->first_end] += 1;
			num_static_edges_per_vertex[itr->second_end] += 1;
			num_static_edges += 1;
		}

	}

	for(long int i = 1; i < num_vertices; i++)
	{
		temporal_starting_indices[i] = temporal_starting_indices[i - 1] + num_temporal_edges_per_vertex[i - 1]; 	
		static_starting_indices[i] = static_starting_indices[i-1] + num_static_edges_per_vertex[i-1];
	}

	temporal_network_CSR = temporalConvertToCSR(num_vertices, list_of_temporal_edges, temporal_starting_indices);
	static_graph_CSR = staticConvertToCSR(num_vertices, num_static_edges,list_of_temporal_edges,  static_starting_indices);	
	auto end_preprocessing_time = high_resolution_clock::now();;
	duration<double, milli> preprocessing_time = end_preprocessing_time - start_preprocessing_time;
	cout<<"Pre processing time: "<<preprocessing_time.count()<<"\n";

	cout<<"--------- METHOD-1------------ \n";

	auto start_sample_and_count_one = high_resolution_clock::now();
	estimate = subgraph_count_estimate( num_vertices,  num_static_edges, static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_network_CSR, static_graph_CSR, delta);
	auto end_sample_and_count_one = high_resolution_clock::now();
	duration<double, milli> run_time_method_one = end_sample_and_count_one - start_sample_and_count_one;
	cout<<"Estimate: "<<estimate<<"\n";
	cout<<"Total time to sample and count with method 1: "<<run_time_method_one.count()<<"\n";

	cout<<"---------- METHOD-2 ------------ \n";

	auto start_sampling_time = high_resolution_clock::now();
	sampled_paths = sample_3_paths(num_vertices, num_static_edges,static_starting_indices, num_static_edges_per_vertex, temporal_starting_indices, num_temporal_edges_per_vertex, temporal_network_CSR, static_graph_CSR, delta);
	auto end_sampling_time = high_resolution_clock::now();
	duration<double, milli> sampling_time = end_sampling_time - start_sampling_time;
	cout<<"Time spent in collecting all 3-paths: "<<sampling_time.count()<<"\n";

	auto start_counting_time = high_resolution_clock::now();
	alternate_estimate = cumulative_subgraph_estimate( sampled_paths,temporal_starting_indices, num_temporal_edges_per_vertex,temporal_network_CSR,  delta);
	auto end_counting_time = high_resolution_clock::now();
	duration<double, milli> counting_time = end_counting_time - start_counting_time;

	cout<<"Estimate: "<<alternate_estimate<<"\n";
	cout<<"Time spent in counting from  all 3-paths: "<<counting_time.count()<<"\n";
	cout<<"------------------------------ \n";

	return 0;
}
