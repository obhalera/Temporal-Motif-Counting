/*
Author: Omkar Bhalerao

Given a list of edges present in any motif on 4-vertices, the program outputs its canonical labelling 
To do so, the algorithm:
1. Extracts the adjacency matrix A from the list of edges specified by the user
2. For every permutation pi of A, construct a 16 bit string by concatenating the entries of pi(A)
3. Return the lexicographically smallest of such strings.
*/

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
#include<string>
//#include <bits/stdc++>

const int NUM_MOTIF_VERTICES = 4;

struct temporal_edge{
	int source, destination, timestamp;
};

using namespace std;

void display(int **input_matrix, int num_vertices)
{
	cout<<"The matrix is: \n";
	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = 0; j < num_vertices; j++)
		{
			cout<<input_matrix[i][j]<<" ";
		}
		cout<<"\n";
	}
}

/*Returns a matrix, that contains entries of input_matrix, but with rows labelled row_id_one and row_id_two swapped.*/

int ** permute_row_adjacency_matrix(int **input_matrix, int num_vertices, int row_id_one, int row_id_two)
{
	int temp;

	for(int i = 0; i < num_vertices; i++)
	{
		temp = input_matrix[row_id_one][i];
		input_matrix[row_id_one][i] = input_matrix[row_id_two][i];
		input_matrix[row_id_two][i] = temp;
	}

	return input_matrix;
}

int ** permute_column_adjacency_matrix(int **input_matrix, int num_vertices, int col_id_one, int col_id_two)
{
	int temp;

	for(int i = 0; i < num_vertices; i++)
	{
		temp = input_matrix[i][col_id_one];
		input_matrix[i][col_id_one] = input_matrix[i][col_id_two];
		input_matrix[i][col_id_two] = temp;
	}

	return input_matrix;
}

int **pairwise_permute_matrix(int **input_matrix, int num_vertices, int id_one, int id_two)
{
	input_matrix = permute_row_adjacency_matrix(input_matrix, num_vertices, id_one, id_two);
	input_matrix = permute_column_adjacency_matrix(input_matrix, num_vertices, id_one, id_two);
	return input_matrix;
}

string generate_labelling(int **input_matrix,int num_vertices)
{
	string labelling = "";

	for(int i = 0; i < num_vertices; i++)
	{
		for(int j = 0; j < num_vertices; j++)
		{
				labelling = labelling + to_string(input_matrix[i][j]);
		}
	}

	return labelling;
}

int **permute_matrix(int **input_matrix, int *permutation, int num_vertices)
{
	int index;
	int *swap_status = new int[num_vertices];
	int **permuted_matrix = new int*[num_vertices];

	for(int i = 0; i < num_vertices; i++)
	{
		permuted_matrix[i] = new int[num_vertices];
		swap_status[i] = i;
	}

	for(int row = 0; row < num_vertices; row ++)
	{
		for(int col = 0; col < num_vertices; col++)
		{
			permuted_matrix[row][col] = input_matrix[row][col];
		}
	}

	for(int i = 0; i < num_vertices; i++)
	{
		// cout<<"Swap status:";
		// for(int j = 0; j < num_vertices; j++)
		// 	cout<<swap_status[j];
		// cout<<"\n";
		if(swap_status[i] == permutation[i])
			continue;
		else
		{
			auto itr = find(swap_status, swap_status + num_vertices, permutation[i]);
			index = distance(swap_status, itr);
			permuted_matrix = pairwise_permute_matrix(permuted_matrix, num_vertices, i , index);
			swap(swap_status[i], swap_status[index]);
		}
	}

	return permuted_matrix;
}

string generate_motif_canonical_labelling(int **input_matrix, int num_vertices)
{
	string min_string = "", current_string, min_achieving_permutation_string = "";
	int *min_achieving_permutation = new int[num_vertices];
	int *permutation = new int[num_vertices];
	int **permuted_matrix;

	for(int i = 0; i < num_vertices; i++)
	{
		permutation[i] = i;
		min_string = min_string + to_string(1);
	}

	do{
		current_string = "";
		permuted_matrix = permute_matrix(input_matrix, permutation, num_vertices);
		current_string = generate_labelling(permuted_matrix, num_vertices);
		if(strcmp(current_string.c_str(), min_string.c_str()) < 0)
		{
			min_string = current_string;
			for(int i = 0; i < num_vertices; i++)
				min_achieving_permutation[i] = permutation[i];
		}
	}while(next_permutation(permutation, permutation + num_vertices));

	for(int j = 0; j < num_vertices; j++)
	{
		min_achieving_permutation_string = min_achieving_permutation_string + to_string(min_achieving_permutation[j]);
	}

	cout<<"The permutation that achieves minima is: "<<min_achieving_permutation_string<<"\n";

	return min_string;
}


/*Input: List of edges
Output: Adjacency matrix*/

int ** get_adjacency_matrix(list<temporal_edge> list_of_edges, int num_vertices)
{
	int **adjacency_matrix = new int*[num_vertices];

	for(int i = 0; i < num_vertices; i++)
		adjacency_matrix[i] = new int[num_vertices];

	for(int row = 0; row < num_vertices; row++)
	{
		for(int column = 0; column < num_vertices; column ++)
		{
			adjacency_matrix[row][column] = 0;
		}
	}

	for(auto edge_list_itr = list_of_edges.begin(); edge_list_itr != list_of_edges.end(); ++edge_list_itr)
	{
		adjacency_matrix[edge_list_itr->source][edge_list_itr->destination] = 1;
	}

	return adjacency_matrix;
}

/*Compares two canonical labellings and returns true iff the two labellings are same.*/

bool compare_canonical_labellings(string labelling_1, string labelling_2)
{
	return strcmp(labelling_1.c_str(), labelling_2.c_str()) == 0 ? true : false;
}

int main()
{
	int **input_matrix, **path_one_matrix, **path_two_matrix, **permuted_matrix, **clique_adjacency_matrix;
	int permutation[] = {1,0,2,3};
	string canonical_labelling_graph, canonical_labelling_path_one, canonical_labelling_path_two ;
	list<temporal_edge> edge_list, path_list_one, path_list_two;

	input_matrix = new int *[NUM_MOTIF_VERTICES];
	path_one_matrix = new int *[NUM_MOTIF_VERTICES];
	path_two_matrix = new int *[NUM_MOTIF_VERTICES];
	permuted_matrix = new int *[NUM_MOTIF_VERTICES];
	clique_adjacency_matrix = new int *[NUM_MOTIF_VERTICES];

	for(int i = 0; i <NUM_MOTIF_VERTICES; i++)
	{
		input_matrix[i] = new int[NUM_MOTIF_VERTICES];
		path_one_matrix[i] = new int[NUM_MOTIF_VERTICES];
		path_two_matrix[i] = new int[NUM_MOTIF_VERTICES];
		permuted_matrix[i] = new int[NUM_MOTIF_VERTICES];
		clique_adjacency_matrix[i] = new int[NUM_MOTIF_VERTICES];
	}

	/*Motif formed by edges 0->1, 0->2, 1->2, 2->3, 3->0*/

	for(int i = 0; i < NUM_MOTIF_VERTICES - 1; i++)
	{
		for(int j = 0; j < NUM_MOTIF_VERTICES - 1; j++)
		{
			if(j > i)
			{
				temporal_edge edge;
				edge.source = i;
				edge.destination = j;
				edge.timestamp = (i + j);
				edge_list.push_back(edge);
			}
		}
	}

	temporal_edge edge_one, edge_two;
	edge_one.source = NUM_MOTIF_VERTICES - 2;
	edge_one.destination = NUM_MOTIF_VERTICES - 1;
	edge_one.timestamp = 10;
	edge_list.push_back(edge_one);

	edge_two.source = NUM_MOTIF_VERTICES - 1;
	edge_two.destination = 0;
	edge_two.timestamp = 10;
	edge_list.push_back(edge_two);

	/*Path 0->1->2->3*/
	for(int i = 0; i < NUM_MOTIF_VERTICES - 1; i++)
	{
		temporal_edge edge;
		edge.source = i;
		edge.destination = i + 1;
		edge.timestamp = 100;
		path_list_one.push_back(edge);
	}

	/*Path 1->2->3->0*/
	for(int i = 1; i < NUM_MOTIF_VERTICES; i++)
	{
		temporal_edge edge;
		edge.source = i;
		edge.destination = (i  + 1)%NUM_MOTIF_VERTICES;
		edge.timestamp = 100;
		path_list_two.push_back(edge); 
	}

	for(int i = 0; i < NUM_MOTIF_VERTICES; i++)
	{
		for(int j = 0; j <  NUM_MOTIF_VERTICES; j++)
		{
			if(i < j)
				clique_adjacency_matrix[i][j] = 1;
		}
	}

	 input_matrix = get_adjacency_matrix(edge_list, NUM_MOTIF_VERTICES);
	 display(input_matrix, NUM_MOTIF_VERTICES);
	 // cout<<"After permutation 1,0,2,3 \n";
	 // input_matrix = permute_matrix(input_matrix, permutation, NUM_MOTIF_VERTICES);
	 // display(input_matrix, NUM_MOTIF_VERTICES);

	canonical_labelling_graph = generate_motif_canonical_labelling(input_matrix, NUM_MOTIF_VERTICES);
	cout<<"The canonical labelling for the graph with edges 0->1, 0->2, 1->2, 2->3, 3->0 is: "<<canonical_labelling_graph<<"\n";

	path_one_matrix = get_adjacency_matrix(path_list_one, NUM_MOTIF_VERTICES);
	canonical_labelling_path_one = generate_motif_canonical_labelling(path_one_matrix, NUM_MOTIF_VERTICES);
	cout<<"The canonical labelling of path 1(0->1->2->3) is: "<<canonical_labelling_path_one<<"\n";

	path_two_matrix = get_adjacency_matrix(path_list_two, NUM_MOTIF_VERTICES);
	canonical_labelling_path_two = generate_motif_canonical_labelling(path_two_matrix, NUM_MOTIF_VERTICES);
	cout<<"The canonical labelling of path 2(1->2->3->0) is: "<<canonical_labelling_path_two<<"\n";

	cout<<"------After Permutation----------- \n";
	permuted_matrix = pairwise_permute_matrix(input_matrix, NUM_MOTIF_VERTICES, 1,2);
	display(permuted_matrix,NUM_MOTIF_VERTICES);
	canonical_labelling_graph = generate_motif_canonical_labelling(permuted_matrix, NUM_MOTIF_VERTICES);
	cout<<"The canonical labelling for the original graph with rows 1 and 2 swapped: "<<canonical_labelling_graph<<"\n";

	cout<<"----Clique-------------- \n";
	display(clique_adjacency_matrix, NUM_MOTIF_VERTICES);
	canonical_labelling_graph = generate_motif_canonical_labelling(clique_adjacency_matrix, NUM_MOTIF_VERTICES);
	cout<<"The canonical labelling for a directed clique is: "<<canonical_labelling_graph<<"\n";

	return 0;
}