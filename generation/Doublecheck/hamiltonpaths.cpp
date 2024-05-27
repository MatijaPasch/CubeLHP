#include <iostream>
#include <vector>
using namespace std;
 
int hpcnt;

// Data structure to store a graph edge
struct Edge {
    int src, dest;
};
 
// A class to represent a graph object
class Graph
{
public:
    // a vector of vectors to represent an adjacency list
    vector<vector<int>> adjList;
 
    // Constructor
    Graph(vector<Edge> const &edges, int n)
    {
        // resize the vector to hold `n` elements of type `vector<int>`
        adjList.resize(n);
 
        // add edges to the undirected graph
        for (Edge edge: edges)
        {
            int src = edge.src;
            int dest = edge.dest;
 
            adjList[src].push_back(dest);
            adjList[dest].push_back(src);
        }
    }
};
 
// Utility function to print a path
void printPath(vector<int> const &path)
{
    //for (int i: path) {
    //    cout << i << ' ';
    //}
    //cout << endl;
    hpcnt++;
}
 
void hamiltonianPaths(Graph const &graph, int v, vector<bool> &visited,
                    vector<int> &path, int n)
{
    // if all the vertices are visited, then the Hamiltonian path exists
    if (path.size() == n)
    {
        // print the Hamiltonian path
        printPath(path);
        return;
    }
 
    // Check if every edge starting from vertex `v` leads to a solution or not
    for (int w: graph.adjList[v])
    {
        // process only unvisited vertices as the Hamiltonian
        // path visit each vertex exactly once
        if (!visited[w])
        {
            visited[w] = true;
            path.push_back(w);
 
            // check if adding vertex `w` to the path leads to the solution or not
            hamiltonianPaths(graph, w, visited, path, n);
 
            // backtrack
            visited[w] = false;
            path.pop_back();
        }
    }
}
 
void findHamiltonianPaths(Graph const &graph, int n)
{
    // start with every node
    for (int start = 0; start < n; start++)
    {
        // add starting node to the path
        vector<int> path;
        path.push_back(start);
 
        // mark the start node as visited
        vector<bool> visited(n);
        visited[start] = true;
 
        hamiltonianPaths(graph, start, visited, path, n);
    }
}
 
int main()
{
    // hypercube {0,1}^4
    vector<Edge> edges = {
        {0, 1}, {0, 2}, {0, 4}, {0, 8},
        {1, 3}, {1, 5}, {1, 9},
        {2, 3}, {2, 6}, {2, 10},
        {3, 7}, {3, 11},
        {4, 5}, {4, 6}, {4, 12},
        {5, 7}, {5, 13},
        {6, 7}, {6, 14},
        {7, 15},
        {8, 9}, {8, 10}, {8, 12},
        {9, 11}, {9, 13},
        {10, 11}, {10, 14},
        {11, 15},
        {12, 13}, {12, 14},
        {13, 15},
        {14, 15}
    };
 
    // total number of nodes in the graph (labelled from 0 to 3)
    int n = 16;
    hpcnt=0;
 
    // build a graph from the given edges
    Graph graph(edges, n);
 
    findHamiltonianPaths(graph, n);
    cout << hpcnt;
 
    return 0;
}
