#Summary for Test 2

##5.1 Greedy Method
The greedy method is applied to **optimization problems** --- that is: problems that involve searching through a set of configurations to find one that minimizes or maximizes a given objective function (score assigned to different configurations). 

It works best for problems that possess the **greedy choice property**. This is the property that a globally optimal configuration can be found by a series of locally optimal choices, starting from a well-defined configuration. 

###5.1.1 Fractional Knapsack Problem
```python
fractionalKnapsackAlgorithm(S,W):
#       input:      set S such that each item i in S has a positive benefit i.b and a positive weight i.w; positive maximum total weight W
#       output:     amount i.x of each item i in S that maximizes the total benefit while not exceeding the maximum total weight W
for i in S:
    i.x = 0
    i.v = i.b/i.w
w = 0 # cumulative weight thus far
while w < W:
    i = minValueIndex(S)
    i.x = min(i.w, W - w)
    w += i.x
```

The Fractional Knapsack algorithm can be implemented in *O(nlogn)* time. To achieve this a heap-based priority queue is used where the key of each item in the heap is its value index. 

###5.1.2 Task Scheduling
```python
taskSchedulingAlgorithm(T):
#       input:      A set T of tasks such that each task i has a start time i.s and a finish time i.f
#       output:     Non-conflicting schedule of the tasks in T using minimum number of machines
M = [] # set of machines for scheduling
while size(T) > 0:
    i = popEarliestStart(T)
    m = M.getNonConflicting(i)
    if m != null:
        m.schedule(i)
    else:
        M.append(m())
```

Let *i* and *j* be two tasks with start and finish times *i.s*, *j.s* and *i.f*, *j.f*. *i* and *j* are said to be **non-conflicting** if *i.f <= j.s* or *j.f <= i.s*. 

This algorithm produces a schedule of the tasks with the minimum number of machines in *O(nlogn)* time.

##6.1 The Graph Abstract Data Type
###Definitions
- Graphs are composed of **vertices** and **edges**. 
- Vertices of a graph *G* are represented by the set *V* and the set of edges of *G*, *E*, are represented by pairs of vertices from *V*. 
- When a pair of vertices that form an edge are ordered, we say the edge is **directed**, otherwise it is **undirected**. 
- A **directed graph** or **digraph** is a graph containing only directed edges. 
- Two vertices are said to be **adjacent** if they are endpoints of the same edge. 
- An edge is said to be **incident** on a vertex if the vertex is one of the edge's endpoints. 
- The **degree** of a vertex is the number of incident edges of that vertex. 
- **Simple** graphs are graphs without **parallel edges** or **self-loops**.
- A **path** in a graph is a sequence of alternating vertices and edges starting at a vertex and ending at a vertex, such that each edge is incident to its predecessor and successor vertex. A path is **simple** if every vertex in the path is distinct.
- A **cycle** in a graph is a path with the same start and end vertex. A cycle is **simple** if every vertex in the cycle is distinct, except the start and end vertex.
- A **subgraph** *H* of *G* is graph whose vertices and edges are subsets of the vertices and edges of G, respectively.
- A **spanning subgraph** of *G* is a graph that contains all the vertices of *G*.
- A graph is **connected** if there is a path between any two vertices in the graph.
- If a graph *G* is not connected, its maximal connected subgraphs are called the **connected components** of *G*.
- A **forest** is a graph without cycles.
- A **tree** is a connected forest (trees in terms of graphs do not have roots).
- A **spanning tree** of a graph is a spanning subgraph that is a tree.

###Theorem 6.6
If *G* is a graph with *m* edges then: *sumDegree(G) == 2m*.

###Theorem 6.7
If *G* is a directed graph with *m* edges then: *sumInDegree(G) == sumOutDegree(G) == m*.

###Theorem 6.8
Let *G* is a simple graph with *n* vertices and *m* edges. If *G* is undirected, then *m <= n(n-1)/2* and if *G* is directed then *m <= n(n-1)*.

###Theorem 6.11
Let *G* be an undirected graph with *n* vertices and *m* edges. Then we have the following: 
- If *G* is connected, them *m >= n - 1*
- If *G* is a tree, then *m == n - 1*
- If *G* is a forest, then *m <= n - 1*

##6.2 Data Structures for Graphs
We will be covering the **edge list**, **adjacency list**, and **adjacency matrix** structures. For a graph with *n* vertices and *m* edges, an edge list or adjacency list representation uses *O(n+m)* space, whereas an adjacency matrix uses *O(n<sup>2</sup>)* space.

###6.2.1 The Edge List Structure
In the edge list structure, objects representing the vertices and edges of the graph are stored in lists. The main feature of the edge list structure is that it provides direct access from the edges to the vertices they are incident upon. This allows constant time algorithms to be defined for certain edge-based methods of the graph ADT. Unfortunately the vertices do not have direct access to their incident edges or adjacent vertices, therefore any methods that require that information will run in *O(m)* time as they must traverse the entire list of edges first.

###6.2.2 The Adjacency List Structure
The adjacency list structure extends the edge list structure by adding extra information that supports direct access to the incident edges of each vertex. This additional information takes the form of an **incidence container** that contains references to the edges incident upon each vertex, partitioned according to their direction or lack thereof. The space used by the incidence container of a vertex *v* is *O(dev(v))*. Methods returning iterators of the incident edges of a vertex *v* can run in *O(deg(v))* time. 

###6.2.3 The Adjacency Matrix Structure
The adjacency matrix structure also extends the edge list structure by the addition of a matrix that allows the determination of adjacency between pairs of vertices in constant time. The performance achievement is traded off by an increase in the space usage, however, which is now *O(n<sup>2</sup>)*, and in the running time of other methods. For example: the *adjacentVertices* methods needs to go through an entire row or column of the matrix, which takes *O(n)* time.

The adjacency list structure is superior to the adjacency matrix in space, and is superior in time for all methods except for the *areAdjacent* method.

##6.3 Graph Traversal
###6.3.1 Depth-First Search (DFS)
```python
DFS(G,v):
#       input:      a graph G and a vertex v of G
#       output:     a labeling of the edges in the connected component of v as discovery edges and back edges
v.explored = true
for e in G.incidentEdges(v):
    if not e.explored:
        e.explored = true
        w = G.oppositeVertex(v,e)
        if not w.explored:
            e.label = "discovery"
            DFS(G,w)
        else:
            e.label = "back"
```
####Theorem 6.12
Let *G* be an undirected graph on which a DFS traversal starting at a vertex *s* has been performed. Then the traversal visits all the vertices in the connected component of *s*, and the discovery edges form a spanning tree of the connected component of *s*.
####Theorem 6.13
Let *G* be a graph with *n* vertices and *m* edges represented with the adjacency list structure. A DFS traversal of *G* can be performed in *O(n+m)* time. Also, there exist *O(n+m)*-time algorithms based on DFS for the following problems:
- Testing whether *G* is connected
- Computing a spanning forest of *G*
- Computing a path between two vertices of *G*, or reporting that no such path exists
- Computing a cycle in *G*, or reporting that *G* has no cycles

###6.3.3 Breadth-First Search (BFS)
```python
BFS(G,s):
#       input:      a graph G and a vertex s of G
#       output:     a labeling of the edges in the connected component of s as discovery edges and cross edges
L = [] # list of containers
L.append([s]) # append s into an empty container in the list of containers
i = 0
while len(L[i]) > 0:
    L.append([]) # add an empty container
    for v in L[i]:
        for e in G.incidentEdges(v):
            if not e.explored:
                e.explored = true
                w = G.oppositeVertex(v,e)
                if not w.explored:
                    e.label = "discovery"
                    L[i+1].append(w)
                else
                    e.label = "cross"
    i++
```
####Theorem 6.18
Let *G* be an undirected graph on which a BFS traversal starting at vertex *s* has been performed. Then:
- The traversal visits all the vertices in the connected component of *s*
- The discovery edges form a spanning tree *T* of the connected component of *s*
- For each vertex *v* at level *i*, the path of tree *T* between *s* and *v* has *i* edges, and any other path of *G* between *s* and *v* has at least *i* edges
- If *(u,v)* is a cross edge, then the level numbers of *u* and *v* differ by at most 1

####Theorem 6.19
Let *G* be a graph with *n* vertices and *m* edges represented with the adjacency list structure. A BFS traversal of *G* takes *O(n+m)* time. Also, there exist *O(n+m)*-time algorithms based on BFS for the following problems:
- Testing whether *G* is connected
- Computing a spanning forest of *G*
- Computing the connected components of *G*
- Given a start vertex *s* of *G*, computing, for every vertex *v* of *G*, a path with the minimum number of edges between *s* and *v*, or reporting that no such path exists
- Computing a cycle in *G*, or reporting that *G* has no cycles

###Comparing DFS and BFS
BFS can do everything that DFS can do. BFS is better for finding shortest paths in a graph. It produces a spanning tree such that all the non-tree edges are cross edges. The DFS is better for answering complex connectivity questions, such as determining if every pair of vertices in a graph can be connected by two disjoint paths. It produces a spanning tree where all the non-tree edges are back edges. **These properties only hold for undirected graphs.**

##6.4 Directed Graphs/Digraphs
Given a directed graph *G* with arbitrary vertices *u* and *v*:
- We say that *v* is **reachable** from *u* if *G* has a directed path from *u* to *v*.
- The vertex *v* **reaches the edge** *(w,z)* if *v* reaches the origin vertex *w* of the edge.
- *G* is **strongly connected** if *u* reaches *v* and *v* reaches *u*.
- A **directed cycle** of *G* is a cycle where all the edges are traversed according to their directions.
- G is **acyclic** if it has no directed cycles.
- The **transitive closure** of *G* is the digraph *C* such that the vertices of *C* are the same as the vertices of *G*, and *C* has an edge *(u,v)*, whenever *C* has a directed path from *u* to *v*.
###6.4.1 Traversing a Digraph
###6.4.2 Transitive Closure
```python 
FloydWarshall(G):
#       input:      a digraph G with n vertices
#       output:     the transitive closure C of G
```
###6.4.3 DFS and Garbage Collection (Keep it short)
###6.4.4 Directed Acyclic Graphs (DAGs)
```python
TopologicalSort(G):
#       input:      a digraph G with n vertices
#       output:     a topological ordering V of G or an indication that G has a directed cycle
```

##7.1 Single-Source Shortest Paths
##7.2 All-Pairs Shortest Paths
##7.3 Minimum Spanning Trees