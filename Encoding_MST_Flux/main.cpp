#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <queue>
#include <utility>
#include <unordered_set>
#include <stack>
#include <list>
#include <stdlib.h>
#include <limits.h>
#include <cstring>

#define INF 1e7
#define NIL -1

using namespace std;

struct Node {
	int degree = 0;
	int data = 0;
	int parent = 0;
};

struct Edge {
	int flow, cap;
	int u, v;
};

struct Vertex {
	int h, e;
};

void DSU(int V, vector<int>& parent, vector<int>& rank) {
	parent.resize(V);
	rank.resize(V);

	for (int i = 0; i < V; i++) {
		parent[i] = -1;
		rank[i] = 1;
	}
}

int findSet(int u, vector<int>& parent) {
	if (parent[u] == -1) {
		return u;
	}

	return parent[u] = findSet(parent[u], parent);
}

void unionSet(int u, int v, vector<int>& parent, vector<int>& rank) {
	int set1 = findSet(u, parent);
	int set2 = findSet(v, parent);

	if (set1 != set2) {
		if (rank[set1] < rank[set2]) {
			parent[set1] = set2;
		}
		else if (rank[set2] < rank[set1]) {
			parent[set2] = set1;
		}
		else {
			parent[set2] = set1;
			rank[set1]++;
		}
	}
}

class Graph {
public:
	string input, output;

	Graph(string i, string o) : input{ i }, output{ o } {}

	// coding algorithms
	void PruferEncoding() {
		vector<Node> tree;
		int size = 0;

		ifstream fin(this->input);
		fin >> size;

		tree.resize(size);

		for (auto i = 0; i < size; i++) {
			int parent;
			fin >> parent;
			tree[i].parent = parent;
			tree[i].data = i;

			if (parent != -1) {
				tree[parent].degree++;
			}
		}

		vector<int> code;

		while (size != 1) {
			int mini = 1e7;

			for (auto i = 0; i < tree.size(); i++) {
				if (tree[i].degree == 0 && tree[i].data != -1 && tree[i].data < mini) {
					mini = tree[i].data;
				}
			}

			tree[mini].data = -1;
			tree[tree[mini].parent].degree--;
			code.push_back(tree[mini].parent);
			size--;
		}

		ofstream fout(this->output);
		fout << code.size() << "\n";
		for (int i = 0; i < code.size(); i++) {
			fout << code[i] << " ";
		}

		fin.close();
		fout.close();
	}
	void PruferDecoding() {
		int size = 0;

		ifstream fin(this->input);
		fin >> size;

		vector<int> tree(size + 1, -1);

		vector<int> code;

		for (auto i = 0; i < size; i++) {
			int data = 0;
			fin >> data;
			code.push_back(data);
		}

		vector<int> exists(size + 1, 0);
		for (auto i = 0; i < size; i++) {
			exists[code[i]] = 1;
		}

		for (auto i = 0; i < size; i++) {
			int current = code.at(0);

			int mini = 1e7;
			
			for (auto j = 0; j < size + 1; j++) {
				if (mini > j && exists[j] == 0) {
					mini = j;
				}
			}
			tree[mini] = current;
			code.erase(code.begin());
			code.push_back(mini);
			
			exists = vector<int>(size + 1, 0);
			for (auto i = 0; i < size; i++) {
				exists[code[i]] = 1;
			}
		}

		ofstream fout(this->output);
		fout << tree.size() << "\n";

		for (auto i = 0; i < tree.size(); i++) {
			fout << tree[i] << " ";
		}

		fin.close();
		fout.close();
	}
	void HuffmannEncoding() {
		struct freqPair {
			char letter = '\0';
			int freq = 0;

			bool operator<(const freqPair& ot) const {
				return freq < ot.freq || (freq == ot.freq && letter < ot.letter);
			}
		};
		 

		ifstream fin(this->input);
		string text;

		getline(fin, text, '\n');

		if (text == "") {
			return;
		}

		vector<int> freq(256, 0);
		int letters = 0;

		for (auto& chr : text) {
			freq[int(chr)]++;
			if (freq[int(chr)] == 1) {
				++letters;
			}
		}

		// litera + codul din secventa
		multimap<freqPair, string> map;
		
		ofstream fout(this->output);
		fout << letters << "\n";

		for (int i = 0; i < 256; i++) {
			if (freq[i]) {
				fout << char(i) << " " << freq[i] << "\n";
			}
		}

		for (int i = 0; i < 256; i++) {
			if (freq[i]) {
				map.insert({ {char(i), freq[i]} , string(1, i) });
			}
		}

		string code[256];
		while (--letters) {
			auto x = *map.begin();

			for (auto& letter : x.second) {
				code[int(letter)] = "0" + code[int(letter)];
			}

			map.erase(map.begin());

			auto y = *map.begin();
			for (auto& letter : y.second) {
				code[int(letter)] = "1" + code[int(letter)];
			}

			map.erase(map.begin());
			map.insert({ {min(x.first.letter, y.first.letter), x.first.freq + y.first.freq}, x.second + y.second });
		}

		for (auto& letter : text) {
			fout << code[int(letter)];
		}

		fin.close();
		fout.close();
	}
	void  HuffmannDecoding() {
		struct freqPair {
			char letter = '\0';
			int freq = 0;

			bool operator<(const freqPair& ot) const {
				return freq < ot.freq || (freq == ot.freq && letter < ot.letter);
			}
		};

		ifstream fin(this->input);

		multimap<freqPair, string> Q;

		vector<int> freq(256, 0);
		int letters = 0;

		string line;
		getline(fin, line);
		letters = stoi(line);

		for (int i = 0; i < letters; i++) {
			char letter;
			getline(fin, line);

			letter = line[0];

			freq[int(letter)] = stoi(line.substr(2));

			Q.insert({ {letter, freq[int(letter)]}, string(1, letter) });
		}

		string encoded;
		fin >> encoded;

		string code[256];

		while (--letters) {
			auto x = *Q.begin();

			for (auto& letter : x.second) {
				code[int(letter)] = "0" + code[int(letter)];
			}

			Q.erase(Q.begin());

			auto y = *Q.begin();
			
			for (auto& letter : y.second) {
				code[int(letter)] = "1" + code[int(letter)];
			}

			Q.erase(Q.begin());

			Q.insert({ {min(x.first.letter, y.first.letter), x.first.freq + y.first.freq}, x.second + y.second });
		}

		map<string, char> codes;
		for (int i = 0; i < 256; i++) {
			if (code[i] != "\0") {
				codes.insert({ code[i], i });
			}
		}

		ofstream fout(this->output);
		string text;
		string currentCode;
		int ind = 0;
		bool found = false;

		while (!encoded.empty()) {
			found = false;
			currentCode += encoded[ind];

			for (auto& code : codes) {
				if (currentCode == code.first) {
					text.append(string(1, code.second));

					encoded.erase(encoded.find(currentCode), currentCode.length());
					currentCode.clear();
					ind = 0;
					found = true;
					break;
				}
			}
			
			if (found == false) {
				ind++;
			}
		}
		fout << text;
		fin.close();
		fout.close();

	}

	// MST algorithms
	void KruskalMST() {
		struct Edge {
			int src = 0;
			int dest = 0;
			int weight = 0;
		};

		ifstream fin(this->input);

		int V, E;
		fin >> V >> E;

		vector<Edge> edges;

		for (int i = 0; i < E; i++) {
			Edge edge;

			fin >> edge.src >> edge.dest >> edge.weight;
			edges.push_back(edge);
		}

		// sorting the edges in increasing order by their weights
		sort(edges.begin(), edges.end(), [](const Edge& e1, const Edge& e2) {
			return e1.weight < e2.weight;
			});

		vector<int> parent;
		vector<int> rank;

		// disjoint set union
		DSU(V, parent, rank);

		vector<Edge> tree;
		int minCost = 0;

		for (auto& edge : edges) {
			int x = edge.src;
			int y = edge.dest;
			int w = edge.weight;

			if (findSet(x, parent) != findSet(y, parent)) {
				unionSet(x, y, parent, rank);
				minCost += w;
				tree.push_back(edge);
			}
		}

		ofstream fout(this->output);
		fout << minCost << "\n";
		fout << tree.size() << "\n";

		sort(tree.begin(), tree.end(), [](const Edge& e1, const Edge& e2) {
			if (e1.src == e2.src) {
				return e1.dest < e2.dest;
				}
			return e1.src < e2.src;
			});

		for (auto& edge : tree) {
			fout << edge.src << " " << edge.dest << "\n";
		}

	}
	void PrimMST() {

		int root = 0;
		ifstream fin(this->input);

		int V, E;
		fin >> V >> E;

		vector<pair<int, int>>* adj = new vector<pair<int, int>>[V];

		for (int i = 0; i < E; i++) {
			int u, v, w;
			fin >> u >> v >> w;
			adj[u].push_back({ v, w });
			adj[v].push_back({ u, w });
		}

		int edges = -1;

		priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;

		vector<int> key(V);
		vector<int> parent(V);
		vector<bool> visited(V);

		for (int i = 0; i < V; i++) {
			key[i] = INF;
			parent[i] = NIL;
			visited[i] = false;
		}

		pq.push({ 0, root });
		key[root] = 0;

		while (!pq.empty()) {
			int current = pq.top().second;
			pq.pop();

			if (visited[current]) {
				continue;
			}

			visited[current] = true;
			edges++;

			for (auto& u : adj[current]) {
				int v = u.first;
				int w = u.second;
				if (key[v] > w && !visited[v]) {
					key[v] = w;
					pq.push({ w, v });
					parent[v] = current;
				}
			}
		}

		int minCost = 0;
		for (int i = 0; i < V; i++) {
			minCost += key[i];
		}

		ofstream fout(this->output);
		fout << minCost << "\n";
		fout << edges << "\n";

		vector<pair<int, int>> tree;

		for (int i = 0; i < V; i++) {
			if (i != root && visited[i] && visited[parent[i]]) {
				tree.push_back({ parent[i], i });
			}
		}

		sort(tree.begin(), tree.end(), [](const pair<int, int>& e1, const pair<int, int>& e2) {
			if (e1.first == e2.first) {
				return e1.second < e2.second;
			}
			return e1.first < e2.first;
		});

		for (auto& e : tree) {
			fout << e.first << " " << e.second << "\n";
		}
	}
};

class Graph_Flux {
private:
	ifstream fin;
	ofstream fout;
	int V, E;

	vector<Vertex> vertices;
	vector<Edge> edges;
	int* current;

	// graf rezidual
	int** rezG;
	int* parent;

	vector<vector<int>> G;
	unordered_set<int>* adj;
	vector<int> circuit;

	// methods to use on flux algorithms
	bool BFS(int** rezG, int S, int T, int* parent) {
		bool visited[this->V];
		memset(visited, false, sizeof(visited));

		queue<int> q;
		q.push(S);

		visited[S] = true;
		parent[S] = -1;

		while (!q.empty()) {
			int u = q.front();
			q.pop();

			for (int v = 0; v < V; v++) {
				if (visited[v] == false && rezG[u][v] > 0) {
					if (v == T) {
						parent[v] = u;
						return true;
					}

					q.push(v);
					parent[v] = u;
					visited[v] = true;
				}
			}
		}
		return false;
	}
	void pushPreflow(int S) {
		vertices[S].h = vertices.size();

		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].u == S) {
				edges[i].flow = edges[i].cap;
				vertices[edges[i].v].e = edges[i].cap;
				vertices[edges[i].u].e -= edges[i].cap;
				edges.push_back({ -edges[i].flow, 0, edges[i].v, S });
			}
		}
	}
	int overflowVertex(vector<Vertex>& ver) {
		for (int i = 1; i < ver.size() - 1; i++) {
			if (ver[i].e > 0) {
				return i;
			}
		}
		return -1;
	}
	void updateReverseEdgeFlow(int i, int flow) {
		int u = edges[i].v;
		int v = edges[i].u;

		for (int j = 0; j < edges.size(); j++) {
			if (edges[j].v == v && edges[j].u == u) {
				edges[j].flow -= flow;
				return;
			}
		}

		Edge e = { 0, flow, u, v };
		edges.push_back(e);
	}
	bool push(int u) {
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].u == u) {
				if (edges[i].flow == edges[i].cap) {
					continue;
				}
				if (vertices[u].h > vertices[edges[i].v].h) {
					int flow = min(edges[i].cap - edges[i].flow, vertices[u].e);
					vertices[u].e -= flow;
					vertices[edges[i].v].e += flow;
					edges[i].flow += flow;
					updateReverseEdgeFlow(i, flow);
					return true;
				}
			}
		}
		return false;
	}
	void relabel(int u) {
		int minHeight = INT_MAX;
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].u == u) {
				if (edges[i].flow == edges[i].cap) {
					continue;
				}
				if (vertices[edges[i].v].h < minHeight) {
					minHeight = vertices[edges[i].v].h;
					vertices[u].h = minHeight + 1;
				}
			}
		}
	}
	void pushFlow(int u, int v) {
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].u == u && edges[i].v == v) {
				int deltaFlow = min(edges[i].cap - edges[i].flow, vertices[u].e);
				vertices[u].e -= deltaFlow;
				vertices[edges[i].v].e += deltaFlow;
				edges[i].flow += deltaFlow;
				updateReverseEdgeFlow(i, deltaFlow);
			}
		}
	}
	void discharge(int u) {
		while (vertices[u].e > 0) {
			if (current[u] == V) {
				relabel(u);
				current[u] = 0;

				continue;
			}

			int v = current[u];
			bool pushed = false;
			for (int i = 0; i < edges.size() && pushed == false; i++) {
				if (edges[i].u == u && edges[i].v == v) {
					if (edges[i].cap - edges[i].flow > 0 && vertices[u].h == vertices[v].h + 1) {
						pushed = true;
						pushFlow(u, v);
					}
				}
			}
			if (pushed == false) {
				current[u]++;
			}
		}
	}

public:
	

	Graph_Flux(string in, string out) : fin{ in }, fout{ out } {
		fin >> this->V >> this->E;

		this->G = vector<vector<int>>(this->V);
		for (int i = 0; i < V; i++) {
			this->G[i] = vector<int>(this->V, 0);
		}

		this->parent = new int[V];

		this->rezG = new int* [V];

		for (int i = 0; i < this->V; i++) {
			this->parent[i] = -1;
			this->rezG[i] = new int[V];
			for (int j = 0; j < this->V; j++) {
				this->rezG[i][j] = 0;
			}
		}

		this->adj = new unordered_set<int>[this->V];

		for (int i = 0; i < this->V; i++) {
			this->vertices.push_back({ 0, 0 });
		}

		for (int i = 0; i < this->E; i++) {
			int u, v, w;
			fin >> u >> v >> w;

			this->adj[u].insert(v);
			this->adj[v].insert(u);
			this->edges.push_back({ 0, w, u, v });

			this->G[u][v] = this->G[v][u] = 1;
			this->rezG[u][v] = w;
		}

		this->current = new int[V];
		this->fin.close();
	}

	void findEulerianCircuit(vector<int>& L, vector<vector<int>>& G, int u) {
		for (int v = 0; v < this->V; v++) {
			if (this->G[u][v] != 0) {
				this->G[u][v] = this->G[v][u] = 0;
				this->findEulerianCircuit(L, G, v);
			}
		}

		L.emplace_back(u);
	}

	void Fleury() {
		vector<int> L;

		findEulerianCircuit(L, G, 0);
		L.pop_back();

		for (auto& u : L) {
			fout << u << " ";
		}

		fout.close();
	}

	void Hierholzer(int s) {
		stack<int> S;
		S.push(s);

		while (!S.empty()) {
			int u = S.top();
			if (adj[u].size() != 0) {
				int v = *adj[u].begin();
				adj[u].erase(v);
				adj[v].erase(u);
				S.push(v);
				continue;
			}

			S.pop();
			circuit.push_back(u);
		}

		circuit.pop_back();

		for (auto& i : circuit) {
			fout << i << " ";
		}

		fout.close();
	}

	// Ford-Fulkerson using Edmond-Karp to find the maximum flow in a network
	int FordFulkerson() {
		int S = 0;
		int T = this->V - 1;

		int maxFlow = 0;

		while (BFS(rezG, S, T, parent)) {
			int pathFlow = INT_MAX;

			for (int v = T; v != S; v = parent[v]) {
				int u = parent[v];
				if (rezG[u][v] < pathFlow) {
					pathFlow = rezG[u][v];
				}
			}

			for (int v = T; v != S; v = parent[v]) {
				int u = parent[v];
				rezG[u][v] -= pathFlow;
				rezG[v][u] += pathFlow;
			}

			maxFlow += pathFlow;
		}

		this->fout << maxFlow;
		return maxFlow;
	}

	// Preflow-Pump
	int PreflowPump() {
		int S = 0;
		int T = this->V - 1;

		pushPreflow(S);
		while (overflowVertex(vertices) != -1) {
			int u = overflowVertex(vertices);
			if (!push(u)) {
				relabel(u);
			}
		}
		int maxFlow = vertices.back().e;
		this->fout << maxFlow;
		return maxFlow;
	}

	// Topologic preflow
	int TopologicPrelow() {
		int S = 0;
		int T = V - 1;

		pushPreflow(S);
		list<int> L;
		L.clear();

		for (int i = 0; i < vertices.size(); i++) {
			if (i != S && i != T) {
				L.push_back(i);
				current[i] = 0;
			}
		}

		auto vertex = L.begin();
		while (vertex != L.end()) {
			int u = *vertex;
			int oldHeight = vertices[u].h;

			discharge(u);

			if (vertices[u].h > oldHeight) {
				L.splice(L.begin(), L, vertex);
			}
			vertex++;
		}
		int maxFlow = vertices[T].e;
		this->fout << maxFlow;
		return maxFlow;
	}
};

int main(int argc, char** argv) {
	Graph graph(argv[1], argv[2]);
	//graph.PruferEncoding();
	//graph.PruferDecoding();
	//graph.HuffmannEncoding();
	//graph.HuffmannDecoding();
	//graph.KruskalMST();
	//graph.PrimMST();

	Graph_Flux graphFlux(argv[1], argv[2]);
	//graphFlux.Fleury();
	//graphFlux.Hierholzer(0);
	//graphFlux.FordFulkerson();
	graphFlux.PreflowPump();
	//graphFlux.TopologicPrelow();

	return 0;
}