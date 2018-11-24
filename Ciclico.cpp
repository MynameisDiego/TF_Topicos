#include "Ciclico.h"



Ciclico::Ciclico(vector<vector<int>> G)
{
	marked = vector<bool>(G.size());
	onStack = vector<bool>(G.size());
	edgeTo = vector<int>(G.size());
	Cycle = stack<int>();


	for (int v = 0; v < G.size(); v++) {
		if (!marked[v] && Cycle.size() == NULL) {
			dfs(G, v);
		}
	}

}


Ciclico::~Ciclico()
{
}

vector<int> Ciclico::obtenerVecinos(vector<vector<int>> G, int v)
{
	
	vector<int> hijos;

	for (int i = 0; i < G.size(); i++) {
		if (G[i][v] != 0) {
			hijos.push_back(i);
		}
	}

	vector<int> nums(hijos.size());

	for (int i = 0; i < hijos.size(); i++) {
		nums[i] = hijos[i];
	}

	
	
	return nums;
}

void Ciclico::dfs(vector<vector<int>> G, int v)
{

	onStack[v] = true;
	marked[v] = true;

	vector<int> vecino = obtenerVecinos(G, v);


	for (int w : vecino) {
		if (Cycle.size() != NULL) {
			return;
		}
		else if (!marked[w]) {
			edgeTo[w] = v;
			dfs(G, w);
		}
		else if (onStack[w]) {
			for (int x = v; x != w; x = edgeTo[x]) {
				Cycle.push(x);
			}

			Cycle.push(w);
			Cycle.push(v);
		}

	}

	onStack[v] = false;

}

bool Ciclico::hasCycle()
{
	if (Cycle.size() == NULL) {
		return false;
	}

	return true;
}
