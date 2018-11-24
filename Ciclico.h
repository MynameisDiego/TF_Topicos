#pragma once
#include<vector>
#include<iostream>
#include<stack>

using namespace std;

class Ciclico
{

private:
	
	vector<bool> marked;
	vector<int> edgeTo;
	vector<bool> onStack;
	stack<int> Cycle;


public:
	Ciclico(vector<vector<int>> G);
	~Ciclico();
	vector<int> obtenerVecinos(vector<vector<int>> G, int v);
	void dfs(vector<vector<int>> G, int v);
	bool hasCycle();

};

