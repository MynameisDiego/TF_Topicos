#include<iostream>
#include<vector>
#include<string>
#include<math.h>
#include<fstream>
#include<sstream>
#include<stdio.h>
#include <math.h> 
#include <stack>
#include "Ciclico.h"
#include <algorithm>    // std::sort


using namespace std;

int TablaConteoMarginal(int var, vector<vector<int>> data, int val) {
	int prob = 0;
	if (data.size() > 0) {
		for (int i = 0; i < data.size(); i++)
		{
			if (data[i][var] == val) {
				prob++;
			}
		}
	}
	return prob;
}
int TablaConteoConjunta(vector<int> vars, vector<vector<int>>data, vector<int>vals) {
	int contador = 0;

	if (data.size() > 0) {
		bool booleano = false;
		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < vars.size(); j++) {
				if (data[i][vars[j]] != vals[j]) {
					booleano = false;
					break;
				}
				else {
					booleano = true;
				}
			}
			if (booleano) {
				contador++;
			}

		}
	}
	return contador;
}
double CalcularProbMarginal(int var, vector<vector<int>> data, int val, vector<int>card, double alpha) {

	double prob = 0.0;
	if (data.size() > 0) {
		for (int i = 0; i < data.size(); i++)
		{
			if (data[i][var] == val)
				prob++;
		}

		prob = (prob + alpha) / (data.size() + (card[var] * alpha));
	}
	return prob;
}

double CalcularProbMarginalEntropia(vector<int> var, vector<vector<int>> data, vector<int> val, vector<int>card, double alpha) {

	double prob = 0.0;
	if (data.size() > 0) {
		double contador = 0;
		bool boleano = false;
		double cardmul = 1.0;
		for (int i = 0; i < data.size(); i++)
		{
			for (int j = 1; j < var.size(); j++) {
				if (data[i][var[j]] != val[j]) {
					boleano = false;
					break;
				}
				else {
					boleano = true;
				}

			}

			if (boleano) {
				contador++;
			}

		}

		for (int i = 1; i < val.size(); i++) {
			cardmul = cardmul * card[var[i]];
		}

		if (contador == 0) {
			contador = (contador + alpha) / (data.size() + (cardmul * alpha));
			return contador;
		}


		prob = (contador + alpha) / (data.size() + (cardmul * alpha));
	}
	return prob;
}


double CalcularProbMarginalFinal(int var, vector<vector<int>> data, int val, vector<int>card, double alpha) {

	double prob = 0.0;
	if (data.size() > 0) {
		for (int i = 0; i < data.size(); i++)
		{
			if (data[i][var] == val)
				prob++;
		}
	}
	return prob;
}
double CalcularProbConjunta(vector<int> vars, vector<vector<int>>data, vector<int>vals, vector<int>card, double alpha, bool conditional, int cTotal) {

	int cardTotal = cTotal;
	double probabilidad = 0.0;

	if (data.size() > 0) {
		double contador = 0;
		bool booleano = false;
		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < vars.size(); j++) {
				if (data[i][vars[j]] != vals[j]) {
					booleano = false;
					break;
				}
				else {
					booleano = true;
				}
			}
			if (booleano) {
				contador++;
			}

		}

		if (conditional && contador == 0) {
			contador = (contador + alpha) / (data.size() + (cardTotal * alpha));
			return contador;
		}

		probabilidad = (contador + alpha) / ((data.size()) + (cardTotal * alpha));
	}
	return probabilidad;


}
double CalcularProbConjuntaFinal(vector<int> vars, vector<vector<int>>data, vector<int>vals) {

	double probabilidad = 0.0;
	double contador = 0;

	if (data.size() > 0) {
		bool booleano = false;
		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < vars.size(); j++) {
				if (data[i][vars[j]] != vals[j]) {
					booleano = false;
					break;
				}
				else {
					booleano = true;
				}
			}
			if (booleano) {
				contador++;
			}

		}
	}
	return contador;


}
vector<double> CalcularProbCondicional(vector<int> vars, vector<vector<int>> data, vector<vector<int>>vals, vector<int> card, int cardTotal, double alpha) {
	int cTotal = cardTotal;
	vector<double> Probabilities(card[vars[0]]);
	vector<int> Aux;
	vector<int> Aux2;
	double Numerator = 0.0;
	double Denominator = 0.0;

	for (int k = 0; k < vars.size() - 1; k++) {
		Aux.push_back(vars[k + 1]);
	}
	for (int i = 0; i < card[vars[0]]; i++) {
		for (int j = 0; j < vals[i].size() - 1; j++) {
			Aux2.push_back(vals[i][j + 1]);
		}
		Numerator = CalcularProbConjunta(vars, data, vals[i], card, alpha, true, cTotal);
		if (Numerator == 0) {
			Probabilities[i] = Numerator;
			continue;
		}
		if (Aux.size() > 1) {
			Denominator = CalcularProbConjunta(Aux, data, Aux2, card, alpha, true, cTotal);
		}
		else {
			Denominator = CalcularProbMarginal(Aux[0], data, Aux2[0], card, alpha);
		}
		Probabilities[i] = Numerator / Denominator;
	}
	//Normalization
	double Sum = 0.0;
	double Total = 1.0;
	int AuxCounter = 0;
	for (int i = 0; i < Probabilities.size(); i++) {
		if (Probabilities[i] == 0) {
			AuxCounter++;
		}
		Sum += Probabilities[i];
	}
	if (Sum != 1.0) {

		Total = Sum;

		for (int i = 0; i < Probabilities.size(); i++) {
			Probabilities[i] = Probabilities[i] / Total;
		}


		/*Total -= Sum;
		Total /= (AuxCounter);
		for (int i = 0; i < Probabilities.size(); i++) {
			if (Probabilities[i] == 0) {
				Probabilities[i] = Total;
			}
		}*/
	}
	/*Sum = 0;
	for (int i = 0; i < Probabilities.size(); i++) {
		Sum += Probabilities[i];
	}*/
	return Probabilities;

}
vector<vector<int>> FactorConjunta(int N, int i) {
	vector<vector<int>> combs;
	int r = 0;
	vector<int>comb(i);
	int index = 0;
	while (r >= 0) {
		if (index <= (N + (r - i))) {
			comb[r] = index;
			if (r == i - 1) {
				vector<int>aux = comb;
				combs.push_back(aux);
				index++;
			}
			else {
				index = comb[r] + 1;
				r++;
			}
		}
		else {
			r--;
			if (r > 0) {
				index = comb[r] + 1;
			}
			else {
				index = comb[0] + 1;
			}
		}
	}
	return combs;
}
vector<double> CalcularProbs(int iFactor, vector<double>probs, vector<int> vars, vector<int>card, vector<vector<int>>data, vector<string>varnames, vector<vector<string>>varvals, bool condicional, int cardTotal, double alpha) {
	int cTotal = cardTotal;
	int N = probs.size();
	int m = vars.size();
	vector<int>vals(m);
	int cardAnt = 1;
	int iVar;
	cout << "\n";
	//cout << "Factor " << iFactor << endl;
	cout << "\n";
	cout << "i" << "\t";
	int w;
	for (w = 0; w < vars.size() - 1; w++) {
		cout << varnames[vars[w]] + "\t";
	}
	cout << varnames[vars[w]] << endl;
	cout << "\n";
	if (!condicional) {
		for (int i = 0; i < N; i++) {
			cardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = vars[j];
				vals[j] = (int)floor(i / cardAnt) % card[iVar];
				cardAnt *= card[iVar];
			}
			double prob = 0.0;
			if (m > 1) {
				prob = CalcularProbConjunta(vars, data, vals, card, alpha, false, cTotal);
			}
			else {
				prob = CalcularProbMarginal(vars[0], data, vals[0], card, alpha);
			}
			probs[i] = prob;
			cout << i << "\t";
			for (w = 0; w < vals.size() - 1; w++) {
				cout << varvals[vars[w]][vals[w]] + "\t";
			}
			cout << varvals[vars[w]][vals[w]] + "\t" << prob;
			cout << "\n";
		}
	}

	else {

		vector<vector<int>>Auxvals(card[vars[0]], vector<int>(m));

		for (int i = 0; i < N; i += card[vars[0]]) {
			for (int k = 0; k < card[vars[0]]; k++) {
				cardAnt = 1;
				for (int j = 0; j < m; j++) {
					iVar = vars[j];
					vals[j] = (int)((int)floor((i + k) / cardAnt) % card[iVar]);
					cardAnt *= card[iVar];
				}
				Auxvals[k] = vals;
			}
			vector<double> ProbFinal;
			ProbFinal = CalcularProbCondicional(vars, data, Auxvals, card, cTotal, alpha);

			double prob;
			for (int p = 0; p < card[vars[0]]; p++) {

				double Total = 0.0;
				double Sum = 0.0;
				for (int h = 0; h < card[vars[0]]; h++) {
					Total = ProbFinal[h];
					Sum = Sum + Total;
				}
				prob = ProbFinal[p] / Sum;
				cout << (i + p) << "\t";
				for (w = 0; w < vals.size() - 1; w++) {
					cout << varvals[vars[w]][Auxvals[p][w]] + "\t";
				}
				cout << varvals[vars[w]][Auxvals[p][w]] + "\t" << prob / Sum;
				cout << "\n";
			}

		}


	}
	return probs;
}
vector<double> CalcularConteo(int iFactor, vector<double>probs, vector<int> vars, vector<int>card, vector<vector<int>>data, vector<string>varnames, vector<vector<string>>varvals, bool condicional) {
	int N = probs.size();
	int m = vars.size();
	vector<int>vals(m);
	int cardAnt = 1;
	int iVar;
	cout << "\n";
	cout << "Factor " << iFactor << endl;
	cout << "\n";
	cout << "i" << "\t";
	int w;
	for (w = 0; w < vars.size() - 1; w++) {
		cout << varnames[vars[w]] + "\t";
	}
	cout << varnames[vars[w]] << endl;
	cout << "\n";
	if (!condicional) {
		for (int i = 0; i < N; i++) {
			cardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = vars[j];
				vals[j] = (int)floor(i / cardAnt) % card[iVar];
				cardAnt *= card[iVar];
			}
			int prob = 0;
			if (m > 1) {
				prob = TablaConteoConjunta(vars, data, vals);
			}
			else {
				prob = TablaConteoMarginal(vars[0], data, vals[0]);
			}
			probs[i] = prob;
			cout << i << "\t";
			for (w = 0; w < vals.size() - 1; w++) {
				cout << varvals[vars[w]][vals[w]] + "\t";
			}
			cout << varvals[vars[w]][vals[w]] + "\t" << prob;
			cout << "\n";
		}
	}
	return probs;
}
vector<int> HipotsisMatrix(int n, int numOfBits) {

	vector<int> binary(numOfBits);
	for (int i = 0; i < numOfBits; i++, n /= 2) {
		switch (n % 2) {
		case 0:
			binary[i] = 0;
			break;
		case 1:
			binary[i] = 1;
			break;
		}
	}
	return binary;
}

vector<int> getPadres(int i, vector<vector<int>> grafo) {
	vector<int> padres;

	for (int k = 0; k < grafo.size(); k++) {
		if (grafo[i][k] == 1) {
			padres.push_back(k);
		}
	}
	return padres;
}

vector<int>unio(int i, vector<int> padres) {
	int N = padres.size() + 1;
	vector<int> arr(N);
	int k = 0;
	arr[k++] = i;
	for (int x : padres) {
		arr[k++] = x;
	}
	return arr;
}

int Factorial(int Number) {
	int Result = 1;
	if (Number < 0) {
		Result = 1;
	}
	else if (Number == 0) {
		Result = 1;
	}
	else {
		for (int i = 1; i <= Number; i++) {
			Result = Result * i;
		}
	}
	return Result;
}

void Marginalization(vector<vector<int>> Data, vector<int> Vars, vector<int> Card, vector<string>VarNames, vector<vector<string>> VarVals, int size, vector<int>marg) {
	int w;
	int CardAnt = 1;
	int CardNew = 1;
	int iVar;
	int iVarn;
	int k = 0;
	double newP = 0.0;
	cout << "\n";
	cout << "Factor ";
	cout << "\n";
	cout << "i" << "\t";
	int m = Vars.size();
	int v = marg.size();
	vector<int>vals(m);
	vector<int>vars1;
	vector<int>vals1;

	for (w = 0; w < Vars.size() - 1; w++) {
		cout << VarNames[Vars[w]] << "\t";
	}
	cout << VarNames[Vars[w]] << endl;
	int FactorSize = 1;
	for (int k = 0; k < Vars.size(); k++) {
		FactorSize *= Card[Vars[k]];
	}
	for (int i = 0; i < FactorSize; i++) {
		CardAnt = 1;
		for (int j = 0; j < m; j++) {
			iVar = Vars[j];
			vals[j] = (int)((int)floor(i / CardAnt) % Card[iVar]);
			CardAnt *= Card[iVar];
		}
		double Prob = 0.0;

		Prob = CalcularProbConjunta(Vars, Data, vals, Card, 1.0, false, CardAnt);

		cout << i << "\t";
		for (w = 0; w < vals.size() - 1; w++) {
			cout << VarVals[Vars[w]][vals[w]] + "\t";
		}
		cout << VarVals[Vars[w]][vals[w]] + "\t" << Prob;
		cout << "\n";
		k++;

		for (int m = size - 1; m < Vars.size(); m++) {
			vars1.push_back(Vars[m]);

		}

		for (int p = size - 1; p < vals.size(); p++) {
			vals1.push_back(vals[p]);

		}
		CardNew = 1;
		for (int s = 0; s < v; s++) {
			iVarn = marg[s];
			CardNew *= Card[iVarn];
		}

		newP += Prob;
		if (k%CardNew == 0) {
			cout << "Marginalizacion de ";
			for (int g = 0; g < marg.size(); g++) {
				cout << VarNames[marg[g]];
			}
			cout << " : " << newP << "\n";
			newP = 0.0;
		}
	}
}



vector<double> CalcularEntropiaCondicional(vector<int> vars, vector<vector<int>> data, vector<vector<int>>vals, vector<int> card, int cardTotal, int alpha) {
	int cTotal = cardTotal;
	vector<double> EntropiaC(card[vars[0]]);
	vector<int> Aux;
	vector<int> Aux2;
	double Numerator = 0.0;
	double Denominator = 0.0;
	int conteo = 0;

	for (int k = 0; k < vars.size(); k++) {
		Aux.push_back(vars[k]);
	}
	for (int i = 0; i < card[vars[0]]; i++) {
		for (int j = 0; j < vals[i].size(); j++) {
			Aux2.push_back(vals[i][j]);
		}
		//Numerator = CalcularProbConjunta(vars, data, vals[i], card, 1.0, true, cTotal);

		if (Aux.size() > 1) {
			Numerator = CalcularProbConjunta(vars, data, vals[i], card, alpha, true, cTotal);
			Denominator = CalcularProbMarginalEntropia(Aux, data, Aux2, card, alpha);
			//CalcularProbConjunta(Aux, data, Aux2, card, 1.0, true, cTotal);
			EntropiaC[i] = Numerator / Denominator;
			Aux2.clear();

		}
		else {
			Denominator = CalcularProbMarginal(Aux[0], data, Aux2[0], card, alpha);
			EntropiaC[i] = Denominator;
			Aux2.clear();
		}


	}


	double total = 0.0;
	double suma = 0.0;

	for (int h = 0; h < card[vars[0]]; h++) {
		total = EntropiaC[h];
		suma += total;
	}

	for (int i = 0; i < card[vars[0]]; i++) {
		EntropiaC[i] = EntropiaC[i] / suma;
	}


	for (int i = 0; i < EntropiaC.size(); i++) {

		conteo = CalcularProbConjuntaFinal(vars, data, vals[i]);
		EntropiaC[i] = conteo * log10(EntropiaC[i]);

	}

	return EntropiaC;

}

vector<double> CalcularProbabilidadesF2(vector<int> vars, vector<vector<int>> data, vector<vector<int>>vals, vector<int> card, int cardTotal, int alpha) {
	int cTotal = cardTotal;
	vector<double> EntropiaC(card[vars[0]]);
	vector<int> Aux;
	vector<int> Aux2;
	double Numerator = 0.0;
	double Denominator = 0.0;
	int conteo = 0;

	for (int k = 0; k < vars.size(); k++) {
		Aux.push_back(vars[k]);
	}
	for (int i = 0; i < card[vars[0]]; i++) {
		for (int j = 0; j < vals[i].size(); j++) {
			Aux2.push_back(vals[i][j]);
		}
		//Numerator = CalcularProbConjunta(vars, data, vals[i], card, 1.0, true, cTotal);

		if (Aux.size() > 1) {
			Numerator = CalcularProbConjunta(vars, data, vals[i], card, alpha, true, cTotal);
			Denominator = CalcularProbMarginalEntropia(Aux, data, Aux2, card, alpha);
			//CalcularProbConjunta(Aux, data, Aux2, card, 1.0, true, cTotal);
			EntropiaC[i] = Numerator / Denominator;
			Aux2.clear();

		}
		else {
			Denominator = CalcularProbMarginal(Aux[0], data, Aux2[0], card, alpha);
			EntropiaC[i] = Denominator;
			Aux2.clear();
		}


	}


	double total = 0.0;
	double suma = 0.0;

	for (int h = 0; h < card[vars[0]]; h++) {
		total = EntropiaC[h];
		suma += total;
	}

	for (int i = 0; i < card[vars[0]]; i++) {
		EntropiaC[i] = EntropiaC[i] / suma;
	}


	/*for (int i = 0; i < EntropiaC.size(); i++) {

		conteo = CalcularProbConjuntaFinal(vars, data, vals[i]);
		EntropiaC[i] = conteo * log10(EntropiaC[i]);

	}*/

	return EntropiaC;

}





vector<double> calcularProbabilidadesF(vector<double> probs, vector<int> vars, vector<int> card, vector<vector<int>> data, int alpha) {

	int N = probs.size();
	int m = vars.size();
	vector<int> vals(m);
	int cardAnt = 1;
	int iVar;

	vector<vector<int>> auxvals(card[vars[0]], vector<int>(m));

	for (int i = 0; i < N; i += card[vars[0]]) {
		for (int k = 0; k < card[vars[0]]; k++) {
			cardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = vars[j];
				vals[j] = (int)floor((i + k) / cardAnt) % card[iVar];
				cardAnt *= card[iVar];
			}

			auxvals[k] = vals;
		}

		vector<double> prob;
		prob = CalcularProbabilidadesF2(vars, data, auxvals, card, cardAnt, alpha);

		for (int l = 0; l < card[vars[0]]; l++) {
			probs[l + i] = prob[l];
		}

	}

	return probs;
}





vector<double> calcularEntropia(vector<double> probs, vector<int> vars, vector<int> card, vector<vector<int>> data, int alpha) {

	int N = probs.size();
	int m = vars.size();
	vector<int> vals(m);
	int cardAnt = 1;
	int iVar;

	vector<vector<int>> auxvals(card[vars[0]], vector<int>(m));

	for (int i = 0; i < N; i += card[vars[0]]) {
		for (int k = 0; k < card[vars[0]]; k++) {
			cardAnt = 1;
			for (int j = 0; j < m; j++) {
				iVar = vars[j];
				vals[j] = (int)floor((i + k) / cardAnt) % card[iVar];
				cardAnt *= card[iVar];
			}

			auxvals[k] = vals;
		}

		vector<double> prob;
		prob = CalcularEntropiaCondicional(vars, data, auxvals, card, cardAnt, alpha);

		for (int l = 0; l < card[vars[0]]; l++) {
			probs[l + i] = prob[l];
		}

	}

	return probs;
}

double GetK2(vector<vector<int>>G, vector<vector<int>>Data, vector<int>Card, int N, double ndga) {
	double K2 = 0.0;
	int CardAnt = 1;
	int iVar;
	vector<double> kv;
	for (int i = 0; i < N; i++) {
		double x = 0.0;
		vector<int> Padres = getPadres(i, G);
		vector<int> Vars = unio(i, Padres);
		int FactorSize = 1;
		for (int v = 0; v < Vars.size(); v++) {
			FactorSize *= Card[Vars[v]];
		}
		vector<int> Vals(Vars.size());

		vector<double> Prob;

		vector<vector<int>>AuxVals(Card[Vars[0]], vector<int>(Vars.size()));
		for (int i = 0; i < FactorSize; i += Card[Vars[0]]) {
			for (int j = 0; j < Card[Vars[0]]; j++) {
				CardAnt = 1;
				for (int k = 0; k < Vars.size(); k++) {
					iVar = Vars[k];
					Vals[k] = (int)((int)floor((i + j) / CardAnt) % Card[iVar]);
					CardAnt *= Card[iVar];
				}
				AuxVals[j] = Vals;
			}
			for (int m = 0; m < AuxVals.size(); m++) {
				int t = CalcularProbConjuntaFinal(Vars, Data, AuxVals[m]);
				double lg = 0.0;
				for (int w = 1; w <= t; w++) {
					lg += log2(w);
				}
				Prob.push_back(lg);
			}
			if (Padres.size() >= 1) {
				vector<vector<int>> valsA;
				for (int f = 0; f < AuxVals.size(); f++) {
					vector<int>temp;
					for (int g = 0; g < AuxVals[f].size(); g++) {
						if (g != 0) {
							temp.push_back(AuxVals[f][g]);
						}
					}
					valsA.push_back(temp);
				}

				double CountP = CalcularProbConjuntaFinal(Padres, Data, valsA[0]);
				double CardVar = Card[Vars[0]];
				double Rest = CardVar - 1;
				double R = 0.0;
				for (int q = 1; q <= CountP; q++) {
					R += log2(1 / (Rest + q));
				}
				Prob.push_back(R);

			}
			else {
				int counP = 0;
				int a = Card[Vars[0]];
				double R = log2(((Factorial(a - 1)) / (Factorial(a - 1 + counP))));
				Prob.push_back(R);
			}
		}

		for (int t = 0; t < Prob.size(); t++) {
			x += Prob[t];
		}
		kv.push_back(x);
	}

	for (int y = 0; y < kv.size(); y++) {
		K2 += kv[y];
	}
	double pDga = 1 / ndga;
	return K2 + log2(pDga);
}

bool myfunction(double i, double j) {
	return (i < j);
}

struct myclass {
	bool operator() (double i, double j) { return (i < j); }
} myobject;


bool sortcol(const vector<double>& v1,
	const vector<double>& v2) {
	return v1[0] > v2[0];
}

int main()
{

	vector<vector<int>> data;
	vector<vector<int>> grafoCLALL;
	vector<int> v;
	vector<int> card;
	vector<string> varnames;
	vector<vector<string>> varvals;
	vector<string> varvalsTemp;
	string entrada1;
	vector<vector<int>> grafo;
	string dataText;
	int filas = 0;
	int cols = 0;
	long n = 0;
	int nAristas = 0;
	vector<int> temp;
	double alpha = 1.0;

	do
	{
		cout << endl;
		getline(cin, entrada1);
		string delimiter = " ";
		size_t pos = 0;
		string archivo;


		pos = entrada1.find(delimiter);
		dataText = entrada1.substr(0, pos);

		if (dataText == "demo") {

			data = { {1,0,0,0},{0,0,1,1},{0,1,0,0},{1,2,0,0},{0,2,1,1},{2,1,1,1} };
			varnames = { "A","B","C","D" };
			varvals = { {"0","1","2"},{"0","1","2"},{"0","1"},{"0","1"} };
			card = { 3,3,2,2 };
			filas = 6;
			cols = 4;

			cout << "Instancias: " << filas << endl;
			cout << "Variables : " << cols << endl;

			n = (int)pow(2, ((pow(cols, 2)) - cols));
			nAristas = (int)pow(cols, 2) - cols;
			grafo.resize(cols, vector<int>(cols));

			int contador = 0;
			for (int w = 0; w < n; w++) {
				temp = HipotsisMatrix(w, nAristas);
				vector<vector<int>> result(cols, vector<int>(cols));
				int k = 0;
				for (int i = 0; i < cols; i++) {
					for (int j = 0; j < cols; j++) {

						if (i == j) {
							result[i][j] = 0;
						}
						else {
							result[i][j] = temp[k];
							k++;
						}
					}
				}

				Ciclico nuevoCiclico(result);

				if (nuevoCiclico.hasCycle()) {
					continue;
				}
				else
				{
					contador++;
					double e = 0.0;
					double aic = 0.0;
					double k2 = 0.0;
					for (int q = 0; q < cols; q++) {

						vector<int> padres = getPadres(q, result);
						vector<int> vars = unio(q, padres);
						int tamFactor = 1;
						for (int l = 0; l < vars.size(); l++) {
							tamFactor *= card[vars[l]];
						}

						vector<double> probs(tamFactor);
						vector<double> conta(tamFactor);
						vector<double> probsNorm(tamFactor);

						probs = calcularEntropia(probs, vars, card, data, alpha);

						int cardAIC = 0;
						if (padres.size() == 0) {
							cardAIC = 1;
						}
						else {

							for (int i = 0; i < padres.size(); i++) {

								cardAIC += card[padres[i]];

							}

						}

						aic += (card[vars[0]] - 1)*cardAIC;



						for (int b = 0; b < probs.size(); b++) {
							e += probs[b];
						}

						k2 = GetK2(result, data, card, cols, 543);



					}

					cout << "Grafo: " << contador << " MDL: " << e + (aic / 2)*log10(filas) << endl;
					cout << "Grafo: " << contador << " Entropia: " << e << endl;
					cout << "Grafo: " << contador << " AIC: " << e - aic << endl;
					cout << "Grafo: " << contador << " K2: " << k2 << endl;
				}


			}

			cout << "----------------------------------------------- K2 ---------------------------------------------------";

			int contador2 = 0;

			for (int w = 0; w < n; w++) {
				bool esK2 = false;
				temp = HipotsisMatrix(w, nAristas);
				vector<vector<int>> result(cols, vector<int>(cols));
				int k = 0;
				for (int i = 0; i < cols; i++) {
					for (int j = 0; j < cols; j++) {

						if (i == j) {
							result[i][j] = 0;
						}
						else {
							result[i][j] = temp[k];
							k++;

							if (j >= i && result[i][j] == 1) {
								esK2 = true;
							}

						}


					}
				}

				if (esK2 == true) {
					continue;
				}

				Ciclico nuevoCiclico(result);

				if (nuevoCiclico.hasCycle()) {
					continue;
				}
				else
				{
					contador2++;
					double e = 0.0;
					double aic = 0.0;
					double k2 = 0.0;
					for (int q = 0; q < cols; q++) {

						vector<int> padres = getPadres(q, result);
						vector<int> vars = unio(q, padres);
						int tamFactor = 1;
						for (int l = 0; l < vars.size(); l++) {
							tamFactor *= card[vars[l]];
						}

						vector<double> probs(tamFactor);
						vector<double> conta(tamFactor);
						vector<double> probsNorm(tamFactor);

						probs = calcularEntropia(probs, vars, card, data, alpha);

						int cardAIC = 0;
						if (padres.size() == 0) {
							cardAIC = 1;
						}
						else {

							for (int i = 0; i < padres.size(); i++) {

								cardAIC += card[padres[i]];

							}

						}

						aic += (card[vars[0]] - 1)*cardAIC;


						for (int b = 0; b < probs.size(); b++) {
							e += probs[b];
						}

						k2 = GetK2(result, data, card, cols, 64);

					}

					cout << "Grafo: " << contador2 << " Entropia: " << e + (aic / 2)*log10(filas) << endl;
					cout << "Grafo: " << contador2 << " AIC: " << e - aic << endl;
					cout << "Grafo: " << contador2 << " MDL: " << e + (aic / 2)*log10(filas) << endl;
					cout << "Grafo: " << contador2 << " K2: " << k2 << endl;

				}


			}


		}


		if (dataText == "salir") {
			continue;
		}

		if (dataText == "alpha") {

			entrada1.erase(0, pos + delimiter.length());
			alpha = stod(entrada1);

		}

		if (dataText == "data") {
			//filas = 0;
			//cols = 0;
			entrada1.erase(0, pos + delimiter.length());
			archivo = entrada1;

			ifstream file(archivo);
			string line, field;

			while (getline(file, line)) {
				filas++;
				stringstream ss(line);
				cols = 0;
				while (getline(ss, field, ','))
				{
					cols++;
					int num = stoi(field.c_str());
					v.push_back(num);
				}
				data.push_back(v);
				v.clear();
			}

			cout << "Instancias: " << filas << endl;
			cout << "Variables : " << cols << endl;

			n = (int)pow(2, ((pow(cols, 2)) - cols));
			nAristas = (int)pow(cols, 2) - cols;
			grafo.resize(cols, vector<int>(cols));

		}

		if (dataText == "vars") {

			entrada1.erase(0, pos + delimiter.length());
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				string temp = entrada1.substr(0, pos);
				varnames.push_back(temp);
				entrada1.erase(0, pos + delimiter.length());
			}
			varnames.push_back(entrada1);

			cout << endl;
			for (int i = 0; i < varnames.size(); i++) {
				cout << "[" << i + 1 << "] " << varnames[i] << endl;
			}
		}

		if (dataText == "arco") {

			entrada1.erase(0, pos + delimiter.length());
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				string temp = entrada1.substr(0, pos);

				for (int i = 0; i < varnames.size(); i++) {

					if (temp == varnames[i]) {
						entrada1.erase(0, pos + delimiter.length());

						for (int j = 0; j < varnames.size(); j++)
						{
							if (entrada1 == varnames[j]) {
								grafo[i][j] = 1;
								continue;
							}
						}
					}

				}
			}
		}

		if (dataText == "dag") {

			for (int i = 0; i < grafo.size(); i++) {

				for (int j = 0; j < grafo.size(); j++) {

					cout << grafo[i][j] << " ";
				}

				cout << endl;

			}

		}

		if (dataText == "vals") {

			entrada1.erase(0, pos + delimiter.length());
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				string temp = entrada1.substr(0, pos);

				for (int i = 0; i < varnames.size(); i++) {

					if (temp == varnames[i]) {

						while ((pos = entrada1.find(delimiter)) != std::string::npos) {
							entrada1.erase(0, pos + delimiter.length());
							pos = entrada1.find(delimiter);
							string temp = entrada1.substr(0, pos);
							varvalsTemp.push_back(temp);
						}
						varvals.push_back(varvalsTemp);
						varvalsTemp.clear();
						break;

					}

				}

			}


		}


		if (dataText == "card") {

			for (int i = 0; i < varvals.size(); i++)
			{
				cout << "Var(" << i << ") = " << varvals[i].size() << endl;
				card.push_back(varvals[i].size());
			}

		}

		if (dataText == "prob") {

			vector<int>vars;
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				entrada1.erase(0, pos + delimiter.length());
				pos = entrada1.find(delimiter);
				string temp = entrada1.substr(0, pos);

				for (int i = 0; i < varnames.size(); i++) {

					if (temp == varnames[i]) {
						vars.push_back(i);
						vector<double>probs(card[i]);
						CalcularProbs(0, probs, vars, card, data, varnames, varvals, false, card[i], alpha);
					}
				}
			}
		}

		if (dataText == "probConj") {

			vector<double>probs;
			vector<string>temporal;
			vector<int>tempNum;
			vector<int>vars;
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				entrada1.erase(0, pos + delimiter.length());
				pos = entrada1.find(delimiter);
				string temp = entrada1.substr(0, pos);
				temporal.push_back(temp);
			}

			int cardTotal = 1;
			for (int i = 0; i < temporal.size(); i++) {
				for (int j = 0; j < varnames.size(); j++) {
					if (temporal[i] == varnames[j]) {
						tempNum.push_back(j);
						cardTotal = cardTotal * card[j];
						break;
					}
				}
			}

			probs.resize(cardTotal);
			CalcularProbs(0, probs, tempNum, card, data, varnames, varvals, false, cardTotal, alpha);

		}


		if (dataText == "probCond") {

			vector<double>probs;
			vector<string>temporal;
			vector<int>tempNum;
			vector<int>vars;
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				entrada1.erase(0, pos + delimiter.length());
				pos = entrada1.find(delimiter);
				string temp = entrada1.substr(0, pos);
				temporal.push_back(temp);
			}

			int cardTotal = 1;
			for (int i = 0; i < temporal.size(); i++) {
				for (int j = 0; j < varnames.size(); j++) {
					if (temporal[i] == varnames[j]) {
						tempNum.push_back(j);
						cardTotal = cardTotal * card[j];
						break;
					}
				}
			}

			probs.resize(cardTotal);
			CalcularProbs(0, probs, tempNum, card, data, varnames, varvals, true, cardTotal, alpha);

		}

		if (dataText == "MI") {

			vector<double>probs;
			vector<string>temporal;
			vector<int>tempNum;
			vector<int>tempVars;
			vector<int>vars;
			vector<int>tempNum2;
			vector<int>tempNum3;

			double val1, val2, val3, val4, val5, valfinal;
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				entrada1.erase(0, pos + delimiter.length());
				pos = entrada1.find(delimiter);
				string temp = entrada1.substr(0, pos);
				temporal.push_back(temp);
			}


			int cardTotal = 1;
			for (int i = 0; i < temporal.size(); i++) {
				for (int j = 0; j < varnames.size(); j++) {
					if (temporal[i] == varnames[j]) {
						tempNum.push_back(j);
						cardTotal = cardTotal * card[j];
						break;
					}
				}
			}

			tempNum2.push_back(tempNum[0]);
			tempNum3.push_back(tempNum[1]);

			for (int i = 2; i < tempNum.size(); i++) {

				tempVars.push_back(tempNum[i]);
				tempNum2.push_back(tempNum[i]);
				tempNum3.push_back(tempNum[i]);

			}

			double valorFinal = 0.0;
			double valorTerminal = 1.0;
			int cardAnt = 1;
			int iVar;
			vector<int> vals(temporal.size());
			vector<int> vals2(temporal.size() - 1);
			vector<int> vals3(temporal.size() - 1);
			vector<int> valsK(temporal.size() - 2);

			vector<int>numVars(temporal.size() - 2);
			for (int i = 0; i < card[tempNum[0]]; i++) {
				for (int j = 0; j < card[tempNum[1]]; j++) {
					for (int k = 0; k < cardTotal / (card[tempNum[0]] * card[tempNum[1]]); k++) {
						cardAnt = 1;
						vals[0] = i;
						vals[1] = j;
						vals2[0] = i;
						vals3[0] = j;
						//vals2[0] = i;
						//vals3[0] = j;

						if (tempVars.size() == 0) {
							valsK = vals;

							val1 = CalcularProbConjunta(tempNum, data, vals, card, alpha, true, cardTotal);
							val2 = val1;
							val3 = CalcularProbConjunta(tempNum2, data, vals2, card, alpha, true, card[tempNum2[0]] * (cardTotal / (card[tempNum[0]] * card[tempNum[1]])));
							val4 = CalcularProbConjunta(tempNum3, data, vals3, card, alpha, true, card[tempNum3[0]] * (cardTotal / (card[tempNum[0]] * card[tempNum[1]])));

							val5 = log2(val2 / (val3*val4));
							valfinal = val1 * val5;

						}
						else {

							for (int b = 0; b < temporal.size() - 2; b++) {
								iVar = tempVars[b];
								numVars[b] = (int)floor(k / cardAnt) % card[iVar];
								cardAnt *= card[iVar];
								vals[2 + b] = numVars[b];
								vals2[1 + b] = numVars[b];
								vals3[1 + b] = numVars[b];
								valsK[b] = numVars[b];
							}

							val1 = CalcularProbConjunta(tempNum, data, vals, card, alpha, true, cardTotal);
							val2 = val1 / CalcularProbConjunta(tempVars, data, valsK, card, alpha, true, cardTotal / (card[tempNum[0]] * card[tempNum[1]]));
							val3 = CalcularProbConjunta(tempNum2, data, vals2, card, alpha, true, card[tempNum2[0]] * (cardTotal / (card[tempNum[0]] * card[tempNum[1]]))) / CalcularProbConjunta(tempVars, data, valsK, card, alpha, true, cardTotal / (card[tempNum[0]] * card[tempNum[1]]));
							val4 = CalcularProbConjunta(tempNum3, data, vals3, card, alpha, true, card[tempNum3[0]] * (cardTotal / (card[tempNum[0]] * card[tempNum[1]]))) / CalcularProbConjunta(tempVars, data, valsK, card, alpha, true, cardTotal / (card[tempNum[0]] * card[tempNum[1]]));

							val5 = log2(val2 / (val3*val4));
							valfinal = val1 * val5;

						}


						valorFinal += valfinal;

					}


				}

			}

			cout << valorFinal;

		}

		if (dataText == "Pearson") {

			vector<double>probs;
			vector<string>temporal;
			vector<int>tempNum;
			vector<int>tempVars;
			vector<int>vars;
			vector<int>tempNum2;
			vector<int>tempNum3;

			double val1, val2, val3, val4, val5, valfinal;
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				entrada1.erase(0, pos + delimiter.length());
				pos = entrada1.find(delimiter);
				string temp = entrada1.substr(0, pos);
				temporal.push_back(temp);
			}


			int cardTotal = 1;
			for (int i = 0; i < temporal.size(); i++) {
				for (int j = 0; j < varnames.size(); j++) {
					if (temporal[i] == varnames[j]) {
						tempNum.push_back(j);
						cardTotal = cardTotal * card[j];
						break;
					}
				}
			}

			tempNum2.push_back(tempNum[0]);
			tempNum3.push_back(tempNum[1]);

			for (int i = 2; i < tempNum.size(); i++) {

				tempVars.push_back(tempNum[i]);
				tempNum2.push_back(tempNum[i]);
				tempNum3.push_back(tempNum[i]);

			}

			double valorFinal = 0.0;
			double valorTerminal = 1.0;
			int cardAnt = 1;
			int iVar;
			vector<int> vals(temporal.size());
			vector<int> vals2(temporal.size() - 1);
			vector<int> vals3(temporal.size() - 1);
			vector<int> valsK(temporal.size() - 2);

			vector<int>numVars(temporal.size() - 2);
			for (int i = 0; i < card[tempNum[0]]; i++) {
				for (int j = 0; j < card[tempNum[1]]; j++) {
					for (int k = 0; k < cardTotal / (card[tempNum[0]] * card[tempNum[1]]); k++) {
						cardAnt = 1;
						vals[0] = i;
						vals[1] = j;
						vals2[0] = i;
						vals3[0] = j;
						//vals2[0] = i;
						//vals3[0] = j;


						for (int b = 0; b < temporal.size() - 2; b++) {
							iVar = tempVars[b];
							numVars[b] = (int)floor(k / cardAnt) % card[iVar];
							cardAnt *= card[iVar];
							vals[2 + b] = numVars[b];
							vals2[1 + b] = numVars[b];
							vals3[1 + b] = numVars[b];
							valsK[b] = numVars[b];
						}

						val1 = CalcularProbConjuntaFinal(tempNum, data, vals);
						val2 = (CalcularProbConjuntaFinal(tempNum2, data, vals2) * CalcularProbConjuntaFinal(tempNum3, data, vals3));
						val4 = CalcularProbConjuntaFinal(tempVars, data, valsK);

						if (val4 == 0) {
							val5 = 0;
						}
						else {
							val5 = val2 / val4;
						}

						val3 = pow(val1 - val5, 2);

						if (val5 == 0) {
							valfinal = 0;
						}
						else {
							valfinal = val3 / val5;
						}



						valorFinal += valfinal;
					}

				}


			}

			cout << valorFinal;
		}








		if (dataText == "entropiaFB") {

			int contador = 0;
			for (int w = 0; w < n; w++) {
				temp = HipotsisMatrix(w, nAristas);
				vector<vector<int>> result(cols, vector<int>(cols));
				int k = 0;
				for (int i = 0; i < cols; i++) {
					for (int j = 0; j < cols; j++) {

						if (i == j) {
							result[i][j] = 0;
						}
						else {
							result[i][j] = temp[k];
							k++;
						}
					}
				}

				Ciclico nuevoCiclico(result);

				if (nuevoCiclico.hasCycle()) {
					continue;
				}
				else
				{
					contador++;
					double e = 0.0;
					for (int q = 0; q < cols; q++) {

						vector<int> padres = getPadres(q, result);
						vector<int> vars = unio(q, padres);
						int tamFactor = 1;
						for (int l = 0; l < vars.size(); l++) {
							tamFactor *= card[vars[l]];
						}

						vector<double> probs(tamFactor);
						vector<double> conta(tamFactor);
						vector<double> probsNorm(tamFactor);

						probs = calcularEntropia(probs, vars, card, data, alpha);

						for (int b = 0; b < probs.size(); b++) {
							e += probs[b];
						}

					}

					cout << "Grafo: " << contador << " Entropia: " << e << endl;

				}


			}



		}

		if (dataText == "aicFB") {

			int contador = 0;
			for (int w = 0; w < n; w++) {
				temp = HipotsisMatrix(w, nAristas);
				vector<vector<int>> result(cols, vector<int>(cols));
				int k = 0;
				for (int i = 0; i < cols; i++) {
					for (int j = 0; j < cols; j++) {

						if (i == j) {
							result[i][j] = 0;
						}
						else {
							result[i][j] = temp[k];
							k++;
						}
					}
				}

				Ciclico nuevoCiclico(result);

				if (nuevoCiclico.hasCycle()) {
					continue;
				}
				else
				{
					contador++;
					double e = 0.0;
					double aic = 0.0;
					for (int q = 0; q < cols; q++) {

						vector<int> padres = getPadres(q, result);
						vector<int> vars = unio(q, padres);
						int tamFactor = 1;
						for (int l = 0; l < vars.size(); l++) {
							tamFactor *= card[vars[l]];
						}

						vector<double> probs(tamFactor);
						vector<double> conta(tamFactor);
						vector<double> probsNorm(tamFactor);

						probs = calcularEntropia(probs, vars, card, data, alpha);

						int cardAIC = 0;
						if (padres.size() == 0) {
							cardAIC = 1;
						}
						else {

							for (int i = 0; i < padres.size(); i++) {

								cardAIC += card[padres[i]];

							}

						}

						aic += (card[vars[0]] - 1)*cardAIC;


						for (int b = 0; b < probs.size(); b++) {
							e += probs[b];
						}

					}

					cout << "Grafo: " << contador << " AIC: " << e - aic << endl;

				}


			}
		}

		if (dataText == "mdlFB") {

			int contador = 0;
			for (int w = 0; w < n; w++) {
				temp = HipotsisMatrix(w, nAristas);
				vector<vector<int>> result(cols, vector<int>(cols));
				int k = 0;
				for (int i = 0; i < cols; i++) {
					for (int j = 0; j < cols; j++) {

						if (i == j) {
							result[i][j] = 0;
						}
						else {
							result[i][j] = temp[k];
							k++;
						}
					}
				}

				Ciclico nuevoCiclico(result);

				if (nuevoCiclico.hasCycle()) {
					continue;
				}
				else
				{
					contador++;
					double e = 0.0;
					double aic = 0.0;
					for (int q = 0; q < cols; q++) {

						vector<int> padres = getPadres(q, result);
						vector<int> vars = unio(q, padres);
						int tamFactor = 1;
						for (int l = 0; l < vars.size(); l++) {
							tamFactor *= card[vars[l]];
						}

						vector<double> probs(tamFactor);
						vector<double> conta(tamFactor);
						vector<double> probsNorm(tamFactor);

						probs = calcularEntropia(probs, vars, card, data, alpha);

						int cardAIC = 0;
						if (padres.size() == 0) {
							cardAIC = 1;
						}
						else {

							for (int i = 0; i < padres.size(); i++) {

								cardAIC += card[padres[i]];

							}

						}

						aic += (card[vars[0]] - 1)*cardAIC;



						for (int b = 0; b < probs.size(); b++) {
							e += probs[b];
						}

					}

					cout << "Grafo: " << contador << " MDL: " << e + (aic / 2)*log10(filas) << endl;

				}


			}
		}

		if (dataText == "k2FB") {

			int contador = 0;
			//vector<double> k2Entropia;
			//bool K2STOP = false;


			for (int w = 0; w < n; w++) {
				temp = HipotsisMatrix(w, nAristas);
				vector<vector<int>> result(cols, vector<int>(cols));
				int k = 0;
				for (int i = 0; i < cols; i++) {
					for (int j = 0; j < cols; j++) {

						if (i == j) {
							result[i][j] = 0;
						}
						else {
							result[i][j] = temp[k];
							k++;
						}
					}
				}

				Ciclico nuevoCiclico(result);

				if (nuevoCiclico.hasCycle()) {
					continue;
				}
				else
				{

					contador++;
					double k2 = 0.0;
					k2 = GetK2(result, data, card, cols, 543);
					cout << "Grafo: " << contador << " K2: " << k2 << endl;
					//k2Entropia.push_back(k2);

				}


			}
		}


		if (dataText == "entropiaK2") {

			int contador = 0;
			vector<vector<int>> grafo(cols, vector<int>(cols));
			vector<vector<int>> grafoCopia = grafo;
			vector<vector<int>> grafoCopia2 = grafo;
			vector<double> k2Entropia;
			bool esK2 = false;
			int esOnoes = 0;

			double eInicial = 0.0;
			double eCopia = 0.0;
			for (int q = 0; q < cols; q++) {

				vector<int> padres = getPadres(q, grafo);
				vector<int> vars = unio(q, padres);
				int tamFactor = 1;
				for (int l = 0; l < vars.size(); l++) {
					tamFactor *= card[vars[l]];
				}

				vector<double> probs(tamFactor);
				vector<double> conta(tamFactor);
				vector<double> probsNorm(tamFactor);

				probs = calcularEntropia(probs, vars, card, data, alpha);

				for (int b = 0; b < probs.size(); b++) {
					eInicial += probs[b];
				}

			}

			eCopia = eInicial;

			cout << "Grafo 0: Entropia: " << eInicial << endl;


			do {

				esOnoes = 0;

				for (int i = 0; i < grafoCopia.size(); i++) {
					for (int j = 0; j < grafoCopia.size(); j++) {

						if (grafoCopia[i][j] == 0 && !(j >= i)) {

							grafoCopia[i][j] = 1;

							Ciclico nuevoCiclico(grafoCopia);

							if (nuevoCiclico.hasCycle()) {
								continue;
							}
							else
							{
								contador++;
								double e = 0.0;
								for (int q = 0; q < cols; q++) {

									vector<int> padres = getPadres(q, grafoCopia);
									vector<int> vars = unio(q, padres);
									int tamFactor = 1;
									for (int l = 0; l < vars.size(); l++) {
										tamFactor *= card[vars[l]];
									}

									vector<double> probs(tamFactor);
									vector<double> conta(tamFactor);
									vector<double> probsNorm(tamFactor);

									probs = calcularEntropia(probs, vars, card, data, alpha);

									for (int b = 0; b < probs.size(); b++) {
										e += probs[b];
									}

								}

								cout << "Grafo: " << contador << " Entropia: " << e << endl;
								k2Entropia.push_back(e);

								if (e < eInicial) {
									eInicial = e;
									grafoCopia2 = grafoCopia;
									esOnoes++;
									grafoCopia[i][j] = 0;
								}

								else {
									grafoCopia[i][j] = 0;
								}




							}
						}




					}
				}

				grafoCopia = grafoCopia2;


			} while (esOnoes > 0);

		}


		if (dataText == "aicK2") {

			int contador = 0;
			vector<vector<int>> grafo(cols, vector<int>(cols));
			vector<vector<int>> grafoCopia = grafo;
			vector<vector<int>> grafoCopia2 = grafo;
			bool esK2 = false;
			int esOnoes = 0;

			double eInicial = 0.0;
			double eCopia = 0.0;
			double aicInicial = 0.0;
			double aicCopia = 0.0;
			double AICPRIMERO = 0.0;

			for (int q = 0; q < cols; q++) {

				vector<int> padres = getPadres(q, grafo);
				vector<int> vars = unio(q, padres);
				int tamFactor = 1;
				for (int l = 0; l < vars.size(); l++) {
					tamFactor *= card[vars[l]];
				}

				vector<double> probs(tamFactor);
				vector<double> conta(tamFactor);
				vector<double> probsNorm(tamFactor);

				probs = calcularEntropia(probs, vars, card, data, alpha);


				int cardAIC = 0;
				if (padres.size() == 0) {
					cardAIC = 1;
				}
				else {

					for (int i = 0; i < padres.size(); i++) {

						cardAIC += card[padres[i]];

					}

				}

				aicInicial += (card[vars[0]] - 1)*cardAIC;


				for (int b = 0; b < probs.size(); b++) {
					eInicial += probs[b];
				}

			}

			eCopia = eInicial;
			aicCopia = aicInicial;
			AICPRIMERO = eInicial - aicInicial;

			cout << "Grafo 0: AIC: " << AICPRIMERO << endl;


			do {
				esOnoes = 0;
				for (int i = 0; i < grafoCopia.size(); i++) {
					for (int j = 0; j < grafoCopia.size(); j++) {

						if (grafoCopia[i][j] == 0 && !(j >= i)) {

							grafoCopia[i][j] = 1;

							Ciclico nuevoCiclico(grafoCopia);

							if (nuevoCiclico.hasCycle()) {
								continue;
							}
							else
							{
								contador++;
								double e = 0.0;
								double aic = 0.0;
								for (int q = 0; q < cols; q++) {

									vector<int> padres = getPadres(q, grafoCopia);
									vector<int> vars = unio(q, padres);
									int tamFactor = 1;
									for (int l = 0; l < vars.size(); l++) {
										tamFactor *= card[vars[l]];
									}

									vector<double> probs(tamFactor);
									vector<double> conta(tamFactor);
									vector<double> probsNorm(tamFactor);

									probs = calcularEntropia(probs, vars, card, data, alpha);

									int cardAIC = 0;
									if (padres.size() == 0) {
										cardAIC = 1;
									}
									else {

										for (int i = 0; i < padres.size(); i++) {

											cardAIC += card[padres[i]];

										}

									}

									aic += (card[vars[0]] - 1)*cardAIC;




									for (int b = 0; b < probs.size(); b++) {
										e += probs[b];
									}

								}

								cout << "Grafo: " << contador << " Entropia: " << e - aic << endl;

								if ((e - aic) < (AICPRIMERO)) {
									AICPRIMERO = e - aic;
									esOnoes++;
									grafoCopia2 = grafoCopia;
									grafoCopia[i][j] = 0;
								}

								else {
									grafoCopia[i][j] = 0;
								}




							}
						}




					}
				}

				grafoCopia = grafoCopia2;


			} while (esOnoes > 0);
		}

		if (dataText == "mdlK2") {
			int contador = 0;
			vector<vector<int>> grafo(cols, vector<int>(cols));
			vector<vector<int>> grafoCopia = grafo;
			vector<vector<int>> grafoCopia2 = grafo;
			bool esK2 = false;
			int esOnoes = 0;

			double eInicial = 0.0;
			double eCopia = 0.0;
			double aicInicial = 0.0;
			double aicCopia = 0.0;
			double AICPRIMERO = 0.0;
			double MDLPRIMERO = 0.0;

			for (int q = 0; q < cols; q++) {

				vector<int> padres = getPadres(q, grafo);
				vector<int> vars = unio(q, padres);
				int tamFactor = 1;
				for (int l = 0; l < vars.size(); l++) {
					tamFactor *= card[vars[l]];
				}

				vector<double> probs(tamFactor);
				vector<double> conta(tamFactor);
				vector<double> probsNorm(tamFactor);

				probs = calcularEntropia(probs, vars, card, data, alpha);


				int cardAIC = 0;
				if (padres.size() == 0) {
					cardAIC = 1;
				}
				else {

					for (int i = 0; i < padres.size(); i++) {

						cardAIC += card[padres[i]];

					}

				}

				aicInicial += (card[vars[0]] - 1)*cardAIC;


				for (int b = 0; b < probs.size(); b++) {
					eInicial += probs[b];
				}

			}

			eCopia = eInicial;
			aicCopia = aicInicial;
			AICPRIMERO = eInicial - aicInicial;
			MDLPRIMERO = eInicial + (aicInicial / 2)*log10(filas);

			cout << "Grafo 0: MDL: " << MDLPRIMERO << endl;


			do {
				esOnoes = 0;
				for (int i = 0; i < grafoCopia.size(); i++) {
					for (int j = 0; j < grafoCopia.size(); j++) {

						if (grafoCopia[i][j] == 0 && !(j >= i)) {

							grafoCopia[i][j] = 1;

							Ciclico nuevoCiclico(grafoCopia);

							if (nuevoCiclico.hasCycle()) {
								continue;
							}
							else
							{
								contador++;
								double e = 0.0;
								double aic = 0.0;
								for (int q = 0; q < cols; q++) {

									vector<int> padres = getPadres(q, grafoCopia);
									vector<int> vars = unio(q, padres);
									int tamFactor = 1;
									for (int l = 0; l < vars.size(); l++) {
										tamFactor *= card[vars[l]];
									}

									vector<double> probs(tamFactor);
									vector<double> conta(tamFactor);
									vector<double> probsNorm(tamFactor);

									probs = calcularEntropia(probs, vars, card, data, alpha);

									int cardAIC = 0;
									if (padres.size() == 0) {
										cardAIC = 1;
									}
									else {

										for (int i = 0; i < padres.size(); i++) {

											cardAIC += card[padres[i]];

										}

									}

									aic += (card[vars[0]] - 1)*cardAIC;




									for (int b = 0; b < probs.size(); b++) {
										e += probs[b];
									}

								}

								cout << "Grafo: " << contador << " MDL: " << e + (aic / 2)*log10(filas) << endl;

								if ((e + (aic / 2)*log10(filas)) < (MDLPRIMERO)) {
									MDLPRIMERO = (e + (aic / 2)*log10(filas));
									esOnoes++;
									grafoCopia2 = grafoCopia;
									grafoCopia[i][j] = 0;
								}

								else {
									grafoCopia[i][j] = 0;
								}




							}
						}




					}
				}

				grafoCopia = grafoCopia2;


			} while (esOnoes > 0);
		}

		if (dataText == "k2K2") {

			int contador = 0;
			vector<vector<int>> grafo(cols, vector<int>(cols));
			vector<vector<int>> grafoCopia = grafo;
			vector<vector<int>> grafoCopia2 = grafo;
			vector<double> k2Entropia;
			bool esK2 = false;
			int esOnoes = 0;

			double eInicial = 0.0;
			double eCopia = 0.0;
			double k2Inicial = 0.0;
			k2Inicial = GetK2(grafo, data, card, cols, 64);
			cout << "Grafo: " << contador << " K2: " << k2Inicial << endl;

			do {

				esOnoes = 0;

				for (int i = 0; i < grafoCopia.size(); i++) {
					for (int j = 0; j < grafoCopia.size(); j++) {

						if (grafoCopia[i][j] == 0 && !(j >= i)) {

							grafoCopia[i][j] = 1;

							Ciclico nuevoCiclico(grafoCopia);

							if (nuevoCiclico.hasCycle()) {
								continue;
							}
							else
							{
								contador++;
								double k2 = 0.0;
								k2 = GetK2(grafoCopia, data, card, cols, 64);
								cout << "Grafo: " << contador << " K2: " << k2 << endl;

								if (k2 > k2Inicial) {
									k2Inicial = k2;
									grafoCopia2 = grafoCopia;
									esOnoes++;
									grafoCopia[i][j] = 0;
								}

								else {
									grafoCopia[i][j] = 0;
								}




							}
						}




					}
				}

				grafoCopia = grafoCopia2;


			} while (esOnoes > 0);

		}

		if (dataText == "CL") {

			vector<vector<int>> grafoCL(cols, vector<int>(cols));
			vector<double> ordenarMI;
			vector<vector<double>> ordenarMI2(6);
			int cuenta = 0;

			for (int i = 0; i < grafoCL.size(); i++) {
				for (int j = 0; j < grafoCL.size(); j++) {

					if ((j > i)) {

						int varA = i;
						int varB = j;

						vector<double>probs;
						vector<string>temporal;
						vector<int>tempNum;
						vector<int>tempVars;
						vector<int>vars;
						vector<int>tempNum2;
						vector<int>tempNum3;

						temporal.push_back(varnames[i]);
						temporal.push_back(varnames[j]);




						double val1, val2, val3, val4, val5, valfinal;
						/*while ((pos = entrada1.find(delimiter)) != std::string::npos) {
							entrada1.erase(0, pos + delimiter.length());
							pos = entrada1.find(delimiter);
							string temp = entrada1.substr(0, pos);
							temporal.push_back(temp);
						}*/



						int cardTotal = 1;
						for (int i = 0; i < temporal.size(); i++) {
							for (int j = 0; j < varnames.size(); j++) {
								if (temporal[i] == varnames[j]) {
									tempNum.push_back(j);
									cardTotal = cardTotal * card[j];
									break;
								}
							}
						}

						tempNum2.push_back(tempNum[0]);
						tempNum3.push_back(tempNum[1]);

						for (int i = 2; i < tempNum.size(); i++) {

							tempVars.push_back(tempNum[i]);
							tempNum2.push_back(tempNum[i]);
							tempNum3.push_back(tempNum[i]);

						}

						double valorFinal = 0.0;
						double valorTerminal = 1.0;
						int cardAnt = 1;
						int iVar;
						vector<int> vals(temporal.size());
						vector<int> vals2(temporal.size() - 1);
						vector<int> vals3(temporal.size() - 1);
						vector<int> valsK(temporal.size() - 2);

						vector<int>numVars(temporal.size() - 2);
						for (int i = 0; i < card[tempNum[0]]; i++) {
							for (int j = 0; j < card[tempNum[1]]; j++) {
								for (int k = 0; k < cardTotal / (card[tempNum[0]] * card[tempNum[1]]); k++) {
									cardAnt = 1;
									vals[0] = i;
									vals[1] = j;
									vals2[0] = i;
									vals3[0] = j;
									//vals2[0] = i;
									//vals3[0] = j;

									if (tempVars.size() == 0) {
										valsK = vals;

										val1 = CalcularProbConjunta(tempNum, data, vals, card, alpha, true, cardTotal);
										val2 = CalcularProbMarginal(i, data, tempNum[0], card, alpha);
										val3 = val2;
										val4 = CalcularProbMarginal(j, data, tempNum[1], card, alpha);
										//val3 = CalcularProbConjunta(tempNum2, data, vals2, card, alpha, true, card[tempNum2[0]] * (cardTotal / (card[tempNum[0]] * card[tempNum[1]])));
										//val4 = CalcularProbConjunta(tempNum3, data, vals3, card, alpha, true, card[tempNum3[0]] * (cardTotal / (card[tempNum[0]] * card[tempNum[1]])));

										val5 = (val2 / (val3 * val4));


										val5 = log10(val5);
										valfinal = val1 * val5;

									}

									valorFinal += valfinal;

								}


							}

						}

						//cout << valorFinal << endl;
						ordenarMI.push_back(valorFinal);



						ordenarMI2[cuenta].push_back(valorFinal);
						ordenarMI2[cuenta].push_back(varA);
						ordenarMI2[cuenta].push_back(varB);





						cuenta++;
					}


				}
			}

			/*std::sort(ordenarMI.begin(), ordenarMI.begin() + (ordenarMI.size()/2));
			std::sort(ordenarMI.begin() + (ordenarMI.size() / 2), ordenarMI.end(), myfunction); // 12 32 45 71(26 33 53 80)

																		 // using object as comp
			std::sort(ordenarMI.begin(), ordenarMI.end(), myobject);


			for (std::vector<double>::iterator it = ordenarMI.begin(); it != ordenarMI.end(); ++it) {

				cout << *it << endl;
			}*/


			std::sort(ordenarMI2.begin(), ordenarMI2.end(), sortcol);

			cout << "Finish";

			int num1 = 0;
			int num2 = 0;

			for (int k = 0; k < ordenarMI2.size(); k++)
			{
				bool esMI = false;
				num1 = ordenarMI2[k][1];
				num2 = ordenarMI2[k][2];

				for (int pl = 0; pl < grafoCL.size(); pl++) {

					if (grafoCL[pl][num2] == 1) {
						esMI = true;
					}


				}


				if (esMI == false) {
					grafoCL[num1][num2] = 1;
				}

			}

			grafoCLALL = grafoCL;
			
			cout << "Final Ahora sí";


		}

		if (dataText == "VE") {

			vector<double>probs;
			vector<string>temporal;
			vector<int>tempNum;
			vector<int>vars;
			while ((pos = entrada1.find(delimiter)) != std::string::npos) {
				entrada1.erase(0, pos + delimiter.length());
				pos = entrada1.find(delimiter);
				string temp = entrada1.substr(0, pos);
				temporal.push_back(temp);
			}

			for (int i = 0; i < temporal.size(); i++) {
				for (int j = 0; j < varnames.size(); j++) {
					if (temporal[i] == varnames[j]) {
						tempNum.push_back(j);
						break;
					}
				}
			}


			vector<vector<int>> Matrix = grafoCLALL;

			vector<int> inference = {tempNum[0],tempNum[1]};

			vector<vector<int>> Dist;

			vector<int> NonObs;

			double totalP1 = 0.0;

			double total = 1.0;

			int Evidencia = stoi(temporal[2]);

			for (int q = 0; q < Matrix.size(); q++) {

				vector<int> padres = getPadres(q, Matrix);
				vector<int> vars = unio(q, padres);

				Dist.push_back(vars);

			}

			for (int p = 0; p < Matrix.size(); p++) {
				if (find(inference.begin(), inference.end(), p) != inference.end()) {
					continue;
				}
				else {
					NonObs.push_back(p);
				}
			}

			int count = 0;

			for (int i = 0; i < NonObs.size(); i++) {
				for (int x = 0; x < Dist.size(); x++) {
					for (int y = 0; y < Dist[x].size(); y++) {
						if (Dist[x][y] == NonObs[i]) {
							Dist[x].erase(Dist[x].begin() + y);
						}
					}
				}
			}

			for (int i = 0; i < Dist.size(); i++) {
				if (Dist[i].size() == 0) {
					Dist.erase(Dist.begin() + i);
					i--;
				}
			}

			for (int n = 0; n < Dist.size(); n++) {
				vector<int> Vars = Dist[n];

				int FactorSize = 1;
				for (int v = 0; v < Vars.size(); v++) {
					FactorSize *= card[Vars[v]];
				}

				vector<double> Probs(FactorSize);
				Probs = calcularProbabilidadesF(Probs, Vars, card, data, alpha);


				vector<double> selecionados;

				double vf = 0.0;

				if (Vars.size() > 1 && Vars[1] == 1) {
					int cardPadre = card[Vars[1]];
					for (int p = cardPadre*Evidencia; p < cardPadre*Evidencia + cardPadre; p++) {

						selecionados.push_back(Probs[p]);
					}

					vf = *max_element(std::begin(selecionados), std::end(selecionados));


				}

				else if (Vars.size() > 1 && Vars[1] != 1) {

					int cardHijo = card[Vars[0]];
					int cardPadre = card[Vars[1]];
					for (int p = Evidencia; p < Probs.size(); p = p + (cardHijo)) {
						selecionados.push_back(Probs[p]);
					}

					vf = *max_element(std::begin(selecionados), std::end(selecionados));
				}

				else if (Vars.size() == 1 && Vars[0] != 1) {

					vf = *max_element(std::begin(Probs), std::end(Probs));

				}

				else {

					vf = Probs[Evidencia];

				}


				total *= vf;


			}

			cout << total;







		}


	} while (dataText != "salir");

	int N = varnames.size();

	getchar();
	return 0;
}