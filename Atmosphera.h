#pragma once
#include "Date_Time.h"

#define OMEGA_EARTH 7.29211587e-5

class Atmosphera
{
public:
	double A[3], B[3], C[4], N[2], Fi, D[3], E[7], L[3]; //коэффициенты модели дл€ различных F0
	double S0, SDSol, CDSol, ASol, d1900;  //звездное врем€, sin delta, cos delta, alfa, кол-во дней от начала года (1900)
	int d; //кол-во дней от начала года
	double DRo;   //средн€€ вариаци€ плотности
	double W[3]; // средн€€ вариаци€ ветра
	unsigned NMounth; //номер мес€ца
	double VZV; // скорость звука
	
	//конструкторы и деструкторы
	Atmosphera();
	~Atmosphera();

	//методы класса
	void   SelectKoef(int aF0, double h);                     //процедура выбора коэффициентов из таблицы
	void   equal(unsigned n, double *mass1, double *mass2);  //процедура присваивани€ зн-й массива2 массиву 1, n - размерность массивов
	void   StarTime (Date_Time & Data);                     //процедура рассчета звездного времени на гринвичскую полночь
	void   Sun (double t);                                 //процедура рассчета положени€ —олнца на текущий момент времени (сначала надо вызвать StarTime!)
	double Ad_D  (double days_year) const;                //процедура выбора коэффициента ј(d) полугодовой эффект 
	double Kp_Ap (double App) const;                     //процедура пересчета из Ap в  р
	int    Selection_F0 (int F);                        //процедура выбора ближайшего з-€ F0 из таблицы    
	double LowAtm (double _h);						   //процедура вычислени€ плотности нижней атмосферы
	void   VarAtm (double H, double Fi, double Lam);  //процедура учета вариаций плотности атмосферы и ветра
	void   LinInterpolMatr(unsigned Col, const double Matr[][10], unsigned NumArg, double DArg, double *RezArr); //процедура линейной интерпол€ции дл€ матриц ветра и плотности
	                                                                                                       //(NumArg-номер аргумента после poi, DArg-разница между табл зн-ем реальным Rez - результирующий массив)
};

inline double lin2 (double x1, double x2, double y, double y1, double y2)// процедура лииейной интерпол€ции
{ 
	return ((y - y2)*(x2 - x1)/(y2 - y1)) + x2;
}
double lin2Matr(unsigned nRow, unsigned nCol, double ** matr, unsigned row1, unsigned col1, unsigned row2, unsigned col2, double argRow, double argCol); //линейна€ интерпол€ци€ в двухмерном массиве
																															//(argCol, argRow - значени€ аргументов по строке и столбцу)

double lin3Matr(unsigned n, unsigned m, unsigned l, double *** matr, unsigned i, unsigned j, unsigned k, unsigned i1, 
				unsigned j1, unsigned k1, double di, double dj, double dk); //линейна€ интерпол€ци€ в трехмерном массиве //(argCol, argRow - значени€ аргументов по строке и столбцу)

int poi (int n, const double *mass, double par); //процедура поиска в массиве, упор€доченном по возрастанию //n - число элементов массива -1, *mass - массив(упор€доченный), par - параметр
int poi (int n, const int *mass, const double par); //процедура поиска в массиве, упор€доченном по возрастанию //n - число элементов массива -1, *mass - массив(упор€доченный), par - параметр
