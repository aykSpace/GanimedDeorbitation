#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include "Date_Time.h"
#include <iomanip>
#include "Atmosphera.h"
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
//Al коэффициент сжати€ «емли

#define PI 3.14159265358979
#define OMEGA_EARTH 7.29211587e-5 
#define MU 398600.5
#define PerGradRad 1.74532925e-2
#define PerRadGrad 57.2957795056
#define KoefPressure 0.003352824419 
#define GEarth 0.00980665

using namespace std;
class Vector
{
public:
	unsigned long vitok;
	double x, y, z, Vx, Vy, Vz, Sb, r, ro;
	bool calcVarAtm;

public:
	double lambda, fi;
	static int f, ap;
	Date_Time DatTime;
	//конструкторы класса 
	Vector()
		:vitok(0), t(0), x(0), y(0), z(0), Vx(0), Vy(0), Vz(0), Sb(0), fi(0), lambda(0), r(0), h(0), ro(0), _ax(0), _ay(0), _az(0)
	{};
	Vector(unsigned long _vitok, double _t, double _x, double _y, double _z, double _vx, double _vy, double _vz, double _sb, Date_Time _dt,
		double _fi = 0, double _lambda = 0, double _r = 0, double _h = 0, double _ro = 0);
	virtual ~Vector();

	const static unsigned SizeOfVector = 14;

	//перегрузка операторов
	Vector &       operator = (const Vector &);
	Vector &       operator = (double *);
	friend Vector  operator + (const Vector &, const Vector &);
	friend Vector  operator - (const Vector &, const Vector &);
    friend Vector  operator * (      double  , const Vector &);
    friend Vector  operator * (const Vector &, double);
	const  double  operator [](size_t) const;

	// ф-и члены класса
	void Perevod(); // переводит вектор в сферические координаты (fi, Lambda, r) 
	void PerevodGanimed(); // переводит вектор в сферические координаты (fi, Lambda, r) 
    void Print() const;
	Vector & Prognoz(Date_Time &_datTime); // возвращает вектор на заданное врем€
    vector<Vector> Prognoz (unsigned _NVitka); //возвращает вектора на каждый из витков интегрировани€
	Vector & Prognoz(); //метод возвращает прогноз на 1 виток (необходимо дл€ отбражени€ ProgressBar)
	Vector  Prognoz(double u); // возвращает вектор на заданный аргумент широты периге€ в градусах
    virtual void RKS (double dt, Atmosphera &Atm); //ф-€ интегрировани€ методом –унге -  утта
	Vector  F (Vector &, double dt, Atmosphera &Atm); // ф-€ рассчета правых частей
	virtual inline void Altitude(); // ф-€ рассчета высоты	
	double CalcCornOfSun(const Atmosphera & _atm); // расчет угла направлени€ на —олнце
	double DLEkv (double _l1, double _l2); //расчет межвиткового рассто€ни€ на экваторе
	unsigned nSutCirc (); //расчет номера суточного витка (если нет межвиткового рассто€ни€)
	unsigned nSutCirc (double _dLekv); // расчет номера суточного витка (если еcть межвитковое рассто€ние)
	double Tet();//расчет угла наклона траектории (рад)

	void Polinom(bool & kk, int n, int m, double* ti, double* xi, double t, double* a, double& xa); //аппроксимаци€ ф-ции полиномом n-1 степени по m точкам	
	//xi - массив ф-ции, ti - массив аргумента, t - аргумент, a - вектор к-тов полинома, xa - результат, kk = false - коэф-ты не посчитаны 

	unsigned long GetCircle(){return vitok;}
	double h, t;

	// ф-и друзь€ 
	friend double Sinusl (const Vector &,  int); // ф-€ понижени€ кратности угла
	friend double Cosinusl (const Vector &,  int); // ф-€ понижени€ кратности угла
	friend void gpz (const Vector &, double *);// ф-€ рассчета гравитационного пол€ «емли
	friend void InterpolMatrix (double ** , unsigned OrderMatr, int i, Vector & Vect, unsigned &j, double arg = 0,
		                         bool *Half_Vitok = false, bool *Izm_Vitka = false, bool Progn = false); //ф-€ рассчета матрицы дл€ интерпол€ции Progn=true - запускает расчет смены витка
	friend void Ilagr (unsigned OrderMatr, double **, unsigned NA, double arg, Vector & Vect, unsigned NumOfCols);// ф-€ интерпол€ции Ћагранжа (возвращает вектор по заданному значению
																													//аргумента arg) NA (номер столбца с аргументом)
	friend void Ilagr (unsigned _orderMatrEnt, double **_matrEnter, unsigned _numOfArgCol, double arg ,double *_rezult, unsigned _countOfColsRez); //возвращает массив результата интерпол€ции
																															//matrEnter (пор€док интерпол€ции зависит от пор€дка входной матрицы)
	friend void GetDensity (int F, double Ap, Atmosphera &, Vector & Vect); // процедура вычислени€ плотности атмосферы
	friend void PerevodV(Vector &Vect); // переводит вектор в сферические координаты (fi, Lambda, r) 

	double ** CreateSquareMatrix(unsigned size);
	void DeleteSquareMatrix(unsigned size, double** matrix);

	//////////////////////////////////////////////////////////
/*
	long double factorial(int n) 
	{
		if (n < 0 ) {
			return 0;
		}
		return !n ? 1 : n * factorial(n - 1);
	}*/
/*
private:
	 //коэф-ты √ѕ«
	double C[665];
	double D[665];
	double C0[36];
	double SQ[75];
	void gpz36 (const Vector& Vect, double* Rez); // ф-€ рассчета гравитационного пол€ «емли*/

	
	void RKS20 (double dt, Atmosphera &Atm); //ф-€ интегрировани€ методом –унге -  утта
	Vector F20 (Vector &, double dt, Atmosphera &Atm); // ф-€ рассчета правых частей
	void Predict20(double tend);


	double ModuleV(){return sqrt(Vx*Vx + Vy*Vy + Vz*Vz);}
public:
	double _ax, _ay, _az;


};


template <typename T> int Sign (T val)
{
	return (val > 0) ? (1) : ((val < 0) ? (-1) : (0));
}


