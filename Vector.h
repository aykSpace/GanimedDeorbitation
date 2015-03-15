#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include "Date_Time.h"
#include <iomanip>
#include "Atmosphera.h"
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
//Al ����������� ������ �����

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
	//������������ ������ 
	Vector()
		:vitok(0), t(0), x(0), y(0), z(0), Vx(0), Vy(0), Vz(0), Sb(0), fi(0), lambda(0), r(0), h(0), ro(0), _ax(0), _ay(0), _az(0)
	{};
	Vector(unsigned long _vitok, double _t, double _x, double _y, double _z, double _vx, double _vy, double _vz, double _sb, Date_Time _dt,
		double _fi = 0, double _lambda = 0, double _r = 0, double _h = 0, double _ro = 0);
	virtual ~Vector();

	const static unsigned SizeOfVector = 14;

	//���������� ����������
	Vector &       operator = (const Vector &);
	Vector &       operator = (double *);
	friend Vector  operator + (const Vector &, const Vector &);
	friend Vector  operator - (const Vector &, const Vector &);
    friend Vector  operator * (      double  , const Vector &);
    friend Vector  operator * (const Vector &, double);
	const  double  operator [](size_t) const;

	// �-� ����� ������
	void Perevod(); // ��������� ������ � ����������� ���������� (fi, Lambda, r) 
	void PerevodGanimed(); // ��������� ������ � ����������� ���������� (fi, Lambda, r) 
    void Print() const;
	Vector & Prognoz(Date_Time &_datTime); // ���������� ������ �� �������� �����
    vector<Vector> Prognoz (unsigned _NVitka); //���������� ������� �� ������ �� ������ ��������������
	Vector & Prognoz(); //����� ���������� ������� �� 1 ����� (���������� ��� ���������� ProgressBar)
	Vector  Prognoz(double u); // ���������� ������ �� �������� �������� ������ ������� � ��������
    virtual void RKS (double dt, Atmosphera &Atm); //�-� �������������� ������� ����� - �����
	Vector  F (Vector &, double dt, Atmosphera &Atm); // �-� �������� ������ ������
	virtual inline void Altitude(); // �-� �������� ������	
	double CalcCornOfSun(const Atmosphera & _atm); // ������ ���� ����������� �� ������
	double DLEkv (double _l1, double _l2); //������ ������������ ���������� �� ��������
	unsigned nSutCirc (); //������ ������ ��������� ����� (���� ��� ������������ ����������)
	unsigned nSutCirc (double _dLekv); // ������ ������ ��������� ����� (���� �c�� ����������� ����������)
	double Tet();//������ ���� ������� ���������� (���)

	void Polinom(bool & kk, int n, int m, double* ti, double* xi, double t, double* a, double& xa); //������������� �-��� ��������� n-1 ������� �� m ������	
	//xi - ������ �-���, ti - ������ ���������, t - ��������, a - ������ �-��� ��������, xa - ���������, kk = false - ����-�� �� ��������� 

	unsigned long GetCircle(){return vitok;}
	double h, t;

	// �-� ������ 
	friend double Sinusl (const Vector &,  int); // �-� ��������� ��������� ����
	friend double Cosinusl (const Vector &,  int); // �-� ��������� ��������� ����
	friend void gpz (const Vector &, double *);// �-� �������� ��������������� ���� �����
	friend void InterpolMatrix (double ** , unsigned OrderMatr, int i, Vector & Vect, unsigned &j, double arg = 0,
		                         bool *Half_Vitok = false, bool *Izm_Vitka = false, bool Progn = false); //�-� �������� ������� ��� ������������ Progn=true - ��������� ������ ����� �����
	friend void Ilagr (unsigned OrderMatr, double **, unsigned NA, double arg, Vector & Vect, unsigned NumOfCols);// �-� ������������ �������� (���������� ������ �� ��������� ��������
																													//��������� arg) NA (����� ������� � ����������)
	friend void Ilagr (unsigned _orderMatrEnt, double **_matrEnter, unsigned _numOfArgCol, double arg ,double *_rezult, unsigned _countOfColsRez); //���������� ������ ���������� ������������
																															//matrEnter (������� ������������ ������� �� ������� ������� �������)
	friend void GetDensity (int F, double Ap, Atmosphera &, Vector & Vect); // ��������� ���������� ��������� ���������
	friend void PerevodV(Vector &Vect); // ��������� ������ � ����������� ���������� (fi, Lambda, r) 

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
	 //����-�� ���
	double C[665];
	double D[665];
	double C0[36];
	double SQ[75];
	void gpz36 (const Vector& Vect, double* Rez); // �-� �������� ��������������� ���� �����*/

	
	void RKS20 (double dt, Atmosphera &Atm); //�-� �������������� ������� ����� - �����
	Vector F20 (Vector &, double dt, Atmosphera &Atm); // �-� �������� ������ ������
	void Predict20(double tend);


	double ModuleV(){return sqrt(Vx*Vx + Vy*Vy + Vz*Vz);}
public:
	double _ax, _ay, _az;


};


template <typename T> int Sign (T val)
{
	return (val > 0) ? (1) : ((val < 0) ? (-1) : (0));
}


