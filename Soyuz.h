#pragma once
#include "spacecraft.h"
class Soyuz : public Spacecraft	
{

//protected:
public:
	void Aerodynamic(double max, double & s, double & s1, double & b) override;
private:
	double DCg[3]; // "бок" в момент разделени€
	double _lSa; //длина —ј
	double _hDivision; //высота разделени€
	double _tDivsion;// врем€ разделени€
	double _vs;  //кажуща€с€ скорость (интеграл от продольного а/д ускорени€)
	unsigned _s; // уставка на направление крена (S = 0 начальный крен правый, S = 1 начальный крен левый)
	unsigned _vsr; // уставка на переворот
	double _gamma, _dgamma, _gammaCom; //приборный угол крена, програмный угол крена и углова€ скорость по крену (приборному)
	double _pxMax; //максимальна€ перегрузка
	bool _guidenceStart; //метка начала управлени€
	//коэффициенты дл€ аэродинамики
	double ** GetTdCx();
	double ** GetTdCy();
	double ** GetTdMz();
	double ** GetTdCxa();
	double ** GetTdCya();
	double ** GetTdMza();
	double ** GetTCx2();
	double ** GetTCy2();
	double *** GetTal();	
	//удаление матриц
	void DeleteMatr3(unsigned i, unsigned j, unsigned k, double *** matr);
	void DeleteMatr(unsigned nRow, unsigned nCol, double ** matr);

	double YtTma(); //вычисление ayt
	void AimPointContDeorb(double longitudeEquator, double & fiAim, double & lAim); //расчет координат прицельной точки ј”—, уст-ки S и первого прибл. VSR
	void ShootGuideDeorbit(const double tVkl, const double aimLambda, vector<Spacecraft> &result) const; //ф-€ рассчета точки посадки при управл€емом спуске(прострел) (tvkl=0 - первый прострел)
	void PredictGuideLandingPoint(vector<Spacecraft> &result);//интегрирование уравнений движени€ после активного участка
	void RKS (double dt, Atmosphera &Atm, Soyuz (*f_ptr) (Soyuz & spCr, double t , Atmosphera & atm)); //ф-€ интегрировани€ методом –унге -  утта	
	void RewriteVector(double * vec); //перезаписывает вектор из массива
	void FindPxMax(double & currentPxmax, double px);

public:
	//конструкторы
	Soyuz();
	Soyuz (const Vector &Vect, double _GKa, double _GBo, double _GSa, double _ImpSizeNu, double _Xt, double _Yt, double _Zt, double _Sbal,
		   double _RazpImpSize, double _NKa, unsigned _ROr, unsigned _RStab);
	Soyuz(const Spacecraft& Spcrafrt);
	~Soyuz();

	vector<Spacecraft> ControlledDeorb (const unsigned long _Vitok); //ф-€ рассчета управл€емого спуска  ј
	void SetAltitudeOfDivision(double hDivision){_hDivision = hDivision;}
	void SetLenghtOfDeorbSpCraft(double lDeorbSpCr){_lSa = lDeorbSpCr;}

	typedef int (*func) (int);
	friend int func1(Spacecraft i);
	friend int func2(int i);
	int testRef(func);
	int TestRef();

	//операторы
	friend Soyuz  operator + (const Soyuz &, const Soyuz &);
	friend Soyuz  operator - (const Soyuz &, const Soyuz &);
	friend Soyuz  operator * (          double  , const Soyuz &);
	friend Soyuz  operator * (const Soyuz &, double);

	//функции
	friend	Soyuz FGuideDeorb (Soyuz &, double dt, Atmosphera &Atm); // ф-€ рассчета правых частей дл€ управл€емого спуска
	friend void InterpolMatrix (double **AM, unsigned n, int i,  Soyuz & soyuz, unsigned &j, double arg, double arg1); //матрица дл€ интерпол€ции (может потом добавить проверку на виток)

};

