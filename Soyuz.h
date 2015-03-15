#pragma once
#include "spacecraft.h"
class Soyuz : public Spacecraft	
{

//protected:
public:
	void Aerodynamic(double max, double & s, double & s1, double & b) override;
private:
	double DCg[3]; // "���" � ������ ����������
	double _lSa; //����� ��
	double _hDivision; //������ ����������
	double _tDivsion;// ����� ����������
	double _vs;  //��������� �������� (�������� �� ����������� �/� ���������)
	unsigned _s; // ������� �� ����������� ����� (S = 0 ��������� ���� ������, S = 1 ��������� ���� �����)
	unsigned _vsr; // ������� �� ���������
	double _gamma, _dgamma, _gammaCom; //��������� ���� �����, ���������� ���� ����� � ������� �������� �� ����� (����������)
	double _pxMax; //������������ ����������
	bool _guidenceStart; //����� ������ ����������
	//������������ ��� ������������
	double ** GetTdCx();
	double ** GetTdCy();
	double ** GetTdMz();
	double ** GetTdCxa();
	double ** GetTdCya();
	double ** GetTdMza();
	double ** GetTCx2();
	double ** GetTCy2();
	double *** GetTal();	
	//�������� ������
	void DeleteMatr3(unsigned i, unsigned j, unsigned k, double *** matr);
	void DeleteMatr(unsigned nRow, unsigned nCol, double ** matr);

	double YtTma(); //���������� ayt
	void AimPointContDeorb(double longitudeEquator, double & fiAim, double & lAim); //������ ��������� ���������� ����� ���, ���-�� S � ������� �����. VSR
	void ShootGuideDeorbit(const double tVkl, const double aimLambda, vector<Spacecraft> &result) const; //�-� �������� ����� ������� ��� ����������� ������(��������) (tvkl=0 - ������ ��������)
	void PredictGuideLandingPoint(vector<Spacecraft> &result);//�������������� ��������� �������� ����� ��������� �������
	void RKS (double dt, Atmosphera &Atm, Soyuz (*f_ptr) (Soyuz & spCr, double t , Atmosphera & atm)); //�-� �������������� ������� ����� - �����	
	void RewriteVector(double * vec); //�������������� ������ �� �������
	void FindPxMax(double & currentPxmax, double px);

public:
	//������������
	Soyuz();
	Soyuz (const Vector &Vect, double _GKa, double _GBo, double _GSa, double _ImpSizeNu, double _Xt, double _Yt, double _Zt, double _Sbal,
		   double _RazpImpSize, double _NKa, unsigned _ROr, unsigned _RStab);
	Soyuz(const Spacecraft& Spcrafrt);
	~Soyuz();

	vector<Spacecraft> ControlledDeorb (const unsigned long _Vitok); //�-� �������� ������������ ������ ��
	void SetAltitudeOfDivision(double hDivision){_hDivision = hDivision;}
	void SetLenghtOfDeorbSpCraft(double lDeorbSpCr){_lSa = lDeorbSpCr;}

	typedef int (*func) (int);
	friend int func1(Spacecraft i);
	friend int func2(int i);
	int testRef(func);
	int TestRef();

	//���������
	friend Soyuz  operator + (const Soyuz &, const Soyuz &);
	friend Soyuz  operator - (const Soyuz &, const Soyuz &);
	friend Soyuz  operator * (          double  , const Soyuz &);
	friend Soyuz  operator * (const Soyuz &, double);

	//�������
	friend	Soyuz FGuideDeorb (Soyuz &, double dt, Atmosphera &Atm); // �-� �������� ������ ������ ��� ������������ ������
	friend void InterpolMatrix (double **AM, unsigned n, int i,  Soyuz & soyuz, unsigned &j, double arg, double arg1); //������� ��� ������������ (����� ����� �������� �������� �� �����)

};

