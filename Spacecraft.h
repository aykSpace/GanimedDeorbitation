#pragma once
#include "Vector.h"
#include "Engine.h"
#include <vector>

enum Orientation
{
	IKV = 1, MANUAL
};

enum Stabilization
{
	OSK = 1, ISK
};


class Spacecraft : public Vector
{
public:
	vector<CppFlightDynamicOperations::Engine> Engines; //������ ���������� ��
	double  GKa, GBo, GSa, ImpSizeCur, Xt, Yt, Zt, Sbal;
	double RazpImpSize, ImpSizeNu, NKa,PDu, dm, Mtek, TRabDu, TurnX, TurnY, TurnZ, pnX;
	unsigned ROr, RStab;
	bool GuidDeorb; // ����� ������������ ������
	bool EngineStart; // ����� ������ ��������� �������
	double MGSK[3][3]; //������� ����������
	static double P2[3][3];  //����������� ������� ����� ����������� ���������� (����� ��� ���� ����� �� ������������� ������ ��� ��� ����� � ��������� Orient)
public:
	double TVkl;
	Spacecraft();
	
	Spacecraft (const Vector &Vect, double _GKa, double _GBo, double _GSa, double _ImpSizeNu, double _Xt, double _Yt, double _Zt, double _Sbal,
				double _RazpImpSize, double _NKa, unsigned _ROr, unsigned _RStab);

	~Spacecraft();	

	//��������� ������
	Spacecraft & operator = (const Spacecraft &);
	friend Spacecraft  operator + (const Spacecraft &, const Spacecraft &);
	friend Spacecraft  operator - (const Spacecraft &, const Spacecraft &);
	friend Spacecraft  operator * (          double  , const Spacecraft &);
	friend Spacecraft  operator * (const Spacecraft &, double);

	
protected:
	explicit Spacecraft (const Vector &);

	//������ ������
	virtual void Aerodynamic(double max, double & s, double & s1, double & b);
	void Orient(double t); //�-� �������� ���������� � ������������ �� �� ����� ������ ��������
	virtual void RKS (double dt, Atmosphera &Atm) override; //�-� �������������� ������� ����� - �����	
	virtual Spacecraft F (Spacecraft &, double dt, Atmosphera &Atm); // �-� �������� ������ ������
	void BalSplashDownPoint (vector<Spacecraft> &Rezult); //�-� �������� ����� ������� ��� �������������� ������ ����� ���������� ��������
	void ShootSplDownPoint (const double _TVkl, const double AimLambda, vector<Spacecraft> &Rezult) const;//�-� ��������� �� ��������� ������� ���������
	double MorePreciseTVkl (const double AimLandLongitude,const Spacecraft &SCNU, const vector<Spacecraft> & Rezult, double *TTVkl, double *TLpos, unsigned NIter, double DLPos = 0.0166667); 
	//�-� ��������� ������� ��������� �� (ScNu - ������ �� �������� ����� ������, TTVkl - ������ ������ ��������� ��)
	void  RezultTimeCheck (vector <Spacecraft> & Rezult);//�������� ������������ 86400
	double TrackAngle () const; //������ ���� ������� ���������� 
	double AzimuthOfLanding (const Spacecraft &SCNU,  double TPos, double LPos); //������ ������� ������� � ����� �������
	double NewDTVkl (double *TTVkl, double *TLpos, double Aimlongitude, unsigned NIter); //������ ���� ����� 2-� ��������
	void RKSLowAtm (double dt, Atmosphera &Atm, Spacecraft (*f_ptr) (Spacecraft & spCr, double t , Atmosphera & atm)); //�-� �������������� ������� ����� - �����

public:
    vector<Spacecraft>  BalDeorb (const unsigned long _Vitok); //�-� �������� ��������������� ������ ��	
	unsigned CountOfEngines(){return _countOfEngines;}

	//friend �-�
	friend double AimLongitudeBS(double _lambda); //�-� ������ ���������� ������� ��� ��
	friend double FirstTimeStart (Spacecraft & SC, double LandLatitude); // �-� ���������� ������� ����������� ������� ���������
	friend Spacecraft EngineZone (const Spacecraft & SC, const double AimLambda, vector<Spacecraft> &Rezult, double _Tvkl = 0, unsigned numberOfCurrentEngine = 0); //�-� �������� ��������� �������
	friend void Mper (double V_A[6], double RezMatr[3][3]); // ������� ������� �������� ���->����
	friend void MultMatr(double A[3][3], double B[3][3], double C[3][3]); //��������� ������� 3�3
	friend void ED (double A[3][3]); //������� ��������� �������
	friend bool Razvor (double ug[3], double p[3][3]); //���� ����������� ����������
	friend void Mavg (double T, double T0, double p[3][3]); //��� ������������ ���� ���� �������� �����
	void OrientFisk(double t);

	//inline functions and properties
	void AddEngine(CppFlightDynamicOperations::Engine engine) //��������� ���������
	{
		Engines.push_back(engine);
		_countOfEngines = Engines.size();
	}; 

	void RemoveEngine(unsigned id);

	double GetSpacecraftNumber(){return NKa;}
	void SetSpacecraftNumber(double value){NKa = value;}

	double GetSpacecraftMass(){return GKa;}
	void SetSpacecraftMass(double value){GKa = value;}

	double GetBoMass(){return GBo;}
	void SetBoMass(double value){GBo = value;}

	double GetDeorbitModuleMass(){return GBo;}
	void SetDeorbitModuleMass(double value){GBo = value;}

	double GetSbal(){return Sbal;}
	void SetSBal(double value){Sbal = value;}

	double GetXt(){return Xt;}
	void SetXt(double value){Xt = value;}

	double GetYt(){return Yt;}
	void SetYt(double value){Yt = value;}

	double GetZt(){return Zt;}
	void SetZt(double value){Zt = value;}

	double GetDeorbitImpulse(){return ImpSizeNu;}
	void SetDeorbitImpulse(double value){ImpSizeNu = value;}

	double GetRaspImpulse(){return RazpImpSize;}
	void SetRaspImpulse(double value){RazpImpSize = value;}

	unsigned GetRateOfOrient(){return ROr;}
	void SetRateOfOrient(unsigned value){ROr = value;}

	unsigned GetRateOfStab(){return RStab;}
	void SetRateOfStab(unsigned value){RStab = value;}

	void RemoveEngine(CppFlightDynamicOperations::Engine engine);
	friend Spacecraft FLowAtm (Spacecraft &, double dt, Atmosphera &Atm); // �-� �������� ������ ������ ��� ������ ���������

private:
	unsigned _countOfEngines;
	void CalcOptOrient (double a[3], double RezMatr[3][3]); // ������� ������� �������� ���->����

};

