#include "Spacecraft.h"
#include "BALCONST.H"
#include <algorithm>

#define RAD_ENTRANCE 6471

//статические объекты
double Spacecraft::P2[3][3]; //выделяю память под статическую матрицу 

// конструкторы и деструкторы
Spacecraft::Spacecraft()
{
	GKa = GBo = GSa = ImpSizeNu = Xt = Yt = Zt = Sbal = 0;
	RazpImpSize = NKa = Mtek = TRabDu = TVkl = 0;
	ROr = 1; RStab = 2;
	ImpSizeCur = PDu = dm = TurnX = TurnY = TurnZ =  0;
	GuidDeorb = false;
	EngineStart = false;
	_countOfEngines = 0;
	pnX = 0;
}

Spacecraft::Spacecraft(const Vector &Vect) :Vector(Vect), GKa(0), GBo(0), GSa(0), ImpSizeNu(0), Xt(0), Yt(0),
	Zt(0), Sbal(0), RazpImpSize(0), NKa(0), ROr(0), RStab(0), _countOfEngines(0), EngineStart(false), pnX(0),
	TVkl(0), TurnX(0), TurnY(0), TurnZ(0), TRabDu(0),dm(0), PDu(0), Mtek(0), ImpSizeCur(0) 
{}

Spacecraft::Spacecraft (const Vector &Vect, double _GKa, double _GBo, double _GSa, double _ImpSizeNu, double _Xt, double _Yt, double _Zt, double _Sbal,
						double _RazpImpSize, double _NKa, unsigned _ROr, unsigned _RStab)
					  : Vector(Vect), GKa(_GKa), GBo(_GBo), GSa(_GSa), ImpSizeNu(_ImpSizeNu), Xt(_Xt), Yt(_Yt),
					    Zt(_Zt), Sbal(_Sbal), RazpImpSize(_RazpImpSize), NKa(_NKa), ROr(_ROr), RStab(_RStab), _countOfEngines(0), EngineStart(false), pnX(0),
						TVkl(0), TurnX(0), TurnY(0), TurnZ(0), TRabDu(0), dm(0), PDu(0), Mtek(0), ImpSizeCur(0)
{}

Spacecraft::~Spacecraft()
{}

//операторы

Spacecraft & Spacecraft::operator = (const Spacecraft &SC)
{
	x = SC.x; y = SC.y; z = SC.z; Vx = SC.Vx; Vy = SC.Vy; Vz = SC.Vz; Sb = SC.Sb; vitok = SC.vitok; t = SC.t; lambda = SC.lambda; fi = SC.fi; h = SC.h; ro = SC.ro; r = SC.r;
	ImpSizeCur = SC.ImpSizeCur;  pnX = SC.pnX; DatTime = SC.DatTime; Engines = SC.Engines; _ax = SC._ax; _ay = SC._ay; _az = SC._az; //Vector::operator=(SC);
	return *this;
}


Spacecraft  operator + (const Spacecraft & a, const Spacecraft & b)
{
	Spacecraft c;
	c.x      = a.x  + b.x;
	c.y      = a.y  + b.y;
	c.z      = a.z  + b.z;
	c.Vx     = a.Vx + b.Vx;
	c.Vy     = a.Vy + b.Vy;
	c.Vz     = a.Vz + b.Vz;
	c.ImpSizeCur = a.ImpSizeCur + b.ImpSizeCur;
	c.pnX     = b.pnX;
	c.Sb     = a.Sb;
	c.vitok  = a.vitok;
	c.t      = a.t;
	c.fi     = a.fi;
	c.lambda = a.lambda;
	c.h      = a.h;
	c.ro     = a.ro;
	c.r      = a.r;
	c.ImpSizeCur = a.ImpSizeCur + b.ImpSizeCur;
	c.EngineStart = a.EngineStart;
	c.Mtek    = a.Mtek;
	c.PDu     = a.PDu;
	c.TVkl    = a.TVkl;
	c.ROr	  = a.ROr;
	c.RStab	  = a.RStab;
	c.Xt      = a.Xt;
	c.Yt      = a.Yt;
	c.Zt      = a.Zt;
	c.GSa     = a.GSa;
	c.DatTime = a.DatTime;
	c.Engines = a.Engines;
	c._ax     = b._ax;
	c._ay     = b._ay;
	c._az     = b._az;
	return(c);
}

Spacecraft  operator - (const Spacecraft & a, const Spacecraft & b)
{
	Spacecraft c;
	c.x  =     a.x  - b.x;
	c.y  =     a.y  - b.y;
	c.z  =     a.z  - b.z;
	c.Vx =     a.Vx - b.Vx;
	c.Vy =     a.Vy - b.Vy;
	c.Vz =     a.Vz - b.Vz;
	c.ImpSizeCur = a.ImpSizeCur - b.ImpSizeCur;
	c.pnX     = a.pnX;
	c.Sb =     a.Sb;
	c.vitok =  a.vitok;
	c.t =      a.t;
	c.fi =     a.fi;
	c.lambda = a.lambda;
	c.h =      a.h;
	c.r =      a.r;
	c.ro =     a.ro;
	c.ImpSizeCur = a.ImpSizeCur - b.ImpSizeCur;
	c.EngineStart = a.EngineStart;
	c.Mtek    = a.Mtek;
	c.PDu     = a.PDu;
	c.TVkl    = a.TVkl;
	c.ROr	  = a.ROr;
	c.RStab	  = a.RStab;
	c.DatTime = a.DatTime;
	c.Engines = a.Engines;
	c._ax    = a._ax;
	c._ay    = a._ay;
	c._az    = a._az;
	return(c);
}

Spacecraft operator * (double a, const Spacecraft & b)
{
	Spacecraft c;
	c.ImpSizeCur = a * b.ImpSizeCur;
	c.x      = a * b.x;
	c.y      = a * b.y;
	c.z      = a * b.z;
	c.Vx     = a * b.Vx;
	c.Vy     = a * b.Vy;
	c.Vz     = a * b.Vz;
	c.Sb     =     b.Sb;
	c.vitok  =     b.vitok;
	c.t      =     b.t;
	c.fi     =     b.fi;
	c.lambda =     b.lambda;
	c.h      =     b.h;
	c.r      =     b.r;
	c.ro     =     b.ro;
	c.ImpSizeCur = a * b.ImpSizeCur;
	c.pnX     = b.pnX;
	c.EngineStart = b.EngineStart;
	c.Mtek    = b.Mtek;
	c.PDu     = b.PDu;
	c.TVkl    = b.TVkl;
	c.ROr	  = b.ROr;
	c.RStab	  = b.RStab;
	c.DatTime = b.DatTime;
	c.Engines = b.Engines;
	c._ax     = b._ax;
	c._ay     = b._ay;
	c._az     = b._az;
	return(c);
}


Spacecraft  operator * (const Spacecraft & b, double a)
{
	Spacecraft c;
	c.x      = a * b.x;
	c.y      = a * b.y;
	c.z      = a * b.z;
	c.Vx     = a * b.Vx;
	c.Vy     = a * b.Vy;
	c.Vz     = a * b.Vz;
	c.ImpSizeCur = a* b.ImpSizeCur;
	c.Sb     =     b.Sb;
	c.vitok  =     b.vitok;
	c.t      =     b.t;
	c.fi     =     b.fi;
	c.lambda =     b.lambda;
	c.h      =     b.h;
	c.r      =     b.r;
	c.ro     =     b.ro;
	c.ImpSizeCur = a* b.ImpSizeCur;
	c.EngineStart = b.EngineStart;
	c.Mtek    = b.Mtek;
	c.PDu     = b.PDu;
	c.TVkl    = b.TVkl;
	c.ROr	  = b.ROr;
	c.RStab	  = b.RStab;
	c.pnX     = b.pnX;
	c.DatTime = b.DatTime;
	c.Engines = b.Engines;
	c._ax     = b._ax;
	c._ay     = b._ay;
	c._az     = b._az;
	return(c);
}

void MultMatr(double A[3][3], double B[3][3], double C[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j =0; j < 3; j++)
		{
			C[i][j] = 0;
			for (int k = 0 ; k < 3; k++)
				C[i][j] = C[i][j] + A[i][k]*B[k][j];
		}
	}
}

void ED (double A[3][3]) 
{
	for (int i = 0; i < 3; i++) 
	{
		for (int j = 0; j < 3; j++)
		{
			A[i][j] = 0;
			if (i ==j)
				A[i][j] = 1;
		}
	}
}

bool Razvor (double ug[3], double p[3][3])
{ 
	//Kv - порядок разворотов
	//ug[3] - углы разворотов (рысканье, крен, тангаж)
	int K[3] = {0,0,0};
	double A0[3][3], A1[3][3], A2[3][3];
	double B[3][3];
	ED(p);
	if (ug[0] != 0)
		K[0] = 1;
	if (ug[1] != 0)
		K[1] = 1;
	if (ug[2] != 0)
	{
		K[2] = 1;
		ED(A2);
		double s = sin(ug[2]);
		double c = cos(ug[2]);
		A2[0][0] = c;
		A2[0][1] = -s;
		A2[1][0] = s;
		A2[1][1] = c;
		MultMatr(A2, p, B);
		for (int i = 0; i < 3; i++)
			for (int j = i; j < 3; j++)
				p[i][j] = B[i][j];
		return true;
	}
	if (K[0] == K[1] && K[0] == K[2] && K[0] == 0)
		return false;
	return true;
}

void Mavg (double T, double T0, double p[3][3])
{
	for (int i = 0; i < 2; i++) 
	{
		p[2][i] = 0;
		p[i][2] = 0;
	}
	p[2][2] = 1;
	double S  = OMEGA_EARTH*(T-T0);
	double SS = sin(S);
	double CS = cos(S);
	p[0][0] = CS;
	p[0][1] = SS;
	p[1][0] = -SS;
	p[1][1] = CS;
}

void Mper (double V_A[6], double RezMatr[3][3])
{
	double C[3];
	double R      = sqrt(V_A[0]*V_A[0] + V_A[1]*V_A[1] + V_A[2]*V_A[2]);
	double RM     = 1/R;
	RezMatr[0][1] = V_A[0]*RM;
	RezMatr[1][1] = V_A[1]*RM;
	RezMatr[2][1] = V_A[2]*RM;
	C[0] = V_A[4]*V_A[2] - V_A[5]*V_A[1];
	C[1] = V_A[5]*V_A[0] - V_A[3]*V_A[2];
	C[2] = V_A[3]*V_A[1] - V_A[4]*V_A[0];
	double A = 1/sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);
	for (unsigned i = 0; i < 3; i++)
		RezMatr[i][2] = C[i]*A;
	RezMatr[0][0] = RezMatr[1][1]*RezMatr[2][2] - RezMatr[2][1]*RezMatr[1][2];
	RezMatr[1][0] = RezMatr[2][1]*RezMatr[0][2] - RezMatr[0][1]*RezMatr[2][2];
	RezMatr[2][0] = RezMatr[0][1]*RezMatr[1][2] - RezMatr[1][1]*RezMatr[0][2];
}

void Spacecraft::Orient(double t)
{
	double V_A[6];         //массив(вектор) в инерциальной СК
	double P1[3][3];      //матрица перехода ОСК - ИСКТ
	double P3[3][3];     //матрица учета вращения Земли для перевода из ИСКФ в ГСК
	double Prazv[3][3]; //матрица  программных разворотов
	double UgRazv[3] = {TurnX, TurnY, TurnZ};
	if (ROr == 1 && RStab == 2 && fabs(t - TVkl) <= 0.01) //первый вход для ориентации ИКВ и стабилизации ИСКТ
	{
		V_A[0] = x;
		V_A[1] = y;
		V_A[2] = z;
		V_A[3] = Vx - OMEGA_EARTH*y;
		V_A[4] = Vy + OMEGA_EARTH*x;
		V_A[5] = Vz;
		Mper(V_A,P1);
		Razvor(UgRazv,Prazv);	
		MultMatr(P1,Prazv, P2);
		Mavg(t,TVkl,P3);
		MultMatr(P3,P2,MGSK);
		return;
	}
	if (ROr == 1 && RStab == 2 && fabs(t - TVkl) > 0.01)
	{
		Mavg(t,TVkl,P3);
		MultMatr(P3,P2,MGSK);
		return;
	}
}

void Spacecraft::RKS(double dt, Atmosphera &Atm)
{
	Spacecraft K1, K2, K3, K4;
	K1 = F(*this, t, Atm) * dt;
	K2 = F(*this + 0.5*K1, t + 0.5*dt, Atm) * dt;
	K3 = F(*this + 0.5*K2, t + 0.5*dt, Atm) * dt;
	K4 = F(*this +     K3, t + dt, Atm) * dt;
	*this = *this + 1.0/6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
	t += dt;
}

void Spacecraft::RKSLowAtm(double dt, Atmosphera &Atm, Spacecraft (*f_ptr) (Spacecraft & spCr, double t , Atmosphera & atm))
{
	Spacecraft K1, K2, K3, K4;
	K1 = f_ptr(*this, t, Atm) * dt;
	K2 = f_ptr(*this + 0.5*K1, t + 0.5*dt, Atm) * dt;
	K3 = f_ptr(*this + 0.5*K2, t + 0.5*dt, Atm) * dt;
	K4 = f_ptr(*this +     K3, t + dt, Atm) * dt;
	*this = *this + 1.0/6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
	t += dt;
}

Spacecraft Spacecraft::F (Spacecraft & SC, double t, Atmosphera &Atm)
{
	Spacecraft Rez;
	double Koef1, Koef2, Koef3;
	GetDensity(Vector::f,Vector::ap,Atm,SC); //рассчет плотности атмосферы (потом надо убрать F and Ap)
	double VB[3] = {SC.Vx - Atm.W[0], SC.Vy - Atm.W[1], SC.Vz - Atm.W[2]};
	double V = sqrt(VB[0]*VB[0] + VB[1]*VB[1] + VB[2]*VB[2]);
	Koef1 = OMEGA_EARTH * OMEGA_EARTH;
	Koef2 = OMEGA_EARTH * 2.0;
	Koef3 = -SC.Sb * SC.ro * V * 1000;
	double *RezGpz = new double [3];
	gpz(SC, RezGpz); // рассчет коэффициентов гравитационного поля Земли
	double AxDu = 0;//компоненты ускорения от двигателя
	double AyDu = 0;
	double AzDu = 0;
	double AT = 0; //ускорение от двигателя в Ньютонах
	Rez.pnX = -Koef3*V*101.9716213;
	if (SC.EngineStart)
	{
		Rez.pnX = SC.PDu/SC.Mtek; //перегрузка
		AT  = Rez.pnX*0.00980665;
		SC.Orient(t);
		Rez.ImpSizeCur = AT;
		AxDu = MGSK[0][0]*AT;
		AyDu = MGSK[1][0]*AT;
		AzDu = MGSK[2][0]*AT;
	}
	Rez.x = SC.Vx;
	Rez.y = SC.Vy;
	Rez.z = SC.Vz;
	Rez.Vx = Koef3 * VB[0] + Koef2*VB[1] + Koef1*SC.x + RezGpz[0] - AxDu;
	Rez.Vy = Koef3 * VB[1] - Koef2*VB[0] + Koef1*SC.y + RezGpz[1] - AyDu;
	Rez.Vz = Koef3 * VB[2] + RezGpz[2] - AzDu;
	delete [] RezGpz;
	return Rez;
}

Spacecraft FLowAtm(Spacecraft &SC, double dt, Atmosphera &Atm)
{
	const double Rzsr  = 6371.0;
	const double TAI20 = -0.10183014444;
	const double omz2  = 53.174954e-10; 
	const double A100  = 62.5648249885;
	const double domz  = 14.58423171e-5;
	const double Rekv  = 6378.14;
	const double szat  = 0.003352824419;
	Spacecraft Rez;
	Rez.x = SC.Vx;
	Rez.y = SC.Vy;
	Rez.z = SC.Vz;
	double R12 = SC.x*SC.x + SC.y*SC.y;
	double Z2  = SC.z*SC.z;
	double R2  = R12 + Z2;
	double R   = sqrt(R2);
	double RM1 = 1/R;
	double RM2 = RM1*RM1;
	double C6  = Rzsr*RM1;
	double C7  = C6*RM2;
	C6 = TAI20*C6*C7;
	double C8 = 5*Z2;
	double C1 = omz2 - C7*(A100 - C6*(R2 - C8));
	Rez.Vx =  domz*SC.Vy + SC.x*C1;
	Rez.Vy = -domz*SC.Vx + SC.y*C1;
	Rez.Vz = -SC.z*C7*(A100 + C6*(C8 - 3 *R2));
	double Rz = Rekv*(1.0 - szat*Z2*RM2);
	SC.h = R - Rz;
	GetDensity(110,16,Atm,SC); //рассчет плотности атмосферы (потом надо убрать F and Ap)
	double VB[3] = {SC.Vx - Atm.W[0], SC.Vy - Atm.W[1], SC.Vz - Atm.W[2]};
	double VB2 = VB[0]*VB[0] + VB[1]*VB[1] + VB[2]*VB[2];
	double VB1 = sqrt(VB2);
	double S   = SC.Sbal*1000;
	double C0  = SC.ro*VB1;
	double CS  = -S*C0;
	Rez.Vx = Rez.Vx + CS*VB[0];
	Rez.Vy = Rez.Vy + CS*VB[1];
	Rez.Vz = Rez.Vz + CS*VB[2];
	return Rez;
}

double AimLongitudeBS(double _Lambda)
{
	double LPr;
	double TleKB[25] = {-180, -160, -133, -128, -123, -110.8, -105, -101, -98, -77, -55, -40, -34, -28, -3, 3.5, 17, 
						 30, 90, 104, 109, 115, 130,150, 180};
	double l0[25]    = {-103, -74.4, -134.2, -89, 0.5, 19, 3, 18, 50, 44, 63, 67.33333,
		                 64, 67.33333, 63, 69, 64, 138, 100.4, -171.8, -214, -123, -65.8,-103, 0};
	double AK[25]    = {0, 0.17, -0.28, 0, 0, 0.167, 0, 0, 0.233, 0, 0, 0, 0, 0, 0, 0,
		                0, 0, 0.357, 0.55, 1, 0.2, -0.24, 0, 0};
	int i = poi(25, TleKB, _Lambda);
	LPr = l0[i] + AK[i]*_Lambda;
	return LPr;
}

double FirstTimeStart (Spacecraft & SC, double LandLatitude)
{
	double fA; //внеатмосферная дальность
	const double pPg = 57.29578; //перевод из рад в градусы
	double DV1    = SC.ImpSizeNu/1000;
	double *Vek_a = new double[6];
	Vek_a[0]      = SC.x;
	Vek_a[1]      = SC.y;
	Vek_a[2]      = SC.z;
	Vek_a[3]      = SC.Vx - OMEGA_EARTH*SC.y;
	Vek_a[4]      = SC.Vy + OMEGA_EARTH*SC.x;
	Vek_a[5]      = SC.Vz;
	//double Rnu    = sqrt(SC.x*SC.x + SC.y*SC.y + SC.z*SC.z);
	double VA     = sqrt(Vek_a[3]*Vek_a[3] + Vek_a[4]*Vek_a[4] + Vek_a[5]*Vek_a[5]);
	double TETHy  = asin((Vek_a[0]*Vek_a[3] + Vek_a[1]*Vek_a[4] + Vek_a[2]*Vek_a[5])/(SC.r*VA));
	double VA1    = VA - DV1; //учет отработанного импульса
	double TgTET  = tan(TETHy); //угловая дальность до входа по Кеплеру 
	double Nu0    = SC.r*VA1*VA1/MU;
	double a      = 2*RAD_ENTRANCE*(1 + TgTET*TgTET) - (SC.r + RAD_ENTRANCE)*Nu0;
	double b      = Nu0*RAD_ENTRANCE*TgTET;
	double c      = (SC.r - RAD_ENTRANCE)*Nu0;
	double D1     = sqrt(b*b + a*c);
	double AA     = atan((b +D1)/a);
	double f      = 2*AA;
	double *C     = new double[3];
	C[0]          = SC.y*Vek_a[5] - SC.z*Vek_a[4];
	C[1]          = SC.z*Vek_a[3] - SC.x*Vek_a[5];
	C[2]          = SC.x*Vek_a[4] - SC.y*Vek_a[3];
	double cosi   = C[2]/sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]); //косинус угла наклона орбиты
	SC.GuidDeorb  == true ? fA = 0.54 : fA = 0.45;
	double fp     = f + fA;
	//разница долгот НУ и посадки
	double Tpol   = SC.r*f/(VA*cos(TETHy)) + 450; //время полета(450с - время полета на атмосферы уч-ке в первом прибл)
	double DL     = LandLatitude - SC.lambda; //приведенная долгота
	double fNu = 0;
	if (DL < 0)
		DL += 360;
	DL = fabs(DL + OMEGA_EARTH*Tpol*pPg);
	double D  = DL;
	double KL = 0;
	if (DL > 180)
	{
		D = D - 180;
		KL = 1;
	}
	if (DL < 180)
	{
		if (fabs(D - 90) < 2)
			fNu = 1.570796;
		else if (D < 90)
			fNu = atan(tan(D*0.0174553)/cosi);
		else 
			fNu = PI - atan(tan((180 - D)*0.017455)/cosi);
	}
	if (KL == 1)
		fNu += PI;
	double Df = fNu - fp;
	double TStart0 = SC.t + Df*SC.r/(VA*cos(TETHy));
	delete[] C;
	delete[] Vek_a;
	return TStart0;
}


Spacecraft EngineZone(const Spacecraft &_SC, const double AimLambda, vector<Spacecraft> &Rezult, double _TVkl, unsigned curEngine)
{
	Spacecraft SC (_SC);
	if (SC.CountOfEngines() == 0)
		throw logic_error("Spacecraft has`t engines");
	double dt = 0;
	if (_TVkl == 0)
		SC.TVkl = FirstTimeStart(SC,AimLambda);
	else
		SC.TVkl = _TVkl;
	Atmosphera Atm;
	Atm.NMounth = SC.DatTime.Get_Mounth() - 1;
	int i = 0; 
	unsigned j = 1; //для корректного первого шага в InterpolMatrix (i = 1, вычисляем AM[j-1][3])
	unsigned OrderMatr = 8;
	double **AM;
	AM = new double * [OrderMatr];
	for (unsigned k = 0; k < OrderMatr; ++k)
		AM[k] = new double [OrderMatr];
	SC.Perevod();
//-----------------------------интегрирование до времени включения двигателя------------------------------------------
	if (SC.t > SC.TVkl)  
	{
		while (SC.t > SC.TVkl)
		{
			dt = -20;
			Atm.StarTime(SC.DatTime);// можно пересчитывать при изменении даты!!!
			Atm.Sun(SC.t);
			SC.RKS(dt, Atm);
			i += 1;
			InterpolMatrix(AM, OrderMatr, i, SC, j, 0);
		}
		if (i >= 8)
			Ilagr(8, AM, 0, SC.TVkl, SC, 7);
		else
		{
			while (i < 8)
			{
				Atm.StarTime(SC.DatTime);
				Atm.Sun(SC.t);
				SC.RKS(dt, Atm);
				i += 1;
				InterpolMatrix(AM, OrderMatr, i, SC, j, 0);
			}
			Ilagr(8, AM, 0, SC.TVkl, SC, 7);
		}
	}
	else
	{
		while (SC.t < SC.TVkl)
		{
			dt = 20;
			Atm.StarTime(SC.DatTime); 
			Atm.Sun(SC.t);
			SC.RKS(dt, Atm);
			i += 1;
			InterpolMatrix(AM, OrderMatr, i, SC, j, 0);
		}
		if (i >= 8)
			Ilagr(8, AM, 0, SC.TVkl, SC, 7);
		else
		{
			while (i < 8)
			{
				Atm.StarTime(SC.DatTime);
				Atm.Sun(SC.t);
				SC.RKS(dt, Atm);
				i += 1;
				InterpolMatrix(AM, OrderMatr, i, SC, j, 0);
			}
			Ilagr(8, AM, 0, SC.TVkl, SC, 7);
		}
	}
	SC.Perevod();
	SC.DatTime.SetFrmtTime(SC.t);
	Rezult.push_back(SC);
//-----------------------------интегрирование активного участка------------------------------------------
	/*Vector test (1819, 84067.7205676898, 2544.06821382599, -3777.51194633820, -4993.13380467807, 6.82645770391152, 1.93692220579422, 2.01518873552143, 0.04158185);
	SC = test; //для отладки
	SC.Perevod();
	SC.TVkl = SC.t; //потом убрать!!*/
	/*Vector::f = 110;
	Vector::ap = 16;
	Date_Time DataTime(2011,10,05);
	Vector test (1819, 84031.4821123286, 2294.79587052663, -3845.15907500935, -5061.90586124525, 6.92907024568318, 1.79609255443045, 1.77980917745032, 0.04158185,DataTime);
	SC = test; //для отладки
	SC.Perevod();
	SC.TVkl = SC.t; //потом убрать!! (aus)*/

	SC.EngineStart = true;
	i = 0; j = 1; //сбрасываем счетчики
	//SC.dm = SC.PSkd/SC.PUdSkd; //может быть ДПО!!!
	SC.dm = SC.Engines[curEngine].Trust/SC.Engines[curEngine].SpecificImpulse;
	//SC.PDu = SC.PSkd;//может быть ДПО!!!
	SC.PDu = SC.Engines[curEngine].Trust;
	if (!SC.GuidDeorb)
		SC.GKa = SC.GKa - SC.GBo; //для БС вычитаем массу БО
	dt = 2;
	SC.ImpSizeCur = 0;
	while (SC.ImpSizeCur < SC.ImpSizeNu/1000)
	{
		SC.TRabDu = SC.t - SC.TVkl;
		SC.Mtek   = SC.GKa - SC.dm*SC.TRabDu;
		Atm.StarTime(SC.DatTime);
		Atm.Sun(SC.t);
		SC.RKS(dt, Atm);
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, SC, j, SC.ImpSizeCur);
	}
	Ilagr(8, AM, 7, SC.ImpSizeNu/1000, SC, 7);
	SC.EngineStart = false;
	SC.Perevod();
	SC.DatTime.SetFrmtTime(SC.t);
	Rezult.push_back(SC);
	for (unsigned i = 0; i < OrderMatr; ++i) // удаление матрицы
		delete [] AM [i];
	delete [] AM;
	return SC;
}

void Spacecraft::RezultTimeCheck(vector<Spacecraft> & Rezult)
{
	for (unsigned i = 0; i < Rezult.size(); i++)
	{
		if (Rezult[i].TVkl > 86400)
			Rezult[i].TVkl -= 86400;
	}
}

void Spacecraft::BalSplashDownPoint(vector<Spacecraft> &Rezult)
{
	/*Vector test (1819, 84282.3583662208, 3916.06388296012, -3279.94113335670, -4420.17274327533, 5.89755525155952, 2.68185228922996, 3.29740858594304, 0.04158185);
	*this = test; //для отладки
	this->Perevod();
	this->TVkl = 84067.7205676898; */
	double TRazd = this->TVkl + 725.5 + 540; //540 - длинная метка разделения	
	Atmosphera Atm;
	Atm.NMounth = this->DatTime.Get_Mounth() - 1;
	int i = 0; 
	unsigned j = 1; //для корректного первого шага в InterpolMatrix (i = 1, вычисляем AM[j-1][3])
	unsigned OrderMatr = 8;
	double **AM;
	AM = new double * [OrderMatr];
	for (unsigned k = 0; k < OrderMatr; ++k)
		AM[k] = new double [OrderMatr];
	this->Perevod();
	double dt = 2;
	while (this->t <= TRazd)
	{
		Atm.StarTime(this->DatTime);// можно пересчитывать при изменении даты!!!
		Atm.Sun(this->t);
		this->RKS(dt, Atm);
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, *this, j, 0);
		this->Altitude();
		if (this->h <= 105.0 )
			throw logic_error("h < 105 km but no devision");
	}
	Ilagr(OrderMatr, AM, 0, TRazd, *this, 7); //интерполяция на время разделения	
	//--------------------------------------------------------------------------------------------------------------------------------------------------------
	/*Vector test (1819, 85333.2205676898, 6465.79278060700, 739.568085783398, 1096.09080325654, -1.71174174286032, 4.23680077886240, 5.90805213351881, 0.04158185);
	*this = test; //для отладки
	this->Perevod();*/
	this->Sb = this->Sbal;
	this->Perevod();
	this->DatTime.SetFrmtTime(this->t);
	Rezult.push_back(*this);
	while (this->r >= 6471.0)
	{
		if (this->h > 120)
		{
			Atm.StarTime(this->DatTime);// можно пересчитывать при изменении даты!!!
			Atm.Sun(this->t);
		}		
		this->RKS(dt, Atm);
		i += 1;
		this->Perevod();
		InterpolMatrix(AM, OrderMatr, i, *this, j, this->r);
	}
	Ilagr(OrderMatr, AM, 7, 6471.0, *this, 7); //интерполяция на начало атмосферного участка
	/*Vector test (1819, 85924.2970546186, 4103.91770907551, 2947.07302663006, 4043.07568725117, -6.00746299372389, 2.96782020056082, 3.62219157033226, 0.04158185);
	*this = test; //для отладки
	this->Sb = this->Sbal;*/
	this->Perevod();
	this->r = 6471.0;
	this->DatTime.SetFrmtTime(this->t);
	Rezult.push_back(*this);
	dt = 0.6;
	while (this->h >= 10.7)
	{
		this->RKS(dt,Atm);
		this->Altitude();
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, *this, j, this->h);
	}
	Ilagr(OrderMatr, AM, 7, 10.7e0, *this, 8); //интерполяция на конец атмосферного участка
	this->Perevod();
	this->DatTime.SetFrmtTime(this->t);
	Rezult.push_back(*this);
	RezultTimeCheck(Rezult);
	for (unsigned i = 0; i < OrderMatr; ++i) // удаление матрицы
		delete [] AM [i];
	delete [] AM;

}

double  Spacecraft::TrackAngle() const 
{
	double V   = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
	double Tet = asin((x*Vx + y*Vy + z*Vz)/(r*V));
	return Tet;
}

double Spacecraft::NewDTVkl(double *TTVkl, double *TLpos, double Aimlongitude, unsigned NIter)
{
	if (NIter >= 10)
		throw logic_error("number of iteration greater then 10");
	if (Aimlongitude < 0)
		Aimlongitude += 360;
	for (unsigned i = 0; i < NIter; i++)
	{
		if (fabs(TLpos[i] - Aimlongitude) > 180)
			TLpos[i] += 360;
	}
	double DTvkl = TTVkl[NIter - 1] - TTVkl[NIter];  //посмотреть niter!!!!!
	double DLpos = TLpos[NIter - 1] - TLpos[NIter];
	if (Sign(DTvkl) != Sign(DLpos) && fabs(DLpos) < 180)
		return  DTvkl;
	//double D1       = DTvkl;
	double DTvklRez = DTvkl*(Aimlongitude - TLpos[NIter])/DLpos;
	if (fabs(DTvkl - DTvklRez) < 0.005)
		DTvklRez *= 0.5;
	if (DTvklRez > -2700)
		return  DTvklRez;
	throw logic_error("invalid dt of engine starting");
}

double Spacecraft::AzimuthOfLanding(const Spacecraft &SCNU, double TPos, double LPos)
{
	double DT = TPos - SCNU.t;
	if (DT < 0)
		DT += 86400;
	double LP = LPos*PerGradRad + DT*OMEGA_EARTH;
	double Ex = -sin(LP);
	double Ey =  cos(LP);
	double H1 =  SCNU.y*SCNU.Vz;
	double H2 = -SCNU.x*SCNU.Vz;
	double H3 =  SCNU.x*SCNU.Vy - SCNU.y*SCNU.Vx;
	double C = 1/(sqrt(H1*H1 + H2*H2 + H3*H3));
	double Azimuth = asin((Ex*H1 + Ey*H2)*C);
	Azimuth = 90 + Azimuth/PerGradRad;
	return Azimuth;
}

double Spacecraft::MorePreciseTVkl(const double AimLandLongitude,const Spacecraft &SC_NU, const vector<Spacecraft> & Rezult, double *TTVkl, double *TLpos, unsigned NIter, double DLPos)
{
	NIter -= 1;
	double DT;
	if (fabs(TLpos[NIter] - AimLandLongitude) > DLPos)
	{
		if (NIter < 1)
		{
			double Tet  = Rezult.front().TrackAngle();
			double VVkl = sqrt((Rezult.front().Vx - OMEGA_EARTH*Rezult.front().y)*(Rezult.front().Vx - OMEGA_EARTH*Rezult.front().y) + 
				               (Rezult.front().Vy + OMEGA_EARTH*Rezult.front().x)*(Rezult.front().Vy + OMEGA_EARTH*Rezult.front().x) +
							   (Rezult.front().Vz*(Rezult.front().Vz)));
			double CC   = (Rezult.back().x*Rezult.back().x + Rezult.back().y*Rezult.back().y);
			double CC1  = 1.0/(CC + Rezult.back().z*Rezult.back().z);
			double Sif2 = Rezult.back().z*Rezult.back().z*CC1;
			double Cof2 = CC*CC1;
			double DLp = TLpos[NIter] - AimLandLongitude;
			if (DLp > 180)
				DLp -= 360;
			if (DLp < -180)
				DLp +=360;
			double All = Sign(DLp)*(90.0*PerGradRad - asin(Cof2*cos(PerGradRad*DLp) + Sif2));
			double Azimuth = AzimuthOfLanding(SC_NU,Rezult.back().t,TLpos[NIter]);
			Azimuth = 90.0-Azimuth;
			double DDL = -Rezult.front().r/(VVkl*cos(Tet)*cos(Azimuth*PerGradRad));
			double DTB1 = All*DDL;
			double DTB2 = 0;
			unsigned Iter = 0;
			do
			{
				double DELI = -OMEGA_EARTH*(DTB1-DTB2)*PerRadGrad;
				DLp += DELI;
				All  = Sign(DLp)*(90.0*PerGradRad - asin(Cof2*cos(PerGradRad*DLp) + Sif2));
				DTB2 = DTB1;
				DTB1 = DDL*All;
				Iter += 1;
				if (Iter >=10)
					break;
			}
			while (fabs(DTB2 - DTB1) > 0.001);
			DT = DTB1;
		}
		else
		{
			DT = NewDTVkl(TTVkl, TLpos, AimLandLongitude,NIter);
		}
	}
	else 
		return 0;
	return DT;
}

void Spacecraft::ShootSplDownPoint(const double _TVkl, const double AimLambda, vector<Spacecraft> &Rezult) const
{
	Spacecraft SpCr = EngineZone(*this, AimLambda, Rezult, _TVkl);
	SpCr.BalSplashDownPoint(Rezult);
}

vector<Spacecraft> Spacecraft::BalDeorb(const unsigned long _Vitok)
{
	//проверка _виток > this.vitok!!
	
	this->GuidDeorb = false;
	unsigned NIter = 0;//номер итерации
	vector<Spacecraft> RezultOfCalculations; //(0 - SC на момент вкл ДУ, 1 - SC после импульса, 2 - SC на момент разделения, 3 - SC на начало атм у-ка 4 - SC на h=10.7
	double *TLPos = new double[10]; //10 ограничение по итерацям
	double *TTVkl = new double[10];
	double DTVkl;
	vector<Vector> arxOfPrognoz = this->Prognoz(static_cast<unsigned int>(_Vitok));
	double AimLandLongitude = AimLongitudeBS(this->lambda);
	this->ShootSplDownPoint(0,AimLandLongitude,RezultOfCalculations); // первый прострел
	RezultOfCalculations.front().TVkl > 0 ? TTVkl[NIter] = RezultOfCalculations.front().TVkl : TTVkl[NIter] = RezultOfCalculations.front().TVkl + 86400;
	TLPos[NIter] = RezultOfCalculations.back().lambda;
	NIter += 1;
	do 
	{
		DTVkl = this->MorePreciseTVkl(AimLandLongitude,*this,RezultOfCalculations,TTVkl,TLPos,NIter); //уточнение времени включения
		if (fabs(DTVkl) <= 0.0001)
			break;
		TTVkl[NIter] = TTVkl[NIter-1] + DTVkl;
		TTVkl[NIter] > 0 ? TTVkl[NIter] : TTVkl[NIter] += 86400;
		if (DTVkl < 10000) //проверка!!
		{
			RezultOfCalculations.clear();
			this->ShootSplDownPoint(TTVkl[NIter],AimLandLongitude,RezultOfCalculations);
			TLPos[NIter] = RezultOfCalculations.back().lambda;
			NIter +=1;
		}
		else break;
	} 
	while (NIter < 10);	
	/*Date_Time RDT; //проверка
	RDT.SetFrmtTime(RezultOfCalculations.back().TVkl)*/
	delete[] TLPos; delete[] TTVkl;
	return RezultOfCalculations;
}

void Spacecraft::Aerodynamic(double max, double & s, double & s1, double & b){
	throw logic_error("not implement");
}

void Spacecraft::RemoveEngine(CppFlightDynamicOperations::Engine engine)
{
	if (_countOfEngines > 0)
	{
		vector<CppFlightDynamicOperations::Engine>::iterator it;
		it = find(Engines.begin(), Engines.end(), engine);
		if (it->Id_Engine >= 0)
		{
			Engines.erase(it);
			_countOfEngines -=1;
		}
		else
			throw invalid_argument("engine not found");
	}
	else
		throw invalid_argument("no engines at this spacecraft");
	
}

void Spacecraft::RemoveEngine(unsigned id)
{
	vector<CppFlightDynamicOperations::Engine>::iterator it;
	it = std::find_if(Engines.begin(), Engines.end(), [id](const CppFlightDynamicOperations::Engine &engine ) -> bool 
	{
		return engine.Id_Engine == id;
	});

	Engines.erase(it);
	_countOfEngines -= 1;
}

void Spacecraft::OrientFisk(double t)
{
	double V_A[6];         //массив(вектор) в инерциальной СК
	double P1[3][3];      //матрица перехода ОСК - ИСКТ
	double Prazv[3][3]; //матрица  программных разворотов
	double UgRazv[3] = {TurnX, TurnY, TurnZ};
	double gradV[3];
	if (ROr == 1 && RStab == 2 /*&& fabs(t - TVkl) <= 0.01*/) //первый вход для ориентации ИКВ и стабилизации ИСКТ
	{
		V_A[0] = x;
		V_A[1] = y;
		V_A[2] = z;
		V_A[3] = Vx;
		V_A[4] = Vy;
		V_A[5] = Vz;
		gradV[0] = _ax;
		gradV[1] = _ay;
		gradV[2] = _az;
		Mper(V_A,P1);
		CalcOptOrient(gradV, P1);
		//Razvor(UgRazv,Prazv);	
		//MultMatr(P1,Prazv, P2);
		return;
	}
/*
	if (ROr == 1 && RStab == 2 && fabs(t - TVkl) > 0.01)
	{
		V_A[0] = x;
		V_A[1] = y;
		V_A[2] = z;
		V_A[3] = Vx;
		V_A[4] = Vy;
		V_A[5] = Vz;
		Mper(V_A,P1);
		Razvor(UgRazv,Prazv);
		//CalcOptOrient(V_A, P2);
		return;
	}*/
		
	if (ROr == 1 && RStab == 1)
	{
		V_A[0] = x;
		V_A[1] = y;
		V_A[2] = z;
		V_A[3] = Vx;
		V_A[4] = Vy;
		V_A[5] = Vz;
		Mper(V_A,P1);
		Razvor(UgRazv,Prazv);	
		MultMatr(P1,Prazv, P2);
		//CalcOptOrient(V_A, P2);
		return;
	}
}

void Spacecraft::CalcOptOrient (double gradV[3], double OrientMatr[3][3])
{
	double lamdaOpt = 0.1;
	for (int i = 0; i < 3; i ++)
	{
		for (int j = i; j < 3; j++)
		{
			P2[i][j] = OrientMatr[i][j];
		}
	}
	P2[0][0] = P2[0][0] - lamdaOpt*gradV[0];
	P2[1][1] = P2[1][1] - lamdaOpt*gradV[1];
	P2[2][0] = P2[2][0] - lamdaOpt*gradV[2];
		
}