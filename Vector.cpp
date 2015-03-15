#include "Vector.h"
#include "OrbitElements.h"
//#include "Printer.h"
#include "Matrix.hpp"

using namespace YMatrix;

typedef Matrix<double,DenseMatrix,MathMatrix> MatrixM;
typedef ColVector<double,DenseMatrix,MathMatrix> CVectorM;

//статические члены
int Vector::ap;
int Vector::f;

// реализация конструкторов
Vector::Vector(unsigned long _vitok, double _t, double _x, double _y, double _z, double _vx, double _vy, double _vz, double _sb, Date_Time _dt,
			   double _fi, double _lambda, double _r, double _h, double _ro)
	: vitok(_vitok), t(_t), x(_x), y(_y), z(_z), Vx(_vx), Vy(_vy), Vz(_vz), Sb(_sb), DatTime(_dt), fi(_fi), lambda(_lambda), r(_r), h(_h), ro(_ro), calcVarAtm(true),
	  _ax(0), _ay(0), _az(0)
{}
Vector::~Vector()
{}

//реализация перегрузки операторов

Vector & Vector::operator = (const Vector & Vect)
{
    vitok = Vect.vitok; x = Vect.x; y = Vect.y; z = Vect.z; Vx = Vect.Vx; Vy = Vect.Vy; Vz = Vect.Vz; Sb = Vect.Sb; vitok = Vect.vitok;
	t = Vect.t; lambda = Vect.lambda; fi = Vect.fi; h = Vect.h; ro = Vect.ro; r = Vect.r; _ax = Vect._ax; _ay = Vect._ay; _az = Vect._az;
	DatTime = Vect.DatTime; calcVarAtm = Vect.calcVarAtm; _ax = Vect._ax; _ay = Vect._ay; _az = Vect._az;
	return *this;
}

Vector & Vector::operator = (double *X)
{
	x = X[1]; y = X[2]; z = X[3]; Vx = X[4]; Vy = X[5]; Vz = X[6]; t = X[0];
	return *this;
}

Vector operator + (const Vector & a, const Vector & b)
{
	Vector c;
	c.x      = a.x  + b.x;
	c.y      = a.y  + b.y;
	c.z      = a.z  + b.z;
	c.Vx     = a.Vx + b.Vx;
	c.Vy     = a.Vy + b.Vy;
	c.Vz     = a.Vz + b.Vz;
	c.Sb     = a.Sb;
	c.vitok  = a.vitok;
	c.t      = a.t;
	c.fi     = a.fi;
	c.lambda = a.lambda;
	c.h      = a.h;
	c.ro     = a.ro;
	c.r      = a.r;
	c.DatTime = a.DatTime;
	c.calcVarAtm = a.calcVarAtm;
	c._ax    = a._ax;
	c._ay    = a._ay;
	c._az    = a._az;
	return(c);
}

Vector operator - (const Vector & a, const Vector & b)
{
	Vector c;
	c.x  =     a.x  - b.x;
	c.y  =     a.y  - b.y;
	c.z  =     a.z  - b.z;
	c.Vx =     a.Vx - b.Vx;
	c.Vy =     a.Vy - b.Vy;
	c.Vz =     a.Vz - b.Vz;
	c.Sb =     a.Sb;
	c.vitok =  a.vitok;
	c.t =      a.t;
	c.fi =     a.fi;
	c.lambda = a.lambda;
	c.h =      a.h;
	c.r =      a.r;
	c.ro =     a.ro;
	c.DatTime = a.DatTime;
	return(c);
}

Vector operator * (double a , const Vector & b)
{
	Vector c;
	c.x       = a * b.x;
	c.y       = a * b.y;
	c.z       = a * b.z;
	c.Vx      = a * b.Vx;
	c.Vy      = a * b.Vy;
	c.Vz      = a * b.Vz;
	c.Sb      =     b.Sb;
	c.vitok   =     b.vitok;
	c.t       =     b.t;
	c.fi      =     b.fi;
	c.lambda  =     b.lambda;
	c.h       =     b.h;
	c.r       =     b.r;
	c.ro      =     b.ro;
	c.DatTime = b.DatTime;
	c._ax     = b._ax;
	c._ay     = b._ay;
	c._az     = b._az;
	return(c);
}


Vector operator * (const Vector & b, double a)
{
	Vector c;
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
	c.DatTime = b.DatTime;
	c._ax     = b._ax;
	c._ay     = b._ay;
	c._az     = b._az;
	return(c);
}

const double  Vector::operator [] (size_t i) const
{
	switch(i)
	{
		case 0: return this->t;			
		case 1: return this->x;
		case 2: return this->y;
		case 3: return this->z;
		case 4: return this->Vx;
		case 5: return this->Vy;
		case 6: return this->Vz;
		case 7: return this->Sb;
		case 8: return this->fi;
		case 9: return this->lambda;
		case 10: return this->h;
		case 11: return this->r;
		case 12: return this->ro;
		case 13: return this->vitok;
		default: return 0;
	}
}

void Vector::Print() const
{
	cout << setprecision (13) << setw(20)<< "   (" << vitok << ",  "<< h <<  ")    " << endl; 
}

void Vector::Perevod()
{
	const double ae = 6378.14; const double alfa = 1/298.256;
	const double r2 = x*x + y*y + z*z;
	r = sqrt(r2);
	double koef1 = PerRadGrad;
	lambda = atan2(y,x)*koef1;
	h = r - ae * (1.0 - alfa * z*z/r2);
	fi = atan2(z/sqrt(x*x + y*y), 1 - alfa*alfa)*koef1;
}


void Vector::PerevodGanimed()
{
	const double ag = 2631.2; const double alfa = 0;
	const double r2 = x*x + y*y + z*z;
	r = sqrt(r2);
	double koef1 = PerRadGrad;
	lambda = atan2(y,x)*koef1;
	h = r - ag * (1.0 - alfa * z*z/r2);
	fi = atan2(z/sqrt(x*x + y*y), 1 - alfa*alfa)*koef1;
}

void Vector::Altitude()
{
	const double ae = 6378.14; const double alfa = 1/298.256;
	r = sqrt(x*x + y*y + z*z);
	h = r - ae * (1.0 - alfa * ((z/r) * (z/r)));
}

void Vector::RKS(double dt, Atmosphera &Atm)
{
	Vector K1, K2, K3, K4;
	K1 = F(*this, t, Atm) * dt;
	K2 = F(*this + 0.5*K1, t + 0.5*dt, Atm) * dt;
	K3 = F(*this + 0.5*K2, t + 0.5*dt, Atm) * dt;
	K4 = F(*this +     K3, t + dt, Atm) * dt;
	*this = *this + 1.0/6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
	t += dt;
}

Vector  Vector::F (Vector &Vect, double dt, Atmosphera &Atm)
{
	Vector Rez;
	double Koef1, Koef2, Koef3;
	double V = sqrt(Vect.Vx * Vect.Vx + Vect.Vy * Vect.Vy + Vect.Vz * Vect.Vz);
	Koef1 = OMEGA_EARTH * OMEGA_EARTH;
	Koef2 = OMEGA_EARTH * 2.0;
	GetDensity(f,ap,Atm,Vect); //рассчет плотности атмосферы (потом надо убрать F and Ap)
	Koef3 = -Vect.Sb * Vect.ro * V * 1000;
	double *RezGpz = new double [3];
	gpz(Vect, RezGpz); // рассчет коэффициентов гравитационного поля Земли
	Rez.x = Vect.Vx;
	Rez.y = Vect.Vy;
	Rez.z = Vect.Vz;
	Rez.Vx = Koef3 * Vect.Vx + Koef2*Vect.Vy + Koef1*Vect.x + RezGpz[0];
	Rez.Vy = Koef3 * Vect.Vy - Koef2*Vect.Vx + Koef1*Vect.y + RezGpz[1];
	Rez.Vz = Koef3 * Vect.Vz + RezGpz[2];
	delete [] RezGpz;
	return Rez;
}

Vector & Vector::Prognoz(Date_Time &_datTime)
{
	int i = 0; unsigned j = 1;
	unsigned OrderMatr = 8;
	double **AM;
	AM = new double * [OrderMatr];
	for (unsigned k = 0; k < OrderMatr; ++k)
	{
		AM[k] = new double [OrderMatr];
	}
	AM[0][0] = t;
	AM[0][1] = x;
	AM[0][2] = y;
	AM[0][3] = z;
	AM[0][4] = Vx;
	AM[0][5] = Vy;
	AM[0][6] = Vz;
	AM[0][7] = 0;
	Perevod();
	bool Half_Vitok = false;
	bool Izm_Vitka = false;
	bool Progn = true; // переменная для инициализации вызова прогноза в Interpolmatrix
	Vector VektEkv;
	double dt = 20;
	Atmosphera Atm;
	Atm.NMounth = DatTime.Get_Mounth() - 1;
	Atm.StarTime(DatTime);
	Atm.Sun(t);
	unsigned checkVitok = vitok;
	while (_datTime > DatTime || !(i > 4) )
	{
		RKS(dt, Atm);
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, *this, j, 0, & Half_Vitok, &Izm_Vitka, Progn);
		DatTime.SetFrmtTime(t);
		if (t - COUNT_SEC_OF_DAY > 0)
		{	
			DatTime.Calendar();
			Atm.StarTime(DatTime);
			t -= COUNT_SEC_OF_DAY;
			DatTime.SetFrmtTime(t);
			if (i > 1000 && checkVitok == vitok)
				throw(-1);
		}
		Izm_Vitka = false;
		if (Half_Vitok)
			Atm.Sun(t);
		Half_Vitok = false;
	}
	double argInterpol = _datTime.GetSecOfDay();
	argInterpol == 0 ? argInterpol += COUNT_SEC_OF_DAY : argInterpol;
	if ((argInterpol - 8*t) > COUNT_SEC_OF_DAY/2 )
		this->DatTime.Day -= 1;
	Ilagr(8, AM, 0, argInterpol, *this, 7);
	this->DatTime.SetFrmtTime(t);
	if (t >= COUNT_SEC_OF_DAY)
		t -= COUNT_SEC_OF_DAY;
	for (unsigned i = 0; i < OrderMatr; ++i) // удаление матрицы
	{
		delete [] AM [i];
	}
	delete [] AM;
	Perevod();
	return *this;
}

vector<Vector> Vector::Prognoz(unsigned _Vitok)
{

	/*==================================== research part of using gpz 36*36 check F method!!!=======================================
	*/
	//=================gpz3636
	/*CppService::Printer::ReadFileName = "..\\dpz4-36.txt";
	CppService::Printer::ReadGpzFromFile2(C0,C,D);
	for (int i = 0; i < 75; i++)
		SQ[i] = sqrt(i + 1);*/
	/*CppService::Printer::OutFileName = "..\\altitude.txt";
	CppService::Printer::Precision = 10;
	double massh[1000];
	double masst[1000];*/

	//==============================================


	int i = 0; unsigned j = 1;
	unsigned OrderMatr = 8;
	double **AM;
	AM = new double * [OrderMatr];
	for (unsigned k = 0; k < OrderMatr; ++k)
	{
		AM[k] = new double [OrderMatr];
	}
	AM[0][0] = t;
	AM[0][1] = x;
	AM[0][2] = y;
	AM[0][3] = z;
	AM[0][4] = Vx;
	AM[0][5] = Vy;
	AM[0][6] = Vz;
	AM[0][7] = 0;
	Perevod();
	bool Half_Vitok = false;
	bool Izm_Vitka = false;
	bool Progn = true; // переменная для инициализации вызова прогноза в Interpolmatrix
	Vector VektEkv;
	double dt = 20;
	Atmosphera Atm;
	Atm.NMounth = DatTime.Get_Mounth() - 1;
	Atm.StarTime(DatTime);
	Atm.Sun(t);
	DatTime.SetFrmtTime(t);
	vector<Vector> rezOfProgn;
	unsigned checkVitok = vitok;
	//======================research gpz ================================
	/*massh[0] = h;
	masst[0] = t;*/
	//==============================
	while (_Vitok != vitok)
	{		
		RKS(dt, Atm);
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, *this, j, 0, & Half_Vitok, &Izm_Vitka, Progn);
		if (Izm_Vitka)
		{
			VektEkv = *this;
			Ilagr(8, AM, 3, 0, VektEkv, 7);
			if (VektEkv.t >= COUNT_SEC_OF_DAY)
				VektEkv.t -= COUNT_SEC_OF_DAY;
			DatTime.SetFrmtTime(VektEkv.t);
			PerevodV(VektEkv);
			VektEkv.DatTime = DatTime;
			rezOfProgn.push_back(VektEkv);
			DatTime.Print();
			VektEkv.Print();
		}
		if (t - COUNT_SEC_OF_DAY > 0)
		{	
			DatTime.Calendar();
			Atm.StarTime(DatTime);
			t -= COUNT_SEC_OF_DAY;
			if (i > 1000 && checkVitok == vitok)
				throw(-1);
		}
		Izm_Vitka = false;
		if (Half_Vitok)
			Atm.Sun(t);
		Perevod();

		//======================research gpz ================================
		/*massh[i] = h;
		masst[i] = t;*/
		//==========================================

		Half_Vitok = false;
	}
	*this = VektEkv; 
	for (unsigned i = 0; i < OrderMatr; ++i) // удаление матрицы
	{
		delete [] AM [i];
	}
	delete [] AM;

	/* research part of using gpz 36*36 check F method!!!
	*/
	//CppService::Printer::PrintDoubleToFile(masst, massh, i);

	return rezOfProgn;
}

Vector & Vector::Prognoz()
{
	int i = 0; unsigned j = 1;
	unsigned OrderMatr = 8;
	double **AM;
	AM = new double * [OrderMatr];
	for (unsigned k = 0; k < OrderMatr; ++k)
	{
		AM[k] = new double [OrderMatr];
	}
	this->z > 0 ? this->z : this->z = fabs(this->z); 
	AM[0][0] = t;
	AM[0][1] = x;
	AM[0][2] = y;
	AM[0][3] = z;
	AM[0][4] = Vx;
	AM[0][5] = Vy;
	AM[0][6] = Vz;
	AM[0][7] = 0;
	Perevod();
	bool Half_Vitok = false;
	bool Izm_Vitka = false;
	bool Progn = true; // переменная для инициализации вызова прогноза в Interpolmatrix
	Vector VektEkv;
	double dt = 20;
	Atmosphera Atm;
	Atm.NMounth = DatTime.Get_Mounth() - 1;
	Atm.StarTime(DatTime);
	Atm.Sun(t);
	while (!Izm_Vitka)
	{		
		RKS(dt, Atm);
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, *this, j, 0, & Half_Vitok, &Izm_Vitka, Progn); 
		if (Izm_Vitka)
		{			
			VektEkv = *this;
			Ilagr(8, AM, 3, 0, VektEkv, 7);
			if (VektEkv.t >= COUNT_SEC_OF_DAY)
				VektEkv.t -= COUNT_SEC_OF_DAY;
			DatTime.SetFrmtTime(VektEkv.t);
			PerevodV(VektEkv);
			VektEkv.DatTime = DatTime;
			DatTime.Print();
			VektEkv.Print();									
		}
		if (abs(t- COUNT_SEC_OF_DAY) <= dt/2)
		{	
			DatTime.Calendar();
			Atm.StarTime(DatTime);
		}
		if (t >= 88000)
			t -= COUNT_SEC_OF_DAY;
		if (Half_Vitok)
			Atm.Sun(t);
		Perevod();
		if (i > 1000)
			throw(-1);
	}
	*this = VektEkv; 
	for (unsigned i = 0; i < OrderMatr; ++i) // удаление матрицы
	{
		delete [] AM [i];
	}
	delete [] AM;
	return *this;
}

Vector Vector::Prognoz(double u)
{
	u = u*PI/180;
	if (u > 0 && u <= 2*PI)
	{
		Vector vect(*this);
		OrbitElements orbElements(vect);
		int i = 0; unsigned j = 1;
		unsigned OrderMatr = 8;
		double **AM;
		AM = new double * [OrderMatr];
		for (unsigned k = 0; k < OrderMatr; ++k)
		{
			AM[k] = new double [OrderMatr];
		}
		vect.z > 0 ? vect.z : vect.z = fabs(vect.z); 
		AM[0][0] = t;
		AM[0][1] = x;
		AM[0][2] = y;
		AM[0][3] = z;
		AM[0][4] = Vx;
		AM[0][5] = Vy;
		AM[0][6] = Vz;
		AM[0][7] = 0;
		vect.Perevod();
		bool Half_Vitok = false;
		bool Izm_Vitka = false;
		bool Progn = true; // переменная для инициализации вызова прогноза в Interpolmatrix
		double dt = 20;
		Atmosphera Atm;
		Atm.NMounth = DatTime.Get_Mounth() - 1;
		Atm.StarTime(DatTime);
		Atm.Sun(vect.t);
		bool hasU = false;
		while (!(hasU && (unsigned)i > OrderMatr))
		{
			vect.RKS(dt, Atm);
			i += 1;
			orbElements.SetVector(vect);
			orbElements.GetElements();
			InterpolMatrix(AM, OrderMatr, i, vect, j, orbElements.GetU(), & Half_Vitok, &Izm_Vitka, Progn); 
			if (abs(vect.t- COUNT_SEC_OF_DAY) <= dt/2)
			{	
				DatTime.Calendar();
				Atm.StarTime(DatTime);
				if (i > 100000)
					throw(-1);
			}
			if (vect.t >= 88000)
				vect.t -= COUNT_SEC_OF_DAY;
			if (Half_Vitok)
				Atm.Sun(vect.t);
			vect.Perevod();
			orbElements.SetVector(vect);
			orbElements.GetElements();
			fabs(u - orbElements.GetU()) > 0.01 ? hasU = false: hasU = true;
		}
		Ilagr(8, AM,7,orbElements.GetU(),vect,7);
		for (unsigned i = 0; i < OrderMatr; ++i) // удаление матрицы
		{
			delete [] AM [i];
		}
		delete [] AM;
		return vect;
	}
	throw(-2);//throw("Invalid argument exception(u)");
}

double Vector::CalcCornOfSun (const Atmosphera & _atm)
{
	double rm = 1/sqrt(x*x + y*y + z*z);
	Atmosphera atm = _atm;
	atm.StarTime(DatTime);
	atm.Sun(t);
	double bet = atm.ASol - atm.S0 - OMEGA_EARTH*(t - 10800.);
	double c1  = rm*(z*atm.SDSol + atm.CDSol*(x*cos(bet) + y*sin(bet)));
	if (fabs(c1) > 1)
		c1 = Sign(c1);
	double s1  = sqrt(1.0 - c1*c1);
	double soz = PerRadGrad*asin(s1);
	if (c1 > 0)
		soz = 180.0 - soz;
	double x1 = x + 3.*Vx;
	double y1 = y + 3.*Vy;
	double z1 = z + 3.*Vz;
	double c2 = rm*(z1*atm.SDSol + atm.CDSol*(x1*cos(bet) + y1*sin(bet)));
	if (fabs(c2) > 1)
		c2 = Sign(c2);
	double s2 = sqrt(1.0 - c2*c2);
	double soz1 = PerRadGrad*asin(s2);
	if (c2 > 0)
		soz1 = 180. - soz1;
	double d = soz1 - soz;
	double dsoz = Sign(d); // 
	return soz;
}
void PerevodV(Vector & Vect)
{
	const double ae = 6378.14; const double alfa = 1/298.256;
	Vect.r = sqrt(Vect.x*Vect.x + Vect.y*Vect.y + Vect.z*Vect.z);
	double koef1 = 180/PI;
	Vect.lambda = atan2(Vect.y,Vect.x)*koef1;
	//Vect.h = Vect.r - ae * (1 - alfa * ((Vect.z/Vect.r) * (Vect.z/Vect.r)));
	Vect.fi = atan2(Vect.z/sqrt(Vect.x*Vect.x + Vect.y*Vect.y), 1 - alfa*alfa)*koef1;
}

double Sinusl (const Vector &Vect,int n)
{
	if (n == 0)
		return 0;
	if (n == 1)
	{
		double r1 = sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y));
		return Vect.y/r1;
	}
	if (n == 2)
	{
		double r1=sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y));
		return 2 * Vect.y/r1 * Vect.x/r1;
	}
	if (n > 2)
	{
		double *s = new double [n+1];
		double r1 = sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y));
		double clamda = Vect.x/r1; 
		s[0] = 0;
		s[1] = Vect.y/r1;
		s[2] = 2*s[1]*clamda;
		for (int i = 3; i <= n; ++i)
		{
			s[i] = 2*clamda*s[i-1] - s[i-2];
		}
		double res = s[n];
		delete [] s;
		return res;
	} 
	return 0;
}
double Cosinusl (const Vector &Vect,int n)
{
	if (n == 0)
		return 1;
	if (n == 1)
	{
		double r1 = sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y));
		return Vect.x/r1;
	}
	if (n == 2)
	{
		double r1 = sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y));
		double clambda = Vect.x/r1;
		double slambda = Vect.y/r1;
		return clambda*clambda - slambda*slambda ;
	}
	if (n > 2)
	{
		double *s = new double [n+1];
		double r1 = sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y));
		double clambda = Vect.x/r1;
		double slambda = Vect.y/r1;
		s[0] = 1;
		s[1] = clambda;
		s[2] =  clambda*clambda - slambda*slambda;
		for (int i = 3; i <= n; ++i)
		{
			s[i] = 2*clambda*s[i-1] - s[i-2];
		}
		double res = s[n];
		delete [] s;
		return res;
	} 
	return 0;
}
void gpz (const Vector &Vect, double *Rez)
{
	const size_t grm = 9;
	const double Re = 6378.14;
	double c[9][9] = {
	{0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0},
	{-1.0826101104e-3, 0.0000000000000, 1.5672672608e-6, 0.0000000000,    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
	{ 2.5399212586e-6, 2.1667276402e-6, 3.1765547374e-7, 1.0151474989e-7, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
	{ 1.6080000000e-6,-5.1228898095e-7, 7.8262379212e-8, 5.8506440427e-8,-3.5707767262e-9, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
	{ 2.2553048574e-7,-4.9668232638e-8, 1.0567807134e-7,-1.4535124398e-8,-2.3358829608e-9, 3.3732644632e-10, 0.0000000000, 0.0000000000, 0.0000000000},
	{-5.3362158877e-7,-7.2385212907e-8, 7.3397969350e-9, 4.1467779294e-11,-1.77917232e-10,-2.4534793021e-10, 1.8638380731e-12, 0.0000000000, 0.0000000000},
	{ 3.6018745120e-7, 1.7712586324e-7, 3.3167593909e-8, 3.7327646251e-9,-6.1794721961e-10, 7.0784332144e-13,-2.2696987675e-11,-1.4840428412e-13, 0.0000000000},
	{ 1.9790907003e-7, 1.9241159586e-8, 2.3818926715e-9,-1.2132029804e-10,-3.2107885416e-10,-3.6199695814e-12,-1.9270783444e-12, 1.8356596162e-13,-1.0453061703e-13}};
	double d[9][9] = {
	{0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0},
	{ 0.0000000000, 0.000000000, -8.9078616963e-7, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
	{ 0.0000000000, 2.6463024518e-7,-2.0937936065e-7, 1.9591788955e-7, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
	{ 0.0000000000,-4.2880485072e-7, 1.5026376809e-7,-1.3147514703e-8, 6.4020434796e-9, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
	{ 0.0000000000,-7.6215046633e-8,-5.0330597528e-8,-8.0603871663e-9, 2.0244318993e-10,-1.6201518371e-9, 0.0000000000, 0.0000000000, 0.0000000000},
	{ 0.0000000000, 2.2817077982e-8,-4.7273268395e-8, 1.2647672685e-9,-1.6693936025e-9,-4.2855181232e-10,-5.3818324361e-11, 0.0000000000, 0.0000000000},
	{ 0.0000000000, 7.4656355581e-8, 1.1055864636e-8,-2.7608372321e-9,-1.8050004696e-10,4.9549032500e-12, 1.1799657201e-11, 8.9042570477e-13, 0.0000000000},
	{ 0.0000000000, 4.3979793340e-8, 6.0779330237e-9,-5.7627141573e-10,6.1344333925e-11,9.9549163490e-12, 8.5182448560e-12, 3.2124043284e-13, 8.0310108211e-14}};
	double r1 = sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y));
	double r  = sqrt((Vect.x*Vect.x) + (Vect.y*Vect.y) + (Vect.z*Vect.z));
	double s2fi,sfi,cfi,c2fi,slam,clam;
	sfi=Vect.z/r;
	s2fi=sfi*sfi;
	cfi=r1/r;
	c2fi=cfi*cfi;
	slam=Vect.y/r1;
	clam=Vect.x/r1;
	double p[grm][grm];
	double p1[grm][grm]; //p`
	p [2][0] = (0.5*(3*s2fi-1)); 
	p [2][1] = 3*sfi*cfi;
	p [2][2] = 3*(1-s2fi);
	p [3][0] = 0.5*(5*s2fi*sfi-3*sfi); p[3][1]=1.5*(5*s2fi-1)*cfi; p[3][2]=15*sfi*c2fi; p[3][3]=15.0*c2fi*cfi;
	p1[2][0] = 3*sfi*cfi;
	p1[2][1] = 3*(1-2*s2fi);
	p1[2][2] = -6*sfi*cfi;
	p1[3][0] = 1.5*(5*s2fi*cfi-cfi);
	p1[3][1] = 1.5*sfi*(11-15*s2fi);
	p1[3][2] = 15*cfi*(1-3*s2fi);
	p1[3][3] = -45*c2fi*sfi;
	for (int n = 4; n < grm; n++) 
	{
		for (int m=0; m <= n; m++) 
		{
			if (m <= n-2) 
			{
				p [n][m] = (((2*n - 1)*sfi*p[n - 1][m] - (n-1 + m)*p[n-2][m])/(n-m));
				p1[n][m] = (((2*n-1)*(cfi*p[n-1][m] + sfi*p1[n-1][m]) - (n-1+m)*p1[n-2][m])/(n-m));
			}
			else
			{
				p [n][m] = (2*(m-1)*(sfi/cfi)*p[n][m-1] - (n+m-1)*(n-m+2)*p[n][m-2]);
				p1[n][m] = ((2*(m-1)/c2fi)*(sfi*cfi*p1[n][m-1] + p[n][m-1]) - (n+m-1)*(n-m+2)*p1[n][m-2]);
			}
		}
	}
	double a_r,a_lamda,a_fi,c1,c11,b1,b2,b3,c10;
	double koef2,koef1,koef3,koef4;
	koef2=koef3=koef4=0;
	for (int n = 2; n < grm; n++)                
	{
		b1 = b2 = b3 = 0;
		koef1 = Re/r;
		c10 = koef1;
		koef1 = pow((Re/r),n-1);
		for (int m=0; m <= n; m++) 
		{
			double csl, snl;
			csl = Cosinusl(Vect,m);
			snl = Sinusl(Vect,m);
			c1  = c[n][m]*csl+d[n][m]*snl;
			c11 = c[n][m]*snl-d[n][m]*csl;
			b1 += c1*p[n][m];
			b2 += c1*p1[n][m];
			b3 += c11*p[n][m]*m;
		}
		koef1  = c10*koef1;
		koef2 += koef1*(n+1.0)*b1;
		koef3 += koef1*b2;
		koef4 += koef1*b3;
	}
	koef1   =  MU/(r*r);
	a_r     = -koef1*(1.0 + koef2);
	a_fi    =  koef1*koef3;
	a_lamda = -(koef1/cfi)*koef4;
	Rez[0]  = a_r*cfi*clam - a_fi*sfi*clam - a_lamda*slam;
	Rez[1]  = a_r*cfi*slam - a_fi*sfi*slam + a_lamda*clam;
	Rez[2]  = a_r*sfi + a_fi*cfi;
}

void InterpolMatrix (double **AM, unsigned n, int i,  Vector & Vect, unsigned &j, double arg, bool *Half_Vitok, bool *Izm_vitka, bool Progn)
{
	if (fmod(float(i),float(n)) != 0)
	{
		AM [j][0] = Vect.t;
		AM [j][1] = Vect.x;
		AM [j][2] = Vect.y;
		AM [j][3] = Vect.z;
		AM [j][4] = Vect.Vx;
		AM [j][5] = Vect.Vy;
		AM [j][6] = Vect.Vz;
		AM [j][7] = arg;
		if (Progn && Sign(AM[j][3]) == 1 && Sign(AM[j-1][3]) == -1)
		{
			Vect.vitok +=1;
			*Half_Vitok = true;
			*Izm_vitka  = true;
		}
		if (Progn && Sign(AM[j][3]) == -1 && Sign(AM[j-1][3]) == 1)
			*Half_Vitok = true;
		j += 1;
	}
	else
	{
		j = 0;
		AM [j][0] = Vect.t;
		AM [j][1] = Vect.x;
		AM [j][2] = Vect.y;
		AM [j][3] = Vect.z;
		AM [j][4] = Vect.Vx;
		AM [j][5] = Vect.Vy;
		AM [j][6] = Vect.Vz;
		AM [j][7] = arg;
		if (Progn && Sign(AM[j][3]) == 1 && Sign(AM[n-1][3]) == -1)
		{
			Vect.vitok +=1;
			*Half_Vitok = true;
			*Izm_vitka = true;
		}
		if (Progn && Sign(AM[j][3]) == -1 && Sign(AM[n-1][3]) == 1)
			*Half_Vitok = true;
		j += 1;
	}
}

void Ilagr (unsigned n, double **AM, unsigned NA, double arg, Vector & Vect, unsigned NumOfCols)
{
	double C1;
	double *R = new double [10]; // максимальное число столбцов
	double *X = new double [NumOfCols];
	for (unsigned i = 0; i < n; ++ i)
	{
		C1 = 1;
		for (unsigned j = 0; j < n; ++ j)
		{
			if (i == j) continue;
			C1 *= (arg-AM[j][NA])/(AM[i][NA] - AM[j][NA]);
			R[i] = C1;
		}
	}
	for (unsigned i = 0; i < NumOfCols; ++ i)
	{
		X[i] = 0;
		for (unsigned j = 0; j < n; ++ j)
		{
			X[i] += R[j]*AM[j][i];
		}
	}
	Vect = X;
	delete [] R;
	delete [] X;
}

void Ilagr (unsigned _orderMatrEnt, double **_matrEnter, unsigned NA, double arg ,double *_rezult, unsigned _numOfColsRez)
{
	double C1;
	double *R = new double [10]; // максимальное число столбцов
	for (unsigned i = 0; i < _orderMatrEnt; ++ i)
	{
		C1 = 1;
		for (unsigned j = 0; j < _orderMatrEnt; ++ j)
		{
			if (i == j) continue;
			C1 *= (arg-_matrEnter[j][NA])/(_matrEnter[i][NA] - _matrEnter[j][NA]);
			R[i] = C1;
		}
	}
	for (unsigned i = 0; i < _numOfColsRez; ++ i)
	{
		_rezult[i] = 0;
		for (unsigned j = 0; j < _orderMatrEnt; ++ j)
		{
			_rezult[i] += R[j]*_matrEnter[j][i];
		}
	}
	delete [] R;
}

void GetDensity (int F, double Ap,  Atmosphera & Atm, Vector & Vect)
{
	if (Vect.h >= 120)
	{
		double ron, K0, K1, K2, K3, K4;
		int F0 = Atm.Selection_F0(F);
		double h2, h3, CFi, Cfi, Beta, Ad, Kp, radius;
		radius = 1.0 / Vect.r;
		Atm.SelectKoef(F0,Vect.h);
		h2   =  Vect.h * Vect.h;
		h3   =  h2 * Vect.h;
		ron  = expl(Atm.A[0] - Atm.A[1] * sqrt(Vect.h - Atm.A[2]));
		K0   = 1 + (Atm.L[0] + Atm.L[1]*Vect.h + Atm.L[2]* h2)*(F - F0);
		Beta = Atm.ASol - Atm.S0 - OMEGA_EARTH*(Vect.t - 10800.0) + Atm.Fi;
		CFi  = (Vect.z * Atm.SDSol + Atm.CDSol * (Vect.x * cos(Beta) + Vect.y * sin(Beta)))*radius;
		Cfi  = 0.5 * (1.0 +CFi);
		K1   = 1 + (Atm.C[0] + Atm.C[1]*Vect.h + Atm.C[2]*h2 + Atm.C[3]*h3)*pow(Cfi,(Atm.N[0] + Atm.N[1]*Vect.h)*0.5);
		Ad   = Atm.Ad_D(Atm.d);
		K2   = 1 + (Atm.D[0] + Atm.D[1]*Vect.h + Atm.D[2]*h2)*Ad;
		Kp   = Atm.Kp_Ap(Ap);
		K3   = 1;
		double Kp2 = Kp*Kp;
		K4 = 1 + (Atm.E[0] + Atm.E[1]*Vect.h + Atm.E[2]*h2 + Atm.E[3]*h3)*(Atm.E[4] + Atm.E[5]*Kp + Atm.E[6]*Kp2);
		Vect.ro = ron*K0*K1*K2*K3*K4;
	}
	else 
	{
		Vect.ro     = Atm.LowAtm(Vect.h);
		double R12  = Vect.x*Vect.x + Vect.y*Vect.y;
		double R121 = sqrt(R12);
		Vect.fi     = asin(Vect.z/sqrt((0.9933055926*R121)*(0.9933055926*R121) + Vect.z*Vect.z));
		double CC1  = 1/R121;
		double Sla  = Vect.y*CC1;
		double Cla  = Vect.x*CC1;
		Vect.lambda = asin(Sla);
		if (Cla > 0 && Vect.calcVarAtm)
			Atm.VarAtm(Vect.h, Vect.fi, Vect.lambda);
		else if (Cla <= 0 && Vect.calcVarAtm)
		{
			Vect.lambda < 0 ? Vect.lambda = -PI - Vect.lambda: Vect.lambda = PI - Vect.lambda;
			Atm.VarAtm(Vect.h, Vect.fi, Vect.lambda);
		}
		Vect.ro = Vect.ro*(1.0 + Atm.DRo); 
	}
}

double Vector::DLEkv(double _l1, double _l2)
{
	double dl = _l2 - _l1;
	if (_l2 > 0 && _l1 < 0)
		dl = -360 + dl;
	return dl;
}

unsigned Vector::nSutCirc()
{
	if (fabs(z) > 5.0 && fabs(fi) > 360.0)
		return 0;
	double v2 = (Vx - OMEGA_EARTH*y)*(Vx - OMEGA_EARTH*y) + (Vy + OMEGA_EARTH*x)*(Vy + OMEGA_EARTH*x) + Vz*Vz;
	double ak = r*v2/MU;
	double a  = r/(2.0 - ak);
	double T = 2.0*PI*a*sqrt(a/MU);
	double dl = -OMEGA_EARTH*T*PerRadGrad; //межвитковое расстояние
	double l1 = lambda;
	if (lambda > 20.0)
		l1 = lambda - 360;
	double c = (l1 - 20.0)/dl;
	unsigned nSutCircle = int(c) + 1;
	return nSutCircle;
}

unsigned Vector::nSutCirc(double _dLekv)
{
	if (fabs(_dLekv) < 0.0001)
	{
		if (abs(z) > 5 && abs(fi) > 360)
			return 0;
		double v2 = (Vx - OMEGA_EARTH*y)*(Vx - OMEGA_EARTH*y) + (Vy + OMEGA_EARTH*x)*(Vy + OMEGA_EARTH*x) + Vz*Vz;
		double ak = r*v2/MU;
		double a  = r/(2.0 - ak);
		double T = 2.0*PI*a*sqrt(a/MU);
		_dLekv = -OMEGA_EARTH*T*PerRadGrad; //межвитковое расстояние
	}
	double l1 = lambda;
	if (lambda > 20.0)
		l1 = lambda - 360;
	double c = (l1 - 20.0)/_dLekv;
	unsigned nSutCircle = int(c) + 1;
	return nSutCircle;
}

double Vector::Tet()
{

	double v = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
	double r = sqrt(x*x + y*y + z*z);
	return asin((x*Vx + y*Vy + z*Vz)/(r*v));
}

double** Vector::CreateSquareMatrix(unsigned size)
{
	double ** b;
	b = new double * [size];
	for (unsigned k = 0; k < size; ++k)
	{
		b[k] = new double [size];
	}
	return b;
}

void Vector::DeleteSquareMatrix(unsigned size, double** matrix)
{
	for (unsigned i = 0; i < size; ++i) 
	{
		delete [] matrix [i];
	}
	delete [] matrix;
}

void Vector::Polinom(bool & kk, int n, int m, double* ti, double* xi, double arg, double* a, double& xa)
{
	if (!kk)
	{
		double * b = new double[m];
		double * bt = new double[m];
		for (int i = 0; i < m; ++i)
			b[i] = ti[i];
		int n2 = 2 * n - 1;
		double * c = new double[n2];
		for (int i = 0; i < m; ++i)
			bt[i] = 1.0;
		int l;
		for (int i = 0; i < n2; i++)
		{
			l = n2 - (i + 1);
			c[l] = 0;
			for (int j = 0; j < m; j++)
			{
				c[l]  = c[l]  + bt[j];
				bt[j] = bt[j] * b[j];
			}
		}
		int nc = 1;
		double ** ak = CreateSquareMatrix(n);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				ak[i][j] = c[j + nc - 1];
			}
			nc += 1;
		}		
		MatrixM ak1(n,n);
		ak1 = ak;
		double det = Det<MatrixM>()(ak1);
		if (fabs(det) < 1e-20)
			throw logic_error("det equal zero");
		MatrixM ainv;
		Inverse<MatrixM>()(ak1, ainv);
		double * ct = new double[n];
		for (int i = 0; i < m; ++i)
			bt[i] = 1.0;
		for (int i = 0; i < n; i++)
		{
			l = n - (i + 1);
			ct[l] = 0.0;
			for (int j = 0; j < m; j++)
			{
				ct[l] = ct[l] + xi[j]*bt[j];
				bt[j] = bt[j] * b[j];
			}
		}
		CVectorM cv(n);
		cv = ct;
		CVectorM aRes(n);
		aRes = ainv*cv;
		for (int i = 0; i < n; i++)
			a[i] = aRes[i];
		delete [] b;
		delete [] bt;
		delete [] c;
		delete [] ct;
		DeleteSquareMatrix(n, ak);
		kk = true;
		double r1 = 1.0;
		xa = 0.0;
		for (int i = 0; i < n; i++)
		{
			xa = xa + r1 * a[n - (i + 1)];
			r1 = r1 * arg;
		}
	}
	else
	{
		double r1 = 1.0; 
		xa = 0.0;
		for (int i = 0; i < n; i++)
		{
			xa = xa + r1 * a[n - (i + 1)];
			r1 = r1 * arg;
		}
	}
}

/*
void Vector::gpz36(const Vector& Vect, double* rez)
{

	double c20, c30, C40, C50, C60;                                     //
	double C21, c22, C31, C32, C41, C42, C51, C52, C61, C62;            //гармонические коэффициенты
	double S21, S22, S31, S32, S41, S42, S51, S52, S61, S62; 

	c20 = C0[1]; c30 = C0[2]; C40 = C0[3];
	C50 = C0[4]; C60 = C0[5];

	C21 = C[0]; c22 = C[1];
	C31 = C[2]; C32 = C[3];
	C41 = C[4]; C42 = C[5];
	C51 = C[6]; C52 = C[7];
	C61 = C[8]; C62 = C[9];

	S21 = D[0]; S22 = D[1];
	S31 = D[2]; S32 = D[3];
	S41 = D[4]; S42 = D[5];
	S51 = D[6]; S52 = D[7];
	S61 = D[8]; S62 = D[9];

	int L, LL, I, J, K1, K2, K3;
	double CK2 = 0.0, U, SS1 = 0.0, SS2 = 0.0, SS3 = 0.0, CK;
	double S1, S2, S3 = 0, MU2;
	double RO, CF, ZP, V2, U2;
	double VD = 0.0, UD = 0.0, V1, V, U1;
	double R2, R, SF, CX, CY, AM, RZR;

	R2 = Vect.r * Vect.r;
	MU2 = -MU / R2;
	R2 = sqrt(R2);
	R = 1.0 / R2;
	SF = Vect.z * R;
	CX = Vect.x * R;
	CY = Vect.y * R;

	rez[0] = MU2 * CX;
	rez[1] = MU2 * CY;
	rez[2] = MU2 * SF;

	AM = MU * R;
	RZR = a_e * R; // !!!!!check
	RO = sqrt(Vect.x*Vect.x + Vect.y*Vect.y);
	RO = 1.0 / RO;
	CF = R2 * RO;
	ZP = Vect.z * RO;
	U2 = R;
	V2 = 0.0;
	U1 = SQ[2] * R * SF * RZR;
	V1 = 0.0;
	S1 = 0.0;
	S2 = 0.0;
	LL = -1;
	int MM = 36;
	int NM = MM;
	for (J = 0; J <= MM; J++)
	{
		LL = LL + J;
		L = LL;

		NM = MM;
		if (J == 0)
		{
			CK2 = 1.0 / SQ[2];
			for (I = 1; I <= NM; I++)
			{
				K3 = I + I + 1;
				U2 = SF * U1 - RZR * CK2 * U2;
				CK2 = (I + 1) / (SQ[K3 - 1] * SQ[K3 + 1]);
				U = RZR / CK2 * U2;
				if (I >= 2)
				{
					S1 = S1 + K3 * C0[I - 1] * U2;
					S2 = S2 + (I + 1) * C0[I - 1] * U1;
				}
				U2 = U1;
				U1 = U;
			}	
			SS1 = S1;
			SS2 = S2;
			S1 = 0.0;
			S2 = 0.0;
			S3 = 0.0;
			SS3 = 0.0;
		} 	
		else
		{
			L = LL;
			if (J > 1)
			{
				K1 = J + J;
				CK = RZR * SQ[K1] / SQ[K1 - 1];
				U1 = CK * (CX * UD - CY * VD);
				V1 = CK * (CY * UD + CX * VD);
			}
			else
			{
				CK = RZR * SQ[2] * R;
				U1 = CK * CX;
				V1 = CK * CY;
			}
			UD = U1;
			VD = V1;
			for (I = J; I <= NM; I++)
			{
				K1 = I + J;
				K2 = I - J;
				K3 = I + I + 1;
				if (I > J)
				{
					CK2 = RZR * CK2;
					U2 = SF * U1 - CK2 * U2;
					V2 = SF * V1 - CK2 * V2;
				}
				else
				{
					U2 = SF * U1;
					V2 = SF * V1;
				}
				CK2 = SQ[K1] * SQ[K2] / (SQ[K3 - 1] * SQ[K3 + 1]);
				CK = RZR / CK2;
				U = CK * U2;
				V = CK * V2;

				if (I >= 2)
				{
					S1 = S1 + K3 * (C[L - 1] * U2 + D[L - 1] * V2);
					S2 = S2 + (I + 1) * (C[L - 1] * U1 + D[L - 1] * V1);
					SS3 = SS3 + C[L - 1] * V1 - D[L - 1] * U1;
				}
				U2 = U1;
				U1 = U;
				V2 = V1;
				V1 = V;
				L = L + I;
			}	
			S3 = S3 + SS3 * J;
			SS3 = 0.0;
		}	
	}	

	S1 = SS1 + S1;
	S2 = SS2 + S2;

	ZP = ZP * S1 - CF * S2;
	rez[0] = rez[0] + AM * CF * (CX * ZP + CF * CY * S3);
	rez[1] = rez[1] + AM * CF * (CY * ZP - CF * CX * S3);
	rez[2] = rez[2] - AM * S1;

}*/
Vector Vector::F20(Vector & SC, double dt, Atmosphera &Atm)
{
	const double Rzsr  = 6371.0;
	const double TAI20 = -0.10183014444;
	const double omz2  = 53.174954e-10; 
	const double A100  = 62.5648249885;
	const double domz  = 14.58423171e-5;
	const double Rekv  = 6378.14;
	const double szat  = 0.003352824419;
	Vector Rez;
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
	Rez.Vx =  domz*SC.Vy + SC.x*C1; //*1.5??
	Rez.Vy = -domz*SC.Vx + SC.y*C1;
	Rez.Vz = -SC.z*C7*(A100 + C6*(C8 - 3 *R2));
	return Rez;
}

void Vector::RKS20(double dt, Atmosphera &Atm)
{
	Vector K1, K2, K3, K4;
	K1 = F20(*this, t, Atm) * dt;
	K2 = F20(*this + 0.5*K1, t + 0.5*dt, Atm) * dt;
	K3 = F20(*this + 0.5*K2, t + 0.5*dt, Atm) * dt;
	K4 = F20(*this +     K3, t + dt, Atm) * dt;
	*this = *this + 1.0/6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
	t += dt;
}

void Vector::Predict20(double tend)
{
	Atmosphera Atm;
	Atm.NMounth = DatTime.Get_Mounth() - 1;
	Atm.StarTime(DatTime);
	Atm.Sun(t);
	while (t < tend)
	{
		RKS20(20, Atm);
	}
	Perevod();
}