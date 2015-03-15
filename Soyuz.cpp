#include "Soyuz.h"

Soyuz::Soyuz() : _lSa(2145.0), _hDivision(140.0), _vs(0), _s(0), _vsr(0), _gamma(0), _dgamma(0),
				 _guidenceStart(false), _gammaCom(0), _pxMax(0)
{}

Soyuz::Soyuz(const Spacecraft& Spcrafrt) : Spacecraft(Spcrafrt), _lSa(2145.0), _hDivision(140.0), _vs(0), _s(0), _vsr(0),
										   _gamma(0), _dgamma(0), _gammaCom(0), _guidenceStart(false), _pxMax(0) 
{}

Soyuz::~Soyuz()
{}


Soyuz::Soyuz (const Vector &Vect, double _GKa, double _GBo, double _GSa, double _ImpSizeNu, double _Xt, double _Yt, double _Zt, double _Sbal,
			  double _RazpImpSize, double _NKa, unsigned _ROr, unsigned _RStab)
			  : Spacecraft(Vect), _lSa(2145.0), _hDivision(140.0), _vs(0), _s(0), _vsr(0), _gamma(0), _dgamma(0), _gammaCom(0), _guidenceStart(false), _pxMax(0)
{
	GKa = _GKa;	
	GBo = _GBo; 
	GSa = _GSa;
	ImpSizeNu = _ImpSizeNu;
	ImpSizeCur = 0;
	Xt = _Xt;
	Yt =_Yt;
	Zt = _Zt;
	Sbal = _Sbal;
	RazpImpSize =_RazpImpSize;
	NKa =_NKa;
	ROr =_ROr;
	RStab = _RStab;
}


Soyuz  operator + (const Soyuz & a, const Soyuz & b)
{
	Soyuz c1 = static_cast<Spacecraft>(a) + static_cast<Spacecraft>(b);
	c1._gamma =  a._gamma + b._gamma;
	c1._dgamma = a._dgamma + b._dgamma;
	c1._gammaCom = a._gammaCom;
	c1.DCg[0] =  a.DCg[0];
	c1.DCg[1] =  a.DCg[1];
	c1.DCg[2] =  a.DCg[2];
	c1._vs    =  a._vs;
	c1._vsr   =  a._vsr;
	c1._s     =  a._s;
	c1._guidenceStart = a._guidenceStart;
	c1._tDivsion = a._tDivsion;
	c1._pxMax = a._pxMax;
	return(c1);
}

Soyuz operator - (const Soyuz & a, const Soyuz & b)
{
	Soyuz c1 = static_cast<Spacecraft>(a) - Spacecraft(b);
	return(c1);
}

Soyuz operator * (double a, const Soyuz & b)
{
	Soyuz c1  =  a * Spacecraft(b);
	c1._gamma =  a * b._gamma;
	c1._dgamma = a * b._dgamma;
	c1._gammaCom = b._gammaCom;
	c1._vs    = b._vs;
	c1._vsr   = b._vsr;
	c1._s     = b._s;
	return(c1);
}


Soyuz operator * (const Soyuz & b, double a)
{
	Soyuz c1 = Spacecraft(b) * a;
	c1._gamma =  a * b._gamma;
	c1._dgamma = a * b._dgamma;
	c1._gammaCom = b._gammaCom;
	c1._vs    = b._vs;
	c1._vsr   = b._vsr;
	c1._s     = b._s;
	return(c1);
}


void Soyuz::RKS(double dt, Atmosphera &Atm, Soyuz (*f_ptr) (Soyuz & soyuz, double t , Atmosphera & atm))
{
	Soyuz K1, K2, K3, K4;
	K1 = f_ptr(*this, t, Atm) * dt;
	K2 = f_ptr(*this + 0.5*K1, t + 0.5*dt, Atm) * dt;
	K3 = f_ptr(*this + 0.5*K2, t + 0.5*dt, Atm) * dt;
	K4 = f_ptr(*this +     K3, t + dt, Atm) * dt;
	*this = *this + 1.0/6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
	t += dt;
}

Soyuz FGuideDeorb (Soyuz & SC, double t, Atmosphera &Atm)
{
	Soyuz Rez;
	const double Rzsr  = 6371.0;
	const double TAI20 = -0.10183014444;
	const double omz2  = 53.174954e-10; 
	const double A100  = 62.5648249885;
	const double domz  = 14.58423171e-5;
	const double Rekv  = 6378.14;
	const double szat  = 0.003352824419;
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

	double Koef1, Koef2;
	GetDensity(Vector::f,Vector::ap,Atm,SC); //рассчет плотности атмосферы
	/*Rez.x = SC.Vx;
	Rez.y = SC.Vy;
	Rez.z = SC.Vz;*/
	Koef1 = OMEGA_EARTH * OMEGA_EARTH;
	Koef2 = OMEGA_EARTH * 2.0;
	/*double *RezGpz = new double [3];
	gpz(SC, RezGpz); // рассчет коэффициентов гравитационного поля Земли
	Rez.Vx =  Koef2*SC.Vy + Koef1*SC.x + RezGpz[0];
	Rez.Vy = -Koef2*SC.Vx + Koef1*SC.y + RezGpz[1];
	Rez.Vz =  RezGpz[2];*/
	double VB[3] = {SC.Vx - Atm.W[0], SC.Vy - Atm.W[1], SC.Vz - Atm.W[2]};
	double VB2 = VB[0]*VB[0] + VB[1]*VB[1] + VB[2]*VB[2];
	double V = sqrt(VB2);
	double max;
	SC.h >= 120 ? max = 30.0 : max = V/Atm.VZV;
	double s = 0;
	double s1 = 0;
	double b = 0;
	SC.Aerodynamic(max, s, s1, b);
	double p[3][3];
	Mavg(SC.t, SC._tDivsion, p);
	double DTek[3];
	double BA[3];
	double CA[3];
	//перевод в ГСК из скоростной СК 
	DTek[0] = SC.DCg[0]*p[0][0] + SC.DCg[1]*p[0][1];
	DTek[1] = SC.DCg[0]*p[1][0] + SC.DCg[1]*p[1][1];
	DTek[2] = SC.DCg[2];
	//текущая вертикаль (приборная - отсчитывается от момента разделения)
	BA[0] = DTek[1]*VB[2] - DTek[2]*VB[1];
	BA[1] = DTek[2]*VB[0] - DTek[0]*VB[2];
	BA[2] = DTek[0]*VB[1] - DTek[1]*VB[0];
	double c = 1/sqrt(BA[0]*BA[0] + BA[1]*BA[1] + BA[2]*BA[2]);
	for (int i = 0; i < 3; i++)
		BA[i] *= c;
	//текущий "бок"
	CA[0] = BA[2]*VB[1] - BA[1]*VB[2];
	CA[1] = BA[0]*VB[2] - BA[2]*VB[0];
	CA[2] = BA[1]*VB[0] - BA[0]*VB[1];
	c = 1/sqrt(CA[0]*CA[0] + CA[1]*CA[1] + CA[2]*CA[2]);
	double cb = b * SC.ro*VB2;
	double gr = SC._gamma*PerGradRad;
	double sGr = c*sin(gr);
	double cGr = cos(gr);
	Rez.Vx = Rez.Vx + cb*(cGr*BA[0] + sGr*CA[0]);
	Rez.Vy = Rez.Vy + cb*(cGr*BA[1] + sGr*CA[1]);
	Rez.Vz = Rez.Vz + cb*(cGr*BA[2] + sGr*CA[2]);
	if (SC._guidenceStart) //начало управления
	{
		Rez._gamma = SC._dgamma; //интегрирование угла крена
		double dg = SC._gamma - SC._gammaCom;
		if (fabs(dg) > 26.0)
			dg = 26*Sign(dg);
		double u = dg + 2*SC._dgamma;
		double u1, aac, bac, ugu;		
		if (SC.pnX > 4.2)
		{
			u1  = 8.0;
			aac = 0.18625;
			bac = 0.69;
			ugu = 11.459156;
		}
		else
		{
			u1  = 6.0;
			aac = 0.3725;
			bac = 1.435;
			ugu = 5.729578;
		}
		double et = 1.0;
		if (u <= -u1)
			et = -1;
		if (u <= -4 && u > -u1)
			et = aac*u + bac;
		if (u < 4 && u > -4)
		{
			double dac = 0.01375;
			et = dac*u;
		}
		if (u < u1 && u >= 4)
			et = aac*u - bac;
		Rez._dgamma = -et*ugu;
	}
	double c0 = SC.ro*V;
	double cs = -s*c0;
	Rez.Vx += cs*VB[0];
	Rez.Vy += cs*VB[1];
	Rez.Vz += cs*VB[2];
	Rez.ImpSizeCur = s1*SC.ro*VB2; //интегрируется (ImpCurSize) - кажущаяся скорость
	Rez.pnX = Rez.ImpSizeCur*101.9716213; //не интегрируется (не прописано в операторах)
	//delete [] RezGpz;
	return Rez;
}

void Soyuz::Aerodynamic(double max, double & s, double & s1, double & b)
{
	double yt1 = Yt/cos(atan2(Zt,fabs(Yt)));
	double dh, dcx, dcy, dmz, dcxa, dcya, dmza, ct, cn, alb;		
	if (max > 6.0)
	{
		//case 3 page 33
		if (h > 90.0)
		{
			//table 21
			int * th1 = new int[5];
			th1[0] = 90; th1[1] = 100; th1[2] = 110; th1[3] = 120; th1[4] = 130;
			double * tdcx1 = new double[5];
			tdcx1[0] = 0.068; tdcx1[1] = 0.50; tdcx1[2] = 0.93; tdcx1[3] = 1.01; tdcx1[4] = 1.12;
			double * tdcy1 = new double[5];
			tdcy1[0] = 0.088; tdcy1[1] = 0.259; tdcy1[2] = 0.51; tdcy1[3] = 0.578; tdcy1[4] = 0.67;
			double * tdmz1 = new double[5];
			tdmz1[0] = -0.017; tdmz1[1] = -0.0432; tdmz1[2] =-0.0968; tdmz1[3] = -0.1152; tdmz1[4] = -0.14;
			double * tdcxa1 = new double[5];
			tdcxa1[0] = 0.0026; tdcxa1[1] = -0.0104; tdcxa1[2] =-0.0145; tdcxa1[3] = -0.0132; tdcxa1[4] =-0.016;
			double * tdcya1 = new double[5];
			tdcya1[0] = 0.004; tdcya1[1] = 0.013; tdcya1[2] = 0.016; tdcya1[3] = 0.0195; tdcya1[4] = 0.018;
			double * tdmza1 = new double[5];
			tdmza1[0] = 0.00034; tdmza1[1] = 0.0007; tdmza1[2] = 0.0027; tdmza1[3] = 0.0032; tdmza1[4] = 0.00389;
			double h1 = h;
			if (h1 > 130.0)
				h1 = 130.0;
			int ih = poi(3, th1, h1);
			int ih1 = ih + 1;
			dh = (h1 - th1[ih])/(th1[ih1] - th1[ih]);
			dcx = tdcx1[ih] + (tdcx1[ih1] - tdcx1[ih])*dh;
			dcy = tdcy1[ih] + (tdcy1[ih1] - tdcy1[ih])*dh;
			dmz = tdmz1[ih] + (tdmz1[ih1] - tdmz1[ih])*dh;
			dcxa = tdcxa1[ih] + (tdcxa1[ih1] - tdcxa1[ih])*dh;
			dcya = tdcya1[ih] + (tdcya1[ih1] - tdcya1[ih])*dh;
			dmza = tdmza1[ih] + (tdmza1[ih1] - tdmza1[ih])*dh;
			delete [] th1;
			delete [] tdcx1;
			delete [] tdcy1;
			delete [] tdmz1;
			delete [] tdcxa1;
			delete [] tdcya1;
			delete [] tdmza1;
			th1 = 0; tdcx1 = 0; tdcy1 = 0; tdmz1 = 0; tdcxa1 = tdcya1 = tdmza1 = 0;
		}
		//case 2 page 31
		if (max > 6.0 && h <= 90)
		{
			int * th2 = new int[7];
			th2[0] = 30; th2[1] = 50; th2[2] = 60; th2[3] = 65; th2[4] = 70; th2[5] = 75; th2[6] = 90;
			int * tm1 = new int[5];
			tm1[0] = 6; tm1[1] = 15; tm1[2] = 20; tm1[3] = 25; tm1[4] = 30;

			unsigned nRow = 5; unsigned nCol = 7;
			double ** tdCx = GetTdCx();
			double ** tdCy = GetTdCy();
			double ** tdMz = GetTdMz();

			int ih = poi(6,th2,h);
			int im = poi(4, tm1, max);
			int ih1 = ih + 1;
			int im1 = im + 1;
			dh = (h - th2[ih])/(th2[ih1] - th2[ih]);
			double dm = (max - tm1[im])/(tm1[im1] - tm1[im]);

			dcx = lin2Matr(nRow,nCol,tdCx,im, ih, im1, ih1, dh, dm);
			dcy = lin2Matr(nRow,nCol,tdCy,im, ih, im1, ih1, dh, dm);
			dmz = lin2Matr(nRow,nCol,tdMz,im, ih, im1, ih1, dh, dm);

			double** tdCxa = GetTdCxa();
			double** tdCya = GetTdCya();
			double** tdMza = GetTdMza();

			dcxa = lin2Matr(nRow,nCol,tdCxa, im, ih, im1, ih1, dh, dm);
			dcya = lin2Matr(nRow,nCol,tdCya, im, ih, im1, ih1, dh, dm);
			dmza = lin2Matr(nRow,nCol,tdMza, im, ih, im1, ih1, dh, dm);

			DeleteMatr(nRow, nCol, tdCx);
			DeleteMatr(nRow, nCol, tdCy);
			DeleteMatr(nRow, nCol, tdMz);
			DeleteMatr(nRow, nCol, tdCxa);
			DeleteMatr(nRow, nCol, tdCya);
			DeleteMatr(nRow, nCol, tdMza);

			delete [] tm1; 
			delete [] th2;
		}
		//alfaBal cx, cy
		double ambs = 0.1013 - 0.14*Xt;
		double amna = 0.005228 - 0.00967*Xt;
		double dal  = -(ambs + 0.0023 + dmz + yt1*(1.33 + dcx))/(amna + dmza + yt1*(-0.0155 + dcxa));
		alb  = 20.0 + dal;
		cn   = 0.14 + 0.009*dal + dcy + dcya*dal;
		ct   = 1.33 - 0.0155*dal + dcx + dcxa*dal;
	}
	//case 1 max <= 6
	if (max <= 6.0)
	{
		double * txt = new double[3];
		txt[0] = 0.37; txt[1] = 0.38; txt[2] = 0.39;
		double * tyt = new double[5];
		tyt[0] = -0.04662; tyt[1] = -0.041958; tyt[2] = -0.037296; tyt[3] = -0.034965; tyt[4] = -0.032634;
		double * tm2 = new double[6];
		tm2[0] = 0.6; tm2[1] = 0.8; tm2[2] = 1.1; tm2[3] = 1.7; tm2[4] = 4.0; tm2[5] = 6.0;
		int * ta = new int[8];
		ta[0] = 0; ta[1] = 5; ta[2] = 10; ta[3] = 15; ta[4] = 20; ta[5] = 25; ta[6] = 30; ta[7] = 35;  

		int ix = poi(2, txt, Xt);
		int iy = poi(4, tyt, yt1);
		int im = poi(5, tm2, max);
		int ix1 = ix + 1;
		int iy1 = iy + 1;
		int im1 = im + 1;
		double dx = (Xt - txt[ix])/(txt[ix1] - txt[ix]);
		double dy = (yt1 - tyt[iy])/(tyt[iy1] - tyt[iy]);
		double dm = (max - tm2[im])/(tm2[im1] - tm2[im]);

		double *** tal = GetTal();
		alb = lin3Matr(3, 5, 6, tal,ix, iy, im, ix1, iy1, im1, dx, dy, dm);
		int ia = poi(7, ta, alb);
		int ia1 = ia + 1;
		double da = (alb - ta[ia])/(ta[ia1] - ta[ia]);
		double ** tCx2 = GetTCx2();
		double ** tCy2 = GetTCy2();
		ct = lin2Matr(8, 6, tCx2, ia, im, ia1, im1, dm, da);
		cn = lin2Matr(8, 6, tCy2, ia, im, ia1, im1, dm, da);

		DeleteMatr(8,6,tCy2);
		DeleteMatr(8,6, tCx2);
		DeleteMatr3(3, 5, 6, tal);
		delete[] txt;
		delete[] tyt;
		delete[] tm2;
		delete[] ta;
	}
	double alr = alb * PerGradRad;
	double sa = sin(alr);
	double ca = cos(alr);
	double cx = ct*ca + cn*sa;
	double cy = ct*sa - cn*ca;
	double smid = 0.0000038;
	double cg = 0.004903325e12*smid/GSa; //0.00490... - g/2
	s = cx * cg;
	b = cy * cg;
	s1 = sqrt(cn*cn + ct*ct)*cg;
}

vector<Spacecraft> Soyuz::ControlledDeorb(const unsigned long _Vitok)
{
	GuidDeorb = true;
	vector<Vector> arxOfPrognoz = this->Prognoz(static_cast<unsigned int>(_Vitok));
	double ayt = YtTma();
	double fiAim, lAim;
	AimPointContDeorb(this->lambda, fiAim, lAim);
	vector<Spacecraft> RezultOfCalculations; //(0 - SC на момент вкл ДУ, ... , last index - SC на h=10.7
	ShootGuideDeorbit(0,lAim, RezultOfCalculations); //first iteration
	double lambdasLanding[10];
	double tEngineStarting[10];
	double dtEngineStart;
	unsigned nIter = 0;
	RezultOfCalculations.front().TVkl > 0 ? tEngineStarting[nIter] = RezultOfCalculations.front().TVkl : tEngineStarting[nIter] = RezultOfCalculations.front().TVkl + 86400;
	lambdasLanding[nIter] = RezultOfCalculations.back().lambda;
	nIter += 1;
	double DLPos = 0.2;
	do 
	{
		double filanding = RezultOfCalculations.back().fi;
		dtEngineStart = this->MorePreciseTVkl(lAim, *this, RezultOfCalculations, tEngineStarting, lambdasLanding, nIter, DLPos);
		if (fabs(dtEngineStart) <= 0.0001)
		 {
			double dfi = 0.07;
			double dFiLanding =  filanding - fiAim;
			if (fabs(dFiLanding ) < dfi)
				break;
			double dVSR = (1.0 - 2.0*_s)*dFiLanding*24.7;
			int vsr1 = _vsr + 1*(int)dVSR;
			if (vsr1 < 0)
				vsr1 = 0;
			if (_s == 1 && dFiLanding < 0.0 && _vsr == 15)
				throw logic_error("Side maneuver is small - check aim attitude!");
			if (_s == 0 && dFiLanding > 0.0 && _vsr == 15)
				throw logic_error("Side maneuver is small - check aim attitude!");
			if (_s == 1 && dFiLanding < 0.0 && _vsr > 15)
				vsr1 = 15;
			if (_s == 0 && dFiLanding > 0.0 && _vsr == 15)
				vsr1 = 15;
			_vsr = vsr1;            
		}
		tEngineStarting[nIter] = tEngineStarting[nIter - 1] + dtEngineStart;
		tEngineStarting[nIter] > 0 ? tEngineStarting[nIter] : tEngineStarting[nIter] += 86400;
		RezultOfCalculations.clear();
		ShootGuideDeorbit(tEngineStarting[nIter], lAim, RezultOfCalculations);
		lambdasLanding[nIter] = RezultOfCalculations.back().lambda;
		nIter += 1;
	} 
	while (nIter < 10);
	return RezultOfCalculations;
}

void Soyuz::AimPointContDeorb(double longitudeEquator, double & fiAim, double & lAim)
{
	const double longEquatorTab[19] = {-55.0, -53.0, -51.2, -48.5, -45.9, -42.9, -39.4, -33.8, -27.5, 
		-18.0, -10.0, -1.5, 0.5, 4.0, 6.3, 8.6, 12.7, 16.5, 20.0};
	const double fiTab[18] = {45.6, 46.12, 46.9, 47.83, 47.33, 48.45, 49.92, 50.67, 51.33, 51.33, 51.0, 50.67, 49.92, 50.35, 49.52, 48.45, 47.33, 46.20};
	const double longTab[18] = {65.0, 65.17, 65.60, 65.25, 69.58, 69.17, 66.95, 67.33, 67.58, 67.58, 67.17, 67.33, 66.95, 70.87, 68.73, 69.17, 69.58, 69.50};
	const double hTab[18] = {0.15, 0.20, 0.22, 0.20, 0.35, 0.40, 0.40, 0.35, 0.33, 0.33, 0.30, 0.35, 0.40, 0.45, 0.50, 0.40, 0.35, 0.30};
	const int tub[18][2] = {
		{-6 , 11},
		{-2 , 11},
		{-9 , 11},
		{-7 , 11},
		{-8 , 12},
		{-10, 12},
		{-11, 13},
		{-12, 8 },
		{-1 , 0 },
		{8  , -3},
		{11 ,-11},
		{-2 ,-10},
		{5  ,-10},
		{-1 ,-13},
		{-2 ,-13},
		{13, -12},
		{13, -12},
		{12, -14}
	};
	int i = poi(18, longEquatorTab, longitudeEquator);
	if (fabs(longitudeEquator - longEquatorTab[i]) < 1e-8 && i != 0)
		i -= 1;
	fiAim = fiTab[i];
	lAim = longTab[i];
	//double hPol = hTab[i];
	int i1 = i + 1;
	double dl = (longitudeEquator - longEquatorTab[i])/(longEquatorTab[i1] - longEquatorTab[i]);
	double as = tub[i][0] + (tub[i][1] - tub[i][0])*dl;
	_s = static_cast<unsigned>((1.0 - Sign(as))*0.5);
	double c = fabs(as) + 1e-10;
	_vsr = 1 * static_cast<int>(c);
	if (abs(c - _vsr) > 0.5)
		_vsr += 1;
}

void Soyuz::ShootGuideDeorbit(const double tVkl, const double aimLambda, vector<Spacecraft> &result) const
{
	Soyuz soyuz = EngineZone(*this, aimLambda, result, tVkl);
	soyuz._s = _s;
	soyuz._vsr = _vsr;
	soyuz.PredictGuideLandingPoint(result);
}

void Soyuz::PredictGuideLandingPoint(vector<Spacecraft> &result)
{	
	Atmosphera Atm;
	Atm.NMounth = DatTime.Get_Mounth() - 1;
	int i = 0; 
	unsigned j = 1;
	unsigned OrderMatr = 9;
	double **AM;
	AM = new double * [OrderMatr];
	for (unsigned k = 0; k < OrderMatr; ++k)
		AM[k] = new double [OrderMatr];
	Perevod();
	double dt = 2;
	while (h >= _hDivision)
	{
		Atm.StarTime(DatTime);
		Atm.Sun(t);
		Spacecraft::RKS(dt, Atm);
		Altitude();
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, *this, j, this->h);	
	}
	Ilagr(OrderMatr, AM, 7, _hDivision, *this, 8); //интерполяция на момент разделения
	//debugging
	/*x = 5291.30417146308; y = 2203.82458112624; z = 3093.30752788590; Vx = -4.58666886933886; Vy = 3.63738571875574; Vz = 4.79982630795501; t = 85699.0772687980;
	_tDivsion = t;*/
	Perevod();
	_tDivsion = t;
	result.push_back(*this);
	double pxi = pnX;
//---------------------интегрирование после разделения-------------------------------------------------------
	Sb = Sbal;
	double vxa = Vx - OMEGA_EARTH*y;
	double vya = Vy + OMEGA_EARTH*x;
	DCg[0] = vya*z - Vz*y;
	DCg[1] = Vz*x  - vxa*z;
	DCg[2] = vxa*y - vya*x;
	double c = 1/sqrt(DCg[0]*DCg[0] + DCg[1]*DCg[1] + DCg[2]*DCg[2]);
	for (int i = 0; i < 3; i++)
		DCg[i] *= c;
	while (r >= 6471.0)
	{
		if (h > 120)
		{
			Atm.StarTime(DatTime);
			Atm.Sun(t);
		}
		RKS(dt, Atm, FGuideDeorb);
		_vs = _vs + GEarth*(pnX + pxi);
		pxi = pnX;
		i += 1;
		Perevod();
		InterpolMatrix(AM, OrderMatr, i, *this, j, r, _vs);
	}
	double vsrev = 2.0 + _vsr*0.147; //кажущаяся скорость на момент переворота
	double resultInterpolVector[9];
	//------------------------------------ радиус входа (начало атмосферного участка)------------------------------------------
	Ilagr(OrderMatr, AM, 7, 6471.0, resultInterpolVector, 9); 
	RewriteVector(resultInterpolVector);
	_vs = resultInterpolVector[8];
	Perevod();
	result.push_back(*this);
	dt = 0.2;
	ImpSizeCur = _vs;
	double pxMaxCurrent = 0;
	double trv = t; //время расчета вариаций плотности (для оптимизации вычислений) 
	while (ImpSizeCur <= 0.0256)
	{
		if (t < trv)
			calcVarAtm = false;
		else
		{
			calcVarAtm = true;
			trv = t + 2;
		}
		RKS(dt,Atm, FGuideDeorb);
		i += 1;
		Altitude();
		InterpolMatrix(AM, OrderMatr, i, *this, j, ImpSizeCur, h);
		if (pnX > pxMaxCurrent)
			pxMaxCurrent = pnX;
	}
	//-----------------------------------------начало управления (Ws = 25.6)-----------------------------------------------
	Ilagr(OrderMatr, AM, 7, 0.0256, resultInterpolVector, 9); 
	RewriteVector(resultInterpolVector);
	h = resultInterpolVector[8];
	ImpSizeCur = resultInterpolVector[7];
	Perevod();
	result.push_back(*this);
	_gammaCom = (1.0 - 2.0*_s)*60.0; //расчет командного угла крена
	_guidenceStart = true;
	double tStartGuidance = resultInterpolVector[0];
	double vsi = 0.2; 
	double dvsi = 0.2; // 200m/c дискрет для кажущейся скорости 
	unsigned ivs = 0;
	const unsigned ivsk = 37; //37 дискретов управления
	double tvsi[ivsk];
	bool isGammaCom = false; //выход на gammaCom
	while (ImpSizeCur <= vsrev)
	{
		if (t < trv)
			calcVarAtm = false;
		else
		{
			calcVarAtm = true;
			trv = t + 2;
		}
		RKS(dt,Atm, FGuideDeorb);
		i += 1;
		Altitude();
		InterpolMatrix(AM, OrderMatr, i, *this, j, ImpSizeCur, h);
		if (fabs(_gamma) >= 60 && !isGammaCom)
		{
			isGammaCom = true;
			_gamma = _gammaCom;
			_dgamma = 0;
			_guidenceStart = false;
			dt = 1;
		}
		if (pnX > pxMaxCurrent)
			pxMaxCurrent = pnX; //максимальная перегрузка
		if (ImpSizeCur >= vsi)
		{
			Ilagr(OrderMatr, AM, 7, vsi, resultInterpolVector, 1); //1 - результат только для первого члена (заполняем времена опорной траектории)
			tvsi[ivs] = resultInterpolVector[0] - tStartGuidance;
			Date_Time::CheckDaySeconds(tvsi[ivs]);
			vsi += dvsi;
			ivs += 1;
		}
	}
	//------------------------------------------------момент переворота----------------------------------------------------------------------
	Ilagr(OrderMatr, AM, 7, vsrev, resultInterpolVector, 9);
	RewriteVector(resultInterpolVector);
	h = resultInterpolVector[8];
	ImpSizeCur = resultInterpolVector[7];
	Perevod();
	result.push_back(*this);
	isGammaCom = false;
	_guidenceStart = true;
	_gammaCom = -(1.0 - 2.0*_s)*60.0;
	dt = 0.2;
	while (h >= 10.7)
	{
		if (t < trv)
			calcVarAtm = false;
		else
		{
			calcVarAtm = true;
			trv = t + 2;
		}
		RKS(dt,Atm, FGuideDeorb);
		i += 1;
		Altitude();
		if (ivs <= 37)
			InterpolMatrix(AM, OrderMatr, i, *this, j, ImpSizeCur, h);
		else
			InterpolMatrix(AM, OrderMatr, i, *this, j, h, ImpSizeCur);
		if (_pxMax == 0)
			FindPxMax(pxMaxCurrent, pnX);
		if (fabs(_gamma) >= 60 && !isGammaCom)
		{
			isGammaCom = true;
			_gamma = _gammaCom;
			_dgamma = 0;
			_guidenceStart = false;
			dt = 1;
		}
		if (ImpSizeCur >= vsi && ivs <= ivsk)
		{
			if (ivs != ivsk)
			{
				Ilagr(OrderMatr, AM, 7, vsi, resultInterpolVector, 1); //1 - результат только для первого члена (заполняем времена опорной траектории)
				tvsi[ivs] = resultInterpolVector[0] - tStartGuidance;
				Date_Time::CheckDaySeconds(tvsi[ivs]);
				vsi += dvsi;
				ivs += 1;
			}
			else
			{
				Ilagr(OrderMatr, AM, 7, vsi, resultInterpolVector, 9);
				RewriteVector(resultInterpolVector);
				h = resultInterpolVector[8];
				ImpSizeCur = resultInterpolVector[7];
				ivs += 1;
			}
			
		}
	}
	//--------------------------------высота ввода ОСП----------------------------------------------------------------
	Ilagr(OrderMatr, AM, 7, 10.7, resultInterpolVector, 9);
	RewriteVector(resultInterpolVector);
	Perevod();
	result.push_back(*this);	
}

double Soyuz::YtTma()
{
	const double tyt[8] = {-85.0, -88.3, -91.7, -95.0, -81.7, -78.3, -75.0, -105.0};
	double yt1 = Yt * _lSa;
	double dmin = 1e10;
	int imin = 0;
	for (int i = 0; i < 8; i++)
	{
		double d = fabs(yt1 - tyt[i]);
		if (d < dmin)
		{
			dmin = d;
			imin = i;
		}
	}
	return 1.0*(imin - 1);

}

//table 16
double ** Soyuz::GetTdCx(){
	double **tdCx = new double * [5];
	for (unsigned k = 0; k < 5; ++k)
	{
		tdCx[k] = new double [7];
	}
	tdCx[0][0] = tdCx[0][1] = tdCx[0][2] = tdCx[0][3] = tdCx[0][4] = tdCx[0][5] = tdCx[0][6] = 0.0; 
	tdCx[1][0] = 0.0370; tdCx[1][1] = 0.0350; tdCx[1][2] = 0.033; tdCx[1][3] = 0.033; tdCx[1][4] = 0.0385; tdCx[1][5] = 0.043; tdCx[1][6] = 0.068;
	tdCx[2][0] = 0.0385; tdCx[2][1] = 0.0380; tdCx[2][2] = 0.035; tdCx[2][3] = 0.035; tdCx[2][4] = 0.0395; tdCx[2][5] = 0.044; tdCx[2][6] = 0.068;
	tdCx[3][0] = 0.0318; tdCx[3][1] = 0.0318; tdCx[3][2] = 0.029; tdCx[3][3] = 0.030; tdCx[3][4] = 0.0365; tdCx[3][5] = 0.043; tdCx[3][6] = 0.068;
	tdCx[4][0] = 0.0220; tdCx[4][1] = 0.0210; tdCx[4][2] = 0.024; tdCx[4][3] = 0.028; tdCx[4][4] = 0.0330; tdCx[4][5] = 0.040; tdCx[4][6] = 0.068;
	return tdCx;
}

//table 15
double ** Soyuz::GetTdCy(){
	double **tdCy = new double * [5];
	for (unsigned k = 0; k < 5; ++k)
	{
		tdCy[k] = new double [7];
	}
	tdCy[0][0] = tdCy[0][1] = tdCy[0][2] = tdCy[0][3] = tdCy[0][4] = tdCy[0][5] = tdCy[0][6] = 0.0; 
	tdCy[1][0] = -0.011; tdCy[1][1] = -0.0100; tdCy[1][2] = -0.006; tdCy[1][3] = 0.0000; tdCy[1][4] = 0.01000; tdCy[1][5] = 0.025; tdCy[1][6] = 0.088;
	tdCy[2][0] = -0.019; tdCy[2][1] = -0.0195; tdCy[2][2] = -0.014; tdCy[2][3] = -0.008; tdCy[2][4] = 0.00400; tdCy[2][5] = 0.020; tdCy[2][6] = 0.088;
	tdCy[3][0] = -0.025; tdCy[3][1] = -0.0260; tdCy[3][2] = -0.025; tdCy[3][3] = -0.022; tdCy[3][4] = -0.0080; tdCy[3][5] = 0.012; tdCy[3][6] = 0.088;
	tdCy[4][0] = -0.030; tdCy[4][1] = -0.0325; tdCy[4][2] = -0.034; tdCy[4][3] = -0.031; tdCy[4][4] = -0.021; tdCy[4][5] = -0.002; tdCy[4][6] = 0.088;
	return tdCy;
}

//table 19
double **Soyuz::GetTdCxa()
{
	double **tdCxa = new double *[5];
	for (unsigned k= 0; k<5; ++k)
	{
		tdCxa[k] = new double [7];
	}
	tdCxa[0][0] =           tdCxa[0][1] =           tdCxa[0][2] =           tdCxa[0][3] =           tdCxa[0][4] =           tdCxa[0][5] =           tdCxa[0][6] = 0.0;
	tdCxa[1][0] = -0.00150; tdCxa[1][1] = -0.00150; tdCxa[1][2] = -0.00176; tdCxa[1][3] = -0.00175; tdCxa[1][4] = -0.00155; tdCxa[1][5] = -0.00095; tdCxa[1][6] = 0.0026;
	tdCxa[2][0] = -0.00126; tdCxa[2][1] = -0.00130; tdCxa[2][2] = -0.00155; tdCxa[2][3] = -0.00165; tdCxa[2][4] = -0.00150; tdCxa[2][5] = -0.00095; tdCxa[2][6] = 0.0026;
	tdCxa[3][0] = -0.00175; tdCxa[3][1] = -0.00190; tdCxa[3][2] = -0.00210; tdCxa[3][3] = -0.00195; tdCxa[3][4] = -0.00160; tdCxa[3][5] = -0.00095; tdCxa[3][6] = 0.0026;
	tdCxa[4][0] = -0.00235; tdCxa[4][1] = -0.00240; tdCxa[4][2] = -0.00250; tdCxa[4][3] = -0.00230; tdCxa[4][4] = -0.00180; tdCxa[4][5] = -0.00115; tdCxa[4][6] = 0.0026;    
	return tdCxa;
}

//table 18
double **Soyuz::GetTdCya()
{
	double **tdCya = new double *[5];
	for (unsigned k= 0; k<5; ++k)
	{
		tdCya[k] = new double [7];
	}
	tdCya[0][0] =           tdCya[0][1] =           tdCya[0][2] =           tdCya[0][3] =           tdCya[0][4] =           tdCya[0][5] =           tdCya[0][6] = 0.0;
	tdCya[1][0] = -0.00075; tdCya[1][1] = -0.00075; tdCya[1][2] = -0.00065; tdCya[1][3] = -0.00050; tdCya[1][4] = -0.00025; tdCya[1][5] = -0.00030; tdCya[1][6] = 0.004;
	tdCya[2][0] = -0.00105; tdCya[2][1] = -0.00105; tdCya[2][2] = -0.00092; tdCya[2][3] = -0.00070; tdCya[2][4] = -0.00040; tdCya[2][5] = -0.00015; tdCya[2][6] = 0.004;
	tdCya[3][0] = -0.00128; tdCya[3][1] = -0.00125; tdCya[3][2] = -0.00115; tdCya[3][3] = -0.00105; tdCya[3][4] = -0.00075; tdCya[3][5] = -0.00015; tdCya[3][6] = 0.004;
	tdCya[4][0] = -0.00140; tdCya[4][1] = -0.00135; tdCya[4][2] = -0.00130; tdCya[4][3] = -0.00120; tdCya[4][4] = -0.00100; tdCya[4][5] = -0.00040; tdCya[4][6] = 0.004; 
	return tdCya;
}

//table 20
double **Soyuz::GetTdMza()
{
	double **tdCza = new double *[5];
	for (unsigned k= 0; k<5; ++k)
	{
		tdCza[k] = new double [7];
	}
	tdCza[0][0] =          tdCza[0][1] =          tdCza[0][2] =           tdCza[0][3] =           tdCza[0][4] =           tdCza[0][5] =           tdCza[0][6] = 0.0;
	tdCza[1][0] = 0.000145; tdCza[1][1] = 0.00015; tdCza[1][2] = 0.000155; tdCza[1][3] = 0.000115; tdCza[1][4] = 0.000075; tdCza[1][5] = 0.000055; tdCza[1][6] = 0.00034;
	tdCza[2][0] = 0.000200; tdCza[2][1] = 0.00021; tdCza[2][2] = 0.000220; tdCza[2][3] = 0.000185; tdCza[2][4] = 0.000125; tdCza[2][5] = 0.000085; tdCza[2][6] = 0.00034;
	tdCza[3][0] = 0.000250; tdCza[3][1] = 0.00026; tdCza[3][2] = 0.000270; tdCza[3][3] = 0.000240; tdCza[3][4] = 0.000130; tdCza[3][5] = 0.000130; tdCza[3][6] = 0.00034;
	tdCza[4][0] = 0.000295; tdCza[4][1] = 0.00031; tdCza[4][2] = 0.000320; tdCza[4][3] = 0.000295; tdCza[4][4] = 0.000175; tdCza[4][5] = 0.000170; tdCza[4][6] = 0.00034; 
	return tdCza;
}
// table 8 - 13
double *** Soyuz::GetTal()
{
	double *** tal;
	tal = new double **[3];
	for (int i = 0; i < 3; i++)
		tal[i] = new double *[5];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 5; j++)
			tal[i][j] = new double [6];
	}
	tal[0][0][0] = 19.35; tal[0][0][1] = 19.60; tal[0][0][2] = 26.52; tal[0][0][3] = 24.70; tal[0][0][4] = 24.60; tal[0][0][5] = 23.66;
	tal[0][1][0] = 17.72; tal[0][1][1] = 18.05; tal[0][1][2] = 23.33; tal[0][1][3] = 22.70; tal[0][1][4] = 22.17; tal[0][1][5] = 21.45;
	tal[0][2][0] = 16.05; tal[0][2][1] = 16.47; tal[0][2][2] = 20.48; tal[0][2][3] = 20.62; tal[0][2][4] = 19.67; tal[0][2][5] = 19.09;
	tal[0][3][0] = 15.21; tal[0][3][1] = 15.67; tal[0][3][2] = 19.36; tal[0][3][3] = 19.58; tal[0][3][4] = 18.56; tal[0][3][5] = 17.85;
	tal[0][4][0] = 14.26; tal[0][4][1] = 14.84; tal[0][4][2] = 18.40; tal[0][4][3] = 18.57; tal[0][4][4] = 17.43; tal[0][4][5] = 16.58;

	tal[1][0][0] = 19.04; tal[1][0][1] = 19.50; tal[1][0][2] = 27.61; tal[1][0][3] = 25.86; tal[1][0][4] = 25.46; tal[1][0][5] = 24.32;
	tal[1][1][0] = 17.40; tal[1][1][1] = 17.93; tal[1][1][2] = 23.65; tal[1][1][3] = 23.37; tal[1][1][4] = 22.89; tal[1][1][5] = 22.04;
	tal[1][2][0] = 15.72; tal[1][2][1] = 16.34; tal[1][2][2] = 20.59; tal[1][2][3] = 21.15; tal[1][2][4] = 20.21; tal[1][2][5] = 19.62;
	tal[1][3][0] = 14.85; tal[1][3][1] = 15.52; tal[1][3][2] = 19.39; tal[1][3][3] = 20.01; tal[1][3][4] = 19.02; tal[1][3][5] = 18.34;
	tal[1][4][0] = 13.90; tal[1][4][1] = 14.68; tal[1][4][2] = 18.41; tal[1][4][3] = 18.94; tal[1][4][4] = 17.85; tal[1][4][5] = 17.03;

	tal[2][0][0] = 18.72; tal[2][0][1] = 19.40; tal[2][0][2] = 29.39; tal[2][0][3] = 27.57; tal[2][0][4] = 26.50; tal[2][0][5] = 25.04;
	tal[2][1][0] = 17.07; tal[2][1][1] = 17.81; tal[2][1][2] = 24.03; tal[2][1][3] = 24.16; tal[2][1][4] = 23.67; tal[2][1][5] = 22.68;
	tal[2][2][0] = 15.39; tal[2][2][1] = 16.20; tal[2][2][2] = 20.72; tal[2][2][3] = 21.78; tal[2][2][4] = 20.88; tal[2][2][5] = 20.19;
	tal[2][3][0] = 14.48; tal[2][3][1] = 15.37; tal[2][3][2] = 19.42; tal[2][3][3] = 20.54; tal[2][3][4] = 19.52; tal[2][3][5] = 18.87;
	tal[2][4][0] = 13.54; tal[2][4][1] = 14.52; tal[2][4][2] = 18.42; tal[2][4][3] = 19.36; tal[2][4][4] = 18.31; tal[2][4][5] = 17.52;     
	return tal;

}

//table 14
double ** Soyuz::GetTCx2()
{
	double **tCx2 = new double * [8];
	for (unsigned k = 0; k < 8; ++k)
	{
		tCx2[k] = new double [6];
	}
	tCx2[0][0] = 1.07;  tCx2[0][1] = 1.150; tCx2[0][2] = 1.340; tCx2[0][3] = 1.510; tCx2[0][4] = 1.460; tCx2[0][5] = 1.46; 
	tCx2[1][0] = 1.065; tCx2[1][1] = 1.145; tCx2[1][2] = 1.340; tCx2[1][3] = 1.505; tCx2[1][4] = 1.445; tCx2[1][5] = 1.46;
	tCx2[2][0] = 1.050; tCx2[2][1] = 1.130; tCx2[2][2] = 1.340; tCx2[2][3] = 1.490; tCx2[2][4] = 1.420; tCx2[2][5] = 1.44;
	tCx2[3][0] = 1.030; tCx2[3][1] = 1.110; tCx2[3][2] = 1.330; tCx2[3][3] = 1.450; tCx2[3][4] = 1.380; tCx2[3][5] = 1.40;
	tCx2[4][0] = 1.000; tCx2[4][1] = 1.075; tCx2[4][2] = 1.305; tCx2[4][3] = 1.410; tCx2[4][4] = 1.320; tCx2[4][5] = 1.34;
	tCx2[5][0] = 0.960; tCx2[5][1] = 1.030; tCx2[5][2] = 1.270; tCx2[5][3] = 1.350; tCx2[5][4] = 1.245; tCx2[5][5] = 1.26;
	tCx2[6][0] = 0.910; tCx2[6][1] = 0.975; tCx2[6][2] = 1.220; tCx2[6][3] = 1.280; tCx2[6][4] = 1.160; tCx2[6][5] = 1.17;
	tCx2[7][0] = 0.840; tCx2[7][1] = 0.905; tCx2[7][2] = 1.150; tCx2[7][3] = 1.190; tCx2[7][4] = 1.060; tCx2[7][5] = 1.07;
	return tCx2;
}

//table 17
double ** Soyuz::GetTdMz(){
	double **tdMz = new double * [5];
	for (unsigned k = 0; k < 5; ++k)
	{
		tdMz[k] = new double [7];
	}
	tdMz[0][0] = tdMz[0][1] = tdMz[0][2] = tdMz[0][3] = tdMz[0][4] = tdMz[0][5] = tdMz[0][6] = 0.0; 
	tdMz[1][0] = 0.00410; tdMz[1][1] = 0.00450; tdMz[1][2] = 0.00460; tdMz[1][3] = 0.00390; tdMz[1][4] = 0.00215; tdMz[1][5] = -0.0008; tdMz[1][6] = -0.017;
	tdMz[2][0] = 0.00520; tdMz[2][1] = 0.00560; tdMz[2][2] = 0.00600; tdMz[2][3] = 0.00530; tdMz[2][4] = 0.00355; tdMz[2][5] = 0.00070; tdMz[2][6] = -0.017;
	tdMz[3][0] = 0.00615; tdMz[3][1] = 0.00630; tdMz[3][2] = 0.00715; tdMz[3][3] = 0.00640; tdMz[3][4] = 0.00470; tdMz[3][5] = 0.00170; tdMz[3][6] = -0.017;
	tdMz[4][0] = 0.00667; tdMz[4][1] = 0.00715; tdMz[4][2] = 0.00805; tdMz[4][3] = 0.00725; tdMz[4][4] = 0.00550; tdMz[4][5] = 0.00230; tdMz[4][6] = -0.017;
	return tdMz;
}

//table 14
double ** Soyuz::GetTCy2()
{
	double **tCy2 = new double * [8];
	for (unsigned k = 0; k < 8; ++k)
	{
		tCy2[k] = new double [6];
	}
	tCy2[0][0] = 0.0000; tCy2[0][1] = 0.0000; tCy2[0][2] = 0.0000; tCy2[0][3] = 0.0000; tCy2[0][4] = 0.000; tCy2[0][5] = 0.000; 
	tCy2[1][0] = -0.040; tCy2[1][1] = -0.025; tCy2[1][2] = -0.025; tCy2[1][3] = -0.030; tCy2[1][4] = 0.020; tCy2[1][5] = 0.030;
	tCy2[2][0] = -0.075; tCy2[2][1] = -0.045; tCy2[2][2] = -0.030; tCy2[2][3] = -0.010; tCy2[2][4] = 0.050; tCy2[2][5] = 0.058;
	tCy2[3][0] = -0.095; tCy2[3][1] = -0.048; tCy2[3][2] = -0.022; tCy2[3][3] =  0.050; tCy2[3][4] = 0.085; tCy2[3][5] = 0.095;
	tCy2[4][0] = -0.090; tCy2[4][1] = -0.032; tCy2[4][2] = -0.015; tCy2[4][3] =  0.140; tCy2[4][4] = 0.140; tCy2[4][5] = 0.140;
	tCy2[5][0] = -0.063; tCy2[5][1] = -0.005; tCy2[5][2] =  0.090; tCy2[5][3] =  0.250; tCy2[5][4] = 0.195; tCy2[5][5] = 0.190;
	tCy2[6][0] = -0.020; tCy2[6][1] =  0.045; tCy2[6][2] =  0.240; tCy2[6][3] =  0.360; tCy2[6][4] = 0.270; tCy2[6][5] = 0.250;
	tCy2[7][0] =  0.055; tCy2[7][1] =  0.105; tCy2[7][2] =  0.400; tCy2[7][3] =  0.480; tCy2[7][4] = 0.345; tCy2[7][5] = 0.325; 
	return tCy2;
}

void Soyuz::DeleteMatr(unsigned nRow, unsigned nCol, double ** matr)
{
	for (unsigned i = 0; i < nRow; ++i) // удаление матрицы
	{
		delete [] matr [i];
	}
	delete [] matr;
}

void Soyuz::DeleteMatr3(unsigned i, unsigned j, unsigned k, double *** matr)
{
	for (unsigned i1 = 0; i1 < i; i1++)
	{
		for (unsigned j1 = 0; j1 < j; j1++)
			delete[] matr[i1][j1];
	}
	for (unsigned i1 = 0; i1 < i; i1++)
	{
		delete[] matr[i1];
	}
	delete[] matr;
}



int func1(int i){return i+1;}
int func2(int i){return i+5;}

int Soyuz::testRef(int(*func)(int))
{
	return func(1);
}

int Soyuz::TestRef()
{
	return this->testRef(func2);
}

void InterpolMatrix (double **AM, unsigned n, int i,  Soyuz & soyuz, unsigned &j, double arg, double arg1)
{
	if (fmod(float(i),float(n)) != 0)
	{
		AM [j][0] = soyuz.t;
		AM [j][1] = soyuz.x;
		AM [j][2] = soyuz.y;
		AM [j][3] = soyuz.z;
		AM [j][4] = soyuz.Vx;
		AM [j][5] = soyuz.Vy;
		AM [j][6] = soyuz.Vz;
		AM [j][7] = arg;
		AM [j][8] = arg1;
		j += 1;
	}
	else
	{
		j = 0;
		AM [j][0] = soyuz.t;
		AM [j][1] = soyuz.x;
		AM [j][2] = soyuz.y;
		AM [j][3] = soyuz.z;
		AM [j][4] = soyuz.Vx;
		AM [j][5] = soyuz.Vy;
		AM [j][6] = soyuz.Vz;
		AM [j][7] = arg;
		AM [j][8] = arg1;
		j += 1;
	}
}

void Soyuz::RewriteVector(double * vec)
{
	t = vec[0]; x = vec[1]; y = vec[2]; z = vec[3]; Vx = vec[4];
	Vy = vec[5]; Vz= vec[6];
}

void Soyuz::FindPxMax(double & currentPxmax, double px)
{
	if (px > currentPxmax)
		currentPxmax = px;
	double px1 = currentPxmax - 0.05;
	if (px < px1)
		_pxMax = currentPxmax;
}