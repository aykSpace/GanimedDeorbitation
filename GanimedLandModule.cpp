#include "GanimedLandModule.h"
#include "OrbitElements.h"


GanimedLandModule::GanimedLandModule()
{
	GKa = 1600;
}


GanimedLandModule::~GanimedLandModule()
{}


void GanimedLandModule::RKSGan(Spacecraft & vect, double dt)
{
	Spacecraft K1, K2, K3, K4, result;
	K1 = FGanimed(vect, t) * dt;
	K2 = FGanimed(vect + 0.5*K1, t + 0.5*dt) * dt;
	K3 = FGanimed(vect + 0.5*K2, t + 0.5*dt) * dt;
	K4 = FGanimed(vect +     K3, t + dt) * dt;
	vect = vect + 1.0/6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
	vect.t += dt;
}

void GanimedLandModule::RKSGan(Spacecraft & vect, double dt, Spacecraft (*f_ptr)(const Spacecraft & vect, double dt))
{
	Spacecraft K1, K2, K3, K4, result;
	K1 = f_ptr(vect, t) * dt;
	K2 = f_ptr(vect + 0.5*K1, t + 0.5*dt) * dt;
	K3 = f_ptr(vect + 0.5*K2, t + 0.5*dt) * dt;
	K4 = f_ptr(vect +     K3, t + dt) * dt;
	vect = vect + 1.0/6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
	vect.t += dt;
}


Spacecraft FGanimed(const Spacecraft & vect, double dt)
{
	Spacecraft Rez (vect);
	Date_Time dateTimeEfem(1997, 01, 16, 00, 00, 00);
	OrbitElements oe1 (vect.vitok, dateTimeEfem.GetSecOfDay(), 1070400, 0.0013, 192.417 * PI / 180.0, 0.177 * M_PI / 180, 63.552 * M_PI / 180, 192.417, 1.016e-5, dateTimeEfem);
	Vector nuYup = oe1.KeplerToVector(vect.DatTime);
	double deltax = -nuYup.x - vect.x;
	double deltay = -nuYup.y - vect.y;
	double deltaz = -nuYup.z - vect.z;

	double deltar = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
	double rganSat = sqrt(vect.x*vect.x + vect.y*vect.y + vect.z*vect.z);
	double koef1 = rganSat*rganSat*rganSat;
	double koef2 = deltar*deltar*deltar;

	Rez.x = vect.Vx;
	Rez.y = vect.Vy;
	Rez.z = vect.Vz;

	Rez.Vx = -mu_g*vect.x/koef1 + mu_yu * ((-nuYup.x - vect.x)/koef2 + nuYup.x/koef2);
	Rez.Vy = -mu_g*vect.y/koef1 + mu_yu * ((-nuYup.y - vect.y)/koef2 + nuYup.y/koef2);
	Rez.Vz = -mu_g*vect.z/koef1 + mu_yu * ((-nuYup.z - vect.z)/koef2 + nuYup.z/koef2);

	double AxDu = 0;//компоненты ускорени€ от двигател€
	double AyDu = 0;
	double AzDu = 0;
	double AT = 0; //ускорение от двигател€ в Ќьютонах
	if (vect.EngineStart)
	{
		Rez.pnX = vect.PDu/vect.Mtek; //перегрузка
		AT  = Rez.pnX*0.001428;
		Rez.ImpSizeCur = AT;
		//Rez.P2[0][0] = -Sign(Rez.P2[0][0]); 
		//Rez.P2[2][0] = -Sign(Rez.P2[2][0]); 
		AxDu = Rez.P2[0][0]*AT;
		AyDu = Rez.P2[0][0]*AT;
		AzDu = Rez.P2[2][0]*AT;

		Rez.Vx -= AxDu;
		Rez.Vy -= AyDu;
		Rez.Vz -= AzDu;
	}
	Rez._ax = Rez.Vx;
	Rez._ay = Rez.Vy;
	Rez._az = Rez.Vz;

	return Rez;
}

Spacecraft FGanimedGSK(const Spacecraft & vect, double dt)
{
	Spacecraft Rez;
	Date_Time dateTimeEfem(1997, 01, 16, 00, 00, 00);
	OrbitElements oe1 (vect.vitok, dateTimeEfem.GetSecOfDay(), 1070400, 0.0013, 192.417 * PI / 180.0, 0.177 * M_PI / 180, 63.552 * M_PI / 180, 192.417, 1.016e-5, dateTimeEfem);
	Vector nuYup = oe1.KeplerToVectorGSK(vect.DatTime);
	double deltax = -nuYup.x - vect.x;
	double deltay = -nuYup.y - vect.y;
	double deltaz = -nuYup.z - vect.z;

	double deltar = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
	double rganSat = sqrt(vect.x*vect.x + vect.y*vect.y + vect.z*vect.z);
	double koef1 = rganSat*rganSat*rganSat;
	double koef2 = deltar*deltar*deltar;
	double koef3 = omega_g*omega_g;
	double koef4 = omega_g*2;
	Rez.x = vect.Vx;
	Rez.y = vect.Vy;
	Rez.z = vect.Vz;
	Rez.Vx = koef3*vect.x + koef4*vect.Vy + (-mu_g*vect.x/koef1 + mu_yu * ((-nuYup.x - vect.x)/koef2 + nuYup.x/koef2));
	Rez.Vy = koef3*vect.y - koef4*vect.Vx + (-mu_g*vect.y/koef1 + mu_yu * ((-nuYup.y - vect.y)/koef2 + nuYup.y/koef2));
	Rez.Vz = -mu_g*vect.z/koef1 + mu_yu * ((-nuYup.z - vect.z)/koef2 + nuYup.z/koef2);
	return Rez;
}


void GanimedLandModule::EngineZoneGanimed(Spacecraft & SC, vector<Spacecraft> &Rezult, double _Tvkl, unsigned curEngine)
{
	int i = 0; 
	unsigned j = 1; //дл€ корректного первого шага в InterpolMatrix (i = 1, вычисл€ем AM[j-1][3])
	unsigned OrderMatr = 8;
	double **AM;
	AM = new double * [OrderMatr];
	for (unsigned k = 0; k < OrderMatr; ++k)
		AM[k] = new double [OrderMatr];
	SC.EngineStart = true;
	SC.dm = SC.Engines[curEngine].Trust/SC.Engines[curEngine].SpecificImpulse;
	SC.PDu = SC.Engines[curEngine].Trust;
	double dt = 2;
	SC.ImpSizeCur = 0;
	SC.TVkl = _Tvkl;
	double vv;
	//SC.TurnZ = 60*PI/180;
	while (SC.ImpSizeCur < SC.ImpSizeNu/1000)
	{
		SC.OrientFisk(SC.t);
		SC.TRabDu = SC.t - SC.TVkl;
		SC.Mtek   = SC.GKa - SC.dm*SC.TRabDu;
		RKSGan(SC, dt, FGanimed);
		SC.PerevodGanimed();
		cout << SC.Vx <<"  " << SC.Vy <<"  " << SC.ModuleV() <<"  " << SC.h  << endl;
		 vv = SC.ModuleV();
		i += 1;
		InterpolMatrix(AM, OrderMatr, i, SC, j, SC.ImpSizeCur);
	}
	Ilagr(8, AM, 7, SC.ImpSizeNu/1000, SC, 7);
	SC.EngineStart = false;
	SC.DatTime.SetFrmtTime(SC.t);
	Rezult.push_back(SC);
	SC.EngineStart = false;
	for (unsigned i = 0; i < OrderMatr; ++i) // удаление матрицы
		delete [] AM [i];
	delete [] AM;
}

vector<Spacecraft> GanimedLandModule::BalDeorb()
{
	vector<Spacecraft> result(2);

	double M0 = 0;
	//double e = 0.0095;
	//double a = 2671.2;
	double e = 0.000045;
	double a = 2651.2;
	double omega = 0 * M_PI / 180.0;
	double i = 0 * M_PI / 2; //equator
	double omega1 = 0 * M_PI / 180.0;
	double n = sqrt(mu_g/(a*a*a));
	Date_Time dateTime(2015, 02, 20, 00, 00, 00);
	Date_Time dateTimeStart(1997, 1, 16, 00, 00, 00);
	Date_Time dateTimeFinish(2015, 8, 20, 00, 00, 00);

	OrbitElements oeGanimedSat (1, dateTimeStart.GetSecOfDay(), a, e, omega, i, omega1, M0, n,dateTimeStart);
	Vector vectorNuSat = oeGanimedSat.KeplerToVector(dateTime);
	GanimedLandModule nuSat(*this);
	nuSat.SetCurrentVector(vectorNuSat);
	oeGanimedSat = OrbitElements::Kepler(nuSat);
	nuSat.PerevodGanimed();

	nuSat.RStab = Stabilization::ISK;
	nuSat.SetDeorbitImpulse(1780.2);
	nuSat.DatTime.SetFrmtTime(nuSat.t);
	//EngineZoneGanimed(nuSat, result, nuSat.t);
	//nuSat.GKa = nuSat.Mtek;
	oeGanimedSat = OrbitElements::Kepler(nuSat);
	double tApog = nuSat.t + oeGanimedSat.T()/2;
	nuSat.PerevodGanimed();
	double hPred;

	do 
	{		
		hPred = nuSat.h;
		RKSGan(nuSat, 20, FGanimed);
		if (nuSat.t - COUNT_SEC_OF_DAY > 0)
		{	
			nuSat.DatTime.Calendar();
			nuSat.t -= COUNT_SEC_OF_DAY;
			nuSat.DatTime.SetFrmtTime(nuSat.t);
		}
		nuSat.PerevodGanimed();
		if ((oeGanimedSat.Hmax() - nuSat.h) < 0.5)
		{
			EngineZoneGanimed(nuSat, result, nuSat.t);
		}
		/*if ((nuSat.h - oeGanimedSat.Hmin()) < 1)
		{
			nuSat.SetDeorbitImpulse(60);
			nuSat.DatTime.SetFrmtTime(nuSat.t);
			nuSat.TurnZ = 0;
			EngineZoneGanimed(nuSat, result, nuSat.t);
		}*/
		 
		double vv = nuSat.ModuleV();
		oeGanimedSat = OrbitElements::Kepler(nuSat);
		nuSat.DatTime.PrintDate();
		cout << oeGanimedSat.Hmin() << "  " << nuSat.Vx << "  " << nuSat.Vy << "  " << nuSat.Vz << endl;
	} while (oeGanimedSat.Hmin() > 2);

	/*Vector vectorNuSat = oeGanimedSat.KeplerToVectorGSK(dateTime);
	GanimedLandModule nuSat(*this);
	nuSat.SetCurrentVector(vectorNuSat);
	nuSat.DatTime.SetFrmtTime(nuSat.t);
	do 
	{
		RKSGan(nuSat, 260, FGanimedGSK);
		nuSat.DatTime.SetFrmtTime(nuSat.t);
		if (nuSat.t - COUNT_SEC_OF_DAY > 0)
		{	
			nuSat.DatTime.Calendar();
			nuSat.t -= COUNT_SEC_OF_DAY;
			nuSat.DatTime.SetFrmtTime(nuSat.t);
		}
		nuSat.PerevodGanimed();
		nuSat.DatTime.PrintDate();
		cout << nuSat.h << endl;
	} while (nuSat.h > 20);*/
	return result;
}


void GanimedLandModule::SetCurrentVector(const Vector & vect)
{
	DatTime = vect.DatTime;
	x  = vect.x;
	y  = vect.y;
	z  = vect.z;
	Vx = vect.Vx;
	Vy = vect.Vy;
	Vz = vect.Vz;
	Sb = vect.Sb;
	vitok = vect.vitok;
	t = vect.t;

}