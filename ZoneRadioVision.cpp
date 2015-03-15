#include "ZoneRadioVision.h"
using namespace std;

ZoneRadioVision::ZoneRadioVision():sizeOfmassives(0), visible(false), nSutCircle(0)
{}

ZoneRadioVision::~ZoneRadioVision()
{
	if (this->sizeOfmassives != 0)
	{
		delete[] datTimes;
		delete[] distances;
		delete[] azimuths;
		delete[] cornersOfLocation;
		delete[] fis;
		delete[] lambdas;
		delete[] altitudes;
		delete[] cornersOfSun;
	}
}

ZoneRadioVision::ZoneRadioVision (Date_Time & _datTimes, double & _distances, double & _azimuths, double & _cornersOfLocation,
								  double & _fis, double & _lamdbas, double & _altitudes, double & _cornersOfSun)
								 :datTimes(&_datTimes), distances(&_distances), azimuths(&_azimuths), cornersOfLocation(&_cornersOfLocation),
								 fis(&_fis), lambdas(&_lamdbas), altitudes(&_altitudes), cornersOfSun(&_cornersOfSun)
{}

ZoneRadioVision::ZoneRadioVision(ZoneRadioVision const& _copy)
	:visible(_copy.visible), vitok(_copy.vitok), sizeOfmassives(_copy.sizeOfmassives), nSutCircle(_copy.nSutCircle)
{
	if (sizeOfmassives != 0)
	{
		datTimes = new Date_Time[sizeOfmassives];
		distances = new double[sizeOfmassives];
		azimuths = new double[sizeOfmassives];
		cornersOfLocation = new double[sizeOfmassives];
		fis = new double[sizeOfmassives];
		lambdas = new double[sizeOfmassives];
		altitudes = new double[sizeOfmassives];
		cornersOfSun = new double[sizeOfmassives];
		stdext::checked_array_iterator<Date_Time*> chd_dateTime (datTimes, _copy.sizeOfmassives);
		std::copy(_copy.datTimes, _copy.datTimes + _copy.sizeOfmassives,chd_dateTime);
		std::copy(_copy.distances, _copy.distances + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (distances,_copy.sizeOfmassives));
		std::copy(_copy.azimuths, _copy.azimuths + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (azimuths,_copy.sizeOfmassives));
		std::copy(_copy.cornersOfLocation, _copy.cornersOfLocation + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (cornersOfLocation,_copy.sizeOfmassives));
		std::copy(_copy.fis, _copy.fis + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (fis,_copy.sizeOfmassives));
		std::copy(_copy.lambdas, _copy.lambdas + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (lambdas,_copy.sizeOfmassives));
		std::copy(_copy.altitudes, _copy.altitudes + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (altitudes,_copy.sizeOfmassives));
		std::copy(_copy.cornersOfSun, _copy.cornersOfSun + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (cornersOfSun,_copy.sizeOfmassives));
	}
}

ZoneRadioVision & ZoneRadioVision::operator= (ZoneRadioVision & _other)
{
	Swap(_other);
	return *this;
}

void ZoneRadioVision::Swap(ZoneRadioVision & _other)
{
	visible = _other.visible;
	std::swap(vitok,_other.vitok);
	std::swap(sizeOfmassives,_other.sizeOfmassives);
	std::swap(nSutCircle,_other.nSutCircle);
	std::swap(datTimes,_other.datTimes);
	std::swap(distances,_other.distances);
	std::swap(azimuths,_other.azimuths);
	std::swap(cornersOfLocation,_other.cornersOfLocation);
	std::swap(fis,_other.fis);
	std::swap(lambdas,_other.lambdas);
	std::swap(altitudes,_other.altitudes);
	std::swap(cornersOfSun,_other.cornersOfSun);
}

vector<ZoneRadioVision> ZoneRadioVision::Calculate(const GroundStation & _grSt, const Vector &_vect, unsigned _numOfCirc, double _gamma1) //нужно в конце процедуры поменять sizeofmassives
{
	double fiRad  = _grSt.Fi * PerGradRad;
	double lamRad = _grSt.Lambda * PerGradRad;
	vector<ZoneRadioVision> zrvs;
	try
	{
		double * vectGrSt = new double[3];
		double ** matr = new double *[3];
		for (unsigned i = 0; i < 3; ++i)
			matr[i] = new double[3];
		double gamma0 = CalcVectorGrStMatrGamma(matr, fiRad, lamRad, _grSt, vectGrSt);
		double **am = new double *[8];
		for (unsigned i = 0; i < 8; i++)
			am[i] = new double[8];
		double ** matrTimDisAz = new double *[4];
		for (unsigned i = 0; i < 4; ++i)
			matrTimDisAz[i] = new double[4];
		Vector vect = _vect;
		int endIndication = _numOfCirc;
		ZoneRadioVision zrv;
		unsigned counter = 0;
		while (endIndication > 0)
		{
			zrv = Integration(matr,am,matrTimDisAz,vectGrSt,gamma0,vect,20, _gamma1, counter);
			endIndication -= 1;
			zrvs.push_back(zrv);
		}
		//-------------test zrv-------
		/*auto ttdt0 = zrvs[5].datTimes[0];
		auto cof0  = zrvs[5].cornersOfLocation[0];
		auto ttdt1 = zrvs[5].datTimes[1];
		auto cof1  = zrvs[5].cornersOfLocation[1];
		auto ttdt2 = zrvs[5].datTimes[2];
		auto cof2  = zrvs[5].cornersOfLocation[2];*/
		/*vector<Date_Time> tt = zrvs.back().GetDatTime();
		vector<double> tt1 = zrvs.back().GetDistance();
		vector<double> tt2 = zrvs.back().GetAzimuts();
		vector<double> tt3 = zrvs.back().GetCornersOfLocation();
		vector<double> tt4 = zrvs.back().GetFis();
		vector<double> tt5 = zrvs.back().GetLambdas();
		vector<double> tt6 = zrvs.back().GetAltitudes();
		vector<double> tt7 = zrvs.back().GetCornersOfSun();
		unsigned tt8 = zrvs.back().GetVitok();
		unsigned tt9 = zrvs.back().GetNumSutCirc();
		int tt10 = zrvs.back().GetSizeOfMassives();
		bool tt11 = zrvs.back().GetVisible();*/

		delete [] vectGrSt;
		for (unsigned i = 0; i < 3; ++i)
			delete[] matr[i];
		delete[] matr;
		for (unsigned i = 0; i < 8; ++i)
			delete[] am[i];
		delete[] am;
		for (unsigned i = 0; i < 4; ++i)
			delete[] matrTimDisAz[i];
		delete[] matrTimDisAz;
		return zrvs;
	}
	catch (std::bad_alloc err)
	{
		this->~ZoneRadioVision();
		throw err;
	}
	catch (...)
	{
		this->~ZoneRadioVision();
		throw;
	}
}

void ZoneRadioVision::VectorOfGrSt(double** _matrPer, double _fiGrSt, double _lamGrSt, double _hGrSt, double *_outVectGrSt)
{
	const double Rz = 6378.14;
	double sf = sin(_fiGrSt);
	double c  = Rz/sqrt(1.0 - KoefPressure*(2.0 - KoefPressure)*sf*sf);
	double cf = cos(_fiGrSt);
	double sl = sin(_lamGrSt);
	double cl = cos(_lamGrSt);
	double c1 = c + _hGrSt;
	_outVectGrSt[0] = c1*cf*cl;
	_outVectGrSt[1] = c1*cf*sl;
	_outVectGrSt[2] = (c*(1.0 - KoefPressure)*(1.0 - KoefPressure) + _hGrSt)*sf;
	_matrPer[0][0] = -sf*cl;
	_matrPer[0][1] = -sf*sl;
	_matrPer[0][2] = cf;
	_matrPer[1][0] = cf*cl;
	_matrPer[1][1] = cf*sl;
	_matrPer[1][2] = sf;
	_matrPer[2][0] = -sl;
	_matrPer[2][1] = cl;
	_matrPer[2][2] = 0.0;
}

void ZoneRadioVision::DisAziCor(double** _matrPer, const Vector &_vect, double *_vecGrSt, double &_distance, double &_azimuth, double &_cornerLoc)
{
	double  differenceOfVectors[3];
	for (int i = 0; i < 3; i++)
		differenceOfVectors[i] = _vect[i+1] - _vecGrSt[i];
	double xp[3]; //вектор разницы в системе координат НИП
	for (int i = 0; i < 3; i++)
	{
		xp[i] = 0;
		for (int j = 0; j < 3; j++)
			xp[i] += _matrPer[i][j]*differenceOfVectors[j];
	}
	_distance  = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2]);
	_cornerLoc = PerRadGrad*asin(xp[1]/_distance);
	_azimuth   = PerRadGrad*atan(xp[2]/xp[0]);
	if (Sign(xp[0]) < 0)
		_azimuth += 180.;
	if(_azimuth < 0)
		_azimuth += 360.;
}

double ZoneRadioVision::CalcVectorGrStMatrGamma( double ** _matrPer, double fiRad, double lamRad, const GroundStation &_grSt, double * vectGrSt )
{
	this->VectorOfGrSt(_matrPer,fiRad,lamRad,_grSt.HUnderSea,vectGrSt); 
	double rZ = 6378.14*(1 - KoefPressure*sin(fiRad)*sin(fiRad)) + _grSt.HUnderEarth; //радиус Земли с учетом сжатия
	double rN = sqrt(vectGrSt[0]*vectGrSt[0] + vectGrSt[1]*vectGrSt[1] + vectGrSt[2]*vectGrSt[2]);
	double gam0;
	fabs(_grSt.HUnderSea - _grSt.HUnderEarth) < 0.05 ? gam0 = 0 : gam0 = PerRadGrad*asin(rZ/rN) - 90.0; //определяем откуда начинается ЗРВ
	return gam0;
}

ZoneRadioVision ZoneRadioVision::Integration(double **_matrPer, double **_am, double **_matrTimeDisAzCor, double *_vectGrSt, double _gamma0, Vector &_vect, int _dt, double _gamma1, unsigned & _counter)
{
	ZoneRadioVision rezult;
	ZoneRadioVision resultCurent;
	resultCurent.InitializationMassives(3);
	resultCurent.vitok = static_cast<unsigned>(_vect[13]);
	double distance, azimuth, cornerOfLocation;
	this->DisAziCor(_matrPer,_vect, _vectGrSt,distance,azimuth,cornerOfLocation);
	_vect.Perevod();
	if (_counter != 0) resultCurent.nSutCircle = _vect.nSutCirc(); // если _vect на экваторе (не первое вхождение)
	bool _halfVitok = false;
	bool izmVitka = false;
	bool isProgn = true;
	int i = 0; unsigned j = 1;
	int jzv = 0;
	unsigned orderMatr = 8;
	_am[0][0] = _vect[0];
	_am[0][1] = _vect[1];
	_am[0][2] = _vect[2];
	_am[0][3] = _vect[3];
	_am[0][4] = _vect[4];
	_am[0][5] = _vect[5];
	_am[0][6] = _vect[6];
	_am[0][7] = 0;
	_matrTimeDisAzCor[0][2] = azimuth;
	Atmosphera atm;
	atm.NMounth = _vect.DatTime.Get_Mounth() - 1;
	atm.StarTime(_vect.DatTime);
	atm.Sun(_am[0][0]);
	double * timeDisAzCorner = new double[4];
	Vector vectZrv;
	unsigned counter = 0;
	visible = false;
	if (cornerOfLocation > _gamma0)  // если НУ на начало ЗРВ
	{
		counter = 2;
		for (unsigned i = 0; i < counter; i++)
		{
			resultCurent.distances[i]  = 0.0;
			resultCurent.datTimes [i]  = Date_Time(2000,1,1);
			resultCurent.azimuths [i]  = 0.0;
			resultCurent.cornersOfLocation [i] = 0.0;
			resultCurent.cornersOfSun [i] = 0.0;
			resultCurent.fis [i] = 0.0;
			resultCurent.lambdas [i] = 0.0;
			resultCurent.altitudes [i] = 0.0;
		}
	}
	auto startVisible = false;
	auto gammaCorrent = _gamma0;
	auto cornerCurrent = cornerOfLocation;
	while (!visible || !izmVitka)
	{
		double h0 = _vect.h;
		_vect.RKS(_dt, atm);
		i += 1;
		InterpolMatrix(_am, orderMatr, i, _vect, j, 0, & _halfVitok, &izmVitka, isProgn);
		this->DisAziCor(_matrPer,_vect,_vectGrSt, distance, azimuth, cornerOfLocation);
		_matrTimeDisAzCor[jzv][0] = _vect.t;
		_matrTimeDisAzCor[jzv][1] = distance;
		_matrTimeDisAzCor[jzv][2] = azimuth;
		_matrTimeDisAzCor[jzv][3] = cornerOfLocation;
		jzv +=1;
		jzv > 3 ? jzv = 0 : jzv;
		counter >= 1  ? cornerCurrent = -1.0*cornerOfLocation : cornerCurrent = cornerOfLocation;
		if (cornerCurrent >= _gamma0 && i >= 4 && !visible)
		{
			startVisible = true;
			resultCurent.visible = true;
			vectZrv = _vect;
			vectZrv.DatTime = _vect.DatTime;
			counter +=1;
			double cornerOfSun = vectZrv.CalcCornOfSun(atm);
			if (counter > 1)
			{
				GetResultZrv(resultCurent,_am, orderMatr, _matrTimeDisAzCor, vectZrv, counter, _gamma0, cornerOfSun);
				vectZrv.DatTime = _vect.DatTime;
			}
			if (counter == 1) //первое вхождение
			{
				vectZrv.DatTime = _vect.DatTime;
				GetResultZrv(resultCurent,_am, orderMatr, _matrTimeDisAzCor, vectZrv, counter, _gamma0, cornerOfSun);
				gammaCorrent = _gamma1;
			}
			if (counter >= 2)
			{
				visible = true;
				startVisible = false;
			}
		}
		if (cornerOfLocation >= _gamma1 && counter == 1)
		{
			vectZrv = _vect;
			vectZrv.DatTime = _vect.DatTime;
			counter +=1;
			double cornerOfSun = vectZrv.CalcCornOfSun(atm);
			GetResultZrv(resultCurent,_am, orderMatr, _matrTimeDisAzCor, vectZrv, counter, _gamma1, cornerOfSun);
		}
		if (azimuth < 90.0 || _matrTimeDisAzCor[jzv][2] > 270.0)
		{
			for (int i = 0; i < 4; i++)
				_matrTimeDisAzCor[i][2] -= 360;
		}
		if (izmVitka && i >=8 && !visible && !startVisible)
			break;
		if (_vect.t - COUNT_SEC_OF_DAY > 0)
		{	
			_vect.DatTime.Calendar();
			atm.StarTime(_vect.DatTime);
			_vect.t -= COUNT_SEC_OF_DAY;
		}
		if (_halfVitok)
		{
			atm.Sun(_vect.t);
			_halfVitok = false;
		}
	}
	_vect.Perevod();
	
	if (_counter == 0)
	{
		unsigned nSutc = _vect.nSutCirc() - 1;
		nSutc == 0 ? nSutc += 1 : nSutc;
		resultCurent.nSutCircle = nSutc;
	}
	_counter += 1;
	fromCourrentToResult(resultCurent,rezult,counter);
	delete[] timeDisAzCorner;
	return rezult;
}

void ZoneRadioVision::InitializationMassives(unsigned _sizeOfMassives)
{
	this->distances = new double[_sizeOfMassives];
	this->azimuths = new double[_sizeOfMassives];
	this->cornersOfLocation = new double[_sizeOfMassives];
	this->fis = new double[_sizeOfMassives];
	this->lambdas = new double[_sizeOfMassives];
	this->altitudes = new double[_sizeOfMassives];
	this->cornersOfSun = new double[_sizeOfMassives];
	this->datTimes = new Date_Time[_sizeOfMassives];
	sizeOfmassives = _sizeOfMassives;
}


//-----------------------------------getters------------------------------------

vector<Date_Time> ZoneRadioVision::GetDatTime()
{
	vector<Date_Time> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->datTimes[i]);
	}
	return rezult;
}

vector<double> ZoneRadioVision::GetDistance()
{
	vector<double> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->distances[i]);
	}
	return rezult;
}

vector<double> ZoneRadioVision::GetAzimuts()
{
	vector<double> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->azimuths[i]);
	}
	return rezult;
}

vector<double> ZoneRadioVision::GetCornersOfLocation()
{
	vector<double> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->cornersOfLocation[i]);
	}
	return rezult;
}

vector<double> ZoneRadioVision::GetFis()
{
	vector<double> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->fis[i]);
	}
	return rezult;
}

vector<double> ZoneRadioVision::GetLambdas()
{
	vector<double> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->lambdas[i]);
	}
	return rezult;
}

vector<double> ZoneRadioVision::GetAltitudes()
{
	vector<double> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->altitudes[i]);
	}
	return rezult;
}

vector<double> ZoneRadioVision::GetCornersOfSun()
{
	vector<double> rezult;
	if (sizeOfmassives == 0)
		return rezult;
	for (int i = 0; i < sizeOfmassives; ++i)
	{
		rezult.push_back(this->cornersOfSun[i]);
	}
	return rezult;
}

unsigned ZoneRadioVision::GetVitok()
{
	return vitok;
}

unsigned ZoneRadioVision::GetNumSutCirc()
{
	return nSutCircle;
}

int ZoneRadioVision::GetSizeOfMassives()
{
	return sizeOfmassives;
}

bool ZoneRadioVision::GetVisible()
{
	return visible;
}

void ZoneRadioVision::GetResultZrv(ZoneRadioVision & res, double ** am, unsigned orderAm, double** matrInterpol, Vector & vectRes, unsigned counter, double gammaInterpol, double cornerOfSun)
{
	if (res.sizeOfmassives != 0)
	{
		double * timeDisAzCorner = new double[4];
		Ilagr(4,matrInterpol,3,gammaInterpol,timeDisAzCorner,4); // timeDisAzCorner на момент достижения _gamma0
		if (timeDisAzCorner[2] < 0)
			timeDisAzCorner[2] += 360;
		Ilagr(orderAm, am, 0, timeDisAzCorner[0], vectRes,7); //вектор на  ЗРВ по времени
		vectRes.Perevod();
		vectRes.DatTime.SetFrmtTime(timeDisAzCorner[0]);
		res.distances[counter - 1]  = timeDisAzCorner[1];
		res.datTimes [counter - 1]  = vectRes.DatTime;
		res.azimuths [counter - 1]  = timeDisAzCorner[2];
		res.cornersOfSun [counter - 1] = cornerOfSun;
		res.fis [counter - 1] = vectRes[8];
		res.lambdas [counter - 1] = vectRes[9];
		res.altitudes [counter - 1] = vectRes[10];
		res.cornersOfLocation [counter - 1] = gammaInterpol;
		delete [] timeDisAzCorner;
	}
}

void ZoneRadioVision::fromCourrentToResult (ZoneRadioVision & current, ZoneRadioVision & result, unsigned counter)
{
	result.vitok = current.vitok;
	result.nSutCircle = current.nSutCircle;
	result.visible = current.visible;
	if (counter != 0)
	{
		result.InitializationMassives(counter);
		for (unsigned i = 0; i < counter; i++)
		{
			result.altitudes[i] = current.altitudes[i];
			result.azimuths[i] = current.azimuths[i];
			result.cornersOfLocation[i] = current.cornersOfLocation[i];
			result.cornersOfSun[i] = current.cornersOfSun[i];
			result.datTimes[i] = current.datTimes[i];
			result.distances[i] = current.distances[i];
			result.fis[i] = current.fis[i];
			result.lambdas[i] = current.lambdas[i];
		}
	}
}