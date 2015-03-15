#include <iostream>
#include <cmath>
#include "OrbitElements.h"
#include "BALCONST.H"

using namespace std;

OrbitElements :: OrbitElements()
	:_vitok(0), _date(0), _t(0), _a(0), _e(0), _omega(0), _i(0), _omega1(0), _u(0), _M(0), _E(0), _S(0)
{}
OrbitElements :: OrbitElements(unsigned long vitok, double date, double t, double a, double e, double omega, double i, double omega1, double u, double M, double E, double S)
	:_vitok(vitok),_date(date),_t(t),_a(a),_e(e),_omega(omega),_i(i),_omega1(omega1),_u(u),_M(M),_E(E),_S(S)
{}
OrbitElements :: OrbitElements(unsigned long vitok, double t, double a, double e, double omega, double i, double omega1, double M, double n, Date_Time dateTime)
	:_vitok(vitok), _t(t), _a(a), _e(e), _omega(omega), _i(i), _omega1(omega1), _M(M), _n(n), _dateTime(dateTime) 
{}
OrbitElements :: OrbitElements(Vector vect)
{
	_vect=vect;
}


OrbitElements :: ~OrbitElements(){}


 double & OrbitElements :: operator [] (int i)
	{
		switch (i)
		{
		case 0: return  (double&)_vitok;
		case 1: return  _date;
		case 2: return  _t;
		case 3: return  _a;
		case 4: return  _e;
		case 5: return  _omega;
		case 6: return  _i;
		case 7: return  _omega1;
		case 8: return  _u;
		case 9: return  _M;
		case 10: return _E;
		case 11: return _S;
		default:			
			throw std::logic_error("Index not found");
		}
	}

void	 OrbitElements :: Rdat(int mdat[3])
	// Распаковка даты в формате ЧЧММГГГГ в дату-массив: число, месяц, год
{
	
	mdat[2]=static_cast<int>(fmod(_vect.DatTime.DateFormat(),10000.));
	mdat[1]=static_cast<int>(fmod(_vect.DatTime.DateFormat(),1000000.));
	mdat[0]=Iround((_vect.DatTime.DateFormat()-mdat[1])/1000000.);
	mdat[1]=Iround((mdat[1]-mdat[2])/10000.);
	return;
}
int		 OrbitElements :: Dat84(const int id[3], const int id1[3], int &cyt)
	// Функция счета числа суток между двумя исходными календарными датами
	// (старый вариант)
{
	int jt,jt1;
	double ft;
	Dj(jt,ft,id[2],id[1],id[0],0.);
	Dj(jt1,ft,id1[2],id1[1],id1[0],0.);
	cyt=jt1-jt;
	return cyt;
}
int		 OrbitElements :: Dat84(const int id[3], const int id1[3])
	// Функция счета числа суток между двумя исходными календарными датами
{
	int jt,jt1,cyt;
	double ft;
	Dj(jt,ft,id[2],id[1],id[0],0.);
	Dj(jt1,ft,id1[2],id1[1],id1[0],0.);
	cyt=jt1-jt;
	return cyt;
}
double	 OrbitElements :: Rtime(double tzip)
	// Преобразование времени в формате ЧЧММСС.ДДД в число секунд
{
	double h,m;
	h=floor(tzip/10000.);
	tzip-=h*10000.;
	m=floor(tzip/100.);
	tzip-=m*100.;
	return (h*60.+m)*60.+tzip;
}
double	 OrbitElements :: RtimeZip (double t)
	// Преобразование времени в секундах в формат ЧЧММСС.ДДД
{
	double h,m,s;
	s=t-(h=floor(t/3600.))*3600.;
	s-=(m=floor(s/60.))*60.;
	return h*10000.+m*100.+s;
}
double	 OrbitElements :: RdatZip(int mdat[3])
{
	return mdat[0]*1000000+mdat[1]*10000+mdat[2];	 
}

void	 OrbitElements :: Ucopy(double a[], double b[], int n)
	// Пересылка элементов из одного массива в другой
{
	int i;
	for (i=0;i<n;i++)
		b[i]=a[i];
	return;
}
int		 OrbitElements :: Iround(double _x)
	// Округление
{
	int r;
	r=(int)((_x-floor(_x)<0.5)?floor(_x):floor(_x)+1);
	return r;
}
double	 OrbitElements :: Angle(double _x, double _y)
{
	double a;
	a=atan2(_x,_y);
	if (a<0.) a+=M_PI*2.;
	return a;
}

double	 OrbitElements :: StarTime(double sut, double t0)
	// Звездное время
{
	double t,omega,eps,time,stime,ret;
	t=sut/36525.;
	omega=(933059.79-6962911.23*t+7.48*t*t+0.0080*t*t*t)*4.848136811e-6;
	eps=-17.23*sin(omega);
	time=23925.836+8640184.542*t+0.0929*t*t+0.061164*eps;
	time=fmod(time,86400.);
	stime=time*7.27220522e-5;
	ret=stime+(t0-10800.)*omega_g;
	ret=fmod(ret,M_PI*2.);
	if (ret==0.) ret=M_PI*2.;
	return ret;
}

void	 OrbitElements :: Dj(int &jt, double &ft, int iy, int im, int id, double d)
	// Вход: iy, im, id - год, месяц, день; d - доли суток от 0 часов.
	// Выход: jt, ft - количество дней юлианского периода и доля суток.
{
	int j1dn,it;
	j1dn=static_cast<int>((im*3057)/100+id+((5-im/3)/5)*(2-(4-fmod(iy,4.))/4
		+(100-fmod(iy,100.))/100-(400-fmod(iy,400.))/400)
		+1721028+(iy*365+iy/4-iy/100+iy/400));
	jt=j1dn;
	ft=d-0.5;
	it=static_cast<int>(ft);
	if (ft<0.) it--;
	jt+=it;
	ft-=it;
	return;
}
void	 OrbitElements :: Jd(int jt, double ft, int &iy, int &im, int &id, double &d)
	// Новая, непроверенная процедура !!! Перевод с Фортрана (Glyba)
	// Вход: jt, ft - количество дней юлианского периода и доля суток.
	// Выход: iy, im, id - год, месяц, день; d - доли суток от 0 часов.
{
	int ii,j1d,j1,j2,j3;
	ii=static_cast<int>(ft+0.5);
	if (ft<-0.5)
		ii--;
	d=ft+0.5-ii;
	j1d=jt+ii-1721119;
	j1=(4*j1d-1)/146097;
	j2=(4*j1d-1-146097*j1)/4;
	j1d=(4*j2+3)/1461;
	j2=(4*j2+7-1461*j1d)/4;
	j3=(5*j2-3)/153;
	id=(5*j2+2-153*j3)/5;
	iy=100*j1+j1d;
	if (j3<10)
		im=j3+3;
	else
		im=j3-9;
	iy++;
	return;
}

void	 OrbitElements :: Nua( double nu2[6], double s, char pp)
	// Перевод вектора состояния из относительной системы координат
	// в абсолютную и обратно
	// s - звездное время
{
	double sn,cs;
	sn=sin(s);
	cs=cos(s);
	switch (pp)
	{
	case 1: // из отн. в абс.
			//nu2[3]=nu1[3]*cs-nu1[4]*sn;
		nu2[3]=_vect[1]*cs - _vect[2]*sn;
		nu2[4]=_vect[1]*sn + _vect[2]*cs;
		nu2[5]=_vect[3];
		nu2[0]=_vect[4]*cs - _vect[5]*sn - nu2[4]*omega_e;
		nu2[1]=_vect[4]*sn + _vect[5]*cs + nu2[3]*omega_e;
		nu2[2]=_vect[6];
		break;
	case 2: // из абс. в отн.
		nu2[3]=_vect[1]*cs + _vect[2]*sn;
		nu2[4]=-_vect[1]*sn + _vect[2]*cs;
		nu2[5]=_vect[3];
		cout <<"nu2[5]="<<nu2[5]<<endl;
		nu2[0]=_vect[4]*cs + _vect[5]*sn + nu2[4]*omega_e;
		nu2[1]=-_vect[4]*sn + _vect[5]*cs - nu2[3]*omega_e;
		cout <<"nu2[1]="<<nu2[1]<<endl;	
		nu2[2]=_vect[6];
		cout <<"nu2[2]="<<nu2[2]<<endl;
		break;
	default:
		// Ошибка
		break;
	}
	return;
}
void	 OrbitElements :: Nua(double nu1[11],double nu2[11],char pp)
	// Перевод вектора состояния из относительной системы координат
	// в абсолютную и обратно
{
	int dat[3],ct;
	double sut,zvr;
	Rdat(/*nu1[1],*/dat);
	Dat84(n75,dat,ct);
	sut=ct+d75;
	zvr=StarTime(sut,nu1[2]);
	Nua(/*&nu1[3],*/&nu2[3],zvr,pp);
	nu2[0]=nu1[0];
	nu2[1]=nu1[1];
	nu2[2]=nu1[2];
	nu2[9]=nu1[9];
	nu2[10]=nu1[10];
	return;
}
void	 OrbitElements :: Rvosk(double nu2[6])
	// Расчет оскулирующих элементов орбиты
	// по заданному вектору состояния в АСК
{
	double r,v2,cx,cy,cz,c,p,rv1,sn,cs,p1,f,t,s,w1;
	r=sqrt(nu2[3]*nu2[3]+nu2[4]*nu2[4]+nu2[5]*nu2[5]);
	v2=nu2[0]*nu2[0]+nu2[1]*nu2[1]+nu2[2]*nu2[2];
	cx=nu2[4]*nu2[2]-nu2[5]*nu2[1];
	cy=nu2[5]*nu2[0]-nu2[3]*nu2[2];
	cz=nu2[3]*nu2[1]-nu2[4]*nu2[0];
	c=cx*cx+cy*cy+cz*cz;
	p=c/mu_e;
	c=sqrt(c);
	_a = mu_e/(2.*mu_e/r-v2);
	_e = sqrt(1.-p/_a);
	rv1=nu2[3]*nu2[0]+nu2[4]*nu2[1]+nu2[5]*nu2[2];
	_i=Angle(sqrt(cx*cx+cy*cy),cz);
	sn=cx/c;
	cs=-cy/c;
	_omega1=Angle(cx,-cy);
	_u=Angle(nu2[5],sn*nu2[4]+cs*nu2[3]);
	p1=p/r;
	f=Angle(rv1*p1/c,p1-1.);
	t=r/_a;
	s=t*sin(f)/sqrt(p/_a);
	_E=Angle(s,t*cos(f)+_e);
	_M=_E-_e*s;
	w1=_u-f;
	_omega=w1;
	
	if (w1<0.) _omega+=M_PI*2.;


	
	return;
}
void	 OrbitElements :: OskRv(OrbitElements &osk, double nu[6])
	// Переход от оскулирующих элементов орбиты
	// к вектору состояния в относительной ГСК
{
	double p,c,f,sn,cs,r,vr,vu,c1,c2,c3,s1,s2,s3,g1,g2,g3,g4,g5;
	p=osk[3]*(1.-osk[4]*osk[4]);
	c=sqrt(mu_e*p);
	f=osk[8]-osk[5];
	sn=sin(f);
	cs=cos(f);
	r=p/(1.+osk[4]*cs);
	vr=c*osk[4]*sn/p;
	vu=c/r;
	c1=cos(osk[7]);
	c2=cos(osk[8]);
	c3=cos(osk[6]);
	s1=sin(osk[7]);
	s2=sin(osk[8]);
	s3=sin(osk[6]);
	g1=c1*c2-s1*s2*c3;
	g2=s1*c2+c1*s2*c3;
	g3=c1*s2+s1*c2*c3;
	g4=s1*s2-c1*c2*c3;
	g5=s2*s3;
	nu[0]=vr*g1-vu*g3;
	nu[1]=vr*g2-vu*g4;
	nu[2]=vr*g5+vu*c2*s3;
	nu[3]=r*g1;
	nu[4]=r*g2;
	nu[5]=r*g5;
	return;
}
void	 OrbitElements :: OskNu(OrbitElements &os, double ny[11])
	// Перевод оскулирующих элементов в вектор в ГСК
{
	double os1[6],nu[6],zvr,sut;
	int dat1[3],ct;
	OskRv(os,nu);
	Rdat(/*os[1],*/dat1);
	Dat84(n75,dat1,ct);
	sut=ct+d75;
	zvr=StarTime(sut,os[2]);
	Nua(/*nu,*/os1,zvr,2);
	Ucopy(os1,&ny[3],6);
	ny[0]=os[0];
	ny[1]=os[1];
	ny[2]=os[2];
	ny[9]=os[10];
	ny[10]=os[11];
	return;
}
void	 OrbitElements :: GetElements()
	// ГСК отн. -> элементы
{
	int dat[3],ct;
	double sut,zvr,nu2[6];
	Rdat(dat);
	Dat84(n75,dat,ct);
	sut=ct+d75;
	zvr=StarTime(sut, _vect[0]);
	Nua(/*nu1,*/nu2,zvr,1);
	Rvosk(nu2);
	_vitok=static_cast<unsigned long>(_vect[13]);
	_date=_vect.DatTime.DateFormat();
	_t=_vect.t;
	_S=_vect[7];
	return;
}


OrbitElements OrbitElements::Kepler(const Vector& nu)
{
	OrbitElements res;
	res._vitok = nu.vitok;
	double w1, w2, w3, r, v2, rpi, ralf, rab, nn, rv;
	double E;
	//double rom, fi2, lambda;
	double alfa[5];
	double betta[5];
	double c1 = nu.y * nu.Vz - nu.z * nu.Vy;
	double c2 = nu.z * nu.Vx - nu.x * nu.Vz;
	double c3 = nu.x * nu.Vy - nu.y * nu.Vx;
	double c = sqrt(c1*c1 + c2*c2 + c3*c3);

	w1 = nu.Vy*c3 - nu.Vz*c2;
	w2 = nu.Vz*c1 - nu.Vx*c3;
	w3 = nu.Vx*c2 - nu.Vy*c1;
	r = sqrt(nu.x*nu.x + nu.y*nu.y + nu.z*nu.z);
	double f1 = w1 - mu_g * nu.x / r;
	double f2 = w2 - mu_g * nu.y / r;
	double f3 = w3 - mu_g * nu.z / r;
	double f = sqrt(f1*f1 + f2*f2 + f3*f3);

	res._e = f/mu_g;

	v2   = nu.Vx*nu.Vx + nu.Vy*nu.Vy + nu.Vz*nu.Vz;
	res._h   = v2 - 2.0*mu_g/r;
	res._a   = -mu_g/res._h;
	rpi  = res._a * (1.0 - res._e);
	ralf = res._a*(1.0 + res._e);
	res._hMin = rpi - a_g;
	res._hMax = ralf - a_g;

	rab = sqrt(mu_g/res._a);
	nn = rab /res._a;
	rv = nu.x*nu.Vx + nu.y*nu.Vy + nu.z*nu.Vz;
	alfa[0] = rab * (rv / f);
	betta[0] = (r - res._a) / f * res._h;

	E   = OrbitElements::asin2(alfa[0], betta[0]);
	Date_Time datTime = nu.DatTime;
	double timeSec = datTime.GetSecOfDay();
	res._tau = timeSec - 1.0 / nn * (E - res._e * alfa[0]);
	res._T = 2 * PI / nn;

	betta[1] = c3 / c;
	alfa[1] = sqrt(1.0 - betta[1] * betta[1]);
	res._i = OrbitElements::asin2(alfa[1], betta[1]);

	if (res._i <=  PI / 10800.0)
	{
		alfa[2] = 0.0;
		betta[2] = 0.0;
		res._omega1 = 0.0;
	}
	else
	{
		alfa[2] = c1 / (c * alfa[1]);
		betta[2] = -c2 / (c * alfa[1]);
		res._omega1 = OrbitElements::asin2(alfa[2], betta[2]);
	}

	if (fabs(betta[1]) < 0.1 * 10e-6)
	{
		alfa[3] = f3 / f;
		betta[3] = f2 / (f * alfa[2]);
	}
	else
	{
		alfa[3] = f3 / (f * alfa[1]);
		betta[3] = (f1 * betta[2] + f2 * alfa[2]) / f;
	}

	res._omega = OrbitElements::asin2(alfa[3], betta[3]);

	/*fi2 = _e * betta[3];
	lambda = sqrt(1.0 - _e * _e);
	rom = (_a * (1 - _e * _e)) / (1 + fi2);
	alfa[4] = -(alfa[3] * lambda) / (1.0 + fi2);
	betta[4] = (_a - rom) / (_a * _e);
	rab = OrbitElements::asin2(alfa[4], betta[4]);
	ElemOrb["tomega"] = ElemOrb["tau"] + 1.0 / nn * (rab - ElemOrb["e"] * alfa[4]);
	if (nu.Data.TimeOfDay.Seconds < ElemOrb["tomega"])
		ElemOrb["tomega"] = ElemOrb["tomega"] - ElemOrb["T"];*/
	return res;

}

double OrbitElements::asin2(double a, double b)
{
	if (a > 0.0)
	{
		if (b > 0.0)
			return asin(a);
		return PI - asin(a);
	}
	if (b > 0.0)
		return 2 * PI + asin(a);
	return 3 * PI / 2 + asin(a);
}

Vector OrbitElements:: KeplerToVector(const Date_Time& dateTime)
{
	double x, y, z, Vx, Vy, Vz, ksi, eta, r, ksi_pr, eta_pr, alpha, beta, gamma, alpha_pr, beta_pr, gamma_pr;	
	_M += _n * Date_Time::DifferenceDateTimeSec(_dateTime, dateTime);
	_E = _M;

	double dE = 0.0;	
	do
	{
		dE = (_E - _e * sin(_E) - _M) / (1 + _e * cos(_E));
		_E = _E - dE;
	} while (abs(dE) > (M_PI / 648000));
	
	ksi = _a * (cos(_E) - _e);
	eta = _a * sqrt(1 - _e * _e) * sin(_E);
	r = _a * (1 - _e * cos(_E));
	ksi_pr = -_a * _a * _n * sin(_E) / r;
	double tt  = sin(_E);
	eta_pr = _a * _a * _n * sqrt(1 - _e * _e) * cos(_E) / r;
	alpha = cos(_omega) * cos(_omega1) - sin(_omega) * sin(_omega1) * cos(_i);
	beta = cos(_omega1) * sin(_omega) + sin(_omega1) * cos(_omega) * cos(_i);
	gamma = sin(_omega1) * sin(_i);
	alpha_pr = -sin(_omega1) * cos(_omega) - cos(_omega1) *sin(_omega) * cos(_i);
	beta_pr = -sin(_omega1) * sin(_omega) + cos(_omega1) * cos(_omega) * cos(_i);
	gamma_pr = cos(_omega1) * sin(_i);
	x = alpha * ksi + alpha_pr * eta;
	y = beta * ksi + beta_pr * eta;
	z = gamma * ksi + gamma_pr * eta;
	Vx = alpha * ksi_pr + alpha_pr * eta_pr;
	Vy = beta * ksi_pr + beta_pr * eta_pr;
	Vz = gamma * ksi_pr + gamma_pr * eta_pr;
	Date_Time tmpDateTime = dateTime;
	
	Vector vect(_vitok, tmpDateTime.GetSecOfDay(), x, y, z, Vx, Vy, Vz, 0, dateTime);	
	return vect;
}

Vector OrbitElements::KeplerToVectorGSK(const Date_Time& dateTime)
{
	_vect = KeplerToVector(dateTime);
	double zvr,sut;
	int dat1[3],ct;
	Rdat(dat1);
	Dat84(n75,dat1,ct);
	sut=ct+d75;
	zvr=StarTime(sut,_t);
	double sn,cs;
	sn=sin(zvr);
	cs=cos(zvr);
	double x = _vect.x*cs + _vect.y*sn;
	_vect.y = -_vect.x*sn + _vect.y*cs;
	_vect.x = x;
	return _vect;
}