#pragma once
#include "Vector.h"
#include "BALCONST.H"

class OrbitElements
{
private:
	unsigned long _vitok;
	double	_date, _t, _a, _e, _omega, _i, _omega1, _u, _M, _E, _S, _h, _hMin, _hMax, _tau, _T; //omega1 - Omega
	double	CalcN()			{return sqrt(mu_e/(_a*_a*_a));}									// Угловая скорость движения по орбите (рад./с)
	double	CalcT()			{return M_PI*2.*sqrt(_a*_a*_a/mu_e);}
	double	CalcGamma2()	{return -C20*(a_e*a_e)/(_a*_a);}
	double	CalcOmega_p()	{return -1.5*sqrt(mu_e)*CalcGamma2()*cos(_i)/sqrt(_a*_a*_a);}		// Угловая скорость прецессии (рад./с)
	double	CalcH()			{return _a-a_e-0.8;}
	Vector	_vect;
	Date_Time _dateTime;
	void	Rdat(int mdat[3]);														// Распаковка даты в формате ЧЧММГГГГ в дату-массив: число, месяц, год
	int		Dat84(const int id[3], const int id1[3], int &cyt);						// Функция счета числа суток между двумя исходными календарными дата и (старый вариант)
	int		Dat84(const int id[3], const int id1[3]);								// Функция счета числа суток между двумя исходными календарными датами
	double	Rtime(double tzip);														// Преобразование времени в формате ЧЧММСС.ДДД в число секунд
	double	RtimeZip(double t);										    			// Преобразование времени в секундах в формат ЧЧММСС.ДДД
	double	RdatZip(int mdat[3]);
	void	Ucopy(double a[], double b[], int n);									// Пересылка элементов из одного массива в другой
	int		Iround(double x);														// Округление
	double  Angle(double x, double y);

	double	StarTime(double sut, double t0);										// Звездное время

	void	Dj(int &jt, double &ft, int iy, int im, int id, double d);
	void	Jd(int jt, double ft, int &iy, int &im, int &id, double &d);			// Новая, непроверенная процедура !!! Перевод с Фортрана (Glyba)

	void	Nua( double nu2[6], double s, char pp);									// Перевод вектора состояния из относительной системы координат в абсолютную и обратно s - звездное время
	void	Nua(double nu1[11],double nu2[11],char pp);								// Перевод вектора состояния из относительной системы координат в абсолютную и обратно
	void	Rvosk(double nu2[6]);													// Расчет оскулирующих элементов орбиты по заданному вектору состояния в АСК
	void	OskRv(OrbitElements &osk, double nu[6]);								// Переход от оскулирующих элементов орбиты к вектору состояния в относительной ГСК
	void	OskNu(OrbitElements &os, double ny[11]);								// Перевод оскулирующих элементов в вектор в ГСК
	static double asin2(double a, double b);

	double _n;

public:
	 // конструкторы класса
	OrbitElements();	
	OrbitElements(unsigned long vitok, double date, double t, double a, double e, double omega, double i, double omega1, double u, double M, double E, double S);
	OrbitElements(unsigned long vitok, double t, double a, double e, double omega, double i, double omega1, double M, double n, Date_Time dateTime);
	OrbitElements(Vector vect);
	

	~OrbitElements();

	
	//перегрузка оператора
	double & operator[] (int);
	
	double GetU(){return _u;}
	void SetVector(const Vector &  vect){_vect = vect;}
	
	void GetElements();
	static OrbitElements Kepler(const Vector& vectorNu); // перевод в элементы орбиты

	//свойства
	double Hmax(){return _hMax;}
	double Hmin(){return _hMin;}
	double A(){return _a;}
	double I(){return _i;}
	double T(){return _T;}
	double Omega(){return _omega1;}
	double omega(){return _omega;}
	double e(){return _e;}
	double Tau(){return _tau;}

	Vector KeplerToVector(const Date_Time& dateTime);//расчет вектора на определенную дату из элементов орбиты
	Vector KeplerToVectorGSK(const Date_Time& dateTime);//расчет вектора на определенную дату из элементов орбиты
	
};

