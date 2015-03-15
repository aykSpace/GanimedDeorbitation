#include "Date_Time.h"

Date_Time::Date_Time(unsigned _Year, unsigned _Mounth, unsigned _Day, unsigned _Hour, unsigned _Min, unsigned _Sec, double _MSec):
	Year(_Year), Mounth(_Mounth), Day(_Day), Hour(_Hour), Min(_Min), Sec(_Sec), MSec(_MSec)
{}

Date_Time::Date_Time(Date_Time const& copy): Year(copy.Year), Mounth(copy.Mounth), Day(copy.Day), Hour(copy.Hour), Min(copy.Min),
	Sec(copy.Sec), MSec(copy.MSec)
{}

Date_Time::~Date_Time()
{}

void Date_Time::SetFrmtTime(double _t)
{
	if (_t < COUNT_SEC_OF_DAY)
	{
		double A, Remainder, B, Remainder1, DC;
		A          = _t/3600.0;
		Hour       = static_cast<unsigned>(A);
		Remainder  = _t - 3600.0*Hour;
		B          = Remainder/60.0;
		Min        = static_cast<unsigned>(B);
		Remainder1 = Remainder - 60.0 * Min;
		Sec        = static_cast<unsigned>(Remainder1);
		DC         = (_t - static_cast<unsigned>(_t)) * 1000.0;
		MSec       = static_cast<unsigned>(DC);
	}
}

void Date_Time::Print() const
{
	cout << setw(2);
	cout << Day <<'.'<< Mounth << '.' << Year << ' ';
	cout.width(2);cout << Hour;
	cout <<':' ; 
	cout.width(2); cout<<Min;
	cout << ':' ;
	cout.width(2); cout << Sec <<'.' ;
	cout.width(3); cout << MSec;
}

void Date_Time::PrintDate() const
{
	cout << setw(2);
	cout << Day <<'.'<< Mounth << '.' << Year << ' ';
}

int Date_Time::Calendar(bool REZIM_FLAG)
{
	unsigned kdney[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	int CountDayYear;
	if (fmod(Year, 4.0) == 0)
		kdney[1] = 29;
	if (!REZIM_FLAG) // если флаг по умолчанию (false) то провер€ем изменение даты
	{
		if (Day != kdney[Mounth - 1])
			Day += 1;
		else 
		{
			Day = 1;
			if (Mounth == 12)
			{
				Mounth = 1;
				Year  += 1;
			}
			else
				Mounth +=1;
		}
	}
	else // если флаг True  возвращаем число дней от начала года
	{
		if (Mounth == 1)
			return Day;
		else 
		{
			CountDayYear = Day;
			unsigned j = Mounth - 1;
			for (unsigned i = 0 ; i < j; ++i)
			{
				CountDayYear += kdney[i];
			}
			return CountDayYear;
		}
	}
	return 0;
}

double Days (Date_Time &Data, int & DYear)
{
	int Year1,R;
	double D = 27393.5;
	Year1    = Data.Year-1;
	for (int i=1975; i <= Year1; i++)
	{
		fmod(i,4.0)==0 ? R=366 : R=365;
		D+=R;
	}
	DYear = Data.Calendar(true);
	D    += DYear - 1;
	return D;

}


bool Date_Time::operator>(Date_Time & _datTime)
{
	bool rez = false;
	if (Year > _datTime.Year)
	{
		rez = true;
		return rez;
	}
	if (Mounth > _datTime.Mounth)
	{
		rez = true;
		return rez;
	}
	if (Day > _datTime.Day)
	{
		rez = true;
		return rez;
	}
	if (Year == _datTime.Year && Mounth == _datTime.Mounth && Day == _datTime.Day)
	{
		double cSecPerDay = this->GetSecOfDay();
		if (cSecPerDay > _datTime.GetSecOfDay())
		{
			rez = true;
			return rez;
		}
	}
	return rez;
}

bool Date_Time::operator<(Date_Time & _datTime)
{
	if (this->operator>(_datTime))
		return false;
	return true;
}

double Date_Time::GetSecOfDay()
{
	return Hour*3600.0 + Min*60.0 + Sec + MSec/1000.0;
}

Date_Time & Date_Time::operator=(const Date_Time &_datTime)
{
	Year =_datTime.Year;
	Mounth = _datTime.Mounth;
	Day = _datTime.Day;
	Hour = _datTime.Hour;
	Min = _datTime.Min;
	Sec = _datTime.Sec;
	MSec= _datTime.MSec;
	//Swap(_datTime);
	return *this;
}

void Date_Time::Swap(const Date_Time& dateTime)
{
	std::swap(Year, const_cast<unsigned &>(dateTime.Year));
	std::swap(Mounth, const_cast<unsigned &>(dateTime.Mounth));
	std::swap(Day, const_cast<unsigned &>(dateTime.Day));
	std::swap(Hour, const_cast<unsigned &>(dateTime.Hour));
	std::swap(Min, const_cast<unsigned &>(dateTime.Min));
	std::swap(Sec, const_cast<unsigned &>(dateTime.Sec));
	std::swap(MSec, const_cast<double &>(dateTime.MSec));
}

int Date_Time :: DateFormat()
{
	return Day*1000000+Mounth*10000+Year;
}


void Date_Time::CheckDaySeconds(double &t)
{
	t > 86400 ? t -= 86400 : t;
	t < 0 ? t += 86400 : t;
}

double Date_Time::DifferenceDateTimeSec(const Date_Time & d1, const Date_Time  & d2)
{
	Date_Time tmpd1(d1);
	Date_Time tmpd2(d2);
	int dyear;
	double dat1 = Days(tmpd1,dyear);
	double dat2 = Days(tmpd2,dyear);
	return (dat2-dat1)*COUNT_SEC_OF_DAY + (tmpd2.GetSecOfDay() - tmpd1.GetSecOfDay());
}