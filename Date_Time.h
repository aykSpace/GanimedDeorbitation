//����� ��� ������������� ���� � ������� 
#pragma once
#include <iomanip>
#include <iostream>

#define  COUNT_SEC_OF_DAY 86400
using namespace std;

class Date_Time
{
	//������ ������� hh:mm:sec.999,9999
public:
	unsigned Year, Mounth, Day, Hour, Min, Sec;
	double MSec;

	static void CheckDaySeconds(double &t);

	//------------------constructors-------------------------------------
	Date_Time()
		:Year(0), Mounth(0), Day(0), Hour(0), Min(0), Sec(0), MSec(0)
	{};
	Date_Time(unsigned _Year, unsigned _Mounth, unsigned _Day, unsigned _Hour = 0, unsigned _Min = 0, unsigned _Sec = 0, double _MSec = 0);
	Date_Time (Date_Time const& copy);
	~Date_Time();

	//------------------operators-----------------------------------------
	bool operator > (Date_Time & _datTime);
	bool operator < (Date_Time & _datTime); //�� ��������!!
	Date_Time & operator = (const Date_Time &);

	//------------------methods--------------------------------------------

	const unsigned Get_Mounth ()
	{
		return Mounth;
	}
	double GetSecOfDay();
	void Print () const;
	void PrintDate () const;
	void SetFrmtTime (double _t);
	int Calendar(bool REZIM_FLAG = false); // ���� REZIM_FLAG = 1,  �� ������� ���-�� ���� �� ������ ����
	int DateFormat();
	static double DifferenceDateTimeSec(Date_Time const & d1, Date_Time const & d2);// ������� ����� ������ (d2 - d1) � ��������
	friend double Days (Date_Time &, int & DYear); //��������� �������� ���-�� ���� �� 1900 ����, DYear - ���-�� ���� �� ������ ����
private:
	void Swap (const Date_Time& dateTime);
}; 

