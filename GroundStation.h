#pragma once
#include <xstring>
using namespace std;
class GroundStation
{
	public:
		GroundStation(unsigned int _numGrSt, double _fi, double _lam, double _hUnderSea, double _hUnderEarth, unsigned int _militaryNum = 0, string _shortName = "-" ,string _code = "-");
		static GroundStation Shelkovo();
		static GroundStation StPetersburg();
		static GroundStation One();
		static GroundStation Djusaly();
		static GroundStation Twelve();
		static GroundStation UlanUde();
		static GroundStation Barnaul();
		static GroundStation Fifteen();
		static GroundStation Petropavlovsk();
		~GroundStation();

		unsigned int NumGrStation; 
		//GroundStation();
		double Fi, Lambda;
		double HUnderSea;  //������ ������� ��� ������� ���� (km)
		double HUnderEarth; // ������ ����������� ��� ����������� (��� ����� �������� � �������������� ����) (km)
		string ShortName; //���������� �� 3-� ����
		string Code;
		unsigned int MilitaryNum;

};

//����� ���������� ����������� ����� ���������� � ���������� GroundStation ��� � ���� ��� � ��!!!!!!