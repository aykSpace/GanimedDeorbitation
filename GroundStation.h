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
		double HUnderSea;  //высота станции над уровнем моря (km)
		double HUnderEarth; // высота наблюдателя над эллипсойдом (ЗРВ может начаться с отрицательного угла) (km)
		string ShortName; //сокращение из 3-х букв
		string Code;
		unsigned int MilitaryNum;

};

//потом необходимо реализовать метод добавления и считывания GroundStation или в файл или в БД!!!!!!