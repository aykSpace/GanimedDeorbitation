#include "GroundStation.h"


//GroundStation::GroundStation()
//{}

GroundStation::GroundStation(unsigned int _numGrSt, double _fi, double _lam, double _hUnderSea, double _hUnderEarth, unsigned int _militaryNum, string _shortName, string _code)
	:NumGrStation(_numGrSt), Fi(_fi), Lambda(_lam), HUnderSea(_hUnderSea), HUnderEarth(_hUnderEarth), MilitaryNum(_militaryNum), ShortName(_shortName), Code(_code)
{}

GroundStation::~GroundStation()
{}

//------------------------------------статические члены--------------------------------

GroundStation GroundStation::Shelkovo()
{
	GroundStation _gr( 14, 55.9518, 37.9655, 0.1746, 0.1746, 26178, "МСК", "Лотограф");
	return _gr;
}

GroundStation GroundStation::StPetersburg()
{
	GroundStation _gr(9, 59.7138, 30.1953, 0.0948, 0.0948, 14108, "СПБ", "Глазница");
	return _gr;
}

GroundStation GroundStation::One()
{
	GroundStation _gr(1, 45.9099, 63.3324, 0.0789, 0.0789, 13951, "---", "Проточка");
	return _gr;
}

GroundStation GroundStation::Djusaly()
{
	GroundStation _gr(5, 45.7049, 63.3408, 0.135, 0.135, 74828, "ДЖС", "Наутилус");
	return _gr;
}

GroundStation GroundStation::Twelve()
{
	GroundStation _gr(12, 58.333, 82.8864, 0.0824, 0.0824, 14174, "КЛП", "Зимовщик");
	return _gr;
}

GroundStation GroundStation::UlanUde()
{
	GroundStation _gr(13, 51.8707, 107.9464, 0.6193, 0.6193, 14129, "УЛД", "Диванчик");
	return _gr;
}

GroundStation GroundStation::Barnaul()
{
	GroundStation _gr(17, 53.32, 83.36, 0.257, 0.257, 96634, "БРН", "Лесовик");
	return _gr;
}

GroundStation GroundStation::Fifteen()
{
	GroundStation _gr(15, 44.0246, 131.7565, 0.0923, 0.0923, 14038, "ЖК", "Жница");
	return _gr;
}

GroundStation GroundStation::Petropavlovsk()
{
	GroundStation _gr(6, 53.1038, 158.3609, 0.0496, 0.0496, 14086, "ППК", "Выправка");
	return _gr;
}