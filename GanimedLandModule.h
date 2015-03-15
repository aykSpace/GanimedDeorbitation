#pragma once
#include "spacecraft.h"
class GanimedLandModule :
	public Spacecraft
{
public:
	GanimedLandModule();
	~GanimedLandModule();

	void RKSGan (Spacecraft & vect, double dt); //ф-я интегрирования методом Рунге - Кутта
	void RKSGan (Spacecraft & vect, double dt, Spacecraft (*f_ptr)(const Spacecraft & vect, double dt));

	friend Spacecraft FGanimed (const Spacecraft & vect, double dt); // ф-я рассчета правых частей
	friend Spacecraft FGanimedGSK (const Spacecraft & vect, double dt); // ф-я рассчета правых частей в ГСК


	vector<Spacecraft> BalDeorb();
	void EngineZoneGanimed (Spacecraft & SC, vector<Spacecraft> &Rezult, double _Tvkl = 0, unsigned numberOfCurrentEngine = 0); //ф-я рассчета активного участка

private:
	void SetCurrentVector(const Vector &  vect);
};

