#pragma once
#include "spacecraft.h"
class GanimedLandModule :
	public Spacecraft
{
public:
	GanimedLandModule();
	~GanimedLandModule();

	void RKSGan (Spacecraft & vect, double dt); //�-� �������������� ������� ����� - �����
	void RKSGan (Spacecraft & vect, double dt, Spacecraft (*f_ptr)(const Spacecraft & vect, double dt));

	friend Spacecraft FGanimed (const Spacecraft & vect, double dt); // �-� �������� ������ ������
	friend Spacecraft FGanimedGSK (const Spacecraft & vect, double dt); // �-� �������� ������ ������ � ���


	vector<Spacecraft> BalDeorb();
	void EngineZoneGanimed (Spacecraft & SC, vector<Spacecraft> &Rezult, double _Tvkl = 0, unsigned numberOfCurrentEngine = 0); //�-� �������� ��������� �������

private:
	void SetCurrentVector(const Vector &  vect);
};

