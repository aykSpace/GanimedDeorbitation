#pragma once

#include "GroundStation.h"
#include "Vector.h"
#include <algorithm>
#include <iterator>

class ZoneRadioVision
{
public:
	ZoneRadioVision();
	ZoneRadioVision(ZoneRadioVision const& _copy);
	~ZoneRadioVision();

	ZoneRadioVision & operator = (ZoneRadioVision & _other);

	vector<ZoneRadioVision> Calculate(const GroundStation & _grSt, const Vector &_vect, unsigned _numOfCirc, double _gamm1 = 7.0); // ��������� ��� ��� ��������� ���-�� ������ _gamma1 - ���� ������
	vector<Date_Time> GetDatTime();  //getter-�, ����� ������������� � ��-��
	vector<double>    GetDistance();
	vector<double>    GetAzimuts();
	vector<double>    GetCornersOfLocation();
	vector<double>    GetFis();
	vector<double>    GetLambdas();
	vector<double>    GetAltitudes();
	vector<double>    GetCornersOfSun();
	unsigned          GetVitok();
	unsigned          GetNumSutCirc();
	int               GetSizeOfMassives();
	bool              GetVisible();
private:
	ZoneRadioVision(Date_Time & _datTimes, double & _distances, double & _azimuths, double & _cornersOfLocation,
		            double & _fis, double & _lamdbas, double & _altitudes, double & _cornersOfSun); //��������� ����������� ��� Calculate

	 Date_Time *datTimes; //can use <vector> or <list> instead of this massives
	 double *distances;
	 double *azimuths;
	 double *cornersOfLocation;
	 double *fis;
	 double *lambdas;
	 double *altitudes;
	 double *cornersOfSun;	 
	 int sizeOfmassives;
	 unsigned vitok;
	 unsigned nSutCircle;
	 bool visible;
	 void VectorOfGrSt(double** _matrPer, double _fiGrSt, double _lamGrSt, double _hGrSt, double *_outVectGrSt); //������ ������� �������� ������� � ������� �������� � ��� ������� ��������� (�������)
	 void DisAziCor(double** _matrPer, const Vector &_vect, double *_vecGrSt, double &_distance, double &_azimuth, double &_cornerLoc);//������ ��������� ������� � ���� ����� (�/� ��� � ��)
	 double CalcVectorGrStMatrGamma( double ** _matrPer, double fiRad, double lamRad, const GroundStation &_grSt, double *vectGrSt);//������ ���� �� �������� ���������� ��� 
	 void InitializationMassives(unsigned _sizeOfMassives); // �������� ������ ��� ������� � �����������
	 ZoneRadioVision Integration (double **_matrPer, double **_am, double **_matrTimeDisAzCor, double *_vectGrSt, double _gamma0, Vector &_vect, int _dt, double _gamma1, unsigned &_counter); 
	 void fromCourrentToResult (ZoneRadioVision & current, ZoneRadioVision & result, unsigned counter);

	 //�������������� ���, (����������� �� ����� �����)

	 void Swap (ZoneRadioVision & _other);
	// ����� other and this
	 void GetResultZrv (ZoneRadioVision & res, double ** am, unsigned orderAm,  double** matrInerpol, Vector & vectRes, unsigned counter, double gammaInterpol, double cornerOfSun); // ������ ��������������� ���
};

