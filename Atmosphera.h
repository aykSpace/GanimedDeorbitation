#pragma once
#include "Date_Time.h"

#define OMEGA_EARTH 7.29211587e-5

class Atmosphera
{
public:
	double A[3], B[3], C[4], N[2], Fi, D[3], E[7], L[3]; //������������ ������ ��� ��������� F0
	double S0, SDSol, CDSol, ASol, d1900;  //�������� �����, sin delta, cos delta, alfa, ���-�� ���� �� ������ ���� (1900)
	int d; //���-�� ���� �� ������ ����
	double DRo;   //������� �������� ���������
	double W[3]; // ������� �������� �����
	unsigned NMounth; //����� ������
	double VZV; // �������� �����
	
	//������������ � �����������
	Atmosphera();
	~Atmosphera();

	//������ ������
	void   SelectKoef(int aF0, double h);                     //��������� ������ ������������� �� �������
	void   equal(unsigned n, double *mass1, double *mass2);  //��������� ������������ ��-� �������2 ������� 1, n - ����������� ��������
	void   StarTime (Date_Time & Data);                     //��������� �������� ��������� ������� �� ����������� �������
	void   Sun (double t);                                 //��������� �������� ��������� ������ �� ������� ������ ������� (������� ���� ������� StarTime!)
	double Ad_D  (double days_year) const;                //��������� ������ ������������ �(d) ����������� ������ 
	double Kp_Ap (double App) const;                     //��������� ��������� �� Ap � ��
	int    Selection_F0 (int F);                        //��������� ������ ���������� �-� F0 �� �������    
	double LowAtm (double _h);						   //��������� ���������� ��������� ������ ���������
	void   VarAtm (double H, double Fi, double Lam);  //��������� ����� �������� ��������� ��������� � �����
	void   LinInterpolMatr(unsigned Col, const double Matr[][10], unsigned NumArg, double DArg, double *RezArr); //��������� �������� ������������ ��� ������ ����� � ���������
	                                                                                                       //(NumArg-����� ��������� ����� poi, DArg-������� ����� ���� ��-�� �������� Rez - �������������� ������)
};

inline double lin2 (double x1, double x2, double y, double y1, double y2)// ��������� �������� ������������
{ 
	return ((y - y2)*(x2 - x1)/(y2 - y1)) + x2;
}
double lin2Matr(unsigned nRow, unsigned nCol, double ** matr, unsigned row1, unsigned col1, unsigned row2, unsigned col2, double argRow, double argCol); //�������� ������������ � ���������� �������
																															//(argCol, argRow - �������� ���������� �� ������ � �������)

double lin3Matr(unsigned n, unsigned m, unsigned l, double *** matr, unsigned i, unsigned j, unsigned k, unsigned i1, 
				unsigned j1, unsigned k1, double di, double dj, double dk); //�������� ������������ � ���������� ������� //(argCol, argRow - �������� ���������� �� ������ � �������)

int poi (int n, const double *mass, double par); //��������� ������ � �������, ������������� �� ����������� //n - ����� ��������� ������� -1, *mass - ������(�������������), par - ��������
int poi (int n, const int *mass, const double par); //��������� ������ � �������, ������������� �� ����������� //n - ����� ��������� ������� -1, *mass - ������(�������������), par - ��������
