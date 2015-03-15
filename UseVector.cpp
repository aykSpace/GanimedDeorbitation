// отладка класса Vector

#include "Vector.h"
#include "Spacecraft.h"
#include <time.h>
#include "GroundStation.h"
#include "ZoneRadioVision.h"
#include "OrbitElements.h"


int testmass (double &_i, int n);
void testmasss (double ** matr);

int main(){
	
	Date_Time DataTime(2011,10,05);
	Date_Time DataTime1(2011,10,07,23,59,59);
	bool re = DataTime1 > DataTime;
	srand(time(NULL));
	int a = rand()%1000;
	int b = rand()%2000;
	Vector::f = 110;
	Vector::ap = 16;
	Vector test(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Date_Time DataTimetst(2014,24,1);
	Vector tst1(2898, 24840, 424.78293230, -6774.344163, 299.88650640, 4.23755698, 0.54054299, 6.00138935, 0.0459,  DataTimetst);
	//Date_Time testDt(2014, 1, 24, 6, 54, 0, 0);

	//Vector testError(1816, 70500, 1234.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	try
	{
		tst1.Prognoz((unsigned int)3000);
		//test.Prognoz((unsigned int)2060);
	}
	catch (int exception)
	{
		cout << exception; 
	}



	//Vector testDatTimeVector(2898, 24840, 424.78293230, -6774.344163, 299.8865064, 4.23755698, 0.54054299, 6.00138935, 0.0459, testDt);

	/*Date_Time testDtPredict(2014, 2, 1, 6, 54, 0, 0);

	if (testDtPredict > testDt)
	{
		testDatTimeVector.Prognoz(testDtPredict);
	}*/

	//test.Prognoz((unsigned int)1817);
	/*OrbitElements orbEl(test);
	orbEl.GetElements();
	try
	{
		Vector vectForU = test.Prognoz(-1.0);
		orbEl.SetVector(vectForU);
		orbEl.GetElements();
	}
	catch(int exception) //catch(char* exception)
	{
		cout << exception; 
	}*/
	


	/*Vector test1(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Vector test2(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Vector test3(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Vector test4(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Vector test5(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Vector test6(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Vector test7(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	span s(series, _T("Phase1"));
	test1.Prognoz((unsigned int)1819);
	s.~span();
	span s1(series, _T("Phase2"));
	test.Prognoz((unsigned int)2016);
	s1.~span();
	span s2(series, _T("Phase3"));
	test2.Prognoz((unsigned int)2000);
	s2.~span();
	span s3(series, _T("Phase4"));
	test3.Prognoz((unsigned int)1850);
	s3.~span();
	span s4(series, _T("Phase5"));
	test4.Prognoz((unsigned int)1850);
	s4.~span();
	span s5(series, _T("Phase6"));
	test5.Prognoz((unsigned int)1850);
	s5.~span();
	span s6(series, _T("Phase7"));
	test6.Prognoz((unsigned int)1850);
	s6.~span();
	span s7(series, _T("Phase7"));
	test7.Prognoz((unsigned int)1850);
	s7.~span();*/
	//vector<Vector> testing = test.Prognoz(2016);
	//Vector test11 = test;
	/*Vector testdt = test.Prognoz(DataTime1);
	double * mass = new double[3]; // * - указатель
	mass[0] = 0.2344;
	mass[1] = 2./3.;
	mass[2] = 3./3.;
	testmass(*mass,3); // *-оператор преобразования от указателя к значению
	delete [] mass;*/

	//--------------тест зрв---------------------------------------
	/*Date_Time testDT(2013,7,18);
	Vector test1(3950, 21600, 4651.198865, 1791.176601, 4613.552035, -4.795183652, 4.706476698, 2.996089725, 0.0421, testDT);*/
	/*ZoneRadioVision zrv;
	zrv.Calculate(GroundStation::Shelkovo(),test,16,7);*/



	/*Vector test1(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	//vector<Vector> testing = test.Prognoz(2016);
	vector<Vector> testing1;
	unsigned i = 1816;
	while (i <= 2016)
	{
		test1.Prognoz();
		i +=1;
	}* /
	//int i = testing.size();*/
	//Spacecraft Soyuz(test, 6746.0, 1158.0, 2886.0, 300.0, 302.0, 115.20, 0, 0.385080, -0.036830, 0.0, 0.008234, 49.40, 273.00, 200, 6, 1, 2);
	//Soyuz.Prognoz(2016, DataTime);
	//Soyuz.BalDeorb(1819);
	cout << "Gotovo" << endl;
	getchar();
	return 0;
}


int testmass (double &_i, int n) //&-ссылка (можно принять за значение)
{
	double *mass = &_i; //*-оператор преобразования от указателя к значению
	for (int i = 0; i < n; i++)
	{
		cout << mass[i] << endl;
	}
	return 1;
}

void testmasss (double ** _matrPer)
{
	_matrPer[0][0] = 0.1;
	_matrPer[0][1] = 0.2;
	_matrPer[0][3] = 0.3;
	_matrPer[1][0] = 0.4;
	_matrPer[1][1] = 0.5;
	_matrPer[1][2] = 0.7535703753;
	_matrPer[2][0] = 0.7847847;
	_matrPer[2][1] = 0.3829474;
	_matrPer[2][2] = 0.0;
}