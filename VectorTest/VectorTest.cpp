#include "gtest/gtest.h"
#include "Date_Time.h"
#include "Date_Time.cpp"
#include "Atmosphera.h"
#include "Atmosphera.cpp"
#include "OrbitElements.h"
#include "OrbitElements.cpp"
#include "Vector.h"
#include "Vector.cpp"
#include "Spacecraft.h"
#include "Spacecraft.cpp"
#include "Soyuz.h"
#include "Soyuz.cpp"
#include "Engine.h"
#include "Engine.cpp"
#include "Printer.h"
#include "Printer.cpp"
#include "Matrix.hpp"
#include "GanimedLandModule.h"
#include "GanimedLandModule.cpp"
#include "ZoneRadioVision.h"
#include "ZoneRadioVision.cpp"
#include "GroundStation.h"
#include "GroundStation.cpp"
class TestVector : public testing::Test
{
protected:
	void SetUp()
	{
		Vector::f = 110;
		Vector::ap = 16;
		Date_Time DataTime(2011,10,05);
		vect = new Vector(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
		soyuz = new Soyuz(*vect, 6746.0, 1158.0, 2886.0, 115.20, 0.385080, -0.036830, 0.0, 0.008234, 273.0, 6, 1, 2);
		Engine skd(1,"SKD", 300, 302);
		soyuz->AddEngine(skd);
	}
	void TearDown() override
	{
		delete vect;
		delete soyuz;		
	}
	Vector* vect;
	Soyuz* soyuz;

};

/*
TEST_F(TestVector, gpz36x36_is_working)
{
	Vector::f = 110;
	Vector::ap = 16;
	Date_Time DataTime(2014,10,05);
	Vector vectorAsn(2, 47519.0, 4877.474915, 810.551593, 4485.757142, -3.800656006, 5.567955078, 3.124009521, 0.03, DataTime);
	vectorAsn.Prognoz((unsigned int)3);
	ASSERT_NEAR(vect->h, 383.39854, 0.0001);
}*/

/*
TEST_F(TestVector, test_readingGpzAndNorm_is_OK)
{
	Printer::ReadFileName = "..\\dpz_norm4-36.txt";
	double ** c = vect->CreateSquareMatrix(37);
	double ** d = vect->CreateSquareMatrix(37);
	double norm[37][37];
	Printer::ReadGpzFromFile(c,d);
	int delta = 0;
	int k = 0;
	double res [1300];
	long double koef1;
	long double koef2;
	double koef3;
	for (int i = 2; i < 37; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			j == 0 ? delta = 1 : delta = 2;
			koef1 = vect->factorial(i - j);
			koef2 = vect->factorial(i+j);
			norm[i][j] = sqrt((delta*(2 * i + 1)* (koef1))/koef2);
			res[k] = c[i][j];
			k ++ ;
		}
	}
	Printer::OutFileName = "..\\gpzN1.txt";
	Printer::Precision = 10;
	Printer::PrintDoubleToFile(res, k);
	vect->DeleteSquareMatrix(37, c);
	vect->DeleteSquareMatrix(37, d);
	SUCCEED();
}*/

/*
TEST_F(TestVector, test_readingGpzAndNorm_is_OK)
{
	Printer::ReadFileName = "..\\dpz4-36.txt";
	double c[665];
	double d[665];
	double c0[35];
	Printer::ReadGpzFromFile2(c0,c,d);
	auto i = 1;

}*/


/*
TEST_F(TestVector, test_polinom_is_working)
{
	bool kk = false;
	int n = 7;
	double* a = new double[n];
	double res;
	double mass[500];
	double mass1[500];
	int countStr;
	Printer::ReadFileName = "C:\\PROG_ASN\\asn.003";
	Printer::OutFileName = "..\\resultAsn.txt";
	Printer::Precision = 11;
	Printer::ReadDoublesFromFile(1,4,mass, mass1, countStr);
	double * resMass = new double [countStr];
	for (int i = 0; i < countStr; i++)
	{
		vect->Polinom(kk, n, countStr, mass, mass1, mass[i], a, res);
		resMass[i] =  fabs(mass1[i]) - fabs(res);
		//resMass[i] = res;
	}
	Printer::PrintDoubleToFile(mass, resMass, countStr);
	delete [] resMass;
	delete [] a;
	ASSERT_EQ(countStr, 408);
}
*/

TEST_F(TestVector, Kepler_To_Vector)
{
	double M0 = 0;
	double e = 0.00045;
	double a = 2831.2;	
	double omega = 0 * M_PI / 180.0;
	double i = M_PI / 2;
	double omega1 = 0 * M_PI / 180.0;
	double n = sqrt(mu_g/(a*a*a));
	Date_Time dateTime(2015, 02, 20, 00, 00, 00);
	Date_Time dateTimeStart(1997, 01, 16, 00, 00, 00);

	OrbitElements oe (1, dateTimeStart.GetSecOfDay(), a, e, omega, i, omega1, M0, n,dateTimeStart);
	Vector vect = oe.KeplerToVector(dateTime);	

	n = 1.016e-5;
	M0 = 192.417;
	e = 0.0013;
	a = 1070400;
	omega = 192.417 * PI / 180.0;
	i = 0.177 * M_PI / 180;
	omega1 = 63.552 * M_PI / 180;

	OrbitElements oe1 (1, dateTime.GetSecOfDay(), a, e, omega, i, omega1, M0, n, dateTime);
	Vector vect1 = oe1.KeplerToVector(dateTime);	
	int ii = 0;

	ASSERT_NEAR(vect.x, 2829.9211857915025, 0.0000001);
	ASSERT_NEAR(vect.y, 0.00000000000000031835908591289033,  0.000001);
	ASSERT_NEAR(vect.z, 5.1993701472998666,  0.0000001);
	ASSERT_NEAR(vect.Vx, -0.003433531336410129, 0.000001);
	ASSERT_NEAR(vect.Vy, 0.00000000000000011447918295357312, 0.0000001);
	ASSERT_NEAR(vect.Vz, 1.8696486850039278 , 0.000001);
}

TEST_F(TestVector, kepler_is_working)
{
	Date_Time DataTime(2015,02,20);
	Vector nu(1, DataTime.GetSecOfDay(), 2829.9211857915025, 0.00000000000000031835908591289033, 5.1993701472998666,
		-0.003433531336410129, 0.00000000000000011447918295357312, 1.8696486850039278, 0, DataTime);
	OrbitElements oe = OrbitElements::Kepler(nu);
	ASSERT_NEAR(oe.Hmax(), 201.27403999999933, 1e-4);
	ASSERT_NEAR(oe.Hmin(), 198.725960000001, 1e-4);
	ASSERT_NEAR(oe.I(), 1.5707963267948966, 1e-4);
	ASSERT_NEAR(oe.e(), 0.000449999999999862, 1e-4);
	ASSERT_NEAR(oe.A(), 2831.1999999999998, 1e-4);
	ASSERT_NEAR(oe.omega(), 6.2831853071795862, 1e-4);
	ASSERT_NEAR(oe.Omega(), 1.8636941088196963e-35, 1e-4);
	ASSERT_NEAR(oe.T(), 9518.8637263964065, 1e-4);
	ASSERT_NEAR(oe.Tau(), -2.7809311670630898, 1e-4);
}



TEST_F(TestVector, spacecraft_ganimedDeorb)
{
	GanimedLandModule landingModule;
	Engine mainEngine(1,"KTD", 600, 320);
	Engine secondEngine(2,"DMP", 120, 298);
	landingModule.AddEngine(mainEngine);
	landingModule.AddEngine(secondEngine);
	vector<Spacecraft> result = landingModule.BalDeorb();
	ASSERT_NEAR(result[0].h, 2, 1e-1);
	ASSERT_NEAR(result[1].h, 0, 0e-1);

}



TEST_F(TestVector, spacecraft_removeEngines)
{
	Date_Time DataTime(2011,10,05);
	Vector vvk(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	vvk.Perevod();
	vvk.Predict20(86380);
	vvk.Perevod();


	Spacecraft testsSpCr(*vect, 6746.0, 1158.0, 2886.0, 115.20, 0.385080, -0.036830, 0.0, 0.008234, 273.0, 6, 1, 2);
	Engine skd(1,"SKD", 300, 302);
	Engine skd2(2,"SKD2", 300, 302);
	testsSpCr.AddEngine(skd);
	testsSpCr.AddEngine(skd2);
	ASSERT_EQ(testsSpCr.CountOfEngines(), 2);
	testsSpCr.RemoveEngine(2);
	ASSERT_EQ(testsSpCr.CountOfEngines(), 1);
	Engine dpo(2,"DPO", 49, 302);
	testsSpCr.AddEngine(dpo);
	ASSERT_EQ(testsSpCr.CountOfEngines(), 2);
	testsSpCr.RemoveEngine(dpo);
	ASSERT_EQ(testsSpCr.CountOfEngines(), 1);
}



TEST_F(TestVector, soyuz_controlledDeorbitation)
{
	ASSERT_EQ(6, soyuz->TestRef());
	vector<Spacecraft> result = soyuz->ControlledDeorb(1819);
	ASSERT_GT(result.size(), (unsigned int)0);
	/*//active part ending
	ASSERT_NEAR(84290.593, result[1].t, 0.01);
	ASSERT_NEAR(3962.798,  result[1][1], 0.03);
	ASSERT_NEAR(-3257.793, result[1][2], 0.03);
	ASSERT_NEAR(-4392.766, result[1][3], 0.03);
	//division	
	ASSERT_NEAR(85699.0773, result[2].t, 0.06);
	ASSERT_NEAR(5291.3042,  result[2][1], 0.3);
	ASSERT_NEAR(2203.8246,  result[2][2], 0.3);
	ASSERT_NEAR(3093.30753, result[2][3], 0.3);
	//start atmosphere part
	ASSERT_NEAR(85904.6220, result[3].t, 0.14);
	ASSERT_NEAR(4207.9380,  result[3][1], 0.9);
	ASSERT_NEAR(2892.3440,  result[3][2], 0.9);
	ASSERT_NEAR(3975.1020,  result[3][3], 0.9);
	ASSERT_NEAR(-1.45640,   result[3].Tet()*PerRadGrad, 0.0014);
	//start guidance part 
	ASSERT_NEAR(86019.063,  result[4].t, 0.3);
	ASSERT_NEAR(3497.10424, result[4][1], 2.0);
	ASSERT_NEAR(3216.07155, result[4][2], 2.0);
	ASSERT_NEAR(4362.23860, result[4][3], 2.0);
	ASSERT_NEAR(81.6013000, result[4].h,  0.07);
	//rotation
	ASSERT_NEAR(86280.339, result[5].t, 0.25);
	ASSERT_NEAR(1864.2230, result[5][1], 1.6);
	ASSERT_NEAR(3687.5830, result[5][2], 0.9);
	ASSERT_NEAR(4902.4850, result[5][3], 0.9);
	ASSERT_NEAR(45.911000, result[5].h,  0.07);
	//parachute open
	ASSERT_NEAR(86463.770, result[6].t, 0.25);
	ASSERT_NEAR(1541.3570, result[6][1], 1.6);
	ASSERT_NEAR(3742.4620, result[6][2], 0.9);
	ASSERT_NEAR(4926.7130, result[6][3], 0.9);
	ASSERT_NEAR(10.700000, result[6].h,  0.001);*/

	ASSERT_NEAR(50.67, result[6][8], 0.07);
	ASSERT_NEAR(67.33, result[6][9], 0.2);
	ASSERT_NEAR(10.700000, result[6].h,  0.001);

	SUCCEED();
}

TEST_F(TestVector, soyuz_ballisticDeorbitation)
{
	Date_Time DataTime(2011,10,05);
	Vector *vect = new Vector(1816, 70500, -4887.4796260, -1988.8661830, 4242.4879740, -0.536599240, -6.3971388400, -3.612259210, 0.04158185, DataTime);
	Spacecraft testsSpCr(*vect, 6746.0, 1158.0, 2886.0, 115.20, 0.385080, -0.036830, 0.0, 0.008234, 273.0, 6, 1, 2);
	Engine skd(1,"SKD", 300, 302);
	testsSpCr.AddEngine(skd);
	vector<Spacecraft> result = testsSpCr.BalDeorb(1819);
	ASSERT_NEAR(86356, result.back().t, 1.0);
	ASSERT_NEAR(63, result.back()[9], 0.0166667);
    ASSERT_NEAR(1881.38985, result.back()[1], 0.00001);
}

TEST_F(TestVector, z_eq_zero)
{
	vect->Prognoz(static_cast<unsigned int>(1817));
	double z = vect->operator[](3);
	double x = vect->operator[](1);
	double y = vect->operator[](2);
	ASSERT_NEAR(0, z, 0.0001);
	ASSERT_NEAR(4651.9107, x, 0.001);
	ASSERT_NEAR(4906.6039, y, 0.001);
}

TEST_F(TestVector, soyuz_aerodynamic_case_1)
{
	double s = 0;
	double s1 = 0;
	double b = 0;
	soyuz->Perevod();
	soyuz->Aerodynamic(5.0, s, s1, b); // case 3 max > 6 && h > 130
	ASSERT_NEAR(s,8.38046713413349, 0.000001);
	ASSERT_NEAR(s1,8.63614637837593, 0.000001);
	ASSERT_NEAR(b,2.08585591123745, 0.000001);
}

TEST_F(TestVector, soyuz_aerodynamic_case_2)
{
	double s = 0;
	double s1 = 0;
	double b = 0;
	soyuz->Perevod();
	soyuz->h = 65.0;
	soyuz->Aerodynamic(8.5, s, s1, b); // case 2 max > 6 &&  30 <= h < 90
	ASSERT_NEAR(s,8.51175843173862, 0.000001);
	ASSERT_NEAR(s1,8.75920062379871, 0.000001);
	ASSERT_NEAR(b,2.06726001453181, 0.000001);
}

TEST_F(TestVector, soyuz_aerodynamic_case_3)
{
	double s = 0;
	double s1 = 0;
	double b = 0;
	soyuz->Perevod();
	soyuz->Aerodynamic(30.0, s, s1, b); // case 3 max > 6 && h > 130
	ASSERT_NEAR(s,14.2968545546068, 0.000001);
	ASSERT_NEAR(s1,14.3173726619517, 0.000001);
	ASSERT_NEAR(b,0.766230895763406, 0.000001);
}


int maint (int argc, char *argv[])
{
	testing::InitGoogleTest(&argc, argv);
	
	return RUN_ALL_TESTS();
}