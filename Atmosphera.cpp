#include "Atmosphera.h"

//конструкторы и деструкторы
Atmosphera::Atmosphera(): DRo(0) 
{
	W[0] = 0;
	W[1] = 0;
	W[2] = 0;
}
Atmosphera::~Atmosphera()
{}

//методы
void Atmosphera::equal ( unsigned n, double *mass1, double *mass2)
{
	for (unsigned i = 0; i < n; i++)  
		*(mass1+i) = *(mass2+i); 
}

void Atmosphera::StarTime (Date_Time & Data)
{
	double a_s[6] = {1.7399358945, 0.0172027912737, 0.675587865e-5, 0.409319754, -0.835464852e-4, 0.617119333e-5};
	double b_s[3] = {4.52360151,-0.0009242202, -0.00003626794};
	double c_s[3] = {0.196365056, 0.230895722, -0.00005604252};
	double k[3] = {6.12152393, 0.212768711, -0.00002504547};
	double R, F, D, T1, T2, SZ1;
	d1900 = Days (Data, d);
	T1  = d1900/36525;
	T2  = T1*T1;
	R   = b_s[0] + b_s[1] * d1900 + b_s[2] * T2;
	F   = c_s[0] + c_s[1] * d1900 + c_s[2] * T2;
	D   =   k[0] +   k[1] * d1900 +   k[2] * T2;
	S0  = a_s[0] + a_s[1] * d1900 + a_s[2] * T2 + cos(a_s[3]) * (a_s[4] * sin(R) - a_s[5] * sin(2*R + 2*F - 2*D));
	SZ1 = S0/6.283185308; //убираем число полных оборотов   с 1900 года
	R   = int(SZ1);
	S0  = S0 - R * 6.283185308;     //конечный результат зв времени меняется вместе с датой!
}

void Atmosphera::Sun (double t)
{
	double TC, TC2, TC3, t1;
	double aL0, aLamda, H, R, E0, DLH, DLH2, aLamda0, dPsi;
	double SL, CL, CE, SE, S1, AL, DSol;
	// рассчет координат Солнца
	t1=t;
	DSol = d1900+(t1-10800)/86400;  // рассчет координат Солнца для текущей даты!!!
	TC = DSol/36525;
	TC2 = TC*TC;
	TC3 = TC2*TC;
	aL0 = 0.01675104-0.0000418*TC-0.000000126*TC2;
	aLamda = 4.881627933+628.3319507*TC+5.279620987e-6*TC2;
	H = 4.908229468+3.000526417e-2*TC+7.902463001e-6*TC2+5.817764173e-8*TC3;
	R = 4.523601515-33.75714624*TC+3.626406333e-5*TC2+3.87850945e-8*TC3;
	E0 = 0.4093197551-2.271109689e-4*TC-2.86040072e-8*TC2+8.77513e-9*TC3+4.465134e-5*cos(R);
	DLH = aLamda-H;
	DLH2 = 2*DLH;
	aLamda0 = aLamda+2*aL0*sin(DLH)+1.25*aL0*aL0*sin(DLH2);
	dPsi = -17.23*sin(R);
	SL = sin(aLamda0);
	CL = cos(aLamda0);
	CE = cos(E0);
	SE = sin(E0);
	S1 = SL*CE/CL;
	AL = atan(S1);
	if (CL < 0) { AL += 3.141592654;}    // прямое восхождение должно быть в той же четверти что и лямда!!!!
	if (AL < 0) { AL += 6.283185308;}
	ASol = AL+(0.061164*15*dPsi-20.496)*4.8481368e-6;
	DSol = atan(SL*SE/sqrt(CL*CL+SL*SL*CE*CE))-9.936741207e-5*SE*cos(ASol);
	SDSol = sin(DSol);
	CDSol = cos(DSol);

}

double Atmosphera::Kp_Ap (double App) const
{
	double Kpp;
	double kp[28] = {0.0, 0.3333, 0.6667, 1.0, 1.3333, 1.6667, 2.0, 2.3333, 2.6667, 3.0, 3.3333,
		3.6667, 4.0, 4.3333, 4.6667, 5.0, 5.3333, 5.6667, 6.0, 6.3333, 6.6667, 7.0,
		7.3333, 7.6667, 8.0, 8.3333, 8.6667, 9.0};
	double ap[28] = {0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48, 56, 67, 80, 94, 111,
		132, 154, 179, 207, 236, 300, 400};
	if (App <= 400) {
		int i = poi(27,ap,App);
		Kpp = lin2(kp[i],kp[i+1],App,ap[i],ap[i+1]);
	}	 
	else { Kpp = 9.0;}
	return Kpp;
}

double Atmosphera::Ad_D (double days_year) const
{
	double AD;
	double ad[38] = {-0.028, -0.045, -0.047, -0.035, -0.011, 0.022, 0.057, 0.090, 0.114, 0.125,
		0.118, 0.096, 0.060, 0.013, -0.037, -0.086, -0.128, -0.162, -0.185, -0.199,
		-0.202, -0.193, -0.173, -0.140, -0.096, -0.042, 0.015, 0.070, 0.115, 0.144,
		0.155, 0.145, 0.120, 0.084, 0.044, 0.006, -0.023, -0.040};
	double days[38] = {0,10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180,
		190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370};
	int i = poi(37,days,days_year);
	AD = lin2(ad[i],ad[i+1],days_year,days[i],days[i+1]);
	return AD;
}

void Atmosphera::SelectKoef(int aF0, double h)
{
	if (h >= 600) 
	{
		switch (aF0) 
		{
		case 75: 
			{
			double a[3] = {-0.332283e+2, 0.1784, 0.5550636e+3};
			double d[3] = {0.10204e+1, 0.2499e-2, -0.1519e-5};
			double b[3] = {0.7833, 0.2861e-2, -0.1944e-5};
			double c[4] = {-0.44e+1, 0.3024e-1, -0.3283e-4, 0.1012e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5411;
			double e[7] = {-0.38e+1, 0.1972e-1, -0.1833e-4, 0.4938e-8, -0.1200, 0.500e-2, 0.1500e-1};
			double l[3] = {0.1083e-1, 0.6694e-4, -0.4277e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 100: 
			{
			double a[3] = {-0.327731e+2, 0.1899, 0.584245e+3};
			double d[3] = {0.10204e+1, 0.2499e-2, -0.1519e-5};
			double b[3] = {0.7250, 0.26750e-2, -0.1750e-5};
			double c[4] = {-0.44e+1, 0.3024e-1, -0.3283e-4, 0.1012e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5515;
			double e[7] = {-0.37e+1, 0.1783e-1, -0.1506e-4, 0.358e-8, -0.1200, 0.250e-1, 0.7500e-2};
			double l[3] = {0.8317e-2, 0.4837e-4, -0.3039e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 125: 
			{
			double a[3] = {-0.316715e+2, 0.2265, 0.5715408e+3};
			double d[3] = {0.10204e+1, 0.2499e-2, -0.1519e-5};
			double b[3] = {0.6100, 0.2343e-2, -0.1433e-5};
			double c[4] = {-0.8980e+1, 0.4087e-1, -0.395e-4, 0.1123e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.37e+1, 0.175e-1, -0.15e-4, 0.3704e-8, -0.100, 0.2083e-1, 0.6251e-2};
			double l[3] = {0.4667e-2, 0.4606e-4, -0.2722e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 150: 
			{
			double a[3] = {-0.297592e+2, 0.2948, 0.5283389e+3};
			double d[3] = {0.10204e+1, 0.2499e-2, -0.1519e-5};
			double b[3] = {0.9333e-1, 0.3038e-2, -0.1711e-5};
			double c[4] = {-0.8980e+1, 0.4087e-1, -0.395e-4, 0.1123e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.44e+1, 0.1981e-1, -0.1806e-4, 0.4938e-8, -0.100, 0.2750e-1, 0.375e-2};
			double l[3] = {-0.1333e-1, 0.7167e-4, -0.3518e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 175: 
			{
			double a[3] = {-0.288463e+2, 0.314, 0.509723e+2};
			double d[3] = {0.10204e+1, 0.2499e-2, -0.1519e-5};
			double b[3] = {-0.3333, 0.3522e-2, -0.1889e-5};
			double c[4] = {-0.1578e+2, 0.5757e-1, -0.5322e-4, 0.1512e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.36e+1, 0.1694e-1, -0.1556e-4, 0.4321e-8, -0.1200, 0.4116e-1, 0.1433e-2};
			double l[3] = {-0.35e-2, 0.4317e-4, -0.2056e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 200:
			{
			double a[3] = {-0.262994e+2, 0.3817, 0.43422e+3};
			double d[3] = {0.10204e+1, 0.2499e-2, -0.1519e-5};
			double b[3] = {-0.4333, 0.3522e-2, -0.1889e-5};
			double c[4] = {-0.1578e+2, 0.5757e-1, -0.5322e-4, 0.1512e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.36e+1, 0.1653e-1, -0.1528e-4, 0.4321e-8, -0.11, 0.3810e-1, 0.1178e-2};
			double l[3] = {-0.15e-2, 0.3250e-4, -0.1389e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 250:
			{
			double a[3] = {-0.236627e+2, 0.4231, 0.3364318e+3};
			double d[3] = {0.10204e+1, 0.2499e-2, -0.1519e-5};
			double b[3] = {-0.175, 0.2642e-2, -0.1417e-5};
			double c[4] = {0.975e+1, 0.3383e-1, -0.2694e-4, 0.6481e-8};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {0.1, 0.2639e-2, -0.2778e-6, -0.6173e-9, -0.9e-1, 0.3118e-1, 0.9662e-3};
			double l[3] = {-0.8333e-3, 0.2817e-4, -0.1129e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		}
	}
	if (h >= 180 && h < 600) {
		switch (aF0) {
		case 75:
			{
			double a[3] = {-0.155605e+2, 0.8248, 0.769132e+2};
			double d[3] = {-0.1721, 0.5756e-2, -0.3635e-5};
			double b[3] = {-0.8607, 0.7861e-2, -0.5711e-5};
			double c[4] = {0.12791e+1, -0.1576e-1, 0.6499e-4, -0.5145e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5411;
			double e[7] = {-0.2152, 0.4167e-2, 0.1587e-5, -0.1651e-8, -0.1200, 0.500e-2, 0.1500e-1};
			double l[3] = {-0.1698e-1, 0.1448e-3, -0.9535e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 100:
			{
			double a[3] = {-0.156408e+2, 0.7754, 0.679162e+2};
			double d[3] = {-0.1721, 0.5756e-2, -0.3635e-5};
			double b[3] = {-0.7540, 0.6850e-2, -0.4600e-5};
			double c[4] = {0.12791e+1, -0.1576e-1, 0.6499e-4, -0.5145e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5515;
			double e[7] = {-0.2162, 0.4086e-2, 0.1270e-5, -0.1587e-8, -0.1200, 0.250e-1, 0.7500e-2};
			double l[3] = {-0.1249e-1, 0.1111e-3, -0.7706e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 125:
			{
			double a[3] = {-0.152229e+2, 0.7569, 0.558165e+2};
			double d[3] = {-0.1721, 0.5756e-2, -0.3635e-5};
			double b[3] = {-0.5700, 0.5250e-2, -0.300e-5};
			double c[4] = {0.12903e+1, -0.1547e-1, 0.5964e-4, -0.4503e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.1486, 0.3263e-2, 0.3143e-5, -0.3429e-8, -0.100, 0.2083e-1, 0.625e-2};
			double l[3] = {-0.7879e-2, 0.7258e-4, -0.3658e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 150:
			{
			double a[3] = {-0.1697520e+2, 0.6736, 0.854440e+2};
			double d[3] = {-0.1721, 0.5756e-2, -0.3635e-5};
			double b[3] = {-0.4760, 0.4400e-2, -0.2400e-5};
			double c[4] = {0.12903e+1, -0.1547e-1, 0.5964e-4, -0.4503e-7};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.1495, 0.3182e-2, 0.2825e-5, -0.3365e-8, -0.100, 0.2750e-1, 0.375e-2};
			double l[3] = {-0.4882e-2, 0.4692e-4, -0.1742e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 175:
			{
			double a[3] = {-0.1730450e+2, 0.6382, 0.819596e+2};
			double d[3] = {-0.1721, 0.5756e-2, -0.3635e-5};
			double b[3] = {-0.2920, 0.2800e-2, -0.800e-6};
			double c[4] = {0.2057, -0.2912e-2, 0.1739e-4, -0.8565e-8};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.8190e-1, 0.2358e-2, 0.4698e-5, -0.5206e-8, -0.1300, 0.4389e-1, 0.1821e-2};
			double l[3] = {-0.5195e-2, 0.4664e-4, -0.2164e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 200:
			{
			double a[3] = {-0.18266e+2, 0.5797, 0.1009417e+3};
			double d[3] = {-0.1721, 0.5756e-2, -0.3635e-5};
			double b[3] = {-0.3113, 0.2839e-2, -0.1089e-5};
			double c[4] = {0.2057, -0.2911e-2, 0.1739e-4, -0.8565e-8};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.8286e-1, 0.2278e-2, 0.4381e-5, -0.5143e-8, -0.11, 0.3810e-1, 0.1178e-2};
			double l[3] = {-0.5017e-2, 0.4282e-4, -0.2132e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 250:
			{
			double a[3] = {-0.192782e+2, 0.5118, 0.1165792e+3};
			double d[3] = {-0.1721, 0.5756e-2, -0.3635e-5};
			double b[3] = {-0.3307, 0.2878e-2, -0.1378e-5};
			double c[4] = {0.1499e-2, -0.2399e-3, 0.7006e-5, -0.5999e-9};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.2048, 0.3596e-2, -0.1587e-5, 0.3175e-9, -0.9e-1, 0.3117e-1, 0.9662e-3};
			double l[3] = {-0.5455e-2, 0.4273e-4, -0.2273e-7};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		}
	}
	if (h >= 120 && h < 180) {
		switch (aF0) {
		case 75:
			{
			double a[3] = {-0.182991e+2, 0.7009, 0.115343e+3};
			double d[3] = {-0.51019e+1, 0.62581e-1, -0.16721e-3};
			double b[3] = {-0.6828, 0.55762e-2, -0.95238e-6};
			double c[4] = {-0.4384e+1, 0.8063e-1, -0.4925e-3, 0.1042e-5};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5411;
			double e[7] = {-0.1238e+1, 0.1203, -0.6450e-3, 0.1208e-5, -0.1200, 0.500e-2, 0.1500e-1};
			double l[3] = {-0.11975e-1, 0.9983e-4, 0};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 100:
			{
			double a[3] = {-0.181908e+2, 0.7000, 0.1146386e+3};
			double d[3] = {-0.51019e+1, 0.6258e-1, -0.16722e-3};
			double b[3] = {-0.7804, 0.7173e-2, -0.5578e-5};
			double c[4] = {-0.4384e+1, 0.8063e-1, -0.4925e-3, 0.1042e-5};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5515;
			double e[7] = {-0.6683e+1, 0.11014, -0.58375e-3, 0.1083e-5, -0.1200, 0.250e-1, 0.7500e-2};
			double l[3] = {-0.9900e-2, 0.8212e-4, 0.3125e-8};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 125: 
			{
			double a[3] = {-0.185209e+2, 0.6419, 0.1159569e+3};
			double d[3] = {-0.51019e+1, 0.6258e-1, -0.16722e-3};
			double b[3] = {-0.8220, 0.8330e-2, -0.1233e-4};
			double c[4] = {0.9776, -0.257e-1, 0.2027e-3, -0.4708e-6};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.5352e+1, 0.8615e-1, -0.44375e-3, 0.8125e-6, -0.100, 0.2083e-1, 0.6251e-2};
			double l[3] = {-0.7680e-2, 0.6362e-4, 0.3125e-8};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 150: 
			{
			double a[3] = {-0.186522e+2, 0.6124, 0.1164154e+3};
			double d[3] = {-0.51019e+1, 0.6258e-1, -0.16722e-3};
			double b[3] = {-0.7376, 0.7597e-2, -0.1209e-4};
			double c[4] = {0.9776, -0.257e-1, 0.2027e-3, -0.4708e-6};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.4799e+1, 0.77792e-1, -0.4075e-3, 0.77083e-6, -0.100, 0.2750e-1, 0.375e-2};
			double l[3] = {-0.56e-2, 0.4667e-4, 0};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 175:
			{
			double a[3] = {-0.186586e+2, 0.6038, 0.1163531e+3};
			double d[3] = {-0.51019e+1, 0.6258e-1, -0.16722e-3};
			double b[3] = {-0.315, 0.2325e-2, 0.25e-5};
			double c[4] = {-0.5632, 0.5743e-2, -0.925e-5, 0.4167e-8};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.4903e+1, 0.82708e-1, -0.45875e-3, 0.9167e-6, -0.1200, 0.4116e-1, 0.1433e-2};
			double l[3] = {-0.4963e-2, 0.4136e-4, 0};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 200:
			{
			double a[3] = {-0.186495e+2, 0.5974, 0.1162144e+3};
			double d[3] = {-0.51019e+1, 0.6258e-1, -0.16722e-3};
			double b[3] = {-0.5161, 0.5341e-2, -0.8672e-5};
			double c[4] = {0.5632, 0.5743e-2, -0.9250e-5, 0.4167e-8};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.5115e+1, 0.85075e-1, -0.45875e-3, 0.875e-6, -0.11,0.3810e-1, 0.1178e-2};
			double l[3] = {-0.411e-2, 0.3463e-4, -0.3125e-8};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		case 250: 
			{
			double a[3] = {-0.187074e+2, 0.5772, 0.1163395e+3};
			double d[3] = {-0.51019e+1, 0.6258e-1, -0.16722e-3};
			double b[3] = {-0.2531, 0.1929e-2, 0.1495e-5};
			double c[4] = {0.4842, -0.1604e-1, 0.1405e-3, -0.3375e-6};
			double n[2] = {0.1500e+1, 0.6000e-2};
			double fi1  = 0.5585;
			double e[7] = {-0.3137e+1, 0.47742e-1, -0.2275e-3, 0.39583e-6, -0.9e-1, 0.3118e-1, 0.9662e-3};
			double l[3] = {-0.3030e-2, 0.2532e-4, -0.5556e-9};
			equal(3,A,a);
			equal(3,B,b);
			equal(4,C,c);
			equal(2,N,n);
			Fi = fi1;
			equal(3,D,d);
			equal(7,E,e);
			equal(3,L,l);
			break;
			}
		}
	}
}

double Atmosphera::LowAtm(double _h)
{
	const double r = 6356.766;  
	const double A = 0.020046796;
	const double C1 = 0.0341632188;
	const double AX = 28.96442; 
	const double AZ = 28.85;
	const double AY = 0.032622; 
	const double AV = 0.1511;

	double Hst[12]     = {0, 11000, 20000, 32000, 47000, 51000, 71000, 80000, 85000, 94000, 102450, 117770};
	double Tmst[12]    = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 196.65, 186.65, 186.65, 212.0, 380.6};
	double Betamst[12] = {-0.0065, 0.0, 0.001, 0.0028, 0, -0.0028, -0.002, -0.002, 0, 0.003, 0.011, 0.0078};
	double Pst[12]     = {124915236e+3, 37109308e+3, 8977020.69e+3, 1348564490, 145566528, 87858748.9,
					      6547648.79, 1600995.24, 691423.676, 133266.712, 27459.428, 2488.52564 };
	double H, T, Tm, ro;
	int j;
	H = 1000*(r * _h)/(r + _h); //геопотенциальная высота (км)
	H > 117770.0 ? j = 11 : j = poi(10,Hst,H);
	double deltaH = H - Hst[j];
	Tm = Tmst[j] + Betamst[j]*deltaH;
	if (fabs(Betamst[j]) >= 0.000001)
		ro = Pst[j]*exp((1 + C1/Betamst[j])*log(Tmst[j]/Tm));
	else
		ro = Pst[j]*exp(-C1*deltaH/Tmst[j]);
	T = Tm;
	double c2;
	if (_h > 97.5) 
	{
		c2 = AZ - AV*(_h-97.5);
		T = Tm*c2/AX;
		VZV = A*sqrt(T);
		return ro * 1e-12;
	}
	if (_h >= 94 && _h <= 97.5) {
		c2 = AX - AY * (_h-94);
		T = Tm*c2/AX;
		VZV = A*sqrt(T);
		return ro * 1e-12;
	}
	VZV = A*sqrt(T); //скорость звука
	return ro * 1e-12 ;
	// считается скорость звука, еcть предложение когда нужно переделать таким образом, чтобы в функции менялась плотность, а возвращалась скорость звука
}

int Atmosphera::Selection_F0 (int F)
{
	int F0[6] = {75, 100, 125, 150, 200, 250};
	if (F >= 250) { return F = 250;}
	if (F <= 75) { return F = 75;}
	for (int i = 0; i < 6; i++) {
		if (i == 5) {
			F = F0[i];
			break;
		}
		if (abs(F0[i] - F) < abs(F0[i+1] - F)) {
			F = F0[i];
			break;
		}
	}
	return F;
}

void Atmosphera::VarAtm(double H, double Fi, double Lam)
{
	const double C[36][10] = {
		{-0.042e0, 0.001e0, 0.024e0, 0.068e0,-0.006e0, 0.002e0, 0.123e0, 0.059e0, 0.014e0,-0.014e0},
		{-0.04e0 , 0.e0   , 0.e0   , 0.023e0, 0.002e0,-0.003e0, 0.031e0, 0.022e0, 0.003e0,-0.008e0},
		{-0.032e0,-0.002e0, 0.e0   , 0.014e0, 0.002e0,-0.002e0, 0.001e0, 0.009e0, 0.001e0,-0.003e0},
		{-0.019e0,-0.004e0, 0.006e0, 0.006e0, 0.002e0,-0.002e0,-0.019e0,-0.003e0,-0.002e0, 0.002e0},
		{ 0.008e0, 0.001e0, 0.009e0,-0.005e0,-0.001e0,-0.006e0,-0.072e0,-0.014e0,-0.003e0, 0.01e0 },
		{ 0.064e0, 0.002e0,-0.005e0,-0.027e0,-0.002e0,-0.005e0,-0.104e0,-0.015e0,-0.003e0, 0.012e0},
		{ 0.169e0, 0.001e0, 0.005e0,-0.064e0,-0.002e0, 0.009e0,-0.086e0,-0.039e0,-0.002e0, 0.012e0},
		{ 0.093e0, 0.001e0, 0.003e0,-0.066e0, 0.002e0, 0.008e0,-0.069e0,-0.07e0, -0.002e0, 0.008e0},
		{ 0.031e0, 0.002e0, 0.006e0,-0.066e0, 0.004e0, 0.007e0,-0.05e0 ,-0.116e0, 0.004e0, 0.002e0},
		{-0.001e0, 0.003e0,-0.001e0,-0.071e0, 0.004e0, 0.009e0,-0.046e0,-0.154e0, 0.025e0,-0.001e0},
		{-0.009e0, 0.003e0,-0.004e0,-0.076e0, 0.003e0, 0.01e0 ,-0.049e0,-0.172e0, 0.036e0,-0.002e0},
		{-0.003e0, 0.005e0,-0.012e0,-0.08e0 ,-0.009e0, 0.017e0,-0.067e0,-0.211e0, 0.06e0 ,-0.002e0},
		{ 0.019e0, 0.005e0,-0.004e0,-0.094e0,-0.001e0, 0.018e0,-0.07e0, -0.241e0, 0.078e0, 0.001e0},
		{ 0.036e0, 0.002e0, 0.002e0,-0.108e0, 0.007e0, 0.018e0,-0.073e0,-0.267e0, 0.095e0, 0.005e0},
		{ 0.051e0,-0.002e0, 0.006e0,-0.12e0 , 0.014e0, 0.019e0,-0.078e0,-0.293e0, 0.109e0, 0.008e0},
		{ 0.068e0,-0.013e0, 0.009e0,-0.128e0, 0.02e0 , 0.021e0,-0.083e0,-0.316e0, 0.122e0, 0.012e0},
		{ 0.082e0,-0.019e0, 0.009e0,-0.13e0 , 0.026e0, 0.027e0,-0.094e0,-0.337e0, 0.133e0, 0.015e0},
		{ 0.086e0,-0.02e0 , 0.007e0,-0.134e0, 0.03e0 , 0.026e0,-0.099e0,-0.356e0, 0.141e0, 0.019e0},
		{ 0.084e0,-0.019e0, 0.008e0,-0.138e0, 0.037e0, 0.025e0,-0.118e0,-0.391e0, 0.149e0, 0.031e0},
		{ 0.081e0,-0.019e0, 0.005e0,-0.152e0, 0.043e0, 0.024e0,-0.128e0,-0.419e0, 0.171e0, 0.033e0},
		{ 0.087e0,-0.018e0, 0.001e0,-0.165e0, 0.048e0, 0.024e0,-0.125e0,-0.433e0, 0.178e0, 0.035e0},
		{ 0.097e0,-0.022e0,-0.003e0,-0.175e0, 0.052e0, 0.024e0,-0.119e0,-0.444e0, 0.184e0, 0.035e0},
		{ 0.102e0,-0.023e0,-0.013e0,-0.184e0, 0.053e0, 0.029e0,-0.113e0,-0.452e0, 0.187e0, 0.034e0},
		{ 0.105e0,-0.023e0,-0.02e0 ,-0.186e0, 0.054e0, 0.034e0,-0.102e0,-0.453e0, 0.188e0, 0.034e0},
		{ 0.107e0,-0.02e0 ,-0.027e0,-0.189e0, 0.054e0, 0.033e0,-0.088e0,-0.448e0, 0.185e0, 0.033e0},
		{ 0.102e0,-0.011e0,-0.02e0 ,-0.181e0, 0.053e0, 0.032e0,-0.066e0,-0.436e0, 0.172e0, 0.034e0},
		{ 0.059e0,-0.002e0,-0.016e0,-0.016e0, 0.051e0, 0.028e0,-0.045e0,-0.42e0 , 0.158e0, 0.031e0},
		{ 0.019e0, 0.007e0,-0.006e0,-0.125e0, 0.049e0, 0.029e0,-0.023e0,-0.393e0, 0.145e0, 0.026e0},
		{-0.023e0, 0.015e0,-0.004e0,-0.077e0, 0.045e0, 0.015e0, 0.008e0,-0.344e0, 0.129e0, 0.013e0},
		{ 0.008e0, 0.021e0, 0.015e0, 0.008e0, 0.036e0, 0.017e0, 0.039e0,-0.24e0 , 0.11e0 ,-0.02e0 },
		{ 0.051e0,-0.021e0, 0.01e0 , 0.123e0, 0.026e0,-0.023e0, 0.073e0,-0.105e0, 0.091e0,-0.069e0},
		{ 0.075e0,-0.057e0, 0.033e0, 0.214e0, 0.018e0,-0.022e0, 0.153e0,-0.034e0, 0.075e0,-0.094e0},
		{ 0.105e0,-0.07e0 , 0.074e0, 0.24e0 , 0.007e0,-0.021e0, 0.229e0, 0.054e0, 0.053e0,-0.101e0},
		{ 0.137e0,-0.072e0, 0.131e0, 0.243e0,-0.002e0,-0.027e0, 0.272e0, 0.056e0, 0.017e0,-0.07e0 },
		{ 0.179e0,-0.079e0, 0.166e0, 0.207e0,-0.008e0,-0.035e0, 0.265e0, 0.024e0, 0.023e0,-0.046e0},
		{ 0.228e0,-0.09e0 , 0.19e0 , 0.174e0,-0.012e0,-0.041e0, 0.245e0,-0.014e0, 0.024e0,-0.03e0 },
	};

	const double F[29][10] = {
		{  1.5e0  ,- 0.5e0 ,-0.4e0 ,-0.1e0 ,- 0.8e0,- 0.3e0,  4.3e0,  0.5e0,  2.e0 ,  0.3e0},
		{- 0.8e0  ,- 3.e0  , 2.8e0 , 0.5e0 ,- 1.2e0,  0.2e0,  6.3e0, -0.6e0,  1.3e0,  0.3e0},
		{  0.e0   ,- 5.6e0 , 5.8e0 , 0.4e0 ,  0.1e0,- 0.2e0,  8.3e0, -0.8e0,  0.6e0,  0.5e0},
		{  0.8e0  ,- 8.1e0 , 7.9e0 , 0.1e0 ,  0.1e0,  0.e0 , 10.4e0, -1.e0 ,- 0.2e0,  0.6e0},
		{  1.e0   ,- 9.e0  , 15.8e0,-4.e0  ,  1.2e0,  3.2e0,  9.e0 , -1.7e0,  0.1e0,  0.1e0},
		{  0.4e0  ,- 5.4e0 , 12.2e0, 0.5e0 ,  0.e0 ,  0.2e0,  6.1e0, -1.1e0,  1.3e0,- 0.1e0},
		{- 1.6e0  ,- 3.3e0 , 11.e0 , 4.6e0 ,- 1.6e0,- 1.4e0,  5.2e0, -0.2e0,  1.2e0,- 1.e0 },
		{- 4.3e0  ,- 1.8e0 , 10.7e0, 9.8e0 ,- 2.3e0,- 2.1e0,  5.1e0,  1.3e0,  1.3e0,- 1.7e0},
		{- 7.8e0  ,- 1.1e0 , 11.2e0, 13.4e0,- 3.3e0,- 2.6e0,  5.4e0,  3.3e0,  0.8e0,- 2.4e0},
		{-11.4e0  ,- 0.6e0 , 12.4e0, 16.6e0,- 4.3e0,- 4.2e0,  5.9e0,  5.4e0,  0.4e0,- 3.1e0},
		{-20.e0   ,  0.e0  , 13.9e0, 19.6e0,- 4.8e0,- 6.4e0,  6.3e0,  7.5e0,  0.e0 ,- 3.8e0},
		{-17.5e0  ,- 2.5e0 , 15.8e0, 26.9e0,- 6.e0 ,- 6.9e0,  8.5e0, 13.3e0,- 0.4e0,- 4.4e0},
		{-15.e0   ,-15.e0  , 17.9e0, 33.8e0,- 7.4e0,- 7.2e0, 10.8e0, 19.1e0,- 0.8e0,- 5.e0 },
		{-12.5e0  ,-27.5e0 , 19.8e0, 41.2e0,- 8.6e0,- 7.7e0, 13.1e0, 25.e0 ,- 1.3e0,- 5.6e0},
		{- 5.e0   ,-35.e0  , 22.e0 , 48.9e0,-10.e0 ,- 8.e0 , 15.8e0, 30.8e0,- 2.5e0,- 6.7e0},
		{  5.9e0  ,-24.1e0 , 20.5e0, 51.e0 ,-10.e0 ,- 9.5e0, 16.7e0, 32.9e0,- 6.3e0,- 9.6e0},
		{ 16.7e0  ,-13.3e0 , 19.2e0, 53.5e0,-10.e0 ,-10.8e0, 17.5e0, 35.e0 ,-10.e0 ,-12.5e0},
		{ 27.5e0  ,- 2.5e0 , 16.4e0, 53.5e0,- 7.4e0,- 7.9e0, 15.4e0, 32.5e0,-11.7e0,-12.9e0},
		{ 17.5e0  ,  7.5e0 , 10.8e0, 48.5e0,- 5.2e0,- 7.4e0, 13.3e0, 30.e0, -13.3e0,-13.3e0},
		{  7.5e0  ,17.5e0  ,  9.2e0, 36.e0 ,- 2.5e0,- 3.2e0, 10.5e0, 24.e0, -10.e0 ,- 9.5e0},
		{- 1.6e0  , 8.3e0  , 12.5e0, 27.e0 ,- 5.e0 ,  2.5e0, 12.3e0, 20.5e0,-10.e0 ,- 7.8e0},
		{-10.8e0  ,-0.8e0  , 15.5e0, 24.e0 ,- 4.e0 ,  8.5e0, 14.e0 , 17.e0 ,-10.e0 ,- 6.e0 },
		{-10.e0   ,-20.e0  , 19.3e0, 12.1e0,-12.5e0, 16.8e0, 18.8e0, 16.7e0,-10.e0 ,- 4.1e0},
		{  0.e0   ,-30.e0  , 25.8e0,  7.1e0,-14.6e0, 21.2e0, 23.6e0, 16.4e0,-10.e0 ,- 2.1e0},
		{  6.e0   ,-24.e0  , 30.7e0, 10.2e0,-12.2e0, 19.5e0, 28.4e0, 16.1e0,-10.e0 ,- 0.2e0},
		{  5.8e0  ,-11.8e0 , 32.8e0, 19.8e0,- 8.3e0, 14.e0 , 33.1e0, 15.8e0,-10.e0 ,  1.7e0},
		{  5.5e0  ,  0.5e0 , 32.4e0, 28.9e0,-10.9e0, 11.5e0, 37.9e0, 15.5e0,-10.e0 ,  3.6e0},
		{  5.3e0  , 12.8e0 , 29.2e0, 38.4e0,-18.e0 , 11.2e0, 42.7e0, 15.3e0,-10.e0 ,  5.6e0},
		{  9.5e0  , 25.e0  , 26.2e0, 47.5e0,-25.e0 , 11.2e0, 47.5e0, 15.e0 ,-10.e0 ,  7.5e0},
	};
	const double Al[3][3] = {
		{-4.4581, 6.4846,-2.0264},
		{ 4.6442,-8.2564, 3.6122},
		{-1.3141, 2.6281,-1.3141},
	};
	const double Bet[3][3] = {
		{-4.7008, 7.3205,-2.6197},
		{ 5.1122,-9.8690, 4.7568},
		{-1.5165, 3.3262,-1.8097},
	};
	const double MassDHC[36] = {0.e0,3.e0,6.e0,8.e0,10.e0,12.e0,16.e0,20.e0,24.e0,28.e0,30.e0,33.e0,36.e0,39.e0,42.e0,45.e0,48.e0,51.e0,54.e0,57.e0,
					     60.e0,63.e0,66.e0,69.e0,72.e0,75.e0,78.e0,81.e0,85.e0,90.e0,95.e0,100.e0,105.e0,110.e0,115.e0,120.e0};
	const double MassDHF[29] = {0.e0,3.e0,6.e0,9.e0,12.e0,15.e0,18.e0,21.e0,24.e0,27.e0,30.e0,35.e0,40.e0,45.e0,50.e0,55.e0,60.e0,65.e0,70.e0,75.e0,
					     80.e0,85.e0,90.e0,95.e0,100.e0,105.e0,110.e0,115.e0,120.e0};
	unsigned FRow      = 29;
	unsigned CRow      = 36;
	unsigned CFCol     = 10;
	const double Pi3   = 1.0471975;
	const double Pi6   = 0.5235987;	
	double  *RezInterp  = new double[10];
	double   Fi2       = Fi*Fi;
	double   Fi3       = fabs(Fi2*Fi);
	double   Fi4       = Fi2*Fi2;
	double   ArgTrigon = Pi6*NMounth;
	double   SinP6     = sin(ArgTrigon);
	double   CosPi6    = cos(ArgTrigon);
	double   CosPi3    = cos(Pi3*NMounth);
	Fi < 0 ? NMounth   += 6 : NMounth;
	NMounth > 11 ? NMounth -= 12: NMounth;
	unsigned NumH      = poi(FRow,MassDHF,H);
	double   DH        = (H - MassDHF[NumH])/(MassDHF[NumH+1] - MassDHF[NumH]);
	//----------------------------------рассчет средних вариаций ветра------------------------------------------------------------------------------
	LinInterpolMatr(CFCol, F, NumH, DH, RezInterp);
	double  R[3];
	R[0] = RezInterp[0] + RezInterp[1]*CosPi3;
	R[1] = RezInterp[2] + RezInterp[3]*CosPi6 + RezInterp[4]*SinP6 + RezInterp[5]*CosPi3;
	R[2] = RezInterp[6] + RezInterp[7]*CosPi6 + RezInterp[8]*SinP6 + RezInterp[9]*CosPi3;
	double U[3];
	for (unsigned i = 0; i < 3; ++i)
	{
		U[i] = 0;
		for (unsigned j = 0; j < 3; ++j)
			U[i] = U[i] + R[j]*Bet[i][j];
	}
	double WindLatitude = (R[0] + Fi2*U[0] + Fi3*U[1] + Fi4*U[2])*0.001; //широтный ветер
	double WindLongitude = 0; //долготный ветер(пока не считается, нужен при вычислении макимальных вариаций)
	//перевод в ГСК
	double SL = sin(Lam);
	double CL = cos(Lam);
	double SF = sin(Fi);
	double CF = cos(Fi);
	W[0]      = -SL*WindLatitude - SF*CL*WindLongitude;
	W[1]      =  CL*WindLatitude - SF*SL*WindLongitude;
	W[2]      =  CF*WindLongitude;
	//------------------------------------рассчет средней вариации плотности----------------------------------------------------------------------------
	NumH      = poi(CRow,MassDHC,H);
	DH        = (H - MassDHC[NumH])/(MassDHC[NumH+1] - MassDHC[NumH]);
	LinInterpolMatr(CFCol, C, NumH, DH, RezInterp);
	R[0] = RezInterp[0] + RezInterp[1]*CosPi3;
	R[1] = RezInterp[2] + RezInterp[3]*CosPi6 + RezInterp[4]*SinP6 + RezInterp[5]*CosPi3;
	R[2] = RezInterp[6] + RezInterp[7]*CosPi6 + RezInterp[8]*SinP6 + RezInterp[9]*CosPi3;
	for (unsigned i = 0; i < 3; ++i)
	{
		U[i] = 0;
		for (unsigned j = 0; j < 3; ++j)
			U[i] = U[i] + R[j]*Al[i][j];
	}
	DRo = R[0] + Fi2*U[0] + Fi3*U[1] + Fi4*U[2];
	delete [] RezInterp;
}

void Atmosphera::LinInterpolMatr(unsigned Col, const double Matr[][10], unsigned NumArg, double DArg, double *RezArr)
{
	unsigned NumArg1 = NumArg + 1;
	for (unsigned i = 0; i < Col; ++i)
	{
		RezArr[i] = Matr[NumArg][i] + (Matr[NumArg1][i] - Matr[NumArg][i])*DArg;
	}
}

int poi (int n, const double *mass, double par) {  //n - число элементов массива -1, *mass - массив(упорядоченный), par - параметр
	int i = 0;
	if (*mass >= *(mass+1) || fabs(*(mass + n)-par) <= 0.0001)  
		return n;
	while (par > *(mass + i)) 
		i += 1;
	if (i == 0)   
		return i;
	return i - 1;
}



int poi (int n, const int *mass, const double par) {  //n - число элементов массива -1, *mass - массив(упорядоченный), par - параметр
	int i = 0;
	if (*mass >= *(mass+1) || abs(*(mass + n)-par) <= 0.1)  
		return n;
	while (par > *(mass + i)) 
		i += 1;
	if (i == 0)   
		return i;
	return i - 1;
}

double lin2Matr(unsigned nRow, unsigned nCol, double ** matr, unsigned row1, unsigned col1, unsigned row2, unsigned col2, double argRow, double argCol)
{
	if (nRow < row1 || nRow < row2)
		throw logic_error("Index of row greater than matrix dimension");
	if (nCol < col1 || nCol < col2)
		throw logic_error("Index of col greater than matrix dimension");
	double a1 = matr[row1][col1] + (matr[row1][col2] - matr[row1][col1])*argRow;
	double a2 = matr[row2][col1] + (matr[row2][col2] - matr[row2][col1])*argRow;
	return a1 + (a2 - a1)*argCol;
}

double lin3Matr(unsigned n, unsigned m, unsigned l, double *** matr, unsigned i, unsigned j, unsigned k, unsigned i1, unsigned j1, unsigned k1, double di, double dj, double dk)
{
	double x1 = matr[i][j][k] + di*(matr[i1][j][k] - matr[i][j][k]);
	double x2 = matr[i][j1][k] + di*(matr[i1][j1][k] - matr[i][j1][k]);
	double x3 = matr[i][j][k1] + di*(matr[i1][j][k1] - matr[i][j][k1]);
	double x4 = matr[i][j1][k1] + di*(matr[i1][j1][k1] - matr[i][j1][k1]);
	double y1 = x1 + dj*(x2 - x1);
	double y2 = x3 + dj*(x4 - x3);
	return y1 + dk*(y2 - y1);
}