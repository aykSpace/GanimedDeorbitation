#pragma once
#include <xstring>
namespace CppService
{
	class Printer
	{
	public:
		Printer();
		~Printer();

		static std::string OutFileName;
		static unsigned Precision;
		static void PrintDoubleToFile(double arg1, double arg2, double arg3 );
		static void PrintDoubleToFile(double arg1, double arg2);
		static void PrintDoubleToFile(double arg1);
		static void PrintDoubleToFile (double* mass, int lengthOfMassiv);
		static void PrintDoubleToFile (double* mass, double *mass1, int lengthOfMassives);

		static std::string ReadFileName;
		static void ReadDoublesFromFile(int nCol1, int nCol2, double* resCol1, double *resCol2, int & countStr);
		static void ReadGpzFromFile(double** res1, double** res2);
		static void ReadGpzFromFile2(double* cI0, double* cIJ, double * dIJ); // возвращает столбец сI0, cij, dij
	};
}


