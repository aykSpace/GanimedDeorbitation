#include "Printer.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <xstring>
#include <string>
#include <stdlib.h>

using namespace CppService;

Printer::Printer(){}

Printer::~Printer(){}

void Printer::PrintDoubleToFile(double arg1, double arg2, double arg3)
{
	std::ofstream f(OutFileName, std::ios::out);
	f.setf(std::ios::fixed); // ios::app write data to the end of file
	f.precision(Precision);
	f << arg1 << ", "<< arg2 << ", " << arg3 <<  std::endl;
	f.close();
}

void Printer::PrintDoubleToFile(double arg1, double arg2)
{
	std::ofstream f(OutFileName, std::ios::out);
	f.setf(std::ios::fixed);
	f.precision(Precision);
	f << arg1 << ", "<< arg2 << std::endl;
	f.close();
}

void Printer::PrintDoubleToFile(double arg1)
{
	std::ofstream f(OutFileName, std::ios::out);
	f.setf(std::ios::fixed);
	f.precision(Precision);
	f << arg1 << std::endl;
	f.close();
}

void Printer::PrintDoubleToFile(double* mass, double *mass1, int lengthOfMassives)
{
	std::ofstream f(OutFileName, std::ios::out);
	f.setf(std::ios::fixed);
	f.precision(Precision);
	for(int i = 0; i < lengthOfMassives; ++i)
		f << mass[i] << ", "<< mass1[i] << std::endl;
	f.close();
}

void Printer::PrintDoubleToFile(double* mass, int lengthOfMassiv)
{
	std::ofstream f(OutFileName, std::ios::out);
	f.setf(std::ios::scientific);
	f.precision(Precision);
	for(int i = 0; i < lengthOfMassiv; ++i)
		f << mass[i] << std::endl;
	f.close();
}


void Printer::ReadDoublesFromFile(int nCol1, int nCol2, double* resCol1, double *resCol2, int & countStr)
{
	std::ifstream f(ReadFileName);
	std::string outStr;
	int offset1, offset2;
	if (nCol1 == 1)
		offset1 = 2;
	else
		offset1 = 17 + 14*(nCol1 - 2);
	if (nCol2 == 1)
		offset2 = 2;
	else
		offset2 = 17 + 14*(nCol2 - 2);
	if (f)
	{
		for (int i = 0; i < 12; ++i)
			getline(f,outStr);
		int i = 0;
		while (!f.eof())
		{
			try
			{
				getline(f,outStr);
				auto str1 = outStr.substr(offset1,12);
				auto str2 = outStr.substr(offset2,12);
				double val1 = atof(str1.c_str());
				double val2 = atof(str2.c_str());
				resCol1[i] = val1;
				resCol2[i] = val2;
				i ++;
			}
			catch(std::exception ex)
			{
				std::cout << ex.what() << std::endl;
			}
		}
		f.close();
		countStr = i;
	}
	else
		throw std::logic_error("file not found");
}

void Printer::ReadGpzFromFile(double** res, double** res2)
{
	std::ifstream f(ReadFileName);
	std::string outStr;
	if (f)
	{
		try
		{
			for (int i = 2; i < 37; i++)
			{
				for (int j = 0; j <= i; j++ )
				{
					getline(f,outStr);
					double val = atof(outStr.c_str());
					res[i][j] = val;
				}
			}
			getline(f,outStr);
			for (int i = 2; i < 37; i++)
			{
				for (int j = 0; j <= i; j++ )
				{
					getline(f,outStr);
					double val = atof(outStr.c_str());
					res2[i][j] = val;
				}
			}
		}
		catch(std::exception ex)
		{
			std::cout << ex.what() << std::endl;
		}
		f.close();
	}
	else
		throw std::logic_error("file not found");
}

void Printer::ReadGpzFromFile2(double* cI0, double* cIJ, double * dIJ)
{
	std::ifstream f(ReadFileName);
	std::string outStr;
	int ind1, ind2;
	ind1 = ind2 = 0;
	ind1 ++;
	if (f)
	{
		try
		{
			for (int i = 2; i < 37; i++)
			{
				for (int j = 0; j <= i; j++ )
				{
					getline(f,outStr);
					double val = atof(outStr.c_str());
					if (j == 0)
					{
						cI0[ind1] = val;
						ind1 ++;
					}
					else
					{
						cIJ[ind2] = val;
						ind2 ++;
					}
						
				}
			}
			ind2 = 0;
			getline(f,outStr);
			for (int i = 2; i < 37; i++)
			{
				for (int j = 0; j <= i; j++ )
				{
					getline(f,outStr);
					double val = atof(outStr.c_str());
					if(j != 0)
					{
						dIJ[ind2] = val;
						ind2++;
					}
						
				}
			}
		}
		catch(std::exception ex)
		{
			std::cout << ex.what() << std::endl;
		}
		f.close();
	}
	else
		throw std::logic_error("file not found");
}

std::string Printer::OutFileName = "unnamedFile.txt";
std::string Printer::ReadFileName = "undifine.txt";
unsigned Printer::Precision = 6;
