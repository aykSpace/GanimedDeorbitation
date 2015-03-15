#include "MyError.h"


MyError::MyError()
{}

MyError::MyError(char *M, int C)
	:Mess(M), Code(C)
{}

MyError::~MyError()
{
	delete [] Mess;
}

void MyError::PrintMess()
{
	std::cout<<Mess;
}