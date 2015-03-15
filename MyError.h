#pragma once
#include <iostream>
class MyError
{
public:
	explicit MyError();
	MyError(char *M, int C);
	void PrintMess();
	~MyError();
	char *Mess; //сообщение об ошибке
	int Code; // код ошибки
};

