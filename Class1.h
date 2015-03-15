#pragma once
class Class1
{
private: 
	static int a;
public:
	Class1();
	~Class1();
	 static int GetA ()
	{
		return a;
	}
	 int Plus()
	 {
		 a = a + 8;
		 return a;
	 }
};

