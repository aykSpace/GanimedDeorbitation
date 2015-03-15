#pragma once
#include <xstring>

using namespace::std;
namespace CppFlightDynamicOperations
{
	class Engine
	{
	public:

		Engine();
		Engine(int idEngine, string nameEngine, double trust, double specificImpulse, double fuelAmount = 0,
			double maxTimeOfWorking = 0, string typeOfEngine = "Dummy Engine", string comment= "");
		~Engine();

		int Id_Engine;
		string NameEngine;
		double Trust;
		double SpecificImpulse;
		double FuelAmount;
		double MaxTimeOfWorking;
		string TypeOfEngine;
		string Comment;

		bool operator == (const Engine & engine1);

	};
}


