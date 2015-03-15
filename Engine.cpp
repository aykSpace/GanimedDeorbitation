#include "Engine.h"

using namespace::CppFlightDynamicOperations;
Engine::Engine()
	: Id_Engine(0), NameEngine("Dummy"), Trust(0), SpecificImpulse(0), FuelAmount(0), MaxTimeOfWorking(0), TypeOfEngine("Dummy engine"), Comment("")
{}

Engine::Engine(int idEngine, std::string nameEngine, double trust, double specificImpulse, double fuelAmount, double maxTimeOfWorking, std::string typeOfEngine, std::string comment)
			  :Id_Engine(idEngine), NameEngine(nameEngine), Trust(trust), SpecificImpulse(specificImpulse),
			  FuelAmount(fuelAmount), MaxTimeOfWorking(maxTimeOfWorking), TypeOfEngine(typeOfEngine), Comment(comment)
{}


Engine::~Engine()
{}

bool Engine::operator== (const Engine & engine)
{
	if (Id_Engine != engine.Id_Engine)
		return false;
	if (Trust != engine.Trust)
		return false;
	if (SpecificImpulse != engine.SpecificImpulse)
		return false;
	return true;
}