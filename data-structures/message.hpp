#ifndef AMMP_PDEVS_MESSAGE_H
#define AMMP_PDEVS_MESSAGE_H

#include <iostream>
#include <string>

using namespace std;

struct Adress {

	string _organelle;
	string _cytoplasm;
	string _ensyme_set;
	string _compartment;
	string _enzyme;

	Adress(string other_organelle, string other_cytoplasm, string other_ensyme_set, string other_compartment, string other_enzyme)
	: _organelle(other_organelle), _cytoplasm(other_cytoplasm), _ensyme_set(other_ensyme_set), _compartment(other_compartment), _enzyme(other_enzyme) {}

	string at(string model_type){
		
		string result = "";
		if(model_type == "organelle") {
			
			result = _organelle;
		} else if(model_type == "cytoplasm") {
			
			result = _cytoplasm;
		} else if(model_type == "ensyme_set") {
			
			result = _ensyme_set;
		} else if(model_type == "compartment") {
			
			result = _compartment;
		} else if(model_type == "enzyme") {
			
			result = _enzyme;
		}

		return result;
	}

};

ostream& operator<<(ostream& os, Adress to) {
	os << "Organelle: " << to._organelle << endl;
	os << "Cytoplasm: " << to._cytoplasm << endl;
	os << "Ensyme_set: " << to._ensyme_set << endl;
	os << "Compartment: " << to._compartment << endl;
	os << "Enzyme: " << to._enzyme << endl;
	return os;
}

struct Message {
	
	Adress to;
	string specie;
	double amount;

	Message(Adress other_to, string other_specie, double other_amount)
	: to(other_to), specie(other_specie), amount(other_amount) {}
};

ostream& operator<<(ostream& os, Message msg) {
	os << "To: " << endl << msg.to << endl;
	os << "Specie: " << msg.specie << endl;
	os << "Amount: " << msg.amount << endl;
	return os;
}

#endif // AMMP_PDEVS_MESSAGE_H