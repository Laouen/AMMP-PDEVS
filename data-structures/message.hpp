#ifndef AMMP_PDEVS_MESSAGE_H
#define AMMP_PDEVS_MESSAGE_H

#include <iostream>
#include <string>
#include "types.hpp"

using namespace std;

struct Adress {

	string _organelle;
	string _cytoplasm;
	string _ensyme_set;
	string _compartment;
	string _enzyme;

	Adress(string other_organelle, string other_cytoplasm, string other_ensyme_set, string other_compartment, string other_enzyme)
	: _organelle(other_organelle), _cytoplasm(other_cytoplasm), _ensyme_set(other_ensyme_set), _compartment(other_compartment), _enzyme(other_enzyme) {}

	Adress()
	: _organelle(""), _cytoplasm(""), _ensyme_set(""), _compartment(""), _enzyme("") {}

	Adress(const Adress& other)
	: _organelle(other._organelle), _cytoplasm(other._cytoplasm), _ensyme_set(other._ensyme_set), _compartment(other._compartment), _enzyme(other._enzyme) {}
	
	Adress(Adress* other)
	: _organelle(other->_organelle), _cytoplasm(other->_cytoplasm), _ensyme_set(other->_ensyme_set), _compartment(other->_compartment), _enzyme(other->_enzyme) {}



	string atModel(string model_type) const{
		
		string result = "";
		if(model_type == "organelle") {
			
			result = _organelle;
		} else if(model_type == "cytoplasm") {
			
			result = _cytoplasm;
		} else if(model_type == "enzyme set") {
			
			result = _ensyme_set;
		} else if(model_type == "compartment") {
			
			result = _compartment;
		} else if(model_type == "enzyme") {
			
			result = _enzyme;
		}

		return result;
	}

	void setModel(string model_type, string new_value) {

		if(model_type == "organelle") {
			
			_organelle = new_value;
		} else if(model_type == "cytoplasm") {
			
			_cytoplasm = new_value;
		} else if(model_type == "enzyme set") {
			
			_ensyme_set = new_value;
		} else if(model_type == "compartment") {
			
			_compartment = new_value;
		} else if(model_type == "enzyme") {
			
			_enzyme = new_value;
		}
	}

	void clear() {
		_organelle 		= "";
		_cytoplasm 		= "";
		_ensyme_set 	= "";
		_compartment 	= "";
		_enzyme 		= "";
	}

};

ostream& operator<<(ostream& os, Adress to) {
	os << "Organelle: " << to._organelle << endl;
	os << "Cytoplasm: " << to._cytoplasm << endl;
	os << "Enzyme set: " << to._ensyme_set << endl;
	os << "Compartment: " << to._compartment << endl;
	os << "Enzyme: " << to._enzyme << endl;
	return os;
}

struct Message {
	
	Adress to;
	string specie;
	integer amount;

	Message(Adress other_to, string other_specie, integer other_amount)
	: to(other_to), specie(other_specie), amount(other_amount) {}

	Message()
	:to(), specie(""), amount(0) {}

	Message(const Message& other)
	: to(other.to), specie(other.specie), amount(other.amount) {}

	Message(Message* other)
	: to(other->to), specie(other->specie), amount(other->amount) {}

	bool sendTo(string model_type, string new_value) {
		to.setModel(model_type, new_value);

		return (to.atModel(model_type) == new_value);
	}

	void clear() {
		to.clear();
		specie = "";
		amount = 0;
	}
};

ostream& operator<<(ostream& os, Message msg) {
	os << "To: " << endl << msg.to << endl;
	os << "Specie: " << msg.specie << endl;
	os << "Amount: " << msg.amount << endl;
	return os;
}

#endif // AMMP_PDEVS_MESSAGE_H