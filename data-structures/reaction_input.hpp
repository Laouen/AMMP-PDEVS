#ifndef AMMP_PDEVS_REACTION_INPUT_H
#define AMMP_PDEVS_REACTION_INPUT_H

using namespace std;

/******** Basic types ************/
typedef pair<string, double> molecule;

/******** Vector types ************/
typedef vector<molecule> vectorOfMolecules;

struct reaction_input {
  
  string name;
  double rate;
  vectorOfMolecules reactants;
  vectorOfMolecules products;
  vectorOfMolecules enzymes;

  reaction_input(string other_name, double other_rate, vector<string> other_reactants, vector<string> other_products, vector<string> other_enzymes)
  : name(other_name),
  rate(other_rate) {

    for (vector<string>::iterator i = other_reactants.begin(); i != other_reactants.end(); ++i){
      reactants.push_back(make_pair(*i, double(0)));
    }

    for (vector<string>::iterator i = other_enzymes.begin(); i != other_enzymes.end(); ++i){
      enzymes.push_back(make_pair(*i, double(0)));
    }

    for (vector<string>::iterator i = other_products.begin(); i != other_products.end(); ++i){
      products.push_back(make_pair(*i, double(0)));
    }
  }
};

#endif // AMMP_PDEVS_REACTION_INPUT_H