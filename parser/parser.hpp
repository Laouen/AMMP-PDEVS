#ifndef BOOST_SIMULATION_PDEVS_PARSER_H
#define BOOST_SIMULATION_PDEVS_PARSER_H

// STL
#include <list>
#include <map>
#include <utility> /* pair */
#include <memory> /* shared_ptr */
#include <string> /* atof */

// XML parser include
#include "../tinyXML/tinyxml.h"

// Data structures
#include "../data-structures/types.hpp"
#include "../data-structures/unit_definition.hpp"

using namespace std;


class Parser_t {

private:
	TiXmlDocument               _document;
  map<string, TiXmlElement*>  _models;
  string                      _biomass_ID;
  string                      _p, _e, _c;

  // data
  list<UnitDefinition_t>            _units;
  map<string, string>               _comps;
  map<string, map<string, string>>  _speciesByComps;
  map<string, reaction_info_t>      _reactions;
  map<string, enzyme_t>             _enzymes;
  bool                              _loaded;

  // extern parameters
  Integer_t                          _default_amount;
  shared_ptr<map<string, Integer_t>> _amounts;
  shared_ptr<map<string, Integer_t>> _konSTPs;
  shared_ptr<map<string, Integer_t>> _konPTSs;
  Integer_t                          _cell_weight;
  long double                        _normalization;

  // private methods
  Integer_t getBiomassStoichiometryFrom(double);
  Integer_t getStoichiometryFrom(double);
  Address_t getReactionAddress(const SetOfMolecules_t&, const SetOfMolecules_t&, string);
  bool loadFile(const char *);
  void setReactionSctry(TiXmlElement *, SetOfMolecules_t&);
  string getGAValue(TiXmlElement *);
  vector<string> getEnzymesHandlerIDs(TiXmlElement *r);
  vector<string> proccesExpresion(vector<string>& tokens, int start, int end);
  int getSubExpresion(const vector<string>& tokens, int start);
  void getEnzymeCompartments(const enzyme_t& en, vector<string>& compartments);

public:
  // Constructors
  Parser_t(const char *filename);
  
  bool loadFile();
  void loadExternParameters(
    Integer_t other_default_amount,
    shared_ptr<map<string, Integer_t>> other_amounts,
    shared_ptr<map<string, Integer_t>> other_konSTPs,
    shared_ptr<map<string, Integer_t>> other_konPTSs,
    long double other_cw,
    Integer_t other_n,
    string other_p,
    string other_e,
    string other_c,
    string other_bID
  );

  // methods to get the information from the current SBML file.
  string specieComp(string);
  pair<string, string> getCompAndSubComp(const SetOfMolecules_t&, const SetOfMolecules_t&);
  list<UnitDefinition_t>& getUnitDefinitions();
  map<string, string>& getCompartments();
  map<string, map<string, string>>& getSpecieByCompartments();
  map<string, reaction_info_t>& getReactions();
  map<string, enzyme_t>& getEnzymes();
  vector<string> getReactionIDs();
  reaction_info_t getBiomass();
  map<string, map<string, enzyme_t>> getEnzymesByCompartments();
};

#endif // BOOST_SIMULATION_PDEVS_PARSER_H