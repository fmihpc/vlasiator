#include "parameters.h"

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <limits>

#include <iomanip>
#include <string>
#include <fstream>
#include <iostream>

#ifndef NAN
#define NAN 0
#endif

using namespace std;

typedef Parameters P;

// Define static members:
Real P::xmin = NAN;
//Real P::xmax = NAN;
Real P::ymin = NAN;
//Real P::ymax = NAN;
Real P::zmin = NAN;
//Real P::zmax = NAN;
Real P::dx_ini = NAN;
Real P::dy_ini = NAN;
Real P::dz_ini = NAN;

Real P::vxmin = NAN;
Real P::vxmax = NAN;
Real P::vymin = NAN;
Real P::vymax = NAN;
Real P::vzmin = NAN;
Real P::vzmax = NAN;

uint P::xcells_ini = numeric_limits<uint>::max();
uint P::ycells_ini = numeric_limits<uint>::max();
uint P::zcells_ini = numeric_limits<uint>::max();
uint P::vxblocks_ini = numeric_limits<uint>::max();
uint P::vyblocks_ini = numeric_limits<uint>::max();
uint P::vzblocks_ini = numeric_limits<uint>::max();

Real P::q = NAN;
Real P::m = NAN;
Real P::q_per_m = NAN;
Real P::t = 0;
Real P::dt = NAN;
luint P::tstep = 0;
luint P::tstep_min = 0;
luint P::tsteps = 0;
luint P::saveRestartInterval = numeric_limits<uint>::max();
luint P::diagnInterval = numeric_limits<uint>::max();

bool P::save_spatial_grid;
bool P::save_velocity_grid;

std::string P::solar_wind_file;

std::vector<Real> P::save_spatial_cells_x;
std::vector<Real> P::save_spatial_cells_y;
std::vector<Real> P::save_spatial_cells_z;

uint P::transmit = 0;

bool P::recalculateStencils = true;
bool P::propagateVlasov = true;
bool P::propagateField = true;

bool P::periodic_x = false;
bool P::periodic_y = false;
bool P::periodic_z = false;

// Handles parameter processing from the user
namespace PO = boost::program_options;

static bool initialized = false;
static PO::options_description* descriptions = NULL;
static PO::variables_map* variables = NULL;

static map<string,string> options;
static map<string,bool>   isOptionParsed;
static map< string,vector<string> > vectorOptions;
static map<string,bool>   isVectorOptionParsed;

static string global_config_file_name = "";
static string user_config_file_name = "";
static string run_config_file_name = "";



int Readparameters::argc;
char** Readparameters::argv;
int Readparameters::rank;
MPI_Comm Readparameters::comm;

/** Constructor for class ReadParameters. The constructor defines some 
 * default parameters and parses the input files. 
 * @param argc Command line argc.
 * @param argv Command line argv.
 * @param mpicomm Communicator for processes that will access options 
 */
Readparameters::Readparameters(int argc, char* argv[],MPI_Comm mpicomm) {
    Readparameters::argc = argc;
    Readparameters::argv = argv;
    MPI_Comm_dup(mpicomm,&(Readparameters::comm));
    MPI_Comm_rank(Readparameters::comm,&(Readparameters::rank));

    if (Readparameters::rank==0){
        if (initialized == false) {
            descriptions = new PO::options_description("Usage: main [options (options given on the command line override options given everywhere else)], where options are:");
            variables = new PO::variables_map;
            initialized = true;
            addDefaultParameters();
        }
    }
    else{
        descriptions = NULL;
        variables = NULL;
    }
    //send as Int as MPI_BOOL is only in C++ bindings
    int init_int=initialized;
    MPI_Bcast(&init_int,1,MPI_INT,0,Readparameters::comm);
    initialized=(init_int==1);
   
}


/** Add a new input parameter to Readparameters. Note that Readparameters::parse must be called
 * in order for the input file(s) to be re-read. This functions only needs to be called by root process.
 * Other processes can call it but those calls have no effect.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Descrin for the parameter.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Readparameters::add(const string& name,const string& desc,const std::string& defValue) {
    if (initialized == false) return false;
    
    if(rank==0){
        options[name]="";
        isOptionParsed[name]=false;
        descriptions->add_options()(name.c_str(), PO::value<string>(&(options[name]))->default_value(defValue), desc.c_str());
    }
   
    return true;
}


/** Add a new input parameter to Readparameters. Note that Readparameters::parse must be called
 * in order for the input file(s) to be re-read. This functions only needs to be called by root process.
 * Other processes can call it but those calls have no effect.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Descrin for the parameter.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Readparameters::add(const string& name,const string& desc,const bool& defValue) {
    if (initialized == false) return false;
    stringstream ss;
    ss<<defValue;

    if(rank==0){
        options[name]="";
        isOptionParsed[name]=false;
        descriptions->add_options()(name.c_str(), PO::value<string>(&(options[name]))->default_value(ss.str()), desc.c_str());
    }
    return true;
}

/** Add a new input parameter to Readparameters. Note that Readparameters::parse must be called
 * in order for the input file(s) to be re-read. This functions only needs to be called by root process.
 * Other processes can call it but those calls have no effect.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Descrin for the parameter.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Readparameters::add(const string& name,const string& desc,const int& defValue) {
    if (initialized == false) return false;
    stringstream ss;
    ss<<defValue;

    if(rank==0){
        options[name]="";
        isOptionParsed[name]=false;
        descriptions->add_options()(name.c_str(), PO::value<string>(&(options[name]))->default_value(ss.str()), desc.c_str());
    }
    return true;
}


/** Add a new input parameter to Readparameters. Note that Readparameters::parse must be called
 * in order for the input file(s) to be re-read. This functions only needs to be called by root process.
 * Other processes can call it but those calls have no effect.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Descrin for the parameter.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Readparameters::add(const string& name,const string& desc,const unsigned int& defValue) {
    if (initialized == false) return false;
    stringstream ss;
    ss<<defValue;

    if(rank==0){
        options[name]="";
        isOptionParsed[name]=false;
        descriptions->add_options()(name.c_str(), PO::value<string>(&(options[name]))->default_value(ss.str()), desc.c_str());
    }
    return true;
}

/** Add a new input parameter to Readparameters. Note that Readparameters::parse must be called
 * in order for the input file(s) to be re-read. This functions only needs to be called by root process.
 * Other processes can call it but those calls have no effect.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Descrin for the parameter.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Readparameters::add(const string& name,const string& desc,const float& defValue) {
    if (initialized == false) return false;
    stringstream ss;
    //set full precision
    ss<<setprecision(numeric_limits<float>::digits10 + 1) <<defValue;
    if(rank==0){
        options[name]="";
        isOptionParsed[name]=false;
        descriptions->add_options()(name.c_str(), PO::value<string>(&(options[name]))->default_value(ss.str()), desc.c_str());
    }
    return true;
}

/** Add a new input parameter to Readparameters. Note that Readparameters::parse must be called
 * in order for the input file(s) to be re-read. This functions only needs to be called by root process.
 * Other processes can call it but those calls have no effect.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Descrin for the parameter.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Readparameters::add(const string& name,const string& desc,const double& defValue) {
    if (initialized == false) return false;
    stringstream ss;

    //set full precision
    ss<<setprecision(numeric_limits<double>::digits10 + 1) <<defValue;
    if(rank==0){
        options[name]="";
        isOptionParsed[name]=false;
        descriptions->add_options()(name.c_str(), PO::value<string>(&(options[name]))->default_value(ss.str()), desc.c_str());
    }
    return true;
}


/** Add a new composing input parameter to Readparameters. Note that Readparameters::parse must be called
 * in order for the input file(s) to be re-read. This functions only needs to be called by root process.
 * Other processes can call it but those calls have no effect.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Descrin for the parameter.
 * @return If true, the new parameter was added successfully.
 */
bool Readparameters::addComposing(const string& name,const string& desc) {
    if (initialized == false) return false;

    if(rank==0){
        isVectorOptionParsed[name]=false;
        descriptions->add_options()(name.c_str(), PO::value<vector<string> >(&(vectorOptions[name]))->composing(), desc.c_str());
    }
    return true;
}

//add names of input files
bool Readparameters::addDefaultParameters() {
    typedef Parameters P;
    
    if (initialized == false) return false;
    if(rank==0){
        descriptions->add_options()
            ("help", "print this help message");
        
        // Parameters which set the names of the configuration file(s):
        descriptions->add_options()
            ("global_config", PO::value<string>(&global_config_file_name)->default_value(""),"read options from the global configuration file arg (relative to the current working directory). Options given in this file are overridden by options given in the user's and run's configuration files and by options given in environment variables (prefixed with MAIN_) and the command line")
            ("user_config", PO::value<string>(&user_config_file_name)->default_value(""), "read options from the user's configuration file arg (relative to the current working directory). Options given in this file override options given in the global configuration file. Options given in this file are overridden by options given in the run's configuration file and by options given in environment variables (prefixed with MAIN_) and the command line")
            ("run_config", PO::value<string>(&run_config_file_name)->default_value(""), "read options from the run's configuration file arg (relative to the current working directory). Options given in this file override options given in the user's and global configuration files. Options given in this override options given in the user's and global configuration files. Options given in this file are overridden by options given in environment variables (prefixed with MAIN_) and the command line");

        
    }
    
    return true;
}
    


/** Deallocate memory reserved by Parameters.
 * @return If true, class Parameters finalized successfully.
 */
bool Readparameters::finalize() {
    
    if(rank==0){
        delete descriptions;
        delete variables;
        descriptions = NULL;
        variables = NULL;
    }
    initialized = false;
    return true;
   
}



/** Get the value of the given parameter added with addComposing(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,std::vector<std::string>& value) {
    if(vectorOptions.find(name) != vectorOptions.end() ){ //check if it exists
        value = vectorOptions[name];
        return true;
    }
    
    return false;
 }

/** Get the value of the given parameter added with addComposing(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,std::vector<int>& value) {
    vector<string> stringValue;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,stringValue);
    if (ret) {
        for (vector<string>::iterator i = stringValue.begin(); i!=stringValue.end(); ++i) {
            value.push_back(lexical_cast<int>(*i));
        }
    }
    return ret;
}



/** Get the value of the given parameter added with addComposing(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,std::vector<float>& value) {
    vector<string> stringValue;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,stringValue);
    if (ret) {
        for (vector<string>::iterator i = stringValue.begin(); i!=stringValue.end(); ++i) {
            value.push_back(lexical_cast<float>(*i));
        }
    }
    return ret;
}


/** Get the value of the given parameter added with addComposing(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,std::vector<double>& value) {
    vector<string> stringValue;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,stringValue);
    if (ret) {
        for (vector<string>::iterator i = stringValue.begin(); i!=stringValue.end(); ++i) {
            value.push_back(lexical_cast<double>(*i));
        }
    }
    return ret;
}







/** Get the value of the given parameter added with add(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,std::string& value) {
    if(options.find(name) != options.end() ){ //check if it exists
        value = options[name];
        return true;
    }
    return false;
}


/** Get the value of the given parameter added with add(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,bool& value) {
    string sval;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,sval);
    if (ret) {
        value = lexical_cast<bool>(sval);
    }
    return ret;
}

/** Get the value of the given parameter added with add(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,int& value) {
    string sval;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,sval);
    if (ret) {
        value = lexical_cast<int>(sval);
    }
    return ret;
}

/** Get the value of the given parameter added with add(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,unsigned int& value) {
    string sval;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,sval);
    if (ret) {
        value = lexical_cast<unsigned int>(sval);
    }
    return ret;
}

/** Get the value of the given parameter added with add(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,unsigned long& value) {
    string sval;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,sval);
    if (ret) {
        value = lexical_cast<unsigned long>(sval);
    }
    return ret;
}





/** Get the value of the given parameter added with add(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,float& value) {
    string sval;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,sval);
    if (ret) {
        value = lexical_cast<float>(sval);
    }
    return ret;
}


/** Get the value of the given parameter added with add(). This may be called after having called Parse, and it may be called
 * by any process, in any order.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Readparameters::get(const std::string& name,double& value) {
    string sval;
    bool ret;
    using boost::lexical_cast;
    ret=Readparameters::get(name,sval);
    if (ret) {
        value = lexical_cast<double>(sval);
    }
    return ret;
}





/** Write the descriptions of known input options to standard output if 
 * an option called "help" has been read.
 * @return If true, option called "help" was found and descriptions were written 
 * to standard output.
 */
bool Readparameters::helpMessage() {
    if(rank==0){
        if (variables->count("help") > 0) {
            cout << *descriptions <<endl;
            return true;
        }
        return false;
    }
    return true;
}



/** Query is Parameters has initialized successfully.
 * @return If true, Parameters is ready for use.
 */
bool Readparameters::isInitialized() {return initialized;}

/** Request Parameters to reparse input file(s). This function needs 
 * to be called after new options have been added via Parameters:add functions.
 * Otherwise the values of the new options are not read. This is a collective function, all processes have to all it.
 * @return If true, input file(s) were parsed successfully.
 */
bool Readparameters::parse() {
    if (initialized == false) return false;
    // Tell Boost to allow undescribed options (throws exception otherwise)
   
    if(rank==0){
        const bool ALLOW_UNKNOWN = true;
        // Read options from command line:
        PO::store(PO::parse_command_line(argc, argv, *descriptions), *variables);
        PO::notify(*variables);
        // Read options from environment variables:
        PO::store(PO::parse_environment(*descriptions, "MAIN_"), *variables);
        PO::notify(*variables);
        // Read options from run config file:
        if (run_config_file_name.size() > 0) {
            ifstream run_config_file(run_config_file_name.c_str(), fstream::in);
            if (run_config_file.good() == true) {
                PO::store(PO::parse_config_file(run_config_file, *descriptions, ALLOW_UNKNOWN), *variables);
                PO::notify(*variables);
                run_config_file.close();
            } else {
                cerr << "Couldn't open or read run config file " << run_config_file_name << endl;
                abort();
            }
        }
        // Read options from user config file:
        if (user_config_file_name.size() > 0) {
            ifstream user_config_file(user_config_file_name.c_str(), fstream::in);
            if (user_config_file.good() == true) {
                PO::store(PO::parse_config_file(user_config_file, *descriptions, ALLOW_UNKNOWN), *variables);
                PO::notify(*variables);
                user_config_file.close();
            } else {
                cerr << "Couldn't open or read user config file " << user_config_file_name << endl;
                abort();
            }
        }
        // Read options from global config file:
        if (global_config_file_name.size() > 0) {
            ifstream global_config_file(global_config_file_name.c_str(), fstream::in);
            if (global_config_file.good() == true) {
                PO::store(PO::parse_config_file(global_config_file, *descriptions, ALLOW_UNKNOWN), *variables);
                PO::notify(*variables);
                global_config_file.close();
            } else {
                cerr << "Couldn't open or read global config file " << global_config_file_name << endl;
                abort();
            }
        }
        
    }

    //Check if the user has specified --help    
    bool hasHelpOption=helpMessage();
    MPI_Bcast(&hasHelpOption,sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
    if(hasHelpOption){
        MPI_Finalize();
        exit(0);
    }

    //Require that there is a run config file. There are so many options so it is unlikely
    //that one would like to define all on the command line, or only in global/user run files
    //If no arguments are given the program (currently r155) would crash later on with nasty error messages
    bool hasRunConfigFile=(run_config_file_name.size() > 0);
    MPI_Bcast(&hasRunConfigFile,sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);    
    if(!hasRunConfigFile){
        if(Readparameters::rank==0){
            cout << "Run config file required. Use --help to list all options" <<endl;
        }
        MPI_Finalize();
        exit(0);
    }
    
    int nOptionsToBroadcast;
    int vectorSize;
    const int maxStringLength=1024;
    char value[maxStringLength];
    char name[maxStringLength];
    value[maxStringLength-1]='\0';
    name[maxStringLength-1]='\0';


    /*
      loop through options and broadcast all options from root rank to the others
      Separate bcasts not optimal from performance point of view, but parse is normally just
      called a few times so it should not matter
    */
  
    //count number of options not parsed/broarcasted previously
    if(rank==0){
        nOptionsToBroadcast=0;
        for( map<string,bool>::iterator ip=isOptionParsed.begin(); ip!=isOptionParsed.end();++ip){
            if(! ip->second) nOptionsToBroadcast++;
        }
    }

    MPI_Bcast(&nOptionsToBroadcast,1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank==0){
        //iterate through map and bcast cstrings of key/value pairs not parsed before
        for( map<string,string>::iterator p=options.begin(); p!=options.end();++p){
           
            if(! isOptionParsed[p->first]) {
                strncpy(name,p->first.c_str(),maxStringLength-1);
                strncpy(value,p->second.c_str(),maxStringLength-1);
                MPI_Bcast(name,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
                MPI_Bcast(value,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
                isOptionParsed[p->first]=true; 
            }
        }
    }
    else{
        //reveive new options
        for(int p=0;p<nOptionsToBroadcast;p++){
            MPI_Bcast(name,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
            MPI_Bcast(value,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
            string sName(name);
            string sValue(value);
            options[sName]=sValue;
        }
    }
    
 
    /*
      loop through vector options and broadcast all vector options from root rank to the others
      Separate bcasts not optimal from performance point of view, but parse is normally just
      called a few times so it should not matter
    */
  
    //count number of vector options not parsed/broarcasted previously
    if(rank==0){
        nOptionsToBroadcast=0;
        for( map<string,bool>::iterator ip=isVectorOptionParsed.begin(); ip!=isVectorOptionParsed.end();++ip){
            if(! ip->second) nOptionsToBroadcast++;
        }
    }

    //root broadcasts its new vector values
    MPI_Bcast(&nOptionsToBroadcast,1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank==0){
        //iterate through map and bcast cstrings of key/value pairs not parsed before
        for( map< string,vector<string> >::iterator p=vectorOptions.begin(); p!=vectorOptions.end();++p){
            if(! isVectorOptionParsed[p->first]) {
                strncpy(name,p->first.c_str(),maxStringLength-1);
                MPI_Bcast(name,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
                vectorSize=vectorOptions[p->first].size();
                MPI_Bcast(&vectorSize,1,MPI_INT,0,MPI_COMM_WORLD);
                for( vector<string>::iterator v=vectorOptions[p->first].begin(); v!=vectorOptions[p->first].end();++v){
                    strncpy(value,v->c_str(),maxStringLength-1);
                    MPI_Bcast(value,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
                }
                isVectorOptionParsed[p->first]=true; 
            }
        }
    }
    else{
        //others receive new options
        for(int p=0;p<nOptionsToBroadcast;p++){
            MPI_Bcast(name,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
            string sName(name);
            MPI_Bcast(&vectorSize,1,MPI_INT,0,MPI_COMM_WORLD);
            for (int i=0;i<vectorSize;i++){
                MPI_Bcast(value,maxStringLength,MPI_CHAR,0,MPI_COMM_WORLD);
                string sValue(value);
                vectorOptions[sName].push_back(sValue);
            }
        }
    }

   
    return true;
}





