/* 
 * This file is part of Vlasiator.
 * Copyright 2021 University of Helsinki
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>    

struct variables{

   std::map<std::string,std::string> attributes;

};


struct File{
   std::string filename,xmlTag;
   std::fstream vlsv;
   std::vector<variables> vars;
   std::vector<std::pair< uint64_t,uint64_t> > tags;
   std::vector<std::pair< std::string,std::string> > lookup;
   bool verbose;


   File(std::string fname,int argc,char * argv[]){
      getArgs(argc,argv);
      filename=fname;
      open();
      getXmlTag();
      lookup.push_back(std::make_pair("<MESH ","</MESH>"));
      lookup.push_back(std::make_pair("<VARIABLE","</VARIABLE>"));
      for (auto c:lookup){
         vars.clear();
         tags.clear();
         tags=getVariableKind(c);
         vars.resize(tags.size());
         extractTags();
         print();
      }
   }
  
   bool getArgs(int argc,char *argv[]){

      //std::cout<<std::string(argv[0])=="-h"<<std::endl;
      if (std::string(argv[1])=="-h"){
         this->verbose=false;
      }else if(std::string(argv[1])=="-hv"){
         this->verbose=true;
      }else{
         std::cerr<<"Invalid usage"<<std::endl;
         helpMenu();
         throw std::exception();
        
      }
      
   }

   void helpMenu(){
      
      std::cout << "vlsvdump: Prints header information of vlsv vlasiator files. " << std::endl;
      std::cout << "Usage: ./vlsvdump [options] <file> " << std::endl;
      std::cout << "Options:" << std::endl;
      std::cout << "\t -h \t\t Header information. No data." << std::endl;
      std::cout << "\t -hv \t\t Header information with a bit of spice. No data." << std::endl;
   }



   uint64_t convUInt64(const char* const ptr,const bool& swapEndian) {
      if (swapEndian == false) return *(reinterpret_cast<const uint64_t*>(ptr));
      int index = 0;
      uint64_t tmp = 0;
      char* const ptrtmp = reinterpret_cast<char*>(&tmp);
      for (int i=sizeof(uint64_t)-1; i>=0; --i) {
          ptrtmp[index] = ptr[i];
          ++index;
      }
      return tmp;
   }


   bool open(){
      vlsv.open(filename,std::fstream::in | std::fstream::binary);
      if (!vlsv.good()) {
         std::cerr<<"Could not open file..."<<std::endl;
         throw std::exception();
      }
      return true;
   }
 

   void getXmlTag(){
      uint64_t footer;
      char buffer[sizeof(uint64_t)];
      vlsv.seekp(sizeof(uint64_t));
      vlsv.read(buffer,sizeof(uint64_t));
      footer=convUInt64(buffer,0);
      vlsv.seekg(0,vlsv.end);
      uint64_t end=vlsv.tellg();
      uint64_t footerSize=end-footer;
      char *xmlbuffer = new char [footerSize];
      vlsv.seekg(footer);
      vlsv.read(xmlbuffer,footerSize);
      this->xmlTag=xmlbuffer;
      delete [] xmlbuffer; xmlbuffer=NULL;
      vlsv.close();
   }


   std::vector<std::pair<uint64_t,uint64_t> > getVariableKind(std::pair<std::string,std::string> node){
      uint64_t count = 0, pos = 0;
      std::vector< std::pair<uint64_t,uint64_t> > tags;
      std::string kind=node.first;
      std::string complement=node.second;
      
      while ((pos = xmlTag.find(kind, pos)) != std::string::npos) {
          tags.push_back(std::make_pair(pos,0)  );
          //std::cout<<pos<<std::endl;
          pos+= kind.size(); 
      }
  
      pos=0;
      while ((pos = xmlTag.find(complement, pos)) != std::string::npos) {
          tags.at(count).second = pos ;
          //std::cout<<pos<<std::endl;
          pos+= complement.size(); 
          count++;
      }
      return tags;
  }


   void extractTags(){

      int cnt=0;
      for (auto c:tags){
         std::string line=xmlTag.substr(c.first,c.second-c.first);
         
         std::size_t count = 0, pos = 0 ;
         std::string tab(" ");
         std::string equals("=");

         while ((pos = line.find(equals, 0)) != std::string::npos) {
            uint64_t p0=0,p1=0,p2=0;
            p0=line.find(tab,p0);
            p1=line.find(equals,p0);
            p2=line.find(tab,p1);

            std::string line2=line.substr(p0,p2-p0);
 
            uint64_t _p0=0,_p1=0,_p2=line2.size(),_p3=0;
            _p1=line2.find(equals,0);
            std::string name=line2.substr(_p0,_p1-_p0);
            std::string info=line2.substr(_p1+equals.size(),_p2-_p1);
            info.erase(remove( info.begin(), info.end(), '\"' ),info.end());
            name.erase(remove_if(name.begin(), name.end(), isspace), name.end());
            _p3=info.find(">",0);

            if( _p3 !=std::string::npos){
               info.erase(_p3,info.size()-_p3);
            }
            line.erase(p0, p2-p0);
            this->vars.at(cnt).attributes.insert(std::make_pair(name,info));
         }
         cnt++;
      }
   }


   void  print(){

      std::map<std::string,std::string>::const_iterator it;
      for (auto c:vars){
         it = c.attributes.find("name"); 
         if (it != c.attributes.end()){
            std::cout<<it->first<<"="<<it->second<<std::endl;
            };
         it = c.attributes.find("max_refinement_level"); 
         if (it != c.attributes.end()){
            std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
            };
         it = c.attributes.find("xperiodic"); 
         if (it != c.attributes.end()){
            std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
            };
         it = c.attributes.find("yperiodic"); 
         if (it != c.attributes.end()){
            std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
            };
         it = c.attributes.find("zperiodic"); 
         if (it != c.attributes.end()){
            std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
            };
         if (this->verbose){
            it = c.attributes.find("mesh"); 
            if (it != c.attributes.end()){
               std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
               };
            it = c.attributes.find("unit"); 
            if (it != c.attributes.end()){
               std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
               };
            it = c.attributes.find("arraysize"); 
            if (it != c.attributes.end()){
               std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
               };
            it = c.attributes.find("datasize"); 
            if (it != c.attributes.end()){
               std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
               };
            it = c.attributes.find("datatype"); 
            if (it != c.attributes.end()){
               std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
               };
            it = c.attributes.find("vectorsize"); 
            if (it != c.attributes.end()){
               std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
               };
          }
      }
      
      return;
   }

};

void helpMenu(){
   
   std::cout << "vlsvdump: Prints header information of vlsv vlasiator files. " << std::endl;
   std::cout << "Usage: ./vlsvdump [options] <file> " << std::endl;
   std::cout << "Options:" << std::endl;
   std::cout << "\t -h \t\t Header information. No data." << std::endl;
   std::cout << "\t -hv \t\t Header information with a bit of spice. No data." << std::endl;
}


int main(int argc, char *argv[]){

   if (argc<3){
      std::cerr<<"Invalid Usage\n";
      helpMenu();
      return 1;
   }
   
   std::string fname(argv[2]);


   printf("File %s contains:\n",fname.c_str());
   File vlsv(fname,argc,argv);
   return 0;

}
