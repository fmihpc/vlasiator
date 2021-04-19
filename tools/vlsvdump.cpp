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
#include <sstream>
#include <cstring>
#include <map>    
#include "zstd.h"


struct variables{

   std::map<std::string,std::string> attributes;

};


struct File{
   std::string filename,xmlTag;
   std::fstream vlsv;
   std::vector<variables> vars;
   std::vector<variables> extract;
   std::vector<std::pair< uint64_t,uint64_t> > tags;
   std::vector<std::pair< std::string,std::string> > lookup;
   bool verbose=false;
   bool show=false;
   bool write=false;
   bool compress=false;
   bool decompress=false;
   float filesize;
   std::vector<std::string> key;


   File(std::string fname,int argc,char * argv[]){
      getArgs(argc,argv);
      filename=fname;
      if (!this->decompress){
         analyse(this->show);
         if (this->write){
            writeExtractedVars();
         }
      }
      else{
         //lz4_Decompress();
         zDecompress();
      }
   }


   bool analyse(bool doPrint){

      this->filesize=get_file_size(filename);
      open();
      getXmlTag();
      lookup.push_back(std::make_pair("<MESH ","</MESH>"));
      //lookup.push_back(std::make_pair("<BLOCKIDS","</BLOCKIDS>"));
      //lookup.push_back(std::make_pair("<BLOCKSPECELL","</BLOCKSPECELL>"));
      //lookup.push_back(std::make_pair("<BLOCKVARIABLE","</BLOCKVARIABLE>"));
      lookup.push_back(std::make_pair("<VARIABLE","</VARIABLE>"));
      for (auto c:lookup){
         vars.clear();
         tags.clear();
         tags=getVariableKind(c);
         vars.resize(tags.size());
         extractTags();
         
         if (this->write){
            //Let's look for variables to extract
            std::map<std::string,std::string>::const_iterator it;
            for (auto c:vars){
               it = c.attributes.find("name"); 
               if (it!=c.attributes.end()){
                  auto itkey = find(this->key.begin(), this->key.end(), it->second);
                  if(itkey!=this->key.end()){
                     this->extract.push_back(c);
                     this->key.erase(itkey);
                  }
               }
            }
         }

         if (doPrint){
            print();
         }
      }
   }


   bool getArgs(int argc,char *argv[]){

      if (std::string(argv[1])=="-h"){
         this->verbose=false;
         this->show=true;
      }else if(std::string(argv[1])=="-hv"){
         this->verbose=true;
         this->show=true;
      }else if(std::string(argv[1])=="-d"){
         this->decompress=true;
      }else if(std::string(argv[1])=="-e"){
         if (argc<4){
            std::cerr<<"Detected -e but no variables specified"<<std::endl;
            helpMenu();
            throw std::exception();
         }
         std::string temp(argv[3]);
         temp+=",";
         std::stringstream lookup(temp);
         while(lookup.good()) {
            std::string substr;
            getline(lookup, substr,','); 
            this->key.push_back(substr);
         }

         if (this->key.size()>1){
            this->write=true;
         }else{
            std::cerr<<"-e detected but could not read variable list.."<<std::endl;
            helpMenu();
            throw std::exception();
         }
      }else if(std::string(argv[1])=="-ez"){
         if (argc<4){
            std::cerr<<"Detected -e but no variables specified"<<std::endl;
            helpMenu();
            throw std::exception();
         }
         std::string temp(argv[3]);
         temp+=",";
         std::stringstream lookup(temp);
         while(lookup.good()) {
            std::string substr;
            getline(lookup, substr,','); 
            this->key.push_back(substr);
         }
         if (this->key.size()>1){
            this->write=true;
            this->compress=true;
         }else{
            std::cerr<<"-e detected but could not read variable list.."<<std::endl;
            helpMenu();
            throw std::exception();
         }
      }else{
         std::cerr<<"Invalid usage"<<std::endl;
         helpMenu();
         throw std::exception();
      }
      
   }


   void helpMenu(){
      
      std::cout << "vlsvdump: Prints header information of vlsv vlasiator files. " << std::endl;
      std::cout << "Usage: ./vlsvdump [options] <file>  [variablelist] " << std::endl;
      std::cout << "Options:" << std::endl;
      std::cout << "\t -h \t\t Header information. No data." << std::endl;
      std::cout << "\t -hv \t\t Header information with a bit of spice. No data." << std::endl;
      std::cout << "\t -e \t\t Extracts variable(s) to binary file(s). Variable(s) should be delimited by comma"<<std::endl;
      std::cout << "\t -ez \t\t Extracts variable(s) and compresses using zstd to binary file(s). Variable(s) should be delimited by comma"<<std::endl;
      std::cout << "\t -d  \t\t Decompresses file compressed by vlsvdump using zstd"<<std::endl;
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
 

   bool writeExtractedVars(){
      
      for (auto c:extract){
      
         //Get stats
         uint64_t arraysize,datasize,vectorsize,tagLocation;
         std::map<std::string,std::string>::const_iterator it_0,it_1,it_2,it_3,it_4;
         it_0 = c.attributes.find("arraysize"); 
         it_1 = c.attributes.find("datasize"); 
         it_2 = c.attributes.find("vectorsize");
         it_3 = c.attributes.find("tag_location");
         it_4 = c.attributes.find("name");
         bool doKeepGoing= (it_0!=c.attributes.end() && 
                            it_1!=c.attributes.end() && 
                            it_2!=c.attributes.end() && 
                            it_3!=c.attributes.end() &&
                            it_4!=c.attributes.end() 
                            );
         if (!doKeepGoing){return false;}
         arraysize=std::stoul(it_0->second,nullptr,10);
         datasize=std::stoul(it_1->second,nullptr,10);
         vectorsize=std::stoul(it_2->second,nullptr,10);
         tagLocation=std::stoul(it_3->second,nullptr,10);

         //Read in  from original vlsv file
         open();
         uint64_t readIn=arraysize*datasize*vectorsize;
         char *buffer = new char [readIn];
         vlsv.seekg(tagLocation);
         vlsv.read(buffer,readIn);
         vlsv.close();

         //Compress if asked for
         std::string outname= this->filename + "_"  +it_4->second+".bin";
         std::replace(outname.begin(), outname.end(), '/', '_');
         if (this->compress){
            //writeBufferLZ4(buffer,outname,readIn);
            zWrite(buffer,outname,readIn);
         }else {
            writeBuffer(buffer,outname,readIn);
         }
         delete [] buffer; buffer=NULL;
      }
      if(this->key.size()>1){
         std::cerr<<"One/Some variable(s) not extracted. Check your input against variables in this file:"<<std::endl;
         for (auto k:this->key){
            std::cerr<<"\t"<<k<<std::endl;
         }
         
      }
   }



   bool writeBuffer(const char* buffer,std::string outFileName,uint64_t amount){
         //Write 
         std::fstream outfile;
         outfile = std::fstream(outFileName, std::ios::out | std::ios::binary);
         outfile.write(buffer,amount);
         outfile.close();
      return true;
   }

   bool zWrite(const char* data,std::string outFileName,uint64_t buffSize){
            outFileName+=".zst";
            std::fstream outfile;
            outfile = std::fstream(outFileName, std::ios::out | std::ios::binary);
            size_t bound=ZSTD_compressBound(buffSize);
            char *buffer_z = new char[bound];
            size_t res= ZSTD_compress(buffer_z,bound,data,buffSize,5); 
            outfile.write(buffer_z,res);
            outfile.close();
            delete [] buffer_z;
    
      return true;
   }


   void zDecompress(){

      std::ifstream infile( this->filename, std::ifstream::binary );
      infile.seekg (0,infile.end);
      uint64_t buffSize_z = infile.tellg();
      infile.clear();
      infile.seekg(0,std::ios::beg);
      char* data=new char[buffSize_z];
      char* buffer=new char[2*buffSize_z];
      infile.read(data,buffSize_z);
      infile.close();

      //Outfile name
      size_t lastindex = this->filename.find_last_of("."); 
      std:: string outname = this->filename.substr(0, lastindex); 
      std::fstream outfile;
      outfile = std::fstream(outname, std::ios::out | std::ios::binary);
      size_t res=ZSTD_decompress(buffer,2*buffSize_z,data,buffSize_z);
      std::cerr<<ZSTD_getErrorName(res)<<std::endl;
      outfile.write(buffer,res);
      outfile.close();
      delete [] data;
      delete [] buffer;
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
          pos+= kind.size(); 
      }
  
      pos=0;
      while ((pos = xmlTag.find(complement, pos)) != std::string::npos) {
          tags.at(count).second = pos ;
          pos+= complement.size(); 
          count++;
      }
      return tags;
  }

   std::string getSizePretty(uint64_t bytes){

      float thissize=(double)bytes;
      float kilo = 1024.0;
      float mega = (1024.0*1024.0);
      float giga = 1024.0*1024.0*1024.0;
      
      
      float portion=100*(thissize/this->filesize);
      std::string retval;
      if (thissize/giga > 1.0){
         retval= std::to_string(thissize/giga)+"GB"+" [" +std::to_string(portion)+ "%]";
      }else if(thissize/mega > 1.0){
         retval= std::to_string(thissize/mega)+"MB"+" [" +std::to_string(portion)+ "%]";
      }else if(thissize/kilo > 1.0){
         retval= std::to_string(thissize/kilo)+"KB"+" [" +std::to_string(portion)+ "%]";
      }else{
         retval= std::to_string(thissize)+"B" + " [" +std::to_string(portion)+ "%]";
      }

      return retval;

   }
   uint64_t get_file_size(std::string filename) // path to file
   {
      FILE *p_file = NULL;
      p_file = fopen(filename.c_str(),"rb");
      fseek(p_file,0,SEEK_END);
      uint64_t size = ftell(p_file);
      fclose(p_file);
      p_file=NULL;
      return size;
   }

   void extractTags(){

      int cnt=0;
      for (auto c:tags){
         std::string line=xmlTag.substr(c.first,c.second-c.first);
         
         std::size_t count = 0, pos = 0 ;
         std::string tab(" ");
         std::string equals("=");
         std::string memIdr(">");
         
         //Get memory location
         std::size_t memloc = line.rfind(memIdr);
         std::string mem ;
         if (memloc !=std::string::npos){
            mem=line.substr(memloc+memIdr.size(),line.size()-memloc);
         }
         
         while ((pos = line.find(equals, 0)) != std::string::npos) {
            uint64_t p0=0,p1=0,p2=0;
            p0=line.find(tab,p0);
            p1=line.find(equals,p0);
            p2=line.find(tab,p1);
            if(p2==std::string::npos){
               p2=line.find(">",p1);
            }
            bool test=line[p2-1]=='"';
            if (!test){
               p2=line.find(tab,p2+1);
            }
            

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
         //Get location in binary file
         this->vars.at(cnt).attributes.insert(std::make_pair("tag_location",mem));
         //Calculate datasize
         std::map<std::string,std::string>::const_iterator it_0;
         std::map<std::string,std::string>::const_iterator it_1;
         std::map<std::string,std::string>::const_iterator it_2;
         it_0 = this->vars.at(cnt).attributes.find("arraysize"); 
         it_1 = this->vars.at(cnt).attributes.find("datasize"); 
         it_2 = this->vars.at(cnt).attributes.find("vectorsize");
         bool doCalculate= ( it_0!=this->vars.at(cnt).attributes.end() &&
                            it_1!=this->vars.at(cnt).attributes.end() &&
                            it_2!=this->vars.at(cnt).attributes.end());

         if (doCalculate){
            uint64_t arraysize;
            uint64_t datasize;
            uint64_t vectorsize;
            arraysize=std::stoul(it_0->second,nullptr,10);
            datasize=std::stoul(it_1->second,nullptr,10);
            vectorsize=std::stoul(it_2->second,nullptr,10);
            uint64_t dSize=arraysize*vectorsize*datasize;
            std::string humanSize=getSizePretty(dSize);
            this->vars.at(cnt).attributes.insert(std::make_pair("full_size",humanSize));
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
            it = c.attributes.find("tag_location"); 
            if (it != c.attributes.end()){
               std::cout<<"  "<<it->first<<"="<<it->second<<std::endl;
               };
            it = c.attributes.find("full_size"); 
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
   std::cout << "Usage: ./vlsvdump [options] <file>  [variablelist] " << std::endl;
   std::cout << "Options:" << std::endl;
   std::cout << "\t -h \t\t Header information. No data." << std::endl;
   std::cout << "\t -hv \t\t Header information with a bit of spice. No data." << std::endl;
   std::cout << "\t -e \t\t Extracts variable(s) to binary file(s). Variable(s) should be delimited by comma"<<std::endl;
   std::cout << "\t -ez \t\t Extracts variable(s) and compresses using zstd to binary file(s). Variable(s) should be delimited by comma"<<std::endl;
   std::cout << "\t -d  \t\t Decompresses file compressed by vlsvdump using zstd"<<std::endl;
}


int main(int argc, char *argv[]){
   if (argc<3){
      std::cerr<<"Invalid Usage\n";
      helpMenu();
      return 1;
   }
   
   std::string fname(argv[2]);


   File vlsv(fname,argc,argv);
   return 0;

}
