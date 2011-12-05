/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include "profile.h"
//include craypat api headers if compiled with craypat on Cray XT/XE
#ifdef CRAYPAT
#include "pat_api.h"
#endif

#ifndef MPILOGGER
  #define output cout
#else
  #include "mpilogger.h"
  extern MPILogger mpilogger;
  #define output mpilogger
#endif


#include "profile.h"

using namespace std;

#ifdef PROFILE

namespace profile 
{
   namespace 
   {
      struct TimerData{
	 string label;          //print label 
	 string workUnitLabel;   //unit for the counter workUnitCount
	 double workUnits;        // how many units of work have we done. If -1 it is not counted, or printed
	 vector<string> groups; // What user-defined groups does this timer belong to, e.g., "MPI", "IO", etc..

	 int id; // unique id identifying this timer (index for timers)
	 int parentId;  //key of parent (id)
	 vector<int> childIds; //children of this timer

	 int level;  //what hierarchy level
	 int count; //how many times have this been accumulated
	 double time; // total time accumulated
	 double startTime; //Starting time of previous start() call
      };
      
      //vector with timers
      vector<TimerData> timers;
      
      //current position in timer hierarchy 
      int currentId=-1;

      //is this profiler initialized
      bool initialized=false;
      
      //used with MPI reduction operator
      struct doubleRankPair {
	 double val;
	 int rank;
      };

      //defines print-area widths for print() output
      const int indentWidth=2; //how many spaces each level is indented
      const int floatWidth=10; //width of float fields;
      const int intWidth=6;   //width of int fields;
      const int unitWidth=4;  //width of workunit label

      
      //djb2 hash function copied from
      //http://www.cse.yorku.ca/~oz/hash.html
      unsigned long hash(const char *str)
      {
	 unsigned long hash = 5381;
	 int c;
	 while ( (c = *str++) )
	       hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
	 return hash;
      }

      
      //Hash value identifying all labels, groups and workunitlabels.
      //If any strings differ, hash should differ.
      int getTimersHash(){
	 string allLabels;
	 for (vector<TimerData>::const_iterator timer = timers.begin();
	       timer != timers.end(); ++timer ) {
	    allLabels.append(timer->label);
	    allLabels.append(timer->workUnitLabel);
	    for (vector<string>::const_iterator group = timer->groups.begin();
		  group != timer->groups.end(); ++group ) {
	       allLabels.append(*group);
	    }
	 }
	 return (int)hash(allLabels.c_str());
      }
   
      
//construct a new timer for the current level.
      int constructTimer(const string &label,int parentId,const vector<string> groups){
	 TimerData timerData;
	 timerData.label=label;
	 timerData.groups=groups;
	 timerData.id=timers.size(); //will be added last to timers vector
	 timerData.parentId=parentId;
	 
	 timerData.workUnits=-1;
	 timerData.workUnitLabel="";
	 
	 if(parentId!=-1) 
	    timerData.level=timers[parentId].level+1;
	 else //this is the special case when one adds the root timer
	    timerData.level=0;
	 timerData.time=0;
	 timerData.startTime=-1;
	 timerData.count=0;
	 //timerData.work  UnitCount initialized in stop
	 //add timer   
	 timers.push_back(timerData);
	 //add timer to tree
	 if(parentId!=-1)
	    timers[parentId].childIds.push_back(timerData.id);
	 return timerData.id;
      }
      
      //initialize profiler, called by first start/initializeTimer. This adds the root timer
      bool init(){
	 if(!initialized){
	    initialized=true;
	    vector<string> group;
	    group.push_back("Total");
	    //no timer yet        
	    currentId=-1;
	    //mainId will be 0, parent is -1 (does not exist)
	    int id=constructTimer("total",-1,group);
	    //start root timer, is stopped in print.      
	    start(id);
	 }
	 return true;
      }
      
      //this function returns the time in seconds 
      double getTime() {
	 return MPI_Wtime();
      }
      //this function returns the accuracy of the timer     
      double getTick() {
	 return MPI_Wtick();
      }
      
      bool computeStatistics(double value,
			      double *average,double *total,
			      doubleRankPair *max,doubleRankPair *min,
			      int rootRank,MPI_Comm comm){
	 int rank,nProcesses;
	 doubleRankPair in;
	 double tot;
	 
	 MPI_Comm_rank(comm,&rank);
	 MPI_Comm_size(comm,&nProcesses);
	 
	 if(average!=NULL || total !=NULL)
	    MPI_Reduce(&value,&tot,1,MPI_DOUBLE,MPI_SUM,rootRank,comm);

	 if(average!=NULL)
	    *average=tot/nProcesses;
	 if(total!=NULL)
	    *total=tot;

	 
	 in.val=value;
	 in.rank=rank;
	 if(max!=NULL)
	    MPI_Reduce(&in,max,1,MPI_DOUBLE_INT,MPI_MAXLOC,rootRank,comm);            
	 if(min!=NULL)
	    MPI_Reduce(&in,min,1,MPI_DOUBLE_INT,MPI_MINLOC,rootRank,comm);
	 return true;
      }


      bool getGroupInfo(string group,int id, double &time,MPI_Comm comm){
	 bool hasGroup=false;
	 for (vector<string>::const_iterator g = timers[id].groups.begin();
	       g != timers[id].groups.end(); ++g ) {
	    if(group==*g){
	       time=time+timers[id].time;
	       hasGroup=true;
	       break;
	    }
	 }
      
      //recursively collect time data, do not collect if this timer already is ib group, avoid double counting
	 if(!hasGroup) { 
	    for(unsigned int i=0;i<timers[id].childIds.size();i++){
	       getGroupInfo(group,timers[id].childIds[i],time,comm);
	    }
	 }
	 return true;
      }
      
      bool printGroup(string group,const string groupId,double time,double minParentFraction,size_t labelWidth,size_t groupWidth,MPI_Comm comm){
	 int rank,nProcesses;
	 doubleRankPair max;
	 doubleRankPair min;
	 double averageTime,totalTime,totalParentTime;
	 double localGroupTime=time;
	 
	 MPI_Comm_rank(comm,&rank);
	 MPI_Comm_size(comm,&nProcesses);

	 computeStatistics(localGroupTime,
			   &averageTime,&totalTime,&max,&min,0,comm);
	 computeStatistics(timers[0].time,
			   NULL,&totalParentTime,NULL,NULL,0,comm);
	 if(rank==0 && minParentFraction<=totalTime/totalParentTime){
	    output << setw(groupWidth+1) << groupId;
	    output << setw(labelWidth+1) << group;
	    output << setw(floatWidth) << averageTime;
	    output << setw(floatWidth) << 100.0*totalTime/totalParentTime;
	    output << setw(floatWidth) << max.val;
	    output << setw(intWidth)   << max.rank;
	    output << setw(floatWidth) << min.val;
	    output << setw(intWidth)   << min.rank;
	    output << endl;
	 }
	 return true;
      }
      
      //print groups
      bool printGroups(double minParentFraction,size_t labelWidth,size_t groupWidth,int totalWidth,MPI_Comm comm,const map<string, string> &groupIds){
	 //construct map from groups to labels in group
	 int rank;
	 map<string,vector<int> > groups;
	 MPI_Comm_rank(comm,&rank);
	 for(unsigned int id=0;id<timers.size();id++){
	    for (vector<string>::const_iterator group = timers[id].groups.begin();
		  group != timers[id].groups.end(); ++group ) {
	       groups[*group].push_back(id);
	    }
	 }
	 if(rank==0){
	    for(int i=0;i<totalWidth/2 -4;i++) output <<"-";
	    output <<" Groups ";
	    for(int i=0;i<totalWidth/2 -3;i++) output <<"-";
	    output<<endl;
	 }


	 
	 for(map<string, vector<int> >::iterator group=groups.begin();
	    group!=groups.end();++group){
	    double time=0.0;
	    getGroupInfo(group->first,0,time,comm);
	    printGroup(group->first, groupIds.at(group->first), time, minParentFraction, labelWidth,groupWidth,comm);
	 }
	 
	 return true;
      }

      //print a timer, call children recursively
      bool printTimer(int id,double minParentFraction,size_t labelWidth,size_t groupWidth,int totalWidth,MPI_Comm comm,const map<string,string> &groupIds){
	 int rank,nProcesses;
	 doubleRankPair max;
	 doubleRankPair min;
	 double averageTime,averageCount,totalTime,totalParentTime;
	 double hasWorkUnits;
	 double totalUnits;
	 double averageUnits;
	 int isPrinted=1;
	 MPI_Comm_rank(comm,&rank);
	 MPI_Comm_size(comm,&nProcesses);

	 //do not write out root (id=0)
	 if(id!=0){
	    //coumpute time statistics      
	    computeStatistics(timers[id].time,
			      &averageTime,&totalTime,&max,&min,0,comm);
	    computeStatistics(timers[timers[id].parentId].time,
			      NULL,&totalParentTime,NULL,NULL,0,comm);

	    if(totalTime/totalParentTime<minParentFraction)
	       isPrinted=0;
	    MPI_Bcast(&isPrinted,1,MPI_INT,0,comm);
	    
	    //label & time statistics
	    if(isPrinted==1 && rank==0){
	       bool hasNoGroups = true;
	       int indent=(timers[id].level-1)*indentWidth;
	       for (vector<string>::const_iterator group = timers[id].groups.begin();
		    group != timers[id].groups.end(); ++group) {
		  output << setw(1) << groupIds.at(*group);
	          hasNoGroups = false;
	       }
	       if(hasNoGroups) output << setw(groupWidth+1) << "";
	       else output << setw(groupWidth-timers[id].groups.size()+1) << "";
	       output << setw(indent) << "";
	       output << setw(labelWidth+1-indent) << setiosflags(ios::left) << timers[id].label;
	       
	       output << setw(floatWidth) << averageTime;
	       output << setw(floatWidth) << 100.0*totalTime/totalParentTime;
	       output << setw(floatWidth) << max.val;
	       output << setw(intWidth)   << max.rank;
	       output << setw(floatWidth) << min.val;
	       output << setw(intWidth)   << min.rank;
	    }

	    //count statistics
	    computeStatistics((double)timers[id].count,
			      &averageCount,NULL,NULL,NULL,0,comm);
	    if(isPrinted==1 && rank==0){
	       output << setw(floatWidth) << averageCount;
	    }
	    
	    //workunits statistics
	    if (timers[id].workUnits<0)
	       hasWorkUnits=0.0;
	    else
	       hasWorkUnits=1.0;             
	    //collect in min data if everybody has workunits
	    computeStatistics(hasWorkUnits,
			      NULL,NULL,NULL,&min,0,comm);
	    computeStatistics(timers[id].workUnits,
			      &averageUnits,&totalUnits,NULL,NULL,0,comm);
	    //print if rank is zero, and units defined for all processes
	    //note how the averages are computed. This is to avoid one process with little data to     
	    //skew results one way or the other     
	    if(isPrinted==1 && rank==0 && min.val >0.5 ){
	       output << setw(floatWidth) << totalUnits/averageTime;
	       output << setw(floatWidth) << averageUnits/averageTime;
	       output << timers[id].workUnitLabel<<"/s";
	    }

	    //new line
	    if(isPrinted==1 && rank==0){
	       output<<endl;
	    }
	 }
	 
	 //recursively print child labels
	 if(isPrinted==1) { //only write out children, it the label itself is printed
	    for(unsigned int i=0;i<timers[id].childIds.size();i++){
	       printTimer(timers[id].childIds[i],minParentFraction,labelWidth,groupWidth,totalWidth,comm,groupIds);
	    }
	 }
	 return true;
      }
      
      bool printFooter(int totalWidth,MPI_Comm comm){
	 int rank;
	 MPI_Comm_rank(comm,&rank);
	 if(rank==0){
	    for(int i=0;i<totalWidth;i++) output <<"-";
	    output<<endl;
	 }
	 return true;
      }
      
      bool printHeader(double minParentFraction,size_t labelWidth,size_t groupWidth,int totalWidth,MPI_Comm comm){
	 int rank;
	 MPI_Comm_rank(comm,&rank);
	 if(rank==0){
	    for(int i=0;i<totalWidth;i++) output <<"-";
	    output<<endl;
	    output << "Profiler results with parent % larger than " << minParentFraction*100.0;
	    output<<endl;
	    for(int i=0;i<totalWidth;i++) output <<"-";
	    output<<endl;
	    output<<setw(groupWidth+1+labelWidth+1)<< setiosflags(ios::left) << "";
	    output<<setw(4*floatWidth+2*intWidth) <<"Time(s)";
	    output<<setw(floatWidth)<<"Calls";
	    output<<setw(2*floatWidth)<<"Workunit-rate";
	    output<<endl;
	    
	    output<<setw(groupWidth+1)<< "Groups";
	    output<<setw(labelWidth+1)<< "Label";
	    //  time
	    output<<setw(floatWidth) <<"Average";
	    output<<setw(floatWidth) <<"parent %";
	    output<<setw(floatWidth) <<"Maximum";
	    output<<setw(intWidth) << "Rank";
	    output<<setw(floatWidth)<< "Minimum";
	    output<<setw(intWidth) << "Rank";
	    //call count
	    output<<setw(floatWidth) << "Average";
	    // workunit rate    
	    output<<setw(floatWidth) << "Average";
	    output<<setw(floatWidth) << "Per process";
	    output<<endl;
	 }
	 return true;
      }
//print out all timers
      bool printTimers(double minParentFraction,size_t labelWidth,size_t groupWidth,int totalWidth,MPI_Comm comm,const map<string, string> &groupIds){
	 int rank;
	 MPI_Comm_rank(comm,&rank);
	 if(rank==0){
	    for(int i=0;i<totalWidth/2-5;i++) output <<"-";
	    output << " Profile ";
	    for(int i=0;i<totalWidth/2-5;i++) output <<"-";
	    output<<endl;
	 }
	 //recursively print all timers
	 printTimer(0,minParentFraction,labelWidth,groupWidth,totalWidth,comm,groupIds);
	 return true;
      }
   }

// end unnamed namespace
//----------------------------------------------------------------------------------------------
// public functions begin    

   //get id number of a timer, return -1 if it does not exist
   int getId(const string &label){
      //find child with this id
      int childId=-1;
      for(unsigned int i=0;i<timers[currentId].childIds.size();i++)
	 if (timers[timers[currentId].childIds[i]].label==label){
	    childId=timers[currentId].childIds[i];
	    break;
	 }
      return childId;
   }

      
   //initialize a timer, with a particular label belonging to some groups
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const string &label,const vector<string> &groups){
      //check if the global profiler is initialized
      if(!initialized)
	 init();
      int id=getId(label); //check if it exists
      if(id>=0)
	 //do nothing if it exists      
	 return id; 
      else
	 //create new timer if it did not exist
	 return constructTimer(label,currentId,groups);
   }


   //initialize a timer, with a particular label belonging to a group
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const string &label,const string &group){
      //check if the global profiler is initialized
      vector<string> groups;
      groups.push_back(group);
      return initializeTimer(label,groups);
   }

   //initialize a timer, with a particular label belonging to two groups
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const string &label,
		     const string &group1,
		     const string &group2
		     ){
      //check if the global profiler is initialized
      vector<string> groups;
      groups.push_back(group1);
      groups.push_back(group2);
      return initializeTimer(label,groups);
   }

   //initialize a timer, with a particular label belonging to three groups
   //returns id of new  timer. If timer exists, then that id is returned.
   int initializeTimer(const string &label,
		     const string &group1,
		     const string &group2,
		     const string &group3
		     ){
      //check if the global profiler is initialized
      vector<string> groups;
      groups.push_back(group1);
      groups.push_back(group2);
      groups.push_back(group3);
      return initializeTimer(label,groups);
   }


   //initialize a timer, with a particular label belonging to no group
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const string &label){
      //check if the global profiler is initialized
      vector<string> groups; //empty vector
      return initializeTimer(label,groups);
   }
   
   //start timer, with id
   bool start(int id){
      if(currentId!=timers[id].parentId){
	 output << "Starting timer that is not a child of the current profiling region" <<endl;
	 return false;
      }
      currentId=id;      
      //start tuner
      timers[currentId].startTime=getTime();      
#ifdef CRAYPAT
      PAT_region_begin(currentId+1,timers[id].label.c_str());
#endif
      return true;        
   }

   //start timer, with label
   bool start(const string &label){
      //If the timer exists, then initializeTimer just returns its id, otherwise it is constructed.
      //Make the timer the current one
      currentId=initializeTimer(label);
      //start timer
      timers[currentId].startTime=getTime();
#ifdef CRAYPAT
      PAT_region_begin(currentId+1,label.c_str());
#endif
      return true;        
   }
   
   //stop with workunits
   bool stop (const string &label,
	    const double workUnits,
	    const string &workUnitLabel){
      if(label != timers[currentId].label ){
	 output << "ERROR: label missmatch in profile::stop Stopping "<< label <<" at level " << timers[currentId].level << endl;
	 return false;
      }
      stop(currentId,workUnits,workUnitLabel);
      return true;
   }

   //stop with workunits        
   bool stop (int id,
	    double workUnits,
	    const string &workUnitLabel){
      double stopTime=getTime();
      if(id != currentId ){
	 output << "ERROR: id missmatch in profile::stop Stopping "<< id <<" at level " << timers[currentId].level << endl;
	 return false;
      }
      
#ifdef CRAYPAT
      PAT_region_end(currentId+1);
#endif  
      
      if(timers[currentId].count!=0){
	 //if this, or a previous, stop did not include work units then do not add them
	 //work units have to be defined for all stops with a certain (full)label
	 if(workUnits<0 || timers[currentId].workUnits<0){
	    timers[currentId].workUnits=-1;
	 }
	 else{
	    timers[currentId].workUnits+=workUnits;
	 }
      }
      else{
	 //firsttime, initialize workUnit stuff here
	 if(workUnits>=0.0 ){
	    //we have workUnits for this counter
	    timers[currentId].workUnits=workUnits;
	    timers[currentId].workUnitLabel=workUnitLabel;
	 }
	 else{
	    //no workUnits for this counter
	    timers[currentId].workUnits=-1.0;
	 }
      }
      
      timers[currentId].time+=(stopTime-timers[currentId].startTime);
      timers[currentId].count++;
      
      //go down in hierarchy    
      currentId=timers[currentId].parentId;
      return true;
   }

   //print out global timers
   // If any labels differ, then the print cannot proceed. It has to be consistent for all processes in communicator
   bool print(MPI_Comm comm,double minParentFraction){
      int rank,nProcesses;
      int timersHash;
      vector<int> allTimersHash;
      size_t labelWidth=0;    //width of column with timer labels
      size_t groupWidth=6;    //width of column with group codes ("Groups" is 6 letters)
      int totalWidth;       //total width of the table

      //stop main timer, will not be printed, but is used to get total time and to get parent % for level=1 processes
      stop(0);
      
      MPI_Comm_rank(comm,&rank);
      MPI_Comm_size(comm,&nProcesses);

      //Synchronize here. This is to make sure that the MPI_Allgather does not interfere with any other ongoing
      //MPI communication (due to load imbalance).
      
      MPI_Barrier(comm);
      //compute labelWidth
      for (vector<TimerData>::iterator timer=timers.begin();
	    timer!=timers.end();++timer){
	    size_t width=timer->label.length()+(timer->level-1)*indentWidth;
	    labelWidth=max(labelWidth,width);
      }

      
      //check that all processes have identical timers to avoid deadlocks
      timersHash=getTimersHash();
      allTimersHash.resize(nProcesses);
      MPI_Allgather(&timersHash,1,MPI_INT,&(allTimersHash[0]),1,MPI_INT,comm);
      for(int i=0;i<nProcesses;i++){
	    if(timersHash!=allTimersHash[i]){
	       if(rank==0){
		  output <<"Error in profile::print, labels, groups or workunits on different processes do not match" <<endl;
	       }
	       return false;
	    }
      }

      //make sure we use default floats, and not fixed or other format
      output <<resetiosflags( ios::floatfield );
      //set float precision
      output <<setprecision(floatWidth-6); //6 is needed for ".", "e+xx" and a space
      
      // Creating map of group names to group one-letter ID
      map<string, string> groupIds;
      for(unsigned int id=0;id<timers.size();id++) {
	 size_t width = timers[id].groups.size();
	 if(width != 0) groupWidth = max(width, groupWidth);
	 for (vector<string>::const_iterator group = timers[id].groups.begin();
	      group != timers[id].groups.end(); ++group ) {
	    groupIds[*group] = *group;
	 }
      }
      int character = 'A'; // ASCII A, TODO skip after Z to a, see ascii table
      for(map<string, string>::const_iterator group = groupIds.begin();
	  group != groupIds.end(); ++group) {
	 groupIds[group->first] = character++;
      }
      
      totalWidth=groupWidth+1+labelWidth+1+floatWidth*7+intWidth*2+unitWidth;
      
      //print header 
      printHeader(minParentFraction,labelWidth,groupWidth,totalWidth,comm);
      //print out all labels recursively
      printTimers(minParentFraction,labelWidth,groupWidth,totalWidth,comm,groupIds);
      //print groups
      printGroups(minParentFraction,labelWidth,groupWidth,totalWidth,comm,groupIds);
      //print footer  
      printFooter(totalWidth,comm);
      //start root timer again in case we continue and call print several times
      start(0);
      return true;
   }
}

#else
namespace profile 
{
   bool start(int id){return true;}
   bool stop (int id,double workUnits,
	    const string &workUnitLabel){return true;}
   bool start(const string &label){return true;}
   bool stop (const string &label,double workUnits,
	    const string &workUnitLabel){return true;}
   bool print(MPI_Comm comm,double minParentPercent){return true;}
   int getId(const string &label) {return 0;}

   int initializeTimer(const string &label,const vector<string> &groups) { return 0;}
   int initializeTimer(const string &label,const string &group){return 0;}
   int initializeTimer(const string &label){return 0;}
}
#endif

#undef output
