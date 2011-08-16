#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <map>

//include craypat api headers if compiled with craypat on Cray XT/XE
#ifdef CRAYPAT
#include "pat_api.h"
#endif

#ifdef NO_MPILOGGER
    #define mpilogger cerr
    #define write flush
#else
    #include "mpilogger.h"
    extern MPILogger mpilogger;
#endif

using namespace std;

#ifdef PROFILE

namespace profile 
{
    namespace 
    {
        struct TimerData{
            string parent; //key of parent (fullLabel)
            string label; //print label 
            string workUnitLabel; //unit for the counter workUnitCount
            double time; // total time accumulated
            int level;  //what hierarchy level
            int count; //how many times have this been accumulated
            int index; // unique index identifying this timer (usefule for Craypat)
            double workUnits; // how many units of work have we done. If -1 it is not counted, or printed
        };
        //store all timers in map, key is fullLabel name
        std::map<string,TimerData> timers;

        //current position in timer hierarchy & label & startTime for the different levels
        int currentLevel=-1;
        std::vector<string> labels;
        std::vector<double> startTime;
        double overHeadTime=0;
        unsigned int startStopCalls=0;
        struct doubleRankPair {
            double val;
            int rank;
        };


        string constructFullLabel(int maxLevel){
            string label;
            if(maxLevel<0)
                label.append("/");
            for(int i=0;i<=maxLevel;i++){
                label.append("/");
                label.append(labels[i]);
            }
            return label;
        }
    }
    
    bool start(const string &label){
        double t1=MPI_Wtime();
        currentLevel++;
        //resize vectors if needed
        if(int(labels.size())<=currentLevel)
            labels.resize(currentLevel+5);
        if(int(startTime.size())<=currentLevel)
            startTime.resize(currentLevel+5);
        
        labels[currentLevel]=label;

        //if fullname timer does not yet exist, it is constructed using the default constructor
        //it could be constructed also in stop, but here we need the index for craypat so we do it in start
        string fullName=constructFullLabel(currentLevel);
        if(timers.find(fullName) == timers.end()){
            //does not exist, add it
            TimerData timerData;
            timerData.level=currentLevel;
            timerData.time=0;
            timerData.count=0;
            //timerData.workUnitCount initialized in stop
            timerData.parent=constructFullLabel(currentLevel-1);
            timerData.label=label;
            timerData.index=timers.size();  //index will be consecutive starting from 0
            timers[fullName]=timerData;
        }
        

#ifdef CRAYPAT
        PAT_region_begin(timers[fullName].index,label.c_str());
#endif
        startTime[currentLevel]=MPI_Wtime();

        //collect information for overHEad statistics
        startStopCalls++;
        overHeadTime+=MPI_Wtime()-t1;
        return true;        
    }

    

    //stop with workunits
    bool stop (const std::string &label,double workUnits,
               const std::string &workUnitLabel){
        double t1=MPI_Wtime();
        double stopTime=MPI_Wtime();
        string fullName=constructFullLabel(currentLevel);
#ifdef CRAYPAT
        PAT_region_end(timers[fullName].index);
#endif
        if(currentLevel<=-1){
            mpilogger << "ERROR: nothing to stop in profile::stop Stopping "<<
                label <<endl;
            return false;
        }
        if(labels[currentLevel] != label){
            mpilogger << "ERROR: label missmatch in profile::stop Stopping "<<
                label <<" at level " <<fullName << endl;
            return false;
        }
        
        //firsttime, initialize workUnit stuff here
        if(timers[fullName].count==0){
            if(workUnits>=0.0 ){
                //we have workUnits for this counter
                timers[fullName].workUnits=workUnits;
                timers[fullName].workUnitLabel=workUnitLabel;
            }
            else{
                //no workUnits for this counter
                timers[fullName].workUnits=-1.0;
            }
        }
        else {
            //if this, or a previous, stop did not include work units then do not add them
            //work units have to be defined for all stops with a certain (full)label
            
            if(workUnits<0 || timers[fullName].workUnits<0){
                timers[fullName].workUnits=-1;
            }
            else{
                timers[fullName].workUnits+=workUnits;
            }
        }

        timers[fullName].time+=(stopTime-startTime[currentLevel]);
        timers[fullName].count++;
        currentLevel--;

        //collect information for overHead computation
        startStopCalls++;
        overHeadTime+=MPI_Wtime()-t1;
        return true;
    }
    
    
    //print out global timers
    //This will assume that all processes have identical labels, will fail otherwise (FIXME)
    bool print(MPI_Comm comm){
        int rank,nProcesses;
        int nTimers;
        vector<int> nAllTimers;

        const int indentWidth=2; //how many spaces each level is indented
        const int floatWidth=12; //width of float fields;
        const int intWidth=6;   //width of int fields;
        const int unitWidth=4;  //width of workunit label
        size_t labelWidth=0;    //width of first column with timer labels
        int totalWidth;       //total width of the table
        
        //compute labelWidth
        for (std::map<string,TimerData>::iterator timer=timers.begin();
             timer!=timers.end();++timer){
            size_t width=timer->second.label.length()+timer->second.level*indentWidth;
            labelWidth=max(labelWidth,width);
        }
        totalWidth=labelWidth+1+floatWidth*7+intWidth*2+unitWidth;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&nProcesses);
        nAllTimers.resize(nProcesses);

        nTimers=timers.size();
        MPI_Allgather(&nTimers,1,MPI_INT,&(nAllTimers[0]),1,MPI_INT,comm);

        //check that all processes have identical number of timers to avoid deadlocks
        //does not ensure the labels themselves are the same (FIXME with e.g. with hash of labels)
        for(int i=0;i<nProcesses;i++){
            if(nTimers!=nAllTimers[i]){
                if(rank==0){
                    mpilogger <<"Error in profile::print, different number of labels on different processes" <<endl;
                }
                return false;
            }
        }

        //make sure we use default floats, and not fixed or other format
        mpilogger <<resetiosflags( ios::floatfield );
        //set float precision
        mpilogger <<setprecision(floatWidth-6); //6 is needed for ".", "e+xx" and a space
        //print out header
        if(rank==0){
            for(int i=0;i<totalWidth/2-5;i++)mpilogger <<"-";
            mpilogger << " Profile ";
            for(int i=0;i<totalWidth/2-5;i++)mpilogger <<"-";
            mpilogger<<endl;

            mpilogger<<setw(labelWidth+1)<< setiosflags(ios::left) << "";
            mpilogger<<setw(4*floatWidth+2*intWidth) <<"Time(s)";
            mpilogger<<setw(floatWidth)<<"Call-count";
            mpilogger<<setw(2*floatWidth)<<"Workunit-rate";
            mpilogger<<endl;
            
            mpilogger<<setw(labelWidth+1)<< "Label";
            //time
            mpilogger<<setw(floatWidth) <<"Average";
            mpilogger<<setw(floatWidth) <<"% of parent";
            mpilogger<<setw(floatWidth) <<"Maximum";
            mpilogger<<setw(intWidth) << "Rank";
            mpilogger<<setw(floatWidth)<< "Minimum";
            mpilogger<<setw(intWidth) << "Rank";
            //call count
            mpilogger<<setw(floatWidth) << "Average";
            // workunit rate
            mpilogger<<setw(floatWidth) << "Total";
            mpilogger<<setw(floatWidth) << "Average";
            mpilogger<<endl;
            for(int i=0;i<totalWidth;i++) mpilogger <<"-";
            mpilogger<<endl;            
        }

        //sort labels according to index (corresponds to first creation time)
        std::vector<string> listOrder(timers.size(),"");
        for (std::map<string,TimerData>::iterator timer=timers.begin();
             timer!=timers.end();++timer){
            listOrder[timer->second.index]=timer->first;
        }

        //loop over listOrder so that timers are printed in order of creation
        for(unsigned int i=0;i<listOrder.size();i++){
            double sum,parentSum;
            double aveTime;
            int intSum;
            doubleRankPair in;
            doubleRankPair out;
            //get timer that is now to be computed
            std::map<string,TimerData>::iterator timer=timers.find(listOrder[i]);

            //first compute parent sum of times
            if(timer->second.level>0){
                in.val=timers[timer->second.parent].time;
                in.rank=rank;
                MPI_Reduce(&(in.val),&parentSum,1,MPI_DOUBLE,MPI_SUM,0,comm);
            }
            
                
            //then current timer sums
            in.val=timer->second.time;
            in.rank=rank;
            MPI_Reduce(&(in.val),&sum,1,MPI_DOUBLE,MPI_SUM,0,comm);
            
            if(rank==0){
                int indent=timer->second.level*indentWidth;
                mpilogger<<setw(indent)<<"";
                mpilogger<<setw(labelWidth+1-indent)<< setiosflags(ios::left)
                    << timer->second.label;

                aveTime= sum/nProcesses;
                mpilogger << setw(floatWidth) << aveTime;
                if(timer->second.level>0)
                    mpilogger << setw(floatWidth) << 100.0*sum/parentSum;
                else
                    mpilogger << setw(floatWidth) << " ";
                
            }
            MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,0,comm);            
            if(rank==0){
                mpilogger << setw(floatWidth) << out.val;
                mpilogger << setw(intWidth) <<out.rank;
            }
            
            MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,0,comm);
            if(rank==0){
                mpilogger << setw(floatWidth) << out.val;
                mpilogger << setw(intWidth)<<out.rank;
            }
            
            //current count statistics
            MPI_Reduce(&( timer->second.count),&intSum,1,MPI_INT,MPI_SUM,0,comm);
            if(rank==0){
                mpilogger << setw(floatWidth) << (double)intSum/nProcesses;
            }
             
            //workunits/time statistics
            double workRate[2];
            double sums[2];

            workRate[0]=timer->second.workUnits/timer->second.time;
            if (workRate[0]<0)
                workRate[1]=1.0;//use this to check if for any process the units were not defined
            else
                workRate[1]=0.0;
            
            MPI_Reduce(workRate,sums,2,MPI_DOUBLE,MPI_SUM,0,comm);
            //print if rank is zero, and units defined for all processes
            if(rank==0 && sums[1] < 0.5){
                mpilogger << setw(floatWidth) << sums[0];
                mpilogger << setw(floatWidth) << sums[0]/nProcesses;
                mpilogger << timer->second.workUnitLabel<<"/s";
            }
            if(rank==0){
                mpilogger<<endl;
            }
        }
        //Overhead statistics
        double sum;
        double aveOverhead=overHeadTime/startStopCalls;
        MPI_Reduce(&aveOverhead,&sum,1,MPI_DOUBLE,MPI_SUM,0,comm);
        if(rank==0){
            for(int i=0;i<totalWidth;i++) mpilogger <<"-";
            mpilogger<<endl;
            mpilogger << setw(labelWidth+1) << "Overhead per call";
            mpilogger << setw(floatWidth) << sum/nProcesses;
            mpilogger<<endl;
        }
        MPI_Reduce(&overHeadTime,&sum,1,MPI_DOUBLE,MPI_SUM,0,comm);
        if(rank==0){
            mpilogger << setw(labelWidth+1) << "Total profiling overhead";
            mpilogger << setw(floatWidth) << sum/nProcesses;
            mpilogger << endl;
        }
        
        //footer line
        if(rank==0){
            for(int i=0;i<totalWidth;i++) mpilogger <<"-";
            mpilogger<<endl<<write;          
        }
        return true;
    }
}

#else
//Profiling disabled
namespace profile 
{
    bool start(const std::string &label){return true;}
    bool stop (const std::string &label,double workUnits,
               const std::string &workUnitLabel){return true;}
    bool print(MPI_Comm comm){return true;}
}

#ifdef NO_MPILOGGER
    #undef mpilogger
    #undef write
#endif

#endif

