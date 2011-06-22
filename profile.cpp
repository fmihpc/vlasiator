#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#ifdef CRAYPAT
//include craypat api headers if compiled with craypat on Cray XT/XE
#include "pat_api.h"
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
            double time; // total time accumulated
            int level;  //what hierarchy level
            int count; //how many times have this been accumulated
            int index; // unique index identifying this timer (usefule for Craypat)
        };
        //store all timers in map, key is fullLabel name
        std::map<string,TimerData> timers;

        //current position in timer hierarchy & label & startTime for the different levels
        int currentLevel=-1;
        std::vector<string> labels;
        std::vector<double> startTime;

        
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
        currentLevel++;
        //resize vectors if needed
        if(labels.size()<=currentLevel)
            labels.resize(currentLevel+5);
        if(startTime.size()<=currentLevel)
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
            timerData.parent=constructFullLabel(currentLevel-1);
            timerData.label=label;
            timerData.index=timers.size(); 
            timers[fullName]=timerData;
        }
        

#ifdef CRAYPAT
        PAT_region_begin(timers[fullName].index,label.c_str());
#endif
        startTime[currentLevel]=MPI_Wtime();
        return true;        
    }
    
    bool stop(const string &label){
        double stopTime=MPI_Wtime();
        string fullName=constructFullLabel(currentLevel);
#ifdef CRAYPAT
        PAT_region_end(timers[fullName].index);
#endif
        

        
        if(currentLevel<=-1){
            cout << "ERROR: nothing to stop in profile::stop Stopping "<< label <<endl;
            return false;
        }
        
        if(labels[currentLevel] != label){
            cout << "ERROR: label missmatch in profile::stop Stopping "<< label <<" at level " <<fullName << endl;
            return false;
        }
        
        
        
        timers[fullName].time+=(stopTime-startTime[currentLevel]);
        timers[fullName].count++;
        currentLevel--;

        return true;
    }

 
    //print out global timers
    //This will assume that all processes have identical labels, will fail otherwise (FIXME)
    bool print(MPI_Comm comm){
        size_t maxWidth=0;
        int rank,nProcesses;
        int nTimers;
        const int indentWidth=2; //how many spaces each level is indented
        vector<int> nAllTimers;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&nProcesses);
        
        for (std::map<string,TimerData>::iterator timer=timers.begin();
             timer!=timers.end();++timer){
            size_t width=timer->second.label.length()+timer->second.level*indentWidth;
            maxWidth=max(maxWidth,width);
        }
        
        nAllTimers.resize(nProcesses);
        nTimers=timers.size();
        MPI_Allgather(&nTimers,1,MPI_INT,&(nAllTimers[0]),1,MPI_INT,comm);

        //check that all processes have identical number of timers to avoid deadlocks
        //does not ensure the labels themselves are the same.
        for(int i=0;i<nProcesses;i++){
            if(nTimers!=nAllTimers[i]){
                if(rank==0){
                    cout <<"Error in profile::print, different number of labels on different processes" <<endl;
                }
                return false;
            }
        }

        //print out header
        if(rank==0){
            
            cout<<setw(maxWidth+1)<< setiosflags(ios::left) << "Label";
            cout<<setw(68) <<"Time(s)";
            cout<<setw(14)<<"Call count";
            cout<<endl;
            
            cout<<setw(maxWidth+1)<< " ";
            //time
            cout<<setw(14) <<"Average";
            cout<<setw(14) <<" % of parent";
            cout<<setw(14) <<"Maximum";
            cout<<setw(6) << "Rank";
            cout<<setw(14)<< "Minimum";
            cout<<setw(6) << "rank";
            //call count
            cout<<setw(14) << "Average";

            cout<<endl;
        }
        
        for (std::map<string,TimerData>::iterator timer=timers.begin();
             timer!=timers.end();++timer){
            
            double sum,parentSum;
            int intSum;
            doubleRankPair in;
            doubleRankPair out;
            

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
                cout<<setw(indent)<<"";
                cout<<setw(maxWidth+1-indent)<< setiosflags(ios::left)
                    << timer->second.label;

                cout << setw(14) << sum/nProcesses;
                if(timer->second.level>0)
                    cout << setw(14) << 100.0*sum/parentSum;
                else
                    cout << setw(14) << " ";
                
            }
            MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,0,comm);            
            if(rank==0){
                cout << setw(14) << out.val;
                cout << setw(6) <<out.rank;
            }
            
            MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,0,comm);
            if(rank==0){
                cout << setw(14) << out.val;
                cout << setw(7)<<out.rank;
            }
            
            //current count statistics
            MPI_Reduce(&( timer->second.count),&intSum,1,MPI_INT,MPI_SUM,0,comm);
            if(rank==0){
                cout << setw(14) << (double)intSum/nProcesses;
                cout<<endl;
            }

            
        }
         return true;
    }
}

#else
//Profiling disabled
namespace profile 
{
  
    bool start(const std::string &label){return true;}
    bool stop (const std::string &label){return true;}
    bool print(MPI_Comm comm){return true;}
}

#endif
