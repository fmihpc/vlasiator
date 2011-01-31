#include <cstdlib>
#include <iostream>

#include "mpilogger.h"
#include "timer.h"

using namespace std;

extern MPILogger mpilogger;
std::vector<Timer::TimerData> Timer::timers;

/** Constuctor for Timer. The constructor does not do anything.*/
Timer::Timer() { }

/** Destructor for Timer. The destructor does not do anything.*/
Timer::~Timer() { }

#ifdef PROFILE

/** Create a new timer.
 * @param name The name of the timer.
 * @return The ID number of the created timer. This ID should be used 
 * when calling start and stop member functions.
 * @see Timer::start
 * @see Timer::stop
 */
unsigned int Timer::create(const std::string& name) {
   timers.push_back(TimerData());
   timers[timers.size()-1].timeInSeconds = 0.0;
   timers[timers.size()-1].name = name;
   return timers.size()-1;
}

/** Get the value of the given timer.
 * @param timerID The ID of the timer.
 * @return Value of the timer in seconds.
 */
double Timer::getValue(const unsigned int& timerID) {
   return timers[timerID].timeInSeconds;
}

/** Write the name and value of each timer to logfile.
 */
void Timer::print() {
   for (vector<TimerData>::const_iterator it=timers.begin(); it!=timers.end(); ++it) {
      mpilogger << "(TIMER) Value of timer '" << it->name << "' is " << it->timeInSeconds << " s." << endl;
   }
   mpilogger << write;
}

/** Start the given timer.
 * @param timerID The ID number of the timer.
 */
void Timer::start(const unsigned int& timerID) {
   timers[timerID].startClock = clock();
}

/** Stop the given timer. The time elapsed between calls to 
 * Timer::start and Timer::stop for the given timer is accumulated 
 * to the value of the timer, which can be requested with 
 * Timer::getValue member function.
 * @param timerID the ID number of the timer.
 */
void Timer::stop(const unsigned int& timerID) {
   clock_t endClock = clock();
   timers[timerID].timeInSeconds += (1.0*(endClock - timers[timerID].startClock))/CLOCKS_PER_SEC;
}

#else

unsigned int Timer::create(const std::string& name) {return 0;}
double Timer::getValue(const unsigned int& timerID) {return 0.0;}
void Timer::print() { }
void Timer::start(const unsigned int& timerID) { }
void Timer::stop(const unsigned int& timerID) { }

#endif // #ifdef PROFILE
