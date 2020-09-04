#ifndef OMPHELPER_H_
#define OMPHELPER_H_

#define USE_OPENMP

#ifdef USE_OPENMP
#include <omp.h>
#ifdef USE_MSVC
#define OMP_PARALLEL __pragma(omp parallel)
#define OMP_FOR __pragma(omp for)
#define OMP_SINGLE __pragma(omp single)
#define OMP_SECTIONS __pragma(omp sections)
#define OMP_SECTION __pragma(omp section)
#else
#define OMP_PARALLEL _Pragma("omp parallel")
#define OMP_FOR _Pragma("omp for")
#define OMP_SINGLE _Pragma("omp single")
#define OMP_SECTIONS _Pragma("omp sections")
#define OMP_SECTION _Pragma("omp section")
#endif
#else
#include <ctime>
#define OMP_PARALLEL
#define OMP_FOR
#define OMP_SINGLE
#define OMP_SECTIONS
#define OMP_SECTION
#endif

#include <cassert>
#include <vector>



class Timer
{
public:

	typedef int EventID;

	EventID get_time()
	{
		EventID id = time_values_.size();

#ifdef USE_OPENMP
		time_values_.push_back(omp_get_wtime());
#else
		time_values_.push_back(clock());
#endif

		return id;
	}

	double elapsed_time(EventID event1, EventID event2)
	{
		assert(event1 >= 0 && event1 < static_cast<EventID>(time_values_.size()));
		assert(event2 >= 0 && event2 < static_cast<EventID>(time_values_.size()));

#ifdef USE_OPENMP
		return time_values_[event2] - time_values_[event1];
#else
		return double(time_values_[event2] - time_values_[event1]) / CLOCKS_PER_SEC;
#endif
	}

	void reset()
	{
		time_values_.clear();
	}

private:
#ifdef USE_OPENMP
	std::vector<double> time_values_;
#else
	std::vector<clock_t> time_values_;
#endif
};

#endif /* OMPHELPER_H_ */
