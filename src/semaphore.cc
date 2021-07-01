#include <cstdint>
#include "semaphore.h"
#include <iostream>
semaphore::semaphore(std::size_t resources) :
global_resources(resources)
{
	counter = resources;
}

void semaphore::wait(std::size_t to_take)
{
	std::unique_lock<decltype(mtx)> lock(mtx);
	while (counter < to_take)
		cond.wait(lock);

	counter -= to_take;
}

void semaphore::signal(std::size_t to_give)
{
	std::unique_lock<decltype(mtx)> lock(mtx);
	counter += to_give;
	if (counter > global_resources)
		counter = global_resources;
	cond.notify_all();
}

semaphore::~semaphore()
{
	//dtor
}
