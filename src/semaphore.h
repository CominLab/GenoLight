#include <cstdint>
#include  <mutex> 
#include <condition_variable>

class semaphore
{

    public:
        semaphore(std::size_t resources);
        void wait(std::size_t to_take);
        void signal(std::size_t to_give);
        ~semaphore();

    private:
        std::mutex mtx;
        std::condition_variable cond;
        std::size_t counter;
        const std::size_t global_resources;
};

