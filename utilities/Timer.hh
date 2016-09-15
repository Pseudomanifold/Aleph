#ifndef ALEPH_TIMER_HH__
#define ALEPH_TIMER_HH__

#include <chrono>

namespace aleph
{

namespace utilities
{

class Timer
{
public:

  /** Starts the timer */
  Timer()
    : _start( std::chrono::steady_clock::now() )
  {
  }

  /** Restarts the timer */
  void restart()
  {
    _start = std::chrono::steady_clock::now();
  }

  // ---------------------------------------------------------------------

  /** @returns Elapsed time in (fractional) seconds */
  double elapsed_s() const
  {
    return std::chrono::duration<double>( now() - _start ).count();
  }

  /** @returns Elapsed time in (fractional) milliseconds */
  double elapsed_ms() const
  {
    return std::chrono::duration<double, std::milli>( now() - _start ).count();
  }

  /** @returns Elapsed time in (fractional) microseconds */
  double elapsed_mu() const
  {
    return std::chrono::duration<double, std::micro>( now() - _start ).count();
  }

  /** @returns Elapsed time in (fractional) nanoseconds */
  double elapsed_ns() const
  {
    return std::chrono::duration<double, std::nano>( now() - _start ).count();
  }

  // ---------------------------------------------------------------------

  /** @returns Current time */
  std::chrono::steady_clock::time_point now() const
  {
    return std::chrono::steady_clock::now();
  }

private:

  /** Start time; will be set automatically in the constructor */
  std::chrono::steady_clock::time_point _start;
};

}

}

#endif
