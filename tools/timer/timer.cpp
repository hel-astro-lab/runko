// NOTE not defining this will cause the timer to not do anything
#define USE_INTERNAL_TIMER

#include "timer.h"


long Fibonacci(unsigned n)
{
    return n < 2 ? n : Fibonacci(n - 1) + Fibonacci(n - 2);
}


int main()
{

  //const auto start{std::chrono::steady_clock::now()};

  //const auto fb{Fibonacci(42)};

  //const auto end{std::chrono::steady_clock::now()};
  //const std::chrono::duration<double> elapsed_seconds{end - start};

  //std::cout << "Fibonacci(42): " << fb << "\nElapsed time: ";
  //std::cout << elapsed_seconds.count() << "s\n"; // Before C++20
  ////std::cout << elapsed_seconds << '\n'; // C++20's chrono::duration operator<<

  ////-------------------------------------------------- 
  //// new version

  //std::cout << "------------------------------\n";
  //Timer timer;
  //timer.start();

  //const auto fb2{Fibonacci(42)};
  //auto elaps = timer.stop();
  //std::cout << elaps << "s\n"; // Before C++20
  //std::cout << "------------------------------\n";


  Timer timer("fibbo timer");


  for(int cycle=0; cycle<5; cycle++){

    timer.start();

    for(int i = 0; i<10; i++){

      auto t1 = timer.start_comp("fib");
      const auto fb3{Fibonacci(12)};
      timer.stop_comp(t1);


      auto t2 = timer.start_comp("quick");
      const auto fb4{Fibonacci(1)};
      timer.stop_comp(t2);
    }

    timer.stop();
  }
  timer.comp_stats();


  // round 2
  timer.clear();

  timer.start();

  auto t1 = timer.start_comp("fib");
  const auto fb5{Fibonacci(42)};
  timer.stop_comp(t1);

  timer.stop();
  timer.comp_stats();


}

