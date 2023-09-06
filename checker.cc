#include <bits/stdc++.h>
using namespace std;
using namespace std::chrono;
int main(int argc, char* argv[]) {
  if (argc < 5) {
    printf("checker generator.cc brute.cc a.cc caseNum\n");
    return -1;
  }
  int tot = atoi(argv[4]);
  system("g++ -std=c++14 -O2 -o _generator generator.cc");
  system("g++ -std=c++14 -O2 -o _brute brute.cc");
  // system("g++ -std=c++14 -O2 -fsanitize=address -o tocheck a.cc");
  system("g++ -std=c++14 -O2 -o tocheck a.cc");
  for (int cas = 1; cas <= tot; ++cas) {
    system("./_generator > task.in");
    system("./_brute < task.in > task.ans");
    // auto t1 = steady_clock::now();
    if (system("./tocheck < task.in > task.out")) {
      cerr << "Runtime Error" << endl;
      return 2;
    }
    /*
    auto t2 = steady_clock::now();
    double dif = duration_cast<duration<double>>(t2 - t1).count();
    if (dif > 1) {
      cerr << "Time Limit Exceeded" << endl;
      return 3;
    }
    */
    if (system("diff -w task.out task.ans")) {
      cerr << "Wrong Answer" << endl;
      return 1;
    }
    cerr << "Test " << cas << " Passed" << endl;
  }
  cerr << "Accepted" << endl;
  return 0;
}
