#include <iostream>
#include <NTL/ZZ.h>
#include<fstream>


using namespace std;
using namespace NTL;
/*  
sudo apt-get install libntl-dev   
g++  search_p_D.cpp -o example -lntl -lm
*/

int main(int argc, char* argv[]){
  
  ZZ D,q,inv,sqrt;
  ZZ x = ZZ(4);
  
  long a1 = atoi(argv[1]);  
  long a2 = atoi(argv[2]);
  long a3 = atoi(argv[3]);
  long a4 = atoi(argv[4]);
  
  ZZ psize = power2_ZZ(a1); 
  ZZ max_iter = power2_ZZ(a2);
  ZZ min_D = power2_ZZ(a3);
  ZZ max_D = power2_ZZ(a4);

  for (ZZ i = min_D ; i < max_D ; i++){
   D = x *  i  + ZZ(3);
   inv = InvMod(x,D);
   for (ZZ m = psize ; m <psize+max_iter ; m++){
      q = D * m*(m+1) + inv;
      sqrt = SqrRoot(q);
      if (sqrt*sqrt == q){
        if (ProbPrime(sqrt) == 1){
          std::ofstream ofile;
          ofile.open("resultat.txt", std::ios::app); 
          cout << "p : " << sqrt << " D : " << D << "\n";
          ofile << D << " p : " << sqrt <<std::endl; 
          ofile.close();
          return 1;
        }
      }
    }
  }
  return 0;
}


