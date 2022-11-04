//compute_zscore.cpp
#include <Rcpp.h>
#include <fstream>
#include <thread>
#include <unistd.h>

using namespace Rcpp;

NumericMatrix compute_write_zscore(NumericMatrix num, int no_threads) ;
void write_file(NumericMatrix num,std::string output_file_name);
void compute_zscore_row(NumericMatrix::Row row);

//Function to write files in tsv format
// [[Rcpp::export]]
void write_file(NumericMatrix num,std::string output_file_name) 
{
    //Rcout<<"Writing z-scores using C++ Implementation"<<"\n";
    std::ofstream myfile;
    myfile.open(output_file_name);
    CharacterVector ch = colnames(num);
    CharacterVector ch1 = rownames(num);

    myfile <<"Hugo_Symbol"<<"\t";
    for(auto p = ch.begin(); p != ch.end(); ++p){
          myfile << *p<<"\t";
      }
    myfile <<"\n";
    for(size_t i = 0; i < num.nrow();++i){
        myfile << ch1(i)<<"\t";
        for(size_t j = 0; j <  num.ncol();++j ){
            myfile << num(i,j)<<"\t";
          }
        myfile <<"\n";
      }
    myfile.close();
}
// Function to compute zscore
void compute_zscore_row(NumericMatrix::Row row)
{
    float mean =0.0;
    float sd = 0.0;
    float var = 0.0;
    int row_size=row.size();
    for (auto i = row.begin(); i != row.end(); ++i){ // compute mean
        if ( *i != 0)
              *i = std::log2(*i+1);
        mean=mean + *i ;      
    }
    mean=mean/row_size;
    for (auto i = row.begin(); i != row.end(); ++i){ //compute SD
          var=var+ ((*i-mean) * (*i-mean));
      }
    sd=std::sqrt(var/(row_size-1));  
    int round_off = 10000;
    float a;
    for(auto i = row.begin(); i != row.end(); ++i) //compute zscore
      {
          a = std::ceil((*i - mean)*round_off/sd) ;
          *i=a/round_off;
      }  
}

//function to compute zscore and replace num matrix with zscores
// [[Rcpp::export]]
NumericMatrix compute_write_zscore(NumericMatrix num, int no_threads =1) 
{
  std::vector<std::thread*> threads;
  Rcout<<"Utilizing "<<no_threads<<" threads to compute zscore"<<"\n";
  for(int i = 0; i < num.nrow();i+=no_threads){
      if (i>num.nrow()-no_threads) //check for last rows i.e corner case
        {
          for (size_t j=i; j<num.nrow();j++)
            {
            compute_zscore_row(num( j , _ ));
            }
        break;
        }
      for(int jj = 0; jj < no_threads; ++jj) //fire threads
      {
            threads.push_back(new std::thread(compute_zscore_row,num(i+jj, _ )));
            //  Rcout <<i+jj<<"\n";
      }
      for(int jj = 0; jj < no_threads; ++jj) // join threads and clean heap
      {
          threads.at(jj) -> join();
          delete threads.at(jj);
      }    
      threads.clear(); //clear the vector
    }
  //NumericMatrix::Row  col = num( 0 , _ );
  //for (auto i = col.begin(); i != col.end(); ++i)
  //     Rcout << *i << "\n";

  //write_file(num,"written.tsv"); //write the num matrix to tsv format
  return num; //return matrix with zscores in it
}

/* keeping it for future
List compute_zscore_1(List num) {
int n = num.size();
float a =0;
float sum =0;
float st = 0;
float var = 0;
float mean =0;
for(int i = 0; i < n; ++i) {
  a=num[i];
  a = std::log2(a + 1);
  sum=sum+a;
  num[i]=a;
} 
mean=sum/n;
for(int i =0; i < n; ++i){
  a=num[i];
  var=var+(a-mean)*(a-mean);
}
st=std::sqrt(var/(n-1));
int round_off = 10000;
for(int i = 0; i < n; ++i) {
    a=num[i];
    a = std::ceil((a - mean)*round_off/st) ;
    num[i]=a/round_off;
}  
return num;
}
*/