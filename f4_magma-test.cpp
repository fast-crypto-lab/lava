#include "f4_algo.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <omp.h>



using namespace std;


gf_t vec[N_VAR] = {
//0,X,X,X, X,0,X,X, X,X,0,X, X,X,X,X, 0,X,X,X, X,X,0,X, 1,1,0,1 /// gh
//0,X,X,X, X,0,X,X, X,X,0,X, X,X,X,X, 0,X,X,X, X,X,0,X, 1,1,0,1 /// daven
  0,1,1,0, 1,0,0,1, 0,1,0,0, 0//,1,0,1, 0,1,1,1, 0,0,0,1, 1,1,0,1 /// daven
};

int selector[8];
int selector_count = 8;


int main( int argc , char * argv[] )
{

	deque<base_poly_t> ii;
	deque<base_poly_t> oo;

	const unsigned n = N_VAR;
	const unsigned m = N_POLY;


	std::vector<base_poly_t> aa;

	const char * filename = (argc>=1)? argv[1] : "input_28.txt";
	if(argc >= 9)
		for(int i=0; i< 8; i++)
			selector[i] = atoi(argv[i+2]);
	ifstream is(filename,ios::in);

	while( ! is.eof() ) {
	//if( ! is.eof() ) {
		base_poly_t bb;
		bb = base_poly_t::magma_in( is );
		if( bb.is_zero() ) break;

		aa.push_back( bb );
	}
	is.close();

	std::cout << "sol: ";
	for(unsigned i=0;i<N_VAR;i++) std::cout << vec[i] << " ";
	std::cout << "\n";

	std::cout << "input:\n";
	for(unsigned i=0;i<aa.size();i++) std::cout << aa[i] << " = " << aa[i].eval( vec ) << "\n";
	for(int i=0; i< aa.size(); i++ ) ii.push_back( aa[i] );

	//std::cout << "here\n"; exit(0);


	cout<<ii.size()<<" inputs.\n";
	for(unsigned i=0;i<ii.size();i++) cout << ii[i] << "\n";

	F4_algo( oo , ii  );

	cout << "\n" << oo.size() << " outpus.\n";
	for(unsigned i=0;i<oo.size();i++) cout << oo[i] << " = " << oo[i].eval(vec) << "\n";

	return 0;
}



