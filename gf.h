
#ifndef _GF_H_
#define _GF_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>

#include <iostream>

//using namespace std;

template <unsigned size>
struct gf
{
	gf( uint8_t c = 0 ) { v = c; }

	uint8_t v;

	bool is_one() const { return (1==v); }
	bool is_zero() const { return (0==v); }
	friend const gf operator + ( const gf &a , const gf &b ) { gf r; r.v=a.v^b.v; return r; }
	const gf& operator += ( const gf &b ) { v^=b.v; return *this; }
	friend const gf operator - ( const gf &a , const gf &b) { gf r; r.v=a.v^b.v; return r; }
	const gf& operator -= ( const gf & b ) { v^=b.v; return *this; }
//////////////
	const gf inv() const;
	const gf & operator *=( const gf & b );
//////////////
	friend const gf operator * (const gf& a , const gf& b){ gf r=a; r*=b; return r; }
	friend bool operator == ( const gf & a ,const gf & b ) { return (a.v==b.v); }
	friend bool operator != ( const gf & a ,const gf & b ) { return (a.v!=b.v); }

	static const gf rand() { return gf( (uint8_t)(::rand()%size) ); }
	void show() const { printf("%2d", v ); }

//	static string out_tab[16];
	friend  std::ostream& operator << ( std::ostream& out, const gf& a){
		//cout<<(out_tab[a.v]);
		out << (unsigned)a.v;
		return out;
	}
};


template <>
inline
const gf<2> gf<2>::inv() const { return *this; }

template <>
inline
const gf<2> & gf<2>::operator *=( const gf<2> & b ) { v &= b.v; return *this; }



extern uint8_t gf16_inv_tab[16];
template <>
inline
const gf<16> gf<16>::inv() const { return gf<16>(gf16_inv_tab[v]); }

extern uint8_t gf16_mul_tab[256];
template <>
inline
const gf<16> & gf<16>::operator *=( const gf<16> & b ) { v = gf16_mul_tab[(v<<4)+b.v]; return *this; }







#if 0
template <>
uint8_t gf<16>::inv_tab[16] = {0,1,9,14,13,11,7,6,15,2,12,5,10,4,3,8};

template <>
uint8_t gf<16>::mul_tab[256] = {
0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
0,  2,  4,  6,  8, 10, 12, 14,  3,  1,  7,  5, 11,  9, 15, 13,
0,  3,  6,  5, 12, 15, 10,  9, 11,  8, 13, 14,  7,  4,  1,  2,
0,  4,  8, 12,  3,  7, 11, 15,  6,  2, 14, 10,  5,  1, 13,  9,
0,  5, 10, 15,  7,  2, 13,  8, 14, 11,  4,  1,  9, 12,  3,  6,
0,  6, 12, 10, 11, 13,  7,  1,  5,  3,  9, 15, 14,  8,  2,  4,
0,  7, 14,  9, 15,  8,  1,  6, 13, 10,  3,  4,  2,  5, 12, 11,
0,  8,  3, 11,  6, 14,  5, 13, 12,  4, 15,  7, 10,  2,  9,  1,
0,  9,  1,  8,  2, 11,  3, 10,  4, 13,  5, 12,  6, 15,  7, 14,
0, 10,  7, 13, 14,  4,  9,  3, 15,  5,  8,  2,  1, 11,  6, 12,
0, 11,  5, 14, 10,  1, 15,  4,  7, 12,  2,  9, 13,  6,  8,  3,
0, 12, 11,  7,  5,  9, 14,  2, 10,  6,  1, 13, 15,  3,  4,  8,
0, 13,  9,  4,  1, 12,  8,  5,  2, 15, 11,  6,  3, 14, 10,  7,
0, 14, 15,  1, 13,  3,  2, 12,  9,  7,  6,  8,  4, 10, 11,  5,
0, 15, 13,  2,  9,  6,  4, 11,  1, 14, 12,  3,  8,  7,  5, 10 };

template <>
string gf<16>::out_tab[16] = {
"0", "1", "$.1", "$.1^4", "$.1^2", "$.1^8", "$.1^5", "$.1^10", "$.1^3", "$.1^14", "$.1^9", "$.1^7", "$.1^6", "$.1^13", "$.1^11", "$.1^12"};

#endif




#endif
