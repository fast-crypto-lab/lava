#ifndef _MONOMIAL_2_H_
#define _MONOMIAL_2_H_


#include <stdint.h>
#include <iostream>
#include "gf.h"


////////////////////// GF 2 //////////////////////////////////////////


/// number of monomial less EQUAL degree ?
template <unsigned deg,unsigned n_mono>
struct n_mono_le_deg
{
	const static uint64_t v = n_mono_le_deg<deg-1,n_mono-1>::v + n_mono_le_deg<deg,n_mono-1>::v;
};


template <> struct n_mono_le_deg<0,0> { const static uint64_t v=0; };
template <unsigned n_mono> struct n_mono_le_deg<0,n_mono> { const static uint64_t v=1; };
template <unsigned deg> struct n_mono_le_deg<deg,0> { const static uint64_t v=0; };

template <> struct n_mono_le_deg<1,1> { const static uint64_t v=2; };
template <unsigned n_mono> struct n_mono_le_deg<1,n_mono> { const static uint64_t v=n_mono+1; };
template <unsigned deg> struct n_mono_le_deg<deg,1> { const static uint64_t v=2; };

template <> struct n_mono_le_deg<1,0> { const static uint64_t v=0; };
template <> struct n_mono_le_deg<0,1> { const static uint64_t v=1; };



///////////////////////////////////////////////////////////////////


template <unsigned n_mono , unsigned deg, unsigned I=(n_mono+1)*(deg+1)-1>
struct n_mono_le_deg_table : public n_mono_le_deg_table<n_mono , deg , I-1> { static const uint64_t dummy; };


template <unsigned n_mono , unsigned deg >
struct n_mono_le_deg_table<n_mono , deg , 0>
{
	static const uint64_t dummy;
	static uint64_t v[(n_mono+1)*(deg+1)];
	static uint64_t get(unsigned nmono , unsigned ndeg ) { return v[nmono*(deg+1)+ndeg]; }
};


template <unsigned n_mono, unsigned deg ,  unsigned I>
const uint64_t n_mono_le_deg_table<n_mono,deg,I>::dummy
		= n_mono_le_deg_table<n_mono,deg,0>::v[I] = n_mono_le_deg<I%(deg+1),I/(deg+1)>::v
							+ 0*n_mono_le_deg_table<n_mono,deg,I-1>::dummy;


template <unsigned n_mono,unsigned deg>
uint64_t n_mono_le_deg_table<n_mono,deg,0>::v[(n_mono+1)*(deg+1)] = {0};


/// initialize needed
/// for ex:
/// template struct n_mono_le_deg_table<10,12>;  /// explicit initialize



template <unsigned deg, unsigned n_mono>
struct n_mono_table
{
	static const uint64_t *v;
};

template <unsigned deg, unsigned n_mono>
const uint64_t * n_mono_table<deg,n_mono>::v = & n_mono_le_deg_table<n_mono,deg>::v[n_mono*(deg+1)];





///////////////////////////////////////////////////////////////////////////////////


#ifdef _SSSE3_POPCNT_

#include <emmintrin.h>
#include <tmmintrin.h>

static const __m128i _count_bit_tab = _mm_set_epi8(4,3,3,2,3,2,2,1,3,2,2,1,2,1,1,0);
static const __m128i _mask_low_nib = _mm_set_epi8(15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15);
//static const __m128i _mask_high_nib = _mm_set_epi8(0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0);

inline unsigned count_bit_on( uint64_t v64 )
{
	uint64_t vv[2];
	vv[0] = v64;
	vv[1] = v64>>4;
	__m128i v128 = _mm_loadu_si128( (__m128i*) vv );
	__m128i r = _mm_shuffle_epi8( _count_bit_tab , _mask_low_nib&v128 );
	r = _mm_hadd_epi16( r , v128 ); // 16 --> 8
	r = _mm_hadd_epi16( r , v128 ); // 8 --> 4
	r = _mm_hadd_epi16( r , v128 ); // 4 --> 2
	return _mm_extract_epi16( _mm_add_epi8( r , _mm_srli_si128( r , 1 ) ), 0 ) & 0xff;
}

#endif



inline bool is_nth_bit_on( uint64_t chk , uint8_t n )
{
	union{
	uint64_t vv;
	uint8_t v[8];
	};
	vv = chk;
	return (v[n>>3]&(1<<(n&7)));
}

inline void set_nth_bit_on( uint8_t * v , uint8_t n )
{
	v[n>>3] |= (1<<(n&7));
}
inline void flip_nth_bit( uint8_t * v , uint8_t n )
{
	v[n>>3] ^= (1<<(n&7));
}


//#ifdef NO_SSE
#if 1

template <unsigned n_mono,unsigned max_deg>
struct mono_var
{
	typedef n_mono_le_deg_table<n_mono,max_deg> tb;

	union{
	uint64_t vv;
	uint8_t v[8];
	};
	uint8_t deg;

	mono_var(uint64_t idx = 0): vv(0), deg(0)
	{
		if( 0 == idx ) return;
		unsigned d=1;
		while( idx >= tb::get(n_mono,d) ) d++;
		deg = d;

		idx = tb::get(n_mono,d) - idx;
		for(unsigned i=n_mono-1; d != 0; i-- ){
			uint64_t num_mono_var_i_deg_d = tb::get(i,d)-tb::get(i,d-1);
			if( idx > num_mono_var_i_deg_d ) {
				set_nth_bit_on(v,i);
				idx -= num_mono_var_i_deg_d;
				d--;
			}
		}
	}

	uint64_t to_idx() const
	{
		unsigned d = deg;
		uint64_t idx=tb::get(n_mono,d)-1;
		uint64_t tmp_vv = vv;
		while( d != 0 ) {
			unsigned i = 63 - __builtin_clzll( tmp_vv );
			tmp_vv ^= ((uint64_t)1)<<i;
			idx -= tb::get(i,d)-tb::get(i,d-1);
			d--;
		}
		return idx;
	}

	static uint64_t min_idx_with_min_mono( uint8_t deg , uint8_t deg1_idx ) {
		if( deg1_idx>n_mono) deg1_idx=n_mono;
		uint64_t r = tb::get( n_mono , deg ) - 1;
		if( deg > n_mono - deg1_idx ) return r;
		for( unsigned i=0; deg !=0 ; i++ ) {
			r -= tb::get( n_mono-deg1_idx-i , deg ) - tb::get( n_mono-deg1_idx-i , deg-1);
			deg--;
		}
		return r;
	/// scalar extension:  idx >= mono_t::min_idx_with_min_mono(deg,deg1_ext_idx+1);
	/// ext displacement:
	/// mono_t::min_idx_with_min_mono(deg+1,deg1_ext_idx)-mono_t::min_idx_with_min_mono(deg,deg1_ext_idx+1);
	/// no extension:   idx >= mono_t::min_idx_with_min_mono(deg,deg1_ext_idx);
	}

	static uint64_t ext_deg1( uint64_t org_idx , uint8_t deg1_idx  ) {
		if( deg1_idx > n_mono) return 0;   ///assert( n_mono >= deg1_idx );
		if( 0 == deg1_idx ) return org_idx;
		if( 0 == org_idx ) return deg1_idx;
		unsigned d=1;
		while( org_idx >= tb::get(n_mono,d) ) d++;
		//unsigned org_deg = d;
		uint64_t idx = tb::get(n_mono,d) - org_idx;
		uint64_t r_idx = tb::get(n_mono,d+1)-1; /// expect degree increased
		for(unsigned i=n_mono-1; d != 0; i-- ){
			uint64_t num_mono_var_i_deg_d = tb::get(i,d)-tb::get(i,d-1);
			if( i == (n_mono-deg1_idx) ) {
				if( idx > num_mono_var_i_deg_d ) return org_idx; /// repeat variable
				return r_idx - tb::get(i,d+1) + tb::get(i,d) - idx + 1;
				//r_idx -= tb::get(i,d+1)-tb::get(i,d) + idx - 1;
			}
			if( idx > num_mono_var_i_deg_d ) {
				//set_nth_bit_on(v,i);
				r_idx -= tb::get(i,d+1)-tb::get(i,d);
				idx -= num_mono_var_i_deg_d;
				d--;
			}
		}
		/// 0 == d
		r_idx -= tb::get(n_mono-deg1_idx,1)-tb::get(n_mono-deg1_idx,0);
		return r_idx;
	}
#if 0
	uint64_t to_idx() const
	{
		unsigned d=0;
		uint64_t idx=0;
		for(unsigned i=0;i<n_mono;i++){
			if( !is_nth_bit_on(vv,i) ) continue;
			/// if i-th bit is on
			d++;
			idx += tb::get(i,d)-tb::get(i,d-1);
		}
		return tb::get(n_mono,d)-1-idx;
	}
#endif
	unsigned degree() const { return deg; }

	int right_most_idx() const { return __builtin_ffsll(vv)-1; }

	mono_var get_deg1mono() const
	{
		if( 0 == deg ) return mono_var(0);
		mono_var r;
		r.vv = vv ^ (vv&(vv-1));
		r.deg = 1;
		return r;
	}

	inline void count_deg() {
#ifdef _SSSE3_POPCNT_
		deg = count_bit_on(vv);
#else
		deg = __builtin_popcountll( vv );
#endif
	}

	bool is_zero() const { return 0==deg; }

	bool is_divisible_by(const mono_var & a ) const
	{
		return ((vv^a.vv)==(vv-a.vv));
	}

	mono_var& operator *=( const mono_var & a )
	{
		vv |= a.vv;
		count_deg();
		return *this;
	}


	friend const mono_var operator *( const mono_var& a , const mono_var& b )
	{
		mono_var r;
		r = a;
		r *= b;
		return r;
	}

	friend const mono_var operator/( const mono_var& a , const mono_var& b )
	{
		mono_var r;
		r = a;
		r.vv ^= b.vv;
		r.count_deg();
		return r;
	}

	friend const mono_var lcm( const mono_var& a , const mono_var& b )
	{
		mono_var r;
		r = a;
		r.vv |= b.vv;
		r.count_deg();
		return r;
	}

	friend const mono_var gcd( const mono_var& a , const mono_var& b )
	{
		mono_var r;
		r = a;
		r.vv &= b.vv;
		r.count_deg();
		return r;
	}

	friend bool operator == (const mono_var& a , const mono_var& b ) { return (a.vv==b.vv); }
	friend bool operator != (const mono_var& a , const mono_var& b ) { return !(a==b); }

	friend bool operator < ( const mono_var& a , const mono_var& b )
	{
		if( a.deg == b.deg ) return a.vv > b.vv;
		return a.deg < b.deg;
	}
	friend bool operator > ( const mono_var& a , const mono_var& b ) { return b<a; }
	friend bool operator <= ( const mono_var& a , const mono_var& b ) { return !(a>b); }
	friend bool operator >= ( const mono_var& a , const mono_var& b ) { return !(b>a); }

	friend std::ostream& operator << ( std::ostream& out, const mono_var& a){
		out << "[" << a.degree() << "|";
		for(unsigned i =0; i < n_mono; i++){
			out << is_nth_bit_on(a.vv,i);
			if( 3 == (i&3) ) out << ",";
		}
		out << "]";
		return out;
	}

	void magma_out( std::ostream& out ) const {
		if( 0 == vv ) { out << "1"; return; }
		bool st = true;
		for(int i=n_mono;i>=0;i--) {
			if( is_nth_bit_on(vv,(unsigned)i) ) {
				if( !st ) out << "*";
				if( st ) st = false;
				out << "$." << (i+1);
			}
		}
	}

	static int magma_in( mono_var & a , std::istream& in ) {
		in >> std::ws; // white space
		int c = in.peek();

		//if( ' ' == c ) {
		//	while( ' ' == c ) c = in.get();
		//	in.unget();
		//}

		a = mono_var();
		if( '1' == c ) {
			c = in.get(); // extract
			return 0;
		}

		while( '$' == c ){
			c = in.get();

			c = in.peek();
			if( '.' != c ) return -1;
			c = in.get();

			c = in.peek();
			if( !std::isdigit(c) ) return -1;
			int ii = 0;
			in >> ii;
			a *= mono_var(n_mono-ii+1);

			c = in.peek();
			if( '*' != c ) return 0;
			c = in.get(); // extract *
			c = in.peek();
		}
		//if( EOF == c ) return -1;
		// else
		return -1;

	}

	template <typename T>
	const T eval( const T * vec ) const { T r(1); for(unsigned i=0;i<n_mono;i++) if(is_nth_bit_on(vv,i)) r*= vec[i]; return r; }

};




template <unsigned n_mono,unsigned max_deg>
struct mono_2
{
	typedef n_mono_le_deg_table<n_mono,max_deg> tb;

	union{
	uint64_t vv;
	uint8_t v[8];
	};

	mono_2(uint64_t idx = 0): vv(0)
	{
		if( 0 == idx ) return;
		unsigned d=1;
		while( idx >= tb::get(n_mono,d) ) d++;
		//deg = d;

		//std::cout << "test, idx = "<< idx  <<", d = " << d<< " , i start at = "<< n_mono-1<< std::endl;
		idx = tb::get(n_mono,d) - idx;

		for(unsigned i=n_mono-1; d != 0; i-- ){
			uint64_t num_mono_var_i_deg_d = tb::get(i,d)-tb::get(i,d-1);
			if( idx > num_mono_var_i_deg_d ) {
				set_nth_bit_on(v,i);
				idx -= num_mono_var_i_deg_d;
				d--;
			}
		}
	}

	uint64_t to_idx() const
	{
		unsigned d = degree();
		uint64_t idx=tb::get(n_mono,d)-1;
		uint64_t tmp_vv = vv;
		while( d != 0 ) {
			unsigned i = 63 - __builtin_clzll( tmp_vv );
			tmp_vv ^= ((uint64_t)1)<<i;
			idx -= tb::get(i,d)-tb::get(i,d-1);
			d--;
		}
		return idx;
	}

	uint64_t fast_idx() const { uint64_t d=degree(); return (d<<n_mono)-vv; } /// XXX: will fail if n_mono ~64

	unsigned degree() const { return __builtin_popcountll(vv); }

	int right_most_idx() const { return __builtin_ffsll(vv)-1; }

	mono_2 get_deg1mono() const
	{
		if( 0 == vv ) return mono_2();
		mono_2 r;
		r.vv = vv ^ (vv&(vv-1));
		return r;
	}

	unsigned split( mono_2 varlist[n_mono] ) const
	{
		unsigned deg = 0;
		mono_2 m1; m1 = *this;
		mono_2 m = m1.get_deg1mono();
		while( 0 != m.vv ) {
			varlist[deg++] = m;
			m1.vv ^= m.vv;
			m = m1.get_deg1mono();
		}
		return deg;
	}
#if 0
	std::vector<mono_2> split() const
	{
		std::vector<mono_2> varlist;
		mono_2 m1; m1 = *this;
		mono_2 m = m1.get_deg1mono();
		while( 0 != m.vv ) {
			varlist.push_back(m);
			m1.vv ^= m.vv;
			m = m1.get_deg1mono();
		}
#if 0
		for(unsigned i =0; i < n_mono; i++){
			if( is_nth_bit_on(vv,i) ) {//probablu exists some smarter way
				mono_2 m;
				set_nth_bit_on(m.v,i);
				varlist.push_back(m);
			}
		}
#endif
		return varlist;
	}
#endif

	const mono_2& operator = ( const mono_2& m ) { 
		vv = m.vv;
		return *this;
	}


	bool is_zero() const { return 0==vv; }

	bool is_divisible_by(const mono_2 & a ) const
	{
		return 0==((vv&a.vv)^a.vv);
		//return ((vv^a.vv)==(vv-a.vv));
	}

	mono_2& operator *=( const mono_2 & a )
	{
		vv |= a.vv;
		return *this;
	}


	friend const mono_2 operator *( const mono_2& a , const mono_2& b )
	{
		mono_2 r = a;
		r *= b;
		return r;
	}

	friend const mono_2 operator/( const mono_2& a , const mono_2& b )
	{
		mono_2 r = a;
		r.vv ^= b.vv;
		return r;
	}

	friend const mono_2 lcm( const mono_2& a , const mono_2& b )
	{
		mono_2 r = a;
		r.vv |= b.vv;
		return r;
	}

	friend const mono_2 gcd( const mono_2& a , const mono_2& b )
	{
		mono_2 r = a;
		r.vv &= b.vv;
		return r;
	}

	friend bool operator == (const mono_2& a , const mono_2& b ) { return (a.vv==b.vv); }
	friend bool operator != (const mono_2& a , const mono_2& b ) { return !(a==b); }

	friend bool operator < ( const mono_2& a , const mono_2& b )
	{
		unsigned adeg = a.degree();
		unsigned bdeg = b.degree();
		if( adeg == bdeg ) return a.vv > b.vv;
		return adeg < bdeg;
	}
	friend bool operator > ( const mono_2& a , const mono_2& b ) { return b<a; }
	friend bool operator <= ( const mono_2& a , const mono_2& b ) { return !(a>b); }
	friend bool operator >= ( const mono_2& a , const mono_2& b ) { return !(b>a); }

	friend std::ostream& operator << ( std::ostream& out, const mono_2& a){
//		out << "[" << a.degree() << "|";
		for(unsigned i =0; i < n_mono; i++){
#if 1
			if( 0 == a.degree() ) { out << 1; break; }
			if( is_nth_bit_on(a.vv,i) ) {
				out << "X[" << i+1 << "]*";
			}
#else
			out << is_nth_bit_on(a.vv,i);
			if( 3 == (i&3) ) out << ",";
#endif
		}
//		out << "]";
		return out;
	}

	void magma_out( std::ostream& out ) const {
		if( 0 == vv ) { out << "1"; return; }
		bool st = true;
		for(int i=n_mono;i>=0;i--) {
			if( is_nth_bit_on(vv,(unsigned)i) ) {
				if( !st ) out << "*";
				if( st ) st = false;
				out << "$." << (i+1);
			}
		}
	}

	static int magma_in( mono_2 & a , std::istream& in ) {
		in >> std::ws; // white space
		int c = in.peek();

		//if( ' ' == c ) {
		//	while( ' ' == c ) c = in.get();
		//	in.unget();
		//}

		a = mono_2();
		if( '1' == c ) {
			c = in.get(); // extract
			return 0;
		}

		while( '$' == c ){
			c = in.get();

			c = in.peek();
			if( '.' != c ) return -1;
			c = in.get();

			c = in.peek();
			if( !std::isdigit(c) ) return -1;
			int ii = 0;
			in >> ii;
			a *= mono_2(n_mono-ii+1);

			c = in.peek();
			if( '*' != c ) return 0;
			c = in.get(); // extract *
			c = in.peek();
		}
		//if( EOF == c ) return -1;
		// else
		return -1;

	}

	template <typename T>
	const T eval( const T * vec ) const { T r(1); uint64_t qq = 1; for(unsigned i=0;i<n_mono;i++,qq<<=1) if(qq&vv) r*= vec[i]; return r; }

};



#else






#endif



template <unsigned n_mono,unsigned max_deg>
struct mono_idx
{


};




#endif
