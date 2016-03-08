#ifndef _F45_COMPONENT_GF2_H_
#define _F45_COMPONENT_GF2_H_


#include <deque>
#include <map>
#include <set>
#include <utility>
#include <vector>
//#include "gmp.h"
#include <iostream>

#ifdef __HAS_SIG__
#include "../gf16/monomial16.h"
#endif

#include "sys_config.h"


#include "gf2_sp_poly.h"

typedef gf<2> gf_t;

typedef mono_2<N_VAR,MAX_DEG> mono_t;
typedef poly2_t<N_VAR,MAX_DEG> base_poly_t;
typedef list_poly<N_VAR,MAX_DEG> list_poly_t;
typedef dense_vec_t<N_VAR,MAX_DEG> vec_t;

#ifdef __HAS_SIG__
typedef monomial<N_VAR,MAX_DEG> mono16_t;
#endif


#define DUMP(x) (x)

#define _USE_F4_


void F4_algo( std::deque<base_poly_t> & out , std::deque<base_poly_t> & inp );




#ifdef __HAS_SIG__
struct Src
{
	unsigned idx;
	//mono_t t;
	mono16_t t;


	friend bool operator < ( const Src & a , const Src & b ) { if(a.idx==b.idx) return a.t<b.t; return a.idx<b.idx; }
	friend bool operator > ( const Src & a , const Src & b ) { if(a.idx==b.idx) return a.t>b.t; return a.idx>b.idx; }
	friend bool operator == ( const Src & a , const Src & b ) { if(a.idx==b.idx) return a.t==b.t; return false; }
	friend bool operator != ( const Src & a , const Src & b ) { return !(a==b); }

	//const Src & operator *= ( const mono_t & m ) { t*=m; return *this; }
	//const Src operator * ( const mono_t & m ) const { Src r; r.idx=idx; r.t=t*m; return r; }
	const Src & operator *= ( const mono_t & m ) { t*=mono16_t(m); return *this; }
	const Src operator * ( const mono_t & m ) const { Src r; r.idx=idx; r.t=t*mono16_t(m); return r; }

	friend std::ostream & operator << ( std::ostream & os , const Src & s ) { os << "(" << s.t << "," << s.idx << ")"; return os; }
	const Src& operator = ( const Src& s ) { 
		t = s.t;
		idx = s.idx;
		return *this;
	}
};
#endif




//extern gf_t vec[N_VAR];
struct pair_t;//need this declaration for pointer?

template<typename poly_t>
struct Src_poly_t
{
#ifdef __HAS_SIG__
	Src s;
#endif
	poly_t poly;
	bool unnecessary_for_new_pair;
	int unnecessary_by;
	Src_poly_t(){
		unnecessary_for_new_pair = false;
		unnecessary_by = -1;
	}//for initialization

	const mono_t & head_term() const { return poly.head_term(); }
	uint64_t head_idx() const { return poly.head_idx(); }
	friend std::ostream & operator << ( std::ostream & os , const Src_poly_t & p ) { os <<p.poly; return os; }

	const Src_poly_t & operator = ( const Src_poly_t & p ) { 
#ifdef __HAS_SIG__
		s = p.s;
#endif
		poly = p.poly;
		return *this;
	}

	void mark_unnecessary(const mono_t & m, int i){
		unnecessary_for_new_pair = true;
		unnecessary_by = i;
	}
/////////////////////////

	bool is_zero() const { return poly.is_zero(); }
	void set_zero(){poly.set_zero();}
	unsigned degree() const { return head_term().degree(); }

	unsigned n_terms() const { return poly.n_terms(); }
	const mono_t & operator[](unsigned i) const { return poly[i]; }

	int get_position(const mono_t & m) const { return poly.get_position(m); }
	gf_t get_coef(const mono_t & m) const { return poly.get_coef(m); }

	bool CanGenFieldSpoly(mono_t & m) const {return !poly.is_mutiple_of(m);}

	Src_poly_t GenFieldSpoly(mono_t & m) const{ 
		Src_poly_t p; 
#ifdef __HAS_SIG__
		p.s = s*m;//not necessary in f4
#endif
		p.poly = poly;
		p.poly *= m;
		return p;
	}

/////////////////////////

	const Src_poly_t & operator -=( const Src_poly_t & a ) { 
#ifdef __HAS_SIG__
		if( a.s > s ) s = a.s; 
#endif		
		poly -= a.poly; return *this; }

#ifdef _DEPRECATE_
	const Src_poly_t & operator |=( const Src_poly_t & a ) { 
#ifdef __HAS_SIG__
		if( a.s > s ) s = a.s; 
#endif
		poly |= a.poly; return *this; }
#endif
	const Src_poly_t & operator *=( const mono_t & m ) { 
#ifdef __HAS_SIG__
	s *= m; 
#endif
poly *= m; return *this; }

//////////////////////////

};


typedef Src_poly_t<base_poly_t> poly_t;


typedef std::deque< poly_t > poly_sto;

typedef std::vector< poly_t *> poly_ref;

typedef std::map< mono_t , poly_t *> ht_poly_ref;

typedef std::multimap< mono_t , poly_t *> mht_poly_ref;




#endif

