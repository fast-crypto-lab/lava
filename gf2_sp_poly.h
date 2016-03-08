#ifndef _GF2_SP_POLY_H_
#define _GF2_SP_POLY_H_

#include "gf.h"
#include "monomial2.h"
#include "bit_sp_vec.h"




template <unsigned gf_size,unsigned n_mono,unsigned max_deg>
struct sp_poly_t
{



};




#include <vector>
#include <queue>
#include <algorithm>


#define _DEPRECATE_


template <unsigned n_mono,unsigned max_deg>
struct poly2_t
{
//////  typedef and static member /////////////////////////
	typedef poly2_t poly_t;
	typedef mono_2<n_mono,max_deg> mono_t;
	typedef gf<2> gf_t;

	typedef std::vector<mono_t> vec_t;
	static const mono_t zero;
//////////////////////////////////////////

	vec_t vv;

///////////////////////////////

	bool is_zero() const { return vv.empty(); }
	const mono_t & head_term() const { if(is_zero()) return zero; else return vv.front(); }
	uint64_t head_idx() const { return head_term().fast_idx(); }
	unsigned degree() const { return head_term().degree(); }

	const poly_t & set_zero() { vv.clear(); return *this; }
	
/////////////////////////

	unsigned n_terms() const { return vv.size(); }
	const mono_t & operator[](unsigned i) const { return vv[i]; }

	const poly_t& operator = ( const poly_t& p ) {
		vv.resize(p.vv.size());
		vv = p.vv;
		return *this;
	}

/////////////////////////

	void split( poly_t & ph , poly_t & pd , const mono_t & m ) const {
		ph.clear();
		pd.clear();
		unsigned i=0;
		for(;i<vv.size();i++) {
			if( m < vv[i] ) break;
			ph.vv.push_back(vv[i]);
		}
		for(;i<vv.size();i++) pd.vv.push_back(vv[i]);
	}

	bool is_mutiple_of(const mono_t & m) const {
		for(unsigned i=0;i<vv.size();i++)
			if(!vv[i].is_divisible_by(m)) return false;
		return true;
	}

#if 0
	bool CanGenFieldSpoly(mono_t & m) const {
//magma style : yes we can!
//		return true;

		for(int i=0; i< vv.size(); i++){
			if(!vv[i].is_divisible_by(m)){
				return true;
			}
		}
		//std::cout<< "cannot into field equation "<< m <<std::endl;
		return false;

	}

//since Lots of -= will take some computation, I decide to saperate it
	poly_t GenFieldSpoly(mono_t & m) const {
		poly_t p;
		poly_t p2;
			for(int i=0; i< vv.size(); i++){
				if(!vv[i].is_divisible_by(m)){
					p.vv.push_back(vv[i]);
				}else{
					p2.vv.push_back(vv[i]);
				}
			}
		p *= m;
		p -= p2;

		return p;
	}
#endif

	const poly_t & operator -=( const poly_t & a )
	{
		if( a.is_zero() ) return *this;
		if( is_zero() ) { *this = a; return *this; }
		vec_t np;
		np.reserve( a.vv.size() + vv.size() );
		unsigned idxa=0,idx=0;

		while(idxa!=a.vv.size()&&idx!=vv.size()) {
			if( vv[idx] == a.vv[idxa] ) {
				idx++;
				idxa++;
			} else if ( vv[idx] > a.vv[idxa] ) {
				np.push_back( vv[idx] );
				idx++;
			} else {
				np.push_back( a.vv[idxa] );
				idxa++;
			}
		}
		while(idxa!=a.vv.size()){
			np.push_back(a.vv[idxa]);
			idxa++;
		}
		while(idx!=vv.size()){
			np.push_back(vv[idx]);
			idx++;
		}
		vv.swap(np);
		return *this;
	}


	static inline int qdeg( const std::queue<mono_t> qq ) { if(qq.empty()) return -1; else return qq.front().degree(); }

	static void ext_deg1( vec_t & p , const mono_t & m_deg1 ) {
		if( 0 == p.size() ) return;
		std::queue<mono_t> tmp;
		unsigned r_idx = 0;
		unsigned w_idx = 0;

		while( r_idx < p.size() ) {
			mono_t check = p[r_idx++];
			unsigned r_deg = check.degree();
//std::cout << "get " << check << "\n";
			while( !tmp.empty() ){
				if( tmp.front().degree() <= r_deg + 1 ) break;
				p[w_idx++] = tmp.front();
				tmp.pop();
			}
			mono_t mm = gcd( check , m_deg1 );
			if( mm.is_zero() ) { /// do extend
				mono_t attempt = check*m_deg1;
				while( ! tmp.empty() ) {
					if( tmp.front() > attempt ) { /// may be slow
						p[w_idx++] = tmp.front();
						tmp.pop();
					} else break;
				}
				if( tmp.empty() ) {
					p[w_idx++] = attempt;
				} else {
					if( tmp.front() == attempt ) {
						tmp.pop();   /// no write
					} else { /// tmp.front < attempt
						p[w_idx++] = attempt;
					}
				}
			} else { /// no need to extend
				tmp.push( check );
			}
		}
		while( !tmp.empty() ) { p[w_idx++] = tmp.front(); tmp.pop(); }
		p.erase( p.begin()+w_idx , p.end() );
	}

	const poly_t & operator *=( const mono_t & m ) { *this = ext(m); return *this; }
	const poly_t ext( const mono_t & m ) const {
		poly_t r = (*this);
		if(is_zero()) return r;

		mono_t m_deg1 = m.get_deg1mono();
		mono_t rem = m;
//		std::cout << "rem: " << m << "\n";
		while( !rem.is_zero()) {
			ext_deg1( r.vv , m_deg1 );
			rem = rem/m_deg1;
			m_deg1 = rem.get_deg1mono();
//		std::cout << "rem: " << m << " is_zero():" << rem.is_zero() << "\n";
		}
		return r;
	}

//////////////////////////

	int get_position(const mono_t & m) const {
		if( is_zero() ) return -1;
		int st = 0;
		int ed = n_terms()-1;
		while( st <= ed ) {
			int mid = (st+ed)/2;
			if( m == vv[mid] ) return mid;
			else if( m > vv[mid] ) ed = mid-1;
			else st = mid +1; // vv[mid] > m
		}
		return -1;
	}

	gf_t get_coef(const mono_t & m) const {
		int pos = get_position(m);
		if( -1 == pos ) return gf_t(0);
		return gf_t(1);
	}

	bool set_coef(const mono_t & m)
	{
		if( vv.empty() ) { vv.push_back(m); return true; }
		if( vv.back() > m ) { vv.push_back(m); return true; }
		if( vv.back() == m || vv.front() == m ) return true;
		if( vv.front() < m ) { vv.insert( vv.begin(),m); return true; }

		int st = 0;
		int ed = vv.size()-1;
		while( st <= ed ) {
			int mid = (st+ed)/2;
			if( m == vv[mid] ) return true;
			else if( m > vv[mid] ) ed = mid-1;
			else st = mid +1; // vv[mid] > m
		}
		if( m > vv[ed] ) { vv.insert( vv.begin()+ed, m); return true; }
		if( m < vv[st] ) { vv.insert( vv.begin()+st+1, m ); return true; }
		if( m > vv[st] ) { vv.insert( vv.begin()+st, m ); return true; }
		return true;
	}

///////////////////////////////////////////////

	const gf_t eval( const gf_t * vars ) const {
		if(is_zero()) return gf_t(0);

		gf_t r(0);
		for(unsigned i=0;i<vv.size();i++) r += vv[i].eval(vars);
		return r;
	}

	const poly_t eval( const mono_t & m , uint8_t val ) const {
		poly_t r;
		if(is_zero()) { return r; }

		for(unsigned i=0;i<vv.size();i++)
			if( !vv[i].is_divisible_by(m) ) r.vv.push_back( vv[i] );

		if( 0 != val ) { /// val == 1
			for(unsigned i=0;i<vv.size();i++)
				if( vv[i].is_divisible_by(m) ) {
					poly_t tt; tt.set_coef( vv[i]/m );
					r -= tt;
				}
		}
		return r;
	}

	static const poly_t rand_quad( const gf_t * vars ){
		poly_t r;
		for(uint64_t idx=n_mono_le_deg<2,n_mono>::v-1; idx>0; idx--) if(gf_t(1) == gf_t::rand()) r.set_coef( mono_t(idx) );
		if( gf_t(1) == r.eval( vars ) ) r.set_coef( mono_t(0) );
		return r;
	}

//////////////////// for debug ////////////////
	bool santy_check() const {
		for(unsigned i=0;i<vv.size()-1;i++) if(vv[i] <= vv[i+1]) { std::cout << vv[i] << "!>" << vv[i+1]; return false; }
		return true;
	}
	friend std::ostream & operator << (std::ostream & os,const poly_t & a ){
		if( a.is_zero() ) os << "[[0]] 0";
		else if( a.head_term().degree() >= 2 && a.n_terms() > 3 ) {
			os << a.head_term() << " +" <<a.n_terms()-1 << " terms ";
		} else {
			os << a[0];
			for(unsigned i=1;i<a.n_terms();i++) os << "+" << a[i];
		}
		return os;
	}

	static const poly_t magma_in( std::istream & is ) {
		poly_t r;
		mono_t m;
		int mr;
		int c = is.peek();
		char tmp[256];
		while( !is.eof() ) {
			mr = mono_t::magma_in( m , is );

			if( 0 == mr ) {
				poly_t tmp; tmp.set_coef(m);
				r -= tmp;
				c = is.peek();
			}
			else c = is.get();
			if( ',' == c ) { is.getline(tmp,256); break; }
		}

		return r;
	}

	void magma_out( std::ostream & os ) const {
		if(is_zero()) { os << "0,\n"; return; }

		head_term().magma_out( os );

		for( unsigned i=1;i<vv.size();i++ ) {
			os << " + ";
			vv[i].magma_out( os );
		}
		os << ",\n";
	}
};

template <unsigned n_mono,unsigned max_deg>
const mono_2<n_mono,max_deg> poly2_t<n_mono,max_deg>::zero = mono_2<n_mono,max_deg>();



template <unsigned n_mono,unsigned max_deg>
struct dense_vec_t
{
	typedef uint64_t sto_t;
	typedef std::vector<sto_t> vec_t;
	static const unsigned n_bit_sto = sizeof(sto_t)*8;
	static const unsigned log_n_bit = 6;
	static const unsigned mod = 63;

	typedef poly2_t<n_mono,max_deg> poly_t;

/////////////////////
	vec_t vv;
	unsigned st_pos;
	unsigned length;

//////////////////////

	dense_vec_t ( unsigned len ) : vv( (len+mod)>>log_n_bit , 0 ), st_pos(len), length(len) {}

	dense_vec_t ( const poly_t & a , const poly_t & mono_list ) : vv( (mono_list.n_terms()+mod)>>log_n_bit,0), length(mono_list.n_terms()) {
		if( a.is_zero() ) { st_pos = mono_list.n_terms(); return; }
		int pos_i = mono_list.get_position( a.head_term() );
		if( 0 > pos_i ) { st_pos=mono_list.n_terms(); return; }
		else st_pos = pos_i;
		set_idx(st_pos);

		unsigned i = st_pos+1;
		unsigned j=1;
		while( i<mono_list.n_terms() && j<a.n_terms() ) {
			if( a[j] == mono_list[i] ) {
				set_idx(i);
				i++; j++;
			} else i++;  // a[i] < mono_list[j]
		}
	}
	bool is_zero() const { return (st_pos>=length);}

	unsigned head_idx() const { return st_pos; }

	void to_poly( poly_t & a , const poly_t & mono_list ) const {
		a.set_zero();
		if( st_pos >= mono_list.n_terms() ) return;

		for(unsigned i=(st_pos>>log_n_bit);i<vv.size();i++) {
			if(0==vv[i]) continue;
			uint64_t v = vv[i];
			while( v ) {
				unsigned ii = __builtin_ffsll(v)-1;
				a.vv.push_back(mono_list[(i<<log_n_bit)+ii]);
				v &= v-1;
			}
		}
		return;
//		a.vv.push_back(mono_list[st_pos]);
//		for(unsigned i=st_pos+1;i<mono_list.n_terms();i++)
//			if(get_idx(i)) a.vv.push_back(mono_list[i]);
//		return;
	}

	uint8_t set_idx( unsigned idx ) { uint64_t v=1; v<<=(idx&mod); vv[idx>>log_n_bit] |= v; return 1; }

	uint8_t get_idx( unsigned idx ) const {
		uint64_t v=1; v<<=(idx&mod);
		v &= vv[idx>>log_n_bit];
		return (0==v)?0:1;
	}

	const dense_vec_t& operator ^=( const dense_vec_t & a ) {
		for(unsigned i= (a.st_pos>>log_n_bit);i<a.vv.size();i++) vv[i] ^= a.vv[i];
		if( a.st_pos < st_pos ) st_pos = a.st_pos;
		else if( a.st_pos == st_pos ) {
			st_pos = (vv.size())<<log_n_bit;
			for(unsigned i=(a.st_pos>>log_n_bit);i<vv.size();i++) {
				if( 0 != vv[i] ) {
					st_pos = (i<<log_n_bit)+ __builtin_ffsll(vv[i])-1;
					break;
				}
			}
		}

		return *this;
	}
};




#include <list>

template <unsigned n_mono,unsigned max_deg>
struct list_poly
{
	typedef poly2_t<n_mono,max_deg> vpoly_t;
	typedef mono_2<n_mono,max_deg> mono_t;
	typedef gf<2> gf_t;

	typedef std::vector<mono_t> vec_t;
	typedef std::list<mono_t> list_t;

//////////////////////
	list_t lv;

	list_poly() {}

	list_poly( const vpoly_t & p ) { *this = p; }

	void to_vpoly( vpoly_t & p ) const {
		p.vv.clear();
		p.vv.reserve( lv.size() );
		for(typename list_t::const_iterator cit=lv.begin();cit!=lv.end();++cit)
			p.vv.push_back(*cit);
	}

	const list_poly & operator=( const vpoly_t &p ) {
		lv.clear();
		for(unsigned i=0;i<p.n_terms();i++) lv.push_back(p[i]);
		return *this;
	}

	const typename list_t::const_iterator find( const mono_t & m ) const {
		for( typename list_t::const_iterator cit=lv.begin();cit!=lv.end();++cit) {
			if( m == *cit ) return cit;
		}
		return lv.end();
	}

	const list_poly & operator |=( const vpoly_t & a ) { return _or( lv.begin() , a ); }

	typedef typename list_t::iterator iterator;
	const iterator begin() { return lv.begin(); }
	const iterator end() { return lv.end(); }

	const list_poly & _or( iterator st , const vpoly_t & a ) {
		if( a.is_zero() ) return *this;
		if( lv.empty() ) { *this = a; return *this; }

		//typename list_t::iterator it = lv.begin();
		iterator it = st;
		unsigned i=0;

		while( i < a.n_terms() && it != lv.end() ) {
			if( *it == a[i] ) { i++; ++it; }
			else if( *it < a[i] ) lv.insert( it , a[i++] );
			else ++it; /// *it > a[i]
		}
		if( it == lv.end() && i < a.n_terms() ) {
			for(;i<a.n_terms();i++) lv.push_back(a[i]);
		}
		return *this;
	}

};








#endif //#ifndef _GF2_SP_POLY_H_

