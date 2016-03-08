#ifndef _BIT_SP_VEC_H_
#define _BIT_SP_VEC_H_

#include <stdint.h>
#include <iostream>

//#ifdef NO_SSE
#if 1

#include <emmintrin.h>




inline uint16_t get_zero_mask( __m128i v )
{
	return _mm_movemask_epi8( _mm_cmpeq_epi8(v,_mm_setzero_si128()) );
}

inline uint16_t get_nonzero_mask( __m128i v )
{
	return (0xffff ^ _mm_movemask_epi8( _mm_cmpeq_epi8(v,_mm_setzero_si128()) ) );
}




static const __m128i _bit_mask_128[128] = {
_mm_set_epi32(0,0,0,0x1),
_mm_set_epi32(0,0,0,0x3),
_mm_set_epi32(0,0,0,0x7),
_mm_set_epi32(0,0,0,0xf),
_mm_set_epi32(0,0,0,0x1f),
_mm_set_epi32(0,0,0,0x3f),
_mm_set_epi32(0,0,0,0x7f),
_mm_set_epi32(0,0,0,0xff),
_mm_set_epi32(0,0,0,0x1ff),
_mm_set_epi32(0,0,0,0x3ff),
_mm_set_epi32(0,0,0,0x7ff),
_mm_set_epi32(0,0,0,0xfff),
_mm_set_epi32(0,0,0,0x1fff),
_mm_set_epi32(0,0,0,0x3fff),
_mm_set_epi32(0,0,0,0x7fff),
_mm_set_epi32(0,0,0,0xffff),
_mm_set_epi32(0,0,0,0x1ffff),
_mm_set_epi32(0,0,0,0x3ffff),
_mm_set_epi32(0,0,0,0x7ffff),
_mm_set_epi32(0,0,0,0xfffff),
_mm_set_epi32(0,0,0,0x1fffff),
_mm_set_epi32(0,0,0,0x3fffff),
_mm_set_epi32(0,0,0,0x7fffff),
_mm_set_epi32(0,0,0,0xffffff),
_mm_set_epi32(0,0,0,0x1ffffff),
_mm_set_epi32(0,0,0,0x3ffffff),
_mm_set_epi32(0,0,0,0x7ffffff),
_mm_set_epi32(0,0,0,0xfffffff),
_mm_set_epi32(0,0,0,0x1fffffff),
_mm_set_epi32(0,0,0,0x3fffffff),
_mm_set_epi32(0,0,0,0x7fffffff),
_mm_set_epi32(0,0,0,0xffffffff),
_mm_set_epi32(0,0,0x1,0xffffffff),
_mm_set_epi32(0,0,0x3,0xffffffff),
_mm_set_epi32(0,0,0x7,0xffffffff),
_mm_set_epi32(0,0,0xf,0xffffffff),
_mm_set_epi32(0,0,0x1f,0xffffffff),
_mm_set_epi32(0,0,0x3f,0xffffffff),
_mm_set_epi32(0,0,0x7f,0xffffffff),
_mm_set_epi32(0,0,0xff,0xffffffff),
_mm_set_epi32(0,0,0x1ff,0xffffffff),
_mm_set_epi32(0,0,0x3ff,0xffffffff),
_mm_set_epi32(0,0,0x7ff,0xffffffff),
_mm_set_epi32(0,0,0xfff,0xffffffff),
_mm_set_epi32(0,0,0x1fff,0xffffffff),
_mm_set_epi32(0,0,0x3fff,0xffffffff),
_mm_set_epi32(0,0,0x7fff,0xffffffff),
_mm_set_epi32(0,0,0xffff,0xffffffff),
_mm_set_epi32(0,0,0x1ffff,0xffffffff),
_mm_set_epi32(0,0,0x3ffff,0xffffffff),
_mm_set_epi32(0,0,0x7ffff,0xffffffff),
_mm_set_epi32(0,0,0xfffff,0xffffffff),
_mm_set_epi32(0,0,0x1fffff,0xffffffff),
_mm_set_epi32(0,0,0x3fffff,0xffffffff),
_mm_set_epi32(0,0,0x7fffff,0xffffffff),
_mm_set_epi32(0,0,0xffffff,0xffffffff),
_mm_set_epi32(0,0,0x1ffffff,0xffffffff),
_mm_set_epi32(0,0,0x3ffffff,0xffffffff),
_mm_set_epi32(0,0,0x7ffffff,0xffffffff),
_mm_set_epi32(0,0,0xfffffff,0xffffffff),
_mm_set_epi32(0,0,0x1fffffff,0xffffffff),
_mm_set_epi32(0,0,0x3fffffff,0xffffffff),
_mm_set_epi32(0,0,0x7fffffff,0xffffffff),
_mm_set_epi32(0,0,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xf,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1f,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3f,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7f,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1ff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3ff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7ff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xfff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1fff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3fff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7fff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1ffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3ffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7ffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xfffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1fffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3fffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7fffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1ffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3ffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7ffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xfffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x1fffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x3fffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0x7fffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xf,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1f,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3f,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7f,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1ff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3ff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7ff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xfff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1fff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3fff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7fff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1ffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3ffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7ffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xfffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1fffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3fffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7fffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1ffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3ffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7ffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xfffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x1fffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x3fffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0x7fffffff,0xffffffff,0xffffffff,0xffffffff),
_mm_set_epi32(0xffffffff,0xffffffff,0xffffffff,0xffffffff)
};


inline int leftmost_bit_8( uint8_t msk )
{
	return 31-__builtin_clz((uint32_t)msk);
}

inline int leftmost_bit_16( uint16_t msk )
{
	if( 0 == msk ) return -1;
	return 31-__builtin_clz((uint32_t)msk);
}


inline int leftmost_bit_128( __m128i vv )
{
	union{
	__m128i v128;
	uint8_t v8[16];
	};
	v128 = vv;

	int r = leftmost_bit_16( get_nonzero_mask( vv ) );
	if( 0 > r ) return -1;
	return (r<<3)+leftmost_bit_8(v8[r]);
}


inline uint16_t bit_mask_16( unsigned n_bit )
{
	return 0xffff>>(16-n_bit);
}

inline int leftmost_bit_less_idx_16( uint16_t msk , unsigned idx )
{
	msk &= bit_mask_16( idx );
	return leftmost_bit_16( msk );
}

inline uint8_t rightmost_bit_8(uint16_t msk)
{
	return __builtin_ctz((uint32_t)msk);
}

inline int rightmost_bit_16( uint16_t msk )
{
//	if( 0 == msk ) return -1;
//	return __builtin_ctz((uint32_t)msk);
	return __builtin_ffs((uint32_t)msk)-1;
}


inline int rightmost_bit_128( __m128i vv )
{
	union{
	__m128i v128;
	uint8_t v8[16];
	};
	v128 = vv;

	int r = rightmost_bit_16( get_nonzero_mask( vv ) );
	if( 0 > r ) return -1;
	return (r<<3)+rightmost_bit_8(v8[r]);
}






//////////////////////////////////////////////////////////////////////////



struct chunk
{
	uint32_t base_idx;
	union{
	__m128i vv;
	uint8_t v[16];
	};

	chunk(): base_idx(0) { vv = _mm_setzero_si128(); }

	bool is_zero() const { return 0xffff == get_zero_mask(vv); }
	int leftmost_bit() const { return leftmost_bit_128(vv); }
	int rightmost_bit() const{ return rightmost_bit_128(vv); }
	int leftmost_bit_le_idx(uint8_t idx) const { return leftmost_bit_128(vv&_bit_mask_128[idx]); }

	bool set_idx( uint8_t idx ) { v[idx>>3]|=(1<<(idx&7)); return true; }
	bool get_idx( uint8_t idx ) const { return v[idx>>3]&(1<<(idx&7)); }

	bool operator==( const chunk & a ) const {return base_idx==a.base_idx;}
	bool operator<( const chunk & a ) const {return base_idx<a.base_idx;}

	bool offset( chunk & hi , chunk & lo , uint64_t offset ) const {
		unsigned off = offset & 127;
		unsigned base_off = offset >> 7;
		if( 0 == off ) { lo = (*this); lo.base_idx += base_off; return true;  }
		__m128i oh_in_l = _mm_srli_si128( vv , 8 );
		__m128i ol_in_h = _mm_slli_si128( vv , 8 );
		if( off > 64 ) {
			lo.vv = _mm_slli_epi64(ol_in_h,off-64);
			lo.base_idx = base_idx+base_off;
			hi.vv = _mm_srli_epi64(vv,64-off)|_mm_slli_epi64(oh_in_l,off-64);
			hi.base_idx = lo.base_idx+1;
		} else if ( off < 64 ) {
			lo.vv = _mm_slli_epi64(vv,off)|_mm_srli_epi64(ol_in_h,64-off);
			lo.base_idx = base_idx+base_off;
			hi.vv = _mm_srli_epi64(oh_in_l,64-off);
			hi.base_idx = lo.base_idx+1;
		} else { /// 64 == off
			lo.vv = ol_in_h;
			lo.base_idx = base_idx + base_off;
			hi.vv = oh_in_l;
			hi.base_idx = lo.base_idx+1;
		}
		return true;
	}

	friend std::ostream & operator <<( std::ostream & os , const chunk & a ) {
//	        os.setf(std::ios::hex);
		os << "[";
		for(unsigned i=16;i--;) os <<((uint32_t)a.v[i])<<" ";
		os << "]_" << a.base_idx;
//		os.setf(std::ios::dec);
		return os;
	}


};




#else // #ifdef NO_SSE
#endif // #ifdef NO_SSE






/////////////////////////////////////////////////////////////////////////////////////////



#include <list>

struct bit_sp_vec
{
	typedef bit_sp_vec bvec;
	typedef std::list<chunk> list_t;
	typedef std::list<chunk>::iterator it_t;
	typedef std::list<chunk>::const_iterator cit_t;

	typedef std::list<chunk>::reverse_iterator rit_t;
	typedef std::list<chunk>::const_reverse_iterator crit_t;

	int64_t head_idx;
	list_t v;

	bit_sp_vec(): head_idx(-1) {}
	//~bit_sp_vec() { }
	bit_sp_vec( const bit_sp_vec& a ): head_idx(a.head_idx) {
		cit_t it2 = a.v.begin();
		while( a.v.end() != it2 ) { if( !(*it2).is_zero() ) v.insert(v.end(),*it2); it2++; }
	}

	bool is_zero() const { return -1==head_idx; }

	int64_t get_head_idx() const { return head_idx; }


	bvec & operator ^= (const bvec & a) {
		if(a.is_zero()) return *this;
		if(this->is_zero()) { *this=a; return *this; }

		it_t it1=v.begin();
		cit_t it2=a.v.begin();
		while( it1!=v.end() && it2 !=a.v.end() ) {
			/// no erase  chunk !!!!!!?????
			/// leave the erasing work in constructor !!!!!
			if((*it1)==(*it2)) { (*it1).vv ^= (*it2).vv; it1++; it2++; }
			else if((*it1)<(*it2)) { if(!(*it2).is_zero()) v.insert(it1,*it2); it2++; }
			else { it1++; }
		}
		/// it2 != end -->  it1 == end
		while( it2 != a.v.end()) { if(!(*it2).is_zero()) v.insert(v.end(),*it2); it2++; }
		//if( it2 != a.v.end() ) v.insert(v.end(),it2,a.v.end());

		flip_idx( a.head_idx );

		return *this;
	}
	const bvec operator^(const bvec & a ) const { bvec r=(*this); r^=a; return r; }


	bvec & operator |= (const bvec & a) {
		if(a.is_zero()) return *this;
		if(this->is_zero()) { *this=a; return *this; }
		if(a.head_idx>head_idx) { uint64_t oh=head_idx; head_idx=a.head_idx; set_idx(oh); }
		else set_idx(a.head_idx);
		it_t it1=v.begin();
		cit_t it2=a.v.begin();
		while( it1!=v.end() && it2 !=a.v.end() ) {
			if((*it1)==(*it2)) { (*it1).vv |= (*it2).vv; it1++; it2++; }
			else if((*it1)<(*it2)) { if(!(*it2).is_zero()) v.insert(it1,*it2); it2++; }
			else { it1++; }
		}
		/// it2 != end -->  it1 == end
		while( it2 != a.v.end()) { if(!(*it2).is_zero()) v.insert(v.end(),*it2); it2++; }
		//if( it2 != a.v.end() ) v.insert(v.end(),it2,a.v.end());
		return *this;
	}
	const bvec operator|(const bvec & a ) const { bvec r=(*this); r|=a; return r; }


	/// check list only, privately use
	int64_t __next_bit(uint64_t idx) const {
		idx -= 1;
		uint32_t bidx = (idx>>7);
		uint8_t rem = (idx&127);
		for(cit_t it=v.begin();it!=v.end();it++) {
			if(bidx==(*it).base_idx) { if( (!(*it).is_zero())&&((*it).rightmost_bit()<=rem) ) return (((uint64_t)bidx)<<7)+(*it).leftmost_bit_le_idx(rem);   }
			if(bidx>(*it).base_idx) if(!(*it).is_zero()) return (((uint64_t)(*it).base_idx)<<7)+(*it).leftmost_bit();
		}
		return -1;
	}
	/// idx would not be return value except 0
	uint64_t next_bit(uint64_t idx) const {
		if(0>head_idx) return 0;
		if((int64_t)idx>head_idx) return head_idx;
		if(0==idx) return 0;

		int64_t r=__next_bit(idx);
		return (r>0)?r:0;
	}


	/// list only operation
	bool __set_idx(uint64_t idx) {
		chunk i;
		i.base_idx = idx>>7;
		i.set_idx( idx&127 );
		for(it_t it=v.begin();it!=v.end();it++) {
			if(i==(*it)) { (*it).vv |= i.vv; return true; }
			if((*it)<i) { v.insert(it,i); return true; }
		}
		v.insert(v.end(),i);
		return true;
	}
	bool set_idx(uint64_t idx) {
		if( 0 > head_idx ) { head_idx = idx; return true; }
		if( (int64_t)idx > head_idx ) { uint64_t oh=head_idx; head_idx=idx; return set_idx(oh); }
		if( (int64_t)idx == head_idx ) return true;
		return __set_idx(idx);
	}

	int64_t __max_idx() const {
		for( cit_t cit=v.begin(); cit!=v.end(); cit++) {
			int lb = leftmost_bit_128((*cit).vv);
			if( lb >= 0 ) return (((*cit).base_idx)<<7)+lb;
		}
		return -1;
	}

	//bvec & operator ^= (const chunk & a) {
	bool xor_chunk( const chunk & a , uint64_t offset ) {
		__xor_chunk( a , offset );
		if( v.empty() ) return true;
		if( 0 > head_idx  ) {
			int64_t max_list_idx = __max_idx();
			__flip_idx( max_list_idx );
			head_idx = max_list_idx;
			return true;
		}
		if( ((unsigned)(head_idx>>7)) > (a.base_idx+(offset>>7)+1) ) return true;
		int64_t max_list_idx = __max_idx();
		if( max_list_idx < 0 ) return true;
		if( max_list_idx > head_idx ) { /// (include case : 0 > head_idx )
			//if( head_idx >= 0 ) __flip_idx( head_idx );
			__flip_idx( head_idx );
			__flip_idx( max_list_idx );
			head_idx = max_list_idx;
		} else if( max_list_idx == head_idx ) {
			__flip_idx( max_list_idx );
			head_idx = __max_idx();
			if( head_idx > 0 ) __flip_idx( head_idx );
		}
		return true;
	}

	/// list only operation
	bool __xor_chunk(const chunk & a , uint64_t offset ) {
		chunk h,l;
		a.offset( h , l , offset );
		if( !h.is_zero()) __xor_chunk(h);
		if( !l.is_zero()) __xor_chunk(l);
		return true;
	}

	bool xor_chunk( const chunk & a ) {
		__xor_chunk( a );
		if( 0 > head_idx  ) {
			int64_t max_list_idx = __max_idx();
			__flip_idx( max_list_idx );
			head_idx = max_list_idx;
			return true;
		}
		if( (unsigned)(head_idx>>7) > a.base_idx ) return true;
		int64_t max_list_idx = __max_idx();
		if( max_list_idx > head_idx ) {
			__flip_idx( head_idx );
			__flip_idx( max_list_idx );
			head_idx = max_list_idx;
		} else if( max_list_idx == head_idx ) {
			__flip_idx( max_list_idx );
			head_idx = __max_idx();
			if( head_idx > 0 ) __flip_idx( head_idx );
		}
		return true;
	}

	bool __xor_chunk(const chunk & a ) {
		for(it_t it=v.end();it!=v.begin();) {
			it--;
			if(a==(*it)) { (*it).vv^=a.vv; if((*it).is_zero()) v.erase(it); return true; }
			if(a<(*it)) { it++; v.insert(it,a); return true; }
		}
		v.push_front(a);
		return true;
	}

	/// list only operation
	bool __flip_idx(uint64_t idx) {
		chunk ins;
		ins.base_idx=idx>>7;
		ins.set_idx(idx&127);
		for(it_t it1=v.begin();it1!=v.end();it1++) {
			if(ins==(*it1)) { (*it1).vv ^= ins.vv; if((*it1).is_zero()) v.erase(it1); return true;} ///// prune
			if((*it1)<ins) { v.insert(it1,ins); return true; }
		}
		v.insert(v.end(),ins);
		return true;
	}
	bool flip_idx(uint64_t idx) {
		if( 0 > head_idx ) { head_idx = idx; return true; }
		if( (int64_t)idx > head_idx ) { uint64_t oh=head_idx; head_idx=idx; return __flip_idx(oh); }
		if(head_idx==(int64_t)idx){
			if(0==head_idx) {head_idx=-1; return true; }
			head_idx=__next_bit(head_idx); ////////////////////////////////////////
			if(-1==head_idx) { v.clear(); return true; } //// prune
			__flip_idx(head_idx);
			return true;
		}
		__flip_idx(idx);
		return true;
	}


	bool get_idx(uint64_t idx) const {
		if((int64_t)idx==head_idx) return true;
		uint32_t bidx=(idx>>7);
		uint8_t rem=(idx&127);
		for(cit_t it=v.begin();it!=v.end();it++) {
			if(bidx==(*it).base_idx) return (*it).get_idx(rem);
			if(bidx>(*it).base_idx) return false;
		}
		return false;
	}

	friend std::ostream & operator <<( std::ostream & os , const bvec & a ) {
		std::cout << a.head_idx << ",";
		for(cit_t it=a.v.begin();it!=a.v.end();it++) std::cout << (*it) << ",";
		return os;
	}

};



#endif  // _BIT_SP_VEC_H_
