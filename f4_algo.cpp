

#include "f4_algo.h"
using namespace std;


template struct n_mono_le_deg_table<N_VAR,MAX_DEG>;
template struct min_mono_tab_init<N_VAR,MAX_DEG>;

#include <algorithm>
#include <m4ri/m4ri.h>
#include <unordered_map>

static int total_pair = 0;
static int total_elim_size = 0;
static int total_reductor = 0;
static int total_reductee = 0;
static int max_elim_size = 0;
static int iteration_count = 0;
static double genpair_time = 0;
static double genspoly_time = 0;
static double findreductor_time = 0;
static double find_divisor_time = 0;
static double conversion_time = 0;
static double computation_time = 0;
static double addpoly_time = 0;

#include <sys/time.h>
static double get_ms_time(void) {
	struct timeval timev;

	gettimeofday(&timev, NULL);
	return (double) timev.tv_sec * 1000 + (double) timev.tv_usec / 1000;
}

std::vector<std::vector<int> > manual_pair_list;

mono_t X[N_VAR+1];
extern int selector[8];
extern int selector_count;

 
struct MonomialHash {
 std::size_t operator()(const mono_t& k) const
 {
     return k.vv;
 }
};
 
struct MonomialEqual {
 bool operator()(const mono_t& lhs, const mono_t& rhs) const
 {
    return lhs == rhs;
 }
};

struct pair_t
{

//common
	mono_t lcm;
	bool is_field_pair;
	int deg;

//field pair only
	mono_t field_pair_mono;
	int field_pair_from;

//normal pair only
	int normal_pair_from_1;
	int normal_pair_from_2;

	pair_t(){
		is_field_pair = false;
		field_pair_from = -1;
		normal_pair_from_1 = -1;
		normal_pair_from_2 = -1;
	}
	bool is_generated_from(int i, int j) const {
		if(i == normal_pair_from_1 && j == normal_pair_from_2)
			return true;
		if(i == normal_pair_from_2 && j == normal_pair_from_1)
			return true;
		return false;
	}

	inline void set_pair_origin(int i, int j) {
		normal_pair_from_1 = i;
		normal_pair_from_2 = j;
	}

	friend std::ostream & operator << ( std::ostream & os , const pair_t & pp ) {
		if(pp.is_field_pair)
			os << "Field pair :(" << pp.field_pair_from <<" , "<< pp.field_pair_mono <<" ) ";
		else
			os << "Normal pair :( " << pp.normal_pair_from_1 <<" , " << pp.normal_pair_from_2 << " )"; 
		return os;
	}
};



struct extension
{
	unsigned idx;
	mono_t t;

	bool operator==( const extension & a ) const { return (idx==a.idx)&&(t==a.t); }
	bool operator!=( const extension & a ) const { return !(*this==a); }

	friend bool operator > ( const extension & a , const extension & b )
	{
		if( a.idx == b.idx) return a.t > b.t;
		return a.idx > b.idx;
	}
	friend bool operator < ( const extension & a , const extension & b ) { return b>a; }
	friend bool operator >= ( const extension & a , const extension & b ) { return !(a<b); }
	friend bool operator <= ( const extension & a , const extension & b ) { return !(a>b); }
};



template <typename ctnr>
void dump_g( const ctnr & g )
{
	std::cout << "[G]: " << g.size() << " polys.\n";

	unsigned i=0;
	for( typename ctnr::const_iterator cit=g.begin();cit!=g.end();cit++)
		std::cout << "[" << i++ << "]:" << *cit->second << "\n";
}

template <typename ctnr>
void dump_g1( const ctnr & g )
{
	std::cout << "[G]: " << g.size() << " polys.\n";
	unsigned i=0;
	for( typename ctnr::const_iterator cit=g.begin();cit!=g.end();cit++)
		std::cout << "[" << i++ << "]:" << **cit << "\n";
}



struct map_mono_t{
	typedef map<mono_t,int> map_t;
	map_t lv;
	void add_poly(const base_poly_t & a ) {
		for(unsigned i=0;i<a.n_terms();i++) 
			lv.insert(pair<mono_t,int>(a[i],0));
	}
	unsigned size() const { return lv.size(); };
	struct iterator{
		iterator(map_t::reverse_iterator it):_it(it){}
		map_t::reverse_iterator _it;
		iterator& operator++(int dummy){
			++_it;
			return *this;
		}
		const mono_t& operator*(){
			return _it->first;
		}
		bool operator!=(const iterator& it)const{
			return _it!=it._it;
		}
	};
	iterator begin(){ return iterator(lv.rbegin()); }
	iterator end(){return iterator(lv.rend()); }
};

struct f4_gb
{
	poly_ref polys;
	ht_poly_ref ht_polys;
	std::map<mono_t,unsigned> ht_idx;
#ifdef __HAS_SIG__
	bool in_syz( const Src & s ) const {
		mono_t t;
		s.t.to_mono2(t);
		for(ht_poly_ref::const_iterator it=ht_polys.begin();it!=ht_polys.end();it++) {
			if( it->first > t ) break;
			if( ! t.is_divisible_by( it->first ) ) continue;
			mono_t ext = t / it->first;
			Src s2 = it->second->s;
			s2 *= ext;
			if( s2 < s ) return true;
		}
		return false;
	}
#endif
	unsigned size() const { return polys.size(); }

	poly_t & operator[](unsigned i) const { return *polys[i]; }

	const mono_t min_ht() const { if(0==size()) return mono_t(); else return ht_polys.begin()->first; }

	int find_ht( const mono_t & m ) const {
		std::map<mono_t,unsigned>::const_iterator it = ht_idx.find(m);
		if( ht_idx.end() == it ) return -1;
		else if (polys[it->second]->unnecessary_for_new_pair) return -1;
		else return it->second;
	}

	int find_divisor( const mono_t & m ) const {
		int r1 = find_ht(m);
		if( r1 >= 0 ) return r1;
		unsigned deg = m.degree();

		if( deg <= 5 ){
			mono_t tmp[N_VAR];
			deg = m.split( tmp );
			for(unsigned k=1;k<=deg;++k){
				for(unsigned i=(1<<deg)-1;i>=1;--i){
					mono_t tmp_mono;
					unsigned counter=0;
					for(int j=0;j<deg;++j){
						if((i>>j)&1){
							++counter;
							tmp_mono = tmp_mono * tmp[j];
							//set_nth_bit_on(tmp_mono.v,tmp[j]);
						}
					}
					if(counter!=k) continue;
					std::map<mono_t,unsigned>::const_iterator it = ht_idx.find(tmp_mono);
					if( ht_idx.end() != it ) 
						return it->second;
					//cout<<tmp_mono<<" ";
				}	
			}
		} else { /// case: deg > 5
			for( std::map<mono_t,unsigned>::const_iterator it=ht_idx.begin();it!=ht_idx.end();it++) {
				if(it->first.degree() >= m.degree() ) continue;
				if( m.is_divisible_by(it->first) ) { return it->second; }
			}
		}
		return -1;
	}

	//I don't want a same ht
	int find_strict_divisor( const mono_t & m , int i ) const {
		for( std::map<mono_t,unsigned>::const_iterator it=ht_idx.begin();it!=ht_idx.end();it++) {
			if(it->first.degree() >= m.degree() ) continue;
			if( m.is_divisible_by(it->first) ) { 
				//if(it->second > i ) continue;
				return it->second; 
			}
		}
		return -1;
	}

//if a new poly with the same term appears, it replace the older one
	bool add_poly( poly_t * p) {
		mono_t ht = p->head_term();
		if( ht_polys.end() != ht_polys.find( ht ) ) {
			DUMP(std::cout << "[GB.add_poly()] repeat ht: " << ht << "\n");
//			return false;
		}

		polys.push_back( p );
		ht_polys[ht] = p;
		ht_idx[ht] = polys.size()-1;

		return true;
	}

};



struct relation_t
{
	std::vector< std::set<unsigned> > relations;
	std::vector< std::set<unsigned> > golden_relations;
	//golden relation is relation between ht divisibles, which can NOT be used in buch 2 
	//(since I'm killing every unnecessary person in spoly, it's naturally killed)
	//maybe remove relation will also work, and cleaner, but I'm afraid of the unawared consequences

	bool is_golden_relation( unsigned i,unsigned j) const {
		if( i >= golden_relations.size() ) return false;
		if(golden_relations[i].find(j) != golden_relations[i].end())	
			return true;
		return false;
	}

	// buch criteria 2
	bool is_unnecessary( unsigned i,unsigned j, const f4_gb & gb, const mono_t & lcm, vector<pair_t>& temp_wl) const {
		if( i >= relations.size() ) return false;
		if( j >= relations.size() ) return false;

//if this is a golden relation, it can't be killed. and cannot kill other people
		if(is_golden_relation(i, j))
			return false;

		std::set<unsigned> int_set;
		std::set_intersection ( relations[i].begin(), relations[i].end()
			, relations[j].begin() , relations[j].end()
			,std::inserter(int_set,int_set.begin()) );
		for(std::set<unsigned>::iterator it=int_set.begin(); it!=int_set.end(); it++){
			int k = *it;
			if(lcm.is_divisible_by(gb[k].head_term())&&
				! is_golden_relation(i, k) && 
				! is_golden_relation(j, k) ){
				return true;
			}
		}

//start killing other only when not killed(will result in a far from optimal result)
		for(std::set<unsigned>::iterator it=int_set.begin(); it!=int_set.end(); it++){
			int k = *it;
			/*
				1. find in wl that is either ik or jk (only one will exist)
				2. check it's lcm can be devided by j  or i
				3. if so, remove it.
				4. actually it's possible that there is another jk or ik in the big pair storage. But to search that, O(wl) is needed, yielding a total of O(wl^2)
				5. or the whole program should be rewritten
			*/
			int t;
			for(t=0; t<temp_wl.size(); t++){
				if (temp_wl[t].is_generated_from(i, k)){
					if(temp_wl[t].lcm.is_divisible_by(gb[j].head_term())){
						temp_wl.erase(temp_wl.begin()+t);
						break;
					}
				}
				if (temp_wl[t].is_generated_from(j, k)){
					if(temp_wl[t].lcm.is_divisible_by(gb[i].head_term())){
						temp_wl.erase(temp_wl.begin()+t);
						break;
					}
				}
			}
		}

		return false;
	}

	bool add_relation( unsigned size , unsigned i,unsigned j) {
		//if( is_unnecessary(size,i,j) ) return false;
		if( size > relations.size() ) relations.resize( size );
		if( i >= relations.size() ) return false;
		if( j >= relations.size() ) return false;
		relations[i].insert(j);
		relations[j].insert(i);
		return true;
	}
	bool add_golden_relation( unsigned size , unsigned i,unsigned j) {
		if( size > golden_relations.size() ) golden_relations.resize( size );
		if( i >= golden_relations.size() ) return false;
		if( j >= golden_relations.size() ) return false;
		golden_relations[i].insert(j);
		golden_relations[j].insert(i);
		return true;
	}
};




using namespace std;



//////////////////////////////////////////////////////////////////////////////////////////////
static unsigned removed_by_syz_check = 0;
static unsigned temp_counter = 0;
static unsigned reduced_to_zero = 0;
bool check_answer(ht_poly_ref& sorter);
bool check_rule(poly_t* p, int n);

static
bool is_cri_pair( pair_t & pp , relation_t & rels, const f4_gb & gb , unsigned i , unsigned j , vector<pair_t>& temp_wl)
{
	mono_t gcd_r = gcd( gb[i].head_term() , gb[j].head_term() );
	mono_t lcm_r = gb[i].head_term() * gb[j].head_term();
	if(gb[i].unnecessary_for_new_pair || gb[j].unnecessary_for_new_pair) { return false; }

	// buch criteria 1
	if( gcd_r.is_zero()  ) {
		//rels.add_relation((i>j)?i:j,i,j); /// no need; have to consider with "unnecessary polynomial"
		return false;
	}

	if(gb[j].head_term().is_divisible_by(gb[i].head_term()) ||
		gb[i].head_term().is_divisible_by(gb[j].head_term())){
		rels.add_golden_relation( gb.size() , i , j );
	} else if( rels.is_unnecessary(i,j, gb, lcm_r, temp_wl) ) { // buch criteria 2
		rels.add_relation( (i>j)?i:j ,i,j);
		return false;
	}

	pp.deg = lcm_r.degree();
	pp.lcm = lcm_r;
	pp.set_pair_origin(i, j);
//	if( pp.lcm.degree() > MAX_SRC_DEG )  {forchecked_large_degree++; return false;}
//	if( pp.p.head_term().degree() > MAX_DEG )  { return false; }

	return true;
}


//typedef std::deque<pair_t> worklist_t;
//typedef std::multimap<mono_t,pair_t> worklist_t;
typedef std::multimap<int,pair_t> worklist_t;

static
void gen_pair( worklist_t & wl , relation_t & rels , const f4_gb & gb , unsigned st , unsigned ed )
{
	genpair_time -= get_ms_time();
	for(unsigned i=st;i<ed;i++) {

		cout << endl;
		cout << "generating pairs" <<endl;
		cout << "polynomial "<< i <<" : " << gb[i] << endl;

		int field_gen_counter = 0;
		int pair_gen_counter = 0;
		int total_gen = wl.size();
		mono_t this_ht = gb[i].head_term();
		//vector<mono_t> varlist = this_ht.split();
		mono_t varlist[N_VAR];
		unsigned varlist_size = this_ht.split(varlist);
		
		vector<pair_t> temp_wl;//this is for a good buchburger's crit 2
		// check if THIS is unnecessary
		int r = gb.find_strict_divisor(this_ht, i);
		if(r != -1){
			pair_t pp;
			rels.add_golden_relation( gb.size() , i , r );
			pp.lcm = gb[i].head_term() * gb[r].head_term();
			pp.deg = pp.lcm.degree();
			pp.set_pair_origin(i, r);
			total_pair++;
			pair_gen_counter++;
			wl.insert( std::pair<int,pair_t>(pp.deg,pp) );
			rels.add_relation( gb.size() , i , r );
			cout << " mark me as unnecessary "<< endl;
			gb[i].mark_unnecessary(gb[r].head_term(), r);
		}else{
//			for(int j=i-1;j>=0;j--) {
//much to my surprise, going from 0 to i yields far less pair than the other way, I don't know why
//the surplus will get elimed later though, so this only affects pair-finding speed
			for(int j=0;j<i;j++) {
				pair_t pp;
				if( ! is_cri_pair( pp , rels , gb , i , j, temp_wl ) ) continue;
			
				if(gb[j].head_term().is_divisible_by(gb[i].head_term())) {
					cout << " mark " <<j << " as unnecessary "<< endl;
					gb[j].mark_unnecessary(gb[i].head_term(), i);
				}
				total_pair++;
				pair_gen_counter++;
				temp_wl.push_back( pp );
				rels.add_relation( gb.size() , i , j );
			}
			for(int j=0; j<temp_wl.size(); j++ )
				wl.insert( std::pair<int,pair_t>(temp_wl[j].deg,temp_wl[j]) );

			//only find field spoly when necessary
			for(int j = 0; j< (int)varlist_size; j++){
				if(varlist_size == 1)//I'm not sure about this
					break;
				if(gb[i].CanGenFieldSpoly(varlist[j])){
					pair_t pp;
					pp.deg = gb[i].degree()+1;
					pp.is_field_pair = true;
					pp.field_pair_from = i;
					pp.field_pair_mono = varlist[j];
					wl.insert( std::pair<int,pair_t>(pp.deg,pp) );
					total_pair++;
					field_gen_counter++;
				}
			}
		}

		total_gen = wl.size()-total_gen;
		cout << total_gen<< " pair generated." << endl;
//		cout << field_gen_counter<< " field equation spolys." << endl;
//		cout << total_gen - field_gen_counter << " normal spolys." << endl;//since some might be later deleted

	}

	genpair_time += get_ms_time();
}



static
bool gen_spoly( mht_poly_ref & new_g , poly_sto & all_polys , const f4_gb & gb , const worklist_t& wl, const relation_t & rels)
{

	if(wl.empty()) return false;
	int total_pair_count = 0;
	int total_field_pair_count = 0;
	int eliminated_field_pair_count = 0;
	int eliminated_normal_pair_count = 0;
	int mindeg = wl.begin()->first;//since this is a multimap, this is true
	//worklist_t::iterator it;
	worklist_t::const_iterator it;

	for(it=wl.begin();it!=wl.end();it++) {
		if(it->first > mindeg){
			break;
		}
		poly_t p;
		if(it->second.is_field_pair){
			total_field_pair_count++;
			if(gb[it->second.field_pair_from].unnecessary_for_new_pair){
				eliminated_field_pair_count++;
				continue;
			}
			p = gb[it->second.field_pair_from].GenFieldSpoly(it->second.field_pair_mono);	

		}else{
			int i = it->second.normal_pair_from_1;
			int j = it->second.normal_pair_from_2;
			//this is ugly, maybe I can use golden relation to improve it
			if(gb[i].unnecessary_for_new_pair || gb[j].unnecessary_for_new_pair){
					if(!rels.is_golden_relation(i, j)){
						eliminated_normal_pair_count++;
						continue;
				}
			}
			p = gb[i];
			poly_t p2 = gb[j];
			p *= it->second.lcm / gb[i].head_term();
			p2 *= it->second.lcm / gb[j].head_term();
			p -= p2;

		}
		all_polys.push_back( p );
		new_g.insert( make_pair(p.head_term(), &all_polys.back()) );
	}
	/// wl.erase(wl.begin(), it);
	
	cout << new_g.size()<<" pairs chosen, degree: " << mindeg <<endl;
	cout << total_field_pair_count << " field pairs, "<< eliminated_field_pair_count <<" eliminated." <<endl;
	cout << new_g.size() + eliminated_normal_pair_count - (total_field_pair_count - eliminated_field_pair_count) << " normal pairs, "<<eliminated_normal_pair_count<<" eliminated." <<endl;
//	cout << wl.size() << " pair left." <<endl;
	return true;
}

static
void guass_elim( poly_ref & g )
{
	for(unsigned i=0;i<g.size();i++) {
		for(unsigned j=0;j<g.size();j++) {
			if(i==j) continue;
			if(! g[j]->get_coef(g[i]->head_term()).is_zero() ) *g[j] -= *g[i];
		}
	}
}

static
bool prepare_reductors(  map_mono_t & list_monos ,poly_ref& selected_reductor,  unordered_map<mono_t, poly_t*, MonomialHash, MonomialEqual> & reductor_finder , poly_sto & reductor_storage
	, const mono_t & min_mono , const f4_gb & gb )
{

	map_mono_t::iterator lit = list_monos.begin();
	while( lit != list_monos.end() ) {
		mono_t m = *lit;
		lit++;
		if(reductor_finder.find(m) != reductor_finder.end()){
			selected_reductor.push_back(reductor_finder[m]);
			list_monos.add_poly(reductor_finder[m]->poly);
			continue;
		}

		find_divisor_time-=get_ms_time();
		int r = gb.find_divisor( m );
		find_divisor_time+=get_ms_time();
		if( 0 > r ) {
			continue;
		}

		mono_t e = m / gb[r].head_term();
		if( e.is_zero() ) {
			//I choose not to add a copy into reductor storage.
			//The reason is I'd like to update some gb in reduction, since it doesn't change ht, nothing should go wrong.
			//I hope so.
			selected_reductor.push_back( &gb[r] );
			list_monos.add_poly(  gb[r].poly );
		} else {
			poly_t p = gb[r];
			p *= e;
			reductor_storage.push_back( p );
			reductor_finder[m] = &reductor_storage.back();
			selected_reductor.push_back( &reductor_storage.back() );
			list_monos.add_poly(  p.poly );
		}
	}
	return selected_reductor.size() > 0;

}

static
bool f4_reduce( mht_poly_ref & new_gb , const mono_t & min_mono, const f4_gb & gb )
{
	if( 0 == new_gb.size() ) return false;

	findreductor_time -= get_ms_time();

	mht_poly_ref::reverse_iterator nrit = new_gb.rbegin();

	//mono_t min_ht = gb.min_ht();
	DUMP(std::cout << "min_mono: " << min_mono << "\n");
	DUMP(std::cout << "composing matrix...\n" );

	map_mono_t list_monos;
	for(;nrit!=new_gb.rend();nrit++) { 
		list_monos.add_poly( (*nrit->second).poly); 
	}

	cout << "Number of pair polynomials: "<< new_gb.size()<<", at "<< list_monos.size()<<" column(s)"<< endl;

	static poly_sto reductor_storage;
	static unordered_map<mono_t, poly_t*, MonomialHash, MonomialEqual> reductor_finder;
	poly_ref selected_reductor;

	prepare_reductors( list_monos , selected_reductor, reductor_finder , reductor_storage , min_mono , gb );

	cout << "All monos including reductor: "<< list_monos.size()<< endl;

	findreductor_time += get_ms_time();

	int total_reductee_size = 0;
	int total_reductor_size = 0;

	for( int i=0; i< selected_reductor.size(); i++) {
		total_reductor_size += selected_reductor[i]->n_terms();
	}
	for(mht_poly_ref::reverse_iterator it=new_gb.rbegin();it!=new_gb.rend();it++) {
		total_reductee_size += it->second->n_terms();
	}

	cout << "Average length for reductees: "<< (double)total_reductee_size/ (double)new_gb.size();
	cout << "["<< new_gb.size() <<"], reductors: "<< (double)total_reductor_size/ (double)selected_reductor.size();
	cout << "["<< selected_reductor.size()<< "]"<<endl;

	total_reductee += new_gb.size();
	total_reductor += selected_reductor.size();



conversion_time -= get_ms_time();
	unordered_map<mono_t, int, MonomialHash, MonomialEqual> mono_to_int;
	unordered_map<int, mono_t> int_to_mono;
	int count = 0;
	int nreduc = selected_reductor.size();
	int nnew = new_gb.size();
	
	for(int i=0; i< selected_reductor.size() ; i++){
		int_to_mono[count] = selected_reductor[i]->head_term();
		mono_to_int[selected_reductor[i]->head_term()] = count;
		count++;
	}


	map_mono_t::iterator lit = list_monos.begin();
	while( lit != list_monos.end()){
		if(mono_to_int.find(*lit) != mono_to_int.end()){
			lit++;
			continue;
		}
		mono_to_int[*lit] = count;
		int_to_mono[count] = *lit;
		count++;
		lit++;
	}
	cout << "total #mono : " << mono_to_int.size() <<endl;

	if(mono_to_int.size() == nreduc){
		for(mht_poly_ref::reverse_iterator it=new_gb.rbegin(); it != new_gb.rend(); it++)
			it->second->set_zero();
		conversion_time += get_ms_time();
		return true;
	}	

	mzd_t* B = mzd_init(nreduc, mono_to_int.size() - nreduc);
	mzd_t* D = mzd_init(nnew , mono_to_int.size() - nreduc);
	vector<vector<int> > A_sparse;
	vector<vector<int> > C_sparse;
	for(int i=0; i< nreduc; i++){
		A_sparse.push_back(vector<int>());
		for(int j=0; j< selected_reductor[i]->n_terms(); j++){
			int ncolumn = mono_to_int[selected_reductor[i]->poly[j]];
			if(ncolumn < nreduc)
				A_sparse[i].push_back(ncolumn);
			else
				mzd_write_bit(B, i, ncolumn-nreduc, 1);
		}
	}

	count = 0;
	for(mht_poly_ref::iterator it=new_gb.begin();it!=new_gb.end();it++){
		C_sparse.push_back(vector<int>());
		for(int j=0; j< it->second->n_terms(); j++){
			int ncolumn = mono_to_int[it->second->poly[j]];
			if(ncolumn < nreduc)
				C_sparse[count].push_back(ncolumn);
			else
				mzd_write_bit(D, count, ncolumn-nreduc, 1);
		}
		count++;
	}



conversion_time += get_ms_time();
cout<< "conversion 1 done" <<endl<<flush;
computation_time -= get_ms_time();
	if(nreduc != 0){
//this is not a parallelizable code, so I don't like this, maybe m4ri is a better idea
	for(int i=nreduc-1; i>=0; i--){
		for(int j=1; j< A_sparse[i].size(); j++){
			mzd_combine_even_in_place(B, i, 0, B, A_sparse[i][j], 0);
		}
	}

//		mzd_top_echelonize_m4ri(A, 0);
//		mzd_submatrix(B, A, 0, nreduc, nreduc, mono_to_int.size());
//sparse mul, used only once so not fuctionized
		for(int i=0; i< nnew; i++){
			for(int j=0; j< C_sparse[i].size(); j++){
				mzd_combine_even_in_place(D, i, 0, B, C_sparse[i][j], 0);
			}
		}
	}
	mzd_echelonize_m4ri(D, 1, 0);// parameter to be determined
computation_time += get_ms_time();
cout<< "computation done" <<endl<<flush;
conversion_time -= get_ms_time();

	count = 0;
	for(mht_poly_ref::iterator it=new_gb.begin(); it != new_gb.end(); it++){
		it->second->set_zero();
		for(int j=0; j< (mono_to_int.size() - nreduc +63)/64; j++){
			uint64_t s = D->rows[count][j];
			if(s == 0) continue;
			for(int k = 0; k < 64; k++){
				if(1&(s>>k)){
					it->second->poly.set_coef(int_to_mono[j*64+k+nreduc]);
				}
			}
		}
		count++;
	}


	for(int i = 0; i< selected_reductor.size(); i++){
		mono_t head = selected_reductor[i] -> head_term();
		selected_reductor[i]->poly.vv.clear();//alright, brute force here
		selected_reductor[i]->poly.vv.push_back(head);
		for(int j=0; j< (mono_to_int.size()- nreduc + 63)/64; j++){
			uint64_t s = B->rows[i][j];
			if(s == 0) continue;
			for(int k = 0; k < 64; k++){
				if(1&(s>>k)){
					selected_reductor[i]->poly.vv.push_back(int_to_mono[j*64+k+nreduc]);
				}
			}
		}
	}




conversion_time += get_ms_time();
cout<< "conversion 2 done" <<endl<<flush;


	return true;
}




///////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////



static
void rem_red( poly_ref & g )
{
	poly_ref ng;
	for(unsigned i=0;i<g.size();i++){
		bool add = true;
		for(unsigned j=0;j<g.size();j++){
			if(i==j) continue;
			//if( g[j]->head_term().degree() < 1 ) continue; ////////////////////////////////////
			if(g[i]->head_term().is_divisible_by(g[j]->head_term()) ){
				add = false;
				break;
			}
		}
		if(add) ng.push_back(g[i]);
	}
	g = ng;
}




template <typename poly_t>
void initial_mix( std::deque<poly_t> & out , std::deque<poly_t> inp )
{
	for(unsigned i=0;i<inp.size();i++){
		if( inp[i].is_zero() ) continue;
		for(unsigned j=0;j<inp.size();j++) {
			if(inp[j].is_zero()) continue;
			if(i==j) continue;
			if( ! inp[j].get_coef(inp[i].head_term()).is_zero() ) inp[j] -= inp[i];
		}
	}

	for( unsigned i=0;i<inp.size();i++) {
		if( inp[i].is_zero() ) continue;
		out.push_back(inp[i]);
	}
	for( unsigned i=0;i<out.size();i++) {
		for(unsigned j=i+1;j<out.size();j++) {
			if( out[i].head_term() > out[j].head_term() ) std::swap( out[i] , out[j] );
		}
	}
}






void F4_algo( std::deque<base_poly_t> & out , std::deque<base_poly_t> & inp )
{

	if(inp.empty()) return;
	double total_time = get_ms_time();

	std::deque<base_poly_t> ii;
	//ii = inp;
	initial_mix<base_poly_t>( ii , inp );

	poly_sto all_polys;
	f4_gb gb;
	for(unsigned i=0;i<ii.size();i++) {
		poly_t pp;
		pp.poly = ii[i];
#ifdef __HAS_SIG__
		pp.s.idx = i;
#endif
		all_polys.push_back( pp );
		gb.add_poly( &all_polys.back() );
	}
	DUMP( std::cout << "[G] " << gb.size() << " polys.\n" );

	worklist_t wl;
	relation_t rels;
	gen_pair( wl , rels , gb , 0 , gb.size() );
	unsigned checked_polys = gb.size();

	DUMP( std::cout << "[GenPair] " << wl.size() << " pairs.\n\n" );

	poly_sto new_gb_sto;
	mht_poly_ref new_gb;

	int step_count = 0;
	while( !wl.empty() ) {
		iteration_count++;
		new_gb.clear();

		cout << "//////////////////////////// reduction step : "<< iteration_count << " //////////////////////////// "<< endl;
		cout << "Basis length: "<<gb .size() <<", queue length: "<< wl.size() <<endl;
		unsigned step_deg = wl.begin()->first;
		cout << "Step degree: " << step_deg << "\n";

		temp_counter = 0;
		genspoly_time -= get_ms_time();
		new_gb_sto.clear();
		gen_spoly( new_gb , new_gb_sto, gb , wl, rels);
		genspoly_time += get_ms_time();
		//DUMP(std::cout << "[SPOLY]: " << polys_in_g.size() << " org polys.\n" );
		//dump_g(polys_in_g);
		DUMP(std::cout << "[SPOLY]: " << new_gb.size() << " new polys, \n" );
		DUMP(std::cout << "[SPOLY]: " << temp_counter << " removed by syz check.\n");
//	dump_g(new_gb);
	
		mono_t min_mono = gb.min_ht();
		f4_reduce( new_gb , min_mono, gb);

		addpoly_time -= get_ms_time();
		ht_poly_ref sorter;//ht should not duplicate except 0
		for(mht_poly_ref::iterator it=new_gb.begin();it!=new_gb.end();it++){
			if(it->second->is_zero() ) { reduced_to_zero++; continue; }
			sorter[it->second->head_term()] = it->second;
		}

//	dump_g(sorter);

		step_count++;
		//DUMP( std::cout << "[[GB]]: " << gb.size() << "\n");
		//dump_g( gb.ht_polys );

		/// add new_gb. In magna, it is added in ht order, from small to large.
		/// So I added another container to simulate it

		unsigned zero = 0;
		unsigned add_deg = sorter.begin()->first.degree();
		for(ht_poly_ref::iterator it=sorter.begin(); it!=sorter.end();it++){
			if(it->second->is_zero() ) { zero++; reduced_to_zero++; continue; }
			if( add_deg != it->first.degree() ) break;
			poly_t pp = *it->second;
			all_polys.push_back(pp);
			bool r = gb.add_poly( &all_polys.back() );
			if( 1 >= pp.degree() ) {
				DUMP( std::cout << pp << "\n" );
			}
			//bool r = gb.add_poly(it->second);
			if( !r ) {
				//ht_poly_ref::iterator ii = polys_in_g.find( it->second->head_term() );
				//DUMP( std::cout << "in polys_in_g: " << !(polys_in_g.end()==ii) << "\n" );
				
			}
		}
		bool has_mutant = (add_deg < step_deg) && (0 != add_deg);
		DUMP( std::cout << "has mutnat ?" << has_mutant << "\n" );
		if( !has_mutant ) {
			worklist_t::iterator it;
			for(it=wl.begin();it!=wl.end();it++)
				if(it->first > step_deg) break;
			wl.erase(wl.begin(), it);
		}
		addpoly_time += get_ms_time();

		DUMP( std::cout << "[REDUCE] " << new_gb.size()-zero << "new polys, " << zero << " zero.\n" );
		DUMP( std::cout << "[[G]] " << gb.size() << " polys.\n" );

		gen_pair( wl , rels , gb , checked_polys , gb.size() );
		checked_polys = gb.size();


		DUMP( std::cout << "[GenPair] " << wl.size() << " pairs.\n\n" );
	}


	rem_red(gb.polys);

	dump_g1( gb.polys );

	guass_elim(gb.polys);

	ht_poly_ref sorter;//ht should not duplicate except 0
	for(int i = 0; i < gb.polys.size(); i++){
		sorter[gb.polys[i]->head_term()] = gb.polys[i];
	}

	for(ht_poly_ref::reverse_iterator it=sorter.rbegin(); it!=sorter.rend();it++){
		out.push_back(it->second->poly);
	}
	
	total_time = get_ms_time() - total_time;

	cout << "removed by syz check: " << removed_by_syz_check << "\n";
	cout << "total pair: " << total_pair << "\n";
	cout << "reduced to zero: " << reduced_to_zero << "\n";
	cout << "total pair in elim " << total_elim_size << "\n";
	cout << "total pair elimed " << total_pair - total_reductee << "\n";
	cout << "total reductee " << total_reductee << "\n";
	cout << "total reductor  " << total_reductor << "\n";
	cout << "max elim size: " << max_elim_size << "\n";
	cout << "iteration count: " << iteration_count << "\n";
	cout << "average elim mat size: "<< total_elim_size / iteration_count << "\n";
	cout << "genpair time: "<< genpair_time <<"ms\n";
	cout << "genspoly time: "<< genspoly_time <<"ms\n";
	cout << "find reductor time: "<< findreductor_time <<"ms\n";
	cout << "find divisor time: "<< find_divisor_time <<"ms\n";
	cout << "total convertion time: "<< conversion_time <<"ms\n";
	cout << "total computation time: "<< computation_time <<"ms\n";
	cout << "addpoly time: "<< addpoly_time <<"ms\n";
	cout << "total time: "<< total_time <<"ms\n";
}






