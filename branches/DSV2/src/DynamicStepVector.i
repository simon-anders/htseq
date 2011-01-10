/* File: DynamicStepVector.i */
%module DynamicStepVector
%include "std_string.i"
%include "AutoPyObjPtr.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_pair.i"
%{
#define SWIG_FILE_WITH_INIT
#include "DynamicStepVector.hpp"
#include <map>
#include <utility>
#include <limits>
//HTSeq::DSV< long int, int >;
//HTSeq::DSV< long int, bool >;
//HTSeq::DSV< long int, double >;

template<>
struct Add< AutoPyObjPtr >{
    Add( AutoPyObjPtr const & v ) : offset( v ) { };
    AutoPyObjPtr operator() ( AutoPyObjPtr const & val ){
        return val;
    }
    
    AutoPyObjPtr offset;
};

template<>
struct Invert< AutoPyObjPtr >{        
    AutoPyObjPtr operator() ( AutoPyObjPtr const & val ){
        return val;
    }
};

std::ostream& operator<< (std::ostream& stream, const AutoPyObjPtr& ) {
    stream << "< PythonAutoPtrObject >";
    return stream;
}

%}

namespace std {
   %template( strvector ) vector< string >;
   %template( intvector ) vector< int >;
   %template( dblvector ) vector< double >;
   %template( pyvector ) vector< AutoPyObjPtr >;
   %template( intintmap ) map< long int, int >;
   %template( intsteppair ) pair< long int, int >;
   %template( dblsteppair ) pair< long int, double >;
   %template( strsteppair ) pair< long int, string >;
   %template( pysteppair ) pair< long int, AutoPyObjPtr >;
};

template< typename TValue >
std::string display_vector( std::vector< TValue > const & );
   
template< typename TKey, typename TValue >
class DSV{
public:
    typedef std::map< TKey, Value< TValue >* > Map;
    typedef SingleValue< TValue > TSV;
    typedef MultipleValue< TValue > TMV;

    DSV( void );
    ~DSV( void );

    DSV( size_t t );

    DSV( DSV< TKey, TValue > const & other );
    
    void add( TKey const & from, TKey const & to, TValue offset );

    void clear();

    void invert( TKey const & from, TKey const & to, TValue offset );
    
    TValue get( TKey const & key );
    
    std::vector< TValue > get( TKey const & from, TKey const & to );

    template< typename TFkt >
    void apply( TKey const & from, TKey const & to, TFkt & fkt );
    
    void set( TKey const & key, TValue val );
    
    void set( TKey const & from, TKey const & to, TValue val );
    
    void set( DSVIter< TValue > other_it );
    
    void refurbish( TKey const & key );
    
    typename Map::iterator get_iter( TKey const & key );
            
    std::string info();
    
    std::vector< std::string > rinfo();
    
    size_t size();
    
    Map get_map();
    
    Map const & get_steps();
    
    DSVIter< TValue > get_step_iter( TKey from = static_cast< TKey >( 0 ), TKey to = std::numeric_limits< TKey >::max(), bool reverse = false );
    
private:
    size_t threshold;
    std::map< TKey, Value< TValue >* > steps;
};

%template(intDSV) DSV< long int, int >;
//%template(boolDSV) DSV< long int, bool >;
%template(floatDSV) DSV< long int, double >;
%template(strDSV) DSV< long int, std::string >;
%template(pyDSV) DSV< long int, AutoPyObjPtr >;

template< typename TValue >
struct DSVIter{
    DSVIter(void); //default-ctor
    ~DSVIter(void); //default-dtor
    DSVIter( DSVIter< TValue > const & other ); //copy-ctor
    DSVIter( std::map< long int, Value< TValue >* > * m, long int from, long int to, bool rev );
    
    std::pair< long int, TValue > next();
    std::pair< long int, TValue > prev();
    
    bool valid();
    std::string info();
    long int start, stop, pos, step_size;
    std::map< long int, Value< TValue >* >::iterator it;
    std::map< long int, Value< TValue >* > * map;
    bool reverse;
};

%template(intDSVIter) DSVIter< int >;
%template(floatDSVIter) DSVIter< double >;
%template(strDSVIter) DSVIter< std::string >;
%template(pyDSVIter) DSVIter< AutoPyObjPtr >;

%{

%}


