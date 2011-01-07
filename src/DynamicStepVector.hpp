#include <map>
#include <utility>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <utility>
#include <limits>

#define NDEBUG

template< typename TValue >
class Value{
public:
    virtual ~Value() {
//#ifdef NDEBUG
//        std::cout << "Value Base Destructor called!" << std::endl;
//#endif
    };
    virtual std::string info() = 0;
    virtual bool multiple(){
        return false;
    }
    virtual TValue & operator[] ( size_t ) = 0;
    virtual void set( TValue const &  ) = 0;
    virtual size_t size() = 0;

};

template< typename TValue >
std::string display_vector( std::vector< TValue > const & vector ){
    std::ostringstream os;
    os << "Vector < ";
    typename std::vector< TValue >::const_iterator it = vector.begin();
    while( it != vector.end() ){
        os << *it << " ";
        ++it;
    }
    os << " >";
    return os.str();
}

template< typename TValue >
class SingleValue : public Value< TValue > {
public:
    SingleValue() : value() { };
    SingleValue( SingleValue< TValue > const & other ){
        value = other.get();
    };
    SingleValue( TValue const & val ) : value( val ) {};
    
    virtual ~SingleValue() {
//    #ifdef NDEBUG
//        std::cout << this->info() << " Destructor called!" << std::endl;
//    #endif
    };
    
    operator TValue() const {
        return value;
    }
    
    template< typename TFkt >
    void apply( TFkt & fkt ){
        value = fkt( value );
    }
    
    virtual TValue & operator[] ( size_t pos ){
        return value;
    }
    
    virtual void set( TValue const & val ){
        value = val;
    }

    virtual size_t size(){
        return 1;
    }
    
    virtual std::string info(){
        std::ostringstream os;
        os << "SingleValue < " << value << " >";
        return os.str();
    }
    
    virtual bool multiple(){
        return false;
    }
private:
    TValue value;
};

template< typename TValue >
class MultipleValue : public Value< TValue > {
public:

    virtual ~MultipleValue() {
    //#ifdef NDEBUG
    //    std::cout << this->info() << " Destructor called!" << std::endl;
    //#endif
    };

    MultipleValue( TValue const & val, size_t count ) : values( count, val ) { 
        //assert( count > 0 && " attempted to create MV with length <= 0!" );
    };
    
    MultipleValue( std::vector< TValue > const & vals ) : values( vals ) {
        //assert( vals.size() > 0 && " attempted to create MV with length <= 0!" );
    };
    
    MultipleValue( Value< TValue > const & val, size_t count ) : values( count, val.get() ) {
        //assert( count > 0 && " attempted to create MV with length <= 0!" );
    };
    
    MultipleValue( SingleValue< TValue > const & val, size_t count ) : values( count, val.get() ) {
        //assert( count > 0 && " attempted to create MV with length <= 0!" );
    };

    void append( std::vector< TValue > const & other ){
        values.insert( values.end(), other.begin(), other.end() );
    }
    
    template< typename TFkt >
    void apply( TFkt & fkt ){
        transform( values.begin(), values.end(), values.begin(), fkt );
    }
    
    template< typename TFkt >
    void apply( size_t from, size_t to, TFkt & fkt ){
        std::transform( values.begin() + from, values.begin() + to, values.begin() + from, fkt );
    }
    
    template< typename TFkt >
    void apply_from( size_t from, TFkt & fkt ){
        apply( from, values.size(), fkt );
    }
    
    template< typename TFkt >
    void apply_to( size_t to, TFkt & fkt ){
        apply( 0, to, fkt );
    }
    
    std::vector< TValue > const & get_all(){
        return values;
    }
    
    std::vector< TValue > get( size_t from, size_t to ){
        return std::vector< TValue >( values.begin() + from, values.begin() + to );
    }

    void swap( std::vector< TValue > const & other ) {
        values.swap( other );
    }
    
    virtual TValue & operator[] ( size_t pos ){
        return values[ pos ];
    }
    
    virtual void set( TValue const & val ){
        values = std::vector< TValue >( values.size(), val );
    }
    
    void set_from_to( TValue const & val, size_t pos_from, size_t pos_to ){
        typename std::vector< TValue >::iterator it = values.begin() + pos_from;
        while( it != values.begin() + pos_to ){
            *it = val;
            ++it;
        }
    }
    
    void del_from( size_t pos ){
        //assert( pos > 0 && " attempted to delete whole values-vector!" );
        values.erase( values.begin() + pos, values.end() );
    }
    
    void push_back( TValue const & val, size_t count ) {
        for( unsigned int i = 0; i < count; ++i ){
            values.push_back( val );
        }
    }
    
    virtual size_t size(){
        return values.size();
    }
    
    virtual std::string info(){
        std::ostringstream os;
        os << "MultipleValue < ";
        typename std::vector< TValue >::iterator it = values.begin();
        while( it != values.end() ){
            os << *it << " ";
            ++it;
        }
        os << ">";
        return os.str();
    }
    
    virtual bool multiple(){
        return true;
    }

private:
    std::vector< TValue > values;
};

template< typename TValue >
struct Add{
    Add( TValue const & v ) : offset( v ) { };
    TValue operator() ( TValue const & val ){
        return val + offset;
    }
    
    TValue offset;
};

template< typename TValue >
struct Invert{        
    TValue operator() ( TValue const & val ){
        return -val;
    }
};

template<>
struct Invert< std::string >{
    std::string operator() ( std::string const & val ){
        return val;
    }
};

template< typename TValue >
struct DSVIter{

    DSVIter(void) : start(0), stop(0), pos(0), it() { }; //default-ctor
    ~DSVIter(void) { }; //default-dtor
    
    DSVIter( DSVIter< TValue > const & other )
     :  start( other.start ),
        stop( other.stop ),
        pos( other.pos ),
        it( other.it ),
        map( other.map ) { }; //copy-ctor
        
    DSVIter( std::map< long int, Value< TValue >* > * m, long int from, long int to, bool reverse = false )
     :  start( from ),
        stop( to ),
        pos( reverse ? to - 1 : from ),
        it( --( m->upper_bound( pos ) ) ),
        map( m ) {
            if( reverse and !( it->second->multiple() ) ){
                pos = it->first;
            }
        };
    
    std::pair< long int, TValue > next(){
        std::pair< long int, TValue > p = std::make_pair( pos, ( *( it->second ) )[ pos - it->first ] );
        ++pos;
        
        if( it != map->end() && pos - it->first >= it->second->size() ){ //end of step
            ++it;
            pos = it->first;
        }
        return p;
    };
    
    std::pair< long int, TValue > prev(){
        std::pair< long int, TValue > p = std::make_pair( pos, ( *( it->second ) )[ pos - it->first ] );
        --pos;
        
        if( it != map->begin() && pos < it->first ){ //begin of step
            --it;
            if( !it->second->multiple() ){
                pos = it->first;
            }
        }
        return p;
    };
    
    bool valid(){
        return it != map->end() && pos < stop && pos >= start;
    }
    
    std::string info() {
        std::ostringstream os;
        os << pos << ":\t" << it->second->info() << " | " << it->second->size();
        return os.str();
    }
    
    long int start, stop, pos, step_size;
    typename std::map< long int, Value< TValue >* >::iterator it;
    std::map< long int, Value< TValue >* > * map;
};

template< typename TKey, typename TValue >
class DSV{
public:
    typedef std::map< TKey, Value< TValue >* > Map;
    typedef SingleValue< TValue > TSV;
    typedef MultipleValue< TValue > TMV;

    DSV() : threshold( 1 ) {
        steps[ static_cast<TKey>(0) ] = new TSV();
    };

    DSV( size_t t ) : threshold( t ) {
        steps[ static_cast<TKey>(0) ] = new TSV();
    };
    
    DSV( DSV< TKey, TValue > const & other ) : threshold( other.get_threshold() ) , steps( other.get_steps() ) { };
    
    ~DSV(){
        typename Map::iterator it = steps.begin();
        while( it != steps.end() ){
            delete it->second;
            ++it;
        };
    };
    
    void add( TKey const & from, TKey const & to, TValue const & offset ){
        Add< TValue > adder( offset ); // The black adder :D
        apply( from, to, adder );
    }

    void clear(){
        typename Map::iterator it = steps.begin();
        while( it != steps.end() ){
            delete it->second;
            ++it;
        };
        steps.clear();
        steps[ static_cast<TKey>(0) ] = new TSV();
    }

    void invert( TKey const & from, TKey const & to, TValue const & offset ){
        Invert< TValue > inv; // Teh 3v0l Invert0r
        apply( from, to, inv );
    }
    
    TValue & get( TKey const & key ){
        typename Map::iterator it = get_iter( key );
        if( it->second->multiple() ){
            return ( *( static_cast< TMV* >( it->second ) ) )[ key - it->first ];
        }else{
            return ( *( static_cast< TSV* >( it->second ) ) )[ 0 ];
        }
    }
    
    std::vector< TValue > get( TKey const & from, TKey const & to ){
        //assert( from < to && " Right-Hand-Interval-Border underrun!" );
        std::vector< TValue > values;
        typename Map::iterator it_from = get_iter( from );
        typename Map::iterator it_to = get_iter( to );
        
        if( it_from == it_to ){
            if( it_from->second->multiple() ){
                return (static_cast<TMV*>(it_from->second))->get( from - it_from->first, to - it_from->first );
            }else{
                return std::vector< TValue >( to - from, (*(static_cast<TSV*>(it_from->second)))[0] );
            }
        }

        typename Map::iterator it_step = it_from;
        ++it_step;

        if( it_from->second->multiple() ){
            values = (static_cast<TMV*>(it_from->second))->get( from - it_from->first, it_step->first - it_from->first);
        }else{
            values = std::vector< TValue >( it_step->first - from, (*(static_cast<TSV*>(it_from->second)))[0] );
        }
        
        ++it_from;
        ++it_step;
        std::vector< TValue > tmp;
        while( it_from != it_to ){
            if( it_from->second->multiple() ){
                tmp = ( (static_cast<TMV*>(it_from->second))->get_all() );
            }else{
                tmp = std::vector< TValue >( it_step->first - it_from->first, (*(static_cast<TSV*>(it_from->second)))[0] );
            }
            values.insert( values.end(), tmp.begin(), tmp.end() );
            ++it_from;
            ++it_step;
        }
        
        if( it_to->second->multiple() ){
            tmp = (static_cast<TMV*>(it_to->second))->get( 0, to - it_to->first );
        }else{
            tmp = std::vector< TValue >( to - it_to->first, (*(static_cast<TSV*>(it_to->second)))[0] );
        }
        values.insert( values.end(), tmp.begin(), tmp.end() );
        
        return values;
    }

    template< typename TFkt >
    void apply( TKey const & from, TKey const & to, TFkt & fkt ){
        typename Map::iterator it_from = get_iter( from );
        typename Map::iterator it_to = get_iter( to );
        
        if( it_from == it_to ){
            if( it_from->second->multiple() ){
                (static_cast<TMV*>(it_from->second))->apply( from - it_from->first, to - it_from->first, fkt );
            }else{
                TValue tmp = (*(static_cast<TSV*>(it_from->second)))[0];
                set( from, fkt( tmp ) );
                set( to, tmp );
            }
            return;
        }
            
        if( it_from->second->multiple() ){
            (static_cast<TMV*>(it_from->second))->apply_from( from - it_from->first, fkt );
            ++it_from;
            apply( it_from->first, to, fkt );
        }else{
            TValue tmp = (*(static_cast<TSV*>(it_from->second)))[0];
            TKey from_new = (++it_from)->first;
            set( from, fkt( tmp ) );
            apply( from_new, to, fkt );
        }

    }
    
    void set( TKey const & key, TValue const & val ) {
        std::cout << "set called on " << this << " with key:val of " << key << ":" << val << std::endl;
        typename Map::iterator it = this->get_iter( key );
        
        if( it->second->multiple() ){
            if( it->first != key ){
                ( *( static_cast<TMV*>( it->second ) ) ).del_from( key - it->first );
                delete steps[ key ];
                steps[ key ] = new TSV( val );
                refurbish( key );
            }else{
                delete steps[key];
                steps.erase( it );
                steps[ key ] = new TSV( val );
            }
//            set( key, val );
        }else{
            if( it->first != key ){
                delete steps[ key ];
                steps[ key ] = new TSV( val );
                refurbish( key );
            }else{
                it->second->set( val );
            }
        }
    }
    
    void set( TKey const & from, TKey const & to, TValue const & val ) {

        std::cout << "set called on " << this << " with <from-to>:val of <" << from << "-" << to << ">:" << val << std::endl;

        typename Map::iterator it_from = this->get_iter( from );
        typename Map::iterator it_to = this->get_iter( to );
        
        if( it_to != it_from ){
        
            if( it_from->second->multiple() ){
                if( from != it_from->first ){
                    static_cast<TMV*>( it_from->second )->del_from( from - it_from->first );
                    ++it_from;
                }
            }else{
                ++it_from;
            }
            
            typename Map::iterator it_tmp = it_from;
            while( it_tmp != it_to ){
                delete it_tmp->second;
                ++it_tmp;
            }
            steps.erase( it_from, it_to ); //erase steps in between
        
            if( it_to->second->multiple() ){

                std::vector< TValue > tmp = static_cast<TMV*>( it_to->second )->get( to - it_to->first, it_to->second->size() );
                delete steps[to];
                steps.erase( it_to );
                steps[ to ] = new TMV( tmp );

            }else{

                TValue to_val = ( *( it_to->second ) )[0];
                delete it_to->second;
                steps.erase(it_to);
                set( to, to_val );
                
            }
            
            set( from, val );
            
        }else{
        
            if( it_from->second->multiple() ){
                static_cast<TMV*>( it_from->second )->set_from_to( val, from - it_from->first, to - it_to->first );
            }else{
                TValue to_val = ( *( it_to->second ) )[0];
                set( from, val );
                set( to, to_val );
            }
            
        }
    }
    
    void refurbish( TKey const & key ){
        typename Map::iterator it = steps.find( key ); // get iterator to key's steps
        //assert( it != steps.end() && " attempted to refurbish key that doesn't exist!" );
        
        TKey new_key = (++it)->first;
        --it;
        
        if( it != steps.begin() ){
            
            typename Map::iterator it_low = --( steps.find( key ) ); // get iterator to key's predecessor
            
            if( it->first - it_low->first <= static_cast<TKey>(threshold) ){ // it_low must be made MultipleValued and may be merged with it
                if( !it_low->second->multiple() ){
                    std::vector< TValue > values( it->first - it_low->first, (*(it_low->second))[0] );    
                    delete it_low->second;
                    it_low->second = new TMV( values );
                }
                if( it_low != steps.end() ){
                    if( it->second->multiple() && it != steps.end() ){
                        (static_cast<TMV*>(it_low->second))->append( (static_cast<TMV*>(it->second))->get_all() );
                        delete it->second;
                        steps.erase( it );
                        it = steps.find( new_key );
                    }
                }
                if( it_low != steps.begin() ){
                    typename Map::iterator it_lower = it_low;
                    --it_lower;
                    if( it_lower->second->multiple() ){
                        (static_cast<TMV*>(it_lower->second))->append( (static_cast<TMV*>(it_low->second))->get_all() );
                        delete it_low->second;
                        steps.erase( it_low );
                    }
                }
            }else{
                return;
            }
            
        }
        ++it;
        if( it != steps.end() ){
            refurbish( it->first ); // check to see if we have to merge the next step also
        }
    }
    
    typename Map::iterator get_iter( TKey const & key ){
        typename Map::iterator it = steps.upper_bound( key );
        //assert( it != steps.begin() && " key with upper_bound == steps.begin() requested, this should not have happened!" );
        //assert( it->first != key && "Upper bound Misbehaves!" );
        return --( it );
    }
            
    std::string info(){
        std::ostringstream os;
        typename Map::iterator it = steps.begin();
        while( it != steps.end() ){
            os << it->first << ": " << it->second->info() << "\n";
            ++it;
        }
        return os.str();
    }
    
    std::vector< std::string > rinfo(){
        std::vector< std::string > ret;
        std::ostringstream os;
        
        typename Map::iterator it = steps.begin();

        while( it != steps.end() ){
            
            os << it->first << ": " << it->second->info();

            ret.push_back( os.str() );
            
            os.str("");
            
            ++it;
            
        }

        return ret;
    }
    
    size_t size(){
        return steps.size();
    }
    
    Map get_map(){
        return steps;
    }
    
    Map const & get_steps() const{
        return steps;
    }
    
    size_t get_threshold() const{
        return threshold;
    }
    
    DSVIter< TValue > get_step_iter( TKey from = static_cast< TKey >( 0 ), TKey to = std::numeric_limits< TKey >::max(), bool reverse = false ){
        return DSVIter< TValue >( &steps, from, to, reverse );
    }
    
private:
    size_t threshold;
    std::map< TKey, Value< TValue >* > steps;
};

