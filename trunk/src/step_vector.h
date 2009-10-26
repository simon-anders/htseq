#ifndef _STEP_VECTOR_H_
#define _STEP_VECTOR_H_

#include <map>
#include <stdexcept>
#include <limits>
#include <iostream>  //for now only

template< class T >
class step_vector {
  protected:
   std::map< long int, T > m;
  public: 
   long int min_index;
   long int max_index;
   typedef typename std::map< long int, T >::const_iterator const_iterator;
   step_vector( long int length, long int min_index_=0 );
   const T operator[]( long int i ) const;
   void set_value( long int from, long int to, T value );
   void add_value( long int from, long int to, T value );
   void apply_to_values( long int from, long int to, void (*func)( T & val ) );
   const_iterator get_values( long int from ) const;
   const_iterator begin( ) const;
   const_iterator end( ) const;
};

template< class T >
step_vector<T>::step_vector( long int length, long int min_index_ ) 
 : min_index( min_index_ ),
   max_index( min_index_ + length - 1 )
{
   m[ min_index ] =  T();
}

template< class T >
const T step_vector<T>::operator[]( long int i ) const
{
   if( i > max_index ) 
      throw std::out_of_range( "Index too large in step_vector." );
   if( i < min_index ) 
      throw std::out_of_range( "Index too small in step_vector." );
   const_iterator it = m.upper_bound( i );
   it--;
   return it->second;
}

template< class T >
void step_vector<T>::set_value( long int from, long int to, T value )
{
   if( from > to )
      throw std::out_of_range( "Indices reversed in step_vector." );
   if( to > max_index )
      throw std::out_of_range( "Index too large in step_vector." );
   if( from < min_index )
      throw std::out_of_range( "Index too small in step_vector." );
   if( to < max_index ) {
      T next_value = (*this)[to+1];
      m[ to + 1 ] = next_value;
   }
   typename std::map< long int, T>::iterator it = m.lower_bound( from );
   if( it->first <= to ) {
      m.erase( it, m.upper_bound( to ) );
   }
   m[ from ] = value;
}

template< class T >
void step_vector<T>::add_value( long int from, long int to, T value )
{
   if( from > to )
      throw std::out_of_range( "Indices reversed in step_vector." );
   if( to > max_index )
      throw std::out_of_range( "Index too large in step_vector." );
   if( from < min_index )
      throw std::out_of_range( "Index too small in step_vector." );
   
   if( to < max_index ) {
      T next_value = (*this)[to+1];
      m[ to + 1 ] = next_value;
   }

   typename std::map< long int, T>::iterator it = m.upper_bound( from );
   it--;
   bool need_to_insert_step_at_from = it->first < from;
   T old_val_at_from;
   if( need_to_insert_step_at_from ) {
      old_val_at_from = it->second;
      it++;
   }
   // Now, it points to the first element with it->first >= from
   
   for( ; it != m.end() && it->first <= to; it++ ) 
      it->second += value;
   
   if( need_to_insert_step_at_from )
      m[ from ] = old_val_at_from + value;
}

template< class T >
void step_vector<T>::apply_to_values( long int from, long int to, 
      void (*func)( T & val ) )
{
   if( from > to )
      throw std::out_of_range( "Indices reversed in step_vector." );
   if( to > max_index )
      throw std::out_of_range( "Index too large in step_vector." );
   if( from < min_index )
      throw std::out_of_range( "Index too small in step_vector." );
   
   if( to < max_index ) {
      T next_value = (*this)[to+1];
      m[ to + 1 ] = next_value;
   }

   typename std::map< long int, T>::iterator it = m.upper_bound( from );
   it--;
   bool need_to_insert_step_at_from = it->first < from;
   T old_val_at_from;
   if( need_to_insert_step_at_from ) {
      old_val_at_from = it->second;
      it++;
   }
   // Now, it points to the first element with it->first >= from
   
   for( ; it != m.end() && it->first <= to; it++ ) 
      func( it->second );
   
   if( need_to_insert_step_at_from ) {
      func( old_val_at_from );
      m[ from ] = old_val_at_from;
   }
}


template< class T >
typename step_vector<T>::const_iterator step_vector<T>::get_values( long int from ) const
{
   return --m.upper_bound( from );
}

template< class T >
typename step_vector<T>::const_iterator step_vector<T>::begin( ) const
{
   return m.begin();
}   

template< class T >
typename step_vector<T>::const_iterator step_vector<T>::end( ) const
{
   return m.end();
}   


//==================== extra stuff for Python =============



      
                  
/*
void test()
{
   step_vector_for_python<double> sv( 100 );
   sv.get_all_values_pystyle();
   step_vector_for_python<double>::const_iterator it1 = sv.begin();
   step_vector_for_python<double>::const_iterator it2 = sv.end();
   step_vector_pystyle_iterator<double>( it1, it2 );
}
*/

#endif //_STEP_VECTOR_H_
