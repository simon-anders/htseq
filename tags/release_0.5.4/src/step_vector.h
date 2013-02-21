#ifndef _STEP_VECTOR_H_
#define _STEP_VECTOR_H_

#include <map>
#include <stdexcept>
#include <climits>
#include <iostream>  //for now only

template< class T >
class step_vector {
  protected:
   std::map< long int, T > m;
  public: 
   static const long int min_index;
   static const long int max_index;
   typedef typename std::map< long int, T >::const_iterator const_iterator;
   step_vector( );
   const T operator[]( long int i ) const;
   void set_value( long int from, long int to, T value );
   void add_value( long int from, long int to, T value );
   void apply_to_values( long int from, long int to, void (*func)( T & val ) );
   const_iterator get_values( long int from ) const;
   const_iterator begin( ) const;
   const_iterator end( ) const;
};

template< class T >
const long int step_vector<T>::min_index = LONG_MIN;

template< class T >
const long int step_vector<T>::max_index = LONG_MAX;

template< class T >
step_vector<T>::step_vector( ) 
{
   m[ min_index ] =  T();
}

template< class T >
const T step_vector<T>::operator[]( long int i ) const
{
   const_iterator it = m.upper_bound( i );
   it--;
   return it->second;
}

template< class T >
void step_vector<T>::set_value( long int from, long int to, T value )
{
   if( from > to )
      throw std::out_of_range( "Indices reversed in step_vector." );

   // Unless the new step extends to the end, we need to insert a new
   // value afterwards unless the step to the right has the same value
   if( to < max_index ) {
      T next_value = (*this)[to+1];
      if( !( next_value == value ) )
         m[ to + 1 ] = next_value;
   }

   // Find the left step, i.e., the step whose start is smaller or equal
   // to 'from':
   typename std::map< long int, T>::iterator left = m.upper_bound( from );
   left--;
   assert( left->first <= from );
                     
   // Get rid of the steps present between from and to
   typename std::map< long int, T>::iterator it = m.lower_bound( from );
   if( it->first == from )
      it++;
   assert( it->first > from );
   if( it->first <= to ) {
      m.erase( it, m.upper_bound( to ) );
   }

   if( ! (left->second == value ) ) {
      if( left->first != from )
         // Insert a new step
         m[ from ] = value;
      else {
         // We have from == left->first, so the step is already present.
         // Would changing m[from] to value make it equal to its left
         // neighbor?
         if( left == m.begin() )      
            // no, there is no left neighbor
            m[ from ] = value;
         else {
            typename std::map< long int, T>::iterator leftleft = left;
            leftleft--;
            if( ! ( leftleft->second == value ) )
               // ok, change the value
               m[ from ] = value;
            else
               // no, rather delete the step
               m.erase( left );
         }
      }
   }
}

template< class T >
void step_vector<T>::add_value( long int from, long int to, T value )
{
   if( from > to )
      throw std::out_of_range( "Indices reversed in step_vector." );
   
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


#endif //_STEP_VECTOR_H_
