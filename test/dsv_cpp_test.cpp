#include <string>
#include <iostream>
#include "DynamicStepVector.hpp"

int main()
{
    DSV< int, std::string > dsv;
    dsv.set( 3, 6, std::string( "Foo" ) );
    dsv.set( 10, 15, std::string( "Bar" ) );
    dsv.set( 12, 13, std::string( "Wakka" ) );
    dsv.set( 12, 13, std::string( "Quacko" ) );
    for( int i = 0; i < 20; i++ )
       std::cout << i << ": " << dsv.get(i) << std::endl;
    std::cout << "Overwriting with Zeroes!" << std::endl;
    dsv.set( 0, 100, std::string( "000" ) );
    for( int i = 0; i < 20; i++ )
       std::cout << i << ": " << dsv.get(i) << std::endl;
    dsv.clear();
}
