#include <iostream>
#include <cstdio>

void Swap(const char * & left, const char * & right )
//void Swap(char* left, char* right )
{
  const char * const temp = left;
  left = right;
  right = temp;
}

void Show(const char*& arg)
{
  std::cout << arg << std::endl;
}

void ShowAll( const char**& args )
{
  std::cout << "OK" << std::endl;
  std::cout << args[0] << std::endl;
  std::cout << args[1] << std::endl;
}

int main(int argc, char* argv[])
{
  const char* x = "Hello"; // Declare and initialize a constant pointer variable `x' to the <char> type as "Hello"
  // This syntax makes the contents (Hello) to be a constant not the pointer so that you can change where the pointer pointing.
  // Therefore x = y sentence at line 33 is just OK.
  const char* y = "World"; // Declare and initialize a constant pointer variable `y' to the <char> type as "World"
  std::cout << "x = " << x << ", y = " << y << std::endl; // Print out to console "x = Hello, y = World"
  Swap(x, y); // Call the Swap() function
  std::cout << "x = " << x << ", y = " << y << std::endl; // You can check x and y are switched.
  x = y; // Now pointer x points where y pointing.
  std::cout << "x = " << x << ", y = " << y << std::endl;
  Show((const char*&)argv[0]);
  ShowAll((const char**&)argv);

  return 0;
}
