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

void ShowAll( const char**& args)
{
  std::cout << "OK" << std::endl;
  std::cout << args[0] << std::endl;
  std::cout << args[1] << std::endl;
}

int main(int argc, char* argv[])
{
  const char* x = "Hello";
  const char* y = "World";
  std::cout << "x = " << x << ", y = " << y << std::endl;
  Swap(x, y);
  std::cout << "x = " << x << ", y = " << y << std::endl;
  x = y;
  std::cout << "x = " << x << ", y = " << y << std::endl;
  Show((const char*&)argv[0]);
  ShowAll((const char**&)argv);

  return 0;
}
