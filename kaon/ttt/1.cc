#include "myClass.h"

int j = MyClass::i;

class Foo {
public:
  Foo(int x) { std::cout <<"Foo: " << x << std::endl;}
};

Foo foo(MyClass::i);

int main() {

  std::cout << j << std::endl;
};


