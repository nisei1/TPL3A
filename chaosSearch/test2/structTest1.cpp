#include <iostream>
#include <vector>

using namespace std;

#define VECTOR_INDEX 10

struct Vec
{
    Vec() : vec1(VECTOR_INDEX, 0), vec2(VECTOR_INDEX, 0) {}
    vector<int> vec1;
    vector<int> vec2;
};

int main(int argc, char const *argv[])
{
    Vec test1;

    for (int i = 0; i < test1.vec1.size(); i++)
    {
        test1.vec1[i] = i;
        test1.vec2[i] = i + 10;
    }

    for (int i = 0; i < test1.vec1.size(); i++)
    {
        cout << test1.vec1[i] << endl;
    }
    for (int i = 0; i < test1.vec2.size(); i++)
    {
        cout << test1.vec2[i] << endl;
    }

    return 0;
}
