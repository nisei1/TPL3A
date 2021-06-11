#include <iostream>
#include <vector>
#include <random>
#include <ctime>

// using std::cin;
// using std::cout;
// using std::endl;
// using std::swap;
// using std::vector;

using namespace std;

int main(void)
{
	time_t t;
	int r1, r2;
	int n;

	cout << "Please input data ";
	cin >> n;

	vector<int> vec;

	for (int i = 0; i < n; i++)
	{
		vec.push_back(i);
	}
	//シャッフル
	// std::random_device seed_gen;
	// std::mt19937 engine(seed_gen());
	// std::shuffle(vec.begin(), vec.end(), engine);

	cout << "before" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << vec.at(i) << " ";
	}
	cout << endl;

	t = time(NULL);
	std::random_device rnd;
	std::mt19937 mt(rnd() * t);

	std::uniform_int_distribution<> dist1(0, n - 1);
	std::uniform_int_distribution<> dist2(0, n - 1);

	r1 = dist1(mt);
	r2 = dist2(mt);

	cout << "range" << r1 << " " << r2 << endl;

	while (r1 == r2 || (r1 == 0 && r2 == n - 1) || (r2 == 0 && r1 == n - 1) || abs(r1 - r2) > 1)
	{
		r1 = dist1(mt);
		r2 = dist2(mt);
		cout << "range" << r1 << " " << r2 << endl;
	}

	if (r1 < r2)
	{
		while (r1 < r2)
		{
			swap(vec[r1++], vec[r2--]);
		}
	}
	else
	{
		while (r1 > r2)
		{
			swap(vec[r1--], vec[r2++]);
		}
	}

	cout << "after" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << vec.at(i) << " ";
	}
	cout << endl;

	return 0;
}