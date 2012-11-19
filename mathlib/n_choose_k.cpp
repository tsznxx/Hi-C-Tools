#include <stdexcept>
#include <limits>

using namespace std;

unsigned long long
gcd(unsigned long long x, unsigned long long y)
{
    while (y != 0)
    {
        unsigned long long t = x % y;
        x = y;
        y = t;
    }
    return x;
}

unsigned long long
choose(unsigned long long n, unsigned long long k)
{
    if (k > n)
        throw std::invalid_argument("invalid argument in choose");
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d, --n)
    {
        unsigned long long g = gcd(r, d);
        r /= g;
        unsigned long long t = n / (d / g);
        if (r > std::numeric_limits<unsigned long long>::max() / t)
           throw std::overflow_error("overflow in choose");
        r *= t;
    }
    return r;
}

int main()
{
    return 0;
}
