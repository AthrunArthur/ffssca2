#ifndef UTILS_H
#define UTILS_H
#include <functional>

class scope_guard {
public:
    scope_guard(const scope_guard & s) = delete;
    const scope_guard & operator = (const scope_guard & ) = delete;

    template<class F1, class F2>
    scope_guard(const F1 & ctor, const F2 & dtor)
        : f(dtor) {
        ctor();
    }
    ~scope_guard() {
        f();
    }
private:
    std::function<void(void)> f;

};//end class scope_guard


#endif // UTILS_H
