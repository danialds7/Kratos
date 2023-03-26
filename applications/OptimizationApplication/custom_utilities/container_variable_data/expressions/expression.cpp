//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <functional>
#include <numeric>

// Project includes

// Application includes

// Include base h
#include "expression.h"

namespace Kratos {

std::size_t Expression::GetLocalSize() const
{
    const auto& r_shape = this->GetShape();
    return std::transform_reduce(r_shape.begin(), r_shape.end(), 1UL, std::multiplies{},
                                 [](const auto v) { return v; });
}

void intrusive_ptr_add_ref(const Expression* x)
{
    x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
}

void intrusive_ptr_release(const Expression* x)
{
    if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
        std::atomic_thread_fence(std::memory_order_acquire);
        delete x;
    }
}

} // namespace Kratos