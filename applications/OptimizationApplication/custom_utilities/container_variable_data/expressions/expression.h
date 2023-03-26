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

#pragma once

// System includes
#include <atomic>
#include <ostream>
#include <string>
#include <vector>

// Project includes
#include "includes/define.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the expression types.
 *
 * This is used to represent an expression for arithmatic evaluation
 *
 */
class Expression {
public:
    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<Expression>;

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    virtual ~Expression() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Evalute the expression for the given entity data start index and component index and returns the value
     *
     * @param EntityDataBeginIndex  Index at which entity data starts
     * @param ComponentIndex        Component index
     * @return double               Evaluated expression
     */
    virtual double Evaluate(
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const = 0;

    /**
     * @brief Get the Shape of the expression
     *
     * @return const std::vector<IndexType>     Size of each dimension is in the vector elements.
     */
    virtual const std::vector<IndexType> GetShape() const = 0;

    /**
     * @brief Get the Local Size of the expression
     *
     * @return IndexType
     */
    IndexType GetLocalSize() const;

    ///@}
    ///@name Input and output
    ///@{

    virtual std::string Info() const = 0;

    ///@}

private:
    ///@name Private operations
    ///@{

    //*********************************************
    // this block is needed for refcounting in th eintrusive ptr
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const Expression* x);

    friend void intrusive_ptr_release(const Expression* x);
    //*********************************************

    ///@}
};

/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const Expression& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos