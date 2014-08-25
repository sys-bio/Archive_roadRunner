// == PREAMBLE ================================================

// * Licensed under the Apache License, Version 2.0; see README

/*
 * BaseTypes.cpp
 *
 *  Created on: Aug 21, 2014
 *
 *  Author: JKM
 */

// == INCLUDES ================================================

# include "BaseTypes.hpp"

// == CODE ====================================================


namespace rr
{

namespace rrgpu
{

namespace dom
{

void BaseType::serialize(Serializer& s) const {
    s << val_;
}

void PointerType::serialize(Serializer& s) const {
    s << *getRoot() << "*";
}

// std::string buildRep() {
//     assert(root_);
//
// }

BaseTypes::SelfPtr BaseTypes::self_;
std::once_flag BaseTypes::once_;

} // namespace dom

} // namespace rrgpu

} // namespace rr