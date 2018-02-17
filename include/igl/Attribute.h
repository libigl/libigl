// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ATTRIBUTE_H
#define IGL_ATTRIBUTE_H

#include <typeindex>

namespace igl
{

struct AttributeBase
{
private:
	std::type_index derived_type_;
public:
	AttributeBase(std::type_index t) : derived_type_(t) { }
	virtual ~AttributeBase() = default;
	std::type_index type() const { return derived_type_; }
};

// -----------------------------------------------------------------------------

template<typename T>
struct Attribute : public AttributeBase
{
	// Constructor
	Attribute() : AttributeBase(typeid(T)) { }

	// Data
	T content_;
};

} // namespace igl

#endif // IGL_ATTRIBUTE_H
