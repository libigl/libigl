// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <cassert>
#include <map>
#include <memory>
#include <typeindex>
#include <vector>

namespace igl
{

////////////////////////////////////////////////////////////////////////////////

struct AttributeBase {
private:
	std::type_index m_DerivedType;
public:
	AttributeBase(std::type_index t) : m_DerivedType(t) { }
	virtual ~AttributeBase() = default;
	std::type_index type() const { return m_DerivedType; }
};

// -----------------------------------------------------------------------------

template<typename T>
struct Attribute : public AttributeBase {
	// Constructor
	Attribute() : AttributeBase(typeid(T)) { }

	// Data
	T content_;
};

////////////////////////////////////////////////////////////////////////////////

// A generic class to manage attributes
class AttributeManager {

private:
	std::shared_ptr<AttributeBase> attr_;

public:
	// Retrieve custom attribute
	template<typename T> T & get();
	template<typename T> const T & get() const;

	// Retrieve attribute type
	std::type_index type() const { return attr_->type(); }
};

// -----------------------------------------------------------------------------

// Retrieve custom attribute
template<typename T>
inline T & AttributeManager::get()
{
	if (!attr_)
	{
		attr_ = std::make_shared<Attribute<T>>();
	}
	auto * derived = dynamic_cast<Attribute<T> *>(attr_.get());
	assert(derived && "Incompatible type requested for attribute");
	return derived->content_;
}


// Retrieve custom attribute
template<typename T>
inline const T & AttributeManager::get() const
{
	assert(attr_);
	const auto * derived = dynamic_cast<const Attribute<T> *>(attr_.get());
	assert(derived && "Incompatible type requested for attribute");
	return derived->content_;
}

} // namespace igl

