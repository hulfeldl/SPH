/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * Dictionary.h
 *
 *  Created on: Aug 4, 2014
 *      Author: kustepha
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <string>
#include <fstream>
#include <vector>
#include <cfloat>
#include <memory>

#include "rapidjson/document.h"

#include "Primitives/MyError.h"
#include "Primitives/DataTypes/MyTypes.h"

namespace SPH {

namespace Settings {

  struct DictionarySource {

	const std::string m_dictFileName;

	std::shared_ptr<rapidjson::Document> m_d;

	void read();

	void write(const std::string& settingsDictOutputFileName);

	DictionarySource ( );

	DictionarySource (const std::string& settingsDictFileName);

  };


  class Dictionary {

  public:

	mutable std::shared_ptr<rapidjson::Value> d_ = nullptr;

	typedef rapidjson::Value::ConstValueIterator	DictIt;

	DictIt begin() const { return d_->Begin();	}

	DictIt end() const { return d_->End(); }

	typedef rapidjson::SizeType 	SizeType;
	typedef uint64_t				IndexType;

	template<typename T>
	T get() const;

	template<typename T>
	T get(const std::string& key) const;

	template<typename T>
	T get(const IndexType& key) const;

	template<typename T>
	void get(T target[], const IndexType& first, const IndexType& n = 1) const;

	template<typename T>
	void get(T target[], const std::string& key, const IndexType& first, const IndexType& n = 1) const;

	template<typename T>
	T get(const std::string& key, const T& defaultVal) const;

    template<typename T>
    std::vector<T> getAsVector() const;

    template<typename T>
    std::vector<T> getAsVector(const std::string& key) const;

	template<typename T>
	Dictionary getDictByValue(const std::string& key, const T& value) const;

	bool isEmpty() const {  return d_ == nullptr; }

	bool hasMember(const std::string& key) const;
	bool hasAllMembers(const std::vector<std::string>& keys) const;

	bool isArray() const {
	  return !isEmpty() && d_->IsArray(); }

	bool isObject() const {
	  return !isEmpty() && d_->IsObject(); }

	bool isValue() const {
	  return !isEmpty() && !isObject() && !isArray(); }

	bool isBool() const	{
	  return !isEmpty() && d_->IsBool(); }

	bool isNumber() const {
	  return !isEmpty() && d_->IsNumber(); }

	bool isString() const {
	  return !isEmpty() && d_->IsString(); }

	template<typename T>
	bool isType() const;

	SizeType size() const { return d_->Size(); }

	std::string
	toString()
	const;

	Dictionary ()
	: d_(nullptr)
	{}

	Dictionary (const Dictionary& other)
	: d_(other.d_)
	{}

	Dictionary (const std::shared_ptr<rapidjson::Value> d)
	: d_(d)
	{}

	Dictionary (const DictionarySource& dictSrc)
	: d_(dictSrc.m_d)
	{}


  };

} /* namespace Settings */
} /* namespace SPH */

// IMPLEMENTATION

#include "Dictionary.hpp"

#endif /* DICTIONARY_H_ */
