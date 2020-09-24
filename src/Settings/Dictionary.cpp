/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * Dictionary.cpp
 *
 *  Created on: Aug 4, 2014
 *      Author: kustepha
 */

#include "Dictionary.h"

#include "rapidjson/prettywriter.h"

#include "rapidjson/error/en.h"

namespace SPH {

namespace Settings {

DictionarySource::DictionarySource ( ) {

}

DictionarySource::DictionarySource (const std::string& settingsDictFileName ) :
	m_dictFileName(settingsDictFileName),
	m_d(std::make_shared<rapidjson::Document>() ) {

	read();
}

void
DictionarySource::read() {

	// open the file
	std::ifstream docstrm(m_dictFileName);

	MYASSERT(	docstrm.is_open(),
				std::string("failed to open dictionary file \"") + m_dictFileName + std::string("\"")
	);

	// read the file
	std::string docstr((std::istreambuf_iterator<char>(docstrm)),
					   std::istreambuf_iterator<char>());

	// close the file
	docstrm.close();

	// parse
	m_d->Parse(docstr.c_str());

	if (m_d->HasParseError()) {

		rapidjson::ParseErrorCode pe = m_d->GetParseError();
		std::string msg(rapidjson::GetParseError_En(pe));

		int64_t errpos_min = std::max(int64_t(0),int64_t(m_d->GetErrorOffset()) - int64_t(25));
		int64_t errpos_max = std::min(int64_t(docstr.size()),int64_t(m_d->GetErrorOffset()) + int64_t(25));

		MYASSERT(false,
				 "failed to parse dictionary file \"" + m_dictFileName
				 + "\" at position " + std::to_string(m_d->GetErrorOffset()) + ":\n"
				 + msg + "\n\"...\n" + std::string( docstr.begin() + errpos_min, docstr.begin() + errpos_max) + "\n...\"" );
	}

}

void
DictionarySource::write(const std::string& settingsDictOutputFileName)
{

	// open the file
	std::ofstream docstrm(settingsDictOutputFileName);

	MYASSERT(	docstrm.is_open(),
				std::string("failed to open dictionary output file \"")
				+ settingsDictOutputFileName
				+ std::string("\"")
	);

	// to string
	rapidjson::StringBuffer buffer(0,1024);
	rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
	m_d->Accept(writer);

	// write
	docstrm.write(buffer.GetString(),buffer.GetSize());


	// close the file
	docstrm.close();

}

bool
Dictionary::hasMember(const std::string& key) const {
	return d_->HasMember(key.c_str());  }

	bool Dictionary::hasAllMembers(const std::vector<std::string>& keys) const {
	for (const auto& k : keys) { if(!hasMember(k)) return false; }
	return true;
}



std::string
Dictionary::toString() const {

	// to string
	rapidjson::StringBuffer buffer(0,1024);
	rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
	d_->Accept(writer);

	return {buffer.GetString(),buffer.GetSize()};

}

} /* namespace Settings */
} /* namespace SPH */





















