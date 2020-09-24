/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * Dictionary.hpp
 *
 *  Created on: Nov 25, 2014
 *      Author: kustepha
 */

#ifndef SETTINGS_DICTIONARY_HPP_
#define SETTINGS_DICTIONARY_HPP_

namespace SPH {

namespace Settings {

  template<typename T>
  inline
  bool
  Dictionary::isType() const
  {
    return !isEmpty();	//isNumber();
  }


  template<>
  inline
  bool
  Dictionary::isType<bool>() const
  {
    return isBool();
  }

  template<>
  inline
  bool
  Dictionary::isType<std::string>() const
  {
    return isString();
  }

  template<>
  inline
  double
  Dictionary::get() const
  {
    return d_->GetDouble();
  }

  template<>
  inline
  bool
  Dictionary::get() const
  {
    return d_->GetBool();
  }

  template<>
  inline
  uint8_t
  Dictionary::get() const
  {
    return d_->GetUint();
  }

  template<>
  inline
  uint32_t
  Dictionary::get() const
  {
    return d_->GetUint();
  }

  template<>
  inline
  uint64_t
  Dictionary::get() const
  {
    return d_->GetUint64();
  }

  template<>
  inline
  int64_t
  Dictionary::get() const
  {
    return d_->GetInt64();
  }

  template<>
  inline
  std::string
  Dictionary::get() const
  {
    return d_->GetString();
  }

// Enum Gets
// ---------------------------------------------------------------------
template<>
inline
Format
Dictionary::get() const {
	std::string s = d_->GetString();

	if (  s == "binary" || s == "bin" ) {
		return Format::Binary;
	}
	else if ( s == "tecplot" || s == "plt" ) {
		return Format::Tecplot;
	}

	return Format::NONE;
}

template<>
inline
PrintStats
Dictionary::get() const {
	std::string s = d_->GetString();

	if (  s == "All" || s == "all" ) {
		return PrintStats::ALL;
	}
	else if ( s == "physical" || s == "Physical" ) {
		return PrintStats::sphParticles;
	}
	else if ( s == "virtual" || s == "Virtual" ) {
		return PrintStats::virtParticles;
	}

	return PrintStats::NONE;
}

template<>
inline
ParticleApproximation
Dictionary::get() const {
	std::string s = d_->GetString();

	if ( s == "Algorithm_1" ) {
		return ParticleApproximation::Algorithm_1;
	}
	else if ( s == "Algorithm_2" ) {
		return ParticleApproximation::Algorithm_2;
	}

	return ParticleApproximation::NONE;
}

template<>
inline
NearestNeighbourSearch
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Direct Search" ||
		 s == "DirectSearch") {
		return NearestNeighbourSearch::DirectSearch;
	}
	else if ( s == "Linked List" ||
			  s == "LinkedList") {
		return NearestNeighbourSearch::LinkedList;
	}
	else if ( s == "Tree Search" ||
			  s == "TreeSearch" ) {
		return NearestNeighbourSearch::TreeSearch;
	}

	return NearestNeighbourSearch::NONE;
}


template<>
inline
SmoothingLengthEvolution
Dictionary::get() const {
	std::string s = d_->GetString();

	if ( s == "Density Algebraic" ||
		 s == "DensityAlbebraic" ) {
		return SmoothingLengthEvolution::DensityAlgebraic;
	}
	else if ( s == "Smoothing Length ODE" ||
			  s == "SmoothingLengthODE" ) {
		return SmoothingLengthEvolution::SmoothingLengthODE;
	}
	else if ( s == "Other" ||
			  s == "OtherApproach" ||
			  s == "Other Approach" ) {
		return SmoothingLengthEvolution::OtherApproach;
	}

	return SmoothingLengthEvolution::NONE;

}

template<>
inline
SmoothingKernel
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "CubicSpline" ||
		 s == "Cubic Spline" ) {
		return SmoothingKernel::CubicSpline;
	}
	else if ( s == "Gauss" ) {
		return SmoothingKernel::Gauss;
	}
	else if ( s == "Quintic" ) {
		return SmoothingKernel::Quintic;
	}

	return SmoothingKernel::NONE;
}

template<>
inline
DensityCalculation
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Summation Density" ||
		 s == "SummationDensity" ) {
		return DensityCalculation::SummationDensity;
	}
	else if ( s == "CSPM" ) {
		return DensityCalculation::CSPM;
	}
	else if ( s == "Continuity" ||
			  s == "Continuity Equation" ||
			  s == "ContinuityEquation" ) {
		return DensityCalculation::ContinuityEquation;
	}

	return DensityCalculation::NONE;
}

template<>
inline
VelocityAveraging
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Monaghan" ) {
		return VelocityAveraging::Monaghan;
	}

	return VelocityAveraging::NONE;
}

template<>
inline
InitialConfiguration
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Generate" ) {
		return InitialConfiguration::Generate;
	}
	else if ( s == "Load" ) {
		return InitialConfiguration::Load;
	}

	return InitialConfiguration::NONE;
}

template<>
inline
ParticleType
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Physical" ) {
		return ParticleType::Physical;
	}
	else if ( s == "Virtual" ) {
		return ParticleType::Virtual;
	}

	return ParticleType::NONE;
}

template<>
inline
PDEs
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Euler Equations" ||
		 s == "EulerEquations" ) {
		return PDEs::EulerEquations;
	}
	else if ( s == "NavierStokesEquations" ||
			  s == "Navier Stokes Equations" ||
			  s == "NSE") {
		return PDEs::NavierStokesEquations;
	}

	return PDEs::NONE;
}

template<>
inline
ExternalForce
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Gravity" ) {
		return ExternalForce::Gravity;
	}
	else if ( s == "Boundary Force" ||
			  s == "Boundary" ||
			  s == "BondaryForce") {
		return ExternalForce::BoundaryForce;
	}

	return ExternalForce::NONE;
}

template<>
inline
Artificial
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Viscosity" ) {
		return Artificial::Viscosity;
	}
	else if ( s == "Heat" ) {
		return Artificial::Heat;
	}

	return Artificial::NONE;
}

template<>
inline
Symmetry
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "axis" || s == "Axis" ) {
		return Symmetry::Axis;
	}
	else if ( s == "center" || s == "Center" ) {
		return Symmetry::Center;
	}

	return Symmetry::NONE;
}

template<>
inline
ProblemType
Dictionary::get() const {

	std::string s = d_->GetString();

	if ( s == "Shock Tube" ) {
		return ProblemType::ShockTube;
	}
	else if ( s == "Shear Cavity" ) {
		return ProblemType::ShearCavity;
	}

	return ProblemType::NONE;
}

  // ---------------------------------------------------------------------

  template<typename T>
  inline
  std::vector<T>
  Dictionary::getAsVector() const
  {
    std::vector<T> ret;
    if (d_->IsArray())
      {
        ret.reserve(d_->Size());
        for (uint64_t i = 0; i < d_->Size(); ++i)
          ret.push_back( get<T>(i) );
      }
    else
      {
        MYASSERT( isType<T>(), "requested wrong type from dictionary");
        ret.push_back( get<T>() );
      }
    return ret;
  }

  template<typename T>
  inline
  std::vector<T> Dictionary::getAsVector(const std::string & key) const {
	    MYASSERT(hasMember(key),
	             std::string("couldn't find key \"") +
	             key +
	             std::string("\" in dict:\n") +
	             this->toString());
	    Dictionary d = get<Dictionary>(key);
	    return d.getAsVector<T>();
  }

  template<typename T>
  inline
  T
  Dictionary::get(const std::string& key) const
  {
    MYASSERT(hasMember(key),
             std::string("couldn't find key \"") +
             key +
             std::string("\" in dict:\n") +
             this->toString());
    Dictionary d = get<Dictionary>(key);
    return d.get<T>();
  }

  template<typename T>
  inline
  T
  Dictionary::get(const IndexType& key) const
  {
    Dictionary d = get<Dictionary>(key);
    MYASSERT(d.isType<T>(),"requested wrong type from dictionary");
    return d.get<T>();
  }

  template<>
  inline
  Dictionary
  Dictionary::get<Dictionary>(const std::string& key) const
  {
    MYASSERT(hasMember(key),
             std::string("couldn't find key \"") +
             key +
             std::string("\" in dict:\n") +
             this->toString());
    rapidjson::Value* raw = d_.get();
    return std::shared_ptr<rapidjson::Value>(d_, &((*raw)[key.c_str()]) );
  }

  template<>
  inline
  Dictionary
  Dictionary::get<Dictionary>(const IndexType& key) const
  {
    rapidjson::Value* raw = d_.get();
    return std::shared_ptr<rapidjson::Value>(d_, &((*raw)[key]) );
  }



  template<typename T>
  inline
  void
  Dictionary::get(
      T target[],
      const IndexType& first,
      const IndexType& n
  ) const
  {

    MYASSERT(
        first+n <= this->size(),
        std::string("requested ") + std::to_string(n) + std::string(" elements, starting at ")
    + std::to_string(first)
    + std::string(", but array only has ") + std::to_string(this->size())
    );
    T* t = target;

    for (SizeType i = first; i < first+n; ++i)
      {
        *(t++) = this->get<T>(i);
      }

  }

  template<typename T>
  inline
  void
  Dictionary::get(
      T target[],
      const std::string& key,
      const IndexType& first,
      const IndexType& n
  ) const
  {
    Dictionary dict = this->get<Dictionary>(key);
    dict.get(target,first,n);
  }


  template<typename T>
  inline
  T
  Dictionary::get(
      const std::string& key,
      const T& defaultVal)
  const
  {

    if ( this->isEmpty() || !(this->hasMember(key)) ) return defaultVal;
    return this->get<T>(key);

  }


  template<typename T>
  inline
  Dictionary
  Dictionary::getDictByValue(
      const std::string& key,
      const T& value )
  const
  {

    if (this->isArray())
      {

        for (DictIt it = this->begin(); it != this->end(); ++it)
          {
            Dictionary d( std::shared_ptr<rapidjson::Value>(d_, (rapidjson::Value*)it ) );
            if (d.get<T>(key) == value)
              return d;
          }
      }
    else
      {
        if (this->get<T>(key) == value)
          return *this;
      }


    MYASSERT(false,"dictionary with key \"" + key + "\" == " + std::string(value) + " not found!");

    return Dictionary();

  }


}/* namespace Settings */
}/* namespace SPH */

#endif /* SETTINGS_DICTIONARY_HPP_ */
