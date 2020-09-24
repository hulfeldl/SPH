/*
 * NeighbourSearch.cpp
 *
 *  Created on: Aug 31, 2020
 *      Author: hulfeldl
 */




#include <algorithm>

#include "NeighbourSearch.h"

namespace SPH {

template < 	uint8_t DIM,
ParticleType iType,
ParticleType jType >
bool AbstractSearch::isInteraction(	const Particle<DIM,iType> & iParticle,
											const Particle<DIM,jType> & jParticle,
											vDdf<DIM> & DX,
											double & r ){



	DX 	= dx(iParticle,jParticle);
	double rsq	= DX ^ DX;

	const double mhsml 		= fmax( iParticle.hsml() , jParticle.hsml() );
	const double scale_k	= fmax( iParticle.k() , jParticle.k() );

	if ( rsq < pow( scale_k * mhsml, 2.0 ) ){
		// Neighboring pair list, and totalinteraction number and
		// the interaction number for each particle
		r	= sqrt(rsq);
		return true;
	}
	else {
		return false;
	}
}

void AbstractSearch::print_interactions() const {

	// Statistics for the interaction
	int32_t maxp = -1, minp = -1;

	uint32_t maxiac 	= 0;
	uint32_t miniac 	= 1000;
	uint32_t noiac 	= 0;

	uint32_t i = 0;
	for( const auto iiac : m_niac ) {

		if ( iiac > maxiac ){
			maxiac = iiac;
			maxp = i;
		}

		if ( iiac < miniac ){
			miniac = iiac;
			minp = i;
		}

		if ( iiac == 0 ){
			noiac++;
		}
		i++;
	}

	std::cout << " >> Statistics: interactions per particle:" << std::endl;
	std::cout << "**** Particle		: " << maxp << " maximal interactions : " << maxiac << std::endl;
	std::cout << "**** Particle		: " << minp << " minimal interactions : " << miniac << std::endl;
	std::cout << "**** Average 		: " << static_cast<double>(2*m_niactot) / static_cast<double>(m_niac.size()) << std::endl;
	std::cout << "**** Total pairs 	: " << m_niactot << std::endl;
	std::cout << "**** Particles with no interactions : " << noiac << std::endl;

}

template < 	uint8_t DIM, ParticleType iType>
void DirectSearch<DIM,iType>::find(){

	size_t nParticles = this->m_iParticles.size();

	this->m_Interactions.clear();
	this->m_niactot = 0;
	this->m_niac = vNsu( nParticles , 0);

	for( size_t i = 0; i < nParticles - 1; i++ ) {

		const auto & iParticle = this->m_iParticles[i];

		for( size_t j = i + 1; j < nParticles; j++ ) {

			const auto & jParticle = this->m_iParticles[j];

			vDdf<DIM> dx;
			double r;
			if ( this->isInteraction(iParticle, jParticle, dx, r) ) {

				this->m_niac[i]++;
				this->m_niac[j]++;
				this->m_niactot++;

				auto maxInt = std::max_element( this->m_niac.begin(), this->m_niac.end() );
				MYASSERT ( 	*maxInt < this->m_maxInteractions,
							"Too many interactions!");

				this->m_Interactions.push_back( InteractionPair <DIM,iType,iType> (
											i, j, r, dx, iParticle, jParticle) );
			}

		}

	}

	this->print_interactions();

}

template < 	uint8_t DIM, ParticleType iType, ParticleType jType >
void DirectSearch<DIM,iType,jType>::find(){

	size_t niParticles = this->m_iParticles.size();
	size_t njParticles = this->m_jParticles.size();

	if ( niParticles == 0 || njParticles == 0){
		return;
	}

	this->m_Interactions.clear();
	this->m_niactot = 0;
	this->m_niac = vNsu( niParticles , 0);

	for( size_t i = 0; i < niParticles; i++ ) {

		const auto & iParticle = this->m_iParticles[i];

		for( size_t j = 0; j < njParticles; j++ ) {

			const auto & jParticle = this->m_jParticles[j];

			vDdf<DIM> dx;
			double r;
			if ( this->isInteraction(iParticle, jParticle, dx, r) ) {

				this->m_niac[i]++;
				this->m_niac[j]++;
				this->m_niactot++;

				auto maxInt = std::max_element( this->m_niac.begin(), this->m_niac.end() );
				MYASSERT ( 	*maxInt < this->m_maxInteractions,
							"Too many interactions!");

				this->m_Interactions.push_back( InteractionPair <DIM,iType,jType> (
											i, j, r, dx, iParticle, jParticle) );
			}

		}

	}

	this->print_interactions();
}


// Compile templates
// -----------------------------------------------------------------------------------
template class DirectSearch<1,ParticleType::Physical>;
template class DirectSearch<2,ParticleType::Physical>;
template class DirectSearch<3,ParticleType::Physical>;

template class DirectSearch<1,ParticleType::Physical,ParticleType::Virtual>;
template class DirectSearch<2,ParticleType::Physical,ParticleType::Virtual>;
template class DirectSearch<3,ParticleType::Physical,ParticleType::Virtual>;

} /* namespace SPH */


















