/*
 * NeighbourSearch.h
 *
 *  Created on: Aug 31, 2020
 *      Author: hulfeldl
 */




#ifndef PARTICLES_NEIGHBOURSEARCH_H_
#define PARTICLES_NEIGHBOURSEARCH_H_

#include <vector>
#include <memory>

#include "Primitives/DataTypes/MyTypes.h"
#include "Particles/Particle.h"
#include "Particles/InteractionPair.h"

namespace SPH {


class AbstractSearch {

protected:

	uint32_t m_maxInteractions;
	uint32_t m_niactot;
	vNsu 	 m_niac;

public:

	AbstractSearch( double maxInteractions ) :
						m_maxInteractions(maxInteractions),
						m_niactot(),
						m_niac() {}

	virtual ~AbstractSearch() = default;

	template < 	uint8_t DIM,
	ParticleType iType,
	ParticleType jType >
	bool isInteraction(	const Particle<DIM,iType> & iParticle,
						const Particle<DIM,jType> & jParticle,
						vDdf<DIM> & dx,
						double & r );

	void print_interactions() const;

	virtual void find() = 0;



};


template < uint8_t DIM, ParticleType... Types >
class NeighbourSearch : public AbstractSearch {

};

template < 	uint8_t DIM,
			ParticleType iType>
class NeighbourSearch<DIM,iType> : public AbstractSearch   {

	typedef Particle<DIM,iType> 				Particle_t;
	typedef std::vector<Particle_t>				ParticleVector;

	typedef InteractionPair<DIM,iType,iType> 	Interaction_t;
	typedef std::vector<Interaction_t> 			InteractionVector;

protected:

	ParticleVector & m_iParticles;

	InteractionVector & m_Interactions;

public:

	NeighbourSearch();
	NeighbourSearch(ParticleVector & Particles,
					InteractionVector & Interactions,
					double maxInteractions ) :
						AbstractSearch(maxInteractions),
						m_iParticles(Particles),
						m_Interactions(Interactions) {}

	virtual ~NeighbourSearch() = default;
	NeighbourSearch(const NeighbourSearch &other);
	NeighbourSearch(NeighbourSearch &&other);
	NeighbourSearch& operator=(const NeighbourSearch &other);
	NeighbourSearch& operator=(NeighbourSearch &&other);

	virtual void find() = 0;
};


template < 	uint8_t DIM,
			ParticleType iType,
			ParticleType jType >
class NeighbourSearch<DIM,iType,jType> : public AbstractSearch  {

	template < ParticleType Type >
	using Particle_t = Particle<DIM,Type>;

	template < ParticleType Type >
	using ParticleVector = std::vector<Particle_t<Type>>;

	typedef InteractionPair<DIM,iType,jType> 	Interaction_t;
	typedef std::vector<Interaction_t> 		InteractionVector;

protected:

	ParticleVector<iType> & m_iParticles;
	ParticleVector<jType> & m_jParticles;

	InteractionVector & m_Interactions;

public:

	NeighbourSearch();
	NeighbourSearch(ParticleVector<iType> & iParticles,
					ParticleVector<jType> & jParticles,
					InteractionVector & Interactions,
					double maxInteractions ) :
						AbstractSearch(maxInteractions),
						m_iParticles(iParticles),
						m_jParticles(jParticles),
						m_Interactions(Interactions) {}

	virtual ~NeighbourSearch() = default;
	NeighbourSearch(const NeighbourSearch &other);
	NeighbourSearch(NeighbourSearch &&other);
	NeighbourSearch& operator=(const NeighbourSearch &other);
	NeighbourSearch& operator=(NeighbourSearch &&other);

	virtual void find() = 0;
};


// Direct Search
// -----------------------------------------------------------------------------------
template < uint8_t DIM, ParticleType... Types >
class DirectSearch {

};


template < uint8_t DIM, ParticleType iType >
class DirectSearch<DIM,iType> : public NeighbourSearch<DIM,iType> {

	typedef Particle<DIM,iType> 				Particle_t;
	typedef std::vector<Particle_t>				ParticleVector;

	typedef InteractionPair<DIM,iType,iType> 	Interaction_t;
	typedef std::vector<Interaction_t> 			InteractionVector;

public:

	DirectSearch();
	DirectSearch( 	ParticleVector & Particles,
					InteractionVector & Interactions,
					double maxInteractions ) :
		NeighbourSearch<DIM,iType>( Particles, Interactions, maxInteractions) {}

	virtual ~DirectSearch() = default;
	DirectSearch(const DirectSearch &other);
	DirectSearch(DirectSearch &&other);

	void find() override;
};


template < 	uint8_t DIM, ParticleType iType, ParticleType jType >
class DirectSearch<DIM,iType,jType> : public NeighbourSearch<DIM,iType,jType> {

	template < ParticleType Type >
	using Particle_t = Particle<DIM,Type>;

	template < ParticleType Type >
	using ParticleVector = std::vector<Particle_t<Type>>;

	typedef InteractionPair<DIM,iType,jType> 	Interaction_t;
	typedef std::vector<Interaction_t> 		InteractionVector;

public:

	DirectSearch();
	DirectSearch(ParticleVector<iType> & iParticles,
				ParticleVector<jType> & jParticles,
				InteractionVector & Interactions,
				double maxInteractions ) :
						NeighbourSearch<DIM,iType,jType>(iParticles,jParticles,Interactions,maxInteractions) {}

	virtual ~DirectSearch() = default;

	void find() override;

};

} /* namespace SPH */

#endif /* PARTICLES_NEIGHBOURSEARCH_H_ */














