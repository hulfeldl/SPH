/*
 * Output.cpp
 *
 *  Created on: Aug 6, 2020
 *      Author: hulfeldl
 */




#include "Output.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>


namespace SPH {

template < uint8_t DIM, ParticleType Type >
Output<DIM,Type>::Output( 	const std::vector<Particle<DIM,Type>> & outParticles,
							const std::string & outDir,
							const Format FileFormat ) :
	m_outParticles( outParticles ),
	m_outDir( outDir ),
	m_FileFormat( FileFormat ) {


}

template < uint8_t DIM, ParticleType Type >
void Output<DIM,Type>::writeOutput( const std::string & Filename, const double t ) const {

	std::string OutFilename = Filename;
	while ( OutFilename.back() == '.' ){
		OutFilename.pop_back();
	}

	OutFilename += "_at_time_" + std::to_string(t) + "s";
	if ( !Filename.empty() ){

		switch( m_FileFormat ){

			case Format::Binary:
				writeBinaryOutput(OutFilename);
				break;

			case Format::Tecplot:
				writeTecplot(OutFilename,t);
				break;

			case Format::NONE:
				MYASSERT (	false,
							"No Fileformat specified!, Unable to wite output File!");
		}
	}
	else {
		MYASSERT (	false,
					"No Filename specified!, Unable to wite output File!");
	}
}

template < uint8_t DIM, ParticleType Type >
void Output<DIM,Type>::writeOutput(	const std::string & Filename,
									const double t,
									const Format Fileformat ){

	m_FileFormat	= Fileformat;

	writeOutput(Filename,t);
}

template < uint8_t DIM, ParticleType Type >
void Output<DIM,Type>::writeBinaryOutput( const std::string & Filename ) const {

	std::string FileExt	= ".bin";

	std::ofstream xv_file, state_file, other_file;

	xv_file.open(	m_outDir + Filename + "_xv"		+ FileExt, std::ios::out | std::ios::trunc | std::ios::binary	 );
	state_file.open(m_outDir + Filename + "_state"	+ FileExt, std::ios::out | std::ios::trunc | std::ios::binary	 );
	other_file.open(m_outDir + Filename + "_other"	+ FileExt, std::ios::out | std::ios::trunc | std::ios::binary	 );

	if( xv_file.is_open() && state_file.is_open() && other_file.is_open()  ) {
		int ntotal	= m_outParticles.size();

		const uint8_t dim = DIM;
		xv_file.write( 		reinterpret_cast<char*>( &ntotal ),	sizeof( ntotal ) );
		xv_file.write( 		reinterpret_cast<const char*>( &dim ),	sizeof( dim ) );

		state_file.write( 	reinterpret_cast<char*>( &ntotal ),	sizeof( ntotal ) );
		state_file.write( 	reinterpret_cast<const char*>( &dim ),		sizeof( dim ) );

		other_file.write( 	reinterpret_cast<char*>( &ntotal ),	sizeof( ntotal ) );
		other_file.write( 	reinterpret_cast<const char*>( &dim ),		sizeof( dim ) );

		for( int i = 0; i < ntotal; i++){

			const Particle<DIM,Type> & iParticle = m_outParticles[i];

			// write  x, xv
			const vDdf<DIM> x ( iParticle.m_x );
			const vDdf<DIM> v ( iParticle.m_v );
			xv_file.write( 		reinterpret_cast<char*>( &i ),		sizeof( i ) );
			xv_file.write( 		reinterpret_cast<const char*>( &x ),		sizeof( x ) );
			xv_file.write( 		reinterpret_cast<const char*>( &v ),		sizeof( v ) );

			// write mass, rho, p, u
			const double mass 	= iParticle.m_mass;
			const double rho 	= iParticle.m_rho;
			const double p 		= iParticle.m_p;
			const double u 		= iParticle.m_u;
			state_file.write(	reinterpret_cast<const char*>( &i ),		sizeof( i) );
			state_file.write( 	reinterpret_cast<const char*>( &mass ),	sizeof( mass ) );
			state_file.write( 	reinterpret_cast<const char*>( &rho ),	sizeof( rho ) );
			state_file.write( 	reinterpret_cast<const char*>( &p ),		sizeof( p ) );
			state_file.write( 	reinterpret_cast<const char*>( &u ),		sizeof( u ) );

			// write itype, hsml
			const int itype 	= static_cast<int> ( iParticle.m_fType );
			const double hsml	= iParticle.hsml();
			state_file.write( 	reinterpret_cast<const char*>( &i ),		sizeof( i ) );
			other_file.write( 	reinterpret_cast<const char*>( &itype ),	sizeof( itype ) );
			other_file.write( 	reinterpret_cast<const char*>( &hsml ),		sizeof( hsml ) );
		}

		xv_file.close();
		state_file.close();
		other_file.close();
	}

}

template < typename T >
void writeToFile( std::ofstream & file, const T out  ) {

	file.write( reinterpret_cast<const char*>( &out ),	sizeof( out ) );

}

template < uint8_t DIM, ParticleType Type >
void Output<DIM,Type>::writeTecplot( const std::string & Filename, const double simTime ) const {

	std::string FileExt	= ".plt";


	std::ofstream file;

	file.open(	m_outDir + Filename + FileExt, std::ios::out | std::ios::trunc | std::ios::binary	 );


	if( file.is_open()  ) {

	// I. HEADER SECTION

	// i. Magic number, Version number
	// --------------------------------------------------------------------------

	//	MAGIC = HEADER.('MAGIC')';

		std::string MAGIC = "#!TDV112";
		writeToFile( file, MAGIC  );
	//	MAGIC = fwrite(fileID,MAGIC,'char');

	// ii. Integer value of 1.
	//--------------------------------------------------------------------------
		writeToFile( file, int32_t(1)  );
	//	fwrite(fileID,1,'int32');

	// iii. Title and variable names.
	// --------------------------------------------------------------------------
	//	FILETYPE = HEADER.('FILETYPE');

		int32_t FILETYPE = 1;
		writeToFile( file, FILETYPE  );
	//	FILETYPE = fwrite(fileID,FILETYPE,'int32');

		std::string TITLE = "SPH Particle Collection" + '\0';
	//    HEADER.('TITLE_char') = [TITLE,char(0)];
	//    HEADER.('TITLE') = uint8(HEADER.('TITLE_char'));
	//	TITLE = HEADER.('TITLE');

		std::vector<int32_t> title ( TITLE.begin() , TITLE.end() );
		writeToFile( file, title );
	//	fwrite(fileID,TITLE,'int32');

		int32_t nParticles = m_outParticles.size();
		writeToFile( file, nParticles );
	//	fwrite(fileID,NVARS,'int32');

		std::vector<std::string> xStr;
		std::vector<std::string> vStr;
		for ( uint8_t i = 0; i < DIM; i++ ) {

			switch( DIM ){

			case 1:
				xStr.push_back("X");
				vStr.push_back("vx");
				break;

			case 2:
				xStr.push_back("Y");
				vStr.push_back("vy");
				break;

			case 3:
				xStr.push_back("Z");
				vStr.push_back("vz");
				break;

			}
		}

		std::vector<std::string> OtherStr = { "m", "rho", "p", "u","h", "FluidType" };

		std::vector<std::string> varNames = xStr;
		varNames.insert( varNames.end(), vStr.begin(), vStr.end() );
		varNames.insert( varNames.end(), OtherStr.begin(), OtherStr.end() );

		std::string VARNAMES;
		for ( size_t i = 0; i < varNames.size(); i++ ){
			VARNAMES += varNames[i] + '\0';;
		}

	//    HEADER.('VARNAMES_char') = VARNAMES_char;
	//    HEADER.('VARNAMES') = uint8(VARNAMES_char);
	//	VARNAMES = HEADER.('VARNAMES');

		std::vector<int32_t> varnames ( VARNAMES.begin() , VARNAMES.end() );
		writeToFile( file, varnames );
	//	fwrite(fileID,VARNAMES,'int32');


	// iv. Zones
	// --------------------------------------------------------------------------

		const float ZoneConst = 299.0;
		writeToFile( file, ZoneConst  );
	//	fwrite(fileID,299.0,'float');

		std::string ZONENAME = TITLE;
		std::vector<int32_t> zonename ( ZONENAME.begin() , ZONENAME.end() );
		writeToFile( file, zonename  );
	//	fwrite(fileID,ZONENAME,'int32');

		int32_t PARENTZONE = -1;
		writeToFile( file, PARENTZONE  );
	//	fwrite(fileID,PARENTZONE,'int32');

		int32_t STRANDID = 1;
		writeToFile( file, STRANDID  );
	//	fwrite(fileID,STRANDIDC,'int32');

		writeToFile( file, simTime  );

		writeToFile( file, int32_t(-1)  );
	//	fwrite(fileID,-1,'int32');

		int32_t ZONETYPE = 0; 		// ordered
		writeToFile( file, ZONETYPE );
	//	fwrite(fileID,ZONETYPE,'int32');

		int32_t VARLOCISSPEC = 1;
		writeToFile( file, VARLOCISSPEC  );
	//	fwrite(fileID,VARLOCISSPEC,'int32');

		if (VARLOCISSPEC == 1) {
			int32_t VARLOC = 0;
			for ( size_t i = 0; i < varNames.size(); i++ ){
				writeToFile( file, VARLOC  );
			}
		}

	// Two unused variables in the Spheres File
		writeToFile( file, int32_t(0)  );
		writeToFile( file, int32_t(0)  );
	//	fwrite(fileID,0,'int32');
	//	fwrite(fileID,0,'int32');

		if (ZONETYPE == 0){ 	//ordered
			int32_t IMAX = nParticles;
			int32_t JMAX = 1;
			int32_t KMAX = 1;
			writeToFile( file, IMAX );
			writeToFile( file, JMAX );
			writeToFile( file, KMAX );
	//		fwrite(fileID,IJKMAX,'int32');
		}

		int32_t AUXVALUES = 0;
		writeToFile( file, AUXVALUES );
	//	fwrite(fileID,AUXVALUES,'int32');

		// Write EoHmarker
	//	fwrite(fileID,357.0,'float');
		const float EoHmarker = 357.0;
		writeToFile( file, EoHmarker );

	// II. DATA SECTION


		// i. For both ordered and fe zoßnes:

		// Zone Marker
		writeToFile( file, float(299.0)  );
	//	fwrite(fileID,299.0,'float');

		// VARDATAFORMAT 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
		for ( size_t i = 0; i < DIM; i++ ){
			writeToFile( file, int32_t(2)  ); 	// x
			writeToFile( file, int32_t(2)  ); 	// v
		}

		writeToFile( file, int32_t(2)  ); 	// m
		writeToFile( file, int32_t(2)  ); 	// rho
		writeToFile( file, int32_t(2)  ); 	// p
		writeToFile( file, int32_t(2)  ); 	// u
		writeToFile( file, int32_t(2)  ); 	// h

		writeToFile( file, int32_t(5)  ); 	// FluidType
	//	fwrite(fileID,VARDATAFORMAT,'int32');

		// Write some other Stuff
	//	fwrite(fileID,0,'int32');
	//	fwrite(fileID,0,'int32');
	//	fwrite(fileID,-1,'int32');

		writeToFile( file, int32_t(0)  );
		writeToFile( file, int32_t(0)  );
		writeToFile( file, int32_t(-1)  );

		std::array<vNdf,DIM> x, v;
		vNdf m, rho, p, u, h;
		std::vector<uint8_t> fType;
		for ( const auto & iParticle : m_outParticles ) {

			for ( uint8_t d = 0; d < DIM; d++ ){
				x[d].push_back( iParticle.m_x[d] );
				v[d].push_back( iParticle.m_v[d] );
			}

			m.push_back( iParticle.m_mass );
			rho.push_back( iParticle.m_rho );
			p.push_back( iParticle.m_p );
			u.push_back( iParticle.m_u );
			h.push_back( iParticle.hsml() );

			fType.push_back( static_cast<uint8_t>(iParticle.m_fType) );


		}

		vNdf MINMAX;
		for ( uint8_t d = 0; d < DIM; d++ ){
			const auto [minX, maxX] = std::minmax_element( begin(x[d]), end(x[d]) );
			MINMAX.push_back(*minX);
			MINMAX.push_back(*maxX);
		}

		for ( uint8_t d = 0; d < DIM; d++ ){
			const auto [minV, maxV] = std::minmax_element( begin(v[d]), end(v[d]) );
			MINMAX.push_back(*minV);
			MINMAX.push_back(*maxV);
		}

		const auto [minM, maxM] = std::minmax_element( begin(m), end(m) );
		MINMAX.push_back(*minM);
		MINMAX.push_back(*maxM);

		const auto [minRho, maxRho] = std::minmax_element( begin(rho), end(rho) );
		MINMAX.push_back(*minRho);
		MINMAX.push_back(*maxRho);

		const auto [minP, maxP] = std::minmax_element( begin(p), end(p) );
		MINMAX.push_back(*minP);
		MINMAX.push_back(*maxP);

		const auto [minU, maxU] = std::minmax_element( begin(u), end(u) );
		MINMAX.push_back(*minU);
		MINMAX.push_back(*maxU);

		const auto [minH, maxH] = std::minmax_element( begin(h), end(h) );
		MINMAX.push_back(*minH);
		MINMAX.push_back(*maxH);

		const auto [minFT, maxFT] = std::minmax_element( begin(fType), end(fType) );
		MINMAX.push_back( static_cast<double>(*minFT) );
		MINMAX.push_back( static_cast<double>(*maxFT) );

		writeToFile( file, MINMAX  );
	//	fwrite(fileID,MINMAX,'float64');


		// Write data
		for ( uint8_t d = 0; d < DIM; d++ ){
			writeToFile( file, x[d] );
		}

		for ( uint8_t d = 0; d < DIM; d++ ){
			writeToFile( file, v[d] );
		}

		writeToFile( file, m );
		writeToFile( file, rho );
		writeToFile( file, p );
		writeToFile( file, u );
		writeToFile( file, h );

		writeToFile( file, fType );

	}
}


template < uint8_t DIM >
void writeTecplotHeader( std::ofstream & file, const double time ) {




}

template < uint8_t DIM, ParticleType Type >
Output<DIM,Type>::~Output() = default;


template class Output<1,ParticleType::Physical>;
template class Output<2,ParticleType::Physical>;
template class Output<3,ParticleType::Physical>;

template class Output<1,ParticleType::Virtual>;
template class Output<2,ParticleType::Virtual>;
template class Output<3,ParticleType::Virtual>;

} /* namespace SPH */






















