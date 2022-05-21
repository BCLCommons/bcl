// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "restraint/bcl_restraint_epr_decay_simulation.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically
#include <math.h>

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> EPRDecaySimulation::s_Instance
    (
      GetObjectInstances().AddInstance( new EPRDecaySimulation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EPRDecaySimulation::EPRDecaySimulation() :
      m_SpinLabelingPairs()
    {
    }

    //! @brief construct from list of spin-labeling pairs
    //! @param SL_PAIRS list of spin-labeling pairs
    EPRDecaySimulation::EPRDecaySimulation( const storage::Vector< SLPair> &SL_PAIRS) :
      m_SpinLabelingPairs( SL_PAIRS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new EPRDecaySimulation
    EPRDecaySimulation *EPRDecaySimulation::Clone() const
    {
      return new EPRDecaySimulation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &EPRDecaySimulation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &EPRDecaySimulation::GetAlias() const
    {
      static const std::string s_alias( "EPRDecaySimulation");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRDecaySimulation::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Simulates EPR decay patterns.");
      // serializer.AddInitializer
      // (
      //   "metropolis",
      //   "metropolis criterion to decide which mutates are accepted",
      //   io::Serialization::GetAgent( &m_Metropolis)
      // );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate decay for the provided distance and time
    //! @param DISTANCE spin-spin distance in angstroms
    //! @param TIME
    //! @return decay
    double EPRDecaySimulation::CalculateDecay( double DISTANCE, double TIME)
    {
      // calculate factors that are independent of the integration
      const double epr_constant( 326.08);
      const double inverse_distance( pow( DISTANCE / 10.0, 3));

      // integrate over the integration points
      double decay( 0.0);
      const size_t number_bins( 201);
      for( size_t bin_index( 0); bin_index < number_bins; ++bin_index)
      {
        const double fy( ( 1.0 - 3.0 * pow( 0.005 * bin_index, 2)) * epr_constant * inverse_distance);
        decay += cos( fy * TIME);
      }

      return decay;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
