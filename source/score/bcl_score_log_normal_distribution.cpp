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
#include "score/bcl_score_log_normal_distribution.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // initialize "s_WellDepth"
    const double LogNormalDistribution::s_WellDepth( -1.0);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LogNormalDistribution::s_Instance
    (
      util::Enumerated< RestraintAtomDistanceAssignment>::AddInstance( new LogNormalDistribution( 1.0))
    );

    // score for a restraint with residues/atoms not found in the protein model
    const double LogNormalDistribution::s_DefaultScore( 0.0);

    // effective distance per bond
    const double LogNormalDistribution::s_EffectiveDistancePerBond( 1.0);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LogNormalDistribution::LogNormalDistribution() :
      m_KVariable( util::GetUndefinedDouble()),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief parameter constructor
    //! @param SCHEME the short tag denoting this scoring function
    LogNormalDistribution::LogNormalDistribution
    (
      const double &KVARIABLE
    ) :
      m_KVariable( KVARIABLE),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief Clone function
    //! @return pointer to new LogNormalDistribution
    LogNormalDistribution *LogNormalDistribution::Clone() const
    {
      return new LogNormalDistribution( *this);
    }

    //! @brief virtual destructor
    LogNormalDistribution::~LogNormalDistribution()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LogNormalDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &LogNormalDistribution::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "log_normal_distribution");

      // end
      return s_default_scheme;
    }

    //! @brief returns scheme being used
    //! @return scheme being used
    const std::string &LogNormalDistribution::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator scores protein model
    //! @param RESTRAINT restraint to be scored
    //! @return score
    double LogNormalDistribution::operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const
    {
      // calculate the distance
      const double cb_distance( RESTRAINT.CalculateAtomDistance());

      // if the calculated distance is undefined
      if( !util::IsDefined( cb_distance))
      {
        // return default score
        return s_DefaultScore;
      }

      // get the bond distance
      const size_t bond_distance( GetTotalBondsFromCB( RESTRAINT));

      // if the bond distance is not defined
      if( !util::IsDefined( bond_distance))
      {
        // return default score
        return s_DefaultScore;
      }

      // calculate the value which will shift the graph enough that the scores will be negative within the upper bound
      // then add the negative well depth in order to shift the graph down and make it a bonus function
      return m_KVariable *
        std::pow
        (
          std::log
          (
            cb_distance / ( RESTRAINT.GetDistance() + double( bond_distance) * s_EffectiveDistancePerBond)
          ),
          2.0
        ) + s_WellDepth;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LogNormalDistribution::GetSerializer() const
    {
      io::Serializer serialize;
      serialize.SetClassDescription( "K * Ln^2(CB-Distance/(Restraint + {Restraint-label # bonds from CB})))  - 1");
      serialize.AddInitializer( "K", "magnitude for Log normal", io::Serialization::GetAgent( &m_KVariable), "1");
      return serialize;
    }

  } // namespace score
} // namespace bcl
