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
#include "score/bcl_score_fuzzy_logic_filter.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_gaussian_function.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FuzzyLogicFilter::s_Instance
    (
      util::Enumerated< RestraintAtomDistanceAssignment>::AddInstance( new FuzzyLogicFilter())
    );

    //! score for a restraint with residues/atoms not found in the protein model
    const double FuzzyLogicFilter::s_DefaultScore( 0.0);

    //! effective distance per bond
    const double FuzzyLogicFilter::s_EffectiveDistancePerBond( 1.0);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief parameter constructor
    //! @param SCHEME the short tag denoting this scoring function
    FuzzyLogicFilter::FuzzyLogicFilter( const std::string &SCHEME) :
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new class_name
    FuzzyLogicFilter *FuzzyLogicFilter::Clone() const
    {
      return new FuzzyLogicFilter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FuzzyLogicFilter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &FuzzyLogicFilter::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "fuzzy_logic_filter");

      // end
      return s_default_scheme;
    }

    //! @brief returns scheme being used
    //! @return scheme being used
    const std::string &FuzzyLogicFilter::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator scores protein model
    //! @param RESTRAINT restraint to be scored
    //! @return score
    double FuzzyLogicFilter::operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const
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

      // create a gaussian function to evaluate the difference between the distances
      return -math::GaussianFunction
      (
        0.0,
        std::max( 0.5, RESTRAINT.GetUpperBound() - RESTRAINT.GetDistance())
      )( cb_distance - RESTRAINT.GetDistance() - double( bond_distance) * s_EffectiveDistancePerBond);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FuzzyLogicFilter::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FuzzyLogicFilter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
