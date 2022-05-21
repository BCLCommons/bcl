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

#ifndef BCL_SCORE_FUZZY_LOGIC_FILTER_H_
#define BCL_SCORE_FUZZY_LOGIC_FILTER_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_restraint_nmr_distance_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FuzzyLogicFilter
    //! @brief calculates a score for distance restraints to be used within the bcl
    //! @details uses a Gaussian Function to determine score allowing the equation to be differentiable
    //!        y = -1 * exp( ( -( X - 0.0) ^ 2) / ( 2 * (std_deviation) ^ 2))
    //! when the distance is greater than the maximum tolerance the score will be 0
    //! when the distance is within the standard deviation, the score will be -1
    //! a linear function will be applied for scores between the standard deviation and max tolerance
    //! Please see Meiler et. al. PMAS 2003 "Rapid protein fold determination using unassigned NMR data"
    //!
    //! @see @link example_score_fuzzy_logic_filter.cpp @endlink
    //! @author akinlr
    //! @date Jul 8, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FuzzyLogicFilter :
      public RestraintNMRDistanceInterface
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used
      std::string m_Scheme;

      //! score for a restraint with residues/atoms not found in the protein model
      static const double s_DefaultScore;

      //! effective distance per bond
      static const double s_EffectiveDistancePerBond;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief parameter constructor
      //! @param SCHEME the short tag denoting this scoring function
      FuzzyLogicFilter( const std::string &SCHEME = GetDefaultScheme());

      //! @brief Clone function
      //! @return pointer to new class_name
      virtual FuzzyLogicFilter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator scores protein model
      //! @param RESTRAINT restraint to be scored
      //! @return score
      double operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class score_fuzzy_logic_filter

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_FUZZY_LOGIC_FILTER_H_ 
