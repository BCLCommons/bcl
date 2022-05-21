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

#ifndef BCL_SCORE_ACCESSIBILITY_H_
#define BCL_SCORE_ACCESSIBILITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_protein_model_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Accessibility
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_score_accessibility.cpp @endlink
    //! @author alexanns
    //! @date Apr 14, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Accessibility :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! scoring function to use to evaulate the distances
      util::ShPtr< math::FunctionInterfaceSerializable< restraint::AccessibilityProfileAssignment, double> > m_ScoringFunction;

      //! the accessibility profile
      util::ShPtr< restraint::AccessibilityProfile> m_Restraints;

      //! scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Accessibility();

      //! @brief construct from a scoring function and scheme
      //! @param SCORING_FUNCTION scoring function to be used
      //! @param RESTRAINTS the actual restraints
      //! @param SCHEME scheme to be used
      Accessibility
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< restraint::AccessibilityProfileAssignment, double> > &SCORING_FUNCTION,
        const util::ShPtr< restraint::AccessibilityProfile> RESTRAINTS,
        const std::string SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new Accessibility
      Accessibility *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief () operator scores the protein model and associated restraints using the member scoring function
      //! @param PROTEIN_MODEL protein model to be scored
      //! @return distance restraint score
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

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

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param PROTEIN_MODEL Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      virtual std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM
      ) const;

    }; // class Accessibility

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_ACCESSIBILITY_H_
