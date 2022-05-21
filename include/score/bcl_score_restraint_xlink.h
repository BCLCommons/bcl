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

#ifndef BCL_SCORE_RESTRAINT_XLINK_H_
#define BCL_SCORE_RESTRAINT_XLINK_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "bcl_score_radius_of_gyration.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_piecewise_function.h"
#include "restraint/bcl_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintXlink
    //! @brief This class scores the agreement of a protein model with cross-linking data
    //! @detail The agreement of a given protein model with given cross-linking data is scored by assuming that the
    //! cross-link's path is along the globular surface of the protein. This path is approximated via computing the arc
    //! of a sphere.
    //!
    //! @see @link example_score_restraint_xlink.cpp @endlink
    //! @author fischea
    //! @date June 12, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintXlink :
      public ProteinModel
    {

    //////////
    // data //
    //////////

    private:

      //! scheme of this score
      std::string m_Scheme;

      //! shared pointer to the observed cross-links with the used cross-linker lengths
      util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > m_Restraints;

      //! calculator for the radius of gyration of a protein model
      RadiusOfGyration m_ROGCalculator;

      //! The transition length used to call the scoring function
      double m_TransitionLength;

      //! shared pointer to the function that scores the agreement of the protein model with the cross-linking data
      util::ShPtr< math::PiecewiseFunction> m_ScoringFunction;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief construct from default values
      RestraintXlink();

      //! @brief construct from a restraint set and scheme
      //! @param SP_RESTRAINTS shared pointer to the restraints obtained from a cross-linking experiment
      //! @param TRANSITION_LENGTH length of the transition region of the scoring function
      //! @param SCHEME scheme of this score
      RestraintXlink
      (
        const util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > &SP_RESTRAINTS,
        double TRANSITION_LENGTH = 5.0,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief returns a pointer to a new RestraintXlink
      //! @return pointer to a new RestraintXlink
      RestraintXlink *Clone() const;

    /////////////////
    // data access //
    /////////////////
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the scheme of this score
      //! @return the scheme of this score
      const std::string &GetScheme() const;

      //! @brief returns the default scheme of this score
      //! @return the default scheme of this score
      static const std::string &GetDefaultScheme();

      //! @brief set the restraints
      //! @param RESTRAINTS vector of restraints
      void SetRestraints( const util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > &RESTRAINTS)
      {
        m_Restraints = RESTRAINTS;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    ////////////////
    // operations //
    ////////////////

      //! @brief scores the agreement of a given protein model with data obtained from a cross-linking experiment
      //! @param PROTEIN_MODEL protein model for which to compute the agreement
      //! @return agreement score normalized for each restraint with -1 being the best and 0 being the worst agreement
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from an input stream
      //! @param ISTREAM input stream to read members from
      //! @return the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into an output stream
      //! @param OSTREAM output stream to write members into
      //! @INDENT number of indentations to use
      //! @return the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief computes the agreement of the given protein model with the cross-link between the given atoms
      //! @param PROTEIN_MODEL protein model to compute the agreement for
      //! @param ATOM_FIRST first end point of the cross-link
      //! @param ATOM_SECOND second end point of the cross-link
      //! @param XL_LENGTH length of the cross-linker
      //! @return agreement of the protein model with the given cross-link
      double ComputeAgreement
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const biol::Atom &ATOM_FIRST,
        const biol::Atom &ATOM_SECOND,
        const double XL_LENGTH
      ) const;

      //! @brief creates the function that scores the agreement of a protein model with cross-linking data
      //! @param TRANSITION_LENGTH length of the transition region
      //! @return shared pointer to the scoring function
      static util::ShPtr< math::PiecewiseFunction> CreateScoringFunction( double TRANSITION_LENGTH);

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );
    }; // class RestraintXlink

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_XLINK_H_
