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

#ifndef BCL_SCORE_RESTRAINT_ATOM_DISTANCE_H_
#define BCL_SCORE_RESTRAINT_ATOM_DISTANCE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintAtomDistance
    //! @brief Scores a protein model and its distance restraints
    //! @details Scores a protein model according to the member scoring function.  The distance restraints must be
    //!          stored with the protein model as ProteinModelData.
    //!
    //! @see @link example_score_restraint_atom_distance.cpp @endlink
    //! @author weinerbe
    //! @date Mar 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintAtomDistance :
      public ProteinModel
    {

    public:

      // typedef for the type of information stored
      typedef util::ShPtrVector< restraint::AtomDistance> Container;

    private:

    //////////
    // data //
    //////////

      //! scoring function to use to evaulate the distances
      util::Implementation< RestraintAtomDistanceAssignment> m_ScoringFunction;

      //! type of distance restraints used
      util::ShPtr< Container> m_Restraints;

      //! fraction of total restraints to use for final score
      double m_FinalFraction;

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
      RestraintAtomDistance();

      //! @brief construct from a scoring function and scheme
      //! @param SCORING_FUNCTION scoring function to be used
      //! @param RESTRAINT_TYPE restraint type to be used
      //! @param FRACTION fraction of total restraints to be used for final score calculation
      //! @param SCHEME scheme to be used
      RestraintAtomDistance
      (
        const RestraintAtomDistanceAssignment &SCORING_FUNCTION,
        const double &FRACTION = 1.0,
        const std::string SCHEME = GetDefaultScheme(),
        const util::ShPtr< Container> &RESTRAINTS = util::ShPtr< Container>()
      );

      //! @brief Clone function
      //! @return pointer to new RestraintAtomDistance
      RestraintAtomDistance *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief gets the score for this function class
      //! @return the score for this function class
      const util::Implementation< RestraintAtomDistanceAssignment> &GetScore() const
      {
        return m_ScoringFunction;
      }

      //! @brief gets the restraints used by this restraint atom distance
      //! @return the restraints used by this score
      const Container &GetRestraints() const
      {
        return *m_Restraints;
      }

      //! @brief set the restraints used by this object
      //! @param RESTRAINTS the restraints used by this object
      void SetRestraints( const util::ShPtr< Container> &RESTRAINTS)
      {
        m_Restraints = RESTRAINTS;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator scores the protein model and associated restraints using the member scoring function
      //! @param PROTEIN_MODEL protein model to be scored
      //! @return distance restraint score
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param PROTEIN_MODEL protein model to be evaluated
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM
      ) const;

    }; // class RestraintAtomDistance

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_ATOM_DISTANCE_H_
