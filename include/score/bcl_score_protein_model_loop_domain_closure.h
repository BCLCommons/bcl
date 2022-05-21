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

#ifndef BCL_SCORE_PROTEIN_MODEL_LOOP_DOMAIN_CLOSURE_H_
#define BCL_SCORE_PROTEIN_MODEL_LOOP_DOMAIN_CLOSURE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "coord/bcl_coord_cyclic_coordinate_descent.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelLoopDomainClosure
    //! @brief Scores how close the a loop is to being closed according to the psuedo residue, target residue, and atoms
    //! @details The distance between the atoms which are desired to be superimposed in order for a loop to be closed
    //!          is calculated for the atoms in the psuedo residue and in the target residue (where the loop is trying
    //!          to be closed towards).
    //!
    //! @see @link example_score_protein_model_loop_domain_closure.cpp @endlink
    //! @author alexanns, fischea
    //! @date Jan 4, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelLoopDomainClosure :
      public ProteinModel
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelLoopDomainClosure();

      //! @brief Clone function
      //! @return pointer to new ProteinModelLoopDomainClosure
      ProteinModelLoopDomainClosure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score after replicating the protein model using symmetry
      //! @param PROTEIN_MODEL protein model to be scored
      //! @return the score after replicating the protein model using symmetry
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief calculates the score for the agreement between the target and moving points
      //! @param TARGET_AND_MOVING_POINTS the target and moving points that will be scored
      //! @return double which is the score of the agreement between the target and moving points
      double CalculateScore
      (
        const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> &TARGET_AND_MOVING_POINTS
      ) const;

    }; // class ProteinModelLoopDomainClosure

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_LOOP_DOMAIN_CLOSURE_H_ 
