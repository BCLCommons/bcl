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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_CHAIN_MOVE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_CHAIN_MOVE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    /////////////////////////////////////////////
    //! @class MutateProteinModelChainMove
    //! @brief A mutate for protein-protein docking.
    //! @detail To be completed
    //!
    //! @author lib14
    //! @modified April 20, 2018
    //!
    /////////////////////////////////////////////

    class BCL_API MutateProteinModelChainMove :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! maximum translation in the X and Y direction
      double m_MaxTranslationXY;

      //!
      double m_MaxTranslationY;

      //! maximum translation in the Z direction
      double m_MaxTranslationZ;

      //! maximum rotation angle in randians
      double m_MaxRotationAngle;

      // euler angles
      double m_MaxPhi;

      //
      double m_MaxTheta;

      //
      double m_MaxPsi;

      //!
      double m_MaxInternalRotationAngle;

    public:

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      MutateProteinModelChainMove();

      //! @param MAX_TRANSLATION_XY
      //! @param MAX_TRANSLATION_Z
      //! @param MAX_ROTATION_ANGLE
      //! @param SCHEME
      MutateProteinModelChainMove
      (
        const double MAX_TRANSLATION_XY,
        const double MAX_TRANSLATION_Z,
        const double MAX_ROTATION_ANGLE,
        const double MAX_INTERNAL_ROTATION_ANGLE
      );

      //! @brief Clone function
      //! @return a pointer to a copy of this object
      MutateProteinModelChainMove *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief Get the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      static char GetLigandChainID( const assemble::ProteinModel &COMPLEX);

      //! @brief Generate placement for the given Chain at a random orientation with respect to a located Chain in the model
      //! @param MOVABLE_CHAIN chain to be moved
      //! @param PROTEIN_MODEL the given protein model to which the chain is to be added
      //! @return
      storage::Pair< math::TransformationMatrix3D, bool> PlaceLigand
      (
        const assemble::Chain &LIGAND,
        const assemble::ProteinModel &COMPLEX
      ) const;

      //! @brief
      //! @param COMPLEX
      static void PreDock( assemble::ProteinModel &COMPLEX, bool MEMBRANE = false);

      //! @param PROTEIN_MODEL
      //! @return
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //!
      //!
      //! @param COMPLEX
      //! @return
      //!
      static const util::ShPtr< assemble::Chain> &GetLigand( const assemble::ProteinModel &COMPLEX);

      //!
      //!
      //! @param COMPLEX
      //! @return
      //!
      static util::ShPtrVector< assemble::Chain> GetReceptor( const assemble::ProteinModel &COMPLEX);

      //! @brief computes the shortest distance between the ligand chain and all receptor chain(s)
      //! @param LIGAND the ligand chain
      //! @param RECEPTOR the receptor chain(s)
      //! @return
      static std::pair< double, linal::Vector3D> ComputeReceptorLigandShortestDistance
      (
        const assemble::ProteinModel &RECEPTOR,
        const assemble::ProteinModel &LIGAND
      );
    };
  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_CHAIN_MOVE_H_
