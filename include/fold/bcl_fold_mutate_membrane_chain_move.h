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

#ifndef BCL_FOLD_MUTATE_MEMBRANE_CHAIN_MOVE_H_
#define BCL_FOLD_MUTATE_MEMBRANE_CHAIN_MOVE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateMembraneChainMove
    //! @brief Moves given chains of a membrane protein complex randomly and rotates the chains randomly
    //! @details This mutate is for docking membrane proteins in the membrane plane. It first moves given chains
    //!          (the ligand) of a membrane protein complex randomly on the membrane plane, relative to the receptor.
    //!          This is followed by random rotation of the ligand.
    //!
    //! @see @link example_fold_mutate_membrane_chain_move.cpp @endlink
    //! @author lib14
    //! @date Sep 26, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateMembraneChainMove :
        public math::MutateInterface< assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

    private:

      //! chain IDs of the receptor
      mutable std::string m_ReceptorChainIDs;

      //! chain IDs of the ligand
      mutable std::string m_LigandChainIDs;

      //! maximally allowed single-step translation along the x axis
      double m_MaxTranslationX;

      //! maximally allowed single-step translation along the y axis
      double m_MaxTranslationY;

      //! maximally allowed single-step translation along the z axis
      double m_MaxTranslationZ;

      //! maximally allowed single-step phi rotation angle
      double m_MaxPhi;

      //! maximally allowed single-step theta rotation angle
      double m_MaxTheta;

      //! maximally allowed single-step psi rotation angle
      double m_MaxPsi;

    public:

      //! single instance of this class
      static const util::SiPtr< util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateMembraneChainMove();

      //! @brief explicit constructor from given arguments
      //! @param RECEPTOR_CHAIN_IDS chain IDs of the receptor
      //! @param LIGAND_CHAIN_IDS chain IDs of the ligand
      //! @param MAX_TRANSLATION_X maximally allowed single-step translation along the x axis
      //! @param MAX_TRANSLATION_Y maximally allowed single-step translation along the y axis
      //! @param MAX_TRANSLATION_Z maximally allowed single-step translation along the y axis
      //! @param MAX_PHI maximally allowed single-step phi rotation angle
      //! @param MAX_THETA maximally allowed single-step theta rotation angle
      //! @param MAX_PSI maximally allowed single-step psi rotation angle
      MutateMembraneChainMove
      (
        const std::string &RECEPTOR_CHAIN_IDS,
        const std::string &LIGAND_CHAIN_IDS,
        double MAX_TRANSLATION_X = 1.0,
        double MAX_TRANSLATION_Y = 1.0,
        double MAX_TRANSLATION_Z = 1.0,
        double MAX_PHI = 1.0,
        double MAX_THETA = 1.0,
        double MAX_PSI = 1.0
      );

      //! @brief Clone function
      //! @return a pointer to a copy of this object
      MutateMembraneChainMove *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the chain IDs of the receptor
      //! @return the chain IDs of the receptor
      const std::string &GetReceptorChainIDs() const;

      //! @brief get the chain IDs of the ligand
      //! @return the chain IDs of the ligand
      const std::string &GetLigandChainIDs() const;

      //! @brief get the maximally allowed single-step translation along the x axis
      //! @return the maximally allowed single-step translation along the x axis
      double GetMaxTranslationX() const;

      //! @brief get the maximally allowed single-step translation along the y axis
      //! @return the maximally allowed single-step translation along the y axis
      double GetMaxTranslationY() const;

      //! @brief get the maximally allowed single-step translation along the z axis
      //! @return the maximally allowed single-step translation along the x axis
      double GetMaxTranslationZ() const;

      //! @brief get the maximally allowed single-step phi rotation angle
      //! @return the maximally allowed single-step phi rotation angle
      double GetMaxPhi() const;

      //! @brief get the maximally allowed single-step theta rotation angle
      //! @return the maximally allowed single-step theta rotation angle
      double GetMaxTheta() const;

      //! @brief get the maximally allowed single-step psi rotation angle
      //! @return the maximally allowed single-step psi rotation angle
      double GetMaxPsi() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief operator that takes a ProteinModel and return a mutated ProteinModel
      //! @param PROTEIN_MODEL ProteinModel which will be mutated
      //! @return MutateResult with the mutated ProteinModel
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set the chain IDs for the ligand
      void SetReceptorLigandChainIDs( const assemble::ProteinModel &MODEL) const;

    }; // end of class MutateMembraneChainMove
  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_MEMBRANE_CHAIN_MOVE_H_
