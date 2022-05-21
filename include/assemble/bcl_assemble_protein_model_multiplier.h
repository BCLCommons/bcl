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

#ifndef BCL_ASSEMBLE_PROTEIN_MODEL_MULTIPLIER_H_
#define BCL_ASSEMBLE_PROTEIN_MODEL_MULTIPLIER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_chain_multiplier.h"
#include "coord/bcl_coord_movable_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelMultiplier
    //! @brief generates a protein multimer
    //! @details () operator takes a monomeric protein model and returns a multimer defined by the member chain
    //!          multipliers.
    //!
    //! @see @link example_assemble_protein_model_multiplier.cpp @endlink
    //! @author weinerbe
    //! @date Nov 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelMultiplier :
      public math::FunctionInterfaceSerializable< ProteinModel, ProteinModel>,
      public coord::MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! chain multipliers containing information on how to handle each chain
      storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan> m_ChainMultipliers;

      //! the number of multimers
      size_t m_NumberMultimers;

      //! the orientation of the multiplier
      math::TransformationMatrix3D m_Orientation;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelMultiplier();

      //! @brief construct from triplets of chain ids and transformations
      //! @param TRANSFORMATIONS vector of current chain id, new chain id and transformation
      //! @param PROTEIN_MODEL original protein model
      //! @param CACHE bool whether to cache the SSE transformer
      ProteinModelMultiplier
      (
        const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > &TRANSFORMATIONS,
        const ProteinModel &PROTEIN_MODEL,
        const bool CACHE = false
      );

      //! @brief construct from an axis of symmetry and number of subunits
      //! @param SYMMETRY_AXIS axis of symmetry
      //! @param SUBUNITS number of subunits
      //! @param PROTEIN_MODEL original protein model
      //! @param CACHE bool whether to cache the SSE transformer
      ProteinModelMultiplier
      (
        const linal::Vector3D &SYMMETRY_AXIS,
        const size_t SUBUNITS,
        const ProteinModel &PROTEIN_MODEL,
        const bool CACHE = false
      );

      //! @brief construct from an axis of symmetry and number of subunits
      //! @param SYMMETRY_AXIS axis of symmetry
      //! @param SUBUNITS number of subunits
      //! @param PROTEIN_MODEL original protein model
      //! @param DIHEDRAL_AXIS optional secondary rotation axis for dihedral symmetry
      //! @param CACHE bool whether to cache the SSE transformer
      ProteinModelMultiplier
      (
        const linal::Vector3D &SYMMETRY_AXIS,
        const linal::Vector3D &DIHEDRAL_AXIS,
        const size_t SUBUNITS,
        const ProteinModel &PROTEIN_MODEL,
        const bool CACHE = false
      );

      //! @brief Clone function
      //! @return pointer to new ProteinModelMultiplier
      ProteinModelMultiplier *Clone() const;

      //! @brief make a copy of this multiplier that also copies the ShPtrs in the set of chain multipliers
      //! @return new copy of this
      util::ShPtr< ProteinModelMultiplier> HardCopy() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the chain multipliers
      //! @return the chain multipliers
      const storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan> &GetChainMultipliers() const
      {
        return m_ChainMultipliers;
      }

      //! @brief get number multimers
      //! @return number multimers
      size_t GetNumberMultimers() const
      {
        return m_NumberMultimers;
      }

      //! @brief return the orientation of the object
      //! @return orientation
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const;

      //! @brief return the orientation and Position as TransformationMatrix3D
      //! @return TransformationMatrix3D that defines orientation and position
      const math::TransformationMatrix3D GetOrientation() const
      {
        // return
        return m_Orientation;
      }

      //! @brief gets the transformation matrices for this multiplier
      //! @return the transformation matrices for this multiplier
      storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > GetTransformationMatrices() const;

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief transform the multimer
      //! @param TRANSFORMATION transformation to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION);

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION Translation to be applied
      void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief rotate the object by a given RotationMatrix3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

      //! @brief mapping of chain id mapping from original to the new chain id
      //! @return table containing the matrix number, the original chainid, the chainid it will have in the model
      storage::Table< char> ChainIDMapping() const;

      //! @brief gets a map of original chain id to target chain ids (as a string)
      //! @return map of original chain id to target chain ids (as a string)
      storage::Map< char, std::string> GetTargetChains() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief returns protein model after transforming all chains
      //! @param PROTEIN_MODEL protein model to be replicated
      //! @return protein model after transforming all chains
      ProteinModel operator ()( const ProteinModel &PROTEIN_MODEL) const;

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

    private:

      //! @brief wrap an SSE transformer into a FunctionCached object
      //! @param SSE_TRANSFORMER SSE transformer to be wrapped
      //! @param CACHE bool whether to cache the SSE transformer
      //! @return FunctionCached object
      static util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > > WrapCacheSSETransformer
      (
        const util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > > &SSE_TRANSFORMER,
        const bool CACHE
      );

      //! @brief converts axis, subunit, and chain id information into triplets used by the constructor
      //! @param SYMMETRY_AXIS axis of symmetry
      //! @param SUBUNITS number of subunits
      //! @param CHAIN_IDS chain ids present in original model
      //! @param DIHEDRAL_AXIS secondary rotation axis for dihedral symmetry
      //! @return transformation information for each chain multiplier
      static storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > GetChainData
      (
        const linal::Vector3D &SYMMETRY_AXIS,
        const size_t SUBUNITS,
        const std::string &CHAIN_IDS,
        const linal::Vector3D &DIHEDRAL_AXIS = linal::Vector3D()
      );

      //! @brief builds the member data from the transformations and protein model
      //! @param TRANSFORMATIONS vector of current chain id, new chain id and transformation
      //! @param CACHE bool whether to cache the SSE transformer
      //! @param PROTEIN_MODEL original protein model
      void BuildChainMultipliers
      (
        const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > &TRANSFORMATIONS,
        const ProteinModel &PROTEIN_MODEL,
        const bool CACHE
      );

    }; // class ProteinModelMultiplier

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PROTEIN_MODEL_MULTIPLIER_H_ 
