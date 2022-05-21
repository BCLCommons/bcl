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

#ifndef BCL_FOLD_PROTEIN_GEOMETRY_H_
#define BCL_FOLD_PROTEIN_GEOMETRY_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_loop_parameters.h"
#include "biol/bcl_biol_aa_base.h"
#include "linal/bcl_linal_vector_3d.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinGeometry
    //!
    //! @see @link example_fold_protein_geometry.cpp @endlink
    //! @author fischea
    //! @date Oct 12, 2015
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ProteinGeometry :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      ProteinGeometry();

      //! @brief clone function
      //! @return pointer to a new ProteinGeometry
      ProteinGeometry *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the local coordinate system of the given amino acid
      //! @param AMINO_ACID amino aid to return the local coordinate system for
      //! @return local coordinate system of the given amino acid
      static linal::Matrix3x3< double> GetLocalCoordinateSystem( const biol::AABase &AMINO_ACID);

      //! @brief fits the given loop to the given template and returns a new protein model containing the loop
      //! @param PROTEIN_MODEL protein model to containing the loop
      //! @param LOOP loop to be the given template
      //! @param TEMPLATE template used to fit the loop to
      //! @return new protein model containing the fitted loop
      static util::ShPtr< assemble::ProteinModel> FitToTemplate
      (
        const assemble::ProteinModel &PROTEIN_MODEL, const LoopParameters &LOOP, const LoopParameters &TEMPLATE
      );

      //! @brief combines the two given sequence at the given merging point
      //! @detail the given sequences are merged at the given merging points, which denote the residue index.
      //! merging points 5 and 3 results in the fifth residue of the n-terminal sequence being connected top the
      //! third residue of the c-terminal sequence
      //! @param SEQ_N n-terminal sequence to be merged
      //! @param SEQ_C c-terminal sequence to be merged
      //! @param MERGE_N merging point of the n-terminal sequence
      //! @param MERGE_C merging point of the c-terminal sequence
      static util::ShPtr< biol::AASequence> CombineSequences
      (
        const biol::AASequence &SEQ_N,
        const biol::AASequence &SEQ_C,
        int MERGE_N,
        int MERGE_C
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from a given input stream
      //! @param ISTREAM input stream to read members from
      //! @return input stream which members were read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into a given output stream
      //! @param OSTREAM output stream to write members into
      //! @param INDENT number of indentations
      //! @return output stream into which members were written
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ProteinGeometry

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PROTEIN_GEOMETRY_H_
