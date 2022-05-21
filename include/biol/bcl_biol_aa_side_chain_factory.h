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

#ifndef BCL_BIOL_AA_SIDE_CHAIN_FACTORY_H_
#define BCL_BIOL_AA_SIDE_CHAIN_FACTORY_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_table.h"

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASideChainFactory
    //! @brief class for adding idealized side chain atoms to Amino Acids
    //! @details This class adds idealized side chain atom coordinates to the given amino acids. It also gives choice
    //! whether to add hdyrogen coordinates as well. It should be noted that the given ProteinModel should have amino
    //! acids of type AABackBone or AAComplete and the returned ProteinModel will always have AAComplete amino acids.
    //!
    //! @see @link example_biol_aa_side_chain_factory.cpp @endlink
    //! @author rouvelgh, woetzen
    //! @date Feb 23, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASideChainFactory :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      bool m_IncludeHydrogen;          //!< whether or not to include hydrogen in side chains
      bool m_IncludeBackBone;          //!< whether or not to consider backbone atoms in side chain transformation
      std::string m_TableFilename;     //!< table file name
      storage::Table< double> m_Table; //!< table used to store side chain coordinates

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief filename of default atom coordinate table
      static const std::string &GetDefaultTableFileName();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AASideChainFactory();

      //! @param INCLUDE_HYDROGEN = true, if hydrogen atoms will be included in attached side chains
      //! @param INCLUDE_BACKBONE = use the backbone to align sidechain for first chi angle (important for proline)
      //! @param TABLE_FILE_NAME filename to table which contains ideal side chain conformations
      AASideChainFactory
      (
        const bool INCLUDE_HYDROGEN,
        const bool INCLUDE_BACKBONE = true,
        const std::string &TABLE_FILE_NAME = GetDefaultTableFileName()
      );

      //! @brief virtual copy constructor
      AASideChainFactory *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add side chain to AABase
      //! @param AA_BASE object to which side chain atoms should be added
      void AddSideChain( AABase &AA_BASE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief generate the sidechain atoms for a given aa type
      //! @param AATYPE the aatype that should have its idealized side chain coordinates retrieved
      //! @return ShPtrVector of atoms for the AAType of interest
      util::ShPtrVector< Atom> GenerateAtoms( const AAType &AATYPE) const;

      //! @brief takes a Protein model and returns a new protein model with side chains attached to new amino acids, sses and chains
      //! @param THIS_MODEL bcl protein model with no side chains
      //! @return ShPtr to new model with side chains added
      util::ShPtr< assemble::ProteinModel> ProteinModelWithSideChains( const assemble::ProteinModel &THIS_MODEL) const;

    private:

      //! @brief initialize table with idealized side chain coordinates
      void InitializeTable();

      //! @brief get transformation matrix to superimpose CA-CB vector of amino acid with the CA-CB vector of the atoms
      //! @param AMINO_ACID the AA considered for transformation
      //! @param ATOMS atoms contained in AA type
      //! @return TransformationMatrix3D that will superimpose the coordinates of ATOMS onto the CA-CB vector of the AA
      math::TransformationMatrix3D GetTransformationMatrix
      (
        const AABase &AMINO_ACID,
        const util::ShPtrVector< Atom> &ATOMS
      ) const;

      //! @brief get transformation matrix to superimpose CA-CB vector of amino acid with the CA-CB vector of the atoms
      //! @param AMINO_ACID the AA considered for transformation
      //! @param ATOMS atoms contained in AA type
      //! @return TransformationMatrix3D that will superimpose the coordinates of ATOMS onto the CA-CB vector of the AA
      math::TransformationMatrix3D GetQuadTransformationMatrix
      (
        const AABase &AMINO_ACID,
        const util::ShPtrVector< Atom> &ATOMS
      ) const;

      //! @brief takes a chain and returns a new chain with side chains attached to new amino acids
      //! @param THIS_CHAIN the chain which contains AAs with no side chains
      //! @return ShPtr to new chain with side chains on the AAs
      util::ShPtr< assemble::Chain> NewChain( const assemble::Chain &THIS_CHAIN) const;

      //! @brief takes an SSE and returns a new SSE with side chains attached to new amino acids
      //! @param THIS_SSE the SSE which contains AAs with no side chains
      //! @return ShPtr to new SSE with side chains on the AAs
      util::ShPtr< assemble::SSE> NewSSE( const assemble::SSE &THIS_SSE) const;

      //! @brief takes an AA and returns a new AA with side chain attached
      //! @param THIS_AA an AA with no side chain
      //! @return ShPtr to new AA with side chain attached
      util::ShPtr< AABase> NewAAWithSideChain( const AABase &THIS_AA) const;

      //! @brief generate atom from given atom type with its ideal coordinates
      //! @param ATOM_TYPE the atom used to lookup in the table
      //! @param ROW the row for a certain AA in the table
      //! @return ATOM which is of ATOM_TYPE and has idealized coordinates based on what AA it's in
      static Atom GenerateAtomWithCoords( const AtomType &ATOM_TYPE, const storage::Row< double> &ROW);

      //! @brief tansform a ShPtrVector of atoms
      //! @param TRANSFORMATION the transformation matrix used
      //! @param ATOMS the atoms to be transformed
      //! @return a ShPtrVector of transformed atoms
      static void TransformAtoms
      (
        const math::TransformationMatrix3D &TRANSFORMATION, util::ShPtrVector< Atom> &ATOMS
      );

    }; // class AASideChainFactory

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_SIDE_CHAIN_FACTORY_H_ 
