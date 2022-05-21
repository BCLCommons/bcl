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

#ifndef BCL_ASSEMBLE_COLLECTOR_COMMON_AA_H_
#define BCL_ASSEMBLE_COLLECTOR_COMMON_AA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorCommonAA
    //! @brief class collects all the amino acids that are common to two protein models.
    //!
    //! @see @link example_assemble_collector_common_aa.cpp @endlink
    //! @author alexanns, woetzen
    //! @date 03/27/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorCommonAA :
      public find::CollectorInterface< storage::VectorND< 2, util::SiPtrList< const biol::AABase> >, storage::VectorND< 2, ProteinModel> >
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CollectorCommonAA();

      //! @brief constructor that takes collector parameter
      //! @param COLLECTOR a collector object
      CollectorCommonAA( const CollectorCommonAA &COLLECTOR);

      //! @brief Clone is the virtual Clone constructor
      //! @return a pointer to new LocatorSSEFurthest which is a copy of this
      virtual CollectorCommonAA *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! Collect returns the amino acids which are common between "m_ProteinModel" and the ProteinModel argument
      //! @param PROTEIN_MODELS is the VectorND which holds the ProteinModels for which common amino acids are desired
      //! @return returns Group of the collected unpaired SSEs objects
      virtual
      storage::VectorND< 2, util::SiPtrList< const biol::AABase> >
      Collect( const storage::VectorND< 2, ProteinModel> &PROTEIN_MODELS) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief collect defined AAs in chains
      //! @param PROTEIN_MODEL model containing chains
      //! @return list of defined AAs
      static util::SiPtrList< const biol::AABase> CollectDefinedAAsInChains
      (
        const ProteinModel &PROTEIN_MODEL
      );

      //! @brief collect defined SSEs in chains
      //! @param PROTEIN_MODEL model containing SSEs
      //! @return list of defined SSEs
      static util::SiPtrList< const biol::AABase> CollectDefinedAAsInSSEs
      (
        const ProteinModel &PROTEIN_MODEL
      );

      //! @brief
      //! @param AAS_DEFINED ???
      //! @return
      static storage::VectorND< 2, util::SiPtrList< const biol::AABase> >
      Collect
      (
        const storage::VectorND< 2, util::SiPtr< const util::SiPtrList< const biol::AABase> > > &AAS_DEFINED
      );

      //! @brief CollectCommonCoordinates is for aligning residues and giving corresponding aligned lists of atoms
      //! @param AAS_DEFINED the two list of amino acids which will be aligned to get their coordinates
      //! @param ATOM_TYPES the atom types in the residues for which you want the coordinates
      //! @param NR_COMMON_RESIDUES number of common residues, this value will be overwritten by the function
      //! @return two lists of coordinates which are the aligned coordinates of the common aas and desired atoms
      static storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >
      CollectCommonCoordinates
      (
        const storage::VectorND< 2, util::SiPtr< const util::SiPtrList< const biol::AABase> > > &AAS_DEFINED,
        const storage::Set< biol::AtomType> &ATOM_TYPES,
        size_t &NR_COMMON_RESIDUES
      );

      //! @brief CollectCommonCoordinates is for aligning residues and giving corresponding aligned lists of atoms
      //! @param AAS_DEFINED_A the list of amino acids which will be aligned to get their coordinates from first model
      //! @param AAS_DEFINED_B the list of amino acids which will be aligned to get their coordinates from second model
      //! @param ATOM_TYPES the atom types in the residues for which you want the coordinates
      //! @param NR_COMMON_RESIDUES number of common residues, this value will be overwritten by the function
      //! @return two lists of coordinates which are the aligned coordinates of the common aas and desired atoms
      static storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >
      CollectCommonCoordinates
      (
        const util::SiPtrList< const biol::AABase> &AAS_DEFINED_A,
        const util::SiPtrList< const biol::AABase> &AAS_DEFINED_B,
        const storage::Set< biol::AtomType> &ATOM_TYPES,
        size_t &NR_COMMON_RESIDUES
      );

    }; // class CollectorSSEUnpaired

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_COLLECTOR_COMMON_AA_H_
