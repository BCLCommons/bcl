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

#ifndef BCL_CHEMISTRY_MERGE_FRAGMENT_COMPLETE_H_
#define BCL_CHEMISTRY_MERGE_FRAGMENT_COMPLETE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MergeFragmentComplete
    //! @brief Merges two given fragments based on the atoms that need to be overlapping or atoms that need to be joined
    //!
    //! @details Given two fragments, the class returns one fragment according to atoms that are specified to be connected
    //!          or overlapped
    //! @see @link example_chemistry_merge_fragment_complete.cpp @endlink
    //! @author kothiwsk, geanesar
    //! @date Oct 08, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MergeFragmentComplete :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new MergeFragmentComplete
      MergeFragmentComplete *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief merge two fragment given common vertices for the fragments
      //! @param MOLECULE_A small Molecule to be merged with MOLECULE_B
      //! @param MOLECULE_B small Molecule to be merged with MOLECULE_A
      //! @param COMMON_INDICES atoms that are common between MOLECULE_A (keys) and MOLECULE_B(mapped values)
      //! @param APPENDED_ATOMS container to store atoms of molecule B that have been appended to molecule A
      //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
      static storage::Pair< bool, FragmentComplete> MergeFragments
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::Map< size_t, size_t> &COMMON_INDICES,
        storage::Vector< size_t> &APPENDED_ATOMS
      );

      //! @brief merge two fragment given common vertices for the fragments
      //! @param MOLECULE_A small Molecule to be merged with MOLECULE_B assembled
      //! @param MOLECULE_B small Molecule to be merged with MOLECULE_A assembled
      //! @param COMMON_INDICES atoms that are common between MOLECULE_A (keys) and MOLECULE_B(mapped values)
      //! @param NEIGHBOR_ATOM atom in molecule B which is bonded to an atom in molecule A (which is not in B)
      //! @param EDGE_DATA edge data between neighbor atom and an atom in A to which it is bonded
      //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
      static storage::Pair< bool, FragmentComplete> MergeFragmentParts
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::Map< size_t, size_t> &COMMON_INDICES,
        const AtomComplete &NEIGHBOR_ATOM,
        const size_t EDGE_DATA,
        const bool PERPENDICULAR_RING
      );

      //! @brief connect two fragment given connection points for the fragments
      //! @param MOLECULE_A small Molecule to be connected with MOLECULE_B assembled
      //! @param MOLECULE_B small Molecule to be connected with MOLECULE_A assembled
      //! @param BOND_TYPE the bond that is to be used to connect
      //! @param INDICES_TO_CONNECT atoms of MOLECULE A (keys) that need to be connected to atoms of MOLECULE B(mapped values)
      //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
      static storage::Pair< bool, FragmentComplete> MergeFragments
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const ConfigurationalBondType &BOND_TYPE,
        const storage::Pair< size_t, size_t> &INDICES_TO_CONNECT
      );

//      //! @brief connect two fragment given connection points for the fragments
//      //! @param FUSED_RING Ring which needs to built after merging constituent rings
//      //! @param CONSTITUENT_RINGS constituent rings and their isomorphisms to the fused ring
//      //! @return a pair of bool and fused ring. Bool is true if merge was successful.
//      static storage::Pair< bool, FragmentComplete> MergeRings
//      (
//        const FragmentComplete &FUSED_RING,
//        const storage::List< storage::Pair< FragmentComplete, storage::Vector< size_t> > > &CONSTITUENT_RINGS
//      );

      //! @brief connect two fragment given multiple connection points for the fragments
      //! @param MOLECULE_A small Molecule to be connected with MOLECULE_B assembled
      //! @param MOLECULE_B small Molecule to be connected with MOLECULE_A assembled
      //! @param BOND_TYPES list of bond types used to connect atoms together
      //! @param INDICES_TO_CONNECT atoms of MOLECULE A (keys) that need to be connected to atoms of MOLECULE B(mapped values)
      //! @return a pair of bool and merged framgents. Bool is true if merge was successful.
      static storage::Pair< bool, FragmentComplete> MultipointMergeFragments
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::List< ConfigurationalBondType> &BOND_TYPES,
        const storage::List< storage::Pair< size_t, size_t> > &INDICES_TO_CONNECT
      );

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

    private:

      //! @brief transformed coordinates of molecule B such that it can be merged with molecule A
      //! @param MOLECULE_A small molecule to be merged with MOLECULE_B assembled
      //! @param MOLECULE_B small molecule to be merged with MOLECULE_A assembled
      //! @param MOLECULE_A_INDICES atoms of molecule A that are common with those in molecule B
      //! @param MOLECULE_B_INDICES atoms of molecule B that are common with those in molecule A
      //! @param NEIGHBOR_ATOM atom in molecule B which is bonded to an atom in molecule A (which is not in B)
      //! @param EDGE_DATA edge data between neighbor atom and an atom in A to which it is bonded
      //! @return transformed coordinates of molecule B
      static storage::Vector< linal::Vector3D> GetTransformedCoordinates
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::Vector< size_t> MOLECULE_A_INDICES,
        const storage::Vector< size_t> MOLECULE_B_INDICES,
        const AtomComplete &NEIGHBOR_ATOM,
        const size_t EDGE_DATA,
        const bool PERPENDICULAR_RING
      );

      //! @brief transformed coordinates of molecule B such that it can be connected with molecule A
      //! @param MOLECULE_A small molecule to be merged with MOLECULE_B assembled
      //! @param MOLECULE_B small molecule to be merged with MOLECULE_A assembled
      //! @param BOND_TYPE bond type to be used to connect the molecules
      //! @param VERTICES_TO_CONNECT atoms of MOLECULE A (key) that need to be connected to atoms of MOLECULE A (mapped value)
      //! @return transformed coordinates of molecule B
      static storage::Vector< linal::Vector3D> GetTransformedCoordinates
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const ConfigurationalBondType &BOND_TYPE,
        const storage::Pair< size_t, size_t> &VERTICES_TO_CONNECT
      );

      //! @brief transformed coordinates after application of a given transformation matrix to a given set of coordinates
      //! @param COORDINATES coordinates that need to be transformed
      //! @param TRANSFORMATION_MATRIX transformation matrix that needs to be applied to coordinates of interest
      //! @return transformed coordinates
      static storage::Vector< linal::Vector3D> TransformCoordinates
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const math::TransformationMatrix3D &TRANSFORMATION_MATRIX
      );

      //! @brief transformed coordinates after application of a given transformation matrix to a given set of coordinates
      //! @param COORDINATES coordinates that need to be transformed
      //! @param TRANSFORMATION_MATRIX transformation matrix that needs to be applied to coordinates of interest
      //! @return transformed coordinates
      static storage::Vector< linal::Vector3D> TransformCoordinates
      (
        const storage::Vector< linal::Vector3D> &COORDINATES,
        const math::TransformationMatrix3D &TRANSFORMATION_MATRIX
      );

    }; // class MergeFragmentComplete

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MERGE_FRAGMENT_COMPLETE_H_
