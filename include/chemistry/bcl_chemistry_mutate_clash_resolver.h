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

#ifndef BCL_CHEMISTRY_MUTATE_CLASH_RESOLVER_H_
#define BCL_CHEMISTRY_MUTATE_CLASH_RESOLVER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_mutate_bond_angles.h"
#include "bcl_chemistry_mutate_dihedral_bond.h"
#include "bcl_chemistry_mutate_dihedrals_interface.h"
#include "bcl_chemistry_rotamer_dihedral_bond_data.h"
#include "coord/bcl_coord_line_segment_3d.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "random/bcl_random_histogram_1d_distribution.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateClashResolver
    //! @brief This class is for randomly using different rotamers of a fragment to change conformation of a molecule of interest.
    //! @details
    //!
    //! @see @link example_chemistry_mutate_clash_resolver.cpp @endlink
    //! @author mendenjl
    //! @date Nov 18, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateClashResolver :
      public math::MutateInterface< FragmentComplete>
    {

    private:

    //////////
    // data //
    //////////

      //! Parameters

      //! maximum cycles; times that clashes are rediscovered
      double m_MaxCycles;

      //! Molecule specific data

      //! rotamer data for fragment which this class mutates
      //! Outer key: atom index, inner key, none, ordered.
      storage::Vector< util::ShPtrVector< MutateDihedralBond> > m_DihedralMutators;

      //! Bond angle mutators. Undefined if bond angle modifications are not allowed
      util::ShPtrVector< MutateBondAngles> m_BondAngleMutators;

      //! Graph of the molecule
      graph::ConstGraph< size_t, size_t> m_Graph;

      //! shared clash score
      util::ShPtr< AtomClashScore> m_ClashScore;

      //! number of rotatable dihedrals
      size_t m_NumberDihedrals;

      //! mobile atoms during clash resolution
      storage::Vector< size_t> m_MobileAtoms;

      //! Cached data

      //! cached ways of handling clashes between two clashing atoms
      mutable storage::Map
      <
        std::pair< size_t, size_t>,
        storage::Pair< util::ShPtrVector< MutateDihedralBond>, util::ShPtrVector< MutateBondAngles> >
      > m_Map;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking the member variable parameters
      MutateClashResolver
      (
        const double &MAX_CYCLES = 1.0,
        const storage::Vector< size_t> &MOBILE_ATOMS = storage::Vector< size_t>()
      ) :
        m_MaxCycles( MAX_CYCLES),
        m_NumberDihedrals( 0),
        m_MobileAtoms( MOBILE_ATOMS)
      {
      }

      //! @brief Clone function
      //! @return pointer to new MutateClashResolver
      MutateClashResolver *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetScheme() const;

      //! @brief setup this mutate to handle a new molecule
      void Setup
      (
        const FragmentComplete &CONF,
        const util::ShPtrVector< MutateBondAngles> &MUT_BA,
        const util::ShPtr< AtomClashScore> &CLASH
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a conformation and returns a mutated conformation
      //! @param MOLECULE conformation of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< FragmentComplete> operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief cache the resolution path for a particular clash
      const storage::Pair< util::ShPtrVector< MutateDihedralBond>, util::ShPtrVector< MutateBondAngles> > &GetResolutionPaths
      (
        const size_t &INDEX_A,
        const size_t &INDEX_B
      ) const;

    }; // class MutateClashResolver

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_CLASH_RESOLVER_H_
