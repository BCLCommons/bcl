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

#ifndef BCL_CHEMISTRY_MUTATE_BOND_ANGLES_H_
#define BCL_CHEMISTRY_MUTATE_BOND_ANGLES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_bond_angle_assignment.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_mutate_chirality.h"
#include "bcl_chemistry_mutate_dihedrals_interface.h"
#include "bcl_chemistry_ring_fragment_map.h"
#include "coord/bcl_coord_line_segment_3d.h"
#include "graph/bcl_graph_const_graph.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "random/bcl_random_histogram_1d_distribution.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateBondAngles
    //! @brief This class is for randomly using different rotamers of a fragment to change conformation of a molecule of interest.
    //! @details
    //!
    //! @see @link example_chemistry_mutate_bond_angles.cpp @endlink
    //! @author mendenjl
    //! @date Jun 19, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateBondAngles :
      public MutateDihedralsInterface
    {

    private:

    //////////
    // data //
    //////////

      //! rotamer data for fragment which this class mutates
      util::ShPtr< BondAngleAssignment>                                                     m_FragmentData;

      //! probability distribution for different rotamers of fragemnt
      util::ShPtr< random::Histogram1DDistribution>                                         m_ProbabilityDistribution;

      //! atoms to transform if a dihedral bond is mutated. The keys are the various non-ring bonds in fragment
      //! and value is a vector of atom vertices whose position is changed when the bond (key) is mutated
      storage::Vector< storage::Vector< size_t> >                                           m_ConnectedAtoms;

      bool                                                                                  m_ChangeChirality;

      //! Whether to perform unbiased sampling
      bool                                                                                  m_Unbiased;

      //! original chirality string
      ChiralityEnum                                                                         m_Chirality;

      // Mutate chirality. Only saved if change chirality or isometry is set and molecule is chiral
      util::SiPtr< const MutateChirality>                                                   m_MutateChirality;

      //! Number of ring bonds
      size_t                                                                                m_NumberRingBonds;

      //! dihedral angles for which this atom is a central atom
      storage::Vector< storage::VectorND< 4, size_t> > m_DihedralAngleIndices;

      //! whether to sample bond lengths
      bool m_SampleBondLengths;

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
      MutateBondAngles();

      //! @brief constructor taking the member variable parameters
      //! @param FRAGMENT fragment which the mutate mutates
      //! @param MOLECULE molecule of interest
      MutateBondAngles
      (
        const util::ShPtr< BondAngleAssignment> &FRAGMENT,
        const FragmentComplete &MOLECULE,
        const MutateChirality &MUTATE_CHIRALITY,
        bool CHANGE_CHIRALITY = false,
        bool UNBIASED = false,
        bool SAMPLE_LENGTHS = false
      );

      //! @brief Clone function
      //! @return pointer to new MutateBondAngles
      MutateBondAngles *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetScheme() const;

      //! @brief returns the fragment that the mutate uses
      //! @return the fragment that the mutate uses
      const util::ShPtr< BondAngleAssignment> &GetFragmentData() const;

      //! @brief returns probability distribution for rotamers of this particular fragments
      //! @return probability distribution for rotamers of this particular fragments
      const random::Histogram1DDistribution &GetProbability() const;

      const storage::Vector< size_t> &GetMoleculeDihedralBondIndices() const
      {
        static const storage::Vector< size_t> s_empty;
        return s_empty;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a conformation and returns a mutated conformation
      //! @param MOLECULE conformation of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< FragmentComplete> Mutate( const FragmentComplete &MOLECULE) const;

      //! @brief operator taking a conformation and returns a mutated conformation
      //! @param MOLECULE conformation of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< FragmentComplete> Mutate
      (
        const FragmentComplete &MOLECULE,
        bool SWAP_Z
      ) const;

      //! @brief rotate the molecule of interest at a particular bond at a particular angle
      //! @param ATOM_INFO the atom info that needs to be updated with coordinates after rotating about a particular bond
      //! @param CONNECTED_ATOMS the connected atoms associated with bond of interest
      //! @param NEW_ANGLE the new angle for the bond of interest
      //! @param EXISTING_ANGLE the angle for the bond of interest in the argument that is passed in
      static void RotateBond
      (
        AtomVector< AtomComplete> &ATOM_INFO,
        const storage::Vector< size_t> &CONNECTED_ATOMS,
        double ANGULAR_CHANGE,
        const size_t &ORIGINATING_ATOM
      );

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

    }; // class MutateBondAngles

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_BOND_ANGLES_H_
