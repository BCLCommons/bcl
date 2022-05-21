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

#ifndef BCL_CHEMISTRY_MUTATE_MULTI_FRAGMENT_H_
#define BCL_CHEMISTRY_MUTATE_MULTI_FRAGMENT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
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
    //! @class MutateMultiFragment
    //! @brief This class is for randomly using different rotamers of a fragment to change conformation of a molecule of interest.
    //! @details
    //!
    //! @see @link example_chemistry_mutate_multi_fragment.cpp @endlink
    //! @author mendenjl
    //! @date Jul 17, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateMultiFragment :
      public math::MutateInterface< FragmentComplete>
    {

    private:

    //////////
    // data //
    //////////

      //! rotamer data for fragment which this class mutates
      storage::Vector< util::ShPtrVector< MutateDihedralsInterface> > m_DihedralMutators;

      //! probability distribution for different rotamers of fragemnt
      storage::Vector< storage::Vector< double> >     m_ProbabilityDistributions;
      //! sums for each distribution
      storage::Vector< double>                        m_ProbabilityDistributionSums;
      storage::Vector< size_t>                        m_RingBondIndices;
      storage::Vector< size_t>                        m_ChainBondIndices;

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
      MutateMultiFragment();

      //! @brief constructor taking the member variable parameters
      //! @param FRAGMENT fragment which the mutate mutates
      //! @param GRAPH constitution graph of the molecule of interest
      //! @param MOLECULE molecule of interest
      MutateMultiFragment
      (
        const util::ShPtrVector< MutateDihedralsInterface> &FRAGMENT_MUTATES,
        const util::ShPtrVector< MutateDihedralsInterface> &DIHEDRAL_BOND_MUTATES_RANDOM,
        const double &RANDOM_DIHEDRAL_MUTATE_FREQUENCY,
        const size_t &N_DIHEDRALS
      );

      //! @brief Clone function
      //! @return pointer to new MutateMultiFragment
      MutateMultiFragment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetScheme() const;

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

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateMultiFragment

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_MULTI_FRAGMENT_H_
