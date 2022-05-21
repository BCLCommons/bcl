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

#ifndef BCL_CHEMISTRY_MUTATE_DIHEDRALS_INTERFACE_H_
#define BCL_CHEMISTRY_MUTATE_DIHEDRALS_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_clash_score.h"
#include "bcl_chemistry_fragment_complete.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDihedralsInterface
    //! @brief Interace for things that mutate a molecules dihedral bonds
    //!
    //! @see @link example_chemistry_mutate_dihedrals_interface.cpp @endlink
    //! @author mendenjl
    //! @date Jul 17, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDihedralsInterface :
      public math::MutateInterface< FragmentComplete>
    {
    protected:

      //! Running average of average mutate weight
      mutable math::RunningAverage< double> m_MutateWeight;

      //! Shared point to bool telling the mutate whether it should only mutate to the best (score-wise) term
      util::ShPtr< bool> m_OnlyMutateToBest;

    public:

      //! @brief default constructor
      MutateDihedralsInterface() :
        m_OnlyMutateToBest( new bool( false))
      {
      }

      //! @brief Clone function
      //! @return pointer to new MutateDihedralsInterface
      virtual MutateDihedralsInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the dihedral bond indices referenced by this mutate
      virtual const storage::Vector< size_t> &GetMoleculeDihedralBondIndices() const = 0;

      //! @brief Virtual function to do the actual mutation. Should not check clashes
      virtual math::MutateResult< FragmentComplete> Mutate( const FragmentComplete &FRAGMENT) const = 0;

      //! @brief implementation for testing clash frequency
      virtual math::MutateResult< FragmentComplete> operator()( const FragmentComplete &FRAGMENT) const
      {
        return Mutate( FRAGMENT);
      }
      
      //! @brief get the choose best only flag
      const bool &GetChooseBestOnlyFlag() const
      {
        return *m_OnlyMutateToBest;
      }

      //! @brief get changeable bool that allows us to switch between only selecting the approximate best choices
      util::ShPtr< bool> &GetChooseBestOnlySharedPointer()
      {
        return m_OnlyMutateToBest;
      }

    }; // class MutateDihedralsInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_DIHEDRALS_INTERFACE_H_
