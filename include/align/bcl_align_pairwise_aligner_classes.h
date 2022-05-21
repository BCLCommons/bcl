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

#ifndef BCL_ALIGN_PAIRWISE_ALIGNER_CLASSES_H_
#define BCL_ALIGN_PAIRWISE_ALIGNER_CLASSES_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_aligner_dp.h"
#include "bcl_align_aligner_dynamic_programming.h"
#include "bcl_align_aligner_merge.h"
#include "bcl_align_aligner_wordbased.h"
#include "bcl_align_pairwise_aligner_interface.h"
#include "biol/bcl_biol_aa_base.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PairwiseAlignerClasses
    //! @brief enumerates all pairwise aligners
    //!
    //! @see @link example_align_pairwise_aligner_classes.cpp @endlink
    //! @author heinzes1
    //! @date Nov 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class PairwiseAlignerClasses :
      public util::Enumerate< util::ShPtr< PairwiseAlignerInterface< t_Member> >, PairwiseAlignerClasses< t_Member> >
    {
      friend class util::Enumerate< util::ShPtr< PairwiseAlignerInterface< t_Member> >, PairwiseAlignerClasses< t_Member> >;

    public:

      //! @brief for PairwiseAlignerClass
      typedef util::Enum< util::ShPtr< PairwiseAlignerInterface< t_Member> >, PairwiseAlignerClasses< t_Member> > PairwiseAlignerClass;

    //////////
    // data //
    //////////

      // declare all pairwise aligner classes
      const PairwiseAlignerClass e_AlignerDynamicProgramming;
      const PairwiseAlignerClass e_AlignerDP;
      const PairwiseAlignerClass e_AlignerMerge;
      const PairwiseAlignerClass e_AlignerWordbased;

      using util::Enumerate< util::ShPtr< PairwiseAlignerInterface< t_Member> >, PairwiseAlignerClasses< t_Member> >::AddEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      typedef util::ShPtr< PairwiseAlignerInterface< t_Member> > SPAlignerType;

      //! @brief construct all PairwiseAlignerClasses
      PairwiseAlignerClasses() :
        e_AlignerDynamicProgramming( AddEnum( "DynamicProgramming", SPAlignerType( new AlignerDynamicProgramming< t_Member>()))),
        e_AlignerDP(                 AddEnum( "DP"                , SPAlignerType( new AlignerDP< t_Member>()))),
        e_AlignerMerge(              AddEnum( "Merge"             , SPAlignerType( new AlignerMerge< t_Member>()))),
        e_AlignerWordbased(          AddEnum( "Wordbased"         , SPAlignerType( new AlignerWordbased< t_Member>())))
      {
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    }; // template class PairwiseAlignerClasses

    //! @brief construct access function for all PairwiseAlignerClasses
    //! @return reference to instances of PairwiseAlignerClasses
    template< typename t_Member>
    const PairwiseAlignerClasses< t_Member> &GetPairwiseAlignerClasses()
    {
      return PairwiseAlignerClasses< t_Member>::GetEnums();
    }

  } // namespace align

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< align::PairwiseAlignerInterface< biol::AABase> >, align::PairwiseAlignerClasses< biol::AABase> >;

  } // namespace util
} // namespace bcl

#endif // BCL_ALIGN_PAIRWISE_ALIGNER_CLASSES_H_
