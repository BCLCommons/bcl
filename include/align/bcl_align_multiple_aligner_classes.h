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

#ifndef BCL_ALIGN_MULTIPLE_ALIGNER_CLASSES_H_
#define BCL_ALIGN_MULTIPLE_ALIGNER_CLASSES_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_aligner_dynamic_programming.h"
#include "bcl_align_aligner_merge.h"
#include "bcl_align_aligner_progressive.h"
#include "bcl_align_multiple_aligner_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MultipleAlignerClasses
    //! @brief enumerates all multiple aligners
    //!
    //! @see @link example_align_multiple_aligner_classes.cpp @endlink
    //! @author heinzes1
    //! @date Nov 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class MultipleAlignerClasses :
      public util::Enumerate< util::ShPtr< MultipleAlignerInterface< t_Member> >, MultipleAlignerClasses< t_Member> >
    {
      friend class util::Enumerate< util::ShPtr< MultipleAlignerInterface< t_Member> >, MultipleAlignerClasses< t_Member> >;

    public:

      //! @brief for MultipleAlignerClasses
      typedef util::Enum< util::ShPtr< MultipleAlignerInterface< t_Member> >, MultipleAlignerClasses< t_Member> > MultipleAlignerClass;

    //////////
    // data //
    //////////

      // declare all pairwise aligner classes
      const MultipleAlignerClass e_AlignerMerge;
      const MultipleAlignerClass e_AlignerProgressiveDynamicProgramming;
      const MultipleAlignerClass e_AlignerProgressiveDP;

      using util::Enumerate< util::ShPtr< MultipleAlignerInterface< t_Member> >, MultipleAlignerClasses< t_Member> >::AddEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      typedef util::ShPtr< MultipleAlignerInterface< t_Member> > SPAlignerType;

      //! @brief construct all PairwiseAlignerClasses
      MultipleAlignerClasses() :
        e_AlignerMerge( AddEnum( "Merge", SPAlignerType( new AlignerMerge< t_Member>()))),
        e_AlignerProgressiveDynamicProgramming
        (
          AddEnum
          (
            "ProgressiveDynamicProgramming",
            SPAlignerType
            (
              new AlignerProgressive< t_Member>
              (
                util::ShPtr< PairwiseAlignerInterface< t_Member> >( new AlignerDynamicProgramming< t_Member>())
              )
            )
          )
        ),
        e_AlignerProgressiveDP
        (
          AddEnum
          (
            "ProgressiveDP",
            SPAlignerType
            (
              new AlignerProgressive< t_Member>
              (
                util::ShPtr< PairwiseAlignerInterface< t_Member> >( new AlignerDP< t_Member>())
              )
            )
          )
        )
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

    }; // template class MultipleAlignerClasses

    //! @brief construct access function for all MultipleAlignerClasses
    //! @return reference to instances of MultipleAlignerClasses
    template< typename t_Member>
    const MultipleAlignerClasses< t_Member> &GetMultipleAlignerClasses()
    {
      return MultipleAlignerClasses< t_Member>::GetEnums();
    }

  } // namespace align

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< align::MultipleAlignerInterface< biol::AABase> >, align::MultipleAlignerClasses< biol::AABase> >;

  } // namespace util
} // namespace bcl

#endif // BCL_ALIGN_MULTIPLE_ALIGNER_CLASSES_H_
