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

#ifndef BCL_ALIGN_WORD_GENERATOR_INTERFACE_H_
#define BCL_ALIGN_WORD_GENERATOR_INTERFACE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_word.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WordGeneratorInterface
    //! @brief generates a vector of words from an alignment
    //!
    //! @tparam t_Member the type of object that is used by the assignments in the alignment
    //!
    //! @remarks example unnecessary
    //! @author heinzes1
    //! @date Mar 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class WordGeneratorInterface :
      public util::ObjectInterface
    {

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief generates a vector of words from an alignment
      //! @param QUERY_ALIGNMENT the alignment
      //! @param WORD_LENGTH the length of the words to be generated
      //! @return a vector of words
      virtual storage::Vector< AlignmentWord< t_Member> > Generate
      (
        const util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
        const size_t WORD_LENGTH
      ) const = 0;

    }; // template class WordGeneratorInterface

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_WORD_GENERATOR_INTERFACE_H_ 
