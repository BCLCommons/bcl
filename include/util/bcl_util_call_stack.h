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

#ifndef BCL_UTIL_CALL_STACK_H_
#define BCL_UTIL_CALL_STACK_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <vector>

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CallStack
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_util_call_stack.cpp @endlink
    //! @author woetzen
    //! @date Aug 28, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CallStack
    {

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Entry
      //! @brief TODO: add a brief comment
      //! @details TODO: add an detailed description to this class
      //!
      //! @remarks example unnecessary
      //! @author woetzen
      //! @date Aug 28, 2011
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class BCL_API Entry
      {

      private:

      //////////
      // data //
      //////////

        std::string m_File;       //! file name
        size_t      m_LineNumber; //! line number in file
        std::string m_Function;   //! function signature

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief constructor
        Entry( const std::string &FILE_NAME, const size_t LINE_NUMBER, const std::string &FUNCTION);

      ////////////////
      // operations //
      ////////////////

        //! @brief return entry as string
        std::string String() const;

      }; // class Entry

    private:

    //////////
    // data //
    //////////

      std::vector< Entry> m_Stack; //!< complete stack

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param DISCARD number of stack entries to discard at the top
      CallStack( const size_t DISCARD = 0);

    ////////////////
    // operations //
    ////////////////

      //! @brief convert the stacktrace to a string
      std::string String() const;

    }; // class CallStack

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_CALL_STACK_H_ 
