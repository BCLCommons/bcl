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

#ifndef BCL_ALIGN_SEQUENCE_INTERFACE_H_
#define BCL_ALIGN_SEQUENCE_INTERFACE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SequenceInterface
    //! @brief class is an interface to a sequence
    //! @details This interface defines all sequence specific functions.
    //!
    //! @tparam t_Member type of the elements in the sequence
    //!
    //! @remarks example unnecessary
    //! @author heinzes1
    //! @date May 31, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class SequenceInterface :
      public virtual util::ObjectInterface
    {
    public:

    //////////////
    // typedefs //
    //////////////

      //! @brief typedef to iterate over non-const members of a sequence
      typedef typename util::ShPtrVector< t_Member>::iterator       iterator;
      //! @brief typedef to iterate over const members of a sequence
      typedef typename util::ShPtrVector< t_Member>::const_iterator const_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct a new SequenceInterface from a given interface implementation, sequence id and members
      //! @param ID sequence identifier
      //! @param MEMBERS sequence members
      //! @param CHAIN_ID chain id used to build the sequence
      //! @return pointer to a SequenceInterface
      virtual SequenceInterface *Construct
      (
        const std::string &ID,
        const std::string &MEMBERS,
        const char CHAIN_ID = 'A'
      ) const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the number of member in the sequence
      //! @return the number of members contained in the sequence
      virtual size_t GetSize() const = 0;

      //! @brief returns a vector of members
      //! @return a vector of SiPtrs of type t_Member contained in the sequence
      virtual util::SiPtrVector< const t_Member> GetMembers() const = 0;

      //! @brief returns a SiPtr to the first element
      //! @return SiPtr to the first t_Member
      virtual util::SiPtr< const t_Member> GetFirstMember() const = 0;

      //! @brief returns a SiPtr to the last element
      //! @return SiPtr to the last t_Member
      virtual util::SiPtr< const t_Member> GetLastMember() const = 0;

      //! @brief add a t_Member to the sequence with type based on ONE_LETTER_CODE
      //! @param ONE_LETTER_CODE char encoding the type of the t_Member
      virtual void AddMember( const char &ONE_LETTER_CODE) = 0;

      //! @brief return sequence identifier
      //! @return sequence identifier of the sequence
      virtual std::string GetSequenceId() const = 0;

      //! @brief return iterator to Begin of Sequence
      //! @return iterator to Begin of Sequence
      virtual typename util::ShPtrVector< t_Member>::iterator Begin() = 0;

      //! @brief return const_iterator to Begin of Sequence
      //! @return const_iterator to Begin of Sequence
      virtual typename util::ShPtrVector< t_Member>::const_iterator Begin() const = 0;

      //! @brief return iterator to End of Sequence
      //! @return iterator to End of Sequence
      virtual typename util::ShPtrVector< t_Member>::iterator End() = 0;

      //! @brief return const_iterator to End of Sequence
      //! @return const_iterator to End of Sequence
      virtual typename util::ShPtrVector< t_Member>::const_iterator End() const = 0;

      //! @brief return iterator to reverse begin of Sequence
      //! @return iterator to reverse begin of Sequence
      virtual typename util::ShPtrVector< t_Member>::reverse_iterator ReverseBegin() = 0;

      //! @brief return const_iterator to reverse begin of Sequence
      //! @return const_iterator to reverse begin of Sequence
      virtual typename util::ShPtrVector< t_Member>::const_reverse_iterator ReverseBegin() const = 0;

      //! @brief return iterator to reverse end of Sequence
      //! @return iterator to reverse end of Sequence
      virtual typename util::ShPtrVector< t_Member>::reverse_iterator ReverseEnd() = 0;

      //! @brief return const_iterator to reverse end of Sequence
      //! @return const_iterator to reverse end of Sequence
      virtual typename util::ShPtrVector< t_Member>::const_reverse_iterator ReverseEnd() const = 0;

    }; // template class SequenceInterface

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_SEQUENCE_INTERFACE_H_ 
