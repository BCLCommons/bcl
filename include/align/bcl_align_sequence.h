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

#ifndef BCL_ALIGN_SEQUENCE_H_
#define BCL_ALIGN_SEQUENCE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_sequence_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Sequence
    //! @brief class is a reference implementation of SequenceInterface
    //! @details This interface defines all sequence specific functions.
    //! The data does never make explicit hardcopies of the members
    //!
    //! @tparam t_Member type of the elements in the sequence
    //!
    //! @see @link example_align_sequence.cpp @endlink
    //! @author woetzen
    //! @date 22/06/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class Sequence :
      public SequenceInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      //! vector of member
      util::ShPtrVector< t_Member> m_Data;

      //! identification
      std::string m_Identification;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from ShPtrVector< t_Member>
      //! @param MEMBER the member
      //! @param IDENTIFICATION identification for that sequence
      Sequence( const util::ShPtrVector< t_Member> &MEMBER, const std::string &IDENTIFICATION) :
        m_Data( MEMBER),
        m_Identification( IDENTIFICATION)
      {
      }

      //! @brief construct from sequence interface and identification
      //! @param SEQUENCE the sequence to copy
      //! @param IDENTIFICATION identification for that sequence
      Sequence( const SequenceInterface< t_Member> &SEQEUENCE, const std::string &IDENTIFICATION) :
        m_Data( SEQEUENCE.Begin(), SEQEUENCE.End()),
        m_Identification( IDENTIFICATION)
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to new Sequence
      Sequence< t_Member> *Clone() const
      {
        return new Sequence< t_Member>( *this);
      }

      //! @brief construct a new SequenceInterface from a given interface implementation, sequence id and members
      //! @param ID sequence identifier
      //! @param MEMBERS sequence members
      //! @param CHAIN_ID chain id used to build the sequence
      //! @return pointer to a SequenceInterface
      Sequence< t_Member> *Construct
      (
        const std::string &ID,
        const std::string &MEMBERS,
        const char CHAIN_ID = 'A'
      ) const
      {
        Sequence< t_Member> *seq( new Sequence< t_Member>( util::ShPtrVector< t_Member>(), ID));
        for( std::string::const_iterator itr( MEMBERS.begin()), itr_end( MEMBERS.end()); itr != itr_end; ++itr)
        {
          seq->AddMember( *itr);
        }

        return seq;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the number of member in the sequence
      //! @return the number of members contained in the sequence
      size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      //! @brief returns a vector of members
      //! @return a vector of SiPtrs of type t_Member contained in the sequence
      util::SiPtrVector< const t_Member> GetMembers() const
      {
        return util::SiPtrVector< const t_Member>( m_Data);
      }

      //! @brief returns a SiPtr to the first element
      //! @return SiPtr to the first t_Member
      util::SiPtr< const t_Member> GetFirstMember() const
      {
        return m_Data.FirstElement();
      }

      //! @brief returns a SiPtr to the last element
      //! @return SiPtr to the last t_Member
      util::SiPtr< const t_Member> GetLastMember() const
      {
        return m_Data.LastElement();
      }

      //! @brief add a t_Member to the sequence with type based on ONE_LETTER_CODE
      //! @param ONE_LETTER_CODE char encoding the type of the t_Member
      void AddMember( const char &ONE_LETTER_CODE)
      {
        m_Data.PushBack( util::ShPtr< t_Member>( t_Member::Construct( ONE_LETTER_CODE, m_Data.GetSize() + 1)));
      }

      //! @brief return sequence identifier
      //! @return sequence identifier of the sequence
      std::string GetSequenceId() const
      {
        return m_Identification;
      }

      //! @brief return iterator to Begin of Sequence
      //! @return iterator to Begin of Sequence
      typename util::ShPtrVector< t_Member>::iterator Begin()
      {
        return m_Data.Begin();
      }

      //! @brief return const_iterator to Begin of Sequence
      //! @return const_iterator to Begin of Sequence
      typename util::ShPtrVector< t_Member>::const_iterator Begin() const
      {
        return m_Data.Begin();
      }

      //! @brief return iterator to End of Sequence
      //! @return iterator to End of Sequence
      typename util::ShPtrVector< t_Member>::iterator End()
      {
        return m_Data.End();
      }

      //! @brief return const_iterator to End of Sequence
      //! @return const_iterator to End of Sequence
      typename util::ShPtrVector< t_Member>::const_iterator End() const
      {
        return m_Data.End();
      }

      //! @brief return iterator to reverse begin of Sequence
      //! @return iterator to reverse begin of Sequence
      typename util::ShPtrVector< t_Member>::reverse_iterator ReverseBegin()
      {
        return m_Data.ReverseBegin();
      }

      //! @brief return const_iterator to reverse begin of Sequence
      //! @return const_iterator to reverse begin of Sequence
      typename util::ShPtrVector< t_Member>::const_reverse_iterator ReverseBegin() const
      {
        return m_Data.ReverseBegin();
      }

      //! @brief return iterator to reverse end of Sequence
      //! @return iterator to reverse end of Sequence
      typename util::ShPtrVector< t_Member>::reverse_iterator ReverseEnd()
      {
        return m_Data.ReverseEnd();
      }

      //! @brief return const_iterator to reverse end of Sequence
      //! @return const_iterator to reverse end of Sequence
      typename util::ShPtrVector< t_Member>::const_reverse_iterator ReverseEnd() const
      {
        return m_Data.ReverseEnd();
      }

      //! @brief append a second sequence
      //! @param SEQUENCE the sequence to append
      void Append( const SequenceInterface< t_Member> &SEQUENCE)
      {
        m_Data.Append( util::ShPtrVector< t_Member>( SEQUENCE.Begin(), SEQUENCE.End()));
      }

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_Data          , ISTREAM);
        io::Serialize::Read( m_Identification, ISTREAM);

        //return
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT indentation
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Data          , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Identification, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

    }; // template class Sequence

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_SEQUENCE_H_ 
