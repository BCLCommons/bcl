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

#ifndef BCL_ASSEMBLE_PICK_SSES_RANDOM_H_
#define BCL_ASSEMBLE_PICK_SSES_RANDOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "find/bcl_find_pick_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickSSEsRandom
    //! @brief randomly picks a given number of SSEs of a given type from a SSE list
    //!
    //! @note if no SSE types and number of SSEs to pick to be considered are provided through the constructor, strands
    //! and helices will be considered as valid and one SSE will be picked
    //!
    //! @see @link example_assemble_pick_sses_random.cpp @endlink
    //! @author fischea
    //! @date April 20, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickSSEsRandom :
      virtual public find::PickInterface< util::SiPtrList< const SSE>, util::SiPtrList< const SSE> >
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! types of SSEs which shall be considered for the picking
      storage::Set< biol::SSType> m_SSTypes;

      //! number of SSEs to pick
      size_t m_NumSSEsToPick;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      //! @detail sets the SSE types to be considered to helix and strand and the number of SSE to pick to one
      PickSSEsRandom();

      //! @brief construct from the SSE types to be considered and number of SSEs to pick
      //! @param SS_TYPES types of SSEs which shall be considered for the picking
      //! @param NUM_SSES number of SSEs to pick randomly
      PickSSEsRandom( const storage::Set< biol::SSType> &SS_TYPES, const size_t NUM_SSES);

      //! @brief returns a pointer to a new PickSSEsRandom
      //! @return pointer to a new PickSSEsRandom
      PickSSEsRandom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief randomly picks the given number of SSEs of the given type from the list
      //! @param SSE_LIST SSE list to pick the SSEs from
      //! @return list with the picked SSEs
      util::SiPtrList< const SSE> Pick( const util::SiPtrList< const SSE> &SSE_LIST) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read object from input stream
      //! @param ISTREAM input stream to read object from
      //! @return input stream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write object into  output stream
      //! @param OSTREAM output stream to write object into
      //! @param INDENT number of indentations to separate members
      //! @return output stream object was written into
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PickSSEsRandom

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PICK_SSES_RANDOM_H_
