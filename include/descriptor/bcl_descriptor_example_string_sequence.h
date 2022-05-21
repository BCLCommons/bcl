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

#ifndef BCL_DESCRIPTOR_EXAMPLE_STRING_SEQUENCE_H_
#define BCL_DESCRIPTOR_EXAMPLE_STRING_SEQUENCE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_sequence_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StringSequence
    //! @brief simple class that derives from descriptor::SequenceInterface and can be used to test various descriptors
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Dec 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API StringSequence :
      public SequenceInterface< char>,
      public util::SerializableInterface
    {
    private:

      //! Internally held string
      std::string m_String;

      // single instance of that class
      static util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

      //! @brief constructor from a string
      StringSequence( const std::string &STRING = "");

      //! @brief virtual copy constructor
      StringSequence *Clone() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the alias for this class over the command line
      const std::string &GetAlias() const;

      //! @brief get the iterator for the sequence
      //! @return the iterator for the sequence
      iterate::Generic< const char> GetIterator() const;

      //! @brief get the internally-held string
      const std::string &GetString() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get the iteration result as a string (semi-colons between entries)
      //! @param DESCRIPTOR the descriptor to use
      //! @param STRING the string of interest
      //! @param PRECISION the fixed float point precision # of digits to use
      //! @return the resulting string
      static std::string WriteIterations
      (
        const util::Implementation< Base< char, float> > &DESCRIPTOR,
        const std::string &STRING,
        const size_t &PRECISION
      );

      //! @brief get the iteration result as a string (semi-colons between entries)
      //! @param DESCRIPTOR the descriptor to use
      //! @param STRING the string of interest
      //! @return the resulting string
      static std::string WriteIterations
      (
        const util::Implementation< Base< char, char> > &DESCRIPTOR,
        const std::string &STRING
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    };

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_EXAMPLE_STRING_SEQUENCE_H_

