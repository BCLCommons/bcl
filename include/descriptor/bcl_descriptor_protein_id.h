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

#ifndef BCL_DESCRIPTOR_PROTEIN_ID_H_
#define BCL_DESCRIPTOR_PROTEIN_ID_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinId
    //! @brief Returns the id (up to 5 character) stored on the protein model (usually the PDB id, optionally the chain)
    //!
    //! @see @link example_descriptor_protein_id.cpp @endlink
    //! @author mendenjl
    //! @date Feb 15, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ProteinId :
      public BaseSequence< biol::AABase, char>
    {

    //////////
    // data //
    //////////

      std::string m_Id; //!< ID of the current protein, cached to avoid repeated dynamic casts
      size_t m_NumChar;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ProteinId( const size_t &N_CHAR = 5) :
        m_NumChar( N_CHAR)
      {
      }

      //! @brief Clone function
      //! @return pointer to new BaseElement
      ProteinId *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    /////////////////
    // data access //
    /////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< char> &STORAGE);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

    }; // class AASeqID

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_PROTEIN_ID_H_
