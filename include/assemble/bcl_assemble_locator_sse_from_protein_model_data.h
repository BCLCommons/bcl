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

#ifndef BCL_ASSEMBLE_LOCATOR_SSE_FROM_PROTEIN_MODEL_DATA_H_
#define BCL_ASSEMBLE_LOCATOR_SSE_FROM_PROTEIN_MODEL_DATA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_protein_model_data.h"
#include "find/bcl_find_locator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSSEFromProteinModelData
    //! @brief retrieves locators from protein model data using given key and uses one randomly to locate an SSE
    //! @details retrieves locators from protein model data using given key and uses one randomly to locate an SSE
    //!
    //! @see @link example_assemble_locator_sse_from_protein_model_data.cpp @endlink
    //! @author karakam, alexanns
    //! @date Jun 16, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorSSEFromProteinModelData :
      public find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! type of data
      ProteinModelData::TypeEnum m_Key;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorSSEFromProteinModelData();

      //! @brief constructor from a key
      //! @param KEY key to retrieve the locators from the protein model data
      LocatorSSEFromProteinModelData( const ProteinModelData::Type KEY);

      //! @brief Clone function
      //! @return pointer to new LocatorSSEFromProteinModelData
      LocatorSSEFromProteinModelData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief locate and return an SSE using the locators retrieved from ProteinModelData
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return returns SiPtr to selected SSE
      util::SiPtr< const SSE> Locate( const DomainInterface &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LocatorSSEFromProteinModelData

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_LOCATOR_SSE_FROM_PROTEIN_MODEL_DATA_H_
