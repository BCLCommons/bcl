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

#ifndef BCL_ASSEMBLE_LOCATOR_CHAIN_H_
#define BCL_ASSEMBLE_LOCATOR_CHAIN_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorChain
    //! @brief This class locates a specified chain in the given protein model
    //! @details For a given protein model, this class iterates over its chains and returns the one with the specified
    //! chain id.
    //!
    //! @see @link example_assemble_locator_chain.cpp @endlink
    //! @author alexanns
    //! @date 01/16/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorChain :
      public find::LocatorInterface< util::SiPtr< const Chain>, ProteinModel>
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      char m_Chain_ID; //!< chain identifier

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorChain() :
        m_Chain_ID( 'A')
      {
      }

      //! @brief constructor from chain id
      //! @param CHAIN_ID char which indicates the chain
      explicit LocatorChain( const char CHAIN_ID) :
        m_Chain_ID( CHAIN_ID)
      {
      }

      //! @brief virtual copy constructor
      LocatorChain *Clone() const
      {
        return new LocatorChain( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetChainID gives the LocatorChain
      //! @return returns "m_Chain_ID" character
      const char &GetChainID() const;

      //! @brief SetChainID changes the chain id
      //! @param CHAINID char which indicates the new chain id
      void SetChainID( const char CHAINID);

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate gets the chain of the ProteinModel symbolized by the LocatorChain
      //! @param MODEL ProteinModel for which the chain is wanted
      //! @return returns
      util::SiPtr< const Chain> Locate( const ProteinModel &MODEL) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;
    }; // class LocatorChain

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_CHAIN_H_
