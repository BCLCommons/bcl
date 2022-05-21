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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_locator_chain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LocatorChain::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorChain())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorChain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! GetChainID gives the LocatorChain
    //! @return returns "m_Chain_ID" character
    const char &LocatorChain::GetChainID() const
    {
      return m_Chain_ID;
    }

    //! SetChainID changes the chain id
    //! @param CHAINID char which indicates the new chain id
    void LocatorChain::SetChainID( const char CHAINID)
    {
      // set "m_Chain_ID" to "CHAINID"
      m_Chain_ID = CHAINID;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &LocatorChain::GetAlias() const
    {
      static const std::string s_Name( "LocatorChain");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate gets the chain of the ProteinModel symbolized by the LocatorChain
    //! @param MODEL ProteinModel for which the chain is wanted
    //! @return returns
    util::SiPtr< const Chain> LocatorChain::Locate( const ProteinModel &MODEL) const
    {
      util::SiPtr< const Chain> chain( MODEL.GetChain( m_Chain_ID));
      if( !chain.IsDefined())
      {
        BCL_MessageDbg
        (
          "chain " + util::Format()( m_Chain_ID) + " does not exist in protein model"
        );
      }
      return chain;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorChain::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Locates a specified chain in the given protein model");

      parameters.AddInitializer
      (
        "chain_id",
        "the chain id of the desired chain",
        io::Serialization::GetAgent( &m_Chain_ID),
        "A"
      );

      return parameters;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorChain::Read( std::istream &ISTREAM)
    {
      return util::SerializableInterface::Read( ISTREAM);
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorChain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return util::SerializableInterface::Write( OSTREAM, INDENT);
    }

  } // namespace assemble
} // namespace bcl
