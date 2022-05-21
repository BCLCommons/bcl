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

#ifndef BCL_ASSEMBLE_LOCATOR_DOMAIN_SPECIFIED_H_
#define BCL_ASSEMBLE_LOCATOR_DOMAIN_SPECIFIED_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorDomainSpecified
    //! @brief locates a domain created from a specific list of located sses
    //! @details domain is located through a list of sses that make up the domain
    //!
    //! @see @link example_assemble_locator_domain_specified.cpp @endlink
    //! @author alexanns
    //! @date Apr 20, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorDomainSpecified :
      public find::LocatorInterface< util::ShPtr< Domain>, ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! domain will be located by list of locators taking domains and returning siptr to sse
      util::ShPtrList< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > m_Locators;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorDomainSpecified();

      //! @brief constructor taking member variable
      //! @param LOCATORS locators to specify the sses that make up the domain
      LocatorDomainSpecified
      (
        const util::ShPtrList
        <
          find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
        > &LOCATORS
      );

      //! @brief Clone function
      //! @return pointer to new LocatorDomainSpecified
      LocatorDomainSpecified *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief writes a pymol formatted script to stream that highlights this domain
      //! @param OSTREAM the stream to which the script will be written to
      //! @param PYMOL_NAME the name of the selection should be in pymol
      //! @return ostream that the script was written to
      std::ostream &WritePymolDomainFile( std::ostream &OSTREAM, const std::string &PYMOL_NAME) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief locate function to locate domain from protein model
      //! @param MODEL the model from which the domain will be located
      //! @return shptr to the located domain
      util::ShPtr< Domain> Locate( const ProteinModel &MODEL) const;

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

    }; // class LocatorDomainSpecified

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_LOCATOR_DOMAIN_SPECIFIED_H_
