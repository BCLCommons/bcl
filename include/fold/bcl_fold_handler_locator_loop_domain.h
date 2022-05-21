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

#ifndef BCL_FOLD_HANDLER_LOCATOR_LOOP_DOMAIN_H_
#define BCL_FOLD_HANDLER_LOCATOR_LOOP_DOMAIN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerLocatorLoopDomain
    //! @brief is for reading in definitions of a loop domain locator from a file.
    //! @details  It also can create a loop domain locator for all the loops in a protein.
    //!
    //! @see @link example_fold_handler_locator_loop_domain.cpp @endlink
    //! @author alexanns
    //! @date Sep 7, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HandlerLocatorLoopDomain :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! bool indicating that the phi and psi values held in the pseudo residue should be randomized based on residue
      //! type and ramachandran distrubution even if they are defined by the starting protein model
      bool m_RandomizePseudoResiduePhiPsi;

      //! the string used to indicate the end of a loop domain
      static const std::string s_DomainEndString;

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
      //! @param RANDOMIZE_PHI_PSI bool true if the phi and psi values held in the pseudo residue should be randomized
      HandlerLocatorLoopDomain( const bool RANDOMIZE_PHI_PSI = true);

      //! @brief Clone function
      //! @return pointer to new HandlerLocatorLoopDomain
      HandlerLocatorLoopDomain *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief HandleRead creates a LocatorLoopDomain when from an istream and protein model
      //! @param ISTREAM the file stream which contains the loop domain locator information
      //! @param MODEL the protein model which will be used to help create the LocatorLoopDomain
      //! @return LocatorLoopDomain which has been created from ISTREAM and MODEL
      LocatorLoopDomain HandleRead( std::istream &ISTREAM, const assemble::ProteinModel &MODEL) const;

      //! @brief HandleReadMultiple can create a list of LocatorLoopDomain objects from a file and a protein model
      //! @param ISTREAM the file stream which contains the loop domain locators information
      //! @param MODEL the protein model which will be used to help create the LocatorLoopDomains
      //! @return list of LocatorLoopDomains which has been created from ISTREAM and MODEL
      util::ShPtrList< LocatorLoopDomain> HandleReadMultiple
      (
        std::istream &ISTREAM, const assemble::ProteinModel &MODEL
      ) const;

      //! @brief CreateLocatorLoopDomainsForInteriorCoil creates a list of loop domain locators from a protein model
      //! It creates a loop domain locator for every coil sse that is in the protein model
      //! @param MODEL the model for which a list of loop domain locators will be created for all its coil sses
      //! @return a list of LocatorLoopDomains created from MODEL - one for each coild sse
      const util::ShPtrList< LocatorLoopDomain> CreateLocatorLoopDomainsForInteriorCoil
      (
        const assemble::ProteinModel &MODEL
      ) const;

      //! @brief creates a list of loop domain locators from a protein model with bidirectionally grown loops
      //! every loop is assumed to be made up of two coils with on attached n and the other attached c terminally
      //! @param MODEL the model for which a list of loop domain locators will be created for all its coil sses
      //! @return a list of LocatorLoopDomains created from MODEL
      util::ShPtr< util::ShPtrList< LocatorLoopDomain> > CreateBidirectionalLocatorsForInteriorCoil
      (
        const assemble::ProteinModel &MODEL
      ) const;

      //! @brief creates a locator loop domain that will give an N to C direction loop domain
      //! @param SSE the sse that will be used to make a loop domain locator
      //! @param MODEL the model that will be used to make LocatorLoopDomain
      //! @return ShPtr to a LocatorLoopDomain
      util::ShPtr< LocatorLoopDomain>
      CreateNToCLocator( const assemble::SSE &SSE, const assemble::ProteinModel &MODEL) const;

      //! @brief creates a locator loop domain that will give an C to N direction loop domain
      //! @param SSE the sse that will be used to make a loop domain locator
      //! @param MODEL the model that will be used to make LocatorLoopDomain
      //! @return ShPtr to a LocatorLoopDomain
      util::ShPtr< LocatorLoopDomain>
      CreateCToNLocator( const assemble::SSE &SSE, const assemble::ProteinModel &MODEL) const;

      //! @brief GetFormat gives the format that this handler needs in order to work
      //! @return string which describes the format needed by this handler in order for it to work
      std::string GetFormat() const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief CreateLocatorLoopSegments creates a list of loop segment locators from an istream
      //! @param ISTREAM the file stream which contains the information to create a list of loop segment locators
      //! @return a list of LocatorLoopSegments created from ISTREAM
      const storage::List< LocatorLoopSegment> CreateLocatorLoopSegments( std::istream &ISTREAM) const;

      //! @brief CreateNTerminalSSELocator creates the locator to the nterminal anchor sse of a loop domain
      //! @param LOOP_SEGMENT_LOCATORS list of loop segment locators which make up a loop domain
      //! @param MODEL the protein model which the loop segment locators and the loop domain refer
      //! @return LocatorSSE which can locate the nterminal anchor sse of the loop domain
      const assemble::LocatorSSE CreateNTerminalSSELocator
      (
        const storage::List< LocatorLoopSegment> &LOOP_SEGMENT_LOCATORS, const assemble::ProteinModel &MODEL
      ) const;

    }; // class HandlerLocatorLoopDomain

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_HANDLER_LOCATOR_LOOP_DOMAIN_H_ 
