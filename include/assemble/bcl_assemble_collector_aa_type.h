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

#ifndef BCL_ASSEMBLE_COLLECTOR_AA_TYPE_H_
#define BCL_ASSEMBLE_COLLECTOR_AA_TYPE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "find/bcl_find_locator_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorAAType
    //! @brief collects residues which match any of a set of desired aa types
    //! @details if a residue has a type that matches one of the desired aa types, then it is added to the collected
    //!          residues
    //!
    //! @see @link example_assemble_collector_aa_type.cpp @endlink
    //! @author alexanns
    //! @date Mar 3, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorAAType :
      public find::CollectorInterface< util::SiPtrList< const biol::AABase>, util::SiPtrVector< const biol::AABase> >
    {

    private:

    //////////
    // data //
    //////////

      //! set of aa types that will be collected from the list of residues
      storage::Set< biol::AAType> m_AATypes;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CollectorAAType();

      //! @brief constructor taking member variables as paramters
      //! @param AA_TYPES set of aa types that will be collected from the list of residues
      CollectorAAType( const storage::Set< biol::AAType> &AA_TYPES);

      //! @brief Clone function
      //! @return pointer to new CollectorAAType
      CollectorAAType *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! Collect the residues that match the desired aa types from the provided residues
      //! @param RESIDUES entity that contains the residues to be collected from
      //! @return returns list of the collected residues that match the desired aa types
      util::SiPtrList< const biol::AABase> Collect( const util::SiPtrVector< const biol::AABase> &RESIDUES) const;

      //! @brief locate the residue that matches the desired aa type - if more than one just returns first encountered
      //! @param RESIDUES entity that contains the residues to be collected from
      //! @return const siptr to biol::AABase which is the first instance of the aa type found in RESIDUES
      util::SiPtr< const biol::AABase> Locate( const util::SiPtrVector< const biol::AABase> &RESIDUES) const;

    ///////////////
    // operators //
    ///////////////

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class CollectorAAType

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_AA_TYPE_H_
