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

#ifndef BCL_CHEMISTRY_MOLECULAR_CONFORMATION_SHARED_H_
#define BCL_CHEMISTRY_MOLECULAR_CONFORMATION_SHARED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_fragment_conformation_shared.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MolecularConformationShared
    //! @brief Class that contains molecular configuration data
    //! @details Models stereochemistry, isomeric fragments, chemical adjacency, and aromatic and ring structures
    //!
    //! @see @link example_chemistry_molecular_conformation_shared.cpp @endlink
    //! @author kothiwsk
    //! @date Dec 22, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MolecularConformationShared :
      public FragmentConformationShared
    {

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
      MolecularConformationShared();

      //! @brief constructor given atoms with conformation and molecule configuration
      //! @params CONFIGURATION molecule configuration
      //! @params ATOMCONFORMATION atom conformation objects
      MolecularConformationShared
      (
        const util::ShPtr< MolecularConfigurationShared> &CONFIGURATION,
        const AtomVector< AtomConformationalShared> &ATOMCONFIGURATION
      );

      //! @brief constructor with a conformation interface
      //! @params CONFORMATION conformation interface
      MolecularConformationShared( const ConformationInterface &CONFORMATION);

      //! @brief Clone function
      //! @return pointer to new MolecularConformationShared
      MolecularConformationShared *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

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
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class MolecularConformationShared

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULAR_CONFORMATION_SHARED_H_ 
