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

#ifndef BCL_CHEMISTRY_MOLECULE_COMPLETE_H_
#define BCL_CHEMISTRY_MOLECULE_COMPLETE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeComplete
    //! @brief Class that contains molecular configuration data
    //! @details Models stereochemistry, isomeric fragments, chemical adjacency, and aromatic and ring structures
    //!
    //! @see @link example_chemistry_molecule_complete.cpp @endlink
    //! @author mendenjl
    //! @date Mar 08, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeComplete :
      public FragmentComplete
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
      MoleculeComplete();

      //! @brief pass-through constructor
      //! @param CONSTRUCTOR_OBJECT object that can be used to construct the base class
      template< typename t_DataType>
      MoleculeComplete( const t_DataType &CONSTRUCTOR_OBJECT) :
        FragmentComplete( CONSTRUCTOR_OBJECT)
      {
        // note that we must call the conformation interfaces' implementation of GetNumberValences, since this class
        // overrides GetNumberValences, always returning 0; but here we need to check that this is true!
        if( ConformationInterface::GetNumberValences() != size_t( 0))
        {
          SaturateWithH();
          BCL_Assert( ConformationInterface::GetNumberValences() == size_t( 0), " Molecule should have 0 valences");
        }
      }

      //! @brief Clone function
      //! @return pointer to new MolecularConformationShared
      MoleculeComplete *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the number of valences
      //! @return the number of valences (=0)
      size_t GetNumberValences() const
      {
        return 0;
      }

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

    }; // class MoleculeComplete

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULE_COMPLETE_H_ 
