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

#ifndef BCL_CHEMISTRY_RDKIT_MOL_UTILS_H_
#define BCL_CHEMISTRY_RDKIT_MOL_UTILS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"

// external includes - sorted alphabetically
#include "GraphMol/Atom.h"
#include "GraphMol/AtomIterators.h"
#include "GraphMol/Bond.h"
#include "GraphMol/BondIterators.h"
#include "GraphMol/Conformer.h"
#include "GraphMol/GraphMol.h"
#include "GraphMol/ROMol.h"
#include "Geometry/point.h"
#include "RDGeneral/RDProps.h"
#include "RDGeneral/types.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RdkitMolUtils
    //! @brief This class contains functions to interconvert RDKit and BCL molecules
    //!
    //! @see @link example_chemistry_rdkit_mol_utils.cpp @endlink
    //! @author brownbp1
    //! @date Aug 26, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RdkitMolUtils :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief delete the default constructor
      RdkitMolUtils() = delete;

      //! @brief virtual copy constructor
      RdkitMolUtils *Clone() const;

      //! @brief returns class name
      //! the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief converts an RDKit ROMol into a BCL FragmentComplete
      //! @param RDKIT_MOL the RDKit ROMol to be converted
      //! @param MOL_ID the molecule identifier on the new molecule
      //! @return the FragmentComplete constructed from the RDKit ROMol
      static std::shared_ptr< FragmentComplete> RDKitROMolToFragmentComplete
      (
        const ::RDKit::ROMol &RDKIT_MOL,
        const std::string &MOL_ID = "RDKit-to-BCL-molecule"
      );

      //! @brief converts an RDKit ROMol into a BCL FragmentComplete
      //! @param RDKIT_MOL the RDKit ROMol to be converted
      //! @param MOL_ID the molecule identifier on the new molecule
      //! @return the FragmentComplete constructed from the RDKit ROMol
      static std::shared_ptr< FragmentComplete> RDKitRWMolToFragmentComplete
      (
        const ::RDKit::RWMol &RDKIT_MOL,
        const std::string &MOL_ID = "RDKit-to-BCL-molecule"
      );

      //! @brief converts a BCL FragmentComplete into an RDKit ROMol
      //! @param MOL the BCL FragmentComplete to be converted
      //! @param MOL_ID the molecule identifier on the new molecule
      //! @return the RDKit ROMol constructed from the FragmentComplete
      static std::shared_ptr< ::RDKit::ROMol> FragmentCompleteToRDKitROMol
      (
        const FragmentComplete &MOL,
        const std::string &MOL_ID = "BCL-to-RDKit-molecule"
      );

      //! @brief converts a BCL FragmentComplete into an RDKit RWMol
      //! @param MOL the BCL FragmentComplete to be converted
      //! @param MOL_ID the molecule identifier on the new molecule
      //! @return the RDKit RWMol constructed from the FragmentComplete
      static std::shared_ptr< ::RDKit::RWMol> FragmentCompleteToRDKitRWMol
      (
        const FragmentComplete &MOL,
        const std::string &MOL_ID = "BCL-to-RDKit-molecule"
      );

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

    private:

    }; // class RdkitMolUtils

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_RDKIT_MOL_UTILS_H_
