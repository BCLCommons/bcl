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
#include "chemistry/bcl_chemistry_fragment_ensemble.h"

// external includes - sorted alphabetically
#include "GraphMol/Atom.h"
#include "GraphMol/AtomIterators.h"
#include "GraphMol/Bond.h"
#include "GraphMol/BondIterators.h"
#include "GraphMol/Conformer.h"
#include "GraphMol/GraphMol.h"
#include "GraphMol/MolOps.h"
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
      //! @param ALLOW_UNDEFINED_ATOMS return null pointer if false
      //! @return the FragmentComplete constructed from the RDKit ROMol
      static std::shared_ptr< FragmentComplete> RDKitROMolToFragmentComplete
      (
        const ::RDKit::ROMol &RDKIT_MOL,
        const std::string &MOL_ID = "RDKit-to-BCL-molecule",
        const bool ALLOW_UNDEFINED_ATOMS = true
      );

      //! @brief converts an RDKit ROMol into a BCL FragmentComplete
      //! @param RDKIT_MOL the RDKit ROMol to be converted
      //! @param MOL_ID the molecule identifier on the new molecule
      //! @param ALLOW_UNDEFINED_ATOMS return null pointer if false
      //! @return the FragmentComplete constructed from the RDKit ROMol
      static std::shared_ptr< FragmentComplete> RDKitRWMolToFragmentComplete
      (
        const ::RDKit::RWMol &RDKIT_MOL,
        const std::string &MOL_ID = "RDKit-to-BCL-molecule",
        const bool ALLOW_UNDEFINED_ATOMS = true
      );

      //! @brief converts a BCL FragmentComplete into an RDKit RWMol
      //! @param MOL the BCL FragmentComplete to be converted
      //! @param UNDEFINED_SUBSTITUTE if the BCL atom is undefined (X), replace with this
      //! element type in the RDKit molecule; default is to leave undefined, which
      //! will skip this atom when building the new molecule
      //! @param SANITIZE perform RDKit-style standardization of molecule after an initial
      //! molecule is created from the BCL FragmentComplete
      //! @param SKIP_COORDS do not try to assign 3D coordinates from BCL object as a 3D conformer
      //! on the RDKit molecule; this means that we are effectively creating a molecule topology
      //! object from our BCL molecule; useful when just converting from FragmentComplete to SMILES
      //! @return the RDKit RWMol constructed from the FragmentComplete
      static std::shared_ptr< ::RDKit::RWMol> FragmentCompleteToRDKitRWMol
      (
        const FragmentComplete &MOL,
        const ElementType &UNDEFINED_SUBSTITUTE = GetElementTypes().e_Undefined,
        const bool SANITIZE = true,
        const bool SKIP_COORDS = false
      );

    private:

      //! @brief add coordinates from BCL conformers to RDKit molecule
      //! @details BCL and RDKit are organized differently with respect to
      //! storage of atomic coordinates. In the BCL, coordinates are a component
      //! of the AtomInfo object, which means each atom can independently have
      //! associated coordinate information. Each molecule contains both the
      //! topological and coordinate data, and thus a conformer ensemble is represented
      //! through a list of molecules (e.g., FragmentEnsemble). In contrast, RDKit
      //! molecule objects, such as RWMol, contain topology and property information,
      //! but no coordinate information. Conformers are stored on RDKit molecules as
      //! attributes, which means conformers need to be added to a molecule (rather than
      //! represented as individual molecules). This function will modify an RDKit molecule
      //! directly to add conformers based on an ensemble of BCL conformers.
      //! @param RDKIT_MOL the RDKit molecule to which conformers will be added
      //! @param MOL_ENS the BCL conformers that will be added to the RDKit molecule
      //! @param MAP_ATOMS if true, convert RDKit molecule to BCL molecule and perform
      //! substructure comparison to map atoms; default false assumes atoms are in the same
      //! order already, such as if the RDKit molecule was derived from the BCL molecule, or
      //! vice versa.
      //! TODO: make public, but need to finish the substructure mapping bit
      //! @param ATOM_MAP user-supplied map between RDKit and BCL molecule atoms; key is the
      //! BCL atom index and value is the RDKit atom index
      static void AddBCLConformersToRDKitRWMol
      (
        std::shared_ptr< ::RDKit::RWMol> RDKIT_MOL,
        const FragmentEnsemble &MOL_ENS,
        const bool MAP_ATOMS = false,
        const std::map< size_t, size_t> ATOM_MAP = std::map< size_t, size_t>()
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
