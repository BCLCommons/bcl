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
#include "chemistry/bcl_chemistry_atom_one_four_interaction_score.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_molecule_one_four_clash_score.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeOneFourClashScore
    MoleculeOneFourClashScore *MoleculeOneFourClashScore::Clone() const
    {
      return new MoleculeOneFourClashScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeOneFourClashScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeOneFourClashScore::GetAlias() const
    {
      static const std::string s_name( "MoleculeOneFourClashScore");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeOneFourClashScore::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // make our molecule
      util::SiPtr< const chemistry::ConformationInterface> reference_mol( this->GetCurrentObject());
      chemistry::FragmentComplete mol( *reference_mol);

      // get atom vdw score object
      chemistry::AtomOneFourInteractionScore one_four_clash_scorer;
      STORAGE = float( one_four_clash_scorer( mol));
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeOneFourClashScore::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the sum of one-four clash scores in the molecule"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
