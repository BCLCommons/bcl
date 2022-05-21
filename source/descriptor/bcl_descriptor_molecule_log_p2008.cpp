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
#include "descriptor/bcl_descriptor_molecule_log_p2008.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeLogP2008
    MoleculeLogP2008 *MoleculeLogP2008::Clone() const
    {
      return new MoleculeLogP2008( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeLogP2008::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeLogP2008::GetAlias() const
    {
      static const std::string s_name( "LogP2008");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeLogP2008::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // prepare atom property used during this calculation
      AtomSigmaCharge sigma_charge_calculator;
      sigma_charge_calculator.SetObject( *this->GetCurrentObject());

      // This formula is from http://onlinelibrary.wiley.com/doi/10.1002/jps.21494/abstract 
      const float hetatm_coeff( -0.11); 
      const float carbon_coeff( 0.11); 
      const float const_term( 1.46);

      size_t n_carbons( 0);
      size_t n_heteroatom( 0);

      for
      (
        Iterator< chemistry::AtomConformationalInterface> itr_atoms( this->GetCurrentObject()->GetIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        const chemistry::AtomConformationalInterface &atom( *itr_atoms( 0));
        if( atom.GetElementType() == chemistry::GetElementTypes().e_Carbon)
        {
          ++n_carbons;
        }
        else // if( atom.GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
        {
          ++n_heteroatom;
        }
      }

      STORAGE = float( const_term + carbon_coeff * n_carbons + hetatm_coeff * n_heteroatom);
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeLogP2008::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates logp using the heuristic formula from http://pubs.acs.org/doi/full/10.1021/ci010315d"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
