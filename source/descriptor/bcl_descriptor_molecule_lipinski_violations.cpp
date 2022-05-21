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
#include "descriptor/bcl_descriptor_molecule_lipinski_violations.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_cheminfo_properties.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @brief VARIANT the variant on the lipinski rule of 5 to use
    MoleculeLipinskiViolations::MoleculeLipinskiViolations
    (
      const MoleculeLipinskiViolations::LipinskiVariant &VARIANT
    ) :
      m_Variant( VARIANT)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeLipinskiViolations
    MoleculeLipinskiViolations *MoleculeLipinskiViolations::Clone() const
    {
      return new MoleculeLipinskiViolations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeLipinskiViolations::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeLipinskiViolations::GetAlias() const
    {
      static const std::string s_original( "LipinskiViolations");
      static const std::string s_veber( "LipinskiViolationsVeber");
      return m_Variant == e_Veber ? s_veber : s_original;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeLipinskiViolations::Calculate( linal::VectorReference< float> &STORAGE)
    {

      // the total number of violations
      size_t violations( 0);

      // the Veber variant of the Lipinski rule of 5
      if( m_Variant == e_Veber)
      {
        size_t n_rot( GetCheminfoProperties().calc_NRotBond->SumOverObject( *GetCurrentObject())( 0));
        float tpsa( GetCheminfoProperties().calc_TopologicalPolarSurfaceArea->SumOverObject( *GetCurrentObject())( 0));
        if( n_rot > 10)
        {
          ++violations;
        }
        if( tpsa > 140)
        {
          ++violations;
        }
      }
      else // Original Lipinski
      {
        float weight( GetCheminfoProperties().calc_MolWeight->SumOverObject( *GetCurrentObject())( 0));
        size_t h_bond_donors( GetCheminfoProperties().calc_HbondDonor->SumOverObject( *GetCurrentObject())( 0));
        size_t h_bond_acceptors( GetCheminfoProperties().calc_HbondAcceptor->SumOverObject( *GetCurrentObject())( 0));
        float log_p( GetCheminfoProperties().calc_LogP->SumOverObject( *GetCurrentObject())( 0));

        if( weight > 500)
        {
          ++violations;
        }
        if( h_bond_donors > 5)
        {
          ++violations;
        }
        if( h_bond_acceptors > 10)
        {
          ++violations;
        }
        if( log_p > 5)
        {
          ++violations;
        }
      }
      STORAGE = violations;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeLipinskiViolations::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "calculates how many times a molecule violates Lipinski's Rule of 5");
      if( m_Variant == e_Veber)
      {
        parameters.SetClassDescription
        (
          "calculates how many times a molecule violates Veber's variant of the Lipinski's Rule of 5 "
          "(<10 rotatable bonds, Polar SA < 140 A^2; see J. Med. Chem., 2002, 45 (12), pp 2615–2623)"
        );
      }
      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
