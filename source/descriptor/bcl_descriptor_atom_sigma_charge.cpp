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
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param GET_CHARGE true if sigma charge is desired, false if sigma electronegativity is desired
    AtomSigmaCharge::AtomSigmaCharge( const bool &GET_CHARGE)
    : m_GetChargeOrElectronegativity( GET_CHARGE)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new AtomSigmaCharge
    AtomSigmaCharge *AtomSigmaCharge::Clone() const
    {
      return new AtomSigmaCharge( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomSigmaCharge::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomSigmaCharge::GetAlias() const
    {
      static const std::string s_name( "Atom_SigmaCharge"), s_en_name( "Atom_SigmaEN");
      return m_GetChargeOrElectronegativity ? s_name : s_en_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomSigmaCharge::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( m_Charges.IsEmpty())
      {
        RecalculateCharges();
      }
      if( m_GetChargeOrElectronegativity)
      {
        STORAGE( 0) = m_Charges( ELEMENT.GetPosition());
      }
      else
      {
        // copy the iterator
        STORAGE( 0) = ELEMENT->GetAtomType()->GetSigmaENFromCharge( m_Charges( ELEMENT.GetPosition()));
      }
    }

    //! @brief Recalculates the charges for the current molecule
    void AtomSigmaCharge::RecalculateCharges()
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);
      storage::Vector< float> q( molecule.GetNumberAtoms()), q_prev;

      // q = charges of each atom
      size_t q_id( 0);
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms_init( molecule.GetAtomsIterator());
        itr_atoms_init.NotAtEnd();
        ++itr_atoms_init, ++q_id
      )
      {
        q( q_id) = itr_atoms_init->GetCharge();
      }

      float alpha( 0.5);
      for
      (
        size_t n( 0);    // current iteration number
        n < 6;           // Adriana uses 10 iterations; we use 6 in accordance with the original gasteiger paper
        alpha /= 2.0, ++n
      )
      {
        // We only need access to the last iteration's q
        q_prev = q;

        storage::Vector< float>::const_iterator itr_prev_q( q_prev.Begin());
        storage::Vector< float>::iterator itr_q( q.Begin());

        // for each atom in the molecule
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms( molecule.GetAtomsIterator());
          itr_atoms.NotAtEnd();
          ++itr_atoms, ++itr_prev_q, ++itr_q
        )
        {
          if( itr_atoms->GetAtomType()->GetOrbitalENegPos() == 0.0)
          {
            continue;
          }

          const float chi_a( itr_atoms->GetAtomType()->GetSigmaENFromCharge( *itr_prev_q));

          // calculate q for bonded atoms
          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_bonds( itr_atoms->GetBonds().Begin()), itr_bonds_end( itr_atoms->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds
          )
          {
            const chemistry::AtomConformationalInterface &target_atom( itr_bonds->GetTargetAtom());
            if( target_atom.GetAtomType()->GetOrbitalENegPos() == 0.0)
            {
              continue;
            }

            // calculate chi_b
            const float chi_b
            (
              target_atom.GetAtomType()->GetSigmaENFromCharge( q_prev( molecule.GetAtomIndex( target_atom)))
            );

            // from the original gasteiger paper, the denominator needs to be the electronegativity of the cation
            // species of the atom type with the more electronegative atom
            if( chi_a >= chi_b)
            {
              // add the delta q for this bonded atom to the delta
              *itr_q += alpha * ( chi_b - chi_a) / itr_atoms->GetAtomType()->GetOrbitalENegPos();
            }
            else
            {
              // add the delta q for this bonded atom to the delta
              *itr_q += alpha * ( chi_b - chi_a) / target_atom.GetAtomType()->GetOrbitalENegPos();
            }
          }
        } // for each atom in the molecule
      } // for alpha
      m_Charges = q;
    } // Recalculate

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomSigmaCharge::SetObjectHook()
    {
      m_Charges.Resize( 0);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomSigmaCharge::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses PEOE to determine sigma-orbital "
        + std::string( m_GetChargeOrElectronegativity ? "partial charge" : "electronegativity")
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
