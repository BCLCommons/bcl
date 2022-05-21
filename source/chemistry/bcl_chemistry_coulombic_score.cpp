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
#include "chemistry/bcl_chemistry_coulombic_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CoulombicScore::s_Instance
    (
      GetObjectInstances().AddInstance( new CoulombicScore())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    CoulombicScore::CoulombicScore( const descriptor::CheminfoProperty &ATOM_CHARGE, double WEIGHT)
    : m_Charge( ATOM_CHARGE),
      m_Weight( WEIGHT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CoulombicScore
    CoulombicScore *CoulombicScore::Clone() const
    {
      return new CoulombicScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CoulombicScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return propensity score for observing rotamers that exist in conformation
    double CoulombicScore::operator()( const FragmentComplete &MOLECULE) const
    {
      const linal::Vector< float> charges( m_Charge->CollectValuesOnEachElementOfObject( MOLECULE));

      float force_sum( 0.0);
      // iterate over all possible pairs of atoms
      // iterate properties and surface areas simultaneously
      linal::Vector< float>::const_iterator itr_charge_a( charges.Begin());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms_a( MOLECULE.GetAtomsIterator());
        itr_atoms_a.NotAtEnd();
        ++itr_atoms_a, ++itr_charge_a
      )
      {
        iterate::Generic< const AtomConformationalInterface> itr_atoms_b( itr_atoms_a);
        linal::Vector< float>::const_iterator itr_charge_b( itr_charge_a + 1);
        for( ++itr_atoms_b; itr_atoms_b.NotAtEnd(); ++itr_atoms_b, ++itr_charge_b)
        {
          // store distance between both atoms
          const float sqr_distance( linal::SquareDistance( itr_atoms_a->GetPosition(), itr_atoms_b->GetPosition()));
          if( sqr_distance < 0.1)
          {
            // overlapping atom
            continue;
          }
          force_sum += *itr_charge_b * *itr_charge_a / sqr_distance;
        }
      }
      return force_sum;
    } // Recalculate

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CoulombicScore::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CoulombicScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
