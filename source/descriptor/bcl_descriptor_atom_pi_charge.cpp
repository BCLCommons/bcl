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
#include "descriptor/bcl_descriptor_atom_pi_charge.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"
#include "graph/bcl_graph_connectivity.h"
#include "linal/bcl_linal_symmetric_eigensolver.h"

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
    AtomPiCharge::AtomPiCharge( const bool &GET_CHARGE)
    : m_GetChargeOrElectronegativity( GET_CHARGE)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new AtomPiCharge
    AtomPiCharge *AtomPiCharge::Clone() const
    {
      return new AtomPiCharge( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomPiCharge::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomPiCharge::GetAlias() const
    {
      static const std::string s_name( "Atom_PiCharge"), s_en_name( "Atom_PiEN");
      return m_GetChargeOrElectronegativity ? s_name : s_en_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomPiCharge::Calculate
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
        STORAGE( 0) = ELEMENT->GetAtomType()->GetPiENFromCharge( m_Charges( ELEMENT.GetPosition()));
      }
    }

    //! @brief Recalculates the charges for the current molecule
    void AtomPiCharge::RecalculateCharges()
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);

      storage::Vector< size_t> conjugated_atom_indices;

      // find all the conjugated atoms and store their index, so we can extract the subgraph
      size_t index( 0);
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr( molecule.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr, ++index
      )
      {
        if( itr->GetAtomType()->IsConjugated())
        {
          conjugated_atom_indices.PushBack( index);
        }
      }

      // make the graph of the whole molecule, then get the subgraph containing only the conjugated component
      const graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> conjugated_graph
      (
        chemistry::ConformationGraphConverter::CreateGraphWithAtoms( molecule).GetSubgraph( conjugated_atom_indices)
      );

      // now go through and eliminate any bonds that are not from supported types

      //need a storage map of sigma charges for each atom
      storage::Map< util::SiPtr< const chemistry::AtomConformationalInterface>, float> sigma_charge_map;

      storage::Vector< float> pi_charges( molecule.GetNumberAtoms(), float( 0.0));

      AtomSigmaCharge sigma_charge_calculator;
      sigma_charge_calculator.SetObject( *this->GetCurrentObject());
      size_t count( 0);
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr( molecule.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr, ++count
      )
      {
        Iterator< chemistry::AtomConformationalInterface> itr_atom( itr);
        sigma_charge_map[ *itr] = sigma_charge_calculator( itr_atom)( 0);
      }

      // make a static map from pi charge type string to storage::VectorND< 2, float > containing the DM h/k parameters
      static const storage::Map< std::string, storage::VectorND< 2, float> >
        pi_charge_string_to_parameters_map( MakePiChargeTypeStringToPiChargeParameterMap());

      // make a list of graphs, 1 for each connected component in the conjugated subgraph
      storage::List< graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> >
        conjugated_systems( graph::Connectivity::GetComponentsAsNewGraphs( conjugated_graph));

      linal::Vector< float> pi_en( molecule.GetNumberAtoms(), 0.0);

      //initialize all atoms to have an electron density and pi charge of zero

      size_t conj_number( 0);
      // the first for loop iterates through different conjugated systems of a molecule
      for
      (
        storage::List< graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> >::const_iterator
          outer_list_itr( conjugated_systems.Begin()),
          outer_list_itr_end( conjugated_systems.End());
        outer_list_itr != outer_list_itr_end;
        ++outer_list_itr, ++conj_number
      )
      {
        const size_t conjugated_system_size( ( *outer_list_itr).GetSize());

        // ignore any graphs of size 1, since a single atom does not give rise to a pi-system
        if( conjugated_system_size <= size_t( 1))
        {
          continue;
        }
        //stores how many electrons should be in an atomic orbital for a given atom
        linal::Vector< float> electron_density( conjugated_system_size, 0.0);
        linal::Vector< float> nominal_electrons( conjugated_system_size, 0.0);

        linal::Matrix< float> conjugated_connectivity( conjugated_system_size, conjugated_system_size, 0.0);
        size_t count_a( 0);

        const storage::Vector< util::SiPtr< const chemistry::AtomConformationalInterface> > &atoms_in_system( outer_list_itr->GetVertices());

        // iterate through all atoms within the current conjugated system
        for
        (
          storage::Vector< util::SiPtr< const chemistry::AtomConformationalInterface> >::const_iterator
            atom_list_itr( atoms_in_system.Begin()),
            atom_list_itr_end( atoms_in_system.End());
          atom_list_itr != atom_list_itr_end;
          ++atom_list_itr, ++count_a
        )
        {
          conjugated_connectivity( count_a, count_a) = std::numeric_limits< float>::infinity();

          // nominal electrons should describe the electron distribution in our conjugated system

          // add sigma charge parameter to the diagonal to influence the correct coulomb integral
          const float sigma_charge_parameter( 0.5 * sigma_charge_map[ *atom_list_itr]);

          const chemistry::ElementType element_a( ( *atom_list_itr)->GetElementType());

          nominal_electrons( count_a) = ( *atom_list_itr)->GetAtomType()->GetMaxEContributionToPiSystem();

          if( element_a == chemistry::GetElementTypes().e_Carbon)
          {
            conjugated_connectivity( count_a, count_a) = -sigma_charge_parameter;
          }

          // to iterate through the atoms neighbors we set up a for loop that itr through the current pointer of atom_list _itr
          for
          (
            storage::Vector< size_t>::const_iterator
              atom_list_b_itr( outer_list_itr->GetNeighborIndices( count_a).Begin()),
              atom_list_b_itr_end( outer_list_itr->GetNeighborIndices( count_a).End());
            atom_list_b_itr != atom_list_b_itr_end;
            ++atom_list_b_itr
          )
          {
            // element types and bond types
            const chemistry::ConfigurationalBondTypeData &bond_type( *chemistry::ConfigurationalBondType( outer_list_itr->GetEdgeData( count_a, *atom_list_b_itr)));
            const chemistry::ElementType element_b( atoms_in_system( *atom_list_b_itr)->GetElementType());

            // get the pi string corresponding to the element types and simple bond type
            const std::string pi_string
            (
              GetPiChargeTypeString
              (
                element_a,
                bond_type.GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic),
                element_b
              )
            );

            storage::Map< std::string, storage::VectorND< 2, float> >::const_iterator
              pi_charge_string_itr( pi_charge_string_to_parameters_map.Find( pi_string));

            if( pi_charge_string_itr == pi_charge_string_to_parameters_map.End())
            {
              continue;
            }

            // extract the parameters from the charge string

            // use the resonance integral parameter
            const float resonance_integral( pi_charge_string_itr->second.Second());
            conjugated_connectivity( count_a, *atom_list_b_itr) = -resonance_integral;
            conjugated_connectivity( *atom_list_b_itr, count_a) = -resonance_integral;

            if( element_b == chemistry::GetElementTypes().e_Carbon)
            {
              const float conjugation_value( pi_charge_string_itr->second.First());

              // if we haven't found the nominal electrons for this atom yet or the new number of nominal electrons is lower
              conjugated_connectivity( count_a, count_a)
                = std::min( conjugated_connectivity( count_a, count_a), -conjugation_value - sigma_charge_parameter);
            }
          }

          // if a conjugated connectivity was found, set up the number of nominal electrons
          if( conjugated_connectivity( count_a, count_a) > std::numeric_limits< float>::max())
          {
            conjugated_connectivity( count_a, count_a) = -sigma_charge_parameter;
          }
        }

        // a temporary matrix to hold the eigenvectors

        // return the eigen values
        linal::SymmetricEigenSolver< float> eigensolver( conjugated_connectivity, true);
        linal::Vector< float> eigenvalues( eigensolver.GetSortedEigenvalues());
        linal::Matrix< float> eigenvectors( eigensolver.GetSortedEigenvectors());

        // compute effective epsilon = NAtoms in conjugated system * machine eta
        const float effective_epsilon( conjugated_system_size * std::numeric_limits< float>::epsilon());
        // zero-out nearly 0 values in the eigenvector matrix to prevent numerical artifacts from influencing the
        // ordering of equivalent orbitals / occupancies
        for
        (
          float *itr_eigen( eigenvectors.Begin()), *itr_eigen_end( eigenvectors.End());
          itr_eigen != itr_eigen_end;
          ++itr_eigen
        )
        {
          if( math::Absolute( *itr_eigen) < effective_epsilon)
          {
            *itr_eigen = 0.0;
          }
        }

        // transpose the eigenvectors for easier access
        eigenvectors.Transpose();

        // orbitals used will count the electrons the given atom type will give to the conjugated system
        size_t orbitals_used( nominal_electrons.Sum());

        for
        (
          size_t eigenvalue_rank( 0), n_eigenvalues( ( orbitals_used + 1) / 2),
                 eigenvalue_number( conjugated_system_size - 1);
          eigenvalue_rank < n_eigenvalues;
          ++eigenvalue_rank, orbitals_used -= 2, --eigenvalue_number
        )
        {
          // determine the number of electrons in this orbital
          const float n_electrons( std::min( orbitals_used, size_t( 2)));

          // add electron density (=probability of electron being at an atom) from coefficients of the atom orbitals
          // electron density = sum( c_i^2)

          // TODO When the eigen values present molecular orbitals that are degenerative in energy, the occupation
          // will vary. It will fill each orbital one at a time where the overall may not have 2 electrons in it. Also,
          // in any ion form of a conjugated system the occupancy can be changed based on if it is a cation or anion.
          // These two instances should both be accounted for.
          //
          // electron density function
          linal::VectorConstReference< float> eigenvector( conjugated_system_size, eigenvectors[ eigenvalue_number]);
          float *itr_electron_density( electron_density.Begin());

          for
          (
            const float *itr_coeff( eigenvector.Begin()), *itr_coeff_end( eigenvector.End());
            itr_coeff != itr_coeff_end;
            ++itr_coeff, ++itr_electron_density
          )
          {
            *itr_electron_density += n_electrons * math::Sqr( *itr_coeff);
          }
        }

        // store pi charge for the actual conjugated system
        // from electron density from Gasteiger equation. (nominal electrons - electron density)
        linal::Vector< float> pi_charge_local( nominal_electrons);
        pi_charge_local -= electron_density;

        float *itr_pi_charge_local( pi_charge_local.Begin());
        for
        (
          storage::Vector< util::SiPtr< const chemistry::AtomConformationalInterface> >::const_iterator
            atom_list_itr( atoms_in_system.Begin()),
            atom_list_itr_end( atoms_in_system.End());
          atom_list_itr != atom_list_itr_end;
          ++atom_list_itr, ++itr_pi_charge_local
        )
        {
          // save pi-charges in molecule vector
          pi_charges( molecule.GetAtomIndex( **atom_list_itr)) = *itr_pi_charge_local;
        }
      }

      m_Charges = pi_charges;
    } // Recalculate

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomPiCharge::SetObjectHook()
    {
      m_Charges.Resize( 0);
    }

    //! @return a unique string for the given connectivity information
    std::string AtomPiCharge::GetPiChargeTypeString
    (
      const chemistry::ElementType &ELEMENT_TYPE_1,
      const size_t &BOND_ORDER_AND_AROMATIC,
      const chemistry::ElementType &ELEMENT_TYPE_2
    )
    {
      std::ostringstream combiner;
      // the string will be of the syntax Lighter-element-symbol,Bond type id,HeavierElementSymbol
      if( ELEMENT_TYPE_1->GetAtomicNumber() <= ELEMENT_TYPE_2->GetAtomicNumber())
      {
        combiner << ELEMENT_TYPE_1->GetChemicalSymbol() << ','
                 << BOND_ORDER_AND_AROMATIC << ','
                 << ELEMENT_TYPE_2->GetChemicalSymbol();
      }
      else
      {
        combiner << ELEMENT_TYPE_2->GetChemicalSymbol() << ','
                 << BOND_ORDER_AND_AROMATIC << ','
                 << ELEMENT_TYPE_1->GetChemicalSymbol();
      }

      return combiner.str();
    } // GetPiChargeTypeString

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create the connectivity string to pi-charge parameter map
    storage::Map< std::string, storage::VectorND< 2, float> >
    AtomPiCharge::MakePiChargeTypeStringToPiChargeParameterMap()
    {
      static storage::Map< std::string, storage::VectorND< 2, float> > map;
      if( map.IsEmpty())
      {
        // compose the map using data from Dixon and Jurs, 1992, also
        // Abraham & Smith, J Comp Chem, 1987 (http://onlinelibrary.wiley.com/doi/10.1002/jcc.540090403/pdf)
        chemistry::ElementType c( chemistry::GetElementTypes().e_Carbon), n( chemistry::GetElementTypes().e_Nitrogen), o( chemistry::GetElementTypes().e_Oxygen);
        chemistry::ElementType f( chemistry::GetElementTypes().e_Fluorine);
        chemistry::ElementType cl( chemistry::GetElementTypes().e_Chlorine), br( chemistry::GetElementTypes().e_Bromine);

        // bonds to carbon
        map[ GetPiChargeTypeString( c, 2, c)] = storage::VectorND< 2, float >( 0.00, 1.0);
        map[ GetPiChargeTypeString( c, 4, c)] = storage::VectorND< 2, float >( 0.00, 1.0);
        map[ GetPiChargeTypeString( c, 1, c)] = storage::VectorND< 2, float >( 0.00, 1.0);

        // bonds to nitrogen
        map[ GetPiChargeTypeString( c, 2, n)] = storage::VectorND< 2, float >( 0.33, 0.72);
        map[ GetPiChargeTypeString( c, 4, n)] = storage::VectorND< 2, float >( 0.33, 0.72);
        map[ GetPiChargeTypeString( c, 1, n)] = storage::VectorND< 2, float >( 1.48, 0.71);

        // bonds to oxygen
        map[ GetPiChargeTypeString( c, 2, o)] = storage::VectorND< 2, float >( 0.31, 1.05);
        map[ GetPiChargeTypeString( c, 4, o)] = storage::VectorND< 2, float >( 1.69, 0.59);
        map[ GetPiChargeTypeString( c, 1, o)] = storage::VectorND< 2, float >( 1.69, 0.59);

        // bonds to halogens
        map[ GetPiChargeTypeString( c, 1, f)] = storage::VectorND< 2, float >( 1.67, 0.54);
        map[ GetPiChargeTypeString( c, 1, cl)] = storage::VectorND< 2, float >( 1.03, 0.29);
        map[ GetPiChargeTypeString( c, 1, br)] = storage::VectorND< 2, float >( 1.03, 0.29);
      }

      return map;
    } // MakePiChargeTypeStringToPiChargeParameterMap

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomPiCharge::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses Hueckel matrix to determine pi-orbital "
        + std::string( m_GetChargeOrElectronegativity ? "partial charge" : "electronegativity")
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
