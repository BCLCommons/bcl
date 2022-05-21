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
#include "descriptor/bcl_descriptor_atom_vcharge.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    AtomVcharge *AtomVcharge::Clone() const
    {
      return new AtomVcharge( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomVcharge::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomVcharge::GetAlias() const
    {
      static const std::string s_Name( "Atom_Vcharge"), s_NameV2( "Atom_VchargeV2");
      return m_ExtrapolationVersion ? s_NameV2 : s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    void AtomVcharge::RecalculateCharges()
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);

      // set up the return value
      m_Vcharge.Resize( molecule.GetNumberAtoms(), 0.0);

      // make a map from unique string identifier for each vcharge atom type to the electronegativity and hardness
      static const storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >
        vcharge_type_string_to_vcharge_params( MakeVchargeTypeStringToVchargeParametersMap());

      // compute vcharge parameters (electronegativity & hardness)
      storage::Vector< storage::Pair< float, float> > vcharge_params;
      vcharge_params.AllocateMemory( molecule.GetNumberAtoms());
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms( molecule.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        // tally up the number of bonds of each type
        size_t number_of_single_bonds( itr_atoms->GetNumberValenceBonds());
        size_t number_of_double_bonds( 0);
        size_t number_of_triple_bonds( 0);
        bool is_aromatic( false);
        bool is_in_ring( false);
        for
        (
          storage::Vector< chemistry::BondConformational>::const_iterator
            itr_bonds( itr_atoms->GetBonds().Begin()), itr_bonds_end( itr_atoms->GetBonds().End());
          itr_bonds != itr_bonds_end;
          ++itr_bonds
        )
        {
          const chemistry::ConfigurationalBondType bond_type( itr_bonds->GetBondType());

          // if the bond is aromatic
          if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic)
          {
            is_aromatic = true;
          }
          else if( bond_type->IsBondInRing())
          {
            is_in_ring = true;
          }

          if( bond_type->GetNumberOfElectrons() == 2) // if the bond is nominally a single bond
          {
            ++number_of_single_bonds;
          }
          else if( bond_type->GetNumberOfElectrons() == 4) // if the bond is nominally a float  bond
          {
            ++number_of_double_bonds;
          }
          else if( bond_type->GetNumberOfElectrons() == 6) // if the bond is nominally a triple bond
          {
            ++number_of_triple_bonds;
          }
        }

        // get the Gasteiger atom type
        const chemistry::AtomType atom_type( itr_atoms->GetAtomType());

        // determine the "special feature" character for the atom, 'A' for aromatic ring, 'P' for planar nitrogen with
        // 2-3 single bonds in a ring
        char special_feature( is_aromatic ? 'A' : ' ');

        if // only determine planarity for nitrogen with 2-3 single bonds in a non-aromatic ring
        (
          itr_atoms->GetElementType() == chemistry::GetElementTypes().e_Nitrogen
          && !is_aromatic
          && is_in_ring
          && number_of_double_bonds == size_t( 0)
          && number_of_triple_bonds == size_t( 0)
          && number_of_single_bonds >= size_t( 2)
          && number_of_single_bonds <= size_t( 3)
          && atom_type->GetHybridOrbitalType() != chemistry::GetHybridOrbitalTypes().e_SP3
        )
        {
          special_feature = 'P';
        }

        // make a string encoding all the information we now have about the type
        const std::string vcharge_string
        (
          GetVchargeTypeString( number_of_single_bonds, number_of_double_bonds, number_of_triple_bonds, special_feature, m_ExtrapolationVersion)
        );

        // key for the value in the vcharge parameter map
        std::pair< chemistry::ElementType, std::string> key( std::make_pair( itr_atoms->GetElementType(), vcharge_string));

        // look for the string in the map
        storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >::const_iterator map_itr
        (
          vcharge_type_string_to_vcharge_params.Find( key)
        );

        // if the type is not found, use the average values for the period
        if( map_itr == vcharge_type_string_to_vcharge_params.End() && m_ExtrapolationVersion == size_t( 1))
        {
          key.second[ key.second.size() - 1] = '0';
          map_itr = vcharge_type_string_to_vcharge_params.Find( key);
        }

        // if the type is not found, use the average values for the period
        if( map_itr == vcharge_type_string_to_vcharge_params.End())
        {
          key.first = chemistry::ElementType();
          key.second = util::Format()( itr_atoms->GetElementType()->GetPeriod()) + util::Format()( m_ExtrapolationVersion);
          map_itr = vcharge_type_string_to_vcharge_params.Find( key);
        }

        // if the string was found in the map, set the corresponding place in the vector to this value
        if( map_itr != vcharge_type_string_to_vcharge_params.End())
        {
          vcharge_params.PushBack( map_itr->second);
        }
        else
        {
          BCL_MessageCrt
          (
            "Vcharge atom type not found: " + key.first->GetChemicalSymbol() + vcharge_string
          );
          vcharge_params.PushBack( storage::Pair< float, float>( 0, 0));
        }
      }

      // declaring doubles to store the numerator and denominator of the langrangrian multiplier, which is unique for
      // the molecule
      float  numerator_lagr_sum( 0);
      float  denominator_lagr_sum( 0);

      // declaring vector for the vcharge calculation
      storage::Vector< float> hardness_vec( molecule.GetNumberAtoms());
      storage::Vector< float>::iterator itr_vcharge( m_Vcharge.Begin()), itr_hardness( hardness_vec.Begin());
      storage::Vector< storage::Pair< float, float> >::const_iterator itr_eneg_hardness( vcharge_params.Begin());

      // calculate the new electronegativities for each atom
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms( molecule.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++itr_vcharge, ++itr_eneg_hardness, ++itr_hardness
      )
      {
        // assigning the electronegativity and hardness values to their variables
        const float  electro_old( itr_eneg_hardness->First());
        const float  hardness( itr_eneg_hardness->Second());

        // store the hardness of each atom for use in calculating the vcharges
        *itr_hardness = hardness;

        // initializing the new electronegativity by the old because the old electronegativity is included in the sum of the new
        float electro_new( electro_old);
        // add up the influence of each bond type for the new electronegativity
        for
        (
          storage::Vector< chemistry::BondConformational>::const_iterator
            itr_bonds( itr_atoms->GetBonds().Begin()), itr_bonds_end( itr_atoms->GetBonds().End());
          itr_bonds != itr_bonds_end;
          ++itr_bonds
        )
        {
          // get the bond type for the connected atoms
          const chemistry::ConfigurationalBondType bond_type( itr_bonds->GetBondType());

          // get the vcharge parameters for the connected atom
          const chemistry::AtomConformationalInterface &connected_atom( itr_bonds->GetTargetAtom());
          const storage::Pair< float, float> &eneg_hardness_connected_atom
          (
            vcharge_params( molecule.GetAtomIndex( connected_atom))
          );
          // retrieving the electronegativity of the new atom
          const float electro_atom( eneg_hardness_connected_atom.First());

          // calculating the difference between the electronegativities, take the absolute value and raise to the beta
          float bond_influence = std::abs( electro_old - electro_atom);
          bond_influence = std::pow( bond_influence, float( 1.378));

          // determine the sign of the term from above
          if( electro_old < electro_atom)
          {
            bond_influence *= -1;
          }

          // multiply by the correct coefficient depending on the type of bond the atoms form
          if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic) // if the bond is aromatic
          {
            bond_influence *= 0.86;
          }
          else if( bond_type->GetNumberOfElectrons() == 2) // if the bond is nominally a single bond
          {
            bond_influence *= 1.00;
          }
          else if( bond_type->GetNumberOfElectrons() == 4) // if the bond is nominally a float  bond
          {
            bond_influence *= 1.74;
          }
          else if( bond_type->GetNumberOfElectrons() == 6) // if the bond is nominally a triple bond
          {
            bond_influence *= 1.67;
          }

          // add to the new electronegativity term
          electro_new += bond_influence;

          // account for farther atoms (1-3 bonds)
          // term for electronegativity taking into account the 1-3 bonding relationships
          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_bonds_connected_atom( connected_atom.GetBonds().Begin()),
              itr_bonds_connected_atom_end( connected_atom.GetBonds().End());
            itr_bonds_connected_atom != itr_bonds_connected_atom_end;
            ++itr_bonds_connected_atom
          )
          {
            // get the next connected atom
            const chemistry::AtomConformationalInterface &next_connected_atom( itr_bonds_connected_atom->GetTargetAtom());

            // get the vcharge parameters
            const storage::Pair< float, float> &eneg_hardness_next_connected_atom
            (
              vcharge_params( molecule.GetAtomIndex( next_connected_atom))
            );

            // retrieving the electronegativity of the new atom
            const float electro_atom2( eneg_hardness_next_connected_atom.First());

            // calculating the difference between the electronegativities, take the absolute value and raise to the beta
            float bond_influence2( std::abs( electro_old - electro_atom2));
            bond_influence2 = std::pow( bond_influence2, float( 1.378));

            // determine the sign of the term from above
            if( electro_old < electro_atom2)
            {
              bond_influence2 *= -1;
            }

            // alpha coefficient for atoms that are 2 bonds away
            bond_influence2 *= 0.057;

            electro_new -= bond_influence2;
          }
        }

        // enter the new electronegativity into its vector
        *itr_vcharge = electro_new;

        // Determine the numerator of the lagrange multiplier for each atom
        const float numerator_lagr( electro_new / hardness);

        // Determine the denominator of the lagrange multiplier for each atom
        const float denominator_lagr( 1.0 / hardness);

        // Sum up all the numerators
        numerator_lagr_sum += numerator_lagr;

        // Sum up all the denominators
        denominator_lagr_sum += denominator_lagr;
      }

      // Determine the lagrange multiplier
      const float lagrange( numerator_lagr_sum / denominator_lagr_sum);

      // calculate the actual vcharge for each atom
      iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms( molecule.GetAtomsIterator());
      for
      (
        storage::Vector< float>::iterator
          itr_electro_new( m_Vcharge.Begin()), itr_electro_new_end( m_Vcharge.End()), itr_hardness( hardness_vec.Begin());
        itr_electro_new != itr_electro_new_end;
        ++itr_electro_new, ++itr_hardness, ++itr_atoms
      )
      {
        // note: the factor of two here is not in the original paper due to a misprint.  See
        //! @see https://structbio.vanderbilt.edu/twiki/bin/view/MeilerLab/InformationAboutTheNewVchargeImplementation
        // for details, specifically the emails to/from the author of VCharge: Dr. Gilson
//        const float electo_old( *itr_electro_new);
        *itr_electro_new = ( lagrange - *itr_electro_new) / ( 2.0 * ( *itr_hardness));
//        if( math::Absolute( *itr_electro_new) > 1.0)
//        {
//          BCL_MessageStd
//          (
//            "AtomType: " + itr_atoms->GetAtomType().GetName() + " Version: " + util::Format()( m_ExtrapolationVersion)
//            + " VChg: " + util::Format()( *itr_electro_new)
//            + " Eneg: " + util::Format()( vcharge_params( itr_atoms.GetPosition()).First())
//            + " Hardness: " + util::Format()( vcharge_params( itr_atoms.GetPosition()).Second())
//            + " e_old: " + util::Format()( electo_old)
//            + " h_old: " + util::Format()( *itr_hardness)
//          );
//        }
      }
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomVcharge::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( m_Vcharge.IsEmpty())
      {
        RecalculateCharges();
      }
      STORAGE( 0) = m_Vcharge( ELEMENT.GetPosition());
    } // Recalculate

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomVcharge::SetObjectHook()
    {
      m_Vcharge.Resize( 0);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @return a unique string for the given connectivity information
    //! @param SPECIAL_FEATURE One of the following letters: A: in aromatic ring, P: in planar ring, (space): otherwise
    //! @note: Planar ring feature is only applicable for N with 2-3 single bonds, no float  bonds
    std::string AtomVcharge::GetVchargeTypeString
    (
      const size_t &NUMBER_SINGLE_BONDS,
      const size_t &NUMBER_DOUBLE_BONDS,
      const size_t &NUMBER_TRIPLE_BONDS,
      const char &SPECIAL_FEATURE,
      const size_t &EXTRAPOLATION_VERSION
    )
    {
      std::ostringstream combiner;
      combiner << NUMBER_SINGLE_BONDS << NUMBER_DOUBLE_BONDS << NUMBER_TRIPLE_BONDS << SPECIAL_FEATURE << EXTRAPOLATION_VERSION;
      return combiner.str();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomVcharge::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Partial charges computed using vcharge 2003 algorithm and parameters, see http://pubs.acs.org/doi/full/10.1021/ci034148o"
      );
      return parameters;
    }

    //! @brief estimate the electronegativity and hardness of a type in an aromatic ring given the base values
    //! @param BASE_ELECTRONEGATIVITY_HARDNESS pair of electronegativity and hardness
    //! @return the estimated electronegativity and hardness
    storage::Pair< float, float> EstimateAromaticTypeParams
    (
      const storage::Pair< float, float> &BASE_ELECTRONEGATIVITY_HARDNESS
    )
    {
      // these values are estimated from a linear fit of all types where we have both electronegativity and hardness
      // for the aromatic and non-aromatic types (except N_TrTrTrPi, which is an outlier), the following best fit equations are obtained
      // Electronegativity of aromatic state = 1.06 * Electronegativity of base state - 2.22 (R^2 = 0.8338)
      // Hardness of aromatic state = 2.56 * hardness of aromatic state - 94.7 (R^2 = 0.967)
      return
        storage::Pair< float, float>
        (
          std::max( 1.06 * BASE_ELECTRONEGATIVITY_HARDNESS.First() - 2.22, 0.0),
          std::max( 2.56 * BASE_ELECTRONEGATIVITY_HARDNESS.Second() - 94.7, 0.0)
        );
    }

    //! @brief estimate the electronegativity and hardness of a period 3 atom type from the corresponding period 2 atom type's params
    //! @param BASE_ELECTRONEGATIVITY_HARDNESS pair of electronegativity and hardness
    //! @return the estimated electronegativity and hardness
    storage::Pair< float, float> EstimatePeriod3TypeParamsFromPeriod2Params
    (
      const storage::Pair< float, float> &BASE_ELECTRONEGATIVITY_HARDNESS
    )
    {
      // performing a linear fit of all types where we have both a period 2 element and a period 3
      // element with the same environment, the following best fit equations are obtained
      // Electronegativity of period 3 type = 2.16 * Electronegativity of period 2 type - 60.5 (R^2 = 0.89)
      // Hardness of period 3 type = 1.88 * hardness of aromatic state - 35.7 (R^2 = 0.921)
      return
        storage::Pair< float, float>
        (
          std::max( 2.16 * BASE_ELECTRONEGATIVITY_HARDNESS.First() - 60.5, 0.0),
          std::max( 1.88 * BASE_ELECTRONEGATIVITY_HARDNESS.Second() - 35.7, 0.0)
        );
    }

    //! @brief estimate the electronegativity and hardness of a type in an aromatic ring given the base values
    //! @param BASE_ELECTRONEGATIVITY_HARDNESS pair of electronegativity and hardness
    //! @return the estimated electronegativity and hardness
    storage::Pair< float, float> EstimateAromaticTypeParamsV2
    (
      const storage::Pair< float, float> &BASE_ELECTRONEGATIVITY_HARDNESS,
      const chemistry::ElementType &ELEMENT_TYPE
    )
    {
      // these values are estimated from a linear fit of all types where we have both electronegativity and hardness
      // for the aromatic and non-aromatic types (except N_TrTrTrPi, which is an outlier), the following best fit equations are obtained
      // The mins and maxes are used here to ensure that the results stay bound within the range of observed hardnesses
      // and electronegativities
      return
        storage::Pair< float, float>
        (
          std::min
          (
            std::max // 0.904 correlation
            (
              1.795 * BASE_ELECTRONEGATIVITY_HARDNESS.First()
              - 0.1 * BASE_ELECTRONEGATIVITY_HARDNESS.Second()
              - 4.103 * ELEMENT_TYPE->GetMainGroup()
              - 2.35,
              20.0
            ),
            60.0
          ),
          std::min
          (
            std::max // 0.988 correlation
            (
              2.347 * BASE_ELECTRONEGATIVITY_HARDNESS.First()
              + 1.315 * BASE_ELECTRONEGATIVITY_HARDNESS.Second()
              + 2.899 * ELEMENT_TYPE->GetMainGroup()
              - 147.29,
              10.0
            ),
            190.0
          )
        );
    }

    //! @brief estimate the electronegativity and hardness of a period 3 atom type from the corresponding period 2 atom type's params
    //! @param BASE_ELECTRONEGATIVITY_HARDNESS pair of electronegativity and hardness
    //! @return the estimated electronegativity and hardness
    storage::Pair< float, float> EstimatePeriod3TypeParamsFromPeriod2ParamsV2
    (
      const storage::Pair< float, float> &BASE_ELECTRONEGATIVITY_HARDNESS,
      const chemistry::ElementType &ELEMENT_TYPE
    )
    {
      // performing a linear fit of all types where we have both a period 2 element and a period 3
      // element with the same environment, the following best fit equations are obtained
      return
        storage::Pair< float, float>
        (
          std::min
          (
            std::max // 0.96 correlation
            (
              (
                1.59 * BASE_ELECTRONEGATIVITY_HARDNESS.First()
                - 0.0454 * BASE_ELECTRONEGATIVITY_HARDNESS.Second()
                - 0.766 * ELEMENT_TYPE->GetMainGroup()
                - 26.13
              ),
              20.0
            ),
            60.0
          ),
          std::min
          (
            std::max // 0.971 correlation
            (
              2.02 * BASE_ELECTRONEGATIVITY_HARDNESS.First()
              + 0.5193 * BASE_ELECTRONEGATIVITY_HARDNESS.Second()
              - 9.227 * ELEMENT_TYPE->GetMainGroup()
              - 22.0,
              10.0
            ),
            190.0
          )
        );
    }

    //! @brief estimate the electronegativity and hardness of a period 2 atom type from the corresponding period 3 atom type's params
    //! @param BASE_ELECTRONEGATIVITY_HARDNESS pair of electronegativity and hardness
    //! @return the estimated electronegativity and hardness
    storage::Pair< float, float> EstimatePeriod2TypeParamsFromPeriod3Params
    (
      const storage::Pair< float, float> &BASE_ELECTRONEGATIVITY_HARDNESS
    )
    {
      // This is the inverse equation from EstimatePeriod3TypeParamsFromPeriod2Params
      return
        storage::Pair< float, float>
        (
          ( BASE_ELECTRONEGATIVITY_HARDNESS.First() + 60.5) / 2.16,
          ( BASE_ELECTRONEGATIVITY_HARDNESS.Second() + 35.7) / 1.88
        );
    }

    //! @brief estimate the electronegativity and hardness of a period 2 atom type from the corresponding period 3 atom type's params
    //! @param BASE_ELECTRONEGATIVITY_HARDNESS pair of electronegativity and hardness
    //! @return the estimated electronegativity and hardness
    storage::Pair< float, float> EstimatePeriod2TypeParamsFromPeriod3ParamsV2
    (
      const storage::Pair< float, float> &BASE_ELECTRONEGATIVITY_HARDNESS,
      const chemistry::ElementType &ELEMENT_TYPE
    )
    {
      // This is the inverse equation from EstimatePeriod3TypeParamsFromPeriod2Params
      return
        storage::Pair< float, float>
        (
          std::min
          (
            std::max // 0.96 correlation
            (
              (
                BASE_ELECTRONEGATIVITY_HARDNESS.First()
                + 0.0454 * BASE_ELECTRONEGATIVITY_HARDNESS.Second()
                + 0.766 * ELEMENT_TYPE->GetMainGroup()
                + 26.13
              ) / 1.59,
              20.0
            ),
            60.0
          ),
          std::min
          (
            std::max // 0.971 correlation
            (
              (
                BASE_ELECTRONEGATIVITY_HARDNESS.First()
                - 0.5193 * BASE_ELECTRONEGATIVITY_HARDNESS.Second()
                + 9.227 * ELEMENT_TYPE->GetMainGroup()
                + 22.0
              ) / 2.02,
              10.0
            ),
            190.0
          )
        );
    }

    //! @brief create the map from element type & bond counts to Vcharge Electronegativity and Hardness
    //! @return the actual Map with the Vcharge Electronegativity and Hardness
    storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >
    AtomVcharge::MakeVchargeTypeStringToVchargeParametersMap()
    {
      static storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> > map;
      if( map.IsEmpty())
      {
        // compose the map using data from http://pubs.acs.org/doi/full/10.1021/ci034148o
        //                        Sym  1x 2x 3x  Special*                     Electronegativity,Hardness
        // Special: One of the following letters: A: in aromatic ring, P: in planar ring, (space): otherwise
        //          Planar ring feature is only applicable for N with 2-3 single bonds, no float  bonds
        map[ std::make_pair( chemistry::GetElementTypes().e_Hydrogen,   GetVchargeTypeString( 1, 0, 0, ' ', 0))] = storage::Pair< float , float >( 27.4, 73.9);
        map[ std::make_pair( chemistry::GetElementTypes().e_Carbon,     GetVchargeTypeString( 4, 0, 0, ' ', 0))] = storage::Pair< float , float >( 30.8, 78.4);
        map[ std::make_pair( chemistry::GetElementTypes().e_Carbon,     GetVchargeTypeString( 2, 1, 0, ' ', 0))] = storage::Pair< float , float >( 33.6, 76.4);
        map[ std::make_pair( chemistry::GetElementTypes().e_Carbon,     GetVchargeTypeString( 0, 2, 0, ' ', 0))] = storage::Pair< float , float >( 37, 65.3);
        map[ std::make_pair( chemistry::GetElementTypes().e_Carbon,     GetVchargeTypeString( 1, 0, 1, ' ', 0))] = storage::Pair< float , float >( 40, 98.5);
        map[ std::make_pair( chemistry::GetElementTypes().e_Carbon,     GetVchargeTypeString( 2, 1, 0, 'A', 0))] = storage::Pair< float , float >( 34.6, 84.7);
        map[ std::make_pair( chemistry::GetElementTypes().e_Oxygen,     GetVchargeTypeString( 2, 0, 0, ' ', 0))] = storage::Pair< float , float >( 45.7, 92.6);
        map[ std::make_pair( chemistry::GetElementTypes().e_Oxygen,     GetVchargeTypeString( 0, 1, 0, ' ', 0))] = storage::Pair< float , float >( 49.5, 86.1);
        map[ std::make_pair( chemistry::GetElementTypes().e_Oxygen,     GetVchargeTypeString( 1, 0, 0, ' ', 0))] = storage::Pair< float , float >( 49.3, 25);
        map[ std::make_pair( chemistry::GetElementTypes().e_Oxygen,     GetVchargeTypeString( 2, 0, 0, 'A', 0))] = storage::Pair< float , float >( 45.9, 137);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 3, 0, 0, ' ', 0))] = storage::Pair< float , float >( 44, 87.6);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 3, 0, 0, 'P', 0))] = storage::Pair< float , float >( 43.6, 94.4);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 1, 1, 0, ' ', 0))] = storage::Pair< float , float >( 44, 72.7);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 0, 0, 1, ' ', 0))] = storage::Pair< float , float >( 57, 111);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 4, 0, 0, ' ', 0))] = storage::Pair< float , float >( 42.8, 188);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 2, 1, 0, ' ', 0))] = storage::Pair< float , float >( 37.6, 41.5);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 0, 2, 0, ' ', 0))] = storage::Pair< float , float >( 24, 104);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 1, 0, 1, ' ', 0))] = storage::Pair< float , float >( 39.4, 29.7);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 3, 0, 0, 'A', 0))] = storage::Pair< float , float >( 43.4, 136);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 1, 1, 0, 'A', 0))] = storage::Pair< float , float >( 53, 102);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 2, 1, 0, 'A', 0))] = storage::Pair< float , float >( 38.7, 8.64);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 0, 1, 0, ' ', 0))] = storage::Pair< float , float >( 31.9, 129);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 2, 0, 0, ' ', 0))] = storage::Pair< float , float >( 28.3, 20.9);
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen,   GetVchargeTypeString( 2, 0, 0, 'P', 0))] = storage::Pair< float , float >( 43.6, 0.176);
        map[ std::make_pair( chemistry::GetElementTypes().e_Chlorine,   GetVchargeTypeString( 1, 0, 0, ' ', 0))] = storage::Pair< float , float >( 37.6, 53.5);
        map[ std::make_pair( chemistry::GetElementTypes().e_Fluorine,   GetVchargeTypeString( 1, 0, 0, ' ', 0))] = storage::Pair< float , float >( 45.2, 96.8);
        map[ std::make_pair( chemistry::GetElementTypes().e_Bromine,    GetVchargeTypeString( 1, 0, 0, ' ', 0))] = storage::Pair< float , float >( 40.1, 75.3);
        map[ std::make_pair( chemistry::GetElementTypes().e_Sulfur,     GetVchargeTypeString( 2, 0, 0, ' ', 0))] = storage::Pair< float , float >( 37.4, 69.1);
        map[ std::make_pair( chemistry::GetElementTypes().e_Sulfur,     GetVchargeTypeString( 3, 0, 0, ' ', 0))] = storage::Pair< float , float >( 31.8, 93.9);
        map[ std::make_pair( chemistry::GetElementTypes().e_Sulfur,     GetVchargeTypeString( 2, 1, 0, ' ', 0))] = storage::Pair< float , float >( 35.8, 93.1);
        map[ std::make_pair( chemistry::GetElementTypes().e_Sulfur,     GetVchargeTypeString( 2, 2, 0, ' ', 0))] = storage::Pair< float , float >( 31.7, 83.2);
        map[ std::make_pair( chemistry::GetElementTypes().e_Sulfur,     GetVchargeTypeString( 2, 0, 0, 'A', 0))] = storage::Pair< float , float >( 33.8, 88.9);
        map[ std::make_pair( chemistry::GetElementTypes().e_Sulfur,     GetVchargeTypeString( 0, 1, 0, ' ', 0))] = storage::Pair< float , float >( 47.5, 74.3);
        map[ std::make_pair( chemistry::GetElementTypes().e_Sulfur,     GetVchargeTypeString( 1, 0, 0, ' ', 0))] = storage::Pair< float , float >( 44.5, 24.8);
        map[ std::make_pair( chemistry::GetElementTypes().e_Phosphorus, GetVchargeTypeString( 3, 0, 0, ' ', 0))] = storage::Pair< float , float >( 37.9, 72.5);
        map[ std::make_pair( chemistry::GetElementTypes().e_Phosphorus, GetVchargeTypeString( 4, 0, 0, ' ', 0))] = storage::Pair< float , float >( 29.6, 108.5);
        map[ std::make_pair( chemistry::GetElementTypes().e_Phosphorus, GetVchargeTypeString( 3, 1, 0, ' ', 0))] = storage::Pair< float , float >( 33.0, 86.6);
        map[ std::make_pair( chemistry::GetElementTypes().e_Iodine,     GetVchargeTypeString( 1, 0, 0, ' ', 0))] = storage::Pair< float , float >( 41.3, 109);
        map[ std::make_pair( chemistry::GetElementTypes().e_Iodine,     GetVchargeTypeString( 2, 0, 0, ' ', 0))] = storage::Pair< float , float >( 34.1, 10.8);

        // the above values represent all the values given in the paper
        // However, about 10% of organic molecules have other atom types, so it is necessary to approximate their
        // parameters, using linear correlations, very similarly to how most of the gasteiger partial charge values
        // were calculated

        // Compute average hardness and electronegativity over all periods
        storage::Vector< math::RunningAverage< float > > ave_vch_eneg_by_period( 7);
        storage::Vector< math::RunningAverage< float > > ave_vch_hard_by_period( 7);

        // Compute average hardness and electronegativity values
        for
        (
          storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >::const_iterator
            itr( map.Begin()), itr_end( map.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the zero-indexed period for this element type
          const size_t period( itr->first.first->GetPeriod() - 1);

          // add the electronegativity and hardness to the running averages
          ave_vch_eneg_by_period( period) += itr->second.First();
          ave_vch_hard_by_period( period) += itr->second.Second();
        }

        // for any zero electronegativity, use the ratio between the last two periods, times the last period's value
        for( size_t period( 2), max_period( 7); period < max_period; ++period)
        {
          if( ave_vch_eneg_by_period( period).GetWeight() == float( 0.0))
          {
            float ratio_last_two_periods( ave_vch_eneg_by_period( period - 1) / ave_vch_eneg_by_period( period - 2));
            ave_vch_eneg_by_period( period) += ave_vch_eneg_by_period( period - 1) * ratio_last_two_periods;
          }
          if( ave_vch_hard_by_period( period).GetWeight() == float( 0.0))
          {
            float ratio_last_two_periods( ave_vch_hard_by_period( period - 1) / ave_vch_hard_by_period( period - 2));
            ave_vch_hard_by_period( period) += ave_vch_hard_by_period( period - 1) * ratio_last_two_periods;
          }
        }

        // estimate from period-3 types to period 2 types
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen, GetVchargeTypeString( 3, 1, 0, ' ', 0))] =
          EstimatePeriod2TypeParamsFromPeriod3Params( storage::Pair< float , float >( 33.0, 86.6));
        map[ std::make_pair( chemistry::GetElementTypes().e_Oxygen,   GetVchargeTypeString( 3, 0, 0, ' ', 0))] =
          EstimatePeriod2TypeParamsFromPeriod3Params( storage::Pair< float , float >( 31.8, 93.9));
        map[ std::make_pair( chemistry::GetElementTypes().e_Nitrogen, GetVchargeTypeString( 3, 1, 0, ' ', 1))] =
          EstimatePeriod2TypeParamsFromPeriod3ParamsV2( storage::Pair< float , float >( 33.0, 86.6), chemistry::GetElementTypes().e_Nitrogen);
        map[ std::make_pair( chemistry::GetElementTypes().e_Oxygen,   GetVchargeTypeString( 3, 0, 0, ' ', 1))] =
          EstimatePeriod2TypeParamsFromPeriod3ParamsV2( storage::Pair< float , float >( 31.8, 93.9), chemistry::GetElementTypes().e_Oxygen);

        // estimate period 3 params from period 2 params given in the table
        for
        (
          storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >::const_iterator
          itr( map.Begin()), itr_end( map.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the element type
          chemistry::ElementType element_type( itr->first.first);
          if( element_type->GetPeriod() != size_t( 2))
          {
            continue;
          }

          // move to the corresponding period 3 element by adding 8 (8 elements in both period 2 and 3)
          element_type = chemistry::ElementType( size_t( element_type.GetIndex() + 8));

          std::pair< chemistry::ElementType, std::string> new_key( element_type, itr->first.second);
          // see if the value was already in the map
          if( element_type != itr->first.first && !map.Has( new_key))
          {
            // nope, so estimate the value from the non-aromatic type
            storage::Pair< float, float> eneg_hardness
            (
              EstimatePeriod3TypeParamsFromPeriod2Params( itr->second)
            );

            // only add the values to the map if both are above zero
            if( std::min( eneg_hardness.First(), eneg_hardness.Second()) > float( 0.0))
            {
              map[ new_key] = eneg_hardness;
            }
            new_key.second[ new_key.second.size() - 1] = '1';
            eneg_hardness = EstimatePeriod3TypeParamsFromPeriod2ParamsV2( itr->second, itr->first.first);
            map[ new_key] = eneg_hardness;
          }
        }

        // estimate aromatic values for missing entries in the above tables
        for
        (
          storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >::const_iterator itr( map.Begin()), itr_end( map.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the vcharge type string
          std::string vcharge_string( itr->first.second);

          // check whether the special feature is just a blank
          if( vcharge_string[ vcharge_string.size() - 2] != ' ')
          {
            continue;
          }

          // ensure that this is type has more than one bond by looking for the possible combinations of 1 bonds
          if( vcharge_string.find( "00") != std::string::npos || vcharge_string.find( "010") != std::string::npos)
          {
            continue;
          }

          // change it to aromatic
          vcharge_string[ vcharge_string.size() - 2] = 'A';

          std::pair< chemistry::ElementType, std::string> new_key( itr->first.first, vcharge_string);

          // see if the value was already in the map
          if( !map.Has( new_key))
          {
            // nope, so estimate the value from the non-aromatic type
            storage::Pair< float, float> eneg_hardness( EstimateAromaticTypeParams( itr->second));

            // only add the values to the map if both are above zero
            if( std::min( eneg_hardness.First(), eneg_hardness.Second()) > float( 0.0))
            {
              map[ new_key] = eneg_hardness;
            }
            new_key.second[ vcharge_string.size() - 1] = '1';
            eneg_hardness = EstimateAromaticTypeParamsV2( itr->second, itr->first.first);
            // only add the values to the map if both are above zero
            if( std::min( eneg_hardness.First(), eneg_hardness.Second()) > float( 0.0))
            {
              map[ new_key] = eneg_hardness;
            }
          }
        }

        // Add the averages to the map; for element type use undefined
        for( size_t period( 0), max_period( 7); period < max_period; ++period)
        {
          map[ std::make_pair( chemistry::ElementType(), util::Format()( period + 1) + '0')] =
            storage::Pair< float, float>( ave_vch_eneg_by_period( period), ave_vch_hard_by_period( period));
          map[ std::make_pair( chemistry::ElementType(), util::Format()( period + 1) + '1')] =
            storage::Pair< float, float>( ave_vch_eneg_by_period( std::min( period, size_t( 2))), ave_vch_hard_by_period( std::min( period, size_t( 2))));
        }
      }
      std::ostringstream output;
      for
      (
        storage::Map< std::pair< chemistry::ElementType, std::string>, storage::Pair< float, float> >::const_iterator
          itr( map.Begin()), itr_end( map.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string no_space_str( itr->first.second);
        if( no_space_str[ no_space_str.size() - 2] == ' ')
        {
          no_space_str[ no_space_str.size() - 2] = '_';
        }
        output << itr->first.first->GetChemicalSymbol() << ' ' << no_space_str << ' ' << itr->second.First() << ' ' << itr->second.Second() << '\n';
      }
      // BCL_Debug( output.str());

      return map;
    } // MakeVchargeTypeStringToVchargeParametersMap

  } // namespace descriptor
} // namespace bcl
