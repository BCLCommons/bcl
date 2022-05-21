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
#include "chemistry/bcl_chemistry_bond_dihedral_angles.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_bond_conformational.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    namespace
    {
      // load the dihedral map
      storage::Map< storage::VectorND< 7, size_t>, storage::Pair< double, double> > LoadDihedralMap()
      {
        storage::Map< storage::VectorND< 7, size_t>, storage::Pair< double, double> > loaded_map;
        io::IFStream input;
        io::File::MustOpenIFStream( input, RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "dihedral.sds.txt.gz");
//        io::File::MustOpenIFStream( input, RotamerLibraryFile::GetRotamerFinder().FindFile( "dihedral.sds.txt.gz"));
        util::ChopHeader( input);

        std::string bondnames[ 8] = { "1xChain", "2xChain", "3xChain", "Amide", "Aromatic", "1xRing", "2xRing", "3xRing"};
        std::map< std::string, size_t> string_to_int;
        for( int i( 0); i < 7; ++i)
        {
          string_to_int[ bondnames[ i]] = size_t( i + 1);
        }

        std::string line;
        std::getline( input, line);
        while( input.good() && std::getline( input, line))
        {
          auto vec( util::SplitString( line, " \t"));
          AtomType a( vec( 0)), b( vec( 2)), c( vec( 4)), d( vec( 6));
          std::string ab( vec( 1)), bc( vec( 3)), cd( vec( 5));
          storage::VectorND< 7, size_t> dihedral;
          dihedral( 0) = a;
          dihedral( 2) = b;
          dihedral( 4) = c;
          dihedral( 6) = d;
          dihedral( 1) = string_to_int[ ab];
          dihedral( 3) = string_to_int[ bc];
          dihedral( 5) = string_to_int[ cd];
          double cnt( util::ConvertStringToNumericalValue< double>( vec( 7)));
          if( cnt < 20)
          {
            continue;
          }
          double ave( util::ConvertStringToNumericalValue< double>( vec( 8)));
          double sd( util::ConvertStringToNumericalValue< double>( vec( 9)));
          loaded_map[ dihedral] = storage::Pair< double, double>( sd, ave);
          std::swap( dihedral( 0), dihedral( 6));
          std::swap( dihedral( 1), dihedral( 5));
          std::swap( dihedral( 2), dihedral( 4));
          loaded_map[ dihedral] = storage::Pair< double, double>( sd, ave);
        }
        io::File::CloseClearFStream( input);
        return loaded_map;
      }
    }

    //! @brief estimate the standard deviation and deviation from perfect (0 degree-centered) bond for a dihedral given only the two atoms
    //! @param ATOM_B, ATOM_C Atoms in the conformation that are bonded
    //! @return estimated standard deviation for the dihedral about this bond
    double BondDihedralAngles::GetEstimatedStdForDihedralBondAngleBin
    (
      const AtomConformationalInterface &ATOM_B,
      const AtomConformationalInterface &ATOM_C
    )
    {
      double min_sd( 10.0);
      // get the type of bond between them
      ConfigurationalBondType bond_bc;
      for( auto itr_b( ATOM_B.GetBonds().Begin()), itr_b_end( ATOM_B.GetBonds().End()); itr_b != itr_b_end; ++itr_b)
      {
        if( &itr_b->GetTargetAtom() == &ATOM_C)
        {
          bond_bc = itr_b->GetBondType();
          break;
        }
      }
      if( bond_bc->IsBondInRing())
      {
        // skip ring bonds
        return 0.0;
      }
      if( bond_bc->GetConjugation() == ConstitutionalBondTypeData::e_Amide)
      {
        return 4.5;
      }
      for( auto itr_b( ATOM_B.GetBonds().Begin()), itr_b_end( ATOM_B.GetBonds().End()); itr_b != itr_b_end; ++itr_b)
      {
        if( &itr_b->GetTargetAtom() == &ATOM_C)
        {
          continue;
        }

        const AtomConformationalInterface &atom_a( itr_b->GetTargetAtom());
        const ConfigurationalBondType &bond_ab( itr_b->GetBondType());

        for( auto itr_c( ATOM_C.GetBonds().Begin()), itr_c_end( ATOM_C.GetBonds().End()); itr_c != itr_c_end; ++itr_c)
        {
          if( &itr_c->GetTargetAtom() == &ATOM_B)
          {
            continue;
          }

          const AtomConformationalInterface &atom_d( itr_c->GetTargetAtom());
          const ConfigurationalBondType &bond_cd( itr_c->GetBondType());
          storage::Pair< double, double> ave_sd
          (
            GetAveStdDihedralDeviation
            (
              atom_a.GetAtomType(),
              bond_ab,
              ATOM_B.GetAtomType(),
              bond_bc,
              ATOM_C.GetAtomType(),
              bond_cd,
              atom_d.GetAtomType()
            )
          );
          min_sd = std::min( min_sd, ave_sd.First());
        }
      }
      return min_sd;
    }

    //! @brief Get the average and standard deviation of a given dihedral angle
    //! @param ATOM_TYPE)A,ATOM_TYPE_B,ATOM_TYPE_C, ATOM_TYPE_D the atom types in the bond
    //! @param BOND_TYPE_AB, BOND_TYPE_BC, BOND_TYPE_CD the bond types in the bond
    //! @return return the estimated ave/std dihedral deviation for a centered bin w/ 30 degree window
    //!         Values will be undefined (nan) if unavailable for the diehdral
    storage::Pair< double, double> BondDihedralAngles::GetAveStdDihedralDeviation
    (
      const AtomType &ATOM_TYPE_A,
      const ConfigurationalBondType &BOND_TYPE_AB,
      const AtomType &ATOM_TYPE_B,
      const ConfigurationalBondType &BOND_TYPE_BC,
      const AtomType &ATOM_TYPE_C,
      const ConfigurationalBondType &BOND_TYPE_CD,
      const AtomType &ATOM_TYPE_D
    )
    {
      if( BOND_TYPE_BC->IsBondInRing())
      {
        return storage::Pair< double, double>( 0.0, 0.0);
      }
      static const storage::Map< storage::VectorND< 7, size_t>, storage::Pair< double, double> > s_map( LoadDihedralMap());
      storage::VectorND< 7, size_t> dihedral;
      dihedral( 0) = ATOM_TYPE_A.GetIndex();
      dihedral( 2) = ATOM_TYPE_B.GetIndex();
      dihedral( 4) = ATOM_TYPE_C.GetIndex();
      dihedral( 6) = ATOM_TYPE_D.GetIndex();
      dihedral( 1) = BOND_TYPE_AB->GetBondData( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness);
      dihedral( 3) = BOND_TYPE_BC->GetBondData( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness);
      dihedral( 5) = BOND_TYPE_CD->GetBondData( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness);
      auto itr( s_map.Find( dihedral));
      if( itr == s_map.End())
      {
        // average SDs for single / non-single bonds
        return storage::Pair< double, double>
               (
                 BOND_TYPE_BC->GetNumberOfElectrons() == size_t( 2) ? double( 7.75) : double( 2.5),
                 double( 0.0)
               );
      }
      return itr->second;
    }

  } // namespace chemistry
} // namespace bcl

