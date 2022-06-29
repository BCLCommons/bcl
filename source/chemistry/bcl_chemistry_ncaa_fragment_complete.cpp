// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users: 
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users: 
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University: 
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the 
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University. 
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org 
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include <chemistry/bcl_chemistry_fragment_complete.h>
#include <chemistry/bcl_chemistry_rotamer_library_file.h>
#include <io/bcl_io_file.h>
#include <linal/bcl_linal_matrix3x3.h>
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_ncaa_fragment_complete.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> NCAAFragmentComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new NCAAFragmentComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    NCAAFragmentComplete::NCAAFragmentComplete()
    {
    }

    //! @brief Construct a NCAA from input a input NCAA structure
    //!        Automatically determines the type of the backbone based on the input NCAA
    //!        Should only be used when the backbone of the input NCAA are complete
    NCAAFragmentComplete::NCAAFragmentComplete
    (
      const FragmentComplete& INPUT_NCAA
    )
    {

    }

    //! @brief Construct a NCAA from input a input NCAA structure and a give type of the backbone
    NCAAFragmentComplete::NCAAFragmentComplete
    (
      const FragmentComplete& INPUT_NCAA,
      const Backbone_Type BACKBONE_TYPE
    ) : m_BackboneType( BACKBONE_TYPE)
    {

    }

    //! @brief Clone function
    //! @return pointer to new NCAAFragmentComplete
    NCAAFragmentComplete *NCAAFragmentComplete::Clone() const
    {
      return new NCAAFragmentComplete( *this);
    }
    ///////////////
    // helpers   //
    ///////////////

    //! @brief return the chi1 atom indices
    //! @param NCAA: the atom vector of NCAA
    //! @param BB_CONN_INDICES: the indices of three bb atoms that are closest to Chi 1 Atom
    //! @return the chi1 atom index
    const size_t NCAAFragmentComplete::FindChi1Index
    (
      const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
      storage::Vector< size_t> BB_CONN_INDICES
    ) const
    {
      // iterate over bonds of carbon alpha
      for
      (
          auto bond_itr( NCAA( BB_CONN_INDICES( 0) ).GetBonds().Begin()),
          bond_itr_end( NCAA( BB_CONN_INDICES( 0) ).GetBonds().End());
          bond_itr != bond_itr_end;
          ++bond_itr
      )
      {
        // skip backbone atoms
        const size_t &atom_index( NCAA.GetAtomIndex( bond_itr->GetTargetAtom()));
        if( atom_index == BB_CONN_INDICES( 1) || atom_index == BB_CONN_INDICES( 2) )
        {
          continue;
        }

        // if the next atom is a heavy atom then accept it as chi1;
        // i realize that if there are two heavy atoms then this is not perfect;
        // R/S will be assigned below a
        if( bond_itr->GetTargetAtom().GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
        {
          return atom_index;
        }
      }
      return util::GetUndefinedSize_t();
    }

    ///////////////
    // Alpha AAs //
    ///////////////

    //! @brief load neutral glycine residue from library
    //! @return the neutral glycine as the ncaa base
    const storage::Pair< bool, chemistry::FragmentComplete> NCAAFragmentComplete::ReadGlycineBase() const
    {
      // Begin
      chemistry::FragmentEnsemble glycine;

      // Read in neutral glycine file
      io::IFStream file;
      io::File::MustOpenIFStream
      (
        file,
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/glycine_bb.sdf.gz"
        //chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/glycine_neutral.sdf.gz"
      );
      glycine.ReadMoreFromMdl( file, sdf::e_Maintain);

      // return the glycine residue
      io::File::CloseClearFStream( file);

      // get the ncaa base (which will provide our reference to Rosetta for peptide backbone atoms)
      chemistry::FragmentComplete ncaa_base( glycine.GetMolecules().FirstElement());

      // make sure no one changed the reference file in a harmful way
      const chemistry::AtomVector< chemistry::AtomComplete> atom_v( ncaa_base.GetAtomVector());
      if
      (
          // central amide oxygen atom
          atom_v( 4).GetElementType() != chemistry::GetElementTypes().e_Oxygen ||

          // central amide nitrogen atom
          atom_v( 1).GetElementType() != chemistry::GetElementTypes().e_Nitrogen ||

          // alpha and beta carbon atoms
          atom_v( 3).GetElementType() != chemistry::GetElementTypes().e_Carbon ||
          atom_v( 2).GetElementType() != chemistry::GetElementTypes().e_Carbon ||

          // alpha and beta carbon hydrogen atoms
          atom_v( 6).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 7).GetElementType() != chemistry::GetElementTypes().e_Hydrogen
      )
      {
        return storage::Pair< bool, chemistry::FragmentComplete>
        (
          std::make_pair( bool( false), chemistry::FragmentComplete( ncaa_base))
        );
      }
      return storage::Pair< bool, chemistry::FragmentComplete>
      (
        std::make_pair( bool( true), chemistry::FragmentComplete( ncaa_base))
      );
    }

    //! @brief load neutral glycine dipeptide from library as backbone for alpha AA
    //! @return the neutral glycine dipeptide as the ncaa base
    const storage::Triplet< bool, size_t, chemistry::FragmentComplete> NCAAFragmentComplete::ReadGlycineDipeptideBackbone() const
    {
      // Begin
      chemistry::FragmentEnsemble backbone;
      // Read in neutral glycine file
      io::IFStream file;
      io::File::MustOpenIFStream
      (
        file,
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/glycine_neutral.sdf.gz"
      );
      backbone.ReadMoreFromMdl( file, sdf::e_Maintain);

      // return the glycine residue
      io::File::CloseClearFStream( file);

      // get the ncaa base (which will provide our reference to Rosetta for peptide backbone atoms)
      chemistry::FragmentComplete alpha_bb( backbone.GetMolecules().FirstElement());

      // make sure no one fucked with the reference file in a harmful way
      const chemistry::AtomVector< chemistry::AtomComplete> atom_v( alpha_bb.GetAtomVector());
      if
      (
          //central amide oxygen atom
          atom_v( 5).GetElementType() != chemistry::GetElementTypes().e_Oxygen ||

          // central amide nitrogen atom
          atom_v( 3).GetElementType() != chemistry::GetElementTypes().e_Nitrogen ||

          // alpha and beta carbon atoms
          atom_v( 4).GetElementType() != chemistry::GetElementTypes().e_Carbon ||
          atom_v( 2).GetElementType() != chemistry::GetElementTypes().e_Carbon ||

          // alpha and beta carbon hydrogen atoms
          atom_v( 0).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 1).GetElementType() != chemistry::GetElementTypes().e_Hydrogen
      )
      {
        return storage::Triplet< bool, size_t, chemistry::FragmentComplete>
        (
          bool( false), size_t(), chemistry::FragmentComplete()
        );
      }
      return storage::Triplet< bool, size_t, chemistry::FragmentComplete>
      (
        bool( true), size_t( 4), chemistry::FragmentComplete( alpha_bb)
      );
    }

    //! @brief Find the CA chirality for alpha + Nmethyl alpha NCAAs
    //! @param NCAA: the atom vector of NCAA
    //! @param C_INDEX: the index of backbone C atom
    //! @param N_INDEX: the index of the backbone N atom
    //! @param CHI1_INDEX: the index of the chi1 angle atom
    const std::string NCAAFragmentComplete::FindAlphaCAChirarity
    (
      const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
      const size_t &CA_INDEX,
      const size_t &C_INDEX,
      const size_t &N_INDEX,
      const size_t &CHI1_INDEX
    ) const
    {
      std::string ca_chirarity;
      if
      (
          GetChiralityName( NCAA[ CA_INDEX]->GetChirality()) == "R" ||
          GetChiralityName( NCAA[ CA_INDEX]->GetChirality()) == "S"
      )
      {
        // make a vector of N, C, chi1 atoms.
        // If this is L-AA, the order of those three atoms will be counter clockwise
        // To determine the order of those atoms in space, we can determine the sign of the determinant
        // of the matrix formed by the 3D coordinate of those three points
        // For the explanation, see https://math.stackexchange.com/questions/2400911/ordering-points-in-mathbbr3-using-the-sign-of-determinant
        storage::Vector< size_t> ca_branches( storage::Vector< size_t>::Create( N_INDEX, C_INDEX, CHI1_INDEX));
        // Create a matrix that will hold the x,y, and z coordinates for the three atoms of highest priority.
        linal::Matrix3x3< float> xyz_coordinates( 0.0); // make a matrix of size 3 X 3
        const linal::Vector3D &root_position( NCAA[ CA_INDEX]->GetPosition());

        // put the positions (relative to the root atom) of the 3 highest priority substituents into a matrix
        // the sign of the determinant of the matrix will give us clock-wise (D-AA) or counter-clockwise (L-AA)
        for( size_t index( 0), size( 3); index < size; ++index)
        {
          // get the atom position out of the vector, which has been sorted by priority
          const linal::Vector3D &position( NCAA[ ca_branches( index)]->GetPosition());

          // put the positions of this atom relative to the root atom into the rows of the matrix
          xyz_coordinates( index, 0) = position( 0) - root_position( 0);
          xyz_coordinates( index, 1) = position( 1) - root_position( 1);
          xyz_coordinates( index, 2) = position( 2) - root_position( 2);
        }

        // Calculate the determinant of a 3 X 3 matrix that has rows sorted in descending order of ca_branches.
        // Opposite orders will have opposite signs.
        const float determinant( xyz_coordinates.Determinant());
        //BCL_Debug( xyz_coordinates.Determinant());
        if( determinant > 0.0)
        {
          ca_chirarity = "L_AA";
        }
        else if( determinant < 0.0)
        {
          ca_chirarity = "D_AA";
        }
        else
        {
          BCL_MessageStd( "Warning: cannot determine whether the NCAA is L or D stereoisomer");
          ca_chirarity = "ACHIRAL_BACKBONE";
        }
      }
      else
      {
        BCL_MessageStd( "Warning: cannot determine whether the NCAA is L or D stereoisomer");
        ca_chirarity = "ACHIRAL_BACKBONE";
      }
      //BCL_Debug( ca_chirarity);
      return ca_chirarity;
    }

    ///////////////
    // peptoid   //
    ///////////////

    //! @brief Find the CA chirality for alpha + Nmethyl alpha NCAAs
    //! @param NCAA: the atom vector of NCAA
    //! @param C_INDEX: the index of backbone C atom
    //! @param N_INDEX: the index of the backbone N atom
    //! @param CHI1_INDEX: the index of the chi1 angle atom
    const std::string NCAAFragmentComplete::FindPeptoidNChirarity
    (
      const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
      const size_t &C_INDEX,
      const size_t &N_INDEX,
      const size_t &CHI1_INDEX
    ) const
    {

    }

    ////////////////////
    // N methyl alpha //
    ////////////////////

    //! @brief load neutral glycine residue from library
    //! @return the neutral glycine as the ncaa base
    const storage::Pair< bool, chemistry::FragmentComplete> NCAAFragmentComplete::ReadNMethylGlycineBase() const
    {
      // Begin
      chemistry::FragmentEnsemble nmethyl_glycine;

      // Read in neutral glycine file
      io::IFStream file;
      io::File::MustOpenIFStream
      (
        file,
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/nmethyl_glycine_bb.sdf.gz"
      );
      nmethyl_glycine.ReadMoreFromMdl( file, sdf::e_Maintain);

      // return the glycine residue
      io::File::CloseClearFStream( file);

      // get the ncaa base (which will provide our reference to Rosetta for peptide backbone atoms)
      chemistry::FragmentComplete ncaa_base( nmethyl_glycine.GetMolecules().FirstElement());

      // make sure no one change the reference file in a harmful way
      const chemistry::AtomVector< chemistry::AtomComplete> atom_v( ncaa_base.GetAtomVector());
      if
      (
          // methyl group
          atom_v( 0).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 2).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 3).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 4).GetElementType() != chemistry::GetElementTypes().e_Carbon ||

          // amide oxygen atom
          atom_v( 1).GetElementType() != chemistry::GetElementTypes().e_Oxygen ||

          // amide nitrogen atom
          atom_v( 8).GetElementType() != chemistry::GetElementTypes().e_Nitrogen ||

          // alpha carbon atom
          atom_v( 9).GetElementType() != chemistry::GetElementTypes().e_Carbon ||

          // alpha carbon hydrogen atoms
          atom_v( 5).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 6).GetElementType() != chemistry::GetElementTypes().e_Hydrogen
      )
      {
        return storage::Pair< bool, chemistry::FragmentComplete>
        (
          std::make_pair( bool( false), chemistry::FragmentComplete( ncaa_base))
        );
      }
      return storage::Pair< bool, chemistry::FragmentComplete>
      (
        std::make_pair( bool( true), chemistry::FragmentComplete( ncaa_base))
      );
    }

    //! @brief load neutral glycine dipeptide from library as backbone for alpha AA
    //! @return the neutral glycine dipeptide as the ncaa base
    const storage::Triplet< bool, size_t, chemistry::FragmentComplete> NCAAFragmentComplete::ReadNMethylGlycineDipeptideBackbone() const
    {
      // Begin
      chemistry::FragmentEnsemble backbone;
      // Read in neutral glycine file
      io::IFStream file;
      io::File::MustOpenIFStream
      (
        file,
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/nmethyl_glycine_neutral.sdf.gz"
      );
      backbone.ReadMoreFromMdl( file, sdf::e_Maintain);

      // return the glycine residue
      io::File::CloseClearFStream( file);

      // get the ncaa base (which will provide our reference to Rosetta for peptide backbone atoms)
      chemistry::FragmentComplete bb( backbone.GetMolecules().FirstElement());

      // make sure no one changed the reference file in a harmful way
      const chemistry::AtomVector< chemistry::AtomComplete> atom_v( bb.GetAtomVector());
      if
      (
          //central amide oxygen atom
          atom_v( 5).GetElementType() != chemistry::GetElementTypes().e_Oxygen ||

          // central amide nitrogen atom
          atom_v( 3).GetElementType() != chemistry::GetElementTypes().e_Nitrogen ||

          // alpha and beta carbon atoms
          atom_v( 4).GetElementType() != chemistry::GetElementTypes().e_Carbon ||
          atom_v( 2).GetElementType() != chemistry::GetElementTypes().e_Carbon ||

          // alpha and beta carbon hydrogen atoms
          atom_v( 0).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 1).GetElementType() != chemistry::GetElementTypes().e_Hydrogen
      )
      {
        return storage::Triplet< bool, size_t, chemistry::FragmentComplete>
        (
          bool( false), size_t(), chemistry::FragmentComplete()
        );
      }
      return storage::Triplet< bool, size_t, chemistry::FragmentComplete>
      (
        bool( true), size_t( 4), chemistry::FragmentComplete( bb)
      );
    }


  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &NCAAFragmentComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NCAAFragmentComplete::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &NCAAFragmentComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
