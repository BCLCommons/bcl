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
#include "biol/bcl_biol_aa_side_chain_factory.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_aa_complete.h"
#include "io/bcl_io_file.h"
#include "quality/bcl_quality_rmsd.h"
#include "score/bcl_score.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASideChainFactory::s_Instance
    (
      GetObjectInstances().AddInstance( new AASideChainFactory())
    );

    //! @brief filename of default atom coordinate table
    const std::string &AASideChainFactory::GetDefaultTableFileName()
    {
      // initialize static table filename
      static const std::string s_filename( "idealized_side_chain_coords.table");

      // end
      return s_filename;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //!@brief default constructor
    AASideChainFactory::AASideChainFactory() :
      m_IncludeHydrogen(),
      m_IncludeBackBone(),
      m_TableFilename()
    {
    }

    //! @brief default constructor
    //! @param INCLUDE_HYDROGEN = true, if hydrogen atoms will be included in attached side chains
    //! @param INCLUDE_BACKBONE = use the backbone to align sidechain for first chi angle (important for proline)
    //! @param TABLE_FILE_NAME filename to table which contains ideal side chain conformations
    AASideChainFactory::AASideChainFactory
    (
      bool INCLUDE_HYDROGEN,
      const bool INCLUDE_BACKBONE,
      const std::string &TABLE_FILE_NAME //= GetDefaultTableFileName()
     ) :
      m_IncludeHydrogen( INCLUDE_HYDROGEN),
      m_IncludeBackBone( INCLUDE_BACKBONE),
      m_TableFilename( TABLE_FILE_NAME)
    {
      InitializeTable();
    }

    //! @brief virtual copy constructor
    AASideChainFactory *AASideChainFactory::Clone() const
    {
      return new AASideChainFactory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASideChainFactory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add side chains to AABase
    //! @param AA_BASE object to which side chain atoms should be added
    void AASideChainFactory::AddSideChain( AABase &AA_BASE) const
    {
      // Generate idealized AA sidechain from AABase's type
      util::ShPtrVector< Atom> side_chain_atoms( GenerateAtoms( AA_BASE.GetType()));

      // determine the transformation matrix to put sidechains in correct place
      math::TransformationMatrix3D transformation;

      if( m_IncludeBackBone)
      {
        transformation = GetQuadTransformationMatrix( AA_BASE, side_chain_atoms);
      }
      else
      {
        transformation = GetTransformationMatrix( AA_BASE, side_chain_atoms);
      }

      // check if transformation is defined
      if( !transformation.IsDefined())
      {
        BCL_MessageStd
        (
          "undefined transformation matrix - cannot attach sidechain atoms with defined coordinates"
        );

        // iterate over all sidechain atoms and render coordiantes undefined (nan)
        for
        (
          util::ShPtrVector< Atom>::iterator
            this_itr( side_chain_atoms.Begin()), end_itr( side_chain_atoms.End());
          this_itr != end_itr; ++this_itr
        )
        {
          ( *this_itr)->SetCoordinates( linal::Vector3D( util::GetUndefined< double>()));
        }

        // set the atoms
        AA_BASE.SetAtoms( side_chain_atoms);

        // end
        return;
      }

      // transform the atoms
      TransformAtoms( transformation, side_chain_atoms);

      if( m_IncludeBackBone)
      {
        // storage for new side chain atoms
        util::ShPtrVector< Atom> side_chain_atoms_new;

        // iterate over the atoms in side_chain_atoms
        for
        (
          util::ShPtrVector< Atom>::iterator
            this_itr( side_chain_atoms.Begin()), end_itr( side_chain_atoms.End());
          this_itr != end_itr; ++this_itr
        )
        {
          // if not a backbone atom type
          if( !( *this_itr)->GetType()->IsBackBone())
          {
            // add it to side chain atoms
            side_chain_atoms_new.PushBack( *this_itr);
          }
        }
        // update the side chain atoms
        side_chain_atoms = side_chain_atoms_new;
      }

      // set the atoms
      AA_BASE.SetAtoms( side_chain_atoms);
    }

    // takes a Protein model and returns a new protein model with side chains attached to new amino acids, sses and chains
    util::ShPtr< assemble::ProteinModel> AASideChainFactory::ProteinModelWithSideChains
    (
      const assemble::ProteinModel &THIS_MODEL
    ) const
    {
      // make ShPtrVector for the new chains
      util::ShPtrVector< assemble::Chain> new_chains;

      // iterate over chains in THIS_MODEL, add side chains to the AAs in it, and push it back into new_chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( THIS_MODEL.GetChains().Begin()), chain_end_itr( THIS_MODEL.GetChains().End());
        chain_itr != chain_end_itr;
        ++chain_itr
      )
      {
        new_chains.PushBack( NewChain( **chain_itr));
      }

      // construct a new ProteinModel and return it
      util::ShPtr< assemble::ProteinModel> new_model( new assemble::ProteinModel( new_chains));
      return new_model;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASideChainFactory::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_IncludeHydrogen, ISTREAM);
      io::Serialize::Read( m_IncludeBackBone, ISTREAM);
      io::Serialize::Read( m_Table, ISTREAM);
      io::Serialize::Read( m_TableFilename, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AASideChainFactory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_IncludeHydrogen, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IncludeBackBone, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Table, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TableFilename, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize table with idealized side chain coordinates
    void AASideChainFactory::InitializeTable()
    {
      io::IFStream read;
      io::File::MustOpenIFStream( read, score::Score::AddHistogramPath( m_TableFilename));
      m_Table.ReadFormatted( read);
      io::File::CloseClearFStream( read);
    }

    //! @brief takes a chain and returns a new chain with side chains attached to new amino acids
    //! @param THIS_CHAIN the chain which contains AAs with no side chains
    //! @return ShPtr to new chain with side chains on the AAs
    util::ShPtr< assemble::Chain> AASideChainFactory::NewChain( const assemble::Chain &THIS_CHAIN) const
    {
      util::ShPtrVector< assemble::SSE> new_sses;
      // iterate over chain and Get SSEs.
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( THIS_CHAIN.GetData().Begin()), sse_end_itr( THIS_CHAIN.GetData().End());
        sse_itr != sse_end_itr; ++sse_itr
      )
      {
        new_sses.PushBack( NewSSE( **sse_itr));
      }
      util::ShPtr< assemble::Chain> new_chain( new assemble::Chain( THIS_CHAIN.GetSequence(), new_sses));
      return new_chain;
    }

    //! @brief similar to assemble::SSE::HardCopy, but adds SideChain
    //! @param THIS_SSE the SSE to be copied and have side chains added
    //! @return ShPtr to new SSE with side chains
    util::ShPtr< assemble::SSE> AASideChainFactory::NewSSE( const assemble::SSE &THIS_SSE) const
    {
      // is the SSE considered a coil?
      //const bool sse_is_coil( THIS_SSE.GetType() == biol::GetSSTypes().COIL);

      // storage for the new amino acids with side chains attached
      util::ShPtrVector< AABase> new_aabase_vector;
      for
      (
        AASequence::const_iterator
          aa_itr( THIS_SSE.Begin()), aa_end_itr( THIS_SSE.End());
        aa_itr != aa_end_itr; ++aa_itr
      )
      {
         // create a new amino acid from this amino acid with side chains
         new_aabase_vector.PushBack( NewAAWithSideChain( **aa_itr));
      }

      AASequence new_sequence( new_aabase_vector, THIS_SSE.GetChainID(), THIS_SSE.GetFastaHeader());
      util::ShPtr< assemble::SSE> new_sse( new assemble::SSE( new_sequence, THIS_SSE.GetType()));
      return new_sse;
    }

    //! @brief adds a side chain to an individual amino acid
    //! @param THIS_AA the amino acid that should have a side chain attached to it
    //! @return ShPtr to an AAComplete with sidechain attached
    util::ShPtr< AABase> AASideChainFactory::NewAAWithSideChain( const AABase &THIS_AA) const
    {
      // make sure that THIS_AA is at least AABackBone or AAComplete
      BCL_Assert
      (
        THIS_AA.GetAAClass() == GetAAClasses().e_AABackBone ||
        THIS_AA.GetAAClass() == GetAAClasses().e_AAComplete,
        "AA Could not have side chain added because it was of type: " +
        util::Format()( THIS_AA.GetAAClass())
      );

      // create new ShPtr to new AAComplete, and add a sidechain to it
      util::ShPtr< AABase> new_aa_shptr( new AAComplete( THIS_AA));
      AddSideChain( *new_aa_shptr);
      return new_aa_shptr;
    }

    //! @brief pull the atoms for a given aa type from the table with their idealized coordinates
    //! @brief basically this is the idealized amino acid copy in space
    //! @param AATYPE the AAType considered
    util::ShPtrVector< Atom> AASideChainFactory::GenerateAtoms( const AAType &AATYPE) const
    {
      // create storage for atoms
      util::ShPtrVector< Atom> atoms;

      // grab the row corresponding to the requested AAType
      const storage::Row< double> &aa_row( m_Table[ AATYPE->GetThreeLetterCode()]);

      // get set of atom types for the AATYPE
      const storage::Set< AtomType> this_atomtypes( AATYPE->GetAllowedAtomTypes());

      // iterate over all the atom types
      for
      (
        storage::Set< AtomType>::const_iterator itr( this_atomtypes.Begin()), itr_end( this_atomtypes.End());
        itr != itr_end;
        ++itr
      )
      {
        // if hydrogens are excluded and the atom type is a hydrogen and the aatype is glycine, still include the HA
        if( AATYPE == GetAATypes().GLY && !m_IncludeHydrogen && ( *itr) == GetAtomTypes().HA)
        {
          //Add the HA for glycine
          Atom this_atom( GenerateAtomWithCoords( *itr, aa_row));

          // insert atom into m_RotamAtoms
          atoms.PushBack( util::ShPtr< Atom>( this_atom.Clone()));

          // go to next atom type
          continue;
        }

        // if hydrogens are not to be included
        if( !m_IncludeHydrogen && ( *itr)->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          continue;
        }

        // exclude backbone if desired
        if( !m_IncludeBackBone && ( *itr)->IsBackBone())
        {
          continue;
        }

        //generate atom for that atomtype using default coordinates from the table
        Atom this_atom( GenerateAtomWithCoords( *itr, aa_row));

        // insert atom into m_RotamAtoms
        atoms.PushBack( util::ShPtr< Atom>( this_atom.Clone()));
      }

      return atoms;
    }

    //! @brief transforms ShPtrVector< Atom> according to transformation matrix
    //! @param TRANSFORMATION matrix used to transform atoms
    //! @param ATOMS atoms to be transformed
    void AASideChainFactory::TransformAtoms
    (
      const math::TransformationMatrix3D &TRANSFORMATION, util::ShPtrVector< Atom> &ATOMS
    )
    {
      // iterate over the atoms in ATOMS transforming each atom
      for( util::ShPtrVector< Atom>::iterator itr( ATOMS.Begin()), end_itr( ATOMS.End()); itr != end_itr; ++itr)
      {
        ( *itr)->Transform( TRANSFORMATION);
      }
    }

    //! @brief get transformation matrix to superimpose CA-CB vector of amino acid with the CA-CB vector of the atoms
    //! @param AMINO_ACID the amino acid that is to be used
    //! @param ATOMS the atoms generated from the table
    //! @return transformationmatrix the transformation matrix used to transform the atoms of the rotamer
    math::TransformationMatrix3D AASideChainFactory::GetTransformationMatrix
    (
      const AABase &AMINO_ACID,
      const util::ShPtrVector< Atom> &ATOMS
    ) const
    {
      // make SiPtrVector of resultant for rotamer Ca-Cb Vector3D
      util::SiPtrVector< linal::Vector3D> rotam_siptr_vec;
      const linal::Vector3D &aa_ca_coords( AMINO_ACID.GetAtom( GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D &aa_cb_coords( AMINO_ACID.GetAtom( GetAtomTypes().CB).GetCoordinates());
      linal::Vector3D rotam_ca_coords;
      linal::Vector3D rotam_cb_coords;

      // iterate over atoms contained in ATOMS (rotamer atoms)
      for( util::ShPtrVector< Atom>::const_iterator itr( ATOMS.Begin()), end_itr( ATOMS.End()); itr != end_itr; ++itr)
      {
        //find the CA and get the coordinates
        if( ( *itr)->GetType() == GetAtomTypes().CA)
        {
          rotam_ca_coords = ( ( *itr)->GetCoordinates());
        }

        //find the CB and get the coordinates
        if( ( *itr)->GetType() == GetAtomTypes().CB)
        {
          rotam_cb_coords = ( ( *itr)->GetCoordinates());
        }
      }
      // define the transformation matrix
      math::TransformationMatrix3D matrix;

      // apply transformation to matrix with the negative of the rotam_ca_coords so that it translates to the origin
      matrix( -rotam_ca_coords);

      // determine the projection angle for Ca-Cb vector alignment
      const double projangle( linal::ProjAngle( rotam_ca_coords, rotam_cb_coords, aa_ca_coords, aa_cb_coords));

      // instantiate the rotation matrix with the axis of rotation (via cross product) and the projection angle
      const math::RotationMatrix3D rotationmatrix
      (
        linal::CrossProduct( aa_ca_coords, aa_cb_coords, rotam_ca_coords, rotam_cb_coords), projangle
      );

      // apply transformation matrix with the rotation matrix from above
      matrix( rotationmatrix);

      // add the aa_ca_coord vector to the transformation matrix so that the Ca-Cb is properly aligned
      // note that this math:Vector3D is implicitly converted to math::TransformationMatrix3D
      matrix( aa_ca_coords);

      // return the matrix that will apply the transformation to a set of un-transformed atoms
      return matrix;
    }

    //! @brief get transformation matrix to superimpose C-CA-CB-N vector of amino acid with the C-CA-CB-N vector of the atoms
    //! @param AMINO_ACID the amino acid that is to be used for superimposition
    //! @param ATOMS the atoms generated from the table
    //! @return transformation matrix the transformation matrix used to transform the atoms of the rotamer
    math::TransformationMatrix3D AASideChainFactory::GetQuadTransformationMatrix
    (
      const AABase &AMINO_ACID,
      const util::ShPtrVector< Atom> &ATOMS
    ) const
    {
      // initialize static set of desired atom types
      static const storage::Set< AtomType> desired_atoms
      (
        storage::Set< AtomType>::Create
        (
          GetAtomTypes().CA, GetAtomTypes().CB, GetAtomTypes().N, GetAtomTypes().C
        )
      );

      // make SiPtrVector of backbone atoms in given amino acid
      util::SiPtrVector< const linal::Vector3D> aa_coords( AMINO_ACID.GetAtomCoordinates( desired_atoms));

      // make SiPtrVector of backbone atoms in rotamer atoms
      util::SiPtrVector< const linal::Vector3D> rotam_siptr_vec;

      // iterate over backbone atom types
      for
      (
        storage::Set< AtomType>::const_iterator
          aatype_itr( desired_atoms.Begin()), aatype_itr_end( desired_atoms.End());
        aatype_itr != aatype_itr_end; ++aatype_itr
      )
      {
        // iterate over atoms contained in ATOMS (rotamer atoms)
        for( util::ShPtrVector< Atom>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End()); itr != itr_end; ++itr)
        {
          if( *aatype_itr == ( *itr)->GetType())
          {
            rotam_siptr_vec.PushBack( util::ToSiPtr( ( *itr)->GetCoordinates()));
            break;
          }
        }
      }
      // superimpose coordiantes and return the superimposition matrix
      return quality::RMSD::SuperimposeCoordinates( aa_coords, rotam_siptr_vec);
    }

    //! @brief a helper function which generates atoms from given atom type with its ideal coordinates from table
    //! @param ATOM_TYPE atom type to be search table for
    //! @param ROW searches the row (predetermined AA type) for the atom name of choice
    //! @return biol::Atom with coordinates from table
    Atom AASideChainFactory::GenerateAtomWithCoords
    (
      const AtomType &ATOM_TYPE,
      const storage::Row< double> &ROW
    )
    {
      // initialize atom names
      const std::string atom_x_name( ATOM_TYPE.GetName() + "_x");
      const std::string atom_y_name( ATOM_TYPE.GetName() + "_y");
      const std::string atom_z_name( ATOM_TYPE.GetName() + "_z");

      //get the coordinates based on atom name given
      const linal::Vector3D coord
      (
        ROW[ atom_x_name],
        ROW[ atom_y_name],
        ROW[ atom_z_name]
      );

      // coordinates should be defined
      BCL_Assert( coord.IsDefined(), "Some coordinates for this atom type are not defined!");

      // construct and return atom
      return Atom( coord, ATOM_TYPE);
    }

  } // namespace biol
} // namespace bcl
