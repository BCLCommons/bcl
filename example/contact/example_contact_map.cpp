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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "contact/bcl_contact_map.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_map.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactMap :
    public ExampleInterface
  {
  public:

    ExampleContactMap *Clone() const
    {
      return new ExampleContactMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // initialize pdb handler
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));

      // read the protein model from pdb handler using pdb factory
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // call default constructor
      BCL_MessageStd( "Creating a contact::Map using default constructor");
      contact::Map map_default;

      // constructor from chain
      BCL_MessageStd( "Creating a contact::Map using a chain");
      contact::Map map_from_chain( protein_model.GetChain( 'A'), 5);

      // constructor from protein model
      BCL_MessageStd( "Creating a contact::Map using a protein model");
      contact::Map map_from_model( protein_model, 5);

      // copy constructor
      BCL_MessageStd( "Creating a contact::Map using copy constructor");
      contact::Map map_from_chain_copy( map_from_chain);

      // clone function
      BCL_MessageStd( "Creating a contact::Map using clone");
      util::ShPtr< contact::Map> sp_map( map_from_chain.Clone());

    /////////////////
    // data access //
    /////////////////

      // the class identifier
      BCL_MessageStd
      (
        "This is the class identifier for this class " + map_default.GetClassIdentifier()
      );
      // check the identifier
      BCL_ExampleCheck( map_default.GetClassIdentifier(), "bcl::contact::Map");

      // get the chains
      BCL_MessageStd
      (
        "Number of chains stored" + util::Format()( map_from_chain.GetChains().GetSize())
      );
      // check the number of chains
      BCL_Example_Check
      (
        map_from_chain.GetChains().GetSize() == 1,
        "There should be only one chain in the map but instead it has " +
          util::Format()( map_from_chain.GetChains().GetSize())
      );

      // check that chain A exists
      BCL_ExampleAssert( map_from_chain.GetChain( 'A').IsDefined(), true);

      // check that chain A exists
      BCL_ExampleCheck( map_from_chain.GetChain( 'A')->GetChainID(), 'A');

      // get a specific chain
      BCL_MessageStd
      (
        "Sequence of the chain A" + map_from_chain.GetChain( 'A')->GetSequence()->Sequence()
      );

      // get the boundary
      BCL_MessageStd
      (
        "The boundary value for the map is " + util::Format()( map_from_chain.GetBoundary())
      )
      // check the boundary value
      BCL_ExampleAssert( map_from_chain.GetBoundary(), 5);

    ////////////////
    // operations //
    ////////////////

      // test the GetContactVector function for various residue couples
      BCL_MessageStd( "Printing out contact vectors for different residue couples :");
      storage::Vector< storage::VectorND< 2, size_t> > indices;
      indices.PushBack( storage::VectorND< 2, size_t>( 7, 196));
      indices.PushBack( storage::VectorND< 2, size_t>( 196, 7));
      indices.PushBack( storage::VectorND< 2, size_t>( 13, 189));
      indices.PushBack( storage::VectorND< 2, size_t>( 189, 13));
      indices.PushBack( storage::VectorND< 2, size_t>( 10, 192));
      indices.PushBack( storage::VectorND< 2, size_t>( 192, 10));

      // iterate over index pairs
      for
      (
        storage::Vector< storage::VectorND< 2, size_t> >::const_iterator index_itr( indices.Begin()),
          index_itr_end( indices.End());
        index_itr != index_itr_end; ++index_itr
      )
      {
        // initialize ptr to amino acids
        const util::SiPtr< const biol::AABase> aa_ptr_a
        (
          protein_model.GetChain( 'A')->GetSequence()->GetAA( index_itr->First())
        );
        const util::SiPtr< const biol::AABase> aa_ptr_b
        (
          protein_model.GetChain( 'A')->GetSequence()->GetAA( index_itr->Second())
        );

        // form aa pair
        const storage::VectorND< 2, util::SiPtr< const biol::AAData> > aa_pair
        (
          aa_ptr_a->GetData(),
          aa_ptr_b->GetData()
        );

        // get the predictions
        const storage::Pair< linal::Vector< size_t>, bool> contact_vector( map_from_chain.GetContactVector( aa_pair));

        // output the predictions
        BCL_MessageStd
        (
          util::Format().W( 5)( aa_ptr_a->GetSeqID()) + " " +
          aa_ptr_a->GetType()->GetOneLetterCode() + " " +
          util::Format().W( 5)( aa_ptr_b->GetSeqID()) + " " +
          aa_ptr_b->GetType()->GetOneLetterCode() + " " +
          util::Format().W( 5)( contact_vector.First()( contact::GetTypes().HELIX_HELIX)) + " " +
          util::Format().W( 5)( contact_vector.First()( contact::GetTypes().HELIX_SHEET)) + " " +
          util::Format().W( 5)( contact_vector.First()( contact::GetTypes().SHEET_HELIX)) + " " +
          util::Format().W( 5)( contact_vector.First()( contact::GetTypes().STRAND_STRAND)) + " " +
          util::Format().W( 5)( contact_vector.First()( contact::GetTypes().SHEET_SHEET)) + " " +
          util::Format().W( 5)( contact_vector.Second())
        );

        // output if they are in contact
        BCL_MessageStd
        (
          "Are these two residues in helix-helix contact "
            + util::Format()( contact::Map::IsInContact( *aa_ptr_a, *aa_ptr_b, contact::GetTypes().HELIX_HELIX))
        );

        // check the value in different contact maps initialized
        BCL_Example_Check
        (
          map_from_chain.GetContactVector( aa_pair) == map_from_model.GetContactVector( aa_pair),
          "The contact values differ between map from chain and map from model\n" +
          util::Format()( map_from_chain.GetContactVector( aa_pair)) + " \nvs\n" +
          util::Format()( map_from_model.GetContactVector( aa_pair))
        );
        BCL_Example_Check
        (
          map_from_chain.GetContactVector( aa_pair) == map_from_chain_copy.GetContactVector( aa_pair),
          "The contact values differ between map from chain and map from chain copy\n" +
          util::Format()( map_from_chain.GetContactVector( aa_pair)) + " \nvs\n" +
          util::Format()( map_from_chain_copy.GetContactVector( aa_pair))
        );
        BCL_Example_Check
        (
          map_from_chain.GetContactVector( aa_pair) == sp_map->GetContactVector( aa_pair),
          "The contact values differ between map from chain and cloned map\n" +
          util::Format()( map_from_chain.GetContactVector( aa_pair)) + " \nvs\n" +
          util::Format()( sp_map->GetContactVector( aa_pair))
        );
      } // end index itr

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // initialize output
      io::OFStream write;
      io::IFStream read;
      const std::string output_filename( AddExampleOutputPathToFilename( map_default, "1IE9.contact_map"));
      // open the output file
      BCL_MessageStd( "writing the contact map to file: " + output_filename);
      BCL_ExampleMustOpenOutputFile( write, output_filename);

      // write the map to the file
      map_from_chain.WriteMap( write, true);
      io::File::CloseClearFStream( write);

      // create a new map
      contact::Map map_read_from_file( protein_model.GetChain( 'A'));

      BCL_ExampleMustOpenInputFile( read, output_filename);

      // read map from the previously written contact map
      map_read_from_file.ReadMap( read);

      // now iterate over previous indices and compare values again
      for
      (
        storage::Vector< storage::VectorND< 2, size_t> >::const_iterator index_itr( indices.Begin()),
          index_itr_end( indices.End());
        index_itr != index_itr_end; ++index_itr
      )
      {
        // initialize ptr to amino acids
        const util::SiPtr< const biol::AABase> aa_ptr_a
        (
          protein_model.GetChain( 'A')->GetSequence()->GetAA( index_itr->First())
        );
        const util::SiPtr< const biol::AABase> aa_ptr_b
        (
          protein_model.GetChain( 'A')->GetSequence()->GetAA( index_itr->Second())
        );

        // form aa pair
        const storage::VectorND< 2, util::SiPtr< const biol::AAData> > aa_pair
        (
          aa_ptr_a->GetData(),
          aa_ptr_b->GetData()
        );

        // get the predictions
        const storage::Pair< linal::Vector< size_t>, bool> contact_vector( map_from_chain.GetContactVector( aa_pair));

        // check the value in different contact maps initialized
        BCL_Example_Check
        (
          map_from_chain.GetContactVector( aa_pair) == map_read_from_file.GetContactVector( aa_pair),
          "The contact values differ between map from chain and map_read_from_file\n" +
          util::Format()( map_from_chain.GetContactVector( aa_pair)) + " \nvs\n" +
          util::Format()( map_read_from_file.GetContactVector( aa_pair))
        );
      }

      // reset read
      io::File::CloseClearFStream( read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleContactMap

  const ExampleClass::EnumType ExampleContactMap::s_Instance
  (
    GetExamples().AddEnum( ExampleContactMap())
  );

} // namespace bcl

