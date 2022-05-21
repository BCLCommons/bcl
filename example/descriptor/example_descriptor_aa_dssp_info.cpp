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
#include "descriptor/bcl_descriptor_aa_dssp_info.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_dssp_info.cpp
  //!
  //! @author mendenjl
  //! @date May 12, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAADSSPInfo :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAADSSPInfo *Clone() const
    {
      return new ExampleDescriptorAADSSPInfo( *this);
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
      // construct an implementation that will read proteins from the biol directory
      const std::string pdb_fasta_dir( AddExampleInputPathToFilename( e_Biology, ""));

      // create an object to read proteins, as either pdbs or fastas, from the directory
      util::Implementation
      <
        assemble::RetrieveProteinModelWithCache
      > retriever_all( "SequenceDirectory( " + pdb_fasta_dir + ")");

      // retrieve 1b43A
      assemble::ProteinModelWithCache protein( *retriever_all->Retrieve( "1b43A"));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////
      // create a descriptor that will retrieve the accessibility
      descriptor::AADSSPInfo accessibility( descriptor::AADSSPInfo::e_Accessibility);
      // call set object with protein
      accessibility.SetObject( protein);

      // create a descriptor that will retrieve the total energy
      descriptor::AADSSPInfo total_energy( descriptor::AADSSPInfo::e_TotalEnergy);
      // call set object with protein
      total_energy.SetObject( protein);

      // create a descriptor that will retrieve the maximum energy of any neighbor
      descriptor::AADSSPInfo max_hbond_potential( descriptor::AADSSPInfo::e_MaxHBondPotential);
      // call set object with protein
      max_hbond_potential.SetObject( protein);

      // create a descriptor that will retrieve the aa seqid offset of the neighbor with the best energetic interaction with this aa
      descriptor::AADSSPInfo max_hbond_neighbor_offset( descriptor::AADSSPInfo::e_MaxHBondNeighborOffset);
      // call set object with protein
      max_hbond_neighbor_offset.SetObject( protein);

      // clone
      util::ShPtr< descriptor::AADSSPInfo> sp_aa_info( accessibility.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( accessibility.GetSizeOfFeatures(), 1);

    ///////////////
    // operators //
    ///////////////

      descriptor::Iterator< biol::AABase> itr( accessibility.GetType(), protein);
      ++itr;

      // operator()
      BCL_ExampleCheck
      (
        accessibility( itr),
        linal::Vector< float>( size_t( 1), 31)
      );
      itr.GotoPosition( 63);
      BCL_ExampleCheck
      (
        accessibility( itr),
        linal::Vector< float>( size_t( 1), 47)
      );

      descriptor::Iterator< biol::AABase> itr_id( total_energy.GetType(), protein);
      ++itr_id;
      // operator()
      BCL_ExampleCheck
      (
        total_energy( itr_id),
        linal::Vector< float>( size_t( 1), -0.4)
      );
      itr_id.GotoPosition( 63);
      BCL_ExampleCheck
      (
        total_energy( itr_id),
        linal::Vector< float>( size_t( 1), -4.5)
      );

      descriptor::Iterator< biol::AABase> itr_pos( max_hbond_potential.GetType(), protein);
      ++itr_pos;
      // operator()
      BCL_ExampleCheck
      (
        max_hbond_potential( itr_pos),
        linal::Vector< float>( size_t( 1), -0.2)
      );
      itr_pos.GotoPosition( 64);
      BCL_ExampleCheck
      (
        max_hbond_potential( itr_pos),
        linal::Vector< float>( size_t( 1), -2.1)
      );

      descriptor::Iterator< biol::AABase> itr_dist( max_hbond_neighbor_offset.GetType(), protein);
      ++++itr_dist;
      // operator()
      BCL_ExampleCheck
      (
        max_hbond_neighbor_offset( itr_dist),
        linal::Vector< float>( size_t( 1), 2)
      );
      itr_dist.GotoPosition( 64);
      BCL_ExampleCheck
      (
        max_hbond_neighbor_offset( itr_dist),
        linal::Vector< float>( size_t( 1), 4)
      );
    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAADSSPInfo

  const ExampleClass::EnumType ExampleDescriptorAADSSPInfo::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAADSSPInfo())
  );

} // namespace bcl
