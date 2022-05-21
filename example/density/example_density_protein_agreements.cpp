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
#include "density/bcl_density_protein_agreements.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_protein_agreements.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensityProteinAgreements :
    public ExampleInterface
  {
  public:

    ExampleDensityProteinAgreements *Clone() const
    { return new ExampleDensityProteinAgreements( *this);}

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

      const double resolution( 6.9);
      const linal::Vector3D grid_spacing( resolution / 3);

      const std::string example_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // preparation
      // create pdb model
      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 0;
      sse_min_size[ biol::GetSSTypes().STRAND] = 0;
      sse_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein( Proteins::GetModel( example_pdb, biol::GetAAClasses().e_AAComplete, sse_min_size));

      // all atoms
      const util::SiPtrVector< const biol::Atom> atoms( protein.GetAtoms());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // no constructors available, since only singleton class

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( density::GetProteinAgreements().GetClassIdentifier(), GetStaticClassName< density::ProteinAgreements>());

    ////////////////
    // operations //
    ////////////////

      // iterator over simulators and generate density maps
      for( density::Simulators::const_iterator sim_map_itr( density::GetSimulators().Begin()), sim_map_itr_end( density::GetSimulators().End()); sim_map_itr != sim_map_itr_end; ++sim_map_itr)
      {
        if( !( *sim_map_itr)->IsDefined())
        {
          continue;
        }

        util::ShPtr< density::Map> sp_map_sim( new density::Map());
        std::string mrc_filename( AddExampleInputPathToFilename( e_Density, "1ubi.pdb"));
        mrc_filename = mrc_filename.substr( 0, mrc_filename.find( ".")) + sim_map_itr->GetName() + ".mrc";

        // try to read the map
        io::IFStream read;
        if( !io::DirectoryEntry( mrc_filename).DoesExist())
        {
          // skip if the example file is missing
          continue;
        }

        BCL_ExampleMustOpenBinaryInputFile( read, mrc_filename);

        // read map
        sp_map_sim->ReadMRC( read);
        io::File::CloseClearFStream( read);

        // iterate over agreement scores
        for( density::ProteinAgreements::const_iterator agr_itr( density::GetProteinAgreements().Begin()), agr_itr_end( density::GetProteinAgreements().End()); agr_itr != agr_itr_end; ++agr_itr)
        {
          // create agreement score
          util::ShPtr< density::ProteinAgreementInterface> sp_agreement( density::GetProteinAgreements().CreateProteinAgreement( *agr_itr, *sim_map_itr, sp_map_sim, resolution));

          if( !sp_agreement.IsDefined())
          {
            continue;
          }

          // calculate agreement
          const double agreement( sp_agreement->operator ()( protein));

          // report agreement
          BCL_MessageStd
          (
            "map simulated " + sim_map_itr->GetName() + " " + agr_itr->GetName() + " agreement:" + util::Format()( agreement)
          );
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

      // protein agreements cannot be written or read - except via the enum read and write mechanism

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDensityProteinAgreements

  const ExampleClass::EnumType ExampleDensityProteinAgreements::s_Instance
  (
    GetExamples().AddEnum( ExampleDensityProteinAgreements())
  );

} // namespace bcl
