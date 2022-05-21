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
#include "assemble/bcl_assemble_protein_model_moment_of_inertia.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_model_moment_of_inertia.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date Apr 16, 2011
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinModelMomentOfInertia :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinModelMomentOfInertia *Clone() const
    {
      return new ExampleAssembleProteinModelMomentOfInertia( *this);
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
      const std::string s_pdbids[] = { "1lgh"};
      const size_t s_number_pdbs( 1);

      storage::Map< std::string, util::ShPtr< assemble::ProteinModelMomentOfInertia> > moment_calculators;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // iterate over transfer free energies
      for( size_t data( biol::AATypeData::e_TransferFreeEnergyWhimleyWhite); data <= biol::AATypeData::e_TransferFreeEnergyPuntaMaritan3D; ++data)
      {
        moment_calculators[ biol::AATypeData::GetPropertyDescriptor( biol::AATypeData::PropertyType( data)) + "_nopositive_nv"] =
          util::ShPtr< assemble::ProteinModelMomentOfInertia>
          (
            new assemble::ProteinModelMomentOfInertia
            (
              biol::AATypeData::PropertyType( data),
              util::ShPtr< assemble::AAExposureInterface>( new assemble::AANeighborVector()),
              true
            )
          );
        moment_calculators[ biol::AATypeData::GetPropertyDescriptor( biol::AATypeData::PropertyType( data)) + "_alsopositive_nv"] =
          util::ShPtr< assemble::ProteinModelMomentOfInertia>
          (
            new assemble::ProteinModelMomentOfInertia
            (
              biol::AATypeData::PropertyType( data),
              util::ShPtr< assemble::AAExposureInterface>( new assemble::AANeighborVector()),
              false
            )
          );
        moment_calculators[ biol::AATypeData::GetPropertyDescriptor( biol::AATypeData::PropertyType( data)) + "_nopositive_no"] =
          util::ShPtr< assemble::ProteinModelMomentOfInertia>
          (
            new assemble::ProteinModelMomentOfInertia
            (
              biol::AATypeData::PropertyType( data),
              util::ShPtr< assemble::AAExposureInterface>(),
              true
            )
          );
        moment_calculators[ biol::AATypeData::GetPropertyDescriptor( biol::AATypeData::PropertyType( data)) + "_alsopositive_no"] =
          util::ShPtr< assemble::ProteinModelMomentOfInertia>
          (
            new assemble::ProteinModelMomentOfInertia
            (
              biol::AATypeData::PropertyType( data),
              util::ShPtr< assemble::AAExposureInterface>(),
              true
            )
          );
      }

      // iterate over all pdb ids
      for( const std::string *pdb_id( s_pdbids), *pdb_id_end( s_pdbids + s_number_pdbs); pdb_id != pdb_id_end; ++pdb_id)
      {
        // iterate over moment_calculators
        for( storage::Map< std::string, util::ShPtr< assemble::ProteinModelMomentOfInertia> >::const_iterator itr( moment_calculators.Begin()), itr_end( moment_calculators.End()); itr != itr_end; ++itr)
        {
          BCL_MessageStd( "processing: " + *pdb_id);
          // construct protein model
          const std::string bio_pdb_filename( AddExampleInputPathToFilename( ExampleInterface::e_Biology, *pdb_id + ".pdb"));
          assemble::ProteinModel model_bio( Proteins::GetModel( bio_pdb_filename));

        /////////////////
        // data access //
        /////////////////

        ///////////////
        // operators //
        ///////////////

        ////////////////
        // operations //
        ////////////////

          // transformation
          const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> transformation_moments( itr->second->TransformationAndMoments( model_bio));
          BCL_MessageStd( "transformation:\t" + *pdb_id + itr->first + "\t" + util::Format()( transformation_moments.First()));
          BCL_MessageStd( "moments:\t" + *pdb_id + itr->first + "\t" + util::Format()( transformation_moments.Second()));

          // transform the model
          model_bio.Transform( transformation_moments.First());

//          // write model to file
//          Proteins::WriteModelToPDB( model_bio, *pdb_id + itr->first + ".pdb");
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinModelMomentOfInertia

  const ExampleClass::EnumType ExampleAssembleProteinModelMomentOfInertia::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinModelMomentOfInertia())
  );

} // namespace bcl
