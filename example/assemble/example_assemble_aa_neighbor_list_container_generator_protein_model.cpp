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
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_aa_neighbor_list_container_generator_protein_model.cpp
  //!
  //! @author karakam
  //! @date Oct 13, 2011
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleAANeighborListContainerGeneratorProteinModel :
    public ExampleInterface
  {
  public:

    ExampleAssembleAANeighborListContainerGeneratorProteinModel *Clone() const
    {
      return new ExampleAssembleAANeighborListContainerGeneratorProteinModel( *this);
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
      // read in model
      assemble::ProteinModel native( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb")));

      const size_t nr_repeats( 2);
      util::Stopwatch timer;

      // initialize expected counts
      const size_t expected_nr_center_aas( 259);
      const size_t expected_full_nr_neighbors( 64770);
      const size_t expected_pruned_nr_neighbors( 4228);
      const size_t expected_ssemove_nr_neighbors( 4228);

    ////////////////////////////////
    // test with full constructor //
    ////////////////////////////////
      BCL_MessageStd( "testing full constructor");
      {
        for( size_t i( 1); i <= nr_repeats; ++i)
        {
          timer.Reset();
          timer.Start();
          assemble::AANeighborListContainer container( native.GetAminoAcids(), 10000, 0, true);
          timer.Stop();
          const size_t nr_aas( container.GetSize());
          const size_t nr_neighbors( container.GetNumberNeighbors());
          BCL_MessageStd
          (
            "#" + util::Format()( i) + " time: " + util::Format()( timer.GetTotalTime().GetTotalMicroseconds()) +
            "\taas: " + util::Format()( nr_aas) + "\tnr_neigh: " + util::Format()( nr_neighbors)
          );
          BCL_ExampleCheck( nr_aas, expected_nr_center_aas);
          BCL_ExampleCheck( nr_neighbors, expected_full_nr_neighbors);
          timer.Reset();
        }
      }

    //////////////////////////////
    // test with full generator //
    //////////////////////////////
      BCL_MessageStd( "testing full generator");
      {
        // generator
        util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > sp_generator
        (
          assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator( 10000, 0, true, false)
        );

        for( size_t i( 1); i <= nr_repeats; ++i)
        {
          timer.Reset();
          timer.Start();
          assemble::AANeighborListContainer container( sp_generator->operator ()( native));
          timer.Stop();
          const size_t nr_aas( container.GetSize());
          const size_t nr_neighbors( container.GetNumberNeighbors());
          BCL_MessageStd
          (
            "#" + util::Format()( i) + " time: " + util::Format()( timer.GetTotalTime().GetTotalMicroseconds()) +
            "\taas: " + util::Format()( nr_aas) + "\tnr_neigh: " + util::Format()( nr_neighbors)
          );
          BCL_ExampleCheck( nr_aas, expected_nr_center_aas);
          BCL_ExampleCheck( nr_neighbors, expected_full_nr_neighbors);
          timer.Reset();
        }
      }

    //////////////////////////////////
    // test with pruned constructor //
    //////////////////////////////////
      BCL_MessageStd( "testing pruned constructor");
      {
        for( size_t i( 1); i <= nr_repeats; ++i)
        {
          timer.Reset();
          timer.Start();
          assemble::AANeighborListContainer container( native.GetAminoAcids(), 12.0, 6, true);
          timer.Stop();
          const size_t nr_aas( container.GetSize());
          const size_t nr_neighbors( container.GetNumberNeighbors());
          BCL_MessageStd
          (
            "#" + util::Format()( i) + " time: " + util::Format()( timer.GetTotalTime().GetTotalMicroseconds()) +
            "\taas: " + util::Format()( nr_aas) + "\tnr_neigh: " + util::Format()( nr_neighbors)
          );
          BCL_ExampleCheck( nr_aas, expected_nr_center_aas);
          BCL_ExampleCheck( nr_neighbors, expected_pruned_nr_neighbors);
          timer.Reset();
        }
      }

      // store time without cache
      size_t time_without_cache( 0);

    ////////////////////////////////
    // test with pruned generator //
    ////////////////////////////////
      BCL_MessageStd( "testing pruned generator");
      {
        // generator
        util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > sp_generator
        (
          assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator( 12.0, 6, true, false)
        );
        for( size_t i( 1); i <= nr_repeats; ++i)
        {
          timer.Reset();
          timer.Start();
          assemble::AANeighborListContainer container( sp_generator->operator ()( native));
          timer.Stop();
          const size_t nr_aas( container.GetSize());
          const size_t nr_neighbors( container.GetNumberNeighbors());
          const size_t this_time( timer.GetTotalTime().GetTotalMicroseconds());
          BCL_MessageStd
          (
            "#" + util::Format()( i) + " time: " + util::Format()( this_time) +
            "\taas: " + util::Format()( nr_aas) + "\tnr_neigh: " + util::Format()( nr_neighbors)
          );
          // store the time
          if( i == 1) { time_without_cache = this_time;}
          BCL_ExampleCheck( nr_aas, expected_nr_center_aas);
          BCL_ExampleCheck( nr_neighbors, expected_pruned_nr_neighbors);
          timer.Reset();

        }
      }

      // define a translation
      linal::Vector3D translation( 0.05, 0.05, 0.0);

    ////////////////////////////////////////////
    // test with constructor after moving SSE //
    ////////////////////////////////////////////
      BCL_MessageStd( "testing constructor after moving SSE");
      {
        // initialize a new model
        util::ShPtr< assemble::ProteinModel> model_ptr( native.HardCopy());
        assemble::ProteinModel &model( *model_ptr);

        for( size_t i( 1); i <= nr_repeats; ++i)
        {
          // make a hardcopy of the SSE that will be moved
          util::ShPtr< assemble::SSE> sp_sse( assemble::LocatorSSE( 'A', 8, 26).Locate( model)->HardCopy());
          //sp_sse->Translate( translation);
          model.Replace( sp_sse);
          timer.Reset();
          timer.Start();
          assemble::AANeighborListContainer container( model.GetAminoAcids(), 12.0, 6, true);
          timer.Stop();
          const size_t nr_aas( container.GetSize());
          const size_t nr_neighbors( container.GetNumberNeighbors());
          const size_t this_time( timer.GetTotalTime().GetTotalMicroseconds());
          BCL_MessageStd
          (
            "#" + util::Format()( i) + " time: " + util::Format()( this_time) +
            "\taas: " + util::Format()( nr_aas) + "\tnr_neigh: " + util::Format()( nr_neighbors)
          );
          // if it's the second run
          BCL_ExampleCheck( nr_aas, expected_nr_center_aas);
          BCL_ExampleCheck( nr_neighbors, expected_ssemove_nr_neighbors);
          timer.Reset();
        }
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleAANeighborListContainerGeneratorProteinModel

  const ExampleClass::EnumType ExampleAssembleAANeighborListContainerGeneratorProteinModel::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleAANeighborListContainerGeneratorProteinModel())
  );

} // namespace bcl
